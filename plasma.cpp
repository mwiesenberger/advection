#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "dg/file/file.h"
#include "common.h"

using Vector = std::array<std::array<dg::HVec,2>,2>;

namespace equations
{

//compute **negative** laplacian of Phi
struct Poisson1d
{
    Poisson1d( dg::Grid1d g, Json::Value js) : m_gh_phi( g.size()+4),
        m_hx( g.h())
    {
        m_bcx = dg::str2bc(js["bc"].get("potential", "PER").asString());
        m_eps_D   = js["physical"].get("epsilon_D", 1e-4).asDouble();
        m_mode = js["physical"].get("type", "original").asString();
    }
    void symv( const dg::HVec& phi, dg::HVec& lapMPhi)
    {
        assign_ghost_cells( phi, m_gh_phi, m_bcx);
        unsigned Nx = phi.size();
        for( unsigned i=0; i<Nx; i++)
        {
            unsigned k=i+2;
            lapMPhi[i] = -m_eps_D*(m_gh_phi[k+1] - 2.*m_gh_phi[k] +
                    m_gh_phi[k-1])/m_hx/m_hx;
        }
        if( "adiabatic" == m_mode)
            for( unsigned i=0; i<Nx; i++)
                lapMPhi[i] += exp(phi[i]);
    }
    void operator()( const dg::HVec& phi, dg::HVec& lapMPhi)
     {
         symv( phi, lapMPhi);

     }

    private:
    dg::HVec m_gh_phi;
    dg::bc m_bcx;
    double m_hx;
    double m_eps_D;
    std::string m_mode;
};

struct PlasmaExplicit
{
    //the rationale for making  it a friend is that we logically would
    //implement the implicit part as a method of this class but the
    //interface of the timestepper mandates a separate class
    //It is as if PlasmaImplicit is an extension to this class
    friend class PlasmaImplicit;
    PlasmaExplicit( dg::Grid1d g, dg::Grid1d vel_g, Json::Value js) :
        m_g(g), m_vel_g(vel_g),
        m_poisson( g, js),
        m_old_phi( 2, dg::evaluate( dg::zero, g))
    {
        m_velocity[0].resize( g.size(), 0.);
        m_velocity[1].resize( g.size(), 0.);
        dg::HVec temp( g.size()+4, 0.);
        m_yg[0].fill( temp);
        m_yg[1].fill( temp);
        m_scheme = js["advection"].get("type", "upwind").asString();
        m_variant = js["advection"].get("variant", "original").asString();
        m_eps_D   = js["physical"].get("epsilon_D", 1e-4).asDouble();
        m_mode = js["physical"].get("type", "original").asString();
        m_tau[0] = -1.;
        m_tau[1] = js["physical"].get("tau", 1.0).asDouble();
        m_mu[0] = js["physical"].get("mu", -0.00027244).asDouble();
        m_mu[1] = 1.;
        m_nu_u[0] = js["physical"]["nu_u"].get( 0u, 0.0).asDouble();
        m_nu_u[1] = js["physical"]["nu_u"].get( 1u, 0.0).asDouble();
        m_nu_n[0] = js["physical"]["nu_n"].get( 0u, 0.0).asDouble();
        m_nu_n[1] = js["physical"]["nu_n"].get( 1u, 0.0).asDouble();
        m_eta     = js["physical"].get("resistivity", 1e-4).asDouble();
        m_bc_n = dg::str2bc(js["bc"].get("density", "PER").asString());
        m_bc_u = dg::str2bc(js["bc"].get("velocity", "PER").asString());
        m_bc_p = dg::str2bc(js["bc"].get("potential", "PER").asString());
        std::cout << "Compute scheme "<<m_scheme<<" Variant " << m_variant<<"\n";
        m_init = js["init"].get( "type", "step").asString();
        if( "mms" == m_init)
        {
            m_n0 = js["init"].get("n_0", 1.0).asDouble();
            m_u0 = js["init"].get("u_0", 1.0).asDouble();
            m_A = js["init"].get("A", 1.0).asDouble();
            m_B = js["init"].get("B", 1.0).asDouble();
            m_k = js["init"].get("k", 1.0).asDouble();
            m_v = js["init"].get("v", 1.0).asDouble();
        }
        m_eps = js["poisson"].get("eps", 1e-6).asDouble();
        m_phi.resize(g.size(), 0.);
        m_ghphi.resize(g.size()+4, 0.);
        m_precond.resize( g.size(), 1./m_eps_D);
        m_rhs.resize( g.size(), 0.);
        m_norm.resize( g.size(), g.h());
        m_source.resize( g.size(), 0.);
        m_method = js["poisson"].get( "type", "gmres").asString();
        if( "gmres" == m_method)
        {
            unsigned max_inner = js["poisson"].get("max_inner", 30).asUInt();
            unsigned max_outer = js["poisson"].get("max_outer", 3).asUInt();
            unsigned restarts = g.size()/max_outer;

            m_lgmres.construct( m_phi, max_inner, max_outer, restarts);
        }
        else if( "cg" == m_method)
        {
            m_pcg.construct( m_phi, m_phi.size());
        }
        else if( "bicgstab" == m_method)
        {
            unsigned l_input = js["poisson"].get( "l_input", 2).asUInt();
            m_bicg.construct( m_phi, m_phi.size(), l_input);
        }
        else if( "anderson" == m_method)
        {
            m_mMax = js["poisson"].get( "mMax", 10).asUInt();
            m_damping = js["poisson"].get("damping", 1e-2).asDouble();
            m_anderson.construct( m_phi, m_mMax);
        }

    }

    const std::array<dg::HVec,2>& velocity() const{return m_velocity;}
    const dg::HVec& potential() const{return m_phi;}
    void solve_poisson( double t, const std::array<dg::HVec,2> & nn, dg::HVec& m_ghphi) {
        //initial guess
        m_old_phi.extrapolate( t, m_phi);
        // a little bit dangerous because what you actually need to check is the rhs
        // the solvers will also realize that the guessed solution is close
        //if( m_old_phi.exists( t))
        //{
        //    assign_ghost_cells( m_phi, m_ghphi, m_bc_p);
        //    //std::cout << "Solution of Poisson equation already available\n";
        //    return;
        //}

        unsigned iter;
        if( m_mode == "adiabatic")
            dg::blas1::copy( nn[1], m_rhs);
        else
            dg::blas1::axpby( 1., nn[1], -1., nn[0], m_rhs);
        if( m_init == "mms")
        {
            dg::HVec tmpN( nn[0]);
            double k = m_k, A = m_A, v= m_v;
            tmpN = dg::evaluate( [=](double x){
                return (A*((-1 + m_eps_D*k*k)*cos(k*(-(t*v) + x)) +
                            sin(k*(-(t*v) + x))));
                }, m_g);
            dg::blas1::axpby( 1., tmpN, 1., m_rhs);
        }
        if( m_method == "gmres")
            iter = m_lgmres.solve( m_poisson, m_phi, m_rhs, m_precond,
                    m_norm, m_eps, m_eps_D);
        else if( m_method == "cg")
            iter = m_pcg( m_poisson, m_phi, m_rhs, m_precond,
                    m_norm, m_eps, m_eps_D);
        else if( m_method == "bicgstab")
            iter = m_bicg.solve( m_poisson, m_phi, m_rhs, m_precond,
                    m_norm, m_eps, m_eps_D);
        else if ( m_method == "anderson")
        {
            if( m_mode == "adiabatic")
                dg::blas1::transform( nn[1], m_phi, dg::LN<double>());
            iter = m_anderson.solve( m_poisson, m_phi, m_rhs, m_norm, m_eps,
                    m_eps*m_eps_D, m_phi.size(), m_damping, m_mMax, false);
        }
        if( iter == m_g.N())
            throw std::runtime_error( "Solution of Poisson equation does not converge!\n");
        m_old_phi.update( t, m_phi);
        std::cout << "Solution of Poisson equation took "<<iter<<" iterations\n";
        assign_ghost_cells( m_phi, m_ghphi, m_bc_p);
    }
    void operator() ( double t, const Vector & y, Vector& yp)
    {
        m_called++;
        // y[0] -> density , y[0][0] electrons, y[0][1] ions
        // y[1] -> velocity, y[1][0] electrons, y[1][0] ions
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y[0][0], m_yg[0][0], m_bc_n);
        assign_ghost_cells( y[0][1], m_yg[0][1], m_bc_n);
        assign_ghost_cells( y[1][0], m_yg[1][0], m_bc_u);
        assign_ghost_cells( y[1][1], m_yg[1][1], m_bc_u);
        solve_poisson( t, y[0], m_ghphi);
        dg::Upwind upwind;
        SlopeLimiter<MinMod> limiter;
        for( unsigned s=0; s<2; s++) // species loop
        {
            if( m_mode == "adiabatic" && s == 0)
                continue;
            // ghost cells are shifted by 2
            if ( m_scheme == "centered")
            {
                const dg::HVec& nn = m_yg[0][s];
                const dg::HVec& uu = m_yg[1][s];

                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    double jpar = m_yg[0][1][k]*m_yg[1][1][k] -
                        m_yg[0][0][k]*m_yg[1][0][k];
                    m_velocity[s][i] = uu[k];
                    yp[0][s][i] = -uu[k]*( nn[k+1] - nn[k-1])/2./hx;
                    yp[0][s][i]+= -nn[k]*( uu[k+1] - uu[k-1])/2./hx;
                    yp[1][s][i] =
                        -uu[k]*(uu[k+1]-uu[k-1])/2./hx
                        - m_tau[s]/m_mu[s]*(nn[k+1] -nn[k-1])/2./hx/nn[k]
                        -1./m_mu[s]*(m_ghphi[k+1] - m_ghphi[k-1])/2./hx
                        -m_eta/m_mu[s]*m_yg[0][1][k]/nn[k]*jpar
                        + m_nu_u[s]/nn[k]*(uu[k+1] - 2.*uu[k] + uu[k-1])/hx/hx;
                }
            }
            else if ( m_scheme == "staggered")
            {
                dg::HVec qST(m_yg[1][s]), q(qST), uh(qST), uST( q), dn(q), du(q);
                const dg::HVec & unST = m_yg[1][s];
                const dg::HVec & nn = m_yg[0][s];
                for( unsigned k=0; k<Nx+3; k++)
                {
                    double nST = (nn[k] + nn[k+1])/2.;
                    uST[k] = unST[k]/nST;
                    dn[k] = nn[k+1]-nn[k];
                }
                for( unsigned k=1; k<Nx+2; k++)
                {
                    qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                    if( m_variant == "slope-limiter-explicit" || m_variant ==
                            "slope-limiter")
                        qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                    qST[k]*= uST[k]; // k + 1/2
                }
                for ( unsigned k=1; k<Nx+3; k++)
                    q [k] = 0.5*(qST[k]+qST[k-1]);
                for( unsigned k=1; k<Nx+4; k++)
                    du[k] = uST[k] - uST[k-1];
                for ( unsigned k=2; k<Nx+3; k++)
                {
                    uh[k] = upwind( q[k], uST[k-1], uST[k]);
                    if( m_variant == "slope-limiter-explicit" || m_variant ==
                            "slope-limiter")
                        uh[k] += limiter(q[k], du[k-1], du[k], du[k+1], 0.5, 0.5);
                }
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    m_velocity[s][i] = 0.5*(unST[k]+unST[k-1])/nn[k];
                    yp[0][s][i] = -( qST[k] - qST[k-1])/hx;
                    //
                    yp[1][s][i] = -(uh[k+1]*q[k+1]-uh[k]*q[k])/hx;
                    yp[1][s][i]+= m_nu_u[s]*(uST[k+1] - 2.*uST[k] + uST[k-1])
                        /hx/hx;
                    //resistivity
                    double niST = (m_yg[0][1][k]+m_yg[0][1][k+1])/2.;
                    yp[1][s][i] += -m_eta*niST*( m_yg[1][1][k] -
                                            m_yg[1][0][k] ) / m_mu[s];
                }
                if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
                {
                    for( unsigned i=0; i<Nx; i++)
                    {
                        unsigned k=i+2;
                        double nST = (nn[k] + nn[k+1])/2.;
                        yp[1][s][i] += -m_tau[s]/m_mu[s]*( nn[k+1] -nn[k])/hx ;
                        yp[1][s][i] += -1./m_mu[s]*nST*( m_ghphi[k+1]
                                -m_ghphi[k])/hx;
                    }
                }
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                yp[0][s][i] +=  m_nu_n[s]*( m_yg[0][s][k+1] - 2.*m_yg[0][s][k] +
                        m_yg[0][s][k-1])/hx/hx;
                yp[0][s][i] += m_source[i];
            }
        } //species loop

        if( "mms" == m_init)
        {
            for( unsigned s=0; s<2; s++)
            {
                double k = m_k, A = m_A, B = m_B, n0 = m_n0, u0 = m_u0, v= m_v;
                double mue = m_mu[0], mui = m_mu[1], taue = m_tau[0];
                double taui = m_tau[1], eta = m_eta, nue = m_nu_u[0];
                double nui = m_nu_u[1];
                //Add sources
                dg::HVec tmpN( yp[0][s]), tmpNST( tmpN), tmpUST(tmpN);
                if( s == 0)
                {
                    tmpN = dg::evaluate( [=](double x){
                        return k*cos(k*(-(t*v) + x))*(B*n0 + A*(u0 - v)
                                + 2*A*B*sin(k*(-(t*v) + x)));
                        }, m_g);
                    tmpNST = dg::evaluate( [=](double x){
                        return k*cos(k*(-(t*v) + x))*(B*n0 + A*(u0 - v)
                                + 2*A*B*sin(k*(-(t*v) + x)));
                        }, m_vel_g);
                    tmpUST = dg::evaluate( [=](double x){
                        return -(B*k*v*cos(k*(-(t*v) + x))) -
                        (A*k*sin(k*(-(t*v) + x)))/mue +
                        (B*k*k*nue*sin(k*(-(t*v) + x)))/(n0 +
                                A*sin(k*(-(t*v) + x))) + B*k*cos(k*(-(t*v) +
                                    x))*(u0 + B*sin(k*(-(t*v) + x))) +
                        (A*k*taue*cos(k*(-(t*v) + x)))/(mue*n0 +
                                A*mue*sin(k*(-(t*v) + x))) + (eta*(n0 +
                                    A*cos(k*(-(t*v) + x)))*(cos(k*(-(t*v) + x))
                                    + sin(k*(t*v - x)))*(B*n0 + A*u0 +
                                    A*B*(cos(k*(-(t*v) + x)) + sin(k*(-(t*v) +
                                                x)))))/ (mue*(n0 +
                                        A*sin(k*(-(t*v) + x))));
                        }, m_vel_g);
                }
                if( s == 1)
                {
                    tmpN = dg::evaluate( [=](double x){
                        return -(k*(B*n0 + A*(u0 - v) + 2*A*B*cos(k*(-(t*v) +
                                            x)))*sin(k*(-(t*v) + x))) ;
                        }, m_g);
                    tmpNST = dg::evaluate( [=](double x){
                        return -(k*(B*n0 + A*(u0 - v) + 2*A*B*cos(k*(-(t*v) +
                                            x)))*sin(k*(-(t*v) + x)));
                        }, m_vel_g);
                    tmpUST = dg::evaluate( [=](double x){
                        return (B*k*k*nui)/A + (eta*(B*n0 +
                                    A*u0)*cos(k*(-(t*v) + x)))/mui +
                        (A*B*eta*cos(2*k*(-(t*v) + x)))/mui + ((B*eta*n0 + A*(k
                                        + eta*u0) + B*k*mui*(u0 -
                                            v))*sin(k*(t*v - x)))/mui -
                        (k*(B*k*mui*n0*nui + A*A*taui*sin(k*(-(t*v) +
                                x))))/(A*mui*(n0 + A*cos(k*(-(t*v) + x))))
                        - (B*B*k*sin(2*k*(-(t*v) + x)))/2. ;
                        }, m_vel_g);
                }
                dg::blas1::axpby( 1., tmpN, 1., yp[0][s]);
                if( m_scheme == "staggered")
                {
                    dg::HVec nST(m_yg[0][s]), uST( nST);
                    const dg::HVec & unST = m_yg[1][s];
                    const dg::HVec & nn = m_yg[0][s];
                    for( unsigned k=0; k<Nx+3; k++)
                    {
                        nST[k] = (nn[k] + nn[k+1])/2.;
                        uST[k] = unST[k]/nST[k];
                    }
                    for( unsigned i=0; i<Nx; i++)
                    {
                        unsigned k=i+2;
                        yp[1][s][i] += tmpUST[i]*nST[k] + uST[k]*tmpNST[i];
                    }
                }
                else
                    dg::blas1::axpby( 1., tmpUST, 1., yp[1][s]);
            }

        }
    }
    unsigned called() const { return m_called;}
    private:
    std::string m_scheme, m_variant, m_init, m_method, m_mode;
    Vector m_yg;
    dg::HVec m_precond, m_norm, m_phi, m_ghphi, m_rhs, m_source;
    std::array<dg::HVec,2> m_velocity; // stores the velocity on non-staggered grid
    dg::Grid1d m_g, m_vel_g;
    dg::LGMRES<dg::HVec> m_lgmres;
    dg::CG<dg::HVec> m_pcg;
    dg::BICGSTABl<dg::HVec> m_bicg;
    dg::AndersonAcceleration<dg::HVec> m_anderson;
    Poisson1d m_poisson;
    dg::Extrapolation<dg::HVec> m_old_phi;
    dg::bc m_bc_n, m_bc_u, m_bc_p;
    std::array<double,2> m_tau, m_mu, m_nu_u, m_nu_n;
    double m_eps, m_eta, m_eps_D;
    double m_n0 = 0, m_u0 = 0, m_A = 0, m_B = 0, m_k = 0, m_v = 0;
    double m_mMax, m_damping;
    unsigned m_called = 0;
};

struct PlasmaImplicit
{
    PlasmaImplicit( PlasmaExplicit& ex) : m_ex ( ex) { }

    void operator() ( double t, const Vector & y, Vector& yp)
    {
        dg::blas1::copy( 0., yp);
        if( m_ex.m_scheme == "staggered")
        {
            // we can access everything in PlasmaExplicit
            unsigned Nx = m_ex.m_g.N();
            double hx = m_ex.m_g.h();
            assign_ghost_cells( y[0][0], m_ex.m_yg[0][0], m_ex.m_bc_n);
            assign_ghost_cells( y[0][1], m_ex.m_yg[0][1], m_ex.m_bc_u);
            // ghost cells are shifted by 2
            if( m_ex.m_variant != "explicit" && m_ex.m_variant !=
                    "slope-limiter-explicit")
            {
                m_ex.solve_poisson( t, y[0], m_ex.m_ghphi);
                dg::HVec& ghphi = m_ex.m_ghphi;
                for( unsigned s=0; s<2; s++)
                {
                    if( m_ex.m_mode == "adiabatic" && s == 0)
                        continue;
                    dg::HVec& nn = m_ex.m_yg[0][s];
                    for( unsigned i=0; i<Nx; i++)
                    {
                        unsigned k=i+2;
                        double nST = (nn[k] + nn[k+1])/2.;
                        yp[1][s][i] += -m_ex.m_tau[s]/m_ex.m_mu[s]*( nn[k+1]
                                -nn[k])/hx ;
                        yp[1][s][i] += -1./m_ex.m_mu[s]*nST*( ghphi[k+1] -
                                ghphi[k])/hx;
                    }
                } // species loop
            }
        }
    }
    private:
    PlasmaExplicit& m_ex;
};

struct PlasmaImplicitSolver
{
    PlasmaImplicitSolver( dg::Grid1d g, Json::Value js) :
        m_tmp( {dg::HVec(g.size(), 0.0), dg::HVec ( g.size(), 0.),
                dg::HVec(g.size(), 0.0), dg::HVec ( g.size(), 0.)}){}
    const Vector& copyable() const{
        return m_tmp;
    }
    // solve (y + alpha I(t,y) = rhs
    void solve( double alpha, PlasmaImplicit& im, double t, Vector& y, const Vector& rhs)
    {
        dg::blas1::copy( rhs[0], y[0]);// I_n = 0
        im( t, y, m_tmp); //ignores y[1], solves Poisson at time t for y[0] and
        // writes 0 in m_tmp[0] and updates m_tmp[1]
        dg::blas1::axpby( 1., rhs[1], -alpha, m_tmp[1], y[1]); // u = rhs_u - alpha I_u
    }
    private:
    Vector m_tmp;

};

////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
struct Variables{
    PlasmaExplicit& f;
    const dg::Grid1d& grid;
    const Vector& y0;
    const double& time;
    Json::Value& js;
    double duration;
    unsigned nfailed;
};

struct Record1d{
    std::string name;
    std::string long_name;
    std::function<double( Variables&)> function;
};

struct Record{
    std::string name;
    std::string long_name;
    std::function<void( dg::HVec&, Variables&)> function;
};

std::vector<Record> diagnostics_list = {
    {"electrons", "Numerical electron density",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.y0[0][0], result);
        }
    },
    {"ions", "Numerical electron density",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.y0[0][1], result);
        }
    },
    {"ue", "Numerical electron velocity",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.velocity()[0], result);
        }
    },
    {"ui", "Numerical ion velocity",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.velocity()[1], result);
        }
    },
    {"potential", "potential",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.potential(), result);
        }
    },
    {"electrons_ana", "Analytical solution to the electron density",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double n_0 = v.js["init"].get("n_0", 1.0).asDouble();
                double A = v.js["init"].get("A", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return n_0 + A*sin( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    },
    {"ions_ana", "Analytical solution to the ion density",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double n_0 = v.js["init"].get("n_0", 1.0).asDouble();
                double A = v.js["init"].get("A", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return n_0 + A*cos( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    },
    {"ue_ana", "Analytical solution to the electron velocity",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double u_0 = v.js["init"].get("u_0", 1.0).asDouble();
                double B = v.js["init"].get("B", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return u_0 + B*sin( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    },
    {"ui_ana", "Analytical solution to the ion velocity",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double u_0 = v.js["init"].get("u_0", 1.0).asDouble();
                double B = v.js["init"].get("B", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return u_0 + B*cos( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    },
    {"potential_ana", "Analytical solution to the potential",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double A = v.js["init"].get("A", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return A*cos( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    },
};
std::vector<Record1d> diagnostics1d_list = {
    {"failed", "Accumulated Number of failed steps",
        []( Variables& v ) {
            return v.nfailed;
        }
    },
    {"duration", "Computation time for the latest output",
        []( Variables& v ) {
            return v.duration;
        }
    },
    {"nsteps", "Accumulated Number of calls to the timestepper (including failed steps)",
        [](Variables& v) {
            return v.f.called();
        }
    }
};
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
} //namespace equations

namespace dg
{

// register Poisson1d with dg library
template<>
struct TensorTraits< equations::Poisson1d >
{
    using value_type      = double;
    using tensor_category = SelfMadeMatrixTag;
};
}
int main( int argc, char* argv[])
{
    ////Parameter initialisation ////////////////////////////////////////////
    Json::Value js;
    if( argc == 1)
        dg::file::file2Json( "input/default.json", js, dg::file::comments::are_discarded);
    else
        dg::file::file2Json( argv[1], js);
    std::cout << js <<std::endl;

    /////////////////////////////////////////////////////////////////
    dg::Grid1d grid = equations::createGrid( js["grid"], dg::PER);
    dg::HVec w1d( dg::create::weights(grid));
    /////////////////////////////////////////////////////////////////
    std::string init = js["init"].get("type", "step").asString();
    std::string scheme = js["advection"].get("type", "staggered").asString();
    std::string mode = js["physical"].get("type", "original").asString();
    dg::Grid1d vel_grid = equations::createGrid( js["grid"], dg::PER);
    if ( "staggered" == scheme )
        vel_grid = equations::createStaggeredGrid( js["grid"], dg::PER);
    Vector y0 = {dg::evaluate( dg::zero, grid), dg::evaluate( dg::zero, grid),
        dg::evaluate( dg::zero, grid), dg::evaluate( dg::zero, grid)};
    if( "step" == init)
    {
        // This is the classical Riemann problem (Dam break for shallow water)
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        double n_l = js["init"].get("n_l", 1.0).asDouble();
        double n_r = js["init"].get("n_r", 0.2).asDouble();
        x_a *= grid.lx();
        y0[0][0] = y0[0][1] = dg::evaluate( [=](double x){ return x < x_a ? n_l
                : n_r;}, grid);
    }
    else if( "soft-step" == init)
    {
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        double n_l = js["init"].get("n_l", 1.0).asDouble();
        double n_r = js["init"].get("n_r", 0.2).asDouble();
        double alpha = js["init"].get("alpha", 0.1).asDouble();
        x_a*= grid.lx();
        alpha*=grid.lx();
        dg::PolynomialHeaviside poly( x_a, alpha/2., -1);
        y0[0][0] = y0[0][1] = dg::evaluate( [&](double x){
                return n_r + (n_l-n_r)*poly(x);}, grid);
    }
    else if( "wave" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double amp = js["init"].get("amp", 1.0).asDouble();
        double k   = js["init"].get("k", 1.0).asDouble();
        double x_0 = js["init"].get("x_0", 1.0).asDouble();
        y0[0][0] = y0[0][1] = dg::evaluate( [=]( double x){ return n_0 +
                amp*sin( k*(x-x_0));}, grid);
    }
    else if ( "mms" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double u_0 = js["init"].get("u_0", 1.0).asDouble();
        double A = js["init"].get("A", 1.0).asDouble();
        double B = js["init"].get("B", 1.0).asDouble();
        double k = js["init"].get("k", 1.0).asDouble();
        y0[0][0] = dg::evaluate( [=](double x){ return n_0 + A*sin( k*x);}, grid);
        y0[0][1] = dg::evaluate( [=](double x){ return n_0 + A*cos( k*x);}, grid);
        y0[1][0] =  dg::evaluate( [=](double x){ return u_0 + B*sin( k*x);}, vel_grid);
        y0[1][1] =  dg::evaluate( [=](double x){ return u_0 + B*cos( k*x);}, vel_grid);
        if( scheme == "staggered")
        {
            y0[1][0] = dg::evaluate( [=](double x){ return
                    (u_0+B*sin(k*x))*(n_0+A*sin(k*x));}, vel_grid);
            y0[1][1] = dg::evaluate( [=](double x){ return
                    (u_0+B*cos(k*x))*(n_0+A*cos(k*x));}, vel_grid);
        }
    }
    if( "adiabatic" == mode)
        y0[0][0] = y0[1][0] = dg::evaluate( dg::zero, grid);

    std::string timestepper = js["timestepper"].get( "type", "ERK").asString();
    std::string tableau= js["timestepper"].get("tableau",
            "ARK-4-2-3").asString();
    double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
    double atol= js["timestepper"].get("atol", 1e-6).asDouble();
    double tend = js["output"].get( "tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = tend/(double)maxout;
    dg::ARKStep<Vector, equations::PlasmaImplicitSolver>
        ark_fixed( dg::EXPLICIT_EULER_1_1, dg::IMPLICIT_EULER_1_1, grid, js);
    dg::RungeKutta<Vector> erk_fixed( "Euler", y0);

    dg::Adaptive<dg::ARKStep<Vector,
        equations::PlasmaImplicitSolver>> ark_adaptive( "ARK-4-2-3", grid, js);
    dg::Adaptive<dg::ERKStep<Vector>> erk_adaptive( "Bogacki-Shampine-4-2-3", y0);
    if( timestepper == "ARK")
        ark_adaptive = dg::Adaptive<dg::ARKStep<Vector,
            equations::PlasmaImplicitSolver>>( tableau, grid, js);
    else if( timestepper == "ERK")
    {
        erk_fixed = dg::RungeKutta<Vector>( tableau, y0);
        erk_adaptive = dg::Adaptive<dg::ERKStep<Vector >>( tableau, y0);
    }
    double dt = 1e-8, time = 0.;
    equations::PlasmaExplicit ex( grid, vel_grid, js);
    equations::PlasmaImplicit im( ex);

    // Set up netcdf
    std::string inputfile = js.toStyledString(); //save input without comments, which is important if netcdf file is later read by another parser
    std::string outputfile;
    if( argc == 1 || argc == 2)
        outputfile = "plasma.nc";
    else
        outputfile = argv[2];
    /// //////////////////////set up netcdf/////////////////////////////////////
    dg::file::NC_Error_Handle err;
    int ncid=-1;
    try{
        err = nc_create( outputfile.c_str(),NC_NETCDF4|NC_CLOBBER, &ncid);
    }catch( std::exception& e)
    {
        std::cerr << "ERROR creating file "<<outputfile<<std::endl;
        std::cerr << e.what()<<std::endl;
        return -1;
    }
    /// Set global attributes
    std::map<std::string, std::string> att;
    att["title"] = "Output file of advection/plasma.cpp";
    att["Conventions"] = "CF-1.7";
    ///Get local time and begin file history
    auto ttt = std::time(nullptr);
    auto tm = *std::localtime(&ttt);

    std::ostringstream oss;
    ///time string  + program-name + args
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    for( int i=0; i<argc; i++) oss << " "<<argv[i];
    att["history"] = oss.str();
    att["comment"] = "Find more info in advection/OneDimensional.ipynb";
    att["source"] = "FELTOR";
    att["references"] = "https://github.com/feltor-dev/feltor";
    att["inputfile"] = inputfile;
    for( auto pair : att)
        err = nc_put_att_text( ncid, NC_GLOBAL,
            pair.first.data(), pair.second.size(), pair.second.data());

    int dim_ids[2], tvarID;
    std::map<std::string, int> id1d, id2d;
    err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid,
            {"time", "x"});
    //Create field IDs
    for( auto& record : equations::diagnostics1d_list)
    {
        std::string name = record.name;
        std::string long_name = record.long_name;
        id1d[name] = 0;
        err = nc_def_var( ncid, name.data(), NC_DOUBLE, 1, &dim_ids[0],
                &id1d.at(name));
        err = nc_put_att_text( ncid, id1d.at(name), "long_name", long_name.size(),
            long_name.data());
    }
    for( auto& record : equations::diagnostics_list)
    {
        std::string name = record.name;
        std::string long_name = record.long_name;
        id2d[name] = 0;
        err = nc_def_var( ncid, name.data(), NC_DOUBLE, 2, dim_ids,
                &id2d.at(name));
        err = nc_put_att_text( ncid, id2d.at(name), "long_name", long_name.size(),
            long_name.data());
    }
    err = nc_enddef(ncid);
    size_t start[2] = {0, 0};
    size_t count[2] = {1, grid.size()};
    dg::HVec result = y0[0][0];
    equations::Variables var = {ex, grid, y0, time, js, 0., 0};
    {   // update internal quantities
        Vector tmp ( y0);
        ex( time, y0, tmp);
    }
    for( auto& record : equations::diagnostics_list)
    {
        record.function( result, var);
        err = nc_put_vara_double( ncid, id2d.at(record.name), start, count,
                result.data());
    }
    for( auto& record : equations::diagnostics1d_list)
    {
        double result = record.function( var);
        nc_put_vara_double( ncid, id1d.at(record.name), start, count, &result);
    }
    err = nc_put_vara_double( ncid, tvarID, start, count, &time);

    //////////////////Main time loop ///////////////////////////////////
    dg::Timer t;
    double t_output = deltaT;
    while( (tend - time) > 1e-14 )
    {
        t.tic();
        while( time < t_output )
        {
            if( time+dt > t_output)
                dt = t_output-time;
            std::cout << "time " <<time <<" "<<dt<<" "<<t_output<<"\n";
            // Compute a step and error
            if( timestepper == "ARK")
                ark_adaptive.step( ex, im, time, y0, time, y0, dt,
                    dg::pid_control, dg::l2norm, rtol, atol);
            else if( timestepper == "ERK")
                erk_adaptive.step( ex, time, y0, time, y0, dt,
                    dg::pid_control, dg::l2norm, rtol, atol);
            if( erk_adaptive.failed() || ark_adaptive.failed())
            {
                var.nfailed++;
                continue;
            }
        }
        t_output += deltaT;
        t.toc();
        var.duration = t.diff();
        /////////////////////////////////output/////////////////////////
        std::cout << "Output time "<<time<<" of "<<tend<<"\n";
        start[0]++;
        for( auto& record : equations::diagnostics_list)
        {
            record.function( result, var);
            err = nc_put_vara_double( ncid, id2d.at(record.name), start, count,
                    result.data());
        }
        for( auto& record : equations::diagnostics1d_list)
        {
            double result = record.function( var);
            nc_put_vara_double( ncid, id1d.at(record.name), start, count, &result);
        }
        err = nc_put_vara_double( ncid, tvarID, start, count, &time);
    }
    err = nc_close(ncid);

    return 0;
}
