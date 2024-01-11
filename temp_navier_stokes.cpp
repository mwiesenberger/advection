#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "dg/file/file.h"
#include "common.h"

using Vector = std::array<dg::HVec,3>;

namespace equations
{

struct NavierStokesExplicit
{
    friend class NavierStokesImplicit;
    NavierStokesExplicit( dg::Grid1d g, dg::Grid1d vel_g, dg::file::WrappedJsonValue js) :
        m_velocity( g.size(), 0.), m_density(m_velocity), m_temperature(m_velocity), m_g(g), m_vel_g(vel_g)
    {
        dg::HVec temp( g.size()+4, 0.);
        m_yg.fill( temp);
        m_scheme = js["advection"].get("type", "upwind").asString();
        m_variant = js["advection"].get("variant", "original").asString();
        m_gamma = js["physical"].get("gamma", 1.0).asDouble();
        m_bc_n = dg::str2bc(js["bc"].get("density", "PER").asString());
        m_bc_u = dg::str2bc(js["bc"].get("velocity", "PER").asString());
        m_bc_t = dg::str2bc(js["bc"].get("temperature", "PER").asString());
        m_nu_n = js["physical"].get( "nu_n", 0.).asDouble();
        m_nu_t = js["physical"].get( "nu_t", 0.).asDouble();
        m_nu_u = js["physical"].get( "nu_u", 1.0).asDouble();
        std::cout << "Compute scheme "<<m_scheme<<" Variant " << m_variant<<"\n";
        m_init = js["init"].get( "type", "step").asString();
        if( "mms" == m_init)
        {
            m_n0 = js["init"].get("n_0", 1.0).asDouble();
            m_u0 = js["init"].get("u_0", 1.0).asDouble();
            m_t0 = js["init"].get("t_0", 1.0).asDouble();
            m_A = js["init"].get("A", 1.0).asDouble();
            m_B = js["init"].get("B", 1.0).asDouble();
            m_C = js["init"].get("C", 1.0).asDouble();
            m_k = js["init"].get("k", 1.0).asDouble();
            m_v = js["init"].get("v", 1.0).asDouble();
        }
    }

    const dg::HVec& density() const{return m_density;}
    const dg::HVec& velocity() const{return m_velocity;}
    const dg::HVec& temperature() const{return m_temperature;}
    void operator() ( double t, const Vector & y, Vector& yp)
    {
        m_called++;
        // y[0] -> density
        // y[1] -> velocity
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y[0], m_yg[0], m_bc_n);
        assign_ghost_cells( y[1], m_yg[1], m_bc_u);
        assign_ghost_cells( y[2], m_yg[2], m_bc_t);
        dg::Upwind upwind;
        SlopeLimiter<MinMod> limiter;
        // ghost cells are shifted by 2
        if ( m_scheme == "staggered")
        {
            dg::HVec qST(m_yg[1]), qTST(qST), q(qST), uh(qST), uST( q), dn(q), dnt(q), du(q);
            const dg::HVec & unST = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            const dg::HVec & nt = m_yg[2];
            for( unsigned k=0; k<Nx+3; k++)
            {
                // don't compute nSTinv here
                double nST = (nn[k] + nn[k+1])/2.;
                uST[k] = unST[k]/nST;
            }

            for( unsigned k=0; k<Nx+3; k++)
            {
                dn[k] = nn[k+1]-nn[k];
                dnt[k] = nt[k+1]-nt[k];
            }
            for( unsigned k=1; k<Nx+4; k++)
                du[k] = uST[k] - uST[k-1];
            for( unsigned k=1; k<Nx+3; k++)
            {
                qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                qTST[k] = upwind( uST[k], nt[k], nt[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                {
                    qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                    qTST[k] += limiter(uST[k], dnt[k-1], dnt[k], dnt[k+1], 0.5, 0.5);
                }
                qST[k]*= uST[k]; // k + 1/2
                qTST[k]*= uST[k]; // k + 1/2
            }
            for ( unsigned k=1; k<Nx+4; k++)
                q [k] = 0.5*(qST[k]+qST[k-1]);
            for ( unsigned k=1; k<Nx+3; k++)
            {
                uh[k] = upwind( q[k], uST[k-1], uST[k]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] += limiter(q[k], du[k-1], du[k], du[k+1], 0.5, 0.5);
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_density[i] = nn[k];
                m_velocity[i] = 0.5*(unST[k]+unST[k-1])/nn[k];
                m_temperature[i] = nt[k]/nn[k];
                yp[0][i] = -( qST[k] - qST[k-1])/hx
                           + m_nu_n*( nn[k+1]-2.*nn[k]+nn[k-1])/hx/hx;
                //
                yp[1][i] = -(uh[k+1]*q[k+1]-uh[k]*q[k])/hx;
                yp[1][i]+= m_nu_u*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx
                          +m_nu_n*( dn[k+1]*uST[k+1] - dn[k-1]*uST[k-1])/hx/hx/2.;

                double tt = nt[k]/ nn[k], tp = nt[k+1]/nn[k+1], tm = nt[k-1]/nn[k-1];
                yp[2][i] = -( qTST[k] - qTST[k-1])/hx
                           - (m_gamma-1.0)*nt[k]*du[k]/hx
                           + m_nu_n*( nt[k+1]-2.*nt[k]+nt[k-1])/hx/hx
                           + (m_gamma-1.)*m_nu_u*du[k]*du[k]/hx/hx
                           + (m_gamma-1.)*m_nu_t*( tp-2.*tt+tm)/hx/hx;
            }
            if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    yp[1][i] += -(nt[k+1] - nt[k])/hx;
                }
            }
        }
        else if ( m_scheme == "velocity-staggered")
        {
            // only velocity is in Lagrangian formulation, we keep pressure equation unchanged
            dg::HVec qST(m_yg[1]), qTST(qST), uu(qST), nST(qST), dn(uu), dnt(uu), du(uu), uh(uu);
            const dg::HVec & uST = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            const dg::HVec & nt = m_yg[2];

            for( unsigned k=0; k<Nx+3; k++)
            {
                nST[k] = 0.5*(nn[k+1]+nn[k]);
                dn[k] = nn[k+1]-nn[k];
                dnt[k] = nt[k+1]-nt[k];
            }
            for ( unsigned k=1; k<Nx+4; k++)
            {
                uu[k] = 0.5*(uST[k]*nST[k]+uST[k-1]*nST[k-1])/nn[k];
            }
            for( unsigned k=1; k<Nx+3; k++)
            {
                qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                qTST[k] = upwind( uST[k], nt[k], nt[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                {
                    qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                    qTST[k] += limiter(uST[k], dnt[k-1], dnt[k], dnt[k+1], 0.5, 0.5);
                }
                qST[k]*= uST[k]; // k + 1/2
                qTST[k]*= uST[k]; // k + 1/2
            }
            for( unsigned k=1; k<Nx+4; k++)
                du[k] = uST[k] - uST[k-1];
            for ( unsigned k=1; k<Nx+3; k++)
            {
                uh[k] = upwind( uu[k], uST[k-1], uST[k]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] += limiter(uu[k], du[k-1], du[k], du[k+1], 0.5, 0.5);
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_density[i] = nn[k];
                m_velocity[i] = 0.5*(uST[k]+uST[k-1]);
                m_temperature[i] = nt[k]/nn[k];
                yp[0][i] = -( qST[k] - qST[k-1])/hx
                           + m_nu_n*( nn[k+1]-2.*nn[k]+nn[k-1])/hx/hx;
                double nSTinv = 2./(nn[k] + nn[k+1]);
                yp[1][i] = -(uu[k+1]*uh[k+1]-uu[k]*uh[k])/2./hx;
                yp[1][i]+= m_nu_u*nSTinv*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx
                          +m_nu_n*nSTinv*(nn[k+1]-nn[k])/hx*(uST[k+1]-uST[k-1])/2./hx;

                double tt = nt[k]/ nn[k], tp = nt[k+1]/nn[k+1], tm = nt[k-1]/nn[k-1];
                yp[2][i] = -( qTST[k] - qTST[k-1])/hx
                           - (m_gamma-1.0)*nt[k]*du[k]/hx
                           + m_nu_n*( nt[k+1]-2.*nt[k]+nt[k-1])/hx/hx
                           + (m_gamma-1.)*m_nu_u*du[k]*du[k]/hx/hx
                           + (m_gamma-1.)*m_nu_t*( tp-2.*tt+tm)/hx/hx;
            }
            if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    double nSTinv = (1./nn[k] + 1./nn[k+1])/2.;
                    yp[1][i] += -nSTinv*(nt[k+1] - nt[k])/hx;
                }
            }
        }
        if( "mms" == m_init)
        {
            //Add sources
            dg::HVec tmpN( yp[0]), tmpNST( tmpN), tmpUST(tmpN), tmpT(tmpN);
            tmpN = dg::evaluate( [=](double x){
                return m_k*(m_A*m_k*m_nu_n*sin(m_k*(-(t*m_v) + x)) +
                        cos(m_k*(-(t*m_v) + x))*(m_B*m_n0 + m_A*m_u0 - m_A*m_v +
                        2*m_A*m_B*sin(m_k*(-(t*m_v) + x))));
                }, m_g);
            tmpNST = dg::evaluate( [=](double x){
                return m_k*(m_A*m_k*m_nu_n*sin(m_k*(-(t*m_v) + x)) +
                        cos(m_k*(-(t*m_v) + x))*(m_B*m_n0 + m_A*m_u0 - m_A*m_v +
                        2*m_A*m_B*sin(m_k*(-(t*m_v) + x))));
                }, m_vel_g);
            tmpUST = dg::evaluate( [=](double x){
                return (m_k*((m_A*(pow(m_B,2) + 4*m_t0) + 4*m_n0*(m_C + m_B*(m_u0 -
                m_v)))*cos(m_k*(-(t*m_v) + x)) - m_B*(4*m_A*m_k*m_nu_n*pow(cos(m_k*(-(t*m_v) +
                x)),2) + m_A*m_B*cos(3*m_k*(-(t*m_v) + x)) - 4*m_k*m_nu_u*sin(m_k*(-(t*m_v) +
                x))) + 2*(2*m_A*m_C + pow(m_B,2)*m_n0 + m_A*m_B*(m_u0 - m_v))*sin(2*m_k*(-(t*m_v)
                + x))))/(4.*(m_n0 + m_A*sin(m_k*(-(t*m_v) + x))));
                }, m_vel_g);
            tmpT = dg::evaluate( [=](double x){
                return (m_k*(-m_k*(2*m_A*m_C*m_nu_n + pow(m_B,2)*(-1 +
                m_gamma)*m_nu_u)*pow(cos(m_k*(-(t*m_v) + x)),2) + cos(m_k*(-(t*m_v) +
                x))*(m_n0 + m_A*sin(m_k*(-(t*m_v) + x)))*(m_B*(-1 + m_gamma)*m_t0 + m_C*(m_u0 -
                m_v) + m_B*m_C*m_gamma*sin(m_k*(-(t*m_v) + x))) + m_C*m_k*sin(m_k*(-(t*m_v) +
                x))*(m_n0*m_nu_n + (-1 + m_gamma)*m_nu_t + m_A*m_nu_n*sin(m_k*(-(t*m_v) +
                x)))))/(m_n0 + m_A*sin(m_k*(-(t*m_v) + x)));
                }, m_g);

            dg::HVec nn( yp[0]), nST(nn), uST(nn), tt(nn);
            nn = dg::evaluate( [=](double x){
                return m_n0 + m_A*sin( m_k*(-(t*m_v) + x));
                }, m_g);
            tt = dg::evaluate( [=](double x){
                return m_t0 + m_C*sin( m_k*(-(t*m_v) + x));
                }, m_g);
            nST = dg::evaluate( [=](double x){
                return m_n0 + m_A*sin( m_k*(-(t*m_v) + x));
                }, m_vel_g);
            uST = dg::evaluate( [=](double x){
                return m_u0 + m_B*sin( m_k*(-(t*m_v) + x));
                }, m_vel_g);
            dg::blas1::axpby( 1., tmpN, 1., yp[0]);
            for( unsigned i=0; i<Nx; i++)
                yp[2][i] += tmpN[i]*tt[i] + nn[i]*tmpT[i];
            if( m_scheme == "staggered")
            {
                for( unsigned i=0; i<Nx; i++)
                    yp[1][i] += tmpUST[i]*nST[i] + uST[i]*tmpNST[i];
            }
            else
                dg::blas1::axpby( 1., tmpUST, 1., yp[1]);

        }
    }
    unsigned called() const { return m_called;}
    private:
    std::string m_scheme, m_variant, m_init;
    Vector m_yg;
    dg::HVec m_velocity; // stores the velocity on non-staggered grid
    dg::HVec m_density;
    dg::HVec m_temperature;
    dg::Grid1d m_g, m_vel_g;
    dg::bc m_bc_n, m_bc_u, m_bc_t;
    double m_nu_u, m_nu_n, m_nu_t;
    unsigned m_called = 0;
    double m_gamma = 0;
    double m_n0 = 0, m_u0 = 0, m_t0 = 0, m_A = 0, m_B = 0, m_C = 0, m_k = 0, m_v = 0;
};

struct NavierStokesImplicit
{
    NavierStokesImplicit( NavierStokesExplicit& ex) : m_ex(ex){}

    void operator() ( double t, const Vector & y, Vector& yp)
    {
        dg::blas1::copy( 0., yp);
        unsigned Nx = m_ex.m_g.N();
        double hx = m_ex.m_g.h();
        assign_ghost_cells( y[0], m_ex.m_yg[0], m_ex.m_bc_n);
        assign_ghost_cells( y[2], m_ex.m_yg[2], m_ex.m_bc_t);
        // ghost cells are shifted by 2
        if( m_ex.m_scheme == "staggered")
        {
            if( m_ex.m_variant != "explicit" && m_ex.m_variant !=
                    "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    yp[1][i] = -(m_ex.m_yg[2][k+1] - m_ex.m_yg[2][k])/hx;
                }
            }
        }
        else if( m_ex.m_scheme == "velocity-staggered")
        {
            if( m_ex.m_variant != "explicit" && m_ex.m_variant !=
                    "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    yp[1][i] += -( m_ex.m_yg[2][k+1] - m_ex.m_yg[2][k])*(
                            1./m_ex.m_yg[0][k+1]+1./m_ex.m_yg[0][k])/2./hx;
                }
            }
        }
    }
    private:
    NavierStokesExplicit& m_ex;
};

struct NavierStokesImplicitSolver
{
    NavierStokesImplicitSolver( dg::Grid1d g, dg::file::WrappedJsonValue js, NavierStokesImplicit& im) :
        m_tmp( {dg::HVec(g.size(), 0.0), dg::HVec ( g.size(), 0.), dg::HVec (g.size(), 0.)}), m_im(im){}
    // solve (y - alpha I(t,y) = rhs
    void operator()( double alpha, double t, Vector& y, const Vector& rhs)
    {
        dg::blas1::copy( rhs[0], y[0]);// I_n = 0
        dg::blas1::copy( rhs[2], y[2]);// I_t = 0
        m_im( t, y, m_tmp); //ignores y[1],
        // writes 0 in m_tmp[0] and m_tmp[2] and updates m_tmp[1]
        dg::blas1::axpby( 1., rhs[1], +alpha, m_tmp[1], y[1]); // u = rhs_u + alpha I_u
    }
    private:
    Vector m_tmp;
    NavierStokesImplicit& m_im;

};

////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
struct Variables{
    NavierStokesExplicit& f;
    const dg::Grid1d& grid;
    const Vector& y0;
    const double& time;
    dg::file::WrappedJsonValue& js;
    double duration;
    const unsigned* nfailed;
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
    {"density", "Numerical density",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.density(), result);
        }
    },
    {"velocity", "Numerical velocity",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.velocity(), result);
        }
    },
    {"temperature", "Numerical temperature",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.temperature(), result);
        }
    },
    {"density_ana", "Analytical solution to the density",
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
    {"velocity_ana", "Analytical solution to the velocity",
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
    {"temperature_ana", "Analytical solution to the temperature",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if ( "mms" == init)
            {
                double t_0 = v.js["init"].get("t_0", 1.0).asDouble();
                double C = v.js["init"].get("C", 1.0).asDouble();
                double k = v.js["init"].get("k", 1.0).asDouble();
                double vel = v.js["init"].get("v", 1.0).asDouble();
                result = dg::evaluate( [=](double x){
                        return t_0 + C*sin( k *(x-vel*v.time));
                        }, v.grid);
            }
        }
    }
};
std::vector<Record1d> diagnostics1d_list = {
    {"failed", "Accumulated Number of failed steps",
        []( Variables& v ) {
            return *v.nfailed;
        }
    },
    {"duration", "Computation time for the latest output",
        []( Variables& v ) {
            return v.duration;
        }
    },
    {"nsteps", "Accumulated Number of calls to the RHS functor (including failed steps)",
        [](Variables& v) {
            return v.f.called();
        }
    }
};
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
////////////////////////////////DIAGNOSTICS END/////////////////////////////////
} //namespace equations

int main( int argc, char* argv[])
{
    ////Parameter initialisation ////////////////////////////////////////////
    dg::file::WrappedJsonValue js( dg::file::error::is_warning);
    if( argc == 1)
        dg::file::file2Json( "input/default.json", js.asJson(), dg::file::comments::are_discarded);
    else
        dg::file::file2Json( argv[1], js.asJson());
    std::cout << js.asJson() <<std::endl;

    /////////////////////////////////////////////////////////////////
    dg::Grid1d grid = equations::createGrid( js, dg::PER);
    dg::HVec w1d( dg::create::weights(grid));
    /////////////////////////////////////////////////////////////////
    std::string init = js["init"].get("type", "step").asString();
    std::string scheme = js["advection"].get("type", "staggered").asString();
    dg::Grid1d vel_grid = equations::createGrid( js, dg::PER);
    if ( "staggered" == scheme || "velocity-staggered" == scheme )
        vel_grid = equations::createStaggeredGrid( js, dg::PER);
    Vector y0 = {dg::evaluate( dg::zero, grid), dg::evaluate( dg::zero, grid), dg::evaluate( dg::zero, grid)};
    if( "riemann" == init)
    {
        // This is the classical Riemann problem
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        double n_l = js["init"].get("n_l", 1.0).asDouble();
        double n_r = js["init"].get("n_r", 1.0).asDouble();
        double u_l = js["init"].get("u_l", 0.).asDouble();
        double u_r = js["init"].get("u_r", 0.).asDouble();
        double t_l = js["init"].get("t_l", 0.).asDouble();
        double t_r = js["init"].get("t_r", 0.).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return x < x_a ? n_l : n_r;}, grid);
        y0[2] = dg::evaluate( [=](double x){ return x < x_a ? t_l : t_r;}, grid);
        dg::blas1::pointwiseDot( y0[0], y0[2], y0[2]);
        y0[1] = dg::evaluate( [=](double x){ return x < x_a ? u_l : u_r;}, vel_grid);
        if( scheme == "staggered")
        {
            //y0[1] = dg::evaluate( [=](double x){ return x < x_a ? n_l*u_l :
            //        n_r*u_r;}, vel_grid);
            //This one is what Herbin does:
            for( unsigned i=0; i<grid.N()-1; i++)
                y0[1][i] = (y0[0][i]+y0[0][i+1])*0.5*y0[1][i];
            y0[1][grid.N()-1] = y0[0][grid.N()-1]*y0[1][grid.N()-1];
        }

    }
    else if( "wave" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double u_0 = js["init"].get("u_0", 1.0).asDouble();
        double t_0 = js["init"].get("t_0", 1.0).asDouble();
        double A = js["init"].get("A", 1.0).asDouble();
        double B = js["init"].get("B", 1.0).asDouble();
        double C = js["init"].get("C", 1.0).asDouble();
        double k = js["init"].get("k", 1.0).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return n_0 + A*sin( k*x);}, grid);
        y0[2] = dg::evaluate( [=](double x){ return t_0 + C*sin( k*x);}, grid);
        dg::blas1::pointwiseDot( y0[0], y0[2], y0[2]);
        y0[1] = dg::evaluate( [=](double x){ return u_0 + B*sin( k*x);}, vel_grid);
        if( scheme == "staggered")
            y0[1] = dg::evaluate( [=](double x){ return
                    (u_0+B*sin(k*x))*(n_0+A*sin(k*x));}, vel_grid);
    }
    else if ( "mms" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double u_0 = js["init"].get("u_0", 1.0).asDouble();
        double t_0 = js["init"].get("t_0", 1.0).asDouble();
        double A = js["init"].get("A", 1.0).asDouble();
        double B = js["init"].get("B", 1.0).asDouble();
        double C = js["init"].get("C", 1.0).asDouble();
        double k = js["init"].get("k", 1.0).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return n_0 + A*sin( k*x);}, grid);
        y0[2] = dg::evaluate( [=](double x){ return t_0 + C*sin( k*x);}, grid);
        dg::blas1::pointwiseDot( y0[0], y0[2], y0[2]);
        y0[1] = dg::evaluate( [=](double x){ return u_0 + B*sin( k*x);}, vel_grid);
        if( scheme == "staggered")
            y0[1] = dg::evaluate( [=](double x){ return
                    (u_0+B*sin(k*x))*(n_0+A*sin(k*x));}, vel_grid);
    }

    equations::NavierStokesExplicit ex( grid, vel_grid, js);
    equations::NavierStokesImplicit im( ex);
    equations::NavierStokesImplicitSolver solver( grid, js, im);

    std::string timestepper = js["timestepper"].get( "type", "ERK").asString();
    std::string tableau= js["timestepper"].get("tableau",
            "ARK-4-2-3").asString();
    double tend = js["output"].get( "tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = tend/(double)maxout;

    dg::Adaptive<dg::ARKStep<Vector>> ark_adaptive;
    dg::Adaptive<dg::ERKStep<Vector >> erk_adaptive;
    auto odeint = std::unique_ptr<dg::aTimeloop<Vector>>();
    double time = 0.;
    const unsigned failed = 0;
    const unsigned* nfailed = &failed;
    if( timestepper == "ARK")
    {
        ark_adaptive = { tableau, y0};
        double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
        double atol= js["timestepper"].get("atol", 1e-6).asDouble();
        odeint = std::make_unique<dg::AdaptiveTimeloop<Vector>>(
                    ark_adaptive, std::tie( ex, im, solver), dg::pid_control,
                    dg::l2norm, rtol, atol);
        nfailed = &ark_adaptive.nfailed();
    }
    else if( timestepper == "ERK")
    {
        erk_adaptive = { tableau, y0};
        double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
        double atol= js["timestepper"].get("atol", 1e-6).asDouble();
        odeint = std::make_unique<dg::AdaptiveTimeloop<Vector>>(
                    erk_adaptive, ex, dg::pid_control,
                    dg::l2norm, rtol, atol);
        nfailed = &erk_adaptive.nfailed();
    }
    else if( timestepper == "RK")
    {
        double dt = js["timestepper"].get("dt", 0.1).asDouble();
        odeint = std::make_unique<dg::SinglestepTimeloop<Vector>>(
                    dg::RungeKutta<Vector>(tableau, y0), ex, dt);
    }

    // Set up netcdf
    std::string inputfile = js.asJson().toStyledString(); //save input without comments, which is important if netcdf file is later read by another parser
    std::string outputfile;
    if( argc == 1 || argc == 2)
        outputfile = "navier-stokes.nc";
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
    att["title"] = "Output file of advection/navier_stokes.cpp";
    att["Conventions"] = "CF-1.7";
    ///Get local time and begin file history
    auto ttt = std::time(nullptr);
    auto tm = *std::localtime(&ttt);

    std::ostringstream oss;
    ///time string  + program-name + args
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    for( int i=0; i<argc; i++) oss << " "<<argv[i];
    att["history"] = oss.str();
    att["comment"] = "Find more info in advection/temp-navier-stokes.ipynb";
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
    dg::HVec result = y0[0];
    equations::Variables var = {ex, grid, y0, time, js, 0., nfailed};
    {   // update internal velocity
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
    for( unsigned u=1; u<=maxout; u++)
    {
        t.tic();
        odeint->integrate( time, y0, time+deltaT, y0, dg::to::exact);
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
