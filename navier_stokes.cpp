#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "dg/file/file.h"

namespace equations
{

dg::Grid1d createGrid( Json::Value grid)
{
    unsigned Nx = grid.get( "Nx", 48).asUInt();
    double x0 = grid["x"].get( 0u, 0.).asDouble();
    double x1 = grid["x"].get( 1u, 1.).asDouble();
    return dg::Grid1d( x0, x1, 1, Nx, dg::PER);
}
dg::Grid1d createStaggeredGrid( Json::Value grid)
{
    dg::Grid1d g1d = createGrid(grid);
    unsigned Nx = grid.get( "Nx", 48).asUInt();
    double x0 = g1d.x0() + g1d.h()/2.;
    double x1 = g1d.x1() + g1d.h()/2.;
    return dg::Grid1d( x0, x1, 1, Nx, dg::PER);
}

void assign_ghost_cells( const dg::HVec& y, dg::HVec& yg, dg::bc bcx)
{
    unsigned Nx = y.size();
    for( unsigned i=0; i<Nx; i++)
        yg[i+2] = y[i];
    //assign ghost values
    if( bcx == dg::PER)
    {
        yg[Nx+2] = y[0]; // right boundary
        yg[Nx+3] = y[1];
        yg[1] = y[Nx-1]; // left boundary
        yg[0] = y[Nx-2];
    }
    else
    {
        if ( bcx == dg::NEU || bcx == dg::NEU_DIR)
        {
            yg[Nx+2] = y[Nx-1];
            yg[Nx+3] = y[Nx-2];
        }
        if( bcx == dg::NEU || bcx == dg::DIR_NEU)
        {
            yg[1] = y[0];
            yg[0] = y[1];
        }
        if ( bcx == dg::DIR || bcx == dg::DIR_NEU)
        {
            yg[Nx+2] = -y[Nx-1];
            yg[Nx+3] = -y[Nx-2];
        }
        if( bcx == dg::DIR || bcx == dg::NEU_DIR)
        {
            yg[1] = -y[0];
            yg[0] = -y[1];
        }
    }
}

double minmod( double x, double y)
{
    if( x >= 0 && y >= 0)
        return std::min( x,y);
    else if( x <=0 && y <= 0)
        return std::max(x,y);
    return 0.;
}

struct NavierStokesExplicit
{
    NavierStokesExplicit( dg::Grid1d g, dg::Grid1d vel_g, Json::Value js) :
        m_velocity( g.size(), 0.), m_g(g), m_vel_g(vel_g)
    {
        dg::HVec temp( g.size()+4, 0.);
        m_yg.fill( temp);
        m_scheme = js["advection"].get("type", "upwind").asString();
        m_variant = js["advection"].get("variant", "original").asString();
        m_alpha = js["physical"].get("alpha", 1.0).asDouble();
        m_gamma = js["physical"].get("gamma", 1.0).asDouble();
        m_nu_u = js["physical"].get("nu_u", 1.0).asDouble();
        m_bc_n = dg::str2bc(js["bc"].get("density", "PER").asString());
        m_bc_u = dg::str2bc(js["bc"].get("velocity", "PER").asString());
        m_nu_n = js["physical"].get( "nu_n", 0.).asDouble();
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
    }

    dg::HVec velocity() const{return m_velocity;}
    void operator() ( double t, const std::array<dg::HVec,2> & y, std::array<dg::HVec, 2>& yp)
    {
        // y[0] -> density
        // y[1] -> velocity
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y[0], m_yg[0], m_bc_n);
        assign_ghost_cells( y[1], m_yg[1], m_bc_u);
        // ghost cells are shifted by 2
        if ( m_scheme == "upwind")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                if( m_yg[1][k] > 0.)
                    yp[0][i] =  -m_yg[1][k]*( m_yg[0][k] - m_yg[0][k-1])/ hx;
                else
                    yp[0][i] =  -m_yg[1][k]*( m_yg[0][k+1] - m_yg[0][k])/ hx;
                yp[0][i] += -m_yg[0][k]*( m_yg[1][k+1]-m_yg[1][k-1])/2./hx;
                yp[1][i] =
                    -m_yg[1][k]*(m_yg[1][k+1]-m_yg[1][k-1])/2./hx -
                    m_alpha*(std::pow(m_yg[0][k+1],
                                m_gamma)-std::pow(m_yg[0][k-1],
                                    m_gamma))/2./hx/m_yg[0][k]+
                    m_nu_u/m_yg[0][k]*(m_yg[1][k+1] -
                            2.*m_yg[1][k]+m_yg[1][k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "upwind2")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                if( m_yg[1][k] > 0.)
                    yp[0][i] =  -m_yg[1][k]*( 3*m_yg[0][k] - 4*m_yg[0][k-1] +
                            m_yg[0][k-2])/2./hx;
                else
                    yp[0][i] =  -m_yg[1][k]*( -m_yg[0][k+2] + 4*m_yg[0][k+1] -
                            3*m_yg[0][k])/2./hx;
                yp[0][i] += -m_yg[0][k]*( m_yg[1][k+1]-m_yg[1][k-1])/2./hx;
                yp[1][i] =
                    -m_yg[1][k]*(m_yg[1][k+1]-m_yg[1][k-1])/2./hx -
                    m_alpha*(std::pow(m_yg[0][k+1],
                                m_gamma)-std::pow(m_yg[0][k-1],
                                    m_gamma))/2./hx/m_yg[0][k]+
                    m_nu_u/m_yg[0][k]*(m_yg[1][k+1] -
                            2.*m_yg[1][k]+m_yg[1][k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "centered")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                yp[0][i] =  -m_yg[1][k]*( m_yg[0][k+1]-m_yg[0][k-1])/2./hx
                            -m_yg[0][k]*( m_yg[1][k+1]-m_yg[1][k-1])/2./hx;
                yp[1][i] =
                    -m_yg[1][k]*(m_yg[1][k+1]-m_yg[1][k-1])/2./hx -
                    m_alpha*(std::pow(m_yg[0][k+1],
                                m_gamma)-std::pow(m_yg[0][k-1],
                                    m_gamma))/2./hx/m_yg[0][k]+
                    m_nu_u/m_yg[0][k]*(m_yg[1][k+1] -
                            2.*m_yg[1][k]+m_yg[1][k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "staggered")
        {
            dg::HVec qST(m_yg[1]), q(qST), uh(qST), uST( q);
            const dg::HVec & unST = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            for( unsigned k=0; k<Nx+3; k++)
            {
                double nST = (nn[k] + nn[k+1])/2.;
                uST[k] = unST[k]/nST;
            }
            for( unsigned k=1; k<Nx+2; k++)
            {
                qST[k] = uST[k] >= 0 ? nn[k] : nn[k+1];
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    qST[k] = uST[k] >=0 ? nn[k] + 1./2.*minmod( nn[k+1]-nn[k],
                            nn[k]-nn[k-1]) : nn[k+1] - 1./2.*minmod(
                                nn[k+2]-nn[k+1], nn[k+1]-nn[k]);
                qST[k]*= uST[k]; // k + 1/2
            }
            for ( unsigned k=1; k<Nx+3; k++)
                q [k] = 0.5*(qST[k]+qST[k-1]);
            for ( unsigned k=2; k<Nx+3; k++)
            {
                uh[k] = q[k] >= 0 ? uST[k-1] : uST[k];
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] = q[k] >=0 ? uST[k-1] + 1./2.*minmod(
                            uST[k]-uST[k-1], uST[k-1]-uST[k-2]) : uST[k] -
                        1./2.*minmod( uST[k+1]-uST[k], uST[k]-uST[k-1]);
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_velocity[i] = 0.5*(unST[k]+unST[k-1])/nn[k];
                yp[0][i] = -( qST[k] - qST[k-1])/hx;
                //
                yp[1][i] = -(uh[k+1]*q[k+1]-uh[k]*q[k])/hx;
                yp[1][i]+= m_nu_u*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx;
            }
            if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    //yp[1][i] += -m_alpha*(  std::pow(nn[k+1], m_gamma) -
                    //                        std::pow(nn[k], m_gamma))/hx;
                    yp[1][i] += -m_alpha*m_gamma*std::pow( nn[k], m_gamma-1)*(
                            nn[k+1] -nn[k])/hx;
                    // eventual sources (due to floor level for example) should also
                    // go here
                }
            }
        }
        else if ( m_scheme == "velocity-staggered")
        {
            dg::HVec qST(m_yg[1]), q(qST), uh(qST);
            const dg::HVec & uST = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            for( unsigned k=1; k<Nx+2; k++)
            {
                qST[k] = uST[k] >= 0 ? nn[k] : nn[k+1];
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    qST[k] = uST[k] >=0 ? nn[k] + 1./2.*minmod( nn[k+1]-nn[k],
                            nn[k]-nn[k-1]) : nn[k+1] - 1./2.*minmod(
                                nn[k+2]-nn[k+1], nn[k+1]-nn[k]);
                qST[k]*= uST[k]; // k + 1/2
            }
            for ( unsigned k=1; k<Nx+3; k++)
                q [k] = 0.5*(uST[k]+uST[k-1]);
            for ( unsigned k=2; k<Nx+3; k++)
            {
                uh[k] = q[k] >= 0 ? uST[k-1] : uST[k];
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] = q[k] >=0 ? uST[k-1] + 1./2.*minmod(
                            uST[k]-uST[k-1], uST[k-1]-uST[k-2]) : uST[k] -
                        1./2.*minmod( uST[k+1]-uST[k], uST[k]-uST[k-1]);
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_velocity[i] = 0.5*(uST[k]+uST[k-1]);
                yp[0][i] = -( qST[k] - qST[k-1])/hx;
                //
                yp[1][i] = -(uh[k+1]*q[k+1]-uh[k]*q[k])/2./hx;
                double nST = (nn[k] + nn[k+1])/2.;
                yp[1][i]+= m_nu_u/nST*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx;
            }
            if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    if( m_gamma == 1)
                        yp[1][i] += -m_alpha*(  log(nn[k+1]) - log(nn[k]))/hx;
                    else
                        yp[1][i] += -m_gamma/(m_gamma-1)*m_alpha*(
                                pow(nn[k+1], m_gamma-1) - pow(nn[k],
                                    m_gamma-1))/hx;
                    // eventual sources (due to floor level for example) should also
                    // go here
                }
            }
        }
        for( unsigned i=0; i<Nx; i++)
        {
            unsigned k=i+2;
            yp[0][i] +=  m_nu_n*( m_yg[0][k+1]-2.*m_yg[0][k]+m_yg[0][k-1])/hx/hx;
        }

        if( "mms" == m_init)
        {
            //Add sources
            dg::HVec tmpN( yp[0]), tmpNST( tmpN), tmpUST(tmpN);
            if( m_gamma == 1)
            {
                tmpN = dg::evaluate( [=](double x){
                    return m_k*cos(m_k*(-(t*m_v) + x))*(m_B*m_n0 + m_A*(m_u0 - m_v)
                            + 2*m_A*m_B*sin(m_k*(-(t*m_v) + x)));
                    }, m_g);
                tmpNST = dg::evaluate( [=](double x){
                    return m_k*cos(m_k*(-(t*m_v) + x))*(m_B*m_n0 + m_A*(m_u0 - m_v)
                            + 2*m_A*m_B*sin(m_k*(-(t*m_v) + x)));
                    }, m_vel_g);
                tmpUST = dg::evaluate( [=](double x){
                    return m_k*((m_B*m_nu_u*m_k*sin(m_k*(-(t*m_v) + x)))/(m_n0
                                + m_A*sin(m_k*(-(t*m_v) + x))) +
                            cos(m_k*(-(t*m_v) + x))*(m_B*m_u0 - m_B*m_v +
                                m_B*m_B*sin(m_k*(-(t*m_v) + x)) +
                                (m_A*m_alpha)/(m_n0 + m_A*sin(m_k*(-(t*m_v) +
                                            x)))));
                    }, m_vel_g);
                dg::blas1::axpby( 1., tmpN, 1., yp[0]);
                if( m_scheme == "staggered")
                {
                    dg::HVec nST(m_yg[0]), uST( nST);
                    const dg::HVec & unST = m_yg[1];
                    const dg::HVec & nn = m_yg[0];
                    for( unsigned k=0; k<Nx+3; k++)
                    {
                        nST[k] = (nn[k] + nn[k+1])/2.;
                        uST[k] = unST[k]/nST[k];
                    }
                    for( unsigned i=0; i<Nx; i++)
                    {
                        unsigned k=i+2;
                        yp[1][i] += tmpUST[i]*nST[k] + uST[k]*tmpNST[i];
                    }
                }
                else
                    dg::blas1::axpby( 1., tmpUST, 1., yp[1]);

            }
        }
    }
    private:
    std::string m_scheme, m_variant, m_init;
    std::array<dg::HVec,2> m_yg;
    dg::HVec m_velocity; // stores the velocity on non-staggered grid
    dg::Grid1d m_g, m_vel_g;
    dg::bc m_bc_n, m_bc_u;
    double m_alpha, m_gamma, m_nu_u, m_nu_n;
    double m_n0 = 0, m_u0 = 0, m_A = 0, m_B = 0, m_k = 0, m_v = 0;
};

struct NavierStokesImplicit
{
    NavierStokesImplicit( dg::Grid1d g, Json::Value js) :
        m_g(g)
    {
        dg::HVec temp( g.size()+4, 0.);
        m_yg.fill( temp);
        m_scheme = js["advection"].get("type", "upwind").asString();
        m_variant = js["advection"].get("variant", "original").asString();
        m_alpha = js["physical"].get("alpha", 1.0).asDouble();
        m_gamma = js["physical"].get("gamma", 1.0).asDouble();
        m_bc_n = dg::str2bc(js["bc"].get("density", "PER").asString());
        m_bc_u = dg::str2bc(js["bc"].get("velocity", "PER").asString());
    }

    void operator() ( double t, const std::array<dg::HVec,2> & y, std::array<dg::HVec,2>& yp)
    {
        dg::blas1::copy( 0., yp);
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y[0], m_yg[0], m_bc_n);
        assign_ghost_cells( y[1], m_yg[1], m_bc_u);
        // ghost cells are shifted by 2
        if( m_scheme == "staggered")
        {
            if( m_variant != "explicit" && m_variant != "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    yp[1][i] = -m_alpha*(std::pow(m_yg[0][k+1], m_gamma) -
                                         std::pow(m_yg[0][k], m_gamma))/hx;
                    // eventual sources (due to floor level for example) should also
                    // go here
                }
            }
        }
        else if( m_scheme == "velocity-staggered")
        {
            if( m_variant != "explicit" && m_variant != "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    if( m_gamma == 1)
                        yp[1][i] += -m_alpha*(  log(m_yg[0][k+1]) -
                                log(m_yg[0][k]))/hx;
                    else
                        yp[1][i] += -m_gamma/(m_gamma-1)*m_alpha*(
                                pow(m_yg[0][k+1], m_gamma-1) - pow(m_yg[0][k],
                                    m_gamma-1))/hx;
                    // eventual sources (due to floor level for example) should also
                    // go here
                }
            }
        }
    }
    private:
    double m_alpha, m_gamma;
    std::array<dg::HVec,2> m_yg;
    dg::Grid1d m_g;
    dg::bc m_bc_n, m_bc_u;
    std::string m_scheme, m_variant;
};

struct NavierStokesImplicitSolver
{
    NavierStokesImplicitSolver( dg::Grid1d g, Json::Value js) :
        m_tmp( {dg::HVec(g.size(), 0.0), dg::HVec ( g.size(), 0.)}){}
    std::array<dg::HVec,2> copyable() const{
        return m_tmp;
    }
    // solve (y + alpha I(t,y) = rhs
    void solve( double alpha, NavierStokesImplicit& im, double t, std::array<dg::HVec,2>& y, const std::array<dg::HVec,2>& rhs)
    {
        im( t, y, m_tmp);
        dg::blas1::copy( rhs[0], y[0]);// I_n = 0
        dg::blas1::axpby( 1., rhs[1], -alpha, m_tmp[1], y[1]); // u = rhs_u - alpha I_u
    }
    private:
    std::array<dg::HVec,2> m_tmp;

};

////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
////////////////////////////////DIAGNOSTICS/////////////////////////////////////
struct Variables{
    NavierStokesExplicit& f;
    const dg::Grid1d& grid;
    const std::array<dg::HVec,2>& y0;
    const double& time;
    Json::Value& js;
    double duration;
    unsigned nfailed;
    unsigned nsteps;
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
             dg::blas1::copy(v.y0[0], result);
        }
    },
    {"velocity", "Numerical velocity",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.f.velocity(), result);
        }
    },
    {"density_ana", "Analytical solution to the density",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if( "step" == init)
            {
                double alpha = v.js["physical"].get("alpha",2).asDouble();
                double gamma = v.js["physical"].get("gamma",2).asDouble();
                if( gamma == 2)
                {
                    double x_a = v.js["init"].get("x_a", 0.1).asDouble();
                    double h_l = v.js["init"].get("n_l", 1.0).asDouble();
                    double h_r = v.js["init"].get("n_r", 1.0).asDouble();
                    const double g = alpha*2;
                    double x_A = x_a - v.time*sqrt( 2*alpha*h_l);
                    double cmm = sqrt(g*h_l), cmp = sqrt(g*h_r);
                    auto lambda = [=](double cm) {
                            return -8.*g*h_r*cm*cm*(sqrt(g*h_l) -
                                    cm)*(sqrt(g*h_l)-cm) +
                            (cm*cm-g*h_r)*(cm*cm-g*h_r)*(cm*cm+g*h_r);};
                    dg::bisection1d( lambda, cmm, cmp, 1e-6);
                    double cm = (cmm+cmp)/2.;
                    double x_B = x_a + v.time*(2.*sqrt( g*h_l)-3*cm);
                    double x_C = x_a + v.time*2*cm*cm*(sqrt(g*h_l) -
                                cm)/(cm*cm-g*h_r);
                    result = dg::evaluate( [=](double x){
                            if( x <= x_A)
                                return h_l;
                            if( x <= x_B)
                                return 4./9./g*std::pow(-(x-x_a)/2./v.time + sqrt(g*h_l), 2);
                            if( x <= x_C)
                                return cm*cm/g;
                            return h_r;}, v.grid);
                }
            }
            else if ( "mms" == init)
            {
                double gamma = v.js["physical"].get("gamma",2).asDouble();
                if( gamma == 1 )
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
        }
    },
    {"velocity_ana", "Analytical solution to the velocity",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            result = dg::evaluate( dg::zero, v.grid);
            if( "step" == init)
            {
                double alpha = v.js["physical"].get("alpha",2).asDouble();
                double gamma = v.js["physical"].get("gamma",2).asDouble();
                if( gamma == 2)
                {
                    double x_a = v.js["init"].get("x_a", 0.1).asDouble();
                    double h_l = v.js["init"].get("n_l", 1.0).asDouble();
                    double h_r = v.js["init"].get("n_r", 1.0).asDouble();
                    const double g = alpha*2;
                    double x_A = x_a - v.time*sqrt( 2*alpha*h_l);
                    double cmm = sqrt(g*h_l), cmp = sqrt(g*h_r);
                    auto lambda = [=](double cm) {
                            return -8.*g*h_r*cm*cm*(sqrt(g*h_l) -
                                    cm)*(sqrt(g*h_l)-cm) +
                            (cm*cm-g*h_r)*(cm*cm-g*h_r)*(cm*cm+g*h_r);};
                    dg::bisection1d( lambda, cmm, cmp, 1e-6);
                    double cm = (cmm+cmp)/2.;
                    double x_B = x_a + v.time*(2.*sqrt( g*h_l)-3*cm);
                    double x_C = x_a + v.time*2*cm*cm*(sqrt(g*h_l) -
                                cm)/(cm*cm-g*h_r);
                    result = dg::evaluate( [=](double x){
                            if( x <= x_A)
                                return 0.;
                            if( x <= x_B)
                                return 2./3.*((x-x_a)/v.time + sqrt(g*h_l));
                            if( x <= x_C)
                                return 2.*(sqrt(g*h_l) -cm);
                            return 0.;}, v.grid);
                }
            }
            else if ( "mms" == init)
            {
                double gamma = v.js["physical"].get("gamma",2).asDouble();
                if( gamma == 1 )
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
        }
    }
};
std::vector<Record1d> diagnostics1d_list = {
    {"failed", "Number of failed steps",
        []( Variables& v ) {
            return v.nfailed;
        }
    },
    {"duration", "Computation time for the latest output (i.e. nsteps timesteps)",
        []( Variables& v ) {
            return v.duration;
        }
    },
    {"nsteps", "Number of calls to the timestepper (including failed steps)",
        [](Variables& v) {
            return v.nsteps;
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
    Json::Value js;
    if( argc == 1)
        dg::file::file2Json( "input/default.json", js, dg::file::comments::are_discarded);
    else
        dg::file::file2Json( argv[1], js);
    std::cout << js <<std::endl;

    /////////////////////////////////////////////////////////////////
    dg::Grid1d grid = equations::createGrid( js["grid"]);
    dg::HVec w1d( dg::create::weights(grid));
    /////////////////////////////////////////////////////////////////
    std::string init = js["init"].get("type", "step").asString();
    std::string scheme = js["advection"].get("type", "staggered").asString();
    dg::Grid1d vel_grid = equations::createGrid( js["grid"]);
    if ( "staggered" == scheme || "velocity-staggered" == scheme)
        vel_grid = equations::createStaggeredGrid( js["grid"]);
    std::array<dg::HVec,2> y0 = {dg::evaluate( dg::zero, grid), dg::evaluate( dg::zero, grid)};
    if( "step" == init)
    {
        // This is the classical Riemann problem (Dam break for shallow water)
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        double n_l = js["init"].get("n_l", 1.0).asDouble();
        double n_r = js["init"].get("n_r", 1.0).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return x < x_a ? n_l : n_r;}, grid);
    }
    else if( "riemann" == init)
    {
        // This is the classical Riemann problem (Dam break for shallow water)
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        double n_l = js["init"].get("n_l", 1.0).asDouble();
        double n_r = js["init"].get("n_r", 1.0).asDouble();
        double u_l = js["init"].get("u_l", 0.).asDouble();
        double u_r = js["init"].get("u_r", 0.).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return x < x_a ? n_l : n_r;}, grid);
        y0[1] = dg::evaluate( [=](double x){ return x < x_a ? u_l : u_r;}, vel_grid);
        if( scheme == "staggered")
            y0[1] = dg::evaluate( [=](double x){ return x < x_a ? n_l*u_l :
                    n_r*u_r;}, vel_grid);
    }
    else if( "wave" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double amp = js["init"].get("amp", 1.0).asDouble();
        double k   = js["init"].get("k", 1.0).asDouble();
        double x_0 = js["init"].get("x_0", 1.0).asDouble();
        y0[0] = dg::evaluate( [=]( double x){ return n_0 + amp*sin( k*(x-x_0));},
                grid);
    }
    else if ( "mms" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double u_0 = js["init"].get("u_0", 1.0).asDouble();
        double A = js["init"].get("A", 1.0).asDouble();
        double B = js["init"].get("B", 1.0).asDouble();
        double k = js["init"].get("k", 1.0).asDouble();
        y0[0] = dg::evaluate( [=](double x){ return n_0 + A*sin( k*x);}, grid);
        y0[1] = dg::evaluate( [=](double x){ return u_0 + B*sin( k*x);}, vel_grid);
        if( scheme == "staggered")
            y0[1] = dg::evaluate( [=](double x){ return
                    (u_0+B*sin(k*x))*(n_0+A*sin(k*x));}, vel_grid);
    }

    std::string timestepper = js["timestepper"].get( "type", "ERK").asString();
    std::string tableau= js["timestepper"].get("tableau",
            "ARK-4-2-3").asString();
    double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
    double atol= js["timestepper"].get("atol", 1e-6).asDouble();
    double tend = js["output"].get( "tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = tend/(double)maxout;
    dg::ARKStep<std::array<dg::HVec,2> , equations::NavierStokesImplicitSolver>
        ark_fixed( dg::EXPLICIT_EULER_1_1, dg::IMPLICIT_EULER_1_1, grid, js);
    dg::RungeKutta<std::array<dg::HVec,2> > erk_fixed( "Euler", y0);

    dg::Adaptive<dg::ARKStep<std::array<dg::HVec,2> ,
        equations::NavierStokesImplicitSolver>> ark_adaptive( "ARK-4-2-3", grid, js);
    dg::Adaptive<dg::ERKStep<std::array<dg::HVec,2> >> erk_adaptive( "Bogacki-Shampine-4-2-3", y0);
    if( timestepper == "ARK")
        ark_adaptive = dg::Adaptive<dg::ARKStep<std::array<dg::HVec,2> ,
            equations::NavierStokesImplicitSolver>>( tableau, grid, js);
    else if( timestepper == "ERK")
    {
        erk_fixed = dg::RungeKutta<std::array<dg::HVec,2> >( tableau, y0);
        erk_adaptive = dg::Adaptive<dg::ERKStep<std::array<dg::HVec,2> >>( tableau, y0);
    }
    double dt = 1e-6, time = 0.;
    equations::NavierStokesExplicit exp( grid, vel_grid, js);
    equations::NavierStokesImplicit imp( grid,js);

    // Set up netcdf
    std::string inputfile = js.toStyledString(); //save input without comments, which is important if netcdf file is later read by another parser
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
    dg::HVec result = y0[0];
    equations::Variables var = {exp, grid, y0, time, js, 0., 0, 0};
    {   // update internal velocity
        std::array<dg::HVec , 2> tmp ( y0);
        exp( time, y0, tmp);
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
                ark_adaptive.step( exp, imp, time, y0, time, y0, dt,
                    dg::pid_control, dg::l2norm, rtol, atol);
            else if( timestepper == "ERK")
                erk_adaptive.step( exp, time, y0, time, y0, dt,
                    dg::pid_control, dg::l2norm, rtol, atol);
            var.nsteps++;
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

    //dg::Timer t;
    //unsigned NT = 1000;
    //dt = grid.h()/ 20.;//tend / (double)NT/(double)maxout;
    //std::cout << "dt "<<dt<<"\n";
    //NT = ( tend/dt/(maxout));
    //std::array<dg::HVec,2> delta(y0);
    //for( unsigned i=0; i<maxout; i++)
    //{
    //    t.tic();
    //    if( timestepper == "ARK")
    //        for( unsigned k=0; k<NT; k++)
    //            ark_fixed.step( exp, imp, time, y0, time, y0, dt, delta);
    //    else if( timestepper == "ERK")
    //        for( unsigned k=0; k<NT; k++)
    //            erk_fixed.step( exp, time, y0, time, y0, dt);
    //    var.nsteps += NT;
    //    t.toc();
    //    var.duration = t.diff();
    //    /////////////////////////////////output/////////////////////////
    //    std::cout << "Output time "<<time<<" of "<<tend<<" "<<var.nsteps<<"\n";
    //    start[0]++;
    //    for( auto& record : equations::diagnostics_list)
    //    {
    //        record.function( result, var);
    //        err = nc_put_vara_double( ncid, id2d.at(record.name), start, count,
    //                result.data());
    //    }
    //    for( auto& record : equations::diagnostics1d_list)
    //    {
    //        double result = record.function( var);
    //        nc_put_vara_double( ncid, id1d.at(record.name), start, count, &result);
    //    }
    //    err = nc_put_vara_double( ncid, tvarID, start, count, &time);
    //}
    err = nc_close(ncid);

    return 0;
}
