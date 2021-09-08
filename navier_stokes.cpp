#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "dg/file/file.h"
#include "common.h"

namespace equations
{

struct NavierStokesExplicit
{
    friend class NavierStokesImplicit;
    NavierStokesExplicit( dg::Grid1d g, dg::Grid1d vel_g, Json::Value js) :
        m_velocity( g.size(), 0.), m_density(m_velocity), m_g(g), m_vel_g(vel_g)
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

    const dg::HVec& density() const{return m_density;}
    const dg::HVec& velocity() const{return m_velocity;}
    void operator() ( double t, const std::array<dg::HVec,2> & y, std::array<dg::HVec, 2>& yp)
    {
        m_called++;
        // y[0] -> density
        // y[1] -> velocity
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y[0], m_yg[0], m_bc_n);
        assign_ghost_cells( y[1], m_yg[1], m_bc_u);
        dg::Upwind upwind;
        SlopeLimiter<MinMod> limiter;
        // ghost cells are shifted by 2
        if ( m_scheme == "upwind")
        {
            dg::HVec q(m_yg[0]);
            const dg::HVec & nn = m_yg[0];
            const dg::HVec & uu = m_yg[1];
            for ( unsigned k=1; k<Nx+3; k++)
            {
                q[k] = 0.5*(uu[k]+uu[k-1]) >= 0 ? uu[k-1]*uu[k-1] : uu[k]*uu[k];
                q[k] /= 2.;
            }
            for( unsigned i=0; i<Nx; i++)
            {
                m_density[i] = y[0][i];
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                if( uu[k] > 0.)
                    yp[0][i] =  -uu[k]*( nn[k] - nn[k-1])/ hx;
                else
                    yp[0][i] =  -uu[k]*( nn[k+1] - nn[k])/ hx;
                yp[0][i] += -nn[k]*( uu[k+1]-uu[k-1])/2./hx;
                yp[1][i] =  -(q[k+1]-q[k])/hx -
                    m_alpha*(std::pow(nn[k+1], m_gamma)-std::pow(nn[k-1],
                                m_gamma))/2./hx/nn[k]+
                    m_nu_u/nn[k]*(uu[k+1]-2.*uu[k]+uu[k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "upwind2")
        {
            const dg::HVec & nn = m_yg[0];
            const dg::HVec & uu = m_yg[1];
            for( unsigned i=0; i<Nx; i++)
            {
                m_density[i] = y[0][i];
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                if( m_yg[1][k] > 0.)
                    yp[0][i] =  -uu[k]*( 3*nn[k] - 4*nn[k-1] + nn[k-2])/2./hx;
                else
                    yp[0][i] =  -uu[k]*( -nn[k+2] + 4*nn[k+1] - 3*nn[k])/2./hx;
                yp[0][i] += -nn[k]*( uu[k+1]-uu[k-1])/2./hx;
                yp[1][i] =
                    -uu[k]*(uu[k+1]-uu[k-1])/2./hx -
                    m_alpha*(std::pow(nn[k+1], m_gamma)-std::pow(nn[k-1],
                                m_gamma))/2./hx/nn[k]+
                    m_nu_u/nn[k]*(uu[k+1] - 2.*uu[k]+uu[k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "centered")
        {
            const dg::HVec & nn = m_yg[0];
            const dg::HVec & uu = m_yg[1];
            for( unsigned i=0; i<Nx; i++)
            {
                m_density[i] = y[0][i];
                m_velocity[i] = y[1][i];
                unsigned k=i+2;
                yp[0][i] =  -uu[k]*( nn[k+1]-nn[k-1])/2./hx
                            -nn[k]*( uu[k+1]-uu[k-1])/2./hx;
                yp[1][i] =
                    -uu[k]*(uu[k+1]-uu[k-1])/2./hx -
                    m_alpha*(std::pow(nn[k+1], m_gamma)-std::pow(nn[k-1],
                                m_gamma))/2./hx/nn[k]+
                    m_nu_u/nn[k]*(uu[k+1] - 2.*uu[k]+uu[k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "staggered")
        {
            dg::HVec qST(m_yg[1]), q(qST), uh(qST), uST( q), dn(q), du(q);
            const dg::HVec & unST = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            for( unsigned k=0; k<Nx+3; k++)
            {
                // don't compute nSTinv here
                double nST = (nn[k] + nn[k+1])/2.;
                uST[k] = unST[k]/nST;
            }
            for( unsigned k=0; k<Nx+3; k++)
                dn[k] = nn[k+1]-nn[k];
            for( unsigned k=1; k<Nx+4; k++)
                du[k] = uST[k] - uST[k-1];
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
                m_density[i] = y[0][i];
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
                    yp[1][i] += -m_alpha*(  std::pow(nn[k+1], m_gamma) -
                                            std::pow(nn[k], m_gamma))/hx;
                    // This is a bad idea (it will lead to systematic
                    // unphysical acceleration):
                    //yp[1][i] += -m_alpha*m_gamma*std::pow( nn[k], m_gamma-1)*(
                    //        nn[k+1] -nn[k])/hx;
                    // It is bad because nn[k] is used in the pow function instead
                    // of (nn[k+1]+nn[k])/2;
                }
            }
        }
        else if ( m_scheme == "velocity-staggered")
        {
            dg::HVec qST(m_yg[1]), uu(qST), nST(qST), fh(qST), u2(qST), dn(u2), du(u2), du2(u2), q(u2), uh(u2);
            const dg::HVec & uST = m_yg[1];
            const dg::HVec & nn = m_yg[0];

            for( unsigned k=0; k<Nx+3; k++)
            {
                nST[k] = 0.5*(nn[k+1]+nn[k]);
                dn[k] = nn[k+1]-nn[k];
            }
            for ( unsigned k=1; k<Nx+3; k++)
            {
                // This one works slightly better (needs less timesteps)
                uu[k] = 0.5*(uST[k]*nST[k]+uST[k-1]*nST[k-1])/nn[k];
                //uu[k] = 0.5*(uST[k]+uST[k-1]);
                u2[k] = uST[k]*uST[k]/2.;
            }
            for( unsigned k=1; k<Nx+2; k++)
            {
                qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                qST[k]*= uST[k]; // k + 1/2
            }
            for( unsigned k=1; k<Nx+4; k++)
                du2[k] = u2[k]-u2[k-1];
            for ( unsigned k=2; k<Nx+3; k++)
            {
                fh[k] = upwind( uu[k], u2[k-1], u2[k]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    fh[k] += limiter(uu[k], du2[k-1], du2[k], du2[k+1], 0.5, 0.5);
            }
            for ( unsigned k=1; k<Nx+3; k++)
                q [k] = 0.5*(qST[k]+qST[k-1]);
            for( unsigned k=1; k<Nx+4; k++)
                du[k] = uST[k] - uST[k-1];
            for ( unsigned k=2; k<Nx+3; k++)
            {
                uh[k] = upwind( uu[k], uST[k-1], uST[k]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] += limiter(uu[k], du[k-1], du[k], du[k+1], 0.5, 0.5);
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_density[i] = y[0][i];
                m_velocity[i] = 0.5*(uST[k]+uST[k-1]);
                // MW: factoring this out into two terms n d.v + v.d n makes
                // it much worse
                yp[0][i] = -( qST[k] - qST[k-1])/hx;
                // Does not really matter:
                double nSTinv = 2./(nn[k] + nn[k+1]);
                //double nSTinv = (1./nn[k] + 1./nn[k+1])/2.;

                // This one is far away from good
                //yp[1][i] = -uST[k]*(uh[k+1]-uh[k])/hx;
                // This one has the "gathering velocity" problem
                //yp[1][i] = -uu[k]*(uh[k+1]-uh[k])/hx;
                // This one captures Burger and almost fluid; almost overlays
                // in wave (has overshoots in the shocks, better than fh)
                yp[1][i] = -(uu[k+1]*uh[k+1]-uu[k]*uh[k])/2./hx;
                // This one almost captures the fluid shock but not Burger
                // needs a lot of timesteps
                // especially w/o slope-limiter
                //yp[1][i] = -(q[k+1]+q[k])/2.*(uh[k+1]-uh[k])/hx*nSTinv;
                // This is almost the same as with fh
                //yp[1][i] = -(uh[k+1]*uh[k+1]-uh[k]*uh[k])/2./hx;
                // This is the most straightforward, works all-round but does
                // not capture shock (but does capture Burger's shock)
                //yp[1][i] = -(fh[k+1]-fh[k])/hx;
                yp[1][i]+= m_nu_u*nSTinv*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx;
            }
            if( m_variant == "explicit" || m_variant == "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    //double nSTinv = 2./(nn[k] + nn[k+1]);
                    // This version is much better than the first or second
                    double nSTinv = (1./nn[k] + 1./nn[k+1])/2.;
                    //double nSTinv = 1./sqrt( nn[k]*nn[k+1]);
                    if( m_gamma == 1)
                        // log is worse the nSTinv
                        //yp[1][i] += -m_alpha*(  log(nn[k+1]) - log(nn[k]))/hx;
                        yp[1][i] += -m_alpha*(  nn[k+1] - nn[k])*nSTinv/hx;
                        //yp[1][i] += -m_alpha*(  nn[k+1]*nn[k+1] - nn[k]*nn[k])*(nn[k+1]*nn[k+1]+nn[k]*nn[k])/nn[k+1]/nn[k]/nn[k+1]/nn[k]/4./hx;
                    else
                        //yp[1][i] += -m_gamma/(m_gamma-1)*m_alpha*(
                        //        pow(nn[k+1], m_gamma-1) - pow(nn[k],
                        //            m_gamma-1))/hx;
                        yp[1][i] += -m_alpha*( pow(nn[k+1], m_gamma)
                                - pow(nn[k], m_gamma))*nSTinv/hx;
                    // possible sources (due to floor level for example) should also
                    // go here
                }
            }
        }
        else if ( m_scheme == "velocity-unstaggered")
        {
            // same as velocity-staggered, just without the staggering
            dg::HVec qST(m_yg[1]), uST(qST), dn(qST), du(qST), uh(qST), u2(qST);
            const dg::HVec & uu = m_yg[1];
            const dg::HVec & nn = m_yg[0];
            for ( unsigned k=0; k<Nx+3; k++)
            {
                uST[k] = 0.5*(uu[k+1]+uu[k]); // this is the local shock speed
                u2[k] = uu[k]*uu[k];
            }
            for( unsigned k=0; k<Nx+3; k++)
                dn[k] = nn[k+1]-nn[k];
            for( unsigned k=1; k<Nx+2; k++)
            {
                qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                qST[k]*= uST[k]; // k + 1/2
            }
            for( unsigned k=0; k<Nx+3; k++)
                du[k] = uu[k+1] - uu[k];
            for ( unsigned k=1; k<Nx+2; k++)
            {
                uh[k] = upwind( uST[k], uu[k], uu[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    uh[k] += limiter(uST[k], du[k-1], du[k], du[k+1], 0.5, 0.5);
                uh[k]*= uST[k]; // k + 1/2

            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_density[i] = y[0][i];
                m_velocity[i] = y[1][i];
                yp[0][i] = -( qST[k] - qST[k-1])/hx;

                yp[1][i] = -(uh[k] - uh[k-1])/2./hx
                //yp[1][i] = -uu[k]*(uu[k+1] - uu[k-1])/2./hx
                    - m_alpha*(std::pow(nn[k+1], m_gamma)-std::pow(nn[k-1],
                                m_gamma))/2./hx/nn[k]
                    + m_nu_u/nn[k]*(uu[k+1]-2.*uu[k]+uu[k-1])/hx/hx;
            }
        }
        else if ( m_scheme == "staggered-direct")
        {
            //This is the Stegmeir variant
            // y[0] -> ln n
            // y[1] -> u^st
            dg::HVec qST(m_yg[1]), dlnST(qST), uu(qST), q(qST);
            const dg::HVec & uST = m_yg[1];
            const dg::HVec & lnn = m_yg[0];
            for ( unsigned k=1; k<Nx+3; k++)
                uu[k] = 0.5*(uST[k]+uST[k-1]); // this is the local shock speed
            for( unsigned k=0; k<Nx+3; k++)
                dlnST[k] = lnn[k+1]-lnn[k];
            for( unsigned k=1; k<Nx+2; k++)
                qST[k] = uST[k]*dlnST[k];
            for ( unsigned k=1; k<Nx+3; k++)
                q [k] = 0.5*(qST[k]+qST[k-1]);
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                m_density[i] = exp(y[0][i]);
                m_velocity[i] = uu[k];
                yp[0][i] = -( uST[k]-uST[k-1] )/hx  - q[k]/hx ;
                yp[1][i] = -uST[k]*(uu[k+1]-uu[k])/hx;
                double nSTinv = 2./(exp(lnn[k]) + exp(lnn[k+1]));
                yp[1][i]+= m_nu_u*nSTinv*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx;
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                if( m_gamma == 1)
                    yp[1][i] += -m_alpha*(  lnn[k+1] - lnn[k])/hx;
                else
                {
                    double nSTinv = (1./exp(lnn[k]) + 1./exp(lnn[k+1]))/2.;
                    yp[1][i] += -m_alpha*( pow(exp(lnn[k+1]), m_gamma)
                            - pow(exp(lnn[k]), m_gamma))*nSTinv/hx;
                }
            }
        }
        else if ( m_scheme == "log-staggered")
        {
            dg::HVec qST(m_yg[1]), uu(qST), fh(qST), u2(qST), dlnn(u2), du(u2), du2(u2), uh(u2), dn(u2), nn(u2);
            const dg::HVec & lnn = m_yg[0];
            const dg::HVec & uST = m_yg[1];
            for ( unsigned k=1; k<Nx+3; k++)
            {
                nn[k] = exp( lnn[k]);
                uu[k] = 0.5*(uST[k]+uST[k-1]); // this is the local shock speed
                u2[k] = uST[k]*uST[k]/2.;
            }
            for( unsigned k=0; k<Nx+3; k++)
                dlnn[k] = lnn[k+1]-lnn[k];
                //dn[k] = nn[k+1]-nn[k];
            for( unsigned k=1; k<Nx+2; k++)
            {
                qST[k] = upwind( uST[k], lnn[k], lnn[k+1]);
                if( m_variant == "slope-limiter-explicit" || m_variant ==
                        "slope-limiter")
                    qST[k] += limiter(uST[k], dlnn[k-1], dlnn[k], dlnn[k+1], 0.5, 0.5);
                // MW: this works better in 1d but in 3d does not really work well
                //qST[k] = upwind( uST[k], nn[k], nn[k+1]);
                //if( m_variant == "slope-limiter-explicit" || m_variant ==
                //        "slope-limiter")
                //    qST[k] += limiter(uST[k], dn[k-1], dn[k], dn[k+1], 0.5, 0.5);
                //qST[k]*= uST[k]; // k + 1/2
            }
            for( unsigned k=1; k<Nx+4; k++)
                du[k] = uST[k] - uST[k-1];
            for ( unsigned k=2; k<Nx+3; k++)
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
                m_velocity[i] = uu[k];
                yp[0][i] = -uu[k]*( qST[k] - qST[k-1])/hx - du[k]/hx;
                //yp[0][i] = -( qST[k] - qST[k-1])/nn[k]/hx;
                double nSTinv = 2./(nn[k] + nn[k+1]);
                yp[1][i] = -(uu[k+1]*uh[k+1]-uu[k]*uh[k])/2./hx;
                yp[1][i]+= m_nu_u*nSTinv*(uST[k+1] - 2.*uST[k] + uST[k-1]) /hx/hx;
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                if( m_gamma == 1)
                    yp[1][i] += -m_alpha*(  lnn[k+1] - lnn[k])/hx;
                else
                {
                    double nSTinv = (1./exp(lnn[k]) + 1./exp(lnn[k+1]))/2.;
                    yp[1][i] += -m_alpha*( pow(exp(lnn[k+1]), m_gamma)
                            - pow(exp(lnn[k]), m_gamma))*nSTinv/hx;
                }
            }
        }
        for( unsigned i=0; i<Nx; i++)
        {
            unsigned k=i+2;
            if( m_scheme == "staggered-direct" || m_scheme == "log-staggered")
                yp[0][i] +=  m_nu_n*( exp(m_yg[0][k+1]) - 2.*exp(m_yg[0][k]) +
                        exp(m_yg[0][k-1]))/exp(m_yg[0][k])/hx/hx;
            else
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
                if( m_scheme == "staggered-direct" || m_scheme == "log-staggered")
                {
                    for( unsigned k=0; k<Nx; k++)
                        tmpN[k]/=exp(y[0][k]);
                }
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
        bool burger = true;
        for( unsigned i=1; i<Nx; i++)
            if( fabs(y[0][i] - y[0][0]) > 1e-14){ burger  = false; break;}
        if( burger)
            for( unsigned i=0; i<Nx; i++)
                yp[0][i] = 0.;
    }
    unsigned called() const { return m_called;}
    private:
    std::string m_scheme, m_variant, m_init;
    std::array<dg::HVec,2> m_yg;
    dg::HVec m_velocity; // stores the velocity on non-staggered grid
    dg::HVec m_density;
    dg::Grid1d m_g, m_vel_g;
    dg::bc m_bc_n, m_bc_u;
    double m_alpha, m_gamma, m_nu_u, m_nu_n;
    double m_n0 = 0, m_u0 = 0, m_A = 0, m_B = 0, m_k = 0, m_v = 0;
    unsigned m_called = 0;
};

struct NavierStokesImplicit
{
    NavierStokesImplicit( NavierStokesExplicit& ex) : m_ex(ex){}

    void operator() ( double t, const std::array<dg::HVec,2> & y, std::array<dg::HVec,2>& yp)
    {
        dg::blas1::copy( 0., yp);
        unsigned Nx = m_ex.m_g.N();
        double hx = m_ex.m_g.h();
        assign_ghost_cells( y[0], m_ex.m_yg[0], m_ex.m_bc_n);
        assign_ghost_cells( y[1], m_ex.m_yg[1], m_ex.m_bc_u);
        // ghost cells are shifted by 2
        if( m_ex.m_scheme == "staggered")
        {
            if( m_ex.m_variant != "explicit" && m_ex.m_variant !=
                    "slope-limiter-explicit")
            {
                for( unsigned i=0; i<Nx; i++)
                {
                    unsigned k=i+2;
                    yp[1][i] = -m_ex.m_alpha*(std::pow(m_ex.m_yg[0][k+1],
                                m_ex.m_gamma) - std::pow(m_ex.m_yg[0][k],
                                    m_ex.m_gamma))/hx;
                    // possible sources (due to floor level for example) should also
                    // go here
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
                    double gamma = m_ex.m_gamma;
                    if( gamma == 1)
                        //yp[1][i] += -m_ex.m_alpha*(  log(m_ex.m_yg[0][k+1]) -
                        //        log(m_ex.m_yg[0][k]))/hx;
                        yp[1][i] += -m_ex.m_alpha*(  m_ex.m_yg[0][k+1] - m_ex.m_yg[0][k])*(1./m_ex.m_yg[0][k+1]+1./m_ex.m_yg[0][k])/2./hx;
                    else
                        //yp[1][i] += -gamma/(gamma-1)*m_ex.m_alpha*(
                        //        pow(m_ex.m_yg[0][k+1], gamma-1) - pow(m_ex.m_yg[0][k],
                        //            gamma-1))/hx;
                        yp[1][i] += -m_ex.m_alpha*( pow(m_ex.m_yg[0][k+1],
                                    gamma) - pow(m_ex.m_yg[0][k],
                                        gamma))*(1./m_ex.m_yg[0][k+1]+1./m_ex.m_yg[0][k])/2./hx;
                    // possible sources (due to floor level for example) should also
                    // go here
                }
            }
        }
    }
    private:
    NavierStokesExplicit& m_ex;
};

struct NavierStokesImplicitSolver
{
    NavierStokesImplicitSolver( dg::Grid1d g, Json::Value js) :
        m_tmp( {dg::HVec(g.size(), 0.0), dg::HVec ( g.size(), 0.)}){}
    const std::array<dg::HVec,2>& copyable() const{
        return m_tmp;
    }
    // solve (y + alpha I(t,y) = rhs
    void solve( double alpha, NavierStokesImplicit& im, double t, std::array<dg::HVec,2>& y, const std::array<dg::HVec,2>& rhs)
    {
        dg::blas1::copy( rhs[0], y[0]);// I_n = 0
        im( t, y, m_tmp); //ignores y[1],
        // writes 0 in m_tmp[0] and updates m_tmp[1]
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
            else if ("riemann" == init)
            {
                double alpha = v.js["physical"].get("alpha",2).asDouble();
                if( alpha == 0)
                {
                    double x_a = v.js["init"].get("x_a", 0.1).asDouble();
                    double u_l = v.js["init"].get("u_l", 1.0).asDouble();
                    double u_r = v.js["init"].get("u_r", 1.0).asDouble();
                    result = dg::evaluate( [=](double x){
                            if( x <= x_a + (u_l+u_r)/2.*v.time)
                                return u_l;
                            return u_r; }, v.grid);

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
    dg::Grid1d vel_grid = equations::createGrid( js["grid"], dg::PER);
    if ( "staggered" == scheme || "velocity-staggered" == scheme ||
            "staggered-direct" == scheme || "log-staggered" == scheme)
        vel_grid = equations::createStaggeredGrid( js["grid"], dg::PER);
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
    if( scheme == "staggered-direct" || scheme == "log-staggered")
        dg::blas1::transform( y0[0], y0[0], dg::LN<double>());

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
    equations::NavierStokesExplicit ex( grid, vel_grid, js);
    equations::NavierStokesImplicit im( ex);

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
    equations::Variables var = {ex, grid, y0, time, js, 0., 0};
    {   // update internal velocity
        std::array<dg::HVec , 2> tmp ( y0);
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
            if( dt < 1e-6)
                throw dg::Error(dg::Message(_ping_)<<"Adaptive failed to converge! dt = "<<std::scientific<<dt);
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
    //            ark_fixed.step( ex, im, time, y0, time, y0, dt, delta);
    //    else if( timestepper == "ERK")
    //        for( unsigned k=0; k<NT; k++)
    //            erk_fixed.step( ex, time, y0, time, y0, dt);
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
