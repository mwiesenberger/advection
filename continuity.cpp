#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/file/file.h"
#include "common.h"

namespace equations
{

struct Filter
{
    Filter( dg::Grid1d g, dg::file::WrappedJsonValue js) : m_tmp( g.size())
    {
        m_scheme = js["advection"].get("type", "upwind").asString();
        if( m_scheme == "dg-upwind")
        {
            m_alpha  = js["advection"].get("alpha", 10).asDouble()*g.h()*g.h();
            m_stencil = dg::create::limiter_stencil( g, g.bcx());
        }
        else if( m_scheme == "dg-upwind-swm")
        {
            m_stencil = dg::create::window_stencil(3, g, g.bcx());
            m_alpha  = js["advection"].get("alpha", 10).asDouble()*g.h()*g.h();
            m_iter = js["advection"].get("iter", 4).asUInt();
        }
    }
    void operator()( dg::HVec& y)
    {
        if( m_scheme == "dg-upwind")
        {
            dg::blas2::stencil( dg::CSRSlopeLimiter<double>(m_alpha), m_stencil, y, m_tmp);
            m_tmp.swap(y);
        }
        else if( m_scheme == "dg-upwind-swm")
        {
            double alpha = m_alpha;
            for( unsigned i=0; i<m_iter; i++)
            {
                dg::blas2::stencil( dg::CSRSWMFilter<double>(alpha), m_stencil, y, m_tmp);
                using std::swap;
                swap( m_tmp, y);
                alpha*=0.8;
            }
        }
    }
    private:
    dg::IHMatrix m_stencil;
    dg::HVec m_tmp;
    std::string m_scheme;
    unsigned m_iter;
    double m_alpha;

};

struct Continuity
{
    Continuity( dg::Grid1d g, dg::file::WrappedJsonValue js) :
        m_yg( g.size()+4, 0.), m_tmp0(g.size()), m_tmp1(g.size()), m_g(g)
    {
        m_scheme = js["advection"].get("type", "upwind").asString();
        m_vel = js["physical"].get("velocity", 1.0).asDouble();
        m_nu = js["physical"].get( "nu", 0.).asDouble();
        m_exp = js["physical"].get( "exponent", 1.).asDouble();
        std::cout << "Compute scheme "<<m_scheme<<"\n";
        m_forward = dg::create::dx( g, g.bcx(), dg::forward);
        m_backward = dg::create::dx( g, g.bcx(), dg::backward);
        m_ell = { g};
    }

    void operator() ( double t, const dg::HVec & y, dg::HVec& yp)
    {
        m_called++;
        unsigned Nx = m_g.N();
        double hx = m_g.h();
        assign_ghost_cells( y, m_yg, m_g.bcx());
        // ghost cells are shifted by 2
        dg::blas1::copy( 0., yp);
        if ( m_scheme == "upwind")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                if( m_vel > 0.)
                    yp[i] -= m_vel*( m_yg[k] - m_yg[k-1])/ hx;
                else
                    yp[i] -= m_vel*( m_yg[k+1] - m_yg[k])/ hx;
            }
        }
        else if ( m_scheme == "upwind2")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                if( m_vel > 0.)
                    yp[i] -= m_vel*( 3*m_yg[k] - 4*m_yg[k-1] + m_yg[k-2])/2./hx;
                else
                    yp[i] -= m_vel*( -m_yg[k+2] + 4*m_yg[k+1] - 3*m_yg[k])/2./hx;
            }
        }
        else if( m_scheme == "diesel")
        {
            dg::HVec tmp(m_yg);
            for( unsigned k=1; k<Nx+3; k++)
            {
                if( fabs(m_yg[k+1]) < fabs(m_yg[k]))
                    tmp[k] = m_vel*m_yg[k];
                else if( fabs(m_yg[k+1]) == fabs(m_yg[k]))
                    tmp[k] = 0;
                else
                    tmp[k] = -m_vel*m_yg[k+1];
            }
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                yp[i] -= (tmp[k]-tmp[k-1])/hx;
            }
        }
        else if ( m_scheme == "centered")
        {
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                yp[i] -= m_vel*( m_yg[k+1]-m_yg[k-1])/2./hx;
            }
        }
        else if (m_scheme.find("dg-upwind") != std::string::npos)
        {
            dg::blas2::symv( m_forward, y, m_tmp0);
            dg::blas2::symv( m_backward, y, m_tmp1);
            dg::blas1::evaluate( yp, dg::minus_equals(), dg::Upwind(), m_vel,
                    m_tmp1, m_tmp0);
        }
        if (m_scheme.find("dg-upwind") != std::string::npos)
        {
            for( unsigned i=0; i<Nx; i++)
                m_tmp0[i] = pow( y[i], m_exp);
            dg::blas2::symv( -m_nu, m_ell, m_tmp0, 1., yp);
        }
        else
        {
            for( unsigned i=0; i<Nx; i++)
            {
                unsigned k=i+2;
                yp[i] +=  m_nu*( pow( m_yg[k+1], m_exp) -2.*pow(m_yg[k], m_exp) + pow( m_yg[k-1], m_exp))/hx/hx;
            }
        }
    }
    unsigned called() const { return m_called;}
    private:
    unsigned m_called = 0;
    std::string m_scheme;
    dg::HVec m_yg, m_tmp0, m_tmp1;
    dg::HMatrix m_forward, m_backward;
    dg::Elliptic1d<dg::Grid1d, dg::HMatrix, dg::HVec> m_ell;
    dg::Grid1d m_g;
    double m_vel, m_nu, m_exp;

};

struct Variables{
    Continuity& f;
    const dg::Grid1d& grid;
    const dg::HVec& y0;
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
    {"numerical", "Numerical solution: Advected variable",
        []( dg::HVec& result, Variables& v ) {
             dg::blas1::copy(v.y0, result);
        }
    },
    {"analytical", "Analytical solution",
        []( dg::HVec& result, Variables& v ) {
            std::string init = v.js[ "init"].get( "type", "step").asString();
            double vel = v.js["physical"].get( "velocity", 1.0).asDouble();
            if( "step" == init)
            {
                double x_a = v.js["init"].get("x_a", 0.1).asDouble();
                result = dg::evaluate( dg::Heaviside(
                            x_a*v.grid.lx()+vel*v.time, -1), v.grid);
            }
            else if( "wave" == init)
            {
                double n_0 = v.js["init"].get("n_0", 1.0).asDouble();
                double amp = v.js["init"].get("amp", 1.0).asDouble();
                double k   = v.js["init"].get("k", 1.0).asDouble();
                double x_0 = v.js["init"].get("x_0", 1.0).asDouble();
                result = dg::evaluate( [=]( double x){ return n_0 + amp*sin(
                            k*(x-x_0-vel*v.time));}, v.grid);
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
    dg::bc bcx = dg::str2bc ( js["grid"].get("bc", "NEU").asString());
    dg::Grid1d grid = equations::createGrid( js, bcx);
    dg::HVec w1d( dg::create::weights(grid));
    /////////////////////////////////////////////////////////////////
    std::string init = js["init"].get("type", "step").asString();
    dg::HVec y0 = dg::evaluate( dg::zero, grid);
    if( "step" == init)
    {
        double x_a = js["init"].get("x_a", 0.1).asDouble();
        y0 = dg::evaluate( dg::Heaviside( x_a*grid.lx(), -1), grid);
    }
    else if( "wave" == init)
    {
        double n_0 = js["init"].get("n_0", 1.0).asDouble();
        double amp = js["init"].get("amp", 1.0).asDouble();
        double k   = js["init"].get("k", 1.0).asDouble();
        double x_0 = js["init"].get("x_0", 1.0).asDouble();
        y0 = dg::evaluate( [=]( double x){ return n_0 + amp*sin( k*(x-x_0));},
                grid);
    }
    else if( "gaussian" == init)
    {
        double amp = js["init"].get("amp", 1.0).asDouble();
        double x_0 = js["init"].get("x_0", 1.0).asDouble();
        double sigma = js["init"].get("sigma", 1.0).asDouble();
        y0 = dg::evaluate( dg::GaussianX( x_0*grid.lx(), sigma*grid.lx(), amp), grid);
    }

    std::string tableau= js["timestepper"].get("tableau",
            "Bogacki-Shampine-4-2-3").asString();
    double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
    double atol= js["timestepper"].get("atol", 1e-6).asDouble();
    double tend = js["output"].get( "tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = tend/(double)maxout;
    dg::Adaptive<dg::FilteredERKStep<dg::HVec>> adaptive( tableau,y0);
    double time = 0.;
    equations::Continuity rhs( grid,js);
    equations::Filter filter( grid,js);
    dg::AdaptiveTimeloop<dg::HVec> timeloop( adaptive,
            std::tie(rhs,filter), dg::pid_control, dg::l2norm, rtol, atol);

    // Set up netcdf
    std::string inputfile = js.asJson().toStyledString(); //save input without comments, which is important if netcdf file is later read by another parser
    std::string outputfile;
    if( argc == 1 || argc == 2)
        outputfile = "equations.nc";
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
    att["title"] = "Output file of advection/equations.cpp";
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
    dg::HVec result = y0;
    equations::Variables var = {rhs, grid, y0, time, js, 0., &adaptive.nfailed()};
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
        timeloop.integrate( time, y0, time+deltaT, y0, dg::to::exact);
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
