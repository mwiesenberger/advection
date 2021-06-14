#pragma once

// common functions to continuity.cpp, navier_stokes.cpp and plasma.cpp
namespace equations
{
dg::Grid1d createGrid( Json::Value grid, dg::bc bcx)
{
    unsigned Nx = grid.get( "Nx", 48).asUInt();
    double x0 = grid["x"].get( 0u, 0.).asDouble();
    double x1 = grid["x"].get( 1u, 1.).asDouble();
    return dg::Grid1d( x0, x1, 1, Nx, bcx);
}

// Actually the staggered grid should have one cell more than the collocated grid?
// | x | x | x | ... | x |
dg::Grid1d createStaggeredGrid( Json::Value grid, dg::bc bcx)
{
    dg::Grid1d g1d = createGrid(grid, bcx);
    unsigned Nx = grid.get( "Nx", 48).asUInt();
    double x0 = g1d.x0() + g1d.h()/2.;
    double x1 = g1d.x1() + g1d.h()/2.;
    return dg::Grid1d( x0, x1, 1, Nx, bcx);
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
        if ( bcx == dg::NEU || bcx == dg::DIR_NEU)
        {
            yg[Nx+2] = y[Nx-1];
            yg[Nx+3] = y[Nx-2];
        }
        if ( bcx == dg::DIR || bcx == dg::NEU_DIR)
        {
            yg[Nx+2] = -y[Nx-1];
            yg[Nx+3] = -y[Nx-2];
        }
        if( bcx == dg::NEU || bcx == dg::NEU_DIR)
        {
            yg[1] = y[0];
            yg[0] = y[1];
        }
        if( bcx == dg::DIR || bcx == dg::DIR_NEU)
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
} //namespace equations
