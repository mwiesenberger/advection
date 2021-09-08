#pragma once

// common functions to continuity.cpp, navier_stokes.cpp and plasma.cpp
namespace equations
{
/**
 * @brief
 \f$ f(x_1, x_2, ...) = \begin{cases}
         \min(x_1, x_2, ...) &\text{ for } x_1, x_2, ... >0 \\
         \max(x_1, x_2, ...) &\text{ for } x_1, x_2, ... <0 \\
         0 &\text{ else}
 \end{cases}
 \f$
 *
 * Useful for Slope limiter
 */
struct MinMod
{
    ///@return minmod(x1, x2)
#ifdef __CUDACC__
    template < class T>
    DG_DEVICE T operator()( T x1, T x2) const
    {
        if( x1 > 0 && x2 > 0)
            return min(x1,x2);
        else if( x1 < 0 && x2 < 0)
            return max(x1,x2);
        return 0.;
    }
#else
    template < class T>
    T operator()( T x1, T x2) const
    {
        if( x1 > 0 && x2 > 0)
            return std::min(x1,x2);
        else if( x1 < 0 && x2 < 0)
            return std::max(x1,x2);
        return 0.;
    }
#endif
    ///@return minmod(x1, x2, x3);
    template<class T>
    DG_DEVICE T operator() ( T x1, T x2, T x3)const
    {
        return this-> operator()( this-> operator()( x1, x2), x3);
    }
};


/**
 * @brief \f$ f(x_1,x_2) = 2\begin{cases}
 *  \frac{x_1x_2}{x_1+x_2} &\text{ if } x_1x_2 > 0 \\
 *  0 & \text { else }
 *  \end{cases}
 *  \f$
 *  @note The first case is the harmonic mean between x_1 and x_2
 */
struct VanLeer
{
    template<class T>
    DG_DEVICE T operator()( T x1, T x2) const
    {
        if( x1*x2 <= 0)
            return 0.;
        return 2.*x1*x2/(x1+x2);
    }
};

/**
 * @brief \f$ \text{up}(v, g_m, g_0, g_p, h_m, h_p ) = \begin{cases}  +h_m \Lambda( g_0, g_m) &\text{ if } v \geq 0 \\
 *  -h_p \Lambda( g_p, g_0) &\text{ else}
 *  \end{cases}
 *  \f$
 *
 * @tparam Limiter Any two-dimensional functor
 * @sa VanLeer, MinMod
 */
template<class Limiter>
struct SlopeLimiter
{
    SlopeLimiter() {}
    SlopeLimiter( Limiter l ) : m_l( l){}
    template<class T>
    DG_DEVICE T operator()( T v, T gm, T g0, T gp, T hm, T hp ) const{
        if( v >= 0)
            return +hm*m_l( g0, gm);
        else
            return -hp*m_l( gp, g0);
    }
    private:
    Limiter m_l;
};
/**
 * @brief \f$ \text{up}(v, g_m, g_0, g_p, h_m, h_p ) = v \begin{cases}  +h_m \Lambda( g_0, g_m) &\text{ if } v \geq 0 \\
 *  -h_p \Lambda( g_p, g_0) &\text{ else}
 *  \end{cases}
 *  \f$
 *
 * @tparam Limiter Any two-dimensional functor
 * @sa VanLeer, MinMod
 */
template<class Limiter>
struct SlopeLimiterProduct
{
    SlopeLimiterProduct() {}
    SlopeLimiterProduct( Limiter l ) : m_s( l){}
    template<class T>
    DG_DEVICE T operator()( T v, T gm, T g0, T gp, T hm, T hp ) const{
        return v*m_s(v,gm,g0,gp,hm,hp);
    }
    private:
    SlopeLimiter<Limiter> m_s;
};
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

} //namespace equations
