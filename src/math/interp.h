#include"definitions.h"

// maximum implemented degree d for bspline_* series functions
int_t bspline_basis_max_degree();

/*
    Basis Functions of BSpline(B).

    defined as:
        B[0,x] = x is in [0,1) ? 1 : 0;
        B[d,x] = ((x+d)*B[d-1,x+1]+(1-x)*B[d-1,x])/d;
    for 0<=d<=bspline_basis_max_degree().

    B[d,x] is a C^{d-1} continuous non-negative function of x on Reals,
        it is positive when x is in (-d,1), is 0 when x is not in [-d,1).
*/
double bspline_basis(int_t d,double x);

/*
    The Chebyshev Coefficients(c) of BSpline Basis Functions.
    for real x in [0,1) and integer k in [0,d],
        B[d,x-k] = sum_{j=0}^{d} c(d,k,j)*chebyshev_T(j,2*x-1);
*/
double bspline_basis_chebyshev_coef(int_t d,int_t k,int_t j);

/*
    same as bspline_basis but use Chebyshev expansion for calculation,
    faster but loses relative accuracy when x is near boundaries of [-d,1).
*/
double bspline_basis_chebyshev(int_t d,double x);
