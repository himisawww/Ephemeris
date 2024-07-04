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

/*
    Fit T data[mn+1][N_Channel] with n+d d-th-degree bspline basis functions.
    Coefficients will be stored to T result[n+d][N_Channel].
    need n>=1 && d>0 && n+d<=mn+1
*/
template<typename T,size_t N_Channel>
bool interp_bspline(T *result,const int_t d,const int_t n,const T *data,const int_t mn);
