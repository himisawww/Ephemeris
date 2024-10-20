#include"math/interp.h"
#include"utils/logger.h"

static double max_relative_error=0;

#define TEST_EPSILON_BSPLINE    2e-15
#define ANGULAR_FREQUENCY       0.1
#define TEST_EPSILON_EXPAND     4e-15

struct bspline_fit_tests{
    int_t d,n,mn;
    double reference_err;

    double do_test() const;
};

static constexpr bspline_fit_tests all_tests[]={
    {   1,     1,     1, 1e-2 },
    {   1,     1,     2, 2e-2 },
    {   1,    20,    20, 1e-2 },
    {   2,     1,     2, 1e-3 },
    {   2,    19,    20, 2e-4 },
    {   2,    10,    20, 1e-3 },
    {   3,     1,     3, 5e-5 },
    {   3,    10,    20, 2e-5 },
    {   4,     1,     4, 5e-6 },
    {   4,   160,  1000, 5e-5 },
    {   5,   200,  1000, 2e-6 },
    {   6,   220,  1000, 2e-7 },
    {   7,   240,  1000, 5e-9 },
    {   8,   200,  1000, 3e-9 },
    {   9,   180,  1000, 1e-9 },
    {  10,   160,  1000, 5e-10},
    {  11,  1234, 10000, 1e-9 },
    {  13,  2500, 10000, 1e-11},
    {  16, 12340,100000, 1e-11},
    {  19, 25000,100000, 4e-8 },
};

double bspline_fit_tests::do_test() const{
    std::vector<double> xy;
    xy.reserve(2*(mn+1));
    for(int_t i=0;i<=mn;++i){
        xy.push_back(std::sin(i*ANGULAR_FREQUENCY));
        xy.push_back(std::cos(i*ANGULAR_FREQUENCY));
    }

    bspline_fitter<double,2> bfitter(d,n,xy.data(),mn);
    bspline_fitter<double,2> cfitter(bfitter.get_fitted_data(),d,n,mn);
    cfitter.expand();

    double max_expand_err=0;
    double max_fit_err=0;
    double failed_at;
    bool failed=false;
    auto ltest=[&](double x){
        double f[2],df[2],fd[2],bdf[2];
        for(int k=0;k<2;++k){
            auto &fitter=k==0?bfitter:cfitter;
            fitter(x,f);
            if(k==1){
                checked_maximize(max_expand_err,std::abs(f[0]-fd[0]));
                checked_maximize(max_expand_err,std::abs(f[1]-fd[1]));
                bdf[0]=df[0];
                bdf[1]=df[1];
            }
            fitter(x,fd,df);
            if(f[0]!=fd[0]||f[1]!=fd[1]){
                failed=true;
                failed_at=x;
            }
            if(k==1){
                checked_maximize(max_expand_err,std::abs(bdf[0]-df[0])/ANGULAR_FREQUENCY);
                checked_maximize(max_expand_err,std::abs(bdf[1]-df[1])/ANGULAR_FREQUENCY);
            }
            double fx=std::sin(x*ANGULAR_FREQUENCY);
            double fy=std::cos(x*ANGULAR_FREQUENCY);
            checked_maximize(max_fit_err,std::abs(fx-f[0]));
            checked_maximize(max_fit_err,std::abs(fy-f[1]));
            checked_maximize(max_fit_err,std::abs(fy*ANGULAR_FREQUENCY-df[0]));
            checked_maximize(max_fit_err,std::abs(fx*ANGULAR_FREQUENCY+df[1]));
        }
    };

    for(int_t i=0;i<=mn;++i)ltest(i);
    for(int_t i=0;i<=2*n;++i)ltest(i*mn/double(2*n));

    if(failed){
        LogError(
            "\nInconsistent Result at (%lld, %lld)[%lld](%.16le)",
            d,n,mn,failed_at);
        return 2;
    }
    if(!(max_expand_err<TEST_EPSILON_EXPAND)){
        LogError(
            "\nAt (%lld, %lld)[%lld], "
            "\nMax Expand Error %.16le Too Large",
            d,n,mn,
            max_expand_err);
        return 3;
    }
    if(!(max_fit_err<reference_err)){
        LogError(
            "\nAt (%lld, %lld)[%lld], "
            "\nMax Fit Error %.16le Too Large (Ref. %.1le)",
            d,n,mn,
            max_fit_err,reference_err);
        return 4;
    }

    checked_maximize(max_relative_error,max_expand_err/TEST_EPSILON_EXPAND);
    checked_maximize(max_relative_error,max_fit_err/reference_err);
    return 0;
}

int test_bspline(){
    double max_bspline_err=0;
    for(int_t d=0;d<=bspline_basis_max_degree();++d){
        for(int_t xi=-d-1;xi<=1;++xi){
            double xf=0;
            int_t xfi=-20;
            do{
                double x=xi+xf;
                double b,db,bc,dbc;
                b=bspline_basis(d,x,&db);
                bc=bspline_basis_chebyshev(d,x,&dbc);
                double err_fac=1;
                if(xi>-d&&xi<0)
                    checked_maximize(err_fac,1/std::sqrt(b));
                checked_maximize(max_bspline_err,err_fac*std::abs(b-bc));
                checked_maximize(max_bspline_err,err_fac*std::abs(db-dbc));
                double exfi=std::exp(++xfi<<1);
                xf=exfi/(1+exfi);
            } while(xf!=1);
        }
    }
    if(!(max_bspline_err<TEST_EPSILON_BSPLINE)){
        LogError(
            "\nMax BSpline Error %.16le Too Large",
            max_bspline_err);
        return 1;
    }

    for(const auto &t:all_tests){
        int ret=t.do_test();
        if(ret)return ret;
    }

    LogInfo("\n      Passed(%f), ",max_relative_error);
    return 0;
}
