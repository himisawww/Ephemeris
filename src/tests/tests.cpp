#include"physics/mass.h"
#include"tests/tests.h"
#include"math/random.h"
#include"utils/calctime.h"

struct test_t{
    int (*test)();
    const char *description;
    const char *long_description;
};

int test_all(uint64_t seed){
    std::vector<test_t> tests;
    printf("Starting tests...");
#define DECLARE_TEST(T,S,L) \
int T();                    \
tests.push_back({T,S,L})
    
    if(seed==0)seed=*(uint64_t*)&(const double &)CalcTime();
    seedrandom(seed);

    DECLARE_TEST(test_prepare,      "Preparing",
        "Preparing resources (initials and params of celestial bodies)\n"
        "   for testing procedures.\n"
        );
    DECLARE_TEST(test_geopotential, "Geopotential",
        "Make sure the implementation of geopotential model is correct\n"
        "   by direct comparing to (Associated-) Legendre Polynomials.\n"
    );
    DECLARE_TEST(test_kepler,       "Keplerians",
        "Test convertion between state vectors and orbital keplerian\n"
        "   parameters.\n"
    );
    DECLARE_TEST(test_conservation, "Conservation",
        "Test conservation of momentum and angular momentum.\n"
    );

    int_t N=tests.size();
    for(int_t i=0;i<N;++i){
        test_t &t=tests[i];
        printf("\n[%lld/%lld] Testing %s...",i+1,N,t.description);
        double start_time=CalcTime();
        int result=t.test();
        if(result){
            fprintf(stderr,
                "\n"
                "\n***********************   Test Failed!   ***********************\n"
                "\n   TEST NAME: %s"
                "\n DESCRIPTION:"
                "\n %s"
                "\nFAILURE CODE: %d"
                "\n   TEST SEED: %016llx\n"
                "\n****************************************************************\n"
                "\nPlease contact the author via github for further information.\n\n",
                t.description,t.long_description,result,seed);
            return result;
        }
        else
            printf("Done in %fs",CalcTime()-start_time);
    }
    printf("\nAll Tests Passed.\n");
    return 0;
}
