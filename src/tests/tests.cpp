#include"physics/CelestialSystem.h"
#include"tests/tests.h"

struct test_t{
    int (*test)();
    const char *description;
    const char *long_description;
};

int test_fun(){
    std::vector<test_t> tests;
    printf("Starting tests...");
#define DECLARE_TEST(T,S,L) \
int T();                    \
tests.push_back({T,S,L})

    DECLARE_TEST(test_geopotential, "Geopotential",
        "Make sure the implementation of geopotential model is correct\n"
        "   by direct comparing to (Associated-) Legendre Polynomials.\n"
    );
    int_t N=tests.size();
    for(int_t i=0;i<N;++i){
        test_t &t=tests[i];
        printf("\n[%lld/%lld] Testing %s...",i+1,N,t.description);
        int result=t.test();
        if(result){
            fprintf(stderr,
                "\n"
                "\n***********************   Test Failed!   ***********************\n"
                "\n TEST  NAME: %s"
                "\nDESCRIPTION:\n %s"
                "\n****************************************************************\n",
                t.description,t.long_description);
            return result;
        }
        printf("\nAll Tests Passed.\n");
    }
    return 0;
}
