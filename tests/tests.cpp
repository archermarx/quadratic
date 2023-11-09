#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <stdlib.h>
#include <random>

#include "doctest.h"
#include "quadratic.h"

#include <iostream>

// Approximate equality of floating point numbers, including comparing if both are nan or inf.
template <typename T>
bool isapprox(T x, T y) {
    return  (std::isnan(x) && std::isnan(y)) ||
            (std::isinf(x) && std::isinf(y) && sgn(x) == sgn(y)) ||
            std::abs(x - y) <= std::sqrt(eps<T>) * std::max(abs(x), abs(y));
}

TEST_CASE("constants") {

    // constants

    // Max and minimum exponents
    CHECK(E_MIN<float> == -126);
    CHECK(E_MAX<float> == 127);
    CHECK(E_MIN<double> == -1022);
    CHECK(E_MAX<double> == 1023);
    CHECK(E_MIN<long double> == -16382);
    CHECK(E_MAX<long double> == 16383);
     
    // Number of bits in mantissa
    CHECK(NUM_DIGITS<float> == 24);
    CHECK(NUM_DIGITS<double> == 53);
    CHECK(NUM_DIGITS<long double> == 64);

    // Maximum exponent to avoid overflow
    CHECK(ECP_MAX<double> ==  995);

    // Minimum exponent to avoid underflow
    CHECK(ECP_MIN<double> == -920);
}

// Helper class for quadratic tests
template<typename T>
class quadratic_test : public quadratic<T> {
    public:
        quadratic_return_code retcode;
        T expected_x1;
        T expected_x2;

        quadratic_test(T _a, T _b, T _c, quadratic_return_code _retcode, T _x1, T _x2) {
            this->a = _a; this->b = _b; this->c = _c; retcode = _retcode; expected_x1 = _x1; expected_x2 = _x2;
        }

        quadratic_test(T _a, T _b, T _c, quadratic_return_code _retcode, T _x1) {
            this->a = _a; this->b = _b; this->c = _c; retcode = _retcode; expected_x1 = _x1; expected_x2 = NaN<T>;
        }

        quadratic_test(T _a, T _b, T _c, quadratic_return_code _retcode) {
            this->a = _a; this->b = _b; this->c = _c; retcode = _retcode; expected_x1 = NaN<T>; expected_x2 = NaN<T>;
        }

        void validate() {
            auto ret = this->solve();
            CHECK(ret == retcode);

            if (ret != retcode || !isapprox(this->x1, expected_x1) || !isapprox(this->x2, expected_x2)) {
                std::cerr << "x1, x2 (expected): " << expected_x1 << ", " << expected_x2 << std::endl;
                std::cerr << "x1, x2 (solved): " << this->x1 << ", " << this->x2 << std::endl << std::endl;
            }

            switch(retcode) {
                case QUAD_ONE_SOLUTION:
                    CHECK(isapprox(this->x1, expected_x1));
                    break;
                case QUAD_TWO_SOLUTIONS:
                    CHECK(isapprox(this->x1, expected_x1));
                    CHECK(isapprox(this->x2, expected_x2));
                    break;
                default:
                    break;
            }
            return;
        }   
};

int64_t fib(int n) {
    double phi = 0.5 * (1 + sqrt(5));
    double phi_n = pow(phi, n);
    return round(phi_n / sqrt(5));
}

TEST_CASE("fibonacci") {
    int n_max = 20;
    int F = 0;
    int F_old = 1;
    int F_new;
    for (int n = 0; n < n_max; n++) {
        CHECK(F == fib(n));

        F_new = F + F_old;
        F_old = F;
        F = F_new;
    }
}

template <typename T>
quadratic_test<T> kahan_random_quadratic(int m) {
    // only even n have real solutions
    auto n = 2 * m;

    auto f2 = fib(n-2);
    auto f1 = fib(n-1);
    auto f0 = fib(n);

    auto x1 = T(f1 - 1) / f0;
    auto x2 = T(f1 + 1) / f0;

    T min_R = exp2(NUM_DIGITS<T> - 1);
    T max_R = exp2(NUM_DIGITS<T>) - 1;

    std::default_random_engine rng;
    std::uniform_real_distribution<T> dist(min_R, max_R);
    T R = dist(rng);

    T M = floor(R / f0);
    T a = M * f0;
    T b = -2 * M * f1;
    T c = M * f2;

    quadratic_test<T> test(a, b, c, QUAD_TWO_SOLUTIONS, x1, x2);

    return test;
}

TEST_CASE("kahan randomized fibonacci quadratic stress test") {

    int min_n = 1;
    int max_n_single = 17;
    int max_n_double = 35;
    int num_tests_per_n = 5;

    // singles
    for (int i = 0; i < num_tests_per_n; ++i) {
        for (int n = min_n; n <= max_n_single; ++n) {
            kahan_random_quadratic<float>(n).validate();
        }
        // doubles
        for (int n = min_n; n <= max_n_double; ++n) {
            kahan_random_quadratic<double>(n).validate();
        }
    }
}

TEST_CASE("no real solutions") {
    std::vector<quadratic_test<double>> tests{
        quadratic_test<double>(2.0, 0.0, 3.0, QUAD_NO_SOLUTIONS),
        quadratic_test<double>(1.0, 1.0, 1.0, QUAD_NO_SOLUTIONS),
        quadratic_test<double>(exp2(600.0), 0.0, exp2(600.0), QUAD_NO_SOLUTIONS),
        quadratic_test<double>(-exp2(600.0), 0.0, -exp2(600.0), QUAD_NO_SOLUTIONS)
    };

    for (auto &test: tests) {
        test.validate();
    }
}

TEST_CASE("degenerate inputs") {
    
    std::vector<quadratic_test<double>> tests {
        // a==0, b==0, c==0
        quadratic_test<double>(0., 0., 0., QUAD_NO_SOLUTIONS),
        // a==0, b!=0, c==0
        quadratic_test<double>(0., 1., 0., QUAD_ONE_SOLUTION, 0.0),
        // a==0, b!=0, c!=0
        quadratic_test<double>(0., 1., 2., QUAD_ONE_SOLUTION, -2.0),
        quadratic_test<double>(0., 3., 4., QUAD_ONE_SOLUTION, -4./3.),
        quadratic_test<double>(0., exp2(600), -exp2(600), QUAD_ONE_SOLUTION, 1.0),
        quadratic_test<double>(0., exp2(600), exp2(600), QUAD_ONE_SOLUTION, -1.0),
        quadratic_test<double>(0., exp2(-600), exp2(600), QUAD_ONE_SOLUTION, -inf<double>),
        quadratic_test<double>(0., exp2(600), exp2(-600), QUAD_ONE_SOLUTION, 0.0),
        quadratic_test<double>(0., 2., -1.0e-323, QUAD_ONE_SOLUTION, 5.0e-324),
        // a!=0, b==0, c==0
        quadratic_test<double>(3., 0., 0., QUAD_ONE_SOLUTION, 0.0),
        quadratic_test<double>(1200., 0., 0., QUAD_ONE_SOLUTION, 0.0),
        // a!=0, b==0, c!=0
        quadratic_test<double>(2., 0.,-3., QUAD_TWO_SOLUTIONS, -sqrt(3./2.), sqrt(3./2.)),
        quadratic_test<double>(exp2(600), 0., -exp2(600), QUAD_TWO_SOLUTIONS, -1.0, 1.0),
        // a!=0, b!=0, c==0
        quadratic_test<double>(3., 2., 0., QUAD_TWO_SOLUTIONS, -2.0/3.0, 0.0),
        quadratic_test<double>(exp2(600), exp2(700), 0., QUAD_TWO_SOLUTIONS, -exp2(100), 0.0),
        quadratic_test<double>(exp2(-600), exp2(700), 0.0, QUAD_TWO_SOLUTIONS, -inf<double>, 0.0),
        quadratic_test<double>(exp2(600), exp2(-700), 0.0, QUAD_TWO_SOLUTIONS, -0.0, 0.0)
    };

    for (auto &test : tests) {
        test.validate();
    }
}

TEST_CASE("one solution") {
    std::vector<quadratic_test<double>> tests{
        quadratic_test<double>(exp2(1023), 2.0, exp2(-1023), QUAD_ONE_SOLUTION, -exp2(-1023))
    };

    for (auto &test : tests) {
        test.validate();
    }
}

TEST_CASE("two solutions") {
    auto ret = QUAD_TWO_SOLUTIONS;
    std::vector<quadratic_test<double>> tests{
        quadratic_test<double>(1., -1., -1., ret, -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(1., 1 + exp2(-52), 0.25 + exp2(-53), ret, (-1-exp2(-51))/2, -0.5),
        quadratic_test<double>(1., exp2(-511) + exp2(-563), exp2(-1024), ret, -7.458340888372987e-155,-7.458340574027429e-155),
        quadratic_test<double>(1., exp2(27), 0.75, ret, -134217728.0, -5.587935447692871e-09),
        quadratic_test<double>(1., -1e9, 1., ret, 1e-9, 1000000000.0),
        quadratic_test<double>(1.3407807929942596e154, -1.3407807929942596e154, -1.3407807929942596e154, ret, -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(exp2(600), 0.5, -exp2(-600), ret, -3.086568504549085e-181,1.8816085719976428e-181),
        quadratic_test<double>(exp2(600), 0.5, -exp2(600), ret, -1.0, 1.0),
        quadratic_test<double>(8.0, exp2(800), -exp2(500), ret, -8.335018041099818e+239, 4.909093465297727e-91),
        quadratic_test<double>(1.0, exp2(26), -0.125, ret, -67108864.0, 1.862645149230957e-09),
        quadratic_test<double>(exp2(-1073), -exp2(-1073), -exp2(-1073), ret, -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(exp2(600), -exp(-600), -exp2(-600), ret, -2.409919865102884e-181, 2.409919865102884e-181),    
        quadratic_test<double>(-158114166017., 316227766017., -158113600000., ret, 0.99999642020057874,1.0),
        quadratic_test<double>(-312499999999.0,707106781186.0,-400000000000.0, ret, 1.131369396027,1.131372303775),
        quadratic_test<double>(-67., 134., -65., ret, 0.82722631488372798,1.17277368511627202),
        quadratic_test<double>(0.247260273973, 0.994520547945, -0.138627953316, ret, -4.157030027041105,0.1348693622211607),
        quadratic_test<double>(1., -2300000., 2e11, ret,90518.994979145,2209481.005020854),
        quadratic_test<double>(1.5*exp2(-1026), 0., -exp2(1022), ret, -1.4678102981723264e308, 1.4678102981723264e308)
    };

    for (auto &test : tests) {
        test.validate();
    }
}