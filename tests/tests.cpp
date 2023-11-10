#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <stdlib.h>
#include <random>

#include "doctest.h"
#include "quadratic.h"
#include <chrono>
#include <iostream>

using namespace quadratic;

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
class quadratic_test {
    public:
        T a, b, c, x1, x2;

        quadratic_test() {}
        quadratic_test(T _a, T _b, T _c, T _x1, T _x2) : a(_a), b(_b), c(_c), x1(_x1), x2(_x2) {}
        quadratic_test(T _a, T _b, T _c, T _x1) : a(_a), b(_b), c(_c), x1(_x1), x2(NaN<T>) {}
        quadratic_test(T _a, T _b, T _c) : a(_a), b(_b), c(_c), x1(NaN<T>), x2(NaN<T>) {}

        void validate() {
            auto solution = solve(a, b, c);
            CHECK(isapprox(x1, solution.first));
            CHECK(isapprox(x2, solution.second));
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

    quadratic_test<T> test(a, b, c, x1, x2);

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
        quadratic_test<double>(2.0, 0.0, 3.0),
        quadratic_test<double>(1.0, 1.0, 1.0),
        quadratic_test<double>(exp2(600.0), 0.0, exp2(600.0)),
        quadratic_test<double>(-exp2(600.0), 0.0, -exp2(600.0))
    };

    for (auto &test: tests) {
        test.validate();
    }
}

TEST_CASE("degenerate inputs") {
    
    std::vector<quadratic_test<double>> tests {
        // a==0, b==0, c==0
        quadratic_test<double>(0., 0., 0.),
        // a==0, b!=0, c==0
        quadratic_test<double>(0., 1., 0., 0.),
        // a==0, b!=0, c!=0
        quadratic_test<double>(0., 1., 2., -2.0),
        quadratic_test<double>(0., 3., 4., -4./3.),
        quadratic_test<double>(0., exp2(600), -exp2(600), 1.0),
        quadratic_test<double>(0., exp2(600), exp2(600), -1.0),
        quadratic_test<double>(0., exp2(-600), exp2(600), -inf<double>),
        quadratic_test<double>(0., exp2(600), exp2(-600), 0.0),
        quadratic_test<double>(0., 2., -1.0e-323, 5.0e-324),
        // a!=0, b==0, c==0
        quadratic_test<double>(3., 0., 0., 0.0),
        quadratic_test<double>(1200., 0., 0., 0.0),
        // a!=0, b==0, c!=0
        quadratic_test<double>(2., 0.,-3., -sqrt(3./2.), sqrt(3./2.)),
        quadratic_test<double>(exp2(600), 0., -exp2(600), -1.0, 1.0),
        // a!=0, b!=0, c==0
        quadratic_test<double>(3., 2., 0., -2.0/3.0, 0.0),
        quadratic_test<double>(exp2(600), exp2(700), 0., -exp2(100), 0.0),
        quadratic_test<double>(exp2(-600), exp2(700), 0.0, -inf<double>, 0.0),
        quadratic_test<double>(exp2(600), exp2(-700), 0.0, -0.0, 0.0)
    };

    for (auto &test : tests) {
        test.validate();
    }
}

TEST_CASE("one solution") {
    std::vector<quadratic_test<double>> tests{
        quadratic_test<double>(exp2(1023), 2.0, exp2(-1023), -exp2(-1023))
    };

    for (auto &test : tests) {
        test.validate();
    }
}

TEST_CASE("two solutions") {
    std::vector<quadratic_test<double>> tests{
        quadratic_test<double>(1., -1., -1., -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(1., 1 + exp2(-52), 0.25 + exp2(-53), (-1-exp2(-51))/2, -0.5),
        quadratic_test<double>(1., exp2(-511) + exp2(-563), exp2(-1024), -7.458340888372987e-155,-7.458340574027429e-155),
        quadratic_test<double>(1., exp2(27), 0.75, -134217728.0, -5.587935447692871e-09),
        quadratic_test<double>(1., -1e9, 1., 1e-9, 1000000000.0),
        quadratic_test<double>(1.3407807929942596e154, -1.3407807929942596e154, -1.3407807929942596e154, -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(exp2(600), 0.5, -exp2(-600), -3.086568504549085e-181,1.8816085719976428e-181),
        quadratic_test<double>(exp2(600), 0.5, -exp2(600), -1.0, 1.0),
        quadratic_test<double>(8.0, exp2(800), -exp2(500), -8.335018041099818e+239, 4.909093465297727e-91),
        quadratic_test<double>(1.0, exp2(26), -0.125, -67108864.0, 1.862645149230957e-09),
        quadratic_test<double>(exp2(-1073), -exp2(-1073), -exp2(-1073), -0.6180339887498948, 1.618033988749895),
        quadratic_test<double>(exp2(600), -exp(-600), -exp2(-600), -2.409919865102884e-181, 2.409919865102884e-181),    
        quadratic_test<double>(-158114166017., 316227766017., -158113600000., 0.99999642020057874,1.0),
        quadratic_test<double>(-312499999999.0,707106781186.0,-400000000000.0, 1.131369396027,1.131372303775),
        quadratic_test<double>(-67., 134., -65., 0.82722631488372798,1.17277368511627202),
        quadratic_test<double>(0.247260273973, 0.994520547945, -0.138627953316, -4.157030027041105,0.1348693622211607),
        quadratic_test<double>(1., -2300000., 2e11,90518.994979145,2209481.005020854),
        quadratic_test<double>(1.5*exp2(-1026), 0., -exp2(1022), -1.4678102981723264e308, 1.4678102981723264e308)
    };

    for (auto &test : tests) {
        test.validate();
    }
}

template <typename T>
std::pair<T, T> quad_naive(T a, T b, T c) {
    if (a == 0) {
        if (b ==0) {
            return std::pair(NaN<T>, NaN<T>);
        } else {
            return std::pair(-c / b, NaN<T>);
        }
    } else {
        auto discrim = b*b - 4*a*c;
        if (discrim < 0) {
            return std::pair(NaN<T>, NaN<T>);
        } else if (discrim == 0) {
            return std::pair(-b / (2*a), NaN<T>);
        } else {
            auto sqrt_d = sqrt(discrim); 
            auto x1 = (-b + sqrt_d) / (2*a);
            auto x2 = (-b - sqrt_d) / (2*a);
            return std::pair(x1, x2);
        }
    }
}

template <typename T>
std::vector<std::pair<T, T>> benchmark_quadratic(T a, T b, T c, int N) {
    int max_m = 10;

    std::vector<std::pair<T,T>> results(2 * N);
    std::vector<quadratic_test<T>> tests(N);

    for (int i = 0; i < N; ++i) {
        tests[i] = kahan_random_quadratic<T>(i % max_m + 1);
    }

    auto start = std::chrono::steady_clock::now();
    
    for (int i = 0; i < N; ++i) {
        auto test = tests[i];
        results[i] = quad_naive(test.a, test.b, test.c);
    }

    auto end = std::chrono::steady_clock::now();

    auto naive_time = (end - start);

    auto naive_duration = std::chrono::duration <double, std::nano> (naive_time).count() / N;
    std::cout << "Naive implementation: " << naive_duration << " ns" << std::endl;

    start = std::chrono::steady_clock::now();

    for (int i = N; i < 2 * N; ++i) {
        auto test = tests[i - N];
        results[i] = solve(test.a, test.b, test.c);
    }
    
    end = std::chrono::steady_clock::now();

    auto robust_time = (end - start);

    auto robust_duration = std::chrono::duration <double, std::nano> (robust_time).count() / N;
    std::cout << "Robust implementation: " << robust_duration << " ns" << std::endl;

    std:: cout << "Robust / Naive: " << robust_duration / naive_duration << std::endl;

    return results;
}

TEST_CASE("benchmark"){
    auto N = 1'000'000;
    auto results = benchmark_quadratic(1., -1., -1., N);
    std::cout << results[0].first << ", " << results[0].second << std::endl;
    std::cout << results[N].first << ", " << results[N].second << std::endl;
}