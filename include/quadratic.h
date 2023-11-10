#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <ostream>
#include <iostream>
#include <algorithm>

namespace quadratic
{

//====================================================
//  Floating point constants and utility functions
//====================================================
template <typename T>
constexpr T NaN = std::numeric_limits<T>::quiet_NaN();

template <typename T>
constexpr T inf = std::numeric_limits<T>::infinity();

// Maximum exponent
template <typename T>
constexpr int E_MIN = std::numeric_limits<T>::min_exponent - 1;

// Minimum exponent
template <typename T>
constexpr int E_MAX = std::numeric_limits<T>::max_exponent - 1;

// Number of bits in mantissa
template <typename T>
constexpr int NUM_DIGITS = std::numeric_limits<T>::digits;

// max val of Ea + Ea - 2 Eb to avoid overflow
template<typename T>
constexpr int ECP_MAX = E_MAX<T> - 2 - NUM_DIGITS<T> / 2;

// min val of Ea + Ea - 2 Eb to avoid underflow
template<typename T>
constexpr int ECP_MIN = E_MIN<T> + 2 * NUM_DIGITS<T> - 4;

// Utility functions
template <typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
inline bool samesign(T a, T b) { return a * b > 0; }

// Return a pair (Kout, K2) of exponents such that neither is outside [E_MIN, E_MAX]
// with Kout = K2 - Kin
// float32:  E_MIN = -126,   E_MAX = 127
// float64:  E_MIN = -1022,  E_MAX = 1023
// float128: E_MIN = -16382, E_MAX = 16383
template <typename T>
inline int keep_exponent_in_check(int Kin) { return std::clamp(Kin, E_MIN<T>, E_MAX<T>); }

// Accurate computation of the discriminant with Kahan's method, with an fma instruction
template <typename T>
inline T kahan_discriminant_fma(T a, T b, T c) {
    T q = 4*a*c;
    T dq = std::fma(-c, 4*a, q);
    T dp = std::fma(b, b, -q);
    return dq + dp;
}

// Algorithm based on "The Ins and Outs of Solving Quadratic Equations with Floating-Point Arithmetic (Goualard 2023)"
template <typename T>
std::pair<T, T> solve(T a, T b, T c) {

    // Function only accepts floats
    static_assert(
        std::is_floating_point_v<T>, 
        "Parameters a, b, and c must be float, double, long double, or other floating point type."
    );

    constexpr std::pair<T,T> NO_SOLUTIONS = std::pair(NaN<T>, NaN<T>);
    constexpr T zero = T(0);
    constexpr T one = T(1);

    // Invalid case: any of the parameters are nan or inf
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
        return NO_SOLUTIONS;
    }

    // Degenerate cases (one of the coefficients is zero)
    if (a == 0) {
        if (b == 0) {
            // if a and b are zero, then the quadratic is invalid
            // if c == 0, then all real numbers are solutions, i.e. 0 == 0,
            // and if c /= 0, then no real numbers are solutions, i.e. c == 0
            return NO_SOLUTIONS;
        } else {
            if (c == 0) {
                // If a and c are zero, but b is nonzero, then we have bx = 0, so only
                // one real solution (i.e. x = 0) exists
                return std::pair(zero, NaN<T>);
            } else {
                // If a == 0 but b and c are nonzero, then we have bx + c = 0, which is
                // a linear equation with solution x = -c / b
                return std::pair(-c / b, NaN<T>);
            }
        }
    } else {
        if (b == 0) {
            if (c == 0) {
                // If a /= 0, but b and c == 0, then we have ax^2 = 0, implying x = 0;
                return std::pair(zero, NaN<T>);
            } else {
                // If a and c /= 0, but b == 0. then we have x^2 = -c/a, which has either no
                // real solutions if sign(a) = sign(c), or two real solutions otherwise
                if (samesign(a, c)) {
                    return NO_SOLUTIONS;
                } else {
                    // Split a and c into significant and exponent parts
                    int exp_a, exp_c;
                    T signif_a = frexp(a, &exp_a);
                    T signif_c = frexp(c, &exp_c);

                    int ecp = exp_c - exp_a;

                    int dM = ecp & ~1;  // dM = floor(ecp/2) * 2
                    int M = dM >> 1;    // M = dM / 2
                    int E = ecp & 1;    // E = odd(ecp) ? 1 : 0
                    T S = sqrt(-ldexp(signif_c, E) / signif_a);

                    int M1 = keep_exponent_in_check<T>(M);
                    int M2 = M - M1;
                    T x = ldexp(ldexp(S, M1), M2);
                    return std::pair(-x, x);
                }   
            }
        } else {
            if (c == 0) {
                return std::minmax(zero, -b / a);
            } else {
                // a, b, and c all nonzero

                // Split a, b, and c into significand and exponent
                int exp_a, exp_b, exp_c;
                T signif_a = frexp(a, &exp_a);
                T signif_b = frexp(b, &exp_b);
                T signif_c = frexp(c, &exp_c);
                
                int K = exp_b - exp_a;
                // ecp = exp_c + exp_a - 2 exp_b
                int ecp = exp_c + exp_a - (exp_b << 1);

                if (ecp >= ECP_MIN<T> && ecp < ECP_MAX<T>) {
                    T c2 = ldexp(signif_c, ecp);
                    T delta = kahan_discriminant_fma(signif_a, signif_b, c2);
                    if (delta < 0) {
                        return NO_SOLUTIONS;
                    }
                    int K1 = keep_exponent_in_check<T>(K);
                    int K2 = K - K1;
                    T expval = ldexp(ldexp(one, K1), K2);

                    if (delta > 0) {
                        T B = signif_b + copysign(sqrt(delta), b);
                        T y1 = -(2 * c2) / B;
                        T y2 = -B / (2 * signif_a);
                        T x1 = y1 * expval;
                        T x2 = y2 * expval;
                        return std::minmax(x1, x2);
                    }
                    // delta == 0
                    T x1 = -(signif_b / (2 * signif_a)) * expval;
                    return std::pair(x1, NaN<T>);
                }

                int dM = ecp & ~1;              // dM = floor(ecp/2) * 2
                int M = dM >> 1;                // M = dM / 2
                int E = ecp & 1;                // E = odd(ecp) ? 1 : 0
                T c3 = ldexp(signif_c, E);   // c3 = signif_c * 2^E
                T S = sqrt(abs(c3 / signif_a));

                if (ecp < ECP_MIN<T>) {
                    T y1 = -signif_b / signif_a;
                    T y2 = c3 / (signif_a * y1);
                    int dMK = dM + K;
                    int dMK1 = keep_exponent_in_check<T>(dMK);
                    int dMK2 = dMK - dMK1;
                    int K1 = keep_exponent_in_check<T>(K);
                    int K2 = K - K1;

                    T x1 = ldexp(ldexp(y1,   K1),   K2); // x1 = y1 * 2^K1 * 2^K2
                    T x2 = ldexp(ldexp(y2, dMK1), dMK2); // x2 = y2 * 2^dMK1 * 2^dMK2
                    return std::minmax(x1, x2);
                }

                // ecp >= ECP_MAX<T>
                if (samesign(a, c)) {
                    return NO_SOLUTIONS;
                } else {
                    int MK = M + K;
                    int MK1 = keep_exponent_in_check<T>(MK);
                    int MK2 = MK - MK1;
                    T x1 = ldexp(ldexp(S, MK1), MK2);
                    return std::pair(-x1, x1);
                }
            }                    
        }
    }
}

}

#endif