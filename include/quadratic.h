#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <ostream>
#include <iostream>

using std::cout;
using std::endl;

// Constants
template <typename T>
const T NaN = std::numeric_limits<T>::quiet_NaN();

template <typename T>
const T inf = std::numeric_limits<T>::infinity();

template <typename T>
const T eps = std::numeric_limits<T>::epsilon();

// Utility functions
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
bool isapprox(T x, T y) {
    return  (std::isnan(x) && std::isnan(y)) ||
            (std::isinf(x) && std::isinf(y) && sgn(x) == sgn(y)) ||
            std::abs(x - y) <= std::sqrt(eps<T>) * std::max(abs(x), abs(y));
}

// Floating point constants

// Maximum exponent
template <typename T>
const int E_MIN = std::numeric_limits<T>::min_exponent - 1;

// Minimum exponent
template <typename T>
const int E_MAX = std::numeric_limits<T>::max_exponent - 1;

// Number of bits in mantissa
template <typename T>
const int NUM_DIGITS = std::numeric_limits<T>::digits;

// max val of Ea + Ea - 2 Eb to avoid overflow
template<typename T>
const int ECP_MAX = E_MAX<T> - 2 - NUM_DIGITS<T> / 2;

// min val of Ea + Ea - 2 Eb to avoid underflow
template<typename T>
const int ECP_MIN = E_MIN<T> + 2 * NUM_DIGITS<T> - 4;


enum quadratic_return_code: char{
    QUAD_NOT_RUN,
    QUAD_NO_SOLUTIONS,
    QUAD_ONE_SOLUTION,
    QUAD_TWO_SOLUTIONS
};

// Printing
inline std::ostream& operator<<(std::ostream &out, const quadratic_return_code q) {
    switch(q) {
        case QUAD_NOT_RUN:
            return out << "QUAD_NOT_RUN";
        case QUAD_NO_SOLUTIONS:
            return out << "QUAD_NO_SOLUTIONS";
        case QUAD_ONE_SOLUTION:
            return out << "QUAD_ONE_SOLUTION";
        case QUAD_TWO_SOLUTIONS:
            return out << "QUAD_TWO_SOLUTIONS";
        default:
            return out << "QUAD_INVALID_RETCODE";
    }
}

// Return a pair (Kout, K2) of exponents such that neither is outside [E_MIN, E_MAX]
// with Kout = K2 = Kin
// float32:  E_MIN = -126,   E_MAX = 127
// float64:  E_MIN = -1022,  E_MAX = 1023
// float128: E_MIN = -16382, E_MAX = 16383
template <typename T>
std::tuple<int, int> keep_exponent_in_check(int Kin) {
    int K2, Kout;
    
    if (Kin <= E_MAX<T> && Kin >= E_MIN<T>) {
        Kout = Kin; K2 = 0;
    } else if (Kin < E_MIN<T>) {
        Kout = E_MIN<T>; K2 = Kin - E_MIN<T>;
    } else {
        Kout = E_MAX<T>; K2 = Kin - E_MAX<T>;
    }
    return std::make_tuple(Kout, K2);
}

// Accurate computation of the discriminant with Kahan's method, with an fma instruction
template <typename T>
T kahan_discriminant_fma(T a, T b, T c) {
    auto p = b*b;
    auto q = 4*a*c;
    auto d = p - q;
    if (3 * std::abs(d) >= p + 4*q) {
        // If b^2 and 4ac are different enough
        return d;
    }
    auto dp = std::fma(b,b,-p);
    auto dq = std::fma(4*a,c,-q);
    d = (p - q) + (dp - dq);
    return d;
}

// Algorithm based on "The Ins and Outs of Solving Quadratic Equations with Floating-Point Arithmetic (Goualard 2023)"
template <typename T>
class quadratic  {
    public:
        T a, b, c;
        T x1, x2;

        quadratic() : x1(NaN<T>), x2(NaN<T>) {}
        quadratic(T _a, T _b, T _c) : a(_a), b(_b), c(_c), x1(NaN<T>), x2(NaN<T>) {}

        quadratic_return_code solve() {
            // Degenerate case 1: any of the parameters are nan or inf
            if (std::isnan(a) || std::isnan(b) || std::isnan(c) || 
                std::isinf(a) || std::isinf(b) || std::isinf(c)) {
                return QUAD_NO_SOLUTIONS;
            }

            // More degenerate cases
            if (a == 0) {
                if (b == 0) {
                    // if a and b are zero, then the quadratic is invalid
                    // if c == 0, then all real numbers are solutions, i.e. 0 == 0,
                    // and if c /= 0, then no real numbers are solutions, i.e. c == 0
                    return QUAD_NO_SOLUTIONS;
                } else {
                    if (c == 0) {
                        // If a and c are zero, but b is nonzero, then we have bx = 0, so only
                        // one real solution (i.e. x = 0) exists
                        x1 = 0;
                        return QUAD_ONE_SOLUTION;
                    } else {
                        // If a == 0 but b and c are nonzero, then we have bx + c = 0, which is
                        // a linear equation with solution x = -c / b
                        x1 = -c / b;
                        return QUAD_ONE_SOLUTION;
                    }
                }
            } else {
                if (b == 0) {
                    if (c == 0) {
                        // If a /= 0, but b and c == 0, then we have ax^2 = 0, implying x = 0;
                        x1 = 0;
                        return QUAD_ONE_SOLUTION;
                    } else {
                        // If a and c /= 0, but b == 0. then we have x^2 = -c/a, which has either no
                        // real solutions if sign(a) = sign(c), or two real solutions otherwise
                        if (sgn(a) == sgn(c)) {
                            return QUAD_NO_SOLUTIONS;
                        } else {
                            // Split a and c into significant and exponent parts
                            int exp_a, exp_c;
                            auto signif_a = frexp(a, &exp_a);
                            auto signif_c = frexp(c, &exp_c);

                            int ecp = exp_c - exp_a;

                            int dM = ecp & ~1;  // dM = floor(ecp/2) * 2
                            int M = dM >> 1;    // M = dM / 2
                            int E = ecp & 1;    // E = odd(ecp) ? 1 : 0
                            auto S = sqrt(-signif_c * exp2(E) / signif_a);

                            std::tuple<int, int> Ms = keep_exponent_in_check<T>(M);
                            int M1 = std::get<0>(Ms);
                            int M2 = std::get<1>(Ms);
                            T x = S * exp2(M1) * exp2(M2);
                            x1 = -x;
                            x2 = x;
                            return QUAD_TWO_SOLUTIONS;
                        }   
                    }
                } else {
                    if (c == 0) {
                        // a /= 0, b /= 0, c /= 0
                        if (sgn(a) == sgn(b)) {
                            x1 = -b / a; x2 = 0.0;
                        } else {
                            x1 = 0.0; x2 = -b / a;
                        }
                        return QUAD_TWO_SOLUTIONS;
                    } else {
                        // a, b, and c all nonzero

                        // Split a, b, and c into significand and exponent
                        int exp_a, exp_b, exp_c;
                        auto signif_a = frexp(a, &exp_a);
                        auto signif_b = frexp(b, &exp_b);
                        auto signif_c = frexp(c, &exp_c);
                        
                        auto K = exp_b - exp_a;
                        auto ecp = exp_c + exp_a - 2 * exp_b;

                        if (ecp >= ECP_MIN<T> && ecp < ECP_MAX<T>) {
                            auto c2 = signif_c * exp2(ecp);
                            auto delta = kahan_discriminant_fma(signif_a, signif_b, c2);
                            if (delta < 0) {
                                return QUAD_NO_SOLUTIONS;
                            }
                            auto Ks = keep_exponent_in_check<T>(K);
                            auto K1 = std::get<0>(Ks);
                            auto K2 = std::get<1>(Ks);
                            if (delta > 0) {
                                auto B = signif_b + sgn(b) * sqrt(delta);
                                auto y1 = -(2 * c2) / B;
                                auto y2 = -B / (2 * signif_a);
                                auto _x1 = y1 * exp2(K1) * exp2(K2);
                                auto _x2 = y2 * exp2(K1) * exp2(K2);
                                x1 = std::min(_x1, _x2);
                                x2 = std::max(_x1, _x2);
                                return QUAD_TWO_SOLUTIONS;
                            }
                            // delta == 0
                            x1 = -(signif_b / (2 * signif_a)) * exp2(K1) * exp2(K2);
                            return QUAD_ONE_SOLUTION;
                        }

                        int dM = ecp & ~1;  // dM = floor(ecp/2) * 2
                        int M = dM >> 1;    // M = dM / 2
                        int E = ecp & 1;    // E = odd(ecp) ? 1 : 0
                        auto c3 = signif_c * exp2(E);
                        auto S = sqrt(abs(c3 / signif_a));

                        if (ecp < ECP_MIN<T>) {
                            auto y1 = -signif_b / signif_a;
                            auto y2 = c3 / (signif_a * y1);
                            auto dMK = keep_exponent_in_check<T>(dM + K);
                            auto Ks = keep_exponent_in_check<T>(K);
                            int K1 = std::get<0>(Ks);
                            int K2 = std::get<1>(Ks);
                            int dMK1 = std::get<0>(dMK);
                            int dMK2 = std::get<1>(dMK);

                            auto _x1 = y1 * exp2(K1) * exp2(K2);
                            auto _x2 = y2 * exp2(dMK1) * exp2(dMK2);
                            x1 = std::min(_x1, _x2);
                            x2 = std::max(_x1, _x2);
                            return QUAD_TWO_SOLUTIONS;
                        }

                        // ecp >= ECP_MAX<T>
                        if (sgn(a) == sgn(c)) {
                            return QUAD_NO_SOLUTIONS;
                        } else {
                            auto MK = keep_exponent_in_check<T>(M + K);
                            int MK1 = std::get<0>(MK);
                            int MK2 = std::get<1>(MK);
                            auto _x1 = S * exp2(MK1) * exp2(MK2);
                            auto _x2 = -_x1;
                            x1 = _x2;
                            x2 = _x1;
                            return QUAD_TWO_SOLUTIONS;
                        }

                    }                    
                }
            }
            return QUAD_NOT_RUN;
        }
};

#endif