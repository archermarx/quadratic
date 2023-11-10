# Quadratic

![Tests](https://github.com/archermarx/quadratic/actions/workflows/test.yml/badge.svg)

Robust C++ quadratic equation solver. Solves quadratics equations without undue overflow, underflow, or catastrophic floating point cancellation.

Based on [*The Ins and Outs of Solving Quadratic Equations with Floating-Point Arithmetic* (2023)](https://www.authorea.com/users/627556/articles/648473-the-ins-and-outs-of-solving-quadratic-equations-with-floating-point-arithmetic) by Frédéric Goualard.

Original Julia implementation at https://github.com/goualard-f/QuadraticEquation.jl.

## Installation

Download the [header file](https://github.com/archermarx/quadratic/include/header.h) and compile with your project using your preferred build system. Then, simply `#include "quadratic.h"`

## Usage

To solve a quadratic equation with parameters `a`, `b`, and `c`, call the `solve` function in the `quadratic` namespace.This function returns a ```std::pair<T, T>```, where `T` is the input type. `T` can be any floating point type.

```cpp
  float a = 1.0f;
  float b = 0.0f;
  float c = -1.0f;
  auto [x1, x2] = quadratic::solve(a, b, c); // [-1, 1]
```

If there are no real solutions, then the pair will contain two NaNs. 

```cpp
  long double a = 1.0;
  long double b = 0.0;
  long double c = 1.0;
  auto [x1, x2] = quadratic::solve(a, b, c); // [nan, nan]
```

If there is only one solution, it will be contained in the first element of the returned pair, while the second element of the pair will be a NaN. 

```cpp
  double a = 1.0;
  double b = 2.0;
  double c = 1.0;
  auto [x1, x2] = quadratic::solve(a, b, c);  // [-1, nan]
```

When the equation has two solutions, the first solution will be the smaller of the two, i.e. `x1 < x2`. 

```cpp
  double a = 1.0;
  double b = -1.0;
  double c = 6.0;
  auto [x1, x2] = quadratic::solve(a, b, c) // [-2, 3]
```

## Performance

In my tests, this is typically about 2-3 times slower than a naive quadratic solver. I think more optimizations are certainly possible, and I welcome performance enhancement suggestions or PRs.
