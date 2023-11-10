# Quadratic

![Tests](https://github.com/archermarx/quadratic/actions/workflows/test.yml/badge.svg)

Robust C++ quadratic equation solver. Solves quadratics equations without undue overflow, underflow, or catastrophic floating point cancellation.

Based on [*The Ins and Outs of Solving Quadratic Equations with Floating-Point Arithmetic* (2023)](https://www.authorea.com/users/627556/articles/648473-the-ins-and-outs-of-solving-quadratic-equations-with-floating-point-arithmetic) by Frédéric Goualard.

Original Julia implementation at https://github.com/goualard-f/QuadraticEquation.jl.

## Installation

Download the [header file](https://github.com/archermarx/quadratic/include/header.h) and compile with your project using your preferred build system. Then, simply `#include "quadratic.h"`

## Usage

To solve a quadratic equation with parameters `a`, `b`, and `c`, call the `solve` function in the `quadratic` namespace. This function returns a ```std::pair<T, T>```, where `T` is the input type.
If there are no real solutions, then the pair will contain two NaNs. If there is only one solution, then it will be contained in the first element, while the second element of the pair will be NaN. 

```cpp
  double a = 1.0;
  double b = 1.0;
  double c = 1.0;
  auto [x1, x2] = quadratic::solve(a, b, c);
```
Robust C++ quadratic equation solver. 
When the equation has two solutions, the first solution will be the smaller of the two, i.e. `x1 < x2`.

## Performance

In my tests, this is typically about 2-3 times slower than a naive quadratic solver. I think more optimizations are certainly possible, and I welcome performance enhancement suggestions or PRs.
