# Optimal control problem

We consider the following optimal control problem 

```math
    \left\{ \begin{array}{ll}
    \displaystyle \min_{x,u} \int_{t_0}^{t_f} x(t) ~\mathrm dt \\[1em]
    \text{s.c.}~\dot x(t) = u(t), & t\in [t_0, t_f]~\mathrm{a.e.}, \\[0.5em]
    \phantom{\mathrm{s.c.}~} u(t) \in [-1,1], & t\in [t_0, t_f], \\[0.5em]
    \phantom{\mathrm{s.c.}~} x(t_0) = x_0, \quad x(t_f) = x_f,
    \end{array} \right.
```

with $x_0$, $t_0$, $x_f$ and $t_f$ fixed. This problem is simple, and can be analytically solve without the use of numerical method. However, the goal is to solve this problem by indirect shooting.  

# Dependencies

All the numerical simulations to generate this documentation from `MRI.jl` are performed with 
the following packages.

```@example
using Pkg
Pkg.status()
```
