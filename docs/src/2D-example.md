### Optimal control problem

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

### Indirect method

We introduce the pseudo-Hamiltonian 

```math
    h(x,p,p^0,u) = p^0 x + p u.
```

For the sake of simplicity, we consider in this notebook only the normal case, and we fix $p^0 = -1$. According to the Pontryagin maximum principle, the maximizing control is given by $u(x,p) \to \mathrm{sign}(p)$. This function is non-differentiable, and may lead to numerical issues.  

Let us import the necessary package and define the studied optimal control problem with some fixed initial and final time and state values.

```@example main
using OptimalControl
using Plots
using ForwardDiff
using DifferentialEquations
using MINPACK

t0 = 0
x0 = 0
tf = 5
xf = 0                      # initial and final time and state

@def ocp begin              # problem definition

    t ∈ [ t0, tf ], time
    x ∈ R, state
    u ∈ R, control

    x(t0) == x0
    x(tf) == xf

    ẋ(t) == u(t)      

    ∫( x(t) ) → min

end
nothing # hide
```