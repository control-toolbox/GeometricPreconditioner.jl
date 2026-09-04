The goal of this application is to improve the convergence of indirect shooting method for solving optimal control problems. The indirect shooting method is a powerful technique for solving optimal control problems via the Pontryagin Maximum Principle, but it often suffers from poor convergence properties due to the sensitivity of the shooting function.

The first part focuses on the proper implementation of the indirect shooting method. A critical aspect of the shooting method is the computation of the Jacobian of the shooting function, which is required by Newton-type solvers. Automatic differentiation can be used to compute this Jacobian, but some implementation choices can lead to a false result.

The second part introduces a geometric preconditioning approach to improve the convergence of the indirect shooting method.  The core idea is to transform the boundary value problem into a coordinate system that aligns with the natural structure of the problem. The comparison between the standard and preconditioned methods will be presented, and will highlight the benefit of using the geometric preconditioned shooting function.

## Optimal control problem

We consider in this application the following toy optimal control problem 

```math
    \left\{ \begin{array}{ll}
    \displaystyle \min_{x,u} \int_{t_0}^{t_f} x(t) ~\mathrm dt \\[1em]
    \text{s.c.}~\dot x(t) = u(t), & t\in [t_0, t_f]~\mathrm{a.e.}, \\[0.5em]
    \phantom{\mathrm{s.c.}~} u(t) \in [-1,1], & t\in [t_0, t_f], \\[0.5em]
    \phantom{\mathrm{s.c.}~} x(t_0) = x_0, \quad x(t_f) = x_f,
    \end{array} \right.
```

with ``x_0``, ``t_0``, ``x_f`` and ``t_f`` fixed. This problem is simple, and can be analytically solved without the use of numerical methods. However, the goal is to solve this problem by indirect shooting.  

## References

- [insert article]
- Olivier Cots, Rémy Dutto, Sophie Jan, Serge Laporte (2024). [Geometric preconditioner for indirect shooting and application to hybrid vehicle](https://www.sciencedirect.com/science/article/pii/S2405896324018950). _4th IFAC MICNON Conference_.
- Rémy Dutto (2024). [Méthode à deux niveaux et préconditionnement géométrique en contrôle optimal. Application au problème de répartition de couple des véhicules hybrides électriques](https://hal.science/tel-04792906v1). _Thèse de doctorat, Université de Toulouse_.

## Reproducibility

```@setup main
using Pkg
using InteractiveUtils
using Markdown

# Download links for the benchmark environment
function _downloads_toml(DIR)
    link_manifest = joinpath("assets", DIR, "Manifest.toml")
    link_project = joinpath("assets", DIR, "Project.toml")
    return Markdown.parse("""
    You can download the exact environment used to build this documentation:
    - 📦 [Project.toml]($link_project) - Package dependencies
    - 📋 [Manifest.toml]($link_manifest) - Complete dependency tree with versions
    """)
end
```

```@example main
_downloads_toml(".") # hide
```

```@raw html
<details style="margin-bottom: 0.5em; margin-top: 1em;"><summary>ℹ️ Version info</summary>
```

```@example main
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📦 Package status</summary>
```

```@example main
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📚 Complete manifest</summary>
```

```@example main
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```
