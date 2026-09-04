"""
    plot_shooting(S::Function, S₁::Function, S₂::Function, coord::Symbol)

Dynamically plot the shooting function S and its components S₁ and S₂.

# Arguments
- `S::Function`: The main shooting function taking a 2D vector [p⁰, p]
- `S₁::Function`: Shooting function with the first method of normalization
- `S₂::Function`: Shooting function with the second method of normalization
- `coord::Symbol`: Coordinate system to use, either :S (old) or :T (new)

# Returns
- A combined plot showing the 3D surface of S and 2D plots of S₁ and S₂

# Details
- For :S coordinates: plots S(p⁰, p) with labels ["S", "S₁", "S₂", "p₀", "p⁰"]
- For :T coordinates: plots T(q⁰, q) with labels ["T", "T₁", "T₂", "q₀", "q⁰"]
- The function η(p0) = -√(1 - p0²) defines the boundary curve
"""
function plot_shooting(S::Function, S₁::Function, S₂::Function, coord::Symbol)
    # Set up labels and axis limits based on coordinate system
    if coord == :S
        labels = ["S", "S₁", "S₂", "p₀", "p⁰"]
        xlim_S1 = [-7,2]
        xlim_S2 = [-1,1]
        xlim_S = [0,-1.5]
        ylim_S = [-8, 2]
    elseif coord == :T
        labels = ["T", "T₁", "T₂", "q₀", "q⁰"]
        xlim_S1 = [-3,3]
        xlim_S2 = [-1,1]
        xlim_S = [0,-1.5]
        ylim_S = [-3, 3]
    else
        error("Coordinates must be the old one (:S) or the new one (:T).")
    end

    # Helper functions
    η(p0) = -sqrt.(1 - p0.^2)               # Function η(⋅) = -√(1 - p0²)
    S_(p⁰, p) = S([p⁰, p])                  # Wrapper to call S with unpacked arguments

    # Plot S₁ (first method of normalization)
    plt_S1 = plot(xlabel = labels[5], ylabel = labels[2], xlim = xlim_S1, legend=nothing)
    plot!(plt_S1, S₁, c = 2, lw = 2)

    # Plot S₂ (second method of normalization)
    plt_S2 = plot(xlabel = labels[4], ylabel = labels[3], xlim = [-1,1], legend=nothing)
    plot!(plt_S2, S₂, c = 3, lw = 2)

    # Plot S as a 3D surface with S₁ and S₂ curves overlaid
    plt_S = surface(xlabel = labels[5], ylabel = labels[4], zlabel = labels[1], xflip = true)
    surface!(plt_S, range(xlim_S[1], xlim_S[2], 100), range(ylim_S[1], ylim_S[2], 100), S_, camera = (30,30))
    # Overlay S₁ curve on the surface (at p⁰ = -1)
    plot3d!(plt_S, -1*ones(100), range(ylim_S[1], ylim_S[2], 100), S₁.(range(ylim_S[1], ylim_S[2], 100)), label = labels[2], lw = 2)
    # Overlay S₂ curve on the surface (along the boundary η)
    plot3d!(plt_S, η.(range(-1, 1, 100)), range(-1, 1, 100), S₂.(range(-1, 1, 100)), label = labels[3], lw = 2)

    # Combine all plots into a single figure
    plt_S12 = plot(plt_S1, plt_S2, layout = (1,2))
    plt_total = plot(plt_S, plt_S12, layout = grid(2,1, heights = [2/3, 1/3]), size=(800, 600))

    return plt_total
end

"""
    fit_ellipse(x, y)

Fit an ellipse to a set of 2D points using least squares.

# Arguments
- `x`: Array of x-coordinates
- `y`: Array of y-coordinates

# Returns
- `a`: Semi-major axis length
- `b`: Semi-minor axis length
- `θ`: Rotation angle (in radians)
- `c`: Center coordinates [cx, cy]

# Details
- Fits the general quadratic form: Ax² + Bxy + Cy² + Dx + Ey + F = 0
- Uses algebraic ellipse fitting with constraint F = -1
- Computes geometric parameters (axes, rotation, center) from quadratic coefficients
"""
function fit_ellipse(x, y)
    # Construct design matrix for quadratic form: [x², xy, y², x, y]
    M = hcat(x.^2, x.*y, y.^2, x, y)
    # Solve least squares problem M*p = 1 to get quadratic coefficients
    p = M\ones(length(x))
    A, B, C, D, E = p
    F = -1.0  # Normalization constraint

    # Calculate geometric parameters from quadratic coefficients
    Δ = B^2 - 4*A*C                              # Discriminant (negative for ellipses)
    Λ = (A-C)^2 + B^2                            # Auxiliary parameter for axis calculation
    # Compute semi-minor (b) and semi-major (a) axes
    b, a = [-sqrt(clamp( 2*(A*E^2 + C*D^2 - B*D*E + Δ*F)*
            ( (A+C) + op(sqrt(Λ)) ), 0, Inf)) / Δ   for op in (+, -)]
    θ = atan(-B, C-A)/2                          # Rotation angle
    c = [(2*C*D - B*E)/Δ, (2*A*E - B*D)/Δ]       # Center coordinates [cx, cy]

    return a, b, -θ+Base.π/2, c
end

"""
    plot_sol(sol, size = (800, 600))

Plot the solution of an optimal control problem, showing state, costate, and control variables.

# Arguments
- `sol`: Solution object from DifferentialEquations.jl (ODE solution)
- `size`: Tuple specifying figure size (width, height), default (800, 600)

# Returns
- A combined plot with state and costate variables (top) and control (bottom)

# Details
- Extracts state variables x⁰, x and costate variables p⁰, p from solution
- Control u is computed as sign(p) (bang-bang control)
- Layout: 2x2 grid for state/costate, with control plot below
"""
function plot_sol(sol, size = (800, 600))

    # Extract time and state/costate variables from solution
    t = sol.t                                     # Time points
    x⁰ = [sol.u[i][1] for i in 1:length(sol.u)]   # First state variable
    x  = [sol.u[i][2] for i in 1:length(sol.u)]   # Second state variable
    p⁰ = [sol.u[i][3] for i in 1:length(sol.u)]   # First costate variable
    p  = [sol.u[i][4] for i in 1:length(sol.u)]   # Second costate variable
    u = sign.(p)                                  # Control (bang-bang: ±1)

    # Create individual plots for each variable
    plt_x⁰ = plot(t, x⁰, ylabel = "x⁰", label = nothing, lw = 2, title = "state", titlefontsize = 10)
    plt_x  = plot(t, x , ylabel = "x", label = nothing, lw = 2)
    plt_p⁰ = plot(t, p⁰, ylim=[-1,1], lw = 2, label = nothing, title = "costate", titlefontsize = 10)
    plt_p  = plot(t, p, label = nothing, lw = 2)
    plt_u  = plot(t, u , xlabel = "u", title = "control", titlefontsize = 10, label = nothing, lw = 2)

    # Combine state and costate plots in 2x2 grid
    plt_xp = plot(plt_x⁰, plt_p⁰, plt_x, plt_p, layout=(2, 2))
    # Combine with control plot (2/3 height for state/costate, 1/3 for control)
    return plot(plt_xp, plt_u, layout = grid(2,1, heights = [2/3, 1/3]), size=size)
end