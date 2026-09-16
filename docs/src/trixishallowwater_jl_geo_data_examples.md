# Cliffs of Moher

The [Rhine river](https://trixi-framework.github.io/TrixiBottomTopography.jl/stable/trixishallowwater_jl_examples/)
examples use bottom topography that is already available in the TrixiBottomTopography.jl
format. This section shows two examples that instead use topography obtained for an
arbitrary region of the world with
[GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl),
as described in [Real topography data](@ref).

The region is the coastline at the Cliffs of Moher in Ireland. The domain covers the open
sea in the west and the cliffs, which rise up to roughly 197 m above sea level, in the east.
In both examples a wave travels east and runs up the cliff face, which requires the
wetting and drying capabilities of
[TrixiShallowWater.jl](https://github.com/trixi-framework/TrixiShallowWater.jl).

## One dimensional wave run-up

The underlying example file can be found [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/trixishallowwater_cliffs_1D.jl).

First, all the necessary packages must be included at the beginning of the file.

```@example geo_trixi_1D
# Include packages
using TrixiBottomTopography
using CairoMakie
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
```

In contrast to the Rhine examples,
[OrdinaryDiffEqSSPRK.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/explicit/SSPRK/)
is used here because the strong stability preserving methods it provides accept the
positivity limiter that keeps the water height non-negative in dry regions.

The one dimensional cut through the topography is shipped with the repository, so it can be
loaded directly.

```@example geo_trixi_1D
root_dir = pkgdir(TrixiBottomTopography)
cliffs_data = joinpath(root_dir, "examples", "data", "cliffs_data_1d_10_x.txt")
nothing #hide
```

The data is used to define the B-spline interpolation function as described in
[B-spline structure](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/structure/)
and [B-spline function](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/function/).
No smoothing is applied here because it would flatten the steep cliff face by several meters.

```@example geo_trixi_1D
const spline_struct = CubicBSpline(cliffs_data; end_condition = "not-a-knot")
spline_func(x::Float64) = spline_interpolation(spline_struct, x)
```

Plotting the interpolated topography shows the deep water in the west, the flat shelf that
the SRTM data reports as zero, and the cliff face in the east.

```@example geo_trixi_1D
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], 500))
plot_topography(x_int_pts, spline_func.(x_int_pts); xlabel = "x [m]", ylabel = "z [m]")
```

The topography is given with respect to sea level. A positive background total water height
$H_0$ therefore floods the shallow shelf in front of the cliffs while the cliff face itself
stays dry.

```@example geo_trixi_1D
equations = ShallowWaterEquations1D(gravity = 9.81, H0 = 10.0)
```

At time $t=0$ the water surface west of $x = -350$ is raised to $20.0$ while the rest of the
domain stays at the background water height $10.0$. This step collapses and sends a wave
towards the cliffs. Because part of the domain is dry, the water surface has to be shifted
by the `threshold_limiter` of the equations to keep the water height `h` strictly positive.

```@example geo_trixi_1D
# Defining initial condition of a wave which travels towards the cliffs
function initial_condition_wave(x, t, equations::ShallowWaterEquations1D)
    H = x[1] < -350.0 ? 20.0 : equations.H0
    v = 0.0
    b = spline_func(x[1])

    H = max(H, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall
boundary_condition = boundary_condition_slip_wall
nothing #hide
```

The upcoming code parts will **not** be covered in full detail. For more information, see
the documentation of [Trixi.jl](https://trixi-framework.github.io/TrixiDocumentation/stable/)
and [TrixiShallowWater.jl](https://trixi-framework.github.io/TrixiShallowWater.jl/stable/).
The essential difference to the Rhine examples is the discretization: wetting and drying
requires the hydrostatic reconstruction of Chen and Noelle together with a shock capturing
volume integral.

```@example geo_trixi_1D
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                    alpha_max = 0.5,
                                                    alpha_min = 0.001,
                                                    alpha_smooth = true,
                                                    variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)
nothing #hide
```

The mesh spans exactly the interval covered by the topography data.

```@example geo_trixi_1D
coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 6,
                periodicity = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)
nothing #hide
```

The positivity limiter is handed to the time integration method as a stage limiter. It cuts
off water heights below the `threshold_limiter` after every Runge-Kutta stage.

```@example geo_trixi_1D
tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 90)

sol = solve(ode, SSPRK43(; stage_limiter!), abstol = 1.0e-6, reltol = 1.0e-6,
            saveat = visnodes)
nothing #hide
```

Finally, the solution is animated. The water surface reaches the cliff face after roughly
50 s, runs up to about 21 m, and is then reflected back towards the open sea.

```@example geo_trixi_1D
j = Observable(1)
time = Observable(0.0)

pd_list = [PlotData1D(sol.u[i], semi) for i in 1:length(sol.t)]
f = Figure()
title_string = lift(t -> "time t = $(round(t, digits=3))", time)
ax = Axis(f[1, 1], xlabel = "x [m]", ylabel = "z [m]", title = title_string)

height = lift(i -> pd_list[i].data[:, 1], j)
bottom = lift(i -> pd_list[i].data[:, 3], j)
lines!(ax, pd_list[1].x, height)
lines!(ax, pd_list[1].x, bottom)
ylims!(ax, -40, 60)

record(f, "animation_cliffs_1d.gif", 1:length(pd_list)) do tt
    j[] = tt
    time[] = sol.t[tt]
end
nothing #hide
```

![simCliffs1D](animation_cliffs_1d.gif)

## Two dimensional wave run-up

The underlying example file can be found [here](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/trixishallowwater_cliffs_2D.jl).

The two dimensional example uses the same scenario on the full topography. Since the
solution is post-processed with
[Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl) instead of Makie.jl, that
package is loaded as well.

```@example geo_trixi_2D
# Include packages
using TrixiBottomTopography
using CairoMakie
using OrdinaryDiffEqSSPRK
using Trixi
using Trixi2Vtk
using TrixiShallowWater
```

The two dimensional data is interpolated with a bicubic B-spline.

```@example geo_trixi_2D
root_dir = pkgdir(TrixiBottomTopography)
cliffs_data = joinpath(root_dir, "examples", "data", "cliffs_data_2d_10.txt")

const spline_struct = BicubicBSpline(cliffs_data; end_condition = "not-a-knot")
spline_func(x::Float64, y::Float64) = spline_interpolation(spline_struct, x, y)
```

Sampling the interpolation function on a finer set of nodes gives a three dimensional view
of the coastline. Note that `Makie.surface` expects the values as `z[x_index, y_index]`
whereas `evaluate_two_dimensional_interpolant` returns them as `z[y_index, x_index]`.

```@example geo_trixi_2D
n = 200
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

z_int_pts = evaluate_two_dimensional_interpolant(spline_func, x_int_pts, y_int_pts)

plot_topography(x_int_pts, y_int_pts, permutedims(z_int_pts);
                xlabel = "x\n [m]", ylabel = "y\n [m]", zlabel = "z\n [m]",
                azimuth_angle = -120 * pi / 180, elevation_angle = 20 * pi / 180)
```

The equations and the initial condition are the direct two dimensional analogue of the one
dimensional case above.

```@example geo_trixi_2D
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 10.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    H = x[1] < -350.0 ? 20.0 : equations.H0
    v1 = 0.0
    v2 = 0.0
    b = spline_func(x[1], x[2])

    H = max(H, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_wave
boundary_condition = boundary_condition_slip_wall
nothing #hide
```

The approximation space is set up in the same way as in the one dimensional case.

```@example geo_trixi_2D
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                    alpha_max = 0.5,
                                                    alpha_min = 0.001,
                                                    alpha_smooth = true,
                                                    variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)
nothing #hide
```

Here a [`P4estMesh`](https://trixi-framework.github.io/TrixiDocumentation/stable/meshes/p4est_mesh/)
is used instead of a `TreeMesh` because it supports the non-conforming adaptive mesh
refinement that resolves the moving wave front.

```@example geo_trixi_2D
coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = P4estMesh((1, 1);
                 polydeg = 1,
                 coordinates_min = coordinates_min,
                 coordinates_max = coordinates_max,
                 initial_refinement_level = 4,
                 periodicity = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

tspan = (0.0, 60.0)
ode = semidiscretize(semi, tspan)
nothing #hide
```

The solution is written to an output directory from which it is later converted to VTK
files. The positivity limiter is passed both to the AMR callback, so that it is applied
after every mesh adaptation, and to the time integration method.

```@example geo_trixi_2D
output_dir = mktempdir()

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

amr_indicator = IndicatorLöhner(semi, variable = first)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 4,
                                      med_level = 5, med_threshold = 0.1,
                                      max_level = 6, max_threshold = 0.5)
amr_callback = AMRCallback(semi, amr_controller,
                           interval = 1,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true,
                           limiter! = stage_limiter!)

stepsize_callback = StepsizeCallback(cfl = 0.2)

save_solution = SaveSolutionCallback(dt = 0.5,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = output_dir,
                                     solution_variables = cons2prim)

callbacks = CallbackSet(amr_callback, stepsize_callback, save_solution)

sol = solve(ode, SSPRK43(; stage_limiter!); dt = 1.0, adaptive = false,
            callback = callbacks)
nothing #hide
```

Finally, the Trixi.jl output files are post-processed with Trixi2Vtk.jl.

```@example geo_trixi_2D
trixi2vtk(joinpath(output_dir, "solution_*.h5"), output_directory = output_dir)
nothing #hide
```

It is possible to open the created `solution_00000.pvd` file with
[ParaView](https://www.paraview.org/) and create a video of the simulation. In ParaView,
one can apply two instances of the *Warp By Scalar* filter to visualize the water height
and the bathymetry in three dimensions. Many additional customizations, e.g., color
scaling and fonts, are available there.
