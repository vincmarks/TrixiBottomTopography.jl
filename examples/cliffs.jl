###########################################################################
# In this example file we will use topography around the cliffs of Moher  #
# to simulate a wave breaking against the cliffs. The topography data is  #
# generated using TrixiBottomTopography and GeophysicalModelGenerator.    #                   
###########################################################################

#load the  packages
using GeophysicalModelGenerator
using GMT
using TrixiBottomTopography
using OrdinaryDiffEqSSPRK
using Trixi2Vtk
using Trixi
using TrixiShallowWater
using CairoMakie

# Here we load the topography data at the coordinates of the cliffs of Moher
Topo, p, Topo_Cart = geo_topo_impression(resolution = "@earth_relief_01s",
                                         lon_min = -9.425224,
                                         lon_max = -9.441139,
                                         lat_min = 52.972167,
                                         lat_max = 52.980888)

# Create a .xyz file of the topogrphy data. Here low_x, high_x, low_y, high_y are choosen
# based on the output of Topo_Cart

df_xyz, Topo_Cart_orth = create_topography_data(low_x = -0.55,
                                                high_x = 0.55,
                                                gridsize_x = 0.001,
                                                low_y = -0.55,
                                                high_y = 0.55,
                                                gridsize_y = 0.001,
                                                write_path = joinpath(@__DIR__, "data"),
                                                dataname = "cliffs.xyz",
                                                Topo = Topo,
                                                p = p)
###############################################################################
# Define file paths

data_dir = joinpath(@__DIR__, "data")
path_src_file = joinpath(data_dir, "cliffs.xyz")

cliffs_2d = joinpath(data_dir, "cliffs.txt")

# Convert data
nx = size(Topo_Cart_orth.x.val, 1)
ny = size(Topo_Cart_orth.y.val, 2)

convert_geo_2d(path_src_file, cliffs_2d, nx = nx, ny = ny; excerpt = 10)

###############################################################################
# visualization of the topography

data_dir = joinpath(@__DIR__, "data")
data = joinpath(data_dir, "cliffs.txt")

# B-spline interpolation of the underlying data
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Define interpolation points
n = 200
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

# Get interpolated matrix
z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)
# Makie expects z[x_index, y_index]; the spline sampler returns z[y_index, x_index].
z_int_pts = permutedims(z_int_pts)
# Plot the topography

#visualize the topography 
plot_topography(x_int_pts,
                y_int_pts,
                z_int_pts;
                xlabel = "x\n [m]",
                ylabel = "y\n [m]",
                zlabel = "z\n [m]",
                azimuth_angle = -120 * pi / 180,
                elevation_angle = 20 * pi / 180,)

####################################################################
# for the simulation we use the shallow water equations in 2D

equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 40.0,)

bathymetry(x::Float64, y::Float64) = spline_interpolation(spline_struct, x, y)

####################################################################
# intial condition

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(0.0, 0.0)
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Water surface elevation above the same datum as the bottom topography.
    H = equations.H0 + (r < 50.0 ? 10.0 : 0.0)
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = spline_func(x1, x2)

    # Keep a thin water film on dry land instead of creating negative depths.
    h = max(equations.threshold_limiter, H - b)
    return SVector(h, h * v1, h * v2, b)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space
basis = LobattoLegendreBasis(7)
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)
indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)
solver = DGSEM(basis,
               surface_flux,
               volume_integral)

###############################################################################
# Get the TreeMesh with wall boundaries

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = false)

mesh = P4estMesh((1, 1);
                   polydeg = 1,
                   coordinates_min = coordinates_min,
                   coordinates_max = coordinates_max,
                   initial_refinement_level = 4,
                   periodicity = false,)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan)


######################

output_dir = "out"

if isdir(output_dir)
    rm(output_dir, recursive = true)  
end
mkpath(output_dir)

#adaptive mesh refinement

amr_indicator = IndicatorLöhner(semi, variable = first)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 4,
                                      med_level = 6, med_threshold = 0.1,
                                      max_level = 7, max_threshold = 0.5)

# positivity limiter for the water height
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 1,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true,
                           limiter! = stage_limiter!)

stepsize_callback = StepsizeCallback(cfl = 0.2)


save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = output_dir,
                                     solution_variables = cons2prim)


callbacks = CallbackSet(amr_callback, stepsize_callback, save_solution);
###############################################################################

sol = solve(ode, SSPRK43(;stage_limiter!); dt = 1.0, callback = callbacks, adaptive = false)

# Create an animation of the solution

# To visualize the solution and bathymetry we post-processing the Trixi.jl output file(s)
# with the Trixi2Vtk.jl functionality and plot them with ParaView.

#Trixi.save_mesh_file(mesh_1, output_dir)
trixi2vtk(joinpath(output_dir, "solution_*.h5"), output_directory = output_dir)

# It is possible to open the created solution_00000.pvd file with ParaView and create a video of the simulation.

# In ParaView, after opening the solution_00000.pvd file, one can apply two instances
# of the Warp By Scalar filter to visualize the water height and bathymetry in three dimensions.
# Many additional customizations, e.g., color scaling, fonts, etc. are available in ParaView.

# A created video of the simulation can be found here: have to fill with a link to the video
# In addition there is also a top view of the solution available here: have to fill with a link to the video
