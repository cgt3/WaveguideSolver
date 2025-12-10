using Plots
using Plots.PlotMeasures
using Printf

include("bessel_solver.jl")


# Set the precision of BigFloat and printing:
setprecision(70, base=10) 
print("\033c") 

# Solver parameters:
print_flag = false
tol_bessel = BigFloat("1e-65")
tol_eigen = BigFloat("1e-60") 
max_itr = 100

# Problem parameters:
BC = (; type="pml", param=BigFloat("800"))
r0_all = [BigFloat("10400"), BigFloat("7800"), BigFloat("5200") ];    # Radius of curvature

a = BigFloat("0.5");      # Distance from centerline to the interior wall of the cladding
b = BigFloat("5.0");      # Distance from centerline to the exterior wall of the cladding

# Material parameters:
n_core = BigFloat("1.4512")  
n_clad = BigFloat("1.45")

# Reference parameters:
ref_length = BigFloat("25.4e-6")
wavelength = BigFloat("1064e-9") / ref_length

k0 = BigFloat("149.993333460866")
w_num_core = k0*n_core
w_num_clad = k0*n_clad

# Precomputing parameters: 
num_terms = 1000 # Max number of iterations and number of precomputed terms for memoization/caching

subfigures = []
subfigure_letter = ["a", "b", "c"]

z_mu_all = zeros(Complex{BigFloat}, 3)
zC_all = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]
zD_all = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]
for i in eachindex(r0_all)
    z_mu_all[i], zC, zD, _, params = solve_bent_waveguide(EVEN2, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
    push!(zC_all, zC)
    push!(zD_all, zD)

    # For plotting
    subdivisions = 256
    dr = LinRange(-b, b, subdivisions + 1)
    zr = r0_all[i] .+ dr 

    # Initialize arrays
    zu_results  = zeros(Complex{BigFloat}, size(dr))
    zdu_results = zeros(Complex{BigFloat}, size(dr))
    for l in eachindex(dr)
        zu_results[l], zdu_results[l] = eval_waveguide_solution( zr[l], z_mu_all[i], zC, zD, params, print_flag=print_flag )
    end

    fig = plot(xlabel="r\n\n($(subfigure_letter[i])) r0 = $(floor(Integer, r0_all[i]))", 
        ylabel="Real part of u",
        framestyle=:box, 
        titlefontsize=32, 
        xlabelfontsize=36, 
        ylabelfontsize=32, 
        top_margin=20mm, 
        left_margin=30mm, 
        xtickfontsize=24, 
        ytickfontsize=24, 
        leg=false, 
        xlim=(-b,b), 
        ylim=(-1.5, 1.5),
        xticks=-b:1:b,
        yticks=-1.5:0.5:1.5,
        fg_color_grid=:grey60, 
        gridalpha=1.0)
    vline!(fig, [-a, a], color=:blue, linewidth=4, linestyle=:dot, label="")
    plot!(fig, dr, real.(zu_results), color=:black, linewidth=8, label="Real Component")
    push!(subfigures, fig)
end

# Results shown in Table 7:
lambda_all  = z_mu_all .* r0_all .^2
beta_all    = r0_all .*sqrt.(z_mu_all)
sqrt_mu_all = sqrt.(z_mu_all)

# Figure 6:
figure6 = plot(subfigures[1], subfigures[2], subfigures[3], layout=(1,3), size=(3000, 1600), plot_title=@sprintf("Second Even Mode with PML Boundary Condition (PML Const=%.1f)\n", BC.param), plot_titlefontsize=48, bottom_margin=30mm)
display(figure6)
png("figures/figure6.png")
