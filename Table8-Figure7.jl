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
BC_imp = (; type="imp", param=BigFloat("1.45"))
BC_pml = (; type="pml", param=BigFloat("800"))
r0_all = [BigFloat("10400"), BigFloat("7800"), BigFloat("5200"), BigFloat("2600")];    # Radius of curvature

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

# Allocate memory for results
z_mu_all_imp = zeros(Complex{BigFloat}, 12)
zC_all_imp = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]
zD_all_imp = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]

z_mu_all_pml = zeros(Complex{BigFloat}, 12)
zC_all_pml = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]
zD_all_pml = @NamedTuple{L::Complex{BigFloat}, core::Complex{BigFloat}, R::Complex{BigFloat}}[]

# First even mode 
i0 = 0
for i in eachindex(r0_all)
    # Run with impedance BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(EVEN1, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_imp, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)

    z_mu_all_imp[i + i0] = success ? z_mu : 0.0
    push!(zC_all_imp, zC)
    push!(zD_all_imp, zD)

    # Run with PML BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(EVEN1, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_pml, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
    
    z_mu_all_pml[i + i0] = success ? z_mu : 0.0
    push!(zC_all_pml, zC)
    push!(zD_all_pml, zD)
end


# First odd mode 
i0 = length(r0_all)
for i in eachindex(r0_all)
    # Run with impedance BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(ODD1, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_imp, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
 
    z_mu_all_imp[i + i0] = success ? z_mu : 0.0
    push!(zC_all_imp, zC)
    push!(zD_all_imp, zD)

    # Run with PML BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(ODD1, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_pml, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
 
    z_mu_all_pml[i + i0] = success ? z_mu : 0.0
    push!(zC_all_pml, zC)
    push!(zD_all_pml, zD)
end

# Second even mode 
i0 = 2*length(r0_all)
for i in eachindex(r0_all)
    # Run with impedance BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(EVEN2, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_imp, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)

    z_mu_all_imp[i + i0] = success ? z_mu : 0.0
    push!(zC_all_imp, zC)
    push!(zD_all_imp, zD)

    # Run with PML BC
    z_mu, zC, zD, success, _ = solve_bent_waveguide(EVEN2, a, b, r0_all[i], wavelength, w_num_core, w_num_clad, 
        BC=BC_pml, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
    
    z_mu_all_pml[i + i0] = success ? z_mu : 0.0
    push!(zC_all_pml, zC)
    push!(zD_all_pml, zD)
end


# Table 8:
beta_all_imp = [r0_all; r0_all; r0_all] .*sqrt.(z_mu_all_imp)
beta_all_pml = [r0_all; r0_all; r0_all] .*sqrt.(z_mu_all_pml)

println("beta - Impedance Boundary Condtion:")
display(beta_all_imp)

println("beta - PML Boundary Condtion:")
display(beta_all_pml)


# Figure 7
subfigures = []
subfigure_titles = ["(a) First even mode", "(b) First odd mode", "(c) Second even mode"]

r0_fig7 = BigFloat("2600")
i = 1
for mode in [EVEN1, ODD1, EVEN2]

    z_mu, zC, zD, success, params = solve_bent_waveguide(mode, a, b, r0_fig7, wavelength, w_num_core, w_num_clad, 
        BC=BC_pml, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)

    subdivisions = 256
    dr = LinRange(-b, b, subdivisions + 1)
    zr = r0_fig7 .+ dr 

    # Initialize arrays
    zu_results  = zeros(Complex{BigFloat}, size(dr))
    zdu_results = zeros(Complex{BigFloat}, size(dr))
    for l in eachindex(dr)
        zu_results[l], zdu_results[l] = eval_waveguide_solution( zr[l], z_mu, zC, zD, params, print_flag=print_flag )
    end

    fig = plot(xlabel="r\n\n$(subfigure_titles[i]) ", 
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
        ylim=i == 2 ? (-0.3, 0.3) : (-1.5, 1.5),
        xticks=-b:1:b,
        yticks= i == 2 ? (-0.3:0.1:0.3) : (-1.5:0.5:1.5),
        gridalpha=0.5)
    vline!(fig, [-a, a], color=:blue, linewidth=4, linestyle=:dot, label="")
    plot!(fig, dr, real.(zu_results), color=:black, linewidth=8, label="Real Component")
    push!(subfigures, fig)

    i+=1
end

figure7 = plot(subfigures[1], subfigures[2], subfigures[3], layout=(1,3), size=(3000, 1600), plot_title=@sprintf("Predicted Modes with PML Boundary Condition (PML Const=%.1f) for r0=%d\n", BC_pml.param, r0_fig7), plot_titlefontsize=40, bottom_margin=30mm)
display(figure7)
png("figures/figure7.png")