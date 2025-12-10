using Plots          # For visualization

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
mode = EVEN1 # Valid options: EVEN1, ODD1, EVEN2
BC_pml = (; type="pml", param=BigFloat(800))
BC_imp = (; type="imp", param=BigFloat("1.45"))

r0 = BigFloat("5200");    # Radius of curvature
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

println("Computing " * mode.type * mode.index * " mode: =================================================================")

z_mu_final, zC, zD, success, params = solve_bent_waveguide(mode, a, b, r0, wavelength, w_num_core, w_num_clad, 
    BC=BC, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
  

# For plotting
subdivisions = 256
dr = LinRange(-b, b, subdivisions + 1)
zr = r0 .+ dr 

# Initialize arrays
zu_results  = zeros(Complex{BigFloat}, size(dr))
zdu_results = zeros(Complex{BigFloat}, size(dr))
for l in eachindex(dr)
    zu_results[l], zdu_results[l] = eval_waveguide_solution( zr[l], z_mu_final, zC, zD, params, print_flag=print_flag )
end


fig1 = plot(
    xtickfontsize=20,
    ytickfontsize=20,
    legendfontsize=28,
    size=(1200, 1000),
    xlim=(-b,b),
    # ylim=(-1.5, 1.5),
    legend=:outerbottom,
    title="u(r - r0) (" * mode.type * mode.index * ")\n",
    titlefontsize=42
)
hline!(fig1, [0.0], color=:black, linewidth=3, label="")
vline!(fig1, [0.0], color=:black, linewidth=3, label="")
vline!(fig1, [-a, a], color=:black, linewidth=4, linestyle=:dot, label="")
plot!(fig1, dr, real.(zu_results), color=:blue, linewidth=8, label="Real Component")
plot!(fig1, dr, imag.(zu_results), color=:red, linewidth=8, label="Complex Component")
display(fig1)



fig2 = plot(
    xtickfontsize=20,
    ytickfontsize=20,
    legendfontsize=28,
    size=(1200, 1000),
    xlim=(-b, b),
    # ylim=(-8,8),
    legend=:outerbottom,
    title="du(r - r0) (" * mode * ")\n",
    titlefontsize=42
)
hline!(fig2, [0.0], color=:black, linewidth=3, label="")
vline!(fig2, [0.0], color=:black, linewidth=3, label="")
vline!(fig2, [-a, a], color=:black, linewidth=4, linestyle=:dot, label="")
plot!(fig2, dr, real.(zdu_results), color=:blue, linewidth=8, label="Real Component")
plot!(fig2, dr, imag.(zdu_results), color=:red, linewidth=8, label="Complex Component")
display(fig2)
