
# Load libraries and helper files: ================================================================
# using BenchmarkTools # For timings (if desired)
# using Polyester      # For parallel for-loops (if desired)

using Plots          # For visualization

include("bessel_solver.jl")

print("\033c") 

MATLAB_DOUBLE_PREC = "matlab"
FORTRAN_QUAD_PREC = "fortran"
ARBITRARY_PREC = "arbitrary"

print_flag = false
precision_mode = ARBITRARY_PREC

# Precision specifier for BigFloat: ===============================================================

if precision_mode == MATLAB_DOUBLE_PREC
    # For testing against MatLab (double precision):
    setprecision(53, base=2) # For comparison with the MatLab code
    tol_bessel = BigFloat("1e-15")   # Convergence tolerance for the Bessel function evaluation
    tol_eigen = BigFloat("1e-13") # MatLab: 1e-13
    max_itr = 50
elseif precision_mode == FORTRAN_QUAD_PREC
    # For testing against Fortran (quadruple precision)
    setprecision(113, base=2) 
    tol_bessel = BigFloat("1e-33")   # Convergence tolerance for the Bessel function evaluation
    tol_eigen = BigFloat("1e-26") 
    max_itr = 50
else
    # High precision:
    setprecision(70, base=10) 
    tol_bessel = BigFloat("1e-65")   # Convergence tolerance for the Bessel function evaluation
    tol_eigen = BigFloat("1e-60") 
    max_itr = 50
end


# Problem parameters: =============================================================================
mode = "even2"           # "even1", "odd1" or "even2"
bctype = "pml"           # "imp" or "pml"
r0 = BigFloat("5200");    # Radius of curvature

a = BigFloat("0.5");      # ?: Distance from centerline to the interior wall of the cladding
b = BigFloat("5.0");      # ?: Distance from centerline to the exterior wall of the cladding

# Material parameters:
n_core = BigFloat("1.4512")  
n_clad = BigFloat("1.45")

# Reference parameters:
ref_length = BigFloat("25.4e-6")
wavelength = BigFloat("1064e-9") / ref_length

k0 = BigFloat("149.993333460866") #2*BigFloat(pi) / wavelength
w_num_core = k0*n_core
w_num_clad = k0*n_clad


# Possibly complex variables:
z_mu = Complex{BigFloat}(BigFloat("0.217655321571770e3")^2, 0)


# Precomputing: ===================================================================================
println("Precomputing array of factors for bessel_frobenius routine...\n")

num_terms = 1000 # Max number of iterations and number of precomputed terms for memoization/caching

factors0, factorsL, factorsR = get_factors(num_terms, r0, a)


# Eigenvalue Solver: ==============================================================================
even1 = 1
even2 = 2
odd1 = 3

println("Computing " * mode * " mode: =================================================================")
if mode == "even1"
    i_mode = even1
    zC_final = Complex{BigFloat}(1.0, 0.0)
    zD_guess = Complex{BigFloat}(0.0, 0.0)

    if precision_mode == MATLAB_DOUBLE_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473785763924115e5"), 0.0) # For MatLab
    else   #if precision_mode == FORTRAN_QUAD_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473785763924114535111724293851737e5"), 0.0)  # For Fortran
    # else
    end

    z_mu_final, zD_final, success = find_even_mode(bctype,z_mu_guess, zD_guess, w_num_core, w_num_clad, r0, a, b, 
            tol_bessel, num_terms, factors0, factorsL, factorsR, max_itr, tol_eigen, print_flag=print_flag)
elseif mode == "even2"
    i_mode = even2
    zC_final = Complex{BigFloat}(1.0, 0.0)
    zD_guess = Complex{BigFloat}(0.0, 0.0)

    if precision_mode == MATLAB_DOUBLE_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473251454095355e5"), 0.0) # For MatLab
    else  #if precision_mode == FORTRAN_QUAD_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473251454095354947292257785388791225e5"), 0.0) # For Fortran
    # else

    end


    z_mu_final, zD_final, success = find_even_mode(bctype,z_mu_guess, zD_guess, w_num_core, w_num_clad, r0, a, b, 
            tol_bessel, num_terms, factors0, factorsL, factorsR, max_itr, tol_eigen, print_flag=print_flag)
elseif mode == "odd1"
    i_mode = odd1
    zC_guess = Complex{BigFloat}(0.0, 0.0)
    zD_final = Complex{BigFloat}(1.0, 0.0)

    if precision_mode == MATLAB_DOUBLE_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473546725588373e5"), 0.0) # For MatLab
    else #if precision_mode == FORTRAN_QUAD_PREC
        z_mu_guess = Complex{BigFloat}(BigFloat("0.473546725588372858352299169480810e5"), 0.0) # For Fortran
    # else

    end

    z_mu_final, zC_final, success = find_odd_mode(bctype,z_mu_guess, zC_guess, w_num_core, w_num_clad, r0, a, b, 
            tol_bessel, num_terms, factors0, factorsL, factorsR, max_itr, tol_eigen, print_flag=print_flag)
else
    @warn "Unknown mode"
end

if success
    println("CONVERGED SOLUTION: " * mode," -------- bctype = " * bctype," -------- r0 = $r0")
    println("    z_mu_final = $z_mu_final")
    println("    lambda_fnl = $(z_mu_final*r0^2)")
    println("    beta_final = $(r0*sqrt(z_mu_final))")
    println("    beta_fnl/r0= $(sqrt(z_mu_final))")
    println("    zC_final   = $zC_final")
    println("    zD_final   = $zD_final")
    println("    loss_final = $(-imag(r0*sqrt(z_mu_final)))")
else
    println("Failed to converge")
end


# Plotting: =======================================================================================

# Get the solutions in the left and right cladding
compute_deriv = false
zCL, zDL, _, _, _, _, _, _ = get_C_D_left(z_mu_final, zC_final, zD_final, w_num_core, r0, a, 
                                            compute_deriv, tol_bessel, num_terms, factors0, print_flag=print_flag)

zCR, zDR, _, _, _, _, _, _ = get_C_D_right(z_mu_final, zC_final, zD_final, w_num_core, r0, a, 
                                            compute_deriv, tol_bessel, num_terms, factors0, print_flag=print_flag)

                                
# For plotting
subdivisions = 256
dr = LinRange(-b, b, subdivisions + 1)
zr = r0 .+ dr 

# Initialize arrays
zu_results  = zeros(Complex{BigFloat}, size(dr))
zdu_results = zeros(Complex{BigFloat}, size(dr))
for l in eachindex(dr)
    zu_results[l], zdu_results[l] = BentThreeLayerWaveguide( zr[l], z_mu_final, zC_final, 
                                        zD_final, zCL, zDL, zCR, zDR, w_num_core, w_num_clad,
                                        r0, a, b, tol_bessel, num_terms, 
                                        factors0, factorsL, factorsR, print_flag=print_flag )
end


fig1 = plot(
    xtickfontsize=20,
    ytickfontsize=20,
    legendfontsize=28,
    size=(1200, 1000),
    xlim=(-b,b),
    # ylim=(-1.5, 1.5),
    legend=:outerbottom,
    title="Julia: u(r - r0) (" * mode * ")\n",
    titlefontsize=42
)
hline!(fig1, [0.0], color=:black, linewidth=3, label="")
vline!(fig1, [0.0], color=:black, linewidth=3, label="")
vline!(fig1, [-a, a], color=:black, linewidth=4, linestyle=:dot, label="")
plot!(fig1, dr, real.(zu_results), color=:blue, linewidth=8, label="Real Component")
plot!(fig1, dr, imag.(zu_results), color=:red, linewidth=8, label="Complex Component")
display(fig1)
savefig("figures/" * mode * "_u.png")



fig2 = plot(
    xtickfontsize=20,
    ytickfontsize=20,
    legendfontsize=28,
    size=(1200, 1000),
    xlim=(-b, b),
    # ylim=(-8,8),
    legend=:outerbottom,
    title="Julia: du(r - r0) (" * mode * ")\n",
    titlefontsize=42
)
hline!(fig2, [0.0], color=:black, linewidth=3, label="")
vline!(fig2, [0.0], color=:black, linewidth=3, label="")
vline!(fig2, [-a, a], color=:black, linewidth=4, linestyle=:dot, label="")
plot!(fig2, dr, real.(zdu_results), color=:blue, linewidth=8, label="Real Component")
plot!(fig2, dr, imag.(zdu_results), color=:red, linewidth=8, label="Complex Component")
display(fig2)
savefig("figures/" * mode * "_du.png")
