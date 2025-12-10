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
BC = (; type="imp", param=BigFloat("1.45"))
r0 = BigFloat("13000");    # Radius of curvature

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


# First even mode:
z_mu_even1, zC_even1, zD_even1, _, _ = solve_bent_waveguide(EVEN1, a, b, r0, wavelength, w_num_core, w_num_clad, 
    BC=BC, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
 
# First odd mode:
z_mu_odd1, zC_odd1, zD_odd1, _, _ = solve_bent_waveguide(ODD1, a, b, r0, wavelength, w_num_core, w_num_clad, 
    BC=BC, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)

# Second even mode:
z_mu_even2, zC_even2, zD_even2, _, _ = solve_bent_waveguide(EVEN2, a, b, r0, wavelength, w_num_core, w_num_clad, 
    BC=BC, tol_bessel=tol_bessel, tol_eigen=tol_eigen, max_itr=max_itr, num_terms=num_terms, print_flag=print_flag)
 

z_mu_all = [z_mu_even1, z_mu_odd1, z_mu_even2]

# Results in Table 5:
lambda_all  = z_mu_all * r0^2
beta_all    = r0*sqrt.(z_mu_all)
sqrt_mu_all = sqrt.(z_mu_all)

println("lambda values:")
display(lambda_all)

println("beta values:")
display(lambda_all)

println("sqrt(mu) values:")
display(lambda_all)


 
