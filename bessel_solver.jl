using LinearAlgebra

function get_factors(num_terms, r0::BigFloat, a::BigFloat)
    n_factors = num_terms - 2
    factors0 = ones(BigFloat, n_factors)
    factorsL = ones(BigFloat, n_factors)
    factorsR = ones(BigFloat, n_factors)

    # Compute the powers and factorials in order to save time
    b0 = BigFloat(2.0) /  r0
    bL = BigFloat(2.0) / (r0 - a)
    bR = BigFloat(2.0) / (r0 + a)
    for l in num_terms - 4:-1:0 # LB should be 0 or 1?
        factors0[l+1] = b0 / BigFloat(num_terms - 3 - l) * factors0[l+2] 
        factorsL[l+1] = bL / BigFloat(num_terms - 3 - l) * factorsL[l+2]
        factorsR[l+1] = bR / BigFloat(num_terms - 3 - l) * factorsR[l+2]
    end

    return factors0, factorsL, factorsR
end

function rfact(N, r0)
    if N == 0 
        return BigFloat(1.0)
    end

    prod = BigFloat(1.0)
    for j in 1:N 
        prod *= BigFloat(2) / BigFloat(j) / BigFloat(r0)
    end

    return prod
end

# Params: ----------------------------------------------------------------------------------------
#        zr - the desired radius to evaluate the Bessel function at
#       z_mu - 
#      w_num - wave number
#     r_base - the radius about which to expand
#       zc0 - base case for coefficent c0 
#       zc1 - base case for coefficent c1 
# compute_du - whether to compute du or not (idec in MatLab code)
#        tol - the tolerance for the expansion
#    factors - precomputed terms
# ------------------------------------------------------------------------------------------------
function bessel_frobenius(zr::Union{BigFloat, Complex{BigFloat}}, 
                          z_mu::Union{BigFloat, Complex{BigFloat}}, 
                          w_num::BigFloat, 
                          r_base::BigFloat, 
                          zc0::Complex{BigFloat}, 
                          zc1::Complex{BigFloat}, 
                          compute_dudmu::Bool, 
                          tol::BigFloat, 
                          num_terms::Integer, 
                          factors::AbstractArray{BigFloat}; 
                          print_flag=false::Bool)

    omega = w_num * w_num

    z_x = r_base * log(zr / r_base)
    z_c = zeros(Complex{BigFloat}, num_terms)

    z_c[1] = zc0
    z_c[2] = zc1

    converged_prev = false
    converged = false

    # Compute f as a power series in the transformed variable
    z_f = zc0 + zc1 * z_x;
  
    z_df = zc1
    z_delta = Complex{BigFloat}(0.0, 0.0)

    for n in 0:num_terms - 3
        z_sum = factors[num_terms-2-n:num_terms - 2]' * z_c[1:n+1]

        z_c[n + 2 + 1] = (z_mu * z_c[n + 1] - omega * z_sum) / BigFloat(n+1) / BigFloat(n+2)
        z_delta = z_c[n + 2 + 1] * z_x^(n+2)

        z_f = z_f + z_delta
        z_df = z_df + z_c[n+2 + 1] * z_x^(n+1) * (n+2)

        if abs(z_delta) < tol 
            converged = true
        else
            converged = false
        end

        if converged && converged_prev
            if print_flag
                println("bessel_frobenius: Converged at n=$n")
            end
            break
        else
            converged_prev = converged
        end
    end

    if !converged_prev || !converged 
        if print_flag
            println("bessel_frobenius: f has not converged in $num_terms terms.")
            println("    Last relative magnitude: $z_delta")
        end
    end

    z_f_mu = Complex{BigFloat}(0.0, 0.0)
    z_df_mu = Complex{BigFloat}(0.0, 0.0)
    converged_prev = false
    if compute_dudmu
        z_b = zeros(Complex{BigFloat}, num_terms)
        z_f_mu = z_b[1] + z_b[2] * z_x
        z_df_mu = z_b[2]

        z_x_pow = z_x
        for n in 0:num_terms - 3
            z_sum = factors[num_terms-2-n:num_terms - 2]' * z_b[1:n+1]

            z_b[n + 2 + 1] = ( z_c[n + 1] + z_mu * z_b[n+1] - omega * z_sum ) / ((n+1)*(n+2))
            z_delta = z_b[n+2+1] * (z_x * z_x_pow)


            z_f_mu = z_f_mu + z_delta
            z_df_mu = z_df_mu + z_b[n+2 + 1] * z_x_pow*BigFloat(n+2) 
            z_x_pow *= z_x

            if abs(z_delta) < tol 
                converged = true
            else
                converged = false
            end

            if converged && converged_prev
                break 
            else
                converged_prev = converged 
            end
        end
    end

    return z_f, z_df, z_f_mu, z_df_mu
end

# Evaluation of solutions and all derivatives at the left boundary
function BC_left(z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, tol_bessel, num_terms, factors0, factorsL; print_flag=false::Bool)

    compute_bessel_deriv = true
    zCL, zDL, zCL_mu, zDL_mu, zCL_C, zDL_C, zCL_D, zDL_D = get_C_D_left(z_mu, zC, zD, w_num_core, r0, a, 
        compute_bessel_deriv, tol_bessel, num_terms, factors0, print_flag=print_flag ) 

    # Evaluate the cladding solution at the left boundary 
    zr       = r0 - b 
    w_num    = w_num_clad 
    r_base   = r0 - a 
    z_mu_loc = z_mu * r0^2 / (r_base^2)

    # Solution 10
    zc0 = Complex{BigFloat}(1.0, 0.0)
    zc1 = Complex{BigFloat}(0.0, 0.0)

    f10, df10, f10_mu, df10_mu = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, 
                                                  compute_bessel_deriv, tol_bessel, num_terms, factorsL, print_flag=print_flag)
    
    # Solution 01
    zc0 = Complex{BigFloat}(0.0, 0.0)
    zc1 = Complex{BigFloat}(1.0, 0.0)

    f01, df01, f01_mu, df01_mu = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, 
                                                  compute_bessel_deriv, tol_bessel, num_terms, factorsL, print_flag=print_flag)
    val    = zCL   * f10 + zDL   * f01
    val_C  = zCL_C * f10 + zDL_C * f01
    val_D  = zCL_D * f10 + zDL_D * f01
    val_mu = zCL * f10_mu + zCL_mu * f10 + zDL * f01_mu + zDL_mu * f01

    dval   = zCL   * df10 + zDL   * df01
    dval_C = zCL_C * df10 + zDL_C * df01
    dval_D = zCL_D * df10 + zDL_D * df01
    dval_mu = zCL * df10_mu + zCL_mu * df10 + zDL * df01_mu + zDL_mu * df01

    return val,    dval, 
           val_mu, dval_mu,
           val_C,  dval_C,
           val_D,  dval_D
end

# Evaluation of solutions and all derivatives at the right boundary
function BC_right(bctype,z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, tol_bessel, num_terms, factors0, factorsR; print_flag=false::Bool)

    # println()
    compute_bessel_deriv = true
    zCR, zDR, zCR_mu, zDR_mu, zCR_C, zDR_C, zCR_D, zDR_D = get_C_D_right(z_mu, zC, zD, w_num_core, r0, a, 
        compute_bessel_deriv, tol_bessel, num_terms, factors0, print_flag=print_flag) 

    # Evaluate the cladding solution at the left boundary
    if bctype == "imp"
        zr       = r0 + b
    elseif bctype == "pml"
        PML_const = BigFloat(800)
        zr       = r0 + b - PML_const / w_num_clad * im 
    end
    
    w_num    = w_num_clad 
    r_base   = r0 + a 
    z_mu_loc = z_mu * r0^2 / (r_base^2)

    # Solution 10
    zc0 = Complex{BigFloat}(1.0, 0.0)
    zc1 = Complex{BigFloat}(0.0, 0.0)

    f10, df10, f10_mu, df10_mu = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, 
                                                  compute_bessel_deriv, tol_bessel, num_terms, factorsR, print_flag=print_flag)
    
    # Solution 01
    zc0 = Complex{BigFloat}(0.0, 0.0)
    zc1 = Complex{BigFloat}(1.0, 0.0)

    f01, df01, f01_mu, df01_mu = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, 
                                                  compute_bessel_deriv, tol_bessel, num_terms, factorsR, print_flag=print_flag)
    
    val    = zCR   * f10 + zDR   * f01
    val_C  = zCR_C * f10 + zDR_C * f01
    val_D  = zCR_D * f10 + zDR_D * f01

    val_mu = zCR * f10_mu + zCR_mu * f10 + zDR * f01_mu + zDR_mu * f01

    dval   = zCR   * df10 + zDR   * df01
    dval_C = zCR_C * df10 + zDR_C * df01
    dval_D = zCR_D * df10 + zDR_D * df01
    dval_mu = zCR * df10_mu + zCR_mu * df10 + zDR * df01_mu + zDR_mu * df01

    return val,    dval, 
           val_mu, dval_mu,
           val_C,  dval_C,
           val_D,  dval_D
end


function check_rel_convergence(a1, a2, denom, rel_tol)
    if abs(a1 - a2) / abs(denom) < rel_tol
        return true
    else
        return false
    end
end

function find_even_mode(bctype::String,
                        z_mu_guess::Complex{BigFloat}, 
                        zD_guess::Complex{BigFloat}, 
                        w_num_core::BigFloat, 
                        w_num_clad::BigFloat, 
                        r0::BigFloat,
                        a::BigFloat, 
                        b::BigFloat, 
                        tol_bessel::BigFloat, 
                        num_terms::Integer,
                        factors0::Vector{BigFloat}, 
                        factorsL::Vector{BigFloat}, 
                        factorsR::Vector{BigFloat}, 
                        max_itr::Integer, 
                        tol_eigen::BigFloat; 
                        print_flag=false::Bool)

    success = false 
    zC = Complex{BigFloat}(1.0, 0.0)

    sol_guess = [ z_mu_guess; zD_guess ]
    ref_norm = norm(sol_guess)

    F = zeros(Complex{BigFloat},2,1)
    Jac = zeros(Complex{BigFloat},2,2)

    sol = sol_guess
    for itr in 1:max_itr
        
        z_mu = sol[1]
        zD = sol[2]

        # Evaluate at the left boundary
        _, dval, _, dval_mu, _, _, _, dval_D = BC_left(z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, 
                                                       tol_bessel, num_terms, factors0, factorsL, print_flag=print_flag)

        # Save value of the first BC and derivatives wrt mu and the unknown coef 
        F[1] = dval 
        Jac[1,:] = [dval_mu dval_D]

        # Evaluate right boundary
        val, dval, val_mu, dval_mu, _, _, val_D, dval_D = BC_right(bctype,z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, 
                                                                   tol_bessel, num_terms, factors0, factorsR, print_flag=print_flag)
        # Save value of the second BC and derivatives wrt mu and the unknown coef 
        if bctype == "imp"
            F[2] = (r0+a)/(r0+b)*dval + 1im*w_num_clad*val
            Jac[2,1] = (r0+a)/(r0+b)*dval_mu + 1im*w_num_clad*val_mu
            Jac[2,2] = (r0+a)/(r0+b)*dval_D + 1im*w_num_clad*val_D
        elseif bctype == "pml"
            # RESTART: Jac[2,1] is off between the two, and it corresponds to val_mu from BC_right
            F[2] = val 
            Jac[2,:] = [val_mu val_D]
        end

        # Increment the solution
        zdet = Jac[1,1] * Jac[2,2] - Jac[2,1]Jac[1,2]
        zdet1 = F[1] * Jac[2,2] - F[2] * Jac[1,2]
        zdet2 = F[2] * Jac[1,1] - F[1] * Jac[2,1]

        delta_sol = -[ zdet1 / zdet ; zdet2 / zdet]
        # delta_sol = -Jac \ F
        sol = sol + delta_sol

        if print_flag
            println("find_even_mode: itr=$itr")
            println("    F[1]     = $(F[1])")
            println("    F[2]     = $(F[2])\n")
            println("    delta mu = $(delta_sol[1])")
            println("    delta_D  = $(delta_sol[2])\n")
            println("    mu       = $(sol[1])")
            println("    D        = $(sol[2])\n")

            println("  Jac:")
            println("    [1,1] = $(Jac[1,1])")
            println("    [1,2] = $(Jac[1,2])")
            println("    [2,1] = $(Jac[2,1])")
            println("    [2,2] = $(Jac[2,2])")

        end

        println("find_even_mode: nonlinear itr=$itr, relative_delta=$(norm(delta_sol)/ref_norm)")

        if norm(delta_sol) / ref_norm < tol_eigen
            success = true
            break
        end
    end

    z_mu_final = sol[1]
    zD_final = sol[2]
    
    return z_mu_final, zD_final, success
end

function find_odd_mode(bctype::String,
                       z_mu_guess::Complex{BigFloat}, 
                       zC_guess::Complex{BigFloat}, 
                       w_num_core::BigFloat, 
                       w_num_clad::BigFloat, 
                       r0::BigFloat, 
                       a::BigFloat, 
                       b::BigFloat, 
                       tol_bessel::BigFloat, 
                       num_terms::Integer, 
                       factors0::Vector{BigFloat}, 
                       factorsL::Vector{BigFloat}, 
                       factorsR::Vector{BigFloat}, 
                       max_itr::Integer, 
                       tol_eigen::BigFloat; 
                       print_flag=false::Bool)
    success = false 
    zD = Complex{BigFloat}(1.0, 0.0)

    sol_guess = [ z_mu_guess; zC_guess ]
    ref_norm = norm(sol_guess)

    F = zeros(Complex{BigFloat},2,1)
    Jac = zeros(Complex{BigFloat},2,2)

    sol = sol_guess
    for itr in 1:max_itr
        z_mu = sol[1]
        zC = sol[2]


        # Evaluate at the left boundary
        _, dval, _, dval_mu, _, dval_C, _, _ = BC_left(z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, 
                                                       tol_bessel, num_terms, factors0, factorsL, print_flag=print_flag)

        # Save value of the first BC and derivatives wrt mu and the unknown coef 
        F[1] = dval 
        Jac[1,:] = [dval_mu dval_C]

        # Evaluate right boundary
        val, dval, val_mu, dval_mu, val_C, dval_C, _, _ = BC_right(bctype,z_mu, zC, zD, w_num_core, w_num_clad, r0, a, b, 
                                                                   tol_bessel, num_terms, factors0, factorsR, print_flag=print_flag)

        # Save value of the second BC and derivatives wrt mu and the unknown coef 
        if bctype == "imp"
            F[2] = (r0+a)/(r0+b)*dval + 1im*w_num_clad*val
            Jac[2,1] = (r0+a)/(r0+b)*dval_mu + 1im*w_num_clad*val_mu
            Jac[2,2] = (r0+a)/(r0+b)*dval_C + 1im*w_num_clad*val_C
        elseif bctype == "pml"
            F[2] = val 
            Jac[2,:] = [val_mu val_C]
        end

        # Increment the solution
        zdet = Jac[1,1] * Jac[2,2] - Jac[2,1]Jac[1,2]
        zdet1 = F[1] * Jac[2,2] - F[2] * Jac[1,2]
        zdet2 = F[2] * Jac[1,1] - F[1] * Jac[2,1]

        delta_sol = -[ zdet1 / zdet ; zdet2 / zdet]

        # delta_sol = -Jac \ F
        sol = sol + delta_sol

        if print_flag
            println("find_odd_mode: itr=$itr")
            println("    F[1]     = $(F[1])")
            println("    F[2]     = $(F[2])\n")
            println("    delta mu = $(delta_sol[1])")
            println("    delta_C  = $(delta_sol[2])\n")
            println("    mu       = $(sol[1])")
            println("    C        = $(sol[2])\n")
            println("\n\n")

            println("  Jac:")
            println("    [1,1] = $(Jac[1,1])")
            println("    [1,2] = $(Jac[1,2])")
            println("    [2,1] = $(Jac[2,1])")
            println("    [2,2] = $(Jac[2,2])")
        end


        println("find_odd_mode: nonlinear itr=$itr, relative_delta=$(norm(delta_sol)/ref_norm)")

        if norm(delta_sol) / ref_norm < tol_eigen
            success = true
            break
        end
    end

    z_mu_final = sol[1]
    zC_final = sol[2]
    
    return z_mu_final, zC_final, success
end

function get_C_D_left(z_mu::Complex{BigFloat}, zC::Complex{BigFloat}, zD::Complex{BigFloat}, w_num_core::BigFloat, r0::BigFloat, a::BigFloat, compute_deriv::Bool, tol_bessel::BigFloat, num_terms::Integer, factors0::Vector{BigFloat}; print_flag=false::Bool)
    zr = r0 - a 
    r_base = r0 
    w_num = w_num_core

    # Core solution 10
    zc0 = Complex{BigFloat}(1.0, 0.0)
    zc1 = Complex{BigFloat}(0.0, 0.0)

    f10, df10, f10_mu, df10_mu = bessel_frobenius(zr, z_mu, w_num, r_base, zc0, zc1, compute_deriv, 
                                                  tol_bessel, num_terms, factors0, print_flag=print_flag)

    # Core solution 01
    zc0 = Complex{BigFloat}(0.0, 0.0)
    zc1 = Complex{BigFloat}(1.0, 0.0)

    f01, df01, f01_mu, df01_mu = bessel_frobenius(zr, z_mu, w_num, r_base, zc0, zc1, compute_deriv, 
                                                  tol_bessel, num_terms, factors0, print_flag=print_flag)


    # Compute the value of u & du/dx at the mid-left interface as a linear comb. of the 10 & 01 solutions
    aux  = zC * f10  + zD * f01 
    daux = zC * df10 + zD * df01

    zCL = aux 
    zDL = daux * r0 / (r0-a)

    if compute_deriv
        zCL_mu = zC * f10_mu + zD * f01_mu
        zDL_mu = (zC * df10_mu + zD * df01_mu) * r0 / (r0-a)

        zCL_C = f10 
        zDL_C = df10 * r0 / zr

        zCL_D = f01 
        zDL_D = df01 * r0 / zr
    else
        zCL_mu = Complex{BigFloat}(0.0, 0.0)
        zDL_mu = Complex{BigFloat}(0.0, 0.0)
        zCL_C  = Complex{BigFloat}(0.0, 0.0)
        zDL_C  = Complex{BigFloat}(0.0, 0.0)
        zCL_D  = Complex{BigFloat}(0.0, 0.0)
        zDL_D  = Complex{BigFloat}(0.0, 0.0)
    end

    return zCL,    zDL, 
           zCL_mu, zDL_mu,
           zCL_C,  zDL_C,
           zCL_D,  zDL_D
end

function get_C_D_right(z_mu, zC, zD, w_num_core, r0, a, compute_deriv, tol_bessel, num_terms, factors0; print_flag=false::Bool)
    zr = r0 + a 
    r_base = r0 
    w_num = w_num_core

    # Core solution 10
    zc0 = Complex{BigFloat}(1.0, 0.0)
    zc1 = Complex{BigFloat}(0.0, 0.0)

    f10, df10, f10_mu, df10_mu = bessel_frobenius(zr, z_mu, w_num, r_base, zc0, zc1, compute_deriv, 
                                                  tol_bessel, num_terms, factors0, print_flag=print_flag)

    # Core solution 01
    zc0 = Complex{BigFloat}(0.0, 0.0)
    zc1 = Complex{BigFloat}(1.0, 0.0)

    f01, df01, f01_mu, df01_mu = bessel_frobenius(zr, z_mu, w_num, r_base, zc0, zc1, compute_deriv, 
                                                  tol_bessel, num_terms, factors0, print_flag=print_flag)
    
    # Compute the value of u & du/dx at the mid-left interface as a linear comb. of the 10 & 01 solutions
    aux  = zC * f10  + zD * f01 
    daux = zC * df10 + zD * df01

    zCR = aux 
    zDR = daux * r0 / (r0+a)

    if compute_deriv
        zCR_mu = zC * f10_mu + zD * f01_mu
        zDR_mu = (zC * df10_mu + zD * df01_mu) * r0 / (r0+a)

        zCR_C = f10 
        zDR_C = df10 * r0 / zr

        zCR_D = f01 
        zDR_D = df01 * r0 / zr
    else
        zCR_mu = Complex{BigFloat}(0.0, 0.0)
        zDR_mu = Complex{BigFloat}(0.0, 0.0)
        zCR_C  = Complex{BigFloat}(0.0, 0.0)
        zDR_C  = Complex{BigFloat}(0.0, 0.0)
        zCR_D  = Complex{BigFloat}(0.0, 0.0)
        zDR_D  = Complex{BigFloat}(0.0, 0.0)
    end

    return zCR,    zDR, 
           zCR_mu, zDR_mu,
           zCR_C,  zDR_C,
           zCR_D,  zDR_D
end


function BentThreeLayerWaveguide(zr, z_mu, zC, zD, zCL, zDL, zCR, zDR, w_num_core, w_num_clad, 
        r0, a, b, tol_bessel, num_terms, factors0, factorsL, factorsR; print_flag=false::Bool )
    
    compute_deriv = false

    if real(zr) < r0 - b || real(zr) > r0 + b
        @warn "BentThreeLayerWaveguide: Point is out of interval [r0 - b, r0 + b]. Returning."
        return
    end

    if real(zr) < r0 - a # In the left cladding
        r_base = r0 - a 
        w_num = w_num_clad 
        z_mu_loc = z_mu * r0^2 / (r_base^2)
        factors = factorsL
        zC = zCL
        zD = zDL
    elseif real(zr) <= r0 + a # In the core
        r_base = r0 
        w_num = w_num_core 
        z_mu_loc = z_mu
        factors = factors0
    else # In the right cladding
        r_base = r0 + a 
        w_num = w_num_clad
        z_mu_loc = z_mu * r0^2 / (r_base^2)
        factors = factorsR
        zC = zCR 
        zD = zDR
    end

    # Solution 10
    zc0 = Complex{BigFloat}(1.0, 0.0)
    zc1 = Complex{BigFloat}(0.0, 0.0)

    f10, df10, _, _ = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, compute_deriv, 
                                        tol_bessel, num_terms, factors, print_flag=print_flag)

    # Solution 10
    zc0 = Complex{BigFloat}(0.0, 0.0)
    zc1 = Complex{BigFloat}(1.0, 0.0)

    f01, df01, _, _ = bessel_frobenius(zr, z_mu_loc, w_num, r_base, zc0, zc1, compute_deriv, 
                                        tol_bessel, num_terms, factors, print_flag=print_flag)

    u = zC * f10 + zD * f01
    du = (zC * df10 + zD * df01) * r_base / zr

    return u, du
end