# Arbitrary Precision Bessel Function and Waveguide Solvers

In this repository we provide:
 1) `bessel_frobenius`: A generalized Bessel function implementation for evaluating a Bessel function (i.e. a solution to the Bessel differential equation with initial value `zc0` and initial derivative `zc1` at point `r0`) of arbitrary complex order `sqrt(mu)` at (possibly complex) point `x` at arbitrary precision and
 2) Non-linear solvers for computing the first three propagating modes of the Bessel eigenvalue problem arising in the study of loss factors in three-layer optical slab waveguides subject to bending. 
 
 These solvers were used to produce the results in the paper, "Bessel Functions and Analysis of Circular Waveguides" by J Mora-Paz, L. Demkowicz, C.G. Taylor, J. Grosek, and S. Henneking, a preprint of which is currently available on ArXiv at: [https://arxiv.org/abs/2512.04348](https://arxiv.org/abs/2512.04348).

## Usage
### Dependencies
The only non-stdlib dependencies of this code is the [Plots.jl](https://docs.juliaplots.org/stable/) library for plotting. The [Manifest.toml](https://github.com/cgt3/WaveguideSolver/blob/main/Manifest.toml) and [Project.toml](https://github.com/cgt3/WaveguideSolver/blob/main/Project.toml) files can be used to replicate the environment used by the authors to generate the results in the paper.

### Setting Arbitrary Precision
Within these codes we make use of Julia's `BigFloat` arbitrary precision data type (see Julia documentation [here](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic)) and its complex equivalent/parameterized type `Complex{BigFloat}`. 

The use of ```BigFloat``` is a necessity for these computations for three reasons: 
 1) the exponentially/factorially fast decay of the coefficients used in the Frobenius expansion in ```bessel_frobenius```, 
 2) the discrepancy in the magnitudes of the real and imaginary components (around 30 orders of magnitude difference in some cases), and 
 3) the extreme sensitivity of the imaginary component of the modes.

The precision needed depends on the magnitude of the input to `bessel_frobenius`. Within the given example codes we use a precision of 70 base-10 digits, but more digits may be necessary for other parameter values.

The precision of Julia's `BigFloat` type is set globally using
```
setprecision(n, base=a)
```
where `n` is the number of digits desired in base `a`. Double and quadruple precision can be duplicated using:
```
setprecision(53,  base=2) # Double precision, approx. 16-17 base-10 digits of accuracy
setprecision(113, base=2) # Quadruple precision, approx 33 base-10 digits of accuracy
```

### Reproducing Paper Results
The files [Table5.jl](https://github.com/cgt3/WaveguideSolver/blob/main/Table5.jl), [Table6-Figure4.jl](https://github.com/cgt3/WaveguideSolver/blob/main/Table6-Figure4.jl), [Table7-Figure6.jl](https://github.com/cgt3/WaveguideSolver/blob/main/Table7-Figure6.jl), and [Table8-Figure7.jl](https://github.com/cgt3/WaveguideSolver/blob/main/Table8-Figure7.jl) will produce the results shown in the denoted entity in the paper. Note that for the figures, the images in the paper were produced in MatLab but are equivalent.

### Testing Other Setups
The file [driver_solver.jl](https://github.com/cgt3/WaveguideSolver/blob/main/driver_solver.jl) can be used to test other choices of parameters, precision, and modes. A guide to the problem parameters is given below.

## Parameters
As discussed in the accompanying paper, the following problem parameters appear in the code.

### Geometry Parameters:
- `r0` the radius of curvature of the bent waveguide
- `a` the distance from the fiber's centerline to the edge of the core
- `b` the distance from the fiber's centerline to the edge of the cladding
- `ref_length` the length scale for the problem of interest (which matches the diameter of the waveguide core).

### Material Parameters:
- `k0` the freespace wavelumber (in non-dimensional units)
- `n_core` the core's refractive index (n_0 in the paper)
- `n_clad` the cladding's refractive index (n_1 in the paper)
- `wavelength` the wavelength of light being studied (in non-dimensional units)

### Boundary Conditions:
`BC = (;type, param)` sets the type of BC to use and any parameters needed for its enforcement (e.g., the PML constant). 

Valid values are: 
 - `BC=(; type="imp")` for impedance BCs
- `BC=(; type="pml", PML_const)` for PML BC with constant `PML_const`

### Modes:
There are two solvers, one for even modes and another for odd modes. Which even/odd mode is recovered depends on the initial guess for `z_mu` (mu in the paper). The initial guesses for recovering the first even, first odd, and second even mode are provided as the constants `EVEN1`, `ODD1`, and `EVEN2` in [bessel_solver.jl](https://github.com/cgt3/WaveguideSolver/blob/main/bessel_solver.jl). 

Other modes can be found by providing a `NamedTuple` with fields:
- `type` valid options are the strings `"even"` or `"odd"`
- `index` this is only used for printing purposes and does not affect the solver in any way
- `z_mu_guess` the initial guess for mu. This decides which even/odd mode is recovered.

### Bessel Function Solver (`bessel_frobenius` ) Parameters:
These parameters should be changed whenever the precision of `BigFloat` is changed.
- `tol_bessel` the tolerance for the Bessel function expansion
- `num_factors` the maximum number of terms to use in the Frobenius/power series expansion of complex ordered Bessel functions

### Non-Linear Solver Parameters:
These parameters should be changed whenever the precision of `BigFloat` is changed.
- `tol_eigen` the tolerance for the non-linear solver for finding the parameters of modes
- `max_itr` the maximum number of iterations for the non-linear solver for finding the paramters of modes

 ## Authors
 This code was developed by Jaime Mora-Paz, Christina G. Taylor, and Leszek Demkowicz for the paper "Bessel Functions and Analysis of Circular Waveguides". If this code is of use to you, please cite this paper.

 ```
 @misc{mora2025bessel-arxiv,
      title={{B}essel {F}unctions and {A}nalysis of {C}ircular {W}aveguides}, 
      author={Jaime Mora-Paz and Leszek Demkowicz and Christina G. Taylor and Jacob Grosek and Stefan Henneking},
      year={2025},
      eprint={2512.04348},
      archivePrefix={arXiv},
      primaryClass={physics.comp-ph},
      url={https://arxiv.org/abs/2512.04348}, 
      doi={10.48550/arXiv.2512.04348}
}
 ```






