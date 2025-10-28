# Arbitrary Precision Bessel Function and Waveguide Solvers

In this repository we provide:
 1) `bessel_frobenius`: A generalized Bessel function implementation for evaluating a Bessel function, `J_mu(x)` of arbitrary complex order `mu` at (possibly complex) point `x` at arbitrary precision and
 2) Non-linear solvers for computing the first three modes of the Bessel eigenvalue problem arising in the study of loss factors in three-layer optical slab waveguides. 
 
 These solvers were used to produce the results in the paper, "Bessel Functions and Analysis of Circular Waveguides" by J Mora-Paz, L. Demkowicz, C.G. Taylor, J. Grosek, and S. Henneking.

## Usage

### Setting Arbitrary Precision
Within these codes we make use of Julia's `BigFloat` arbitrary precision data type (see Julia documentation [here](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic)) and its complex equivalent/parameterized type `Complex{BigFloat}`. 

The use of ```BigFloat``` is a necessity for these computations for three reasons: 
 1) the exponetially/factorially fast decay of the coefficients used in Frobenius expansion used in ```bessel_frobenius```, 
 2) the discrepancy in the magnitudes of the real and imaginary components (up to 30 orders of magnitude difference), and 
 3) the extreme sensistivity of the imaginary component of the modes.

The precision needed depends on the magnitude of the input to `bessel_frobenius`. Within the given example codes we use a precision of 70 base-10 digits. 

To globally set the precision of the `BigFloat` data type use
```
setprecision(n, base=a)
```
where `n` is the number of digits desired in base `a`. To duplicate double and quadruple precision by using the following:
```
setprecision(53,  base=2) # Double precision, approx. 16-17 base-10 digits of accuracy
setprecision(113, base=2) # Quadruple precision, approx 33 base-10 digits of accuracy
```
 
 

 ### Complex-Ordered Bessel Function

 ### Reproducing Paper Results






