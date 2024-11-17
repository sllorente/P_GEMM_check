# test_p?gemm
Program to check the correct functioning of the p?gemm subroutines in the the mkl_scalapack library.

The test compute the following operation using the subroutine p?gemm:

B(1:N2, 1:N2) = C(N1+1:N1+N2, 1:N1)*C(1:N1, N1+1:N1+N2) - B(1:N2, 1:N2)

where all coefficients of matrix C are set to 1, and all coeffs of matrix B are initially set to N1. The expected result is therefore that all coeffs of B are equal to 0 after the computation.


## Compile
```
$ make
```
It is assumed that the intel compilers are installed, and the environment variables are properly set.s

## Run test
```
$ make test_[sdcz]
```
## Results
The test fails for different combinations of the values of the matrix sizes (N1 and N2), the block sizes (mb and nb) and the number of rows and columns in the process grid (nprow and npcol):
  - `mb=4 ; nb=4 ; N1=5*4 ; N2=2**10 + 1 ; nprow=2; npcol=2`
  - `mb=4 ; nb=4 ; N1=5*4 ; N2=2**10 + 1 ; nprow=3; npcol=2`
  - `mb=4 ; nb=4 ; N1=5*4 ; N2=2**9*3 + 1; nprow=3; npcol=4`
  - `mb=4 ; nb=4 ; N1=5*4 ; N2=2**9*4 + 1; nprow=4; npcol=4`
  - `mb=32; nb=32; N1=5*32; N2=2**9*4 + 1; nprow=4; npcol=4`
  - ...

The error occurs for mkl versions 2022.2.0, 2024.0 and 2025.0 (and probably others) when linking with intelmpi, libmkl_blacs_intelmpi_lp64 and libmkl_scalapack_lp64. The error is silent as the code runs and exits normally. 

When linking with a manually compiled version of Scalapack 2.2.0, the subroutine p?gemm returns the expected values and the test program terminate without displaying any erroneous values.

