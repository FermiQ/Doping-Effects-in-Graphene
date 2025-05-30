# `transfermatrix/linalg.f`

## Overview

This Fortran file provides a collection of linear algebra routines for operations on square matrices of double complex numbers. These routines are wrappers around standard BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra Package) routines, specifically their complex versions (prefixed with `Z`). They are used extensively in the transfer matrix calculations for matrix multiplication, inversion, solving linear systems, and singular value decomposition (SVD).

## Key Components

All subroutines operate on square matrices of size `N x N`.

*   **`SQDOT(C, A, B, N)`**: Calculates matrix product `C = A * B`.
    *   Wraps `ZGEMM`.
*   **`SQDOTNH(C, A, B, N)`**: Calculates `C = A * B^H` (B Hermitian conjugate).
    *   Wraps `ZGEMM` with `TRANSA = 'N'`, `TRANSB = 'C'`.
*   **`SQDOTHN(C, A, B, N)`**: Calculates `C = A^H * B` (A Hermitian conjugate).
    *   Wraps `ZGEMM` with `TRANSA = 'C'`, `TRANSB = 'N'`.
*   **`SQUPDAXPY(Y, ALPHA, X, N)`**: Updates matrix `Y` as `Y = Y + ALPHA * X`.
    *   Wraps `ZAXPY` (treats matrices as vectors).
*   **`SQDOTAX(C, ALPHA, A, B, N)`**: Calculates `C = ALPHA * A * B`.
    *   Wraps `ZGEMM`.
*   **`SQDOTUPD(BETA, C, ALPHA, A, B, N)`**: Calculates `C = BETA * C + ALPHA * A * B`.
    *   Wraps `ZGEMM`.
*   **`INVERTMATRIX(MATRIX, N)`**: Inverts `MATRIX` in place.
    *   Uses LAPACK routines `ZGETRF` (LU decomposition) and `ZGETRI` (matrix inversion from LU).
    *   Stops execution if the matrix is non-invertible.
*   **`DNORMDIFF(A, B, N)`**: `DOUBLE PRECISION` function that calculates the Frobenius norm of the difference between matrices A and B, i.e., `||A - B||_F`.
    *   Uses `SQCOPY`, `SQUPDAXPY`, and LAPACK's `ZLANGE` (to compute matrix norm).
*   **`INVSOLVE(A, B, N)`**: Solves the system `A * X = B` for `X` and stores the result `X` back into `A`. So, on output, `A` becomes `inv(A_original) * B`.
    *   Uses LAPACK routines `ZGETRF` (LU decomposition of `A`) and `ZGETRS` (solves system using LU factors).
    *   Stops execution if `A` is singular or the system is non-solvable.
*   **`SQSVDVALS(MATRIX, OUTPUTS, N)`**: Computes the singular values of `MATRIX`.
    *   `OUTPUTS`: `DOUBLE PRECISION` array (output) to store the singular values.
    *   Wraps LAPACK routine `ZGESVD`. `U` and `V` matrices are not computed ('N' option).
    *   Stops execution if SVD fails.
*   **`SQSVDFULL(MATRIX, OUTPUTS, U, V, N)`**: Computes the full SVD: `MATRIX = U * S * V^H`.
    *   `OUTPUTS`: `DOUBLE PRECISION` array (output) for singular values `S`.
    *   `U`: `DOUBLE COMPLEX` array (output), left singular vectors.
    *   `V`: `DOUBLE COMPLEX` array (output), right singular vectors (Hermitian conjugate is taken by ZGESVD convention).
    *   Wraps LAPACK routine `ZGESVD`.
    *   Stops execution if SVD fails.

## Important Variables/Constants

*   `ZEROC /0.0/`, `ONEC /1.0/`: `DOUBLE COMPLEX` parameters for zero and one.
*   `MINUSONE /-1.0/`: `DOUBLE COMPLEX` parameter for -1.0, used in `DNORMDIFF`.
*   `S`: Integer, status code returned by LAPACK routines. A non-zero value indicates an error, and most routines here will stop execution if `S /= 0`.
*   `PIVOT(N,N)`: Integer array used by LU decomposition routines.
*   `WORK(...)`, `RWORK(...)`: `DOUBLE COMPLEX` and `DOUBLE PRECISION` workspace arrays required by various LAPACK routines.

## Usage Examples

These routines are building blocks for higher-level functions in the transfer matrix code, such as those in `tmatrixx.f`, `tmatrixy.f`, and `trupdate.f`. For instance, `INVERTMATRIX` is used in constructing parts of the slice transfer matrices, and `SQSVDVALS` is crucial for extracting transmission eigenvalues.

```fortran
! Conceptual example of using SQSVDVALS
INTEGER, PARAMETER :: N_SIZE = 4
DOUBLE COMPLEX MY_MATRIX(N_SIZE, N_SIZE)
DOUBLE PRECISION SINGULAR_VALUES(N_SIZE)

! ... (Initialize MY_MATRIX with some values) ...

CALL SQSVDVALS(MY_MATRIX, SINGULAR_VALUES, N_SIZE)

PRINT *, "Singular Values:"
PRINT *, (SINGULAR_VALUES(I), I=1, N_SIZE)
```

## Dependencies and Interactions

*   **BLAS**: Basic Linear Algebra Subprograms library (specifically complex double precision versions like `ZGEMM`, `ZAXPY`).
*   **LAPACK**: Linear Algebra Package library (specifically complex double precision versions like `ZGETRF`, `ZGETRI`, `ZGETRS`, `ZGESVD`, `ZLANGE`).
*   **`SQCOPY(FROM, TO, N)`**: Assumed to be an external subroutine (likely in `util.f`) that copies matrix `FROM` to matrix `TO`.

This file is critical for the numerical stability and performance of the entire simulation, leveraging optimized standard libraries for computationally intensive tasks.
```
