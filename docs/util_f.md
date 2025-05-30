# `transfermatrix/util.f`

## Overview

This Fortran file provides a collection of utility subroutines and functions used throughout the transfer matrix project. These include matrix initialization, copying, printing routines for debugging, a function to convert polar coordinates to complex numbers, and a function to calculate conductance from transmission eigenvalues.

## Key Components

*   **`SQCOPY(A, B, LIMX)`**: Copies a square complex matrix `A` to `B`.
    *   `A(LIMX, LIMX)`: Source matrix.
    *   `B(LIMX, LIMX)`: Destination matrix.
    *   `LIMX`: Dimension of the square matrices.
    *   Wraps BLAS routine `ZCOPY`, treating matrices as flat arrays.

*   **`SQZERO(M, LIMX)`**: Initializes a square complex matrix `M` to all zeros.
    *   `M(LIMX, LIMX)`: Matrix to be zeroed.
    *   `LIMX`: Dimension.
    *   Uses LAPACK routine `ZLASET`.

*   **`SQUNIT(M, LIMX)`**: Initializes a square complex matrix `M` to an identity matrix.
    *   `M(LIMX, LIMX)`: Matrix to be initialized.
    *   `LIMX`: Dimension.
    *   Uses LAPACK routine `ZLASET` (sets diagonal to 1.0, off-diagonal to 0.0).

*   **`SQUNITZ(M, ALPHA, LIMX)`**: Initializes a square complex matrix `M` to `ALPHA` times an identity matrix (i.e., `M(i,i) = ALPHA`, `M(i,j) = 0` for `i != j`).
    *   `M(LIMX, LIMX)`: Matrix to be initialized.
    *   `ALPHA`: `DOUBLE COMPLEX` scalar.
    *   `LIMX`: Dimension.
    *   Uses LAPACK routine `ZLASET`.

*   **`PRINTVECTOR(INPUT, LIMX, MNAME)`**: Prints the square of elements of a `DOUBLE PRECISION` vector `INPUT` to standard output.
    *   `INPUT(LIMX)`: Vector to print.
    *   `LIMX`: Dimension of the vector.
    *   `MNAME`: Character string (length 2) used as a label in the output.
    *   Formats output using `ES15.5E2`.

*   **`PRINTT(T, LIMX, MNAME)`**: Prints a square complex matrix `T`, separating real and imaginary parts.
    *   `T(LIMX, LIMX)`: Matrix to print.
    *   `LIMX`: Dimension.
    *   `MNAME`: Character string (length 3) as a label.
    *   Formats output using `ES15.5E3`.

*   **`PRINTM(M, LIMX, MNAME)`**: Prints the real part of a `2*LIMX` by `2*LIMX` square complex matrix `M`.
    *   `M(2*LIMX, 2*LIMX)`: Matrix to print.
    *   `LIMX`: Base dimension (matrix is `2*LIMX x 2*LIMX`).
    *   `MNAME`: Character string (length 3) as a label.
    *   Formats output using `F6.2`.

*   **`ZPRINTM(M, LIMX, MNAME)`**: Prints a square complex matrix `M` in the format `real + imagI`.
    *   `M(LIMX, LIMX)`: Matrix to print.
    *   `LIMX`: Dimension.
    *   `MNAME`: Character string (length 3) as a label.
    *   Formats output using `ES15.5E4`.

*   **`ZPOLAR(ARG, ZCOMPLEX)`**: Converts a number in polar form (magnitude 1, angle `ARG`) to a `DOUBLE COMPLEX` number.
    *   `ARG`: `DOUBLE PRECISION`, the argument (angle in radians).
    *   `ZCOMPLEX`: `DOUBLE COMPLEX` (output), the resulting complex number `cos(ARG) + i*sin(ARG)`.
    *   Uses intrinsic functions `DCOS`, `DSIN`, `DCMPLX`.

*   **`CONDUCTANCE(TVALS, LIMX)`**: `DOUBLE PRECISION` function that calculates the conductance as the sum of squares of transmission eigenvalues.
    *   `TVALS(LIMX)`: `DOUBLE PRECISION` array of transmission eigenvalues (typically, these are singular values of the transmission matrix, so `TVALS(i)*TVALS(i)` gives `T_i`, the transmission probability for channel `i`).
    *   `LIMX`: Number of transmission eigenvalues.
    *   Uses BLAS routine `DDOT` to compute the sum of squares `TVALS(i)^2`. This is equivalent to `sum(T_i)`.
    *   Note: The Landauer formula for conductance is `G = (e^2/h) * Sum(T_n)`. This function calculates `Sum(T_n)`. The prefactor `e^2/h` (the quantum of conductance) is often set to 1 in theoretical calculations.

## Important Variables/Constants

*   `ZERO /0.0/`, `ONE /1.0/`: `DOUBLE COMPLEX` parameters for zero and one used in initialization routines.

## Usage Examples

These utilities are called from various parts of the project. `SQCOPY`, `SQZERO`, `SQUNIT` are used for matrix setup. Printing routines are for debugging. `ZPOLAR` is essential for incorporating Peierls phases from magnetic fields. `CONDUCTANCE` is used at a higher level (e.g., in Python scripts) to get a physical quantity from the calculated transmission eigenvalues.

```fortran
! Example of using ZPOLAR
DOUBLE PRECISION PHASE_ANGLE
DOUBLE COMPLEX PHASE_FACTOR
PHASE_ANGLE = 1.57079632679 ! pi/2

CALL ZPOLAR(PHASE_ANGLE, PHASE_FACTOR)
! PHASE_FACTOR will be approximately (0.0 + 1.0i)

! Example of using CONDUCTANCE
INTEGER, PARAMETER :: N_CHANNELS = 2
DOUBLE PRECISION T_EIGENVALUES(N_CHANNELS)
DOUBLE PRECISION G_VALUE

T_EIGENVALUES(1) = 0.8 ! Corresponds to T_1 = 0.8*0.8 = 0.64
T_EIGENVALUES(2) = 0.6 ! Corresponds to T_2 = 0.6*0.6 = 0.36
! If TVALS are singular values, then TVALS(i)^2 is the transmission probability.
! If TVALS are already transmission probabilities, then DDOT is not the right way.
! Given the context of SVD (SQSVDVALS in linalg.f returns singular values),
! TVALS are likely singular values.

G_VALUE = CONDUCTANCE(T_EIGENVALUES, N_CHANNELS)
! G_VALUE will be 0.64 + 0.36 = 1.0
```

## Dependencies and Interactions

*   **BLAS**: `ZCOPY` (for `SQCOPY`), `DDOT` (for `CONDUCTANCE`).
*   **LAPACK**: `ZLASET` (for `SQZERO`, `SQUNIT`, `SQUNITZ`).
*   Fortran intrinsic functions: `DCOS`, `DSIN`, `DCMPLX`, `REAL`, `AIMAG`, `DIMAG`.

This file centralizes common low-level operations and conversions, promoting code reuse and maintainability.
```
