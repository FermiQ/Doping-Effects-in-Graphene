# `transfermatrix/tmatrixx.f`

## Overview

This Fortran file implements the function `GETTRANSX`, which calculates the transmission properties of a 2D system when the current is flowing in the x-direction. It employs a transfer matrix approach, iteratively building up the total transfer matrix by multiplying individual slice matrices. The final transmission eigenvalues are obtained via Singular Value Decomposition (SVD) of the transmission block of the total matrix.

## Key Components

*   **`GETTRANSX(GAUGE, TVALS, LIMX, NSIZE, E, FLUX)`**: A `DOUBLE PRECISION` function that returns a value from `CHECKUNI2` (presumably a check of unitarity). The primary output, the transmission eigenvalues, is through the `TVALS` array.
    *   `GAUGE`: Character, specifies the gauge choice ('X' or 'Y'). This determines how the magnetic field (Peierls phase) is incorporated into the matrices.
    *   `TVALS`: `DOUBLE PRECISION` array (output), stores the calculated transmission eigenvalues (singular values of the transmission matrix T).
    *   `LIMX`: Integer, the number of slices in the x-direction (length of the system).
    *   `NSIZE`: Integer, half of the system's width in the y-direction (`LIMY/2`). This corresponds to the number of propagating modes or channels at one end.
    *   `E`: `DOUBLE PRECISION`, energy at which the transmission is calculated.
    *   `FLUX`: `DOUBLE PRECISION`, magnetic flux.

*   **`CALCMULTXY(E, FLUX, POS, MULT, NSIZE)`**: Subroutine to calculate the transfer matrix `MULT` for a single slice when the gauge is chosen along the y-direction (`GAUGE = 'Y'`).
    *   `POS`: Integer, current position (slice index) in the x-direction.
    *   `MULT`: `DOUBLE COMPLEX` array (output), the slice transfer matrix.

*   **`CALCMULTXX(E, FLUX, POS, MULT, NSIZE)`**: Subroutine to calculate the transfer matrix `MULT` for a single slice when the gauge is chosen along the x-direction (`GAUGE = 'X'`).
    *   `POS`: Integer, current position (slice index) in the x-direction.
    *   `MULT`: `DOUBLE COMPLEX` array (output), the slice transfer matrix.

*   **`FILLUXY(U, FLUX, LIMX)`**: Subroutine to fill the matrix `U` used in `GENABCD` when `GAUGE = 'Y'`. `U` relates wavefunctions across a slice boundary. For the 'Y' gauge, this seems to be a complex identity matrix. (Note: `LIMX` parameter here is likely a typo and should be `NSIZE` as `U` is `NSIZE x NSIZE`).
    *   `U`: `DOUBLE COMPLEX` array (output).
    *   `FLUX`: `DOUBLE PRECISION`, magnetic flux (not used in this specific 'Y' gauge implementation).

*   **`FILLUXX(U, FLUX, LIMX)`**: Subroutine to fill the matrix `U` when `GAUGE = 'X'`. This matrix incorporates the Peierls phase factors due to the magnetic field. (Note: `LIMX` parameter here is likely a typo and should be `NSIZE`).
    *   `U`: `DOUBLE COMPLEX` array (output).
    *   `FLUX`: `DOUBLE PRECISION`, magnetic flux.

## Important Variables/Constants

*   `MULT(2*NSIZE, 2*NSIZE)`: `DOUBLE COMPLEX` array, stores the transfer matrix for a single slice of the system.
*   `A(NSIZE, NSIZE)`, `B(NSIZE, NSIZE)`, `C(NSIZE, NSIZE)`, `D(NSIZE, NSIZE)`: `DOUBLE COMPLEX` arrays, sub-blocks of the slice transfer matrix `MULT` after transformation with `U`.
*   `T(NSIZE, NSIZE)`, `R(NSIZE, NSIZE)`, `TTILDE(NSIZE, NSIZE)`, `RTILDE(NSIZE, NSIZE)`: `DOUBLE COMPLEX` arrays, representing the cumulative transmission and reflection matrices for the system built up so far.
*   `TINC(NSIZE, NSIZE)`, `RINC(NSIZE, NSIZE)`, `TTILDEINC(NSIZE, NSIZE)`, `RTILDEINC(NSIZE, NSIZE)`: `DOUBLE COMPLEX` arrays, representing the incremental transmission and reflection matrices for the current slice being added.
*   `U(NSIZE, NSIZE)`: `DOUBLE COMPLEX` array, transformation matrix incorporating gauge choice and magnetic field effects.

## Algorithm Outline (within `GETTRANSX`)

1.  Initialize cumulative transmission `T`, `TTILDE` to identity matrices and reflection `R`, `RTILDE` to zero matrices.
2.  Determine the `U` matrix based on the `GAUGE` parameter using `FILLUXX` or `FILLUXY`.
3.  Loop `I` from 1 to `LIMX` (iterate through slices):
    a.  Calculate the slice transfer matrix `MULT` for the current slice `I` using `CALCMULTXX` or `CALCMULTXY` based on `GAUGE`.
    b.  Generate the `A, B, C, D` sub-blocks from `MULT` and `U` using `GENABCD`.
    c.  Calculate the incremental `TINC, RINC, TTILDEINC, RTILDEINC` for the current slice using `GENTANDRINC`.
    d.  Update the cumulative `T, R, TTILDE, RTILDE` by combining them with the incremental matrices using `UPDATETANDR`.
4.  After iterating through all slices, calculate the singular values of the final cumulative transmission matrix `T` using `SQSVDVALS`. These singular values are the transmission eigenvalues and are stored in `TVALS`.
5.  Return the result of `CHECKUNI2` (a unitarity check on the final T, R, Ttilde, Rtilde matrices).

## Usage Examples

This function is not meant to be called directly by an end-user but is a core component of the `GETTRANS` interface, which is then wrapped by `_tmatrix.c` for Python usage.

## Dependencies and Interactions

*   **`SQUNIT(matrix, size)`**: Subroutine (external, likely in `util.f` or `linalg.f`) to initialize `matrix` as an identity matrix.
*   **`SQZERO(matrix, size)`**: Subroutine (external, likely in `util.f` or `linalg.f`) to initialize `matrix` as a zero matrix.
*   **`SQUNITZ(matrix, val, size)`**: Subroutine (external, likely in `util.f` or `linalg.f`) to initialize `matrix` as `val` times an identity matrix.
*   **`INVERTMATRIX(matrix, size)`**: Subroutine (external, likely in `linalg.f`) to invert `matrix`.
*   **`GENABCD(MULT, U, A, B, C, D, NSIZE)`**: Subroutine (external, likely in `trupdate.f` or a related file) to extract sub-blocks A, B, C, D from MULT and U.
*   **`GENTANDRINC(A, B, C, D, TINC, RINC, TTILDEINC, RTILDEINC, NSIZE)`**: Subroutine (external, likely in `trupdate.f` or a related file) to calculate incremental T/R matrices.
*   **`UPDATETANDR(T, R, TTILDE, RTILDE, TINC, RINC, TTILDEINC, RTILDEINC, NSIZE)`**: Subroutine (external, likely in `trupdate.f` or a related file) to update cumulative T/R matrices.
*   **`SQSVDVALS(matrix, values, size)`**: Subroutine (external, likely in `linalg.f`) to compute singular values of `matrix` and store them in `values`.
*   **`CHECKUNI2(T, R, TTILDE, RTILDE, NSIZE)`**: `DOUBLE PRECISION` Function (external, likely in `unitarity.f`) to check the unitarity condition of the T/R matrices.
*   **`ZPOLAR(angle, complex_num)`**: Subroutine (external, system or library intrinsic if not in local files) to create a complex number `complex_num` from polar coordinates (magnitude 1, angle `angle`).

The file assumes that `NSIZE` corresponds to `LIMY/2` where `LIMY` is the y-dimension of the system, and that `LIMY` would be even.
The `CALCMULTXX` and `CALCMULTXY` subroutines construct the slice transfer matrix based on a tight-binding Hamiltonian model, incorporating energy `E` and magnetic flux `FLUX` (via Peierls substitution).
The specific structure of `N2` and `N3` matrices within `CALCMULTXX`/`CALCMULTXY` suggests a nearest-neighbor hopping model on a square lattice.
The `POS` variable in `CALCMULTXX`/`CALCMULTXY` is used to alternate the connectivity pattern for odd/even slices, which is characteristic of some lattice discretizations or models.
```
