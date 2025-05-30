# `transfermatrix/tmatrixy.f`

## Overview

This Fortran file implements the function `GETTRANSY`, which calculates the transmission properties of a 2D system when the current is flowing in the y-direction. Similar to `tmatrixx.f`, it uses a transfer matrix approach. The system is sliced along the y-direction, and the total transfer matrix is built iteratively. Transmission eigenvalues are obtained via SVD of the transmission block of the final matrix. This file also handles options for wrapping in the x-direction (periodic boundary conditions).

## Key Components

*   **`GETTRANSY(GAUGE, TVALS, LIMX, LIMY, E, FLUX, WRAPX)`**: A `DOUBLE PRECISION` function that returns a value from `CHECKUNI2` (a unitarity check). The primary output, transmission eigenvalues, is via the `TVALS` array.
    *   `GAUGE`: Character, specifies the gauge choice ('X' or 'Y').
    *   `TVALS`: `DOUBLE PRECISION` array (output), stores calculated transmission eigenvalues.
    *   `LIMX`: Integer, dimension of the system in the x-direction (width). This is `NSIZE` for y-current.
    *   `LIMY`: Integer, dimension of the system in the y-direction (length, number of slices).
    *   `E`: `DOUBLE PRECISION`, energy.
    *   `FLUX`: `DOUBLE PRECISION`, magnetic flux.
    *   `WRAPX`: Integer, flag for periodic boundary conditions in x (1 for wrap, 0 for no wrap).

*   **`CALCMULTYX(E, FLUX, POS, WRAPX, MULT, LIMX)`**: Subroutine to calculate the slice transfer matrix `MULT` for y-current when `GAUGE = 'X'`.
    *   `POS`: Integer, current slice index in the y-direction.
    *   `WRAPX`: Integer, x-wrapping flag.
    *   `MULT`: `DOUBLE COMPLEX` array (output), the slice transfer matrix.
    *   `LIMX`: Integer, width of the system.

*   **`CALCMULTYY(E, FLUX, POS, WRAPX, MULT, LIMX)`**: Subroutine to calculate the slice transfer matrix `MULT` for y-current when `GAUGE = 'Y'`.
    *   `POS`: Integer, current slice index.
    *   `WRAPX`: Integer, x-wrapping flag.
    *   `MULT`: `DOUBLE COMPLEX` array (output).
    *   `LIMX`: Integer, width of the system.

*   **`FILLUYX(U, FLUX, LIMX)`**: Subroutine to fill matrix `U` when `GAUGE = 'X'`. `U` seems to be a complex identity matrix in this case.
    *   `U`: `DOUBLE COMPLEX` array (output).
    *   `FLUX`: `DOUBLE PRECISION` (not used here).

*   **`FILLUYY(U, FLUX, LIMX)`**: Subroutine to fill matrix `U` when `GAUGE = 'Y'`. Incorporates Peierls phases.
    *   `U`: `DOUBLE COMPLEX` array (output).
    *   `FLUX`: `DOUBLE PRECISION`.

## Important Variables/Constants

*   `MULT(2*LIMX, 2*LIMX)`: `DOUBLE COMPLEX` array, slice transfer matrix.
*   `A(LIMX, LIMX)`, `B(LIMX, LIMX)`, `C(LIMX, LIMX)`, `D(LIMX, LIMX)`: `DOUBLE COMPLEX` arrays, sub-blocks of `MULT` after transformation with `U`.
*   `T(LIMX, LIMX)`, `R(LIMX, LIMX)`, `TTILDE(LIMX, LIMX)`, `RTILDE(LIMX, LIMX)`: Cumulative T/R matrices.
*   `TINC(LIMX, LIMX)`, etc.: Incremental T/R matrices.
*   `U(LIMX, LIMX)`: Transformation matrix for gauge and field effects.
*   `THR/1e-12/`: A threshold parameter, possibly for `CORRUNI` (commented out).

## Algorithm Outline (within `GETTRANSY`)

1.  Initialize cumulative `T`, `TTILDE` to identity, `R`, `RTILDE` to zero.
2.  Determine `U` matrix based on `GAUGE` using `FILLUYX` or `FILLUYY`.
3.  Loop `I` from 1 to `LIMY` (iterate through y-slices):
    a.  Calculate slice transfer matrix `MULT` using `CALCMULTYX` or `CALCMULTYY`.
    b.  Generate `A, B, C, D` from `MULT` and `U` using `GENABCD`.
    c.  Calculate incremental `TINC, RINC, TTILDEINC, RTILDEINC` using `GENTANDRINC`.
    d.  Update cumulative `T, R, TTILDE, RTILDE` using `UPDATETANDR`.
    e.  (Commented out: `CALL CORRUNI` - a unitarity correction step).
4.  Calculate SVD of final `T` matrix using `SQSVDVALS`, store in `TVALS`.
5.  Return result of `CHECKUNI2`.

## Error Handling

*   The code checks if `LIMX` is even when `WRAPX = 1` (periodic boundaries in x). If `LIMX` is odd, it prints an error and stops. This is likely due to the way periodic boundary conditions are implemented for the chosen lattice discretization, requiring an even number of sites for consistency.

## Usage Examples

This function is part of the `GETTRANS` interface, typically called via the Python wrapper.

## Dependencies and Interactions

*   **`SQUNIT(matrix, size)`**: (External) Initialize to identity.
*   **`SQZERO(matrix, size)`**: (External) Initialize to zero.
*   **`SQUNITZ(matrix, val, size)`**: (External) Initialize to `val` * identity.
*   **`ZPOLAR(angle, complex_num)`**: (External/System) Create complex number from polar form.
*   **`GENABCD(...)`**: (External) Extract A,B,C,D sub-blocks.
*   **`GENTANDRINC(...)`**: (External) Calculate incremental T/R matrices.
*   **`UPDATETANDR(...)`**: (External) Update cumulative T/R matrices.
*   **`SQSVDVALS(matrix, values, size)`**: (External) Compute SVD.
*   **`CHECKUNI2(...)`**: (External) Function for unitarity check.
*   **`CORRUNI(...)`**: (External, but commented out) Subroutine for unitarity correction.

The subroutines `CALCMULTYX` and `CALCMULTYY` define the on-slice Hamiltonian and inter-slice hopping terms, including Peierls phases for the magnetic field (`FLUX`) and handling of boundary conditions (`WRAPX`). The logic within these routines for `NEIGH` (neighbor indexing) and phase factors `ZPLUS`, `ZMINUS`, `CNUM` is specific to the tight-binding model on a square lattice and the chosen gauge. The `POS` variable is used to alternate coupling terms, similar to `tmatrixx.f`.
```
