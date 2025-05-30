# `transfermatrix/tmatrix.f`

## Overview

This Fortran file defines the main function `GETTRANS`, which serves as a primary interface for calculating transmission probabilities through a 2D system. It acts as a dispatcher, calling specific functions (`GETTRANSX` or `GETTRANSY`) based on the specified current direction. This function is crucial for simulating quantum transport properties.

## Key Components

*   **`GETTRANS(CURRENT, GAUGE, TVALS, NVALS, LIMX, LIMY, E, FLUX, WRAPX)`**: A `DOUBLE PRECISION` function that returns a value indicating the success or failure of the transmission calculation (though its primary output is via the `TVALS` array).
    *   `CURRENT`: Character, specifies the direction of current flow ('X' or 'Y').
    *   `GAUGE`: Character, specifies the gauge choice for the magnetic field (e.g., Landau gauge).
    *   `TVALS`: `DOUBLE PRECISION` array (output), stores the calculated transmission values (eigenvalues of t't matrix).
    *   `NVALS`: Integer (output), the number of transmission values written to `TVALS`.
    *   `LIMX`, `LIMY`: Integers, dimensions of the system in x and y directions.
    *   `E`: `DOUBLE PRECISION`, energy at which the transmission is calculated.
    *   `FLUX`: `DOUBLE PRECISION`, magnetic flux through the system.
    *   `WRAPX`: Integer, flag indicating whether periodic boundary conditions are applied in the x-direction (1 for wrap, 0 for no wrap). Relevant for Y-current calculations.

## Important Variables/Constants

*   **`MAXSIZE`**: Integer parameter, defines the maximum size of the `TVALS` array (currently 10000). This limits the maximum number of transmission channels that can be handled.
*   **`GETTRANS` (return value)**: Initialized to -1. If `NVALS` remains 0 after checking `CURRENT` direction (i.e., invalid `CURRENT` identifier), the function will effectively return -1 (though this value isn't explicitly used by callers shown in `hofstadter.py` which rely on `NVALS` or `TVALS`). The actual transmission calculation results are passed via the `TVALS` array.

## Usage Examples

This function is typically called from other routines (like the Python wrapper `_tmatrix.c`) that set up the system parameters and then pass them to `GETTRANS` to obtain the transmission eigenvalues.

```fortran
! Example (conceptual)
IMPLICIT NONE
CHARACTER CURRENT_DIR, GAUGE_TYPE
PARAMETER (CURRENT_DIR = 'X', GAUGE_TYPE = 'L') ! Example: Landau gauge for X current
INTEGER, PARAMETER :: SYS_X = 100, SYS_Y = 50
DOUBLE PRECISION, PARAMETER :: ENERGY = 0.5, MAG_FLUX = 0.1
INTEGER, PARAMETER :: MAX_CHANNELS = SYS_Y / 2
DOUBLE PRECISION TRANSMISSION_EIGENVALUES(MAX_CHANNELS)
INTEGER NUM_EIGENVALUES
DOUBLE PRECISION RESULT_STATUS
INTEGER WRAP_FLAG_X

! For X current, LIMY must be even. WRAPX is not used by GETTRANSX.
! For Y current, NVALS will be LIMX.
IF (CURRENT_DIR .EQ. 'X') THEN
    WRAP_FLAG_X = 0 ! Not directly used by GETTRANSX
    NUM_EIGENVALUES = SYS_Y / 2
ELSE IF (CURRENT_DIR .EQ. 'Y') THEN
    WRAP_FLAG_X = 1 ! Example: wrap in x for y-current
    NUM_EIGENVALUES = SYS_X
ELSE
    WRITE(*,*) "Invalid current direction for example"
    STOP
END IF


RESULT_STATUS = GETTRANS(CURRENT_DIR, GAUGE_TYPE,
     +                     TRANSMISSION_EIGENVALUES, NUM_EIGENVALUES,
     +                     SYS_X, SYS_Y,
     +                     ENERGY, MAG_FLUX, WRAP_FLAG_X)

IF (NUM_EIGENVALUES .GT. 0) THEN
    PRINT *, "Calculated transmission eigenvalues:"
    PRINT *, (TRANSMISSION_EIGENVALUES(I), I=1,NUM_EIGENVALUES)
ELSE
    PRINT *, "Failed to get transmission values or invalid setup."
ENDIF
```

## Dependencies and Interactions

*   **`GETTRANSX`**: External `DOUBLE PRECISION` function. Called when `CURRENT = 'X'`. Expected to be defined in `tmatrixx.f`. This function calculates transmission for current in the x-direction.
*   **`GETTRANSY`**: External `DOUBLE PRECISION` function. Called when `CURRENT = 'Y'`. Expected to be defined in `tmatrixy.f`. This function calculates transmission for current in the y-direction.
*   The file expects `LIMY` to be even if `CURRENT = 'X'`.
*   The results (transmission eigenvalues) are stored in the `TVALS` array, and `NVALS` indicates how many such values were computed.
```
