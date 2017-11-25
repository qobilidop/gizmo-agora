#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/* this does the local wind coupling on a per-star basis, as opposed to in the SFR routine */

/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


