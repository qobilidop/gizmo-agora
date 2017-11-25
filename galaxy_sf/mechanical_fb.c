#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc */

/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */



