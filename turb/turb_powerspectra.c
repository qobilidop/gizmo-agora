#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

/*
 *  This code was originally written for GADGET3 by Andreas Bauer; it has been
 *   modified slightly by Phil Hopkins for GIZMO, but is largely intact.
 */

#if defined(TURB_DRIVING_SPECTRUMGRID) && defined(BOX_PERIODIC) && (defined(TURB_DRIVING))


#ifndef USE_FFTW3
#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif
#else /* FFTW3 */
#include "../gravity/myfftw3.h"
#endif

#define  TURB_DRIVING_SPECTRUMGRID2 (2*(TURB_DRIVING_SPECTRUMGRID/2 + 1))

#if (TURB_DRIVING_SPECTRUMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#ifndef USE_FFTW3
static rfftwnd_mpi_plan fft_forward_plan;
static int slabstart_x, nslab_x, slabstart_y, nslab_y;

static int fftsize, maxfftsize;
#else 
//static fftw_plan fft_forward_plan;
static fftw_plan fft_velx_plan, fft_vely_plan, fft_velz_plan;
static fftw_plan fft_svelx_plan, fft_svely_plan, fft_svelz_plan;
static fftw_plan fft_vrhox_plan, fft_vrhoy_plan, fft_vrhoz_plan;
static fftw_plan fft_vortx_plan, fft_vorty_plan, fft_vortz_plan;
static fftw_plan fft_dis1field_plan, fft_dis2field_plan; 
static fftw_plan fft_rand_plan; 
static fftw_plan fft_dens_plan; 

static ptrdiff_t slabstart_x, nslab_x, slabstart_y, nslab_y;

static ptrdiff_t fftsize, maxfftsize;
static MPI_Datatype MPI_TYPE_PTRDIFF; 
#endif
static fftw_real *velfield[3];
#ifdef TURB_DIFF_DYNAMIC
static fftw_real *velbarfield[3];
static fftw_real *velhatfield[3];
#endif
static fftw_real *smoothedvelfield[3];
static fftw_real *vorticityfield[3];
static fftw_real *velrhofield[3];
static fftw_real *dis1field;
static fftw_real *dis2field;
static fftw_real *densityfield;

static fftw_real *randomfield;
static fftw_real *workspace;

static float    *RandomValue;

static fftw_complex *fft_of_field;

static float *powerspec_turb_nearest_distance, *powerspec_turb_nearest_hsml;

#ifndef USE_FFTW3
void powerspec_turb_calc_and_bin_spectrum(fftw_real *field, int flag);
#else 
void powerspec_turb_calc_and_bin_spectrum(fftw_plan plan, fftw_real *field, int flag);
#endif

#ifndef USE_FFTW3
#define cmplx_re(c) ((c).re)
#define cmplx_im(c) ((c).im)
#else 
#define cmplx_re(c) ((c).[0])
#define cmplx_im(c) ((c).[1])
