#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

/* Routines for simple photo-ionization heating feedback model */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#if defined(FLAG_NOT_IN_PUBLIC_CODE) || (defined(FLAG_NOT_IN_PUBLIC_CODE) && defined(GALSF))

double particle_ionizing_luminosity_in_cgs(long i)
{
    double lm_ssp=0;
    if(P[i].Type != 5)
    {
        /* use updated SB99 tracks: including rotation, new mass-loss tracks, etc. */
        double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age >= 0.02) return 0; // skip since old stars don't contribute
        if(star_age < 0.0035) {lm_ssp=500.;} else {
            double log_age=log10(star_age/0.0035);
            lm_ssp=470.*pow(10.,-2.24*log_age-4.2*log_age*log_age) + 60.*pow(10.,-3.6*log_age);
        }
        lm_ssp *= calculate_relative_light_to_mass_ratio_from_imf(i);
#ifdef SINGLE_STAR_FORMATION
        /* use effective temperature as a function of stellar mass and size to get ionizing photon production */
        double l_sol = bh_lum_bol(0,P[i].Mass,i) * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)); // L/Lsun
        double m_sol = P[i].Mass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS)); // M/Msun
        double r_sol = pow(m_sol, 0.738); // R/Rsun
        double T_eff = 5780. * pow(l_sol/(r_sol*r_sol), 0.25); // ZAMS effective temperature
        double x0 = 157800./T_eff; // h*nu/kT for nu>13.6 eV
        double fion = 0.0; // fraction of blackbody emitted above x0
        if(x0 < 30.) {double q=18./(x0*x0) + 1./(8. + x0 + 20.*exp(-x0/10.)); fion = exp(-1./q);} // accurate to <10% for a Planck spectrum to x0>30, well into vanishing flux //
        lm_ssp = fion * l_sol / m_sol; // just needs to be multiplied by the actual stellar luminosity to get luminosity to mass ratio
#endif
    } // (P[i].Type != 5)
    
    
    lm_ssp *= (1.95*P[i].Mass*All.UnitMass_in_g/All.HubbleParam); // convert to luminosity from L/M
    if(lm_ssp <= 0) {lm_ssp=0;} // trap for negative values (shouldnt happen)
    if(isnan(lm_ssp)) {lm_ssp=0;} // trap for nans (if stellar age routine cant evaluate non-zero value)
    return lm_ssp;
}


#endif



