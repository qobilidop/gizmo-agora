#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"

/* Routines for gas equation-of-state terms (collects things like calculation of gas pressure)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* return the pressure of particle i */
double get_pressure(int i)
{
    MyFloat press = GAMMA_MINUS1 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); /* ideal gas EOS (will get over-written it more complex EOS assumed) */
    
#ifdef GALSF_EFFECTIVE_EQS
    /* modify pressure to 'interpolate' between effective EOS and isothermal, with the Springel & Hernquist 2003 'effective' EOS */
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh)
        press = All.FactorForSofterEQS * press + (1 - All.FactorForSofterEQS) * All.cf_afac1 * GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;
#endif
    
    
#ifdef EOS_HELMHOLTZ
    /* pass the necessary quantities to wrappers for the Timms EOS */
    struct eos_input eos_in;
    struct eos_output eos_out;
    eos_in.rho  = SphP[i].Density;
    eos_in.eps  = SphP[i].InternalEnergyPred;
    eos_in.Ye   = SphP[i].Ye;
    eos_in.Abar = SphP[i].Abar;
    eos_in.temp = SphP[i].Temperature;
    int ierr = eos_compute(&eos_in, &eos_out);
    assert(!ierr);
    press              = eos_out.press;
    SphP[i].SoundSpeed = eos_out.csound;
    SphP[i].Temperature= eos_out.temp;
#endif

    
#ifdef EOS_TILLOTSON
    press = calculate_eos_tillotson(i);
#endif
    
#ifdef EOS_ENFORCE_ADIABAT
    press = EOS_ENFORCE_ADIABAT * pow(SphP[i].Density, GAMMA);
#endif
    
    
    
    
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
    /* add an artificial pressure term to suppress fragmentation at/below the explicit resolution scale */
    double h_eff = DMAX(Get_Particle_Size(i), All.ForceSoftening[0]/2.8); /* need to include latter to account for inter-particle spacing << grav soft cases */
    /* standard finite-volume formulation of this (note there is some geometric ambiguity about whether there should be a "pi" in the equation below, but this 
        can be completely folded into the (already arbitrary) definition of NJeans, so we just use the latter parameter */
    double NJeans = 4; // set so that resolution = lambda_Jeans/NJeans -- fragmentation with Jeans/Toomre scales below this will be artificially suppressed now
    double xJeans = (NJeans * NJeans / GAMMA) * All.G * h_eff*h_eff * SphP[i].Density * SphP[i].Density * All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
    SphP[i].SoundSpeed = sqrt(GAMMA * press / Particle_density_for_energy_i(i));
#endif
    
    return press;
}




/* trivial function to check if particle falls below the minimum allowed temperature */
void check_particle_for_temperature_minimum(int i)
{
    if(All.MinEgySpec)
    {
        if(SphP[i].InternalEnergy < All.MinEgySpec)
        {
            SphP[i].InternalEnergy = All.MinEgySpec;
            SphP[i].DtInternalEnergy = 0;
        }
    }
}



double INLINE_FUNC Particle_density_for_energy_i(int i)
{
#ifdef HYDRO_PRESSURE_SPH
    return SphP[i].EgyWtDensity;
#endif
    return SphP[i].Density;
}




double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{
#ifdef EOS_GENERAL
    return SphP[i].SoundSpeed;
#endif
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
}





