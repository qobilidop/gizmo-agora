#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"

#ifdef COOL_GRACKLE
#include <grackle.h>
#define ENDRUNVAL 91234

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return pressure
//     4 == calculate and return gamma (only valid when COOL_GRACKLE_CHEMISTRY>0)
//
double CallGrackle(double u_old, double rho, double dt, double ne_guess, double udot, int target, int mode)
{
    gr_float returnval = 0.0;

    // Create struct for storing grackle field data
    grackle_field_data my_fields;

    // Set grid dimension and size.
    // grid_start and grid_end are used to ignore ghost zones.
    my_fields.grid_rank = 3;
    my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
    my_fields.grid_start = malloc(my_fields.grid_rank * sizeof(int));
    my_fields.grid_end = malloc(my_fields.grid_rank * sizeof(int));
    for (int i = 0; i < my_fields.grid_rank; i++) {
        my_fields.grid_dimension[i] = 1;
        my_fields.grid_start[i] = 0;
        my_fields.grid_end[i] = 0;
    }
    my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

    my_fields.density          = malloc(sizeof(gr_float));
    my_fields.internal_energy  = malloc(sizeof(gr_float));
    my_fields.x_velocity       = malloc(sizeof(gr_float));
    my_fields.y_velocity       = malloc(sizeof(gr_float));
    my_fields.z_velocity       = malloc(sizeof(gr_float));

    *my_fields.density         = rho;
    *my_fields.internal_energy = u_old;
    *my_fields.x_velocity      = SphP[target].VelPred[0];
    *my_fields.x_velocity      = SphP[target].VelPred[0];
    *my_fields.y_velocity      = SphP[target].VelPred[1];
    *my_fields.z_velocity      = SphP[target].VelPred[2];
    // specific heating rate (provide in units [egs s^-1 g^-1]
    my_fields.specific_heating_rate = malloc(sizeof(gr_float));
    *my_fields.specific_heating_rate = udot;

#if (COOL_GRACKLE_CHEMISTRY >=  1)
    my_fields.HI_density       = malloc(sizeof(gr_float));
    my_fields.HII_density      = malloc(sizeof(gr_float));
    my_fields.HeI_density      = malloc(sizeof(gr_float));
    my_fields.HeII_density     = malloc(sizeof(gr_float));
    my_fields.HeIII_density    = malloc(sizeof(gr_float));
    my_fields.e_density        = malloc(sizeof(gr_float));

    *my_fields.HI_density      = rho * SphP[target].grHI;
    *my_fields.HII_density     = rho * SphP[target].grHII;
    *my_fields.HeI_density     = rho * SphP[target].grHeI;
    *my_fields.HeII_density    = rho * SphP[target].grHeII;
    *my_fields.HeIII_density   = rho * SphP[target].grHeIII;
    *my_fields.e_density       = rho * ne_guess;
#endif

#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    my_fields.HM_density       = malloc(sizeof(gr_float));
    my_fields.H2I_density      = malloc(sizeof(gr_float));
    my_fields.H2II_density     = malloc(sizeof(gr_float));

    *my_fields.HM_density      = rho * SphP[target].grHM;
    *my_fields.H2I_density     = rho * SphP[target].grH2I;
    *my_fields.H2II_density    = rho * SphP[target].grH2II;
#endif

#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    my_fields.DI_density       = malloc(sizeof(gr_float));
    my_fields.DII_density      = malloc(sizeof(gr_float));
    my_fields.HDI_density      = malloc(sizeof(gr_float));

    *my_fields.DI_density      = rho * SphP[target].grDI;
    *my_fields.DII_density     = rho * SphP[target].grDII;
    *my_fields.HDI_density     = rho * SphP[target].grHDI;
#endif

#ifdef METALS
    my_fields.metal_density    = malloc(sizeof(gr_float));

    *my_fields.metal_density   = rho * P[target].Metallicity[0];
#endif

    switch(mode) {
        case 0: {  //solve chemistry & update values
            if (solve_chemistry(&All.GrackleUnits, &my_fields, dt) == 0) {
                fprintf(stderr, "Error in solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            // Assign variables back
            gr_float irho = 1 / *my_fields.density;

#if (COOL_GRACKLE_CHEMISTRY >= 1)
            SphP[target].grHI    = irho * *my_fields.HI_density;
            SphP[target].grHII   = irho * *my_fields.HII_density;
            SphP[target].grHeI   = irho * *my_fields.HeI_density;
            SphP[target].grHeII  = irho * *my_fields.HeII_density;
            SphP[target].grHeIII = irho * *my_fields.HeIII_density;
#endif

#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grHM    = irho * *my_fields.HM_density;
            SphP[target].grH2I   = irho * *my_fields.H2I_density;
            SphP[target].grH2II  = irho * *my_fields.H2II_density;
#endif

#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = irho * *my_fields.DI_density;
            SphP[target].grDII   = irho * *my_fields.DII_density;
            SphP[target].grHDI   = irho * *my_fields.HDI_density;
#endif

            returnval = *my_fields.internal_energy;
            break;
        }

        case 1: {  //cooling time
            gr_float *cooling_time;
            cooling_time = malloc(sizeof(gr_float));
            if (calculate_cooling_time(&All.GrackleUnits, &my_fields,
                                       cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *cooling_time;
            break;
        }

        case 2: {  //calculate temperature
            gr_float *temperature;
            temperature = malloc(sizeof(gr_float));
            if (calculate_temperature(&All.GrackleUnits, &my_fields,
                                      temperature) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *temperature;
            break;
        }

        case 3: {  //calculate pressure
            gr_float *pressure;
            pressure = malloc(sizeof(gr_float));
            if (calculate_pressure(&All.GrackleUnits, &my_fields,
                                   pressure) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *pressure;
            break;
        }

        case 4: {  //calculate gamma
            gr_float *gamma;
            gamma = malloc(sizeof(gr_float));
            if (calculate_gamma(&All.GrackleUnits, &my_fields,
                                gamma) == 0) {
                fprintf(stderr, "Error in calculate_gamma.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *gamma;
            break;
        }
    } //end switch

    return returnval;
}

//Initialize Grackle
void InitGrackle(void)
{
    grackle_verbose = 0;
    // Enable output
    if(ThisTask == 0)
        grackle_verbose = 1;

    // First, set up the units system.
    // These are conversions from code units to cgs.
    All.GrackleUnits.comoving_coordinates = 0;
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm / All.HubbleParam;
    All.GrackleUnits.time_units           = All.UnitTime_in_s / All.HubbleParam;
    All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
    All.GrackleUnits.a_units              = 1.0;
    All.GrackleUnits.a_value              = 1.0;
    if (All.ComovingIntegrationOn) {
        // Cosmological case
        All.GrackleUnits.comoving_coordinates = 1;
        All.GrackleUnits.a_value              = All.TimeBegin;
    }

    // Second, create a chemistry object for parameters.
    chemistry_data *my_grackle_data;
    my_grackle_data = malloc(sizeof(chemistry_data));
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
        fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
        exit(ENDRUNVAL);
    }

    // Third, set parameter values for chemistry.

    // Pass parameters
    grackle_data->use_grackle = 1;
#ifndef COOL_GRACKLE_CHEMISTRY
    grackle_data->primordial_chemistry = 0;
#else
    grackle_data->primordial_chemistry = COOL_GRACKLE_CHEMISTRY;
#endif
#ifdef METALS
    grackle_data->metal_cooling = 1;
#else
    grackle_data->metal_cooling = 0;
#endif
    grackle_data->grackle_data_file = All.GrackleDataFile;
    grackle_data->Gamma = GAMMA;
    grackle_data->use_specific_heating_rate = 1;
    grackle_data->SolarMetalFractionByMass = All.SolarAbundance;

    // AGORA defaults
    grackle_data->UVbackground = 1;
    grackle_data-> UVbackground_redshift_on = 15;
    grackle_data->UVbackground_redshift_off = 0;
    grackle_data->UVbackground_redshift_fullon = 15;
    grackle_data->UVbackground_redshift_drop = 0;

    // Grackle defaults
    grackle_data->with_radiative_cooling = 1;
    grackle_data->h2_on_dust = 0;
    grackle_data->cmb_temperature_floor = 1;
    grackle_data->three_body_rate = 0;
    grackle_data->cie_cooling = 0;
    grackle_data->h2_optical_depth_approximation = 0;
    grackle_data->photoelectric_heating = 0;
    grackle_data->photoelectric_heating_rate = 8.5e-26;
    grackle_data->Compton_xray_heating = 0;
    grackle_data->LWbackground_intensity = 0;
    grackle_data->LWbackground_sawtooth_suppression = 0;
    grackle_data->use_volumetric_heating_rate = 0;
    grackle_data->use_radiative_transfer = 0;
    grackle_data->radiative_transfer_coupled_rate_solver = 0;
    grackle_data->radiative_transfer_intermediate_step = 0;
    grackle_data->radiative_transfer_hydrogen_only = 0;
    grackle_data->H2_self_shielding = 0;
    grackle_data->self_shielding_method = 0;

    if(ThisTask == 0)
        printf("Grackle Initialized\n");
}

#endif  //COOL_GRACKLE
