/*! \file blackhole_swallow_and_kick.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clarity by parsing the existing code into
 *   smaller files and routines. Some communication and black hole structures were modified
 *   to reduce memory usage. Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#include "blackhole_local.h"

static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_BH_swallowed;

void blackhole_swallow_and_kick_loop(void)
{
    int i, j, k;
    int ndone_flag, ndone;
    int ngrp, recvTask, place, nexport, nimport, dummy;
    MPI_Status status;
    
    int Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed;
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    N_gas_swallowed = N_star_swallowed = N_dm_swallowed = N_BH_swallowed = 0;
    Ntot_gas_swallowed = Ntot_star_swallowed = Ntot_dm_swallowed = Ntot_BH_swallowed = 0;
    
    i = FirstActiveParticle;	/* first particle for this task */
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            if(P[i].Type == 5)
                if(P[i].SwallowID == 0)     /* this particle not being swallowed */
                    if(blackhole_swallow_and_kick_evaluate(i, 0, &nexport, Send_count) < 0)
                        break;
        
        qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        BlackholeDataGet = (struct blackholedata_in *) mymalloc("BlackholeDataGet", nimport * sizeof(struct blackholedata_in));
        BlackholeDataIn = (struct blackholedata_in *) mymalloc("BlackholeDataIn", nexport * sizeof(struct blackholedata_in));
        
        
        /* populate the struct to be exported */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(FLAG_NOT_IN_PUBLIC_CODE_X) || defined(FLAG_NOT_IN_PUBLIC_CODE)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = BlackholeTempInfo[P[place].IndexMapToTempStruc].Jgas_in_Kernel[k];
#endif
            }
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].ID = P[place].ID;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            memcpy(BlackholeDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_G,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_G, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
        
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult", nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_swallow_and_kick_evaluate(j, 1, &dummy, &dummy);  /* set BlackholeDataResult based on BlackholeDataGet */
        
        
        if(i < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_H,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_H, MPI_COMM_WORLD, &status);
                }
            }
        }
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_Mass += BlackholeDataOut[j].Mass;
            BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_BH_Mass += BlackholeDataOut[j].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BPP(place).BH_Mass_AlphaDisk += BlackholeDataOut[j].BH_Mass_AlphaDisk;
#endif
            for(k = 0; k < 3; k++)
                BlackholeTempInfo[P[place].IndexMapToTempStruc].accreted_momentum[k] += BlackholeDataOut[j].accreted_momentum[k];
#ifdef BH_COUNTPROGS
            BPP(place).BH_CountProgs += BlackholeDataOut[j].BH_CountProgs;
#endif
#ifdef GALSF
            if(P[place].StellarAge > BlackholeDataOut[j].Accreted_Age)
                P[place].StellarAge = BlackholeDataOut[j].Accreted_Age;
#endif
        }
        myfree(BlackholeDataOut);
        myfree(BlackholeDataResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    
    MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_star_swallowed, &Ntot_star_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_dm_swallowed, &Ntot_dm_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if((ThisTask == 0)&&(Ntot_gas_swallowed+Ntot_star_swallowed+Ntot_dm_swallowed+Ntot_BH_swallowed>0))
    {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d BH particles\n",
               Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed);
    }
    
}




int blackhole_swallow_and_kick_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, bin, listindex = 0;
    MyIDType id;
    MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3];
    MyFloat *pos, h_i, bh_mass;
    MyFloat f_accreted=0;
    
    MyFloat dir[3], norm, mom;
    mom=0; norm=0; dir[0]=0;
#if defined(FLAG_NOT_IN_PUBLIC_CODE_X) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    MyFloat *Jgas_in_Kernel;
#endif
#ifdef GALSF
    double accreted_age = 1;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat accreted_BH_mass_alphadisk;   
#endif
    
    int mod_index = 0;
    
    if(mode == 0)
    {
        pos = P[target].Pos;
        h_i = PPP[target].Hsml;
        id = P[target].ID;
        bh_mass = BPP(target).BH_Mass;
#if defined(FLAG_NOT_IN_PUBLIC_CODE_X) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        Jgas_in_Kernel = BlackholeTempInfo[P[target].IndexMapToTempStruc].Jgas_in_Kernel;
#endif
        mod_index = P[target].IndexMapToTempStruc;  /* the index of the BlackholeTempInfo should we modify*/
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        h_i = BlackholeDataGet[target].Hsml;
        id = BlackholeDataGet[target].ID;
        bh_mass = BlackholeDataGet[target].BH_Mass;
#if defined(FLAG_NOT_IN_PUBLIC_CODE_X) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
#endif
    }

    
    accreted_mass = 0;
    accreted_BH_mass = 0;
#ifdef BH_ALPHADISK_ACCRETION
    accreted_BH_mass_alphadisk = 0;
#endif
    accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
#ifdef BH_COUNTPROGS
    int accreted_BH_progs = 0;
#endif
    
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_variable_targeted(pos, h_i, target, &startnode, mode, nexport, nSend_local, BH_NEIGHBOR_BITFLAG); // BH_NEIGHBOR_BITFLAG defines which types of particles we search for
            if(numngb < 0) return -1;
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                
                /* we've found a particle to be swallowed.  This could be a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == id && P[j].Mass > 0)
                {
#ifndef IO_REDUCED_MODE
                    printf("found particle P[j].ID = %llu with P[j].SwallowID = %llu of type P[j].Type = %d nearby id = %llu \n",
                           (unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID, P[j].Type, (unsigned long long) id);
#endif
                    /* this is a BH-BH merger */
                    if(P[j].Type == 5)
                    {
//#ifndef IO_REDUCED_MODE   DAA-IO: BH_OUTPUT_MOREINFO overrides IO_REDUCED_MODE
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhMergerDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f\n",
                              All.Time,  id,bh_mass,pos[0],pos[1],pos[2],  P[j].ID,BPP(j).BH_Mass,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]);
#else
#ifndef IO_REDUCED_MODE
                        fprintf(FdBlackHolesDetails,
                                "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
                                ThisTask, All.Time, id, P[j].ID, bh_mass, BPP(j).BH_Mass);
#endif
#endif

#ifdef BH_INCREASE_DYNAMIC_MASS
                        /* the true dynamical mass of the merging BH is P[j].Mass/BH_INCREASE_DYNAMIC_MASS unless exceeded by physical growth
                         - in the limit BPP(j).BH_Mass > BH_INCREASE_DYNAMIC_MASS x m_b, then bh_mass=P[j].Mass on average and we are good as well  */
                        accreted_mass    += FLT( DMAX(BPP(j).BH_Mass, P[j].Mass/BH_INCREASE_DYNAMIC_MASS) );
#else
                        accreted_mass    += FLT(P[j].Mass);
#endif
                        accreted_BH_mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
                        for(k = 0; k < 3; k++)
                            accreted_momentum[k] += FLT(BPP(j).BH_Mass * P[j].Vel[k]);
#ifdef BH_COUNTPROGS
                        accreted_BH_progs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin;
                        TimeBin_BH_mass[bin] -= BPP(j).BH_Mass;
                        TimeBin_BH_dynamicalmass[bin] -= P[j].Mass;
                        TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0)
                            TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;
                        P[j].Mass = 0;
                        BPP(j).BH_Mass = 0;
                        BPP(j).BH_Mdot = 0;
#ifdef GALSF
                        accreted_age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5)
                    
                    

/* DAA: DM and star particles can only be accreted ifdef BH_GRAVCAPTURE_NONGAS */
#ifdef BH_GRAVCAPTURE_NONGAS

                    /* this is a DM particle:
                     In this case, no kick, so just zero out the mass and 'get rid of' the
                     particle (preferably by putting it somewhere irrelevant) */

                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {
#ifndef IO_REDUCED_MODE
                        printf("BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n",
                               j,P[j].Type,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],pos[0],pos[1],pos[2]);
#endif
                        accreted_mass += FLT(P[j].Mass);
                        accreted_BH_mass += FLT(P[j].Mass);
                        P[j].Mass = 0;		// zero out particle mass.  it has now been fully swallowed.
                        N_dm_swallowed++;
                    }


                    /* this is a star particle:
                     If there is an alpha-disk, we let them go to the disk.
                     If there is no alpha-disk, stars go to the BH directly and won't affect feedback.
                     (Can be simply modified if we need something different.) */
                    if((P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {
                        accreted_mass += FLT(P[j].Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        accreted_BH_mass_alphadisk += FLT(P[j].Mass);
#else 
                        accreted_BH_mass += FLT(P[j].Mass);   /* mass goes directly to the BH, not just the parent particle */
#endif
                        P[j].Mass = 0;          // zero out particle mass.  it has now been fully swallowed.
                        N_star_swallowed++;
                    }
                   
#endif // #ifdef BH_GRAVCAPTURE_NONGAS



                    /* this is a gas particle:
                     DAA: we need to see if the gas particle has to be accreted in full or not, depending on BH_WIND_KICK
                     the only difference with BH_ALPHADISK_ACCRETION should be that the mass goes first to the alphadisk */
                    if(P[j].Type == 0)                    
                    {
                        f_accreted = 1;                           // DAA: no "kick winds" so we need to accrete gas particle in full

                        accreted_mass += FLT(f_accreted*P[j].Mass);

#ifdef BH_GRAVCAPTURE_GAS
#ifdef BH_ALPHADISK_ACCRETION       /* mass goes into the alpha disk, before going into the BH */
                        accreted_BH_mass_alphadisk += FLT(f_accreted*P[j].Mass);
#else                               /* mass goes directly to the BH, not just the parent particle */
                        accreted_BH_mass += FLT(f_accreted*P[j].Mass);
#endif
#endif
                        P[j].Mass *= (1-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue *= (1-f_accreted);
#endif



                        /* BAL kicking operations 
                         NOTE: we have two separate BAL wind models, particle kicking and smooth wind model.
                         This is where we do the particle kicking BAL model
                         DAA: This should also work when there is alpha-disk. */

                        N_gas_swallowed++;

                    }  // if(P[j].Type == 0)

                    /* DAA: make sure it is not accreted (or ejected) by the same BH again if inactive in the next timestep */
                    P[j].SwallowID = 0; 
                
                } // if(P[j].SwallowID == id)
                
                
                
                
                

                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        BlackholeTempInfo[mod_index].accreted_Mass = accreted_mass;
        BlackholeTempInfo[mod_index].accreted_BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        // DAA: could be better to include this in BlackholeTempInfo and update BH_Mass_AlphaDisk only at the end (like Mass and BH_Mass)
        BPP(target).BH_Mass_AlphaDisk += accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++) {BlackholeTempInfo[mod_index].accreted_momentum[k] = accreted_momentum[k];}
#ifdef BH_COUNTPROGS
        BPP(target).BH_CountProgs += accreted_BH_progs;
#endif
#ifdef GALSF
        if(P[target].StellarAge > accreted_age)
            P[target].StellarAge = accreted_age;
#endif
    }
    else
    {
        BlackholeDataResult[target].Mass = accreted_mass;
        BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef BH_ALPHADISK_ACCRETION
        BlackholeDataResult[target].BH_Mass_AlphaDisk = accreted_BH_mass_alphadisk;
#endif
        for(k = 0; k < 3; k++) {BlackholeDataResult[target].accreted_momentum[k] = accreted_momentum[k];}
#ifdef BH_COUNTPROGS
        BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
#ifdef GALSF
        BlackholeDataResult[target].Accreted_Age = accreted_age;
#endif
    }
    
    return 0;
} /* closes bh_evaluate_swallow */



