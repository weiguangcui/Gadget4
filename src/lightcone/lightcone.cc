/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  lightcone.cc
 *
 *  \brief contains code to collect particles on the lightcone
 */

#include "gadgetconfig.h"

#ifdef LIGHTCONE

#include <gsl/gsl_rng.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../io/hdf5_util.h"
#include "../lightcone/lightcone.h"
#include "../lightcone/lightcone_massmap_io.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"

void lightcone::lightcone_init_intposconverter(double linklength)
{
  double maxlength = 0;

#ifdef LIGHTCONE_PARTICLES
  if(ConeGlobComDistStart > maxlength)
    maxlength = ConeGlobComDistStart;
#endif

#ifdef LIGHTCONE_MASSMAPS
  for(int i = 0; i < NumMassMapBoundaries; i++)
    if(MassMapBoundariesComDist[i] > maxlength)
      maxlength = MassMapBoundariesComDist[i];
#endif

  maxlength += linklength;

#ifdef LIGHTCONE_PARTICLES
  Lp->RegionLen     = 2.0 * maxlength;
  Lp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Lp->RegionLen;
  Lp->FacIntToCoord = Lp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);
#endif
}

#ifdef LIGHTCONE_MASSMAPS
int lightcone::lightcone_add_position_massmaps(particle_data *P, double *pos, double ascale)
{
  int buffer_full_flag = 0;

  long ipnest;
  vec2pix_ring(All.LightConeMassMapsNside, pos, &ipnest);

  int task;
  subdivide_evenly_get_bin(Mp->Npix, NTask, ipnest, &task);

  if(Mp->NumPart >= Mp->MaxPart)
    {
      int new_maxpart = std::max<int>(Mp->NumPart + 1, (1 + ALLOC_TOLERANCE) * (1 + ALLOC_TOLERANCE) * Mp->MaxPart);
      if(Mp->NumPart >= new_maxpart)
        Terminate("Mp->NumPart >= new_maxpart");

      Mp->reallocate_memory_maxpart(new_maxpart);

      buffer_full_flag = 1;
    }

  int p = Mp->NumPart++;

  Mp->P[p].Ascale = ascale;
  Mp->P[p].setMass(P->getMass());

  Mp->P[p].PixIndex = ipnest;
  Mp->P[p].Task     = task;

  if(task < 0 || task >= NTask || ipnest < 0 || ipnest >= Mp->Npix)
    Terminate("strange assignment:  task=%d  NTask=%d  pos=(%g|%g|%g)  ipnest=%d\n", task, NTask, pos[0], pos[1], pos[2], (int)ipnest);

  return buffer_full_flag;
}
#endif

#ifdef LIGHTCONE_PARTICLES
void lightcone::lightcone_add_position_particles(particle_data *P, double *pos, double ascale, int oindex)
{
  if(Lp->NumPart >= Lp->MaxPart)
    {
      int new_maxpart = (1 + ALLOC_TOLERANCE) * Lp->MaxPart;
      if(Lp->NumPart >= new_maxpart)
        Terminate("Lp->NumPart >= new_maxpart.  Lp->NumPart=%d Lp->MaxPart=%d  new_maxpart=%d", Lp->NumPart, Lp->MaxPart, new_maxpart);

      Lp->reallocate_memory_maxpart(new_maxpart);
    }

  int q = Lp->NumPart++;

  MyIntPosType intpos[3];

  Lp->pos_to_signedintpos(pos, (MySignedIntPosType *)intpos);

  for(int i = 0; i < 3; i++)
    Lp->P[q].IntPos[i] = intpos[i];

  Lp->P[q].Ascale = ascale;

  /* we change to physical velocity here */
  for(int i = 0; i < 3; i++)
    Lp->P[q].Vel[i] = P->Vel[i] / ascale;

#ifdef LIGHTCONE_OUTPUT_ACCELERATIONS
  for(int i = 0; i < 3; i++)
    Lp->P[q].GravAccel[i] = P->GravAccel[i];
#endif

  Lp->P[q].setType(P->getType());
  Lp->P[q].setMass(P->getMass());
  Lp->P[q].ID.set(P->ID.get());
#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
  Lp->P[q].setSofteningClass(P->getSofteningClass());
#endif

  if(P->ID.is_previously_most_bound())
    Lp->P[q].ID.mark_as_formerly_most_bound();

#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
  Lp->P[q].setFlagSaveDistance();
#endif

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  Lp->P[q].OriginIndex = oindex;
#endif
}
#endif

int lightcone::lightcone_test_for_particle_addition(particle_data *P, integertime time0, integertime time1, double dt_drift)
{
  int buffer_full_flag = 0;

  int flag = 0;

#ifdef LIGHTCONE_PARTICLES
  if(time0 < ConeGlobTime_end && time1 > ConeGlobTime_start) /* we have a potential overlap */
    flag = 1;
#endif

#ifdef LIGHTCONE_MASSMAPS
  if(time0 < MassMapBoundariesTime[NumMassMapBoundaries - 1] && time1 > MassMapBoundariesTime[0])
    flag = 1;
#endif

  if(flag == 0)
    return buffer_full_flag;

  static integertime last_time0 = -1, last_time1 = -1;
  static double last_R0 = -1, last_R1 = -1;

  if(time0 != last_time0 || time1 != last_time1)
    {
      last_R0 = Driftfac.get_comoving_distance(time0);
      last_R1 = Driftfac.get_comoving_distance(time1);

      last_time0 = time0;
      last_time1 = time1;
    }

  double R0 = last_R0;
  double R1 = last_R1;

  double R0_squared = R0 * R0;
  double R1_squared = R1 * R1;

  double R1prime = R1 - (R0 - R1);

  double pos[3];
  Sp->intpos_to_pos(P->IntPos, pos);

#ifdef LIGHTCONE_PARTICLES
  bool previously = P->ID.is_previously_most_bound();
#endif

#ifndef LIGHTCONE_MULTIPLE_ORIGINS
  int oindex = 0;
#else
  for(int oindex = 0; oindex < NlightconeOrigins; oindex++)
#endif
  {
    for(int n = 0; n < BoxOrigin[oindex].NumBoxes; n++)
      {
        if(R0 < BoxOrigin[oindex].BoxList[n].Rmin)
          continue;

        if(R1prime > BoxOrigin[oindex].BoxList[n].Rmax)
          break;

        int i = BoxOrigin[oindex].BoxList[n].i;
        int j = BoxOrigin[oindex].BoxList[n].j;
        int k = BoxOrigin[oindex].BoxList[n].k;

        double PosA[3];

        PosA[0] = pos[0] + i * All.BoxSize;
        PosA[1] = pos[1] + j * All.BoxSize;
        PosA[2] = pos[2] + k * All.BoxSize;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
        PosA[0] -= ConeOrigins[oindex].PosOrigin[0];
        PosA[1] -= ConeOrigins[oindex].PosOrigin[1];
        PosA[2] -= ConeOrigins[oindex].PosOrigin[2];
#endif

        double rA2 = PosA[0] * PosA[0] + PosA[1] * PosA[1] + PosA[2] * PosA[2];

        if(rA2 < R0_squared)
          {
            double PosB[3];
            double diffBminusA[3];

            for(int q = 0; q < 3; q++)
              {
                diffBminusA[q] = P->Vel[q] * dt_drift * Sp->FacIntToCoord;
                PosB[q]        = PosA[q] + diffBminusA[q];
              }

            double rB2 = PosB[0] * PosB[0] + PosB[1] * PosB[1] + PosB[2] * PosB[2];

            if(rB2 > R1_squared)
              {
                /* ok, particle crossed the lightcone. Interpolate the coordinate of the crossing */

                double dr2 = diffBminusA[0] * diffBminusA[0] + diffBminusA[1] * diffBminusA[1] + diffBminusA[2] * diffBminusA[2];

                double a = pow(R1 - R0, 2) - dr2;
                double b = 2 * R0 * (R1 - R0) - 2 * (PosA[0] * diffBminusA[0] + PosA[1] * diffBminusA[1] + PosA[2] * diffBminusA[2]);
                double c = R0 * R0 - rA2;

                double det = b * b - 4 * a * c;

                if(det < 0)
                  Terminate(
                      "det=%g R0=%g  R1=%g  rA=%g rB=%g  dr=%g  dt_drift=%g  dx=(%g|%g|%g) vel=(%g|%g|%g)  "
                      "posA=(%g|%g|%g)  posB=(%g|%g|%g)\n",
                      det, R0, R1, sqrt(rA2), sqrt(rB2), sqrt(dr2), dt_drift, P->Vel[0] * dt_drift * Sp->FacIntToCoord,
                      P->Vel[1] * dt_drift * Sp->FacIntToCoord, P->Vel[2] * dt_drift * Sp->FacIntToCoord, P->Vel[0], P->Vel[1],
                      P->Vel[2], PosA[0], PosA[1], PosA[2], PosB[0], PosB[1], PosB[2]);

                double fac = (-b - sqrt(det)) / (2 * a);

                vector<double> Pos;

                for(int q = 0; q < 3; q++)
                  Pos[q] = PosA[q] + fac * diffBminusA[q];

                double ascale = All.TimeBegin * exp((time0 + (time1 - time0) * fac) * All.Timebase_interval);

                /* now we can add particle at position Pos[] to the lightcone, provided it fits into the angular mask */

                if(fac < 0 || fac > 1)
                  {
                    warn(
                        "ascale=%g  fac=%g  fac-1%g R0=%g  R1=%g  rA=%g rB=%g  dr=%g  dt_drift=%g  dx=(%g|%g|%g) vel=(%g|%g|%g)  "
                        "posA=(%g|%g|%g)  posB=(%g|%g|%g)\n",
                        ascale, fac, fac - 1, R0, R1, sqrt(rA2), sqrt(rB2), sqrt(dr2), dt_drift,
                        P->Vel[0] * dt_drift * Sp->FacIntToCoord, P->Vel[1] * dt_drift * Sp->FacIntToCoord,
                        P->Vel[2] * dt_drift * Sp->FacIntToCoord, P->Vel[0], P->Vel[1], P->Vel[2], PosA[0], PosA[1], PosA[2], PosB[0],
                        PosB[1], PosB[2]);
                  }
                else
                  {
#ifdef LIGHTCONE_PARTICLES
                    if(ascale >= ConeGlobAstart && ascale < ConeGlobAend)
                      for(int cone = 0; cone < Nlightcones; cone++)
#ifdef LIGHTCONE_MULTIPLE_ORIGINS
                        if(oindex == Cones[cone].OriginIndex)
#endif
                          if(lightcone_is_cone_member_basic(ascale, Pos, previously, cone))
                            {
                              /* we only add the particle once if it is at least contained in one of the cones */
                              lightcone_add_position_particles(P, Pos.da, ascale, oindex);
                              break;
                            }
#endif

#ifdef LIGHTCONE_MASSMAPS
                    if(ascale >= MassMapBoundariesAscale[0] && ascale < MassMapBoundariesAscale[NumMassMapBoundaries - 1])
                      buffer_full_flag |= lightcone_add_position_massmaps(P, Pos.da, ascale);
#endif
                  }
              }
          }
      }
  }

  return buffer_full_flag;
}

#ifdef LIGHTCONE_PARTICLES

bool lightcone::lightcone_is_cone_member(int i, int cone)
{
#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
  if(!Lp->P[i].getFlagSaveDistance())
    return false;
#endif

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  if(Lp->P[i].OriginIndex != Cones[cone].OriginIndex)
    return false;
#endif

  vector<double> pos;

  if(i >= Lp->NumPart)
    Terminate("i=%d Lp->NumPart=%d\n", i, Lp->NumPart);

  Lp->signedintpos_to_pos((MySignedIntPosType *)Lp->P[i].IntPos, pos.da);

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  int oindex = Cones[cone].OriginIndex;
  pos[0] -= ConeOrigins[oindex].PosOrigin[0];
  pos[1] -= ConeOrigins[oindex].PosOrigin[1];
  pos[2] -= ConeOrigins[oindex].PosOrigin[2];
#endif

  return lightcone_is_cone_member_basic(Lp->P[i].Ascale, pos, Lp->P[i].ID.is_previously_most_bound(), cone);
}

bool lightcone::lightcone_is_cone_member_basic(double ascale, vector<double> &pos, bool previously, int cone)
{
  if(ascale < Cones[cone].Astart || ascale > Cones[cone].Aend)
    return false;

  if(Cones[cone].OnlyMostBoundFlag == 1 && previously == false)
    return false;

  if(Cones[cone].LcType == LC_TYPE_FULLSKY)
    {
      return true;
    }
  else if(Cones[cone].LcType == LC_TYPE_OCTANT)
    {
      double rad   = pos.norm();
      double phi   = atan2(pos[1], pos[0]);
      double theta = acos(pos[2] / rad);

      if(phi < 0)
        phi += 2 * M_PI;

      int octnr = phi / (0.5 * M_PI);

      if(octnr >= 4)
        octnr = 3;

      if(octnr < 0)
        octnr = 0;

      if(theta > 0.5 * M_PI)
        octnr += 4;

      if(octnr == Cones[cone].OctantNr)
        return true;
      else
        return false;
    }
  else if(Cones[cone].LcType == LC_TYPE_PENCIL)
    {
      double rad = pos.norm();

      // note: PencilDir is already normalized

      double angle = acos((pos * Cones[cone].PencilDir) / rad);

      if(angle < Cones[cone].PencilAngleRad)
        return true;
      else
        return false;
    }
  else if(Cones[cone].LcType == LC_TYPE_DISK)
    {
      double dist = pos * Cones[cone].DiskNormal;

      if(fabs(dist) < 0.5 * Cones[cone].DiskThickness)
        return true;
      else
        return false;
    }
  else if(Cones[cone].LcType == LC_TYPE_SQUAREMAP)
    {
      double x = pos * Cones[cone].SquareMapXdir;
      double y = pos * Cones[cone].SquareMapYdir;
      double z = pos * Cones[cone].SquareMapZdir;

      /* not angle has been converted from degrees to rad */
      if(z > 0 && fabs(x) < Cones[cone].SquareMapAngleRad * z && fabs(y) < Cones[cone].SquareMapAngleRad * z)
        return true;
      else
        return false;
    }
  return false;
}

void lightcone::lightcone_init_geometry(char *fname)
{
  if(!All.ComovingIntegrationOn)
    Terminate("LIGHTCONE_PARTICLES: makes only sense for cosmological simulations with ComovingIntegrationOn enabled\n");

  if(ThisTask == 0)
    {
#ifdef LIGHTCONE_MULTIPLE_ORIGINS

      for(int iter = 0; iter < 2; iter++)
        {
          NlightconeOrigins = 0;
          FILE *fd;

          if(!(fd = fopen(All.LightConeOriginsFile, "r")))
            Terminate("LIGHTCONE_MULTIPLE_ORIGINS: cannot read lightcone origins from file `%s'\n", All.LightConeOriginsFile);

          if(iter == 0)
            {
              while(1)
                {
                  double dummy;
                  if(fscanf(fd, "%lg %lg %lg", &dummy, &dummy, &dummy) != 3)
                    break;

                  NlightconeOrigins++;
                  if(NlightconeOrigins > LIGHTCONE_MAX_NUMBER_ORIGINS)
                    Terminate("LIGHTCONE_MULTIPLE_ORIGINS: Too many entries in file %s  (maximum number of origins set to %d)",
                              All.LightConeOriginsFile, LIGHTCONE_MAX_NUMBER_ORIGINS);
                }

              if(NlightconeOrigins == 0)
                Terminate("LIGHTCONE_MULTIPLE_ORIGINS: No entry in file %s", All.LightConeOriginsFile);

              ConeOrigins = (cone_origin *)Mem.mymalloc("ConeOrigins", (NlightconeOrigins + 1) * sizeof(cone_origin));

              mpi_printf("LIGHTCONE_MULTIPLE_ORIGINS: read specification for %d origins from file `%s'.\n", NlightconeOrigins,
                         All.LightConeOriginsFile);
            }
          else
            {
              while(1)
                {
                  if(fscanf(fd, "%lg %lg %lg", &ConeOrigins[NlightconeOrigins].PosOrigin[0],
                            &ConeOrigins[NlightconeOrigins].PosOrigin[1], &ConeOrigins[NlightconeOrigins].PosOrigin[2]) != 3)
                    break;

                  NlightconeOrigins++;
                };
            }

          fclose(fd);
        }

      for(int n = 0; n < NlightconeOrigins; n++)
        {
          for(int i = 0; i < 0; i++)
            {
              while(ConeOrigins[n].PosOrigin[i] < 0)
                ConeOrigins[n].PosOrigin[i] += All.BoxSize;
              while(ConeOrigins[n].PosOrigin[i] >= All.BoxSize)
                ConeOrigins[n].PosOrigin[i] -= All.BoxSize;
            }
          mpi_printf("LIGHTCONE_MULTIPLE_ORIGINS:    Origin #%03d:   %10g %10g %10g\n", n, ConeOrigins[n].PosOrigin[0],
                     ConeOrigins[n].PosOrigin[1], ConeOrigins[n].PosOrigin[2]);
        }
#endif

      for(int iter = 0; iter < 2; iter++)
        {
          Nlightcones = 0;
          FILE *fd;

          if(!(fd = fopen(fname, "r")))
            Terminate("LIGHTCONE_PARTICLES: cannot read lightcone geometries in file `%s'\n", fname);

          if(iter == 0)
            {
              while(1)
                {
                  int lc_type;
                  if(fscanf(fd, "%d", &lc_type) != 1)
                    break;

                  double dummy;

                  switch(lc_type)
                    {
                      case LC_TYPE_FULLSKY:
                        if(fscanf(fd, "%lg %lg %lg", &dummy, &dummy, &dummy) != 3)
                          Terminate("LIGHTCONE_PARTICLES: incorrect data for lightcone type %d in file '%s'", lc_type, fname);
                        break;

                      case LC_TYPE_OCTANT:
                        if(fscanf(fd, "%lg %lg %lg %lg", &dummy, &dummy, &dummy, &dummy) != 4)
                          Terminate("LIGHTCONE_PARTICLES: incorrect data for lightcone type %d in file '%s'", lc_type, fname);
                        break;

                      case LC_TYPE_PENCIL:
                        if(fscanf(fd, "%lg %lg %lg %lg %lg %lg %lg", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != 7)
                          Terminate("LIGHTCONE_PARTICLES: incorrect data for lightcone type %d in file '%s'", lc_type, fname);
                        break;

                      case LC_TYPE_DISK:
                        if(fscanf(fd, "%lg %lg %lg %lg %lg %lg %lg", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != 7)
                          Terminate("LIGHTCONE_PARTICLES: incorrect data for lightcone type %d in file '%s'", lc_type, fname);
                        break;

                      case LC_TYPE_SQUAREMAP:
                        if(fscanf(fd, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,
                                  &dummy, &dummy, &dummy, &dummy) != 10)
                          Terminate("LIGHTCONE_PARTICLES: incorrect data for squaremap type %d in file '%s'", lc_type, fname);
                        break;

                      default:
                        Terminate("LIGHTCONE_PARTICLES: unknown lightcone type %d in file '%s'", lc_type, fname);
                        break;
                    }

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
                  int lc_origin;
                  if(fscanf(fd, "%d", &lc_origin) != 1)
                    Terminate("LIGHTCONE_PARTICLES: can't read origin identifier in file '%s'", fname);
#endif

                  Nlightcones++;
                }

              Cones = (cone_data *)Mem.mymalloc("Cones", Nlightcones * sizeof(cone_data));

              mpi_printf("LIGHTCONE_PARTICLES: read definitions with %d entries in file `%s'.\n", Nlightcones, fname);
            }
          else if(iter == 1)
            {
              while(1)
                {
                  int lc_type;
                  if(fscanf(fd, "%d", &lc_type) != 1)
                    break;

                  Cones[Nlightcones].LcType = lc_type;

                  switch(lc_type)
                    {
                      case LC_TYPE_FULLSKY:
                        fscanf(fd, "%d %lg %lg", &Cones[Nlightcones].OnlyMostBoundFlag, &Cones[Nlightcones].Astart,
                               &Cones[Nlightcones].Aend);
                        snprintf(Cones[Nlightcones].Tag, MAXLEN_PATH, "Full-sky");
                        break;

                      case LC_TYPE_OCTANT:
                        fscanf(fd, "%d %lg %lg", &Cones[Nlightcones].OnlyMostBoundFlag, &Cones[Nlightcones].Astart,
                               &Cones[Nlightcones].Aend);
                        fscanf(fd, "%d", &Cones[Nlightcones].OctantNr);
                        snprintf(Cones[Nlightcones].Tag, MAXLEN_PATH, "Octant");
                        break;

                      case LC_TYPE_PENCIL:
                        fscanf(fd, "%d %lg %lg", &Cones[Nlightcones].OnlyMostBoundFlag, &Cones[Nlightcones].Astart,
                               &Cones[Nlightcones].Aend);
                        fscanf(fd, "%lg %lg %lg", &Cones[Nlightcones].PencilDir[0], &Cones[Nlightcones].PencilDir[1],
                               &Cones[Nlightcones].PencilDir[2]);
                        fscanf(fd, "%lg", &Cones[Nlightcones].PencilAngle);

                        /* normalize the normal vector in case it is not normalized yet */
                        Cones[Nlightcones].PencilDir *= 1.0 / Cones[Nlightcones].PencilDir.norm();

                        /* convert to rad */
                        Cones[Nlightcones].PencilAngleRad = Cones[Nlightcones].PencilAngle * M_PI / 180.0;

                        snprintf(Cones[Nlightcones].Tag, MAXLEN_PATH, "Pencil-Beam");
                        break;

                      case LC_TYPE_DISK:
                        fscanf(fd, "%d %lg %lg", &Cones[Nlightcones].OnlyMostBoundFlag, &Cones[Nlightcones].Astart,
                               &Cones[Nlightcones].Aend);
                        fscanf(fd, "%lg %lg %lg", &Cones[Nlightcones].DiskNormal[0], &Cones[Nlightcones].DiskNormal[1],
                               &Cones[Nlightcones].DiskNormal[2]);
                        fscanf(fd, "%lg", &Cones[Nlightcones].DiskThickness);

                        /* normalize the normal vector in case it is not normalized yet */
                        Cones[Nlightcones].DiskNormal *= 1.0 / Cones[Nlightcones].DiskNormal.norm();
                        snprintf(Cones[Nlightcones].Tag, MAXLEN_PATH, "Disk (for image)");
                        break;

                      case LC_TYPE_SQUAREMAP:
                        fscanf(fd, "%d %lg %lg", &Cones[Nlightcones].OnlyMostBoundFlag, &Cones[Nlightcones].Astart,
                               &Cones[Nlightcones].Aend);
                        fscanf(fd, "%lg %lg %lg", &Cones[Nlightcones].SquareMapZdir[0], &Cones[Nlightcones].SquareMapZdir[1],
                               &Cones[Nlightcones].SquareMapZdir[2]);
                        fscanf(fd, "%lg %lg %lg", &Cones[Nlightcones].SquareMapXdir[0], &Cones[Nlightcones].SquareMapXdir[1],
                               &Cones[Nlightcones].SquareMapXdir[2]);
                        fscanf(fd, "%lg", &Cones[Nlightcones].SquareMapAngle);

                        /* establish coordinate system */

                        Cones[Nlightcones].SquareMapZdir *= 1.0 / Cones[Nlightcones].SquareMapZdir.norm();

                        // cross product
                        Cones[Nlightcones].SquareMapYdir = Cones[Nlightcones].SquareMapZdir ^ Cones[Nlightcones].SquareMapXdir;

                        Cones[Nlightcones].SquareMapYdir *= 1.0 / Cones[Nlightcones].SquareMapYdir.norm();

                        // cross product
                        Cones[Nlightcones].SquareMapXdir = Cones[Nlightcones].SquareMapYdir ^ Cones[Nlightcones].SquareMapZdir;

                        Cones[Nlightcones].SquareMapAngleRad = Cones[Nlightcones].SquareMapAngle * M_PI / 180.0;

                        mpi_printf("LIGHTCONE_SQUAREMAP: cone=%2d   x-axis  =   %15g %15g %15g\n", Nlightcones,
                                   Cones[Nlightcones].SquareMapXdir[0], Cones[Nlightcones].SquareMapXdir[1],
                                   Cones[Nlightcones].SquareMapXdir[2]);
                        mpi_printf("LIGHTCONE_SQUAREMAP: cone=%2d   y-axis  =   %15g %15g %15g\n", Nlightcones,
                                   Cones[Nlightcones].SquareMapYdir[0], Cones[Nlightcones].SquareMapYdir[1],
                                   Cones[Nlightcones].SquareMapYdir[2]);
                        mpi_printf("LIGHTCONE_SQUAREMAP: cone=%2d   z-axis  =   %15g %15g %15g\n", Nlightcones,
                                   Cones[Nlightcones].SquareMapZdir[0], Cones[Nlightcones].SquareMapZdir[1],
                                   Cones[Nlightcones].SquareMapZdir[2]);
                        snprintf(Cones[Nlightcones].Tag, MAXLEN_PATH, "Square-map");
                        break;
                      default:
                        Terminate("odd");
                    }

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
                  fscanf(fd, "%d", &Cones[Nlightcones].OriginIndex);

                  if(Cones[Nlightcones].OriginIndex < 0 || Cones[Nlightcones].OriginIndex >= NlightconeOrigins)
                    Terminate("lightcone origin '%d' out of range, we have only %d origins", Cones[Nlightcones].OriginIndex,
                              NlightconeOrigins);
#endif
                  Nlightcones++;
                }
            }
          fclose(fd);
        }

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
      for(int i = 0; i < Nlightcones; i++)
        mpi_printf("LIGHTCONE_PARTICLES: lightcone #%2d:  %18s %20s  Astart=%10g  Aend=%10g  Origin: (%8g|%8g|%8g)\n", i, Cones[i].Tag,
                   Cones[i].OnlyMostBoundFlag ? "(only most bound)" : "(all particles)", Cones[i].Astart, Cones[i].Aend,
                   ConeOrigins[Cones[i].OriginIndex].PosOrigin[0], ConeOrigins[Cones[i].OriginIndex].PosOrigin[1],
                   ConeOrigins[Cones[i].OriginIndex].PosOrigin[2]);
#else
      for(int i = 0; i < Nlightcones; i++)
        mpi_printf("LIGHTCONE_PARTICLES: lightcone #%2d:  %18s %20s  Astart=%10g  Aend=%10g\n", i, Cones[i].Tag,
                   Cones[i].OnlyMostBoundFlag ? "(only most bound)" : "(all particles)", Cones[i].Astart, Cones[i].Aend);
#endif
      mpi_printf("\n");
    }

  /* now tell also all other ranks about the lightcones */

  MPI_Bcast(&Nlightcones, 1, MPI_INT, 0, Communicator);

  if(ThisTask != 0)
    Cones = (cone_data *)Mem.mymalloc("Cones", Nlightcones * sizeof(cone_data));

  MPI_Bcast(Cones, Nlightcones * sizeof(cone_data), MPI_BYTE, 0, Communicator);

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  MPI_Bcast(&NlightconeOrigins, 1, MPI_INT, 0, Communicator);

  if(ThisTask != 0)
    ConeOrigins = (cone_origin *)Mem.mymalloc("ConeOrigins", (NlightconeOrigins + 1) * sizeof(cone_origin));

  MPI_Bcast(ConeOrigins, NlightconeOrigins * sizeof(cone_origin), MPI_BYTE, 0, Communicator);
#endif
}

int lightcone::lightcone_init_times(void)
{
  for(int i = 0; i < Nlightcones; i++)
    {
      Cones[i].Time_start = log(Cones[i].Astart / All.TimeBegin) / All.Timebase_interval;
      Cones[i].Time_end   = log(Cones[i].Aend / All.TimeBegin) / All.Timebase_interval;

      Cones[i].ComDistStart = Driftfac.get_comoving_distance(Cones[i].Time_start);
      Cones[i].ComDistEnd   = Driftfac.get_comoving_distance(Cones[i].Time_end);
    }

  double astart = 1.0e10;
  double aend   = 0;

  for(int i = 0; i < Nlightcones; i++)
    {
      if(Cones[i].Astart < astart)
        astart = Cones[i].Astart;

      if(Cones[i].Aend > aend)
        aend = Cones[i].Aend;
    }

  ConeGlobAstart = astart;
  ConeGlobAend   = aend;

  ConeGlobTime_start = log(ConeGlobAstart / All.TimeBegin) / All.Timebase_interval;
  ConeGlobTime_end   = log(ConeGlobAend / All.TimeBegin) / All.Timebase_interval;

  ConeGlobComDistStart = Driftfac.get_comoving_distance(ConeGlobTime_start);
  ConeGlobComDistEnd   = Driftfac.get_comoving_distance(ConeGlobTime_end);

  double fac = (4 * M_PI / 3.0) * pow(ConeGlobComDistStart, 3) / pow(All.BoxSize, 3);

  mpi_printf(
      "LIGHTCONE_PARTICLES:  scale_factor: %10g to %10g    comoving distance: %10g to %10g   covered volume in units of box "
      "volume=%g\n",
      ConeGlobAstart, ConeGlobAend, ConeGlobComDistStart, ConeGlobComDistEnd, fac);

  return 0;
}

#endif

int lightcone::lightcone_init_boxlist(void)
{
  double max_GlobComDistStart = 0;

#ifdef LIGHTCONE_PARTICLES
  if(ConeGlobComDistStart > max_GlobComDistStart)
    max_GlobComDistStart = ConeGlobComDistStart;
#endif

#ifdef LIGHTCONE_MASSMAPS
  if(MassMapBoundariesComDist[0] > max_GlobComDistStart)
    max_GlobComDistStart = MassMapBoundariesComDist[0];
#endif

  int n = ceil(max_GlobComDistStart / All.BoxSize + 1);

#ifndef LIGHTCONE_MULTIPLE_ORIGINS
  int oindex = 0;
#else
  for(int oindex = 0; oindex < NlightconeOrigins; oindex++)
#endif
  {
    for(int rep = 0; rep < 2; rep++)
      {
        if(rep == 1)
          BoxOrigin[oindex].BoxList = (boxlist *)Mem.mymalloc_movable(&BoxOrigin[oindex].BoxList, "BoxOrigin[oindex].BoxList",
                                                                      BoxOrigin[oindex].NumBoxes * sizeof(boxlist));

        BoxOrigin[oindex].NumBoxes = 0;

        for(int i = -n; i <= n; i++)
          for(int j = -n; j <= n; j++)
            for(int k = -n; k <= n; k++)
              {
                double corner[3];

                corner[0] = i * All.BoxSize;
                corner[1] = j * All.BoxSize;
                corner[2] = k * All.BoxSize;

                double Rmin, Rmax;

                if(lightcone_box_at_corner_overlaps_at_least_with_one_cone(corner, Rmin, Rmax, oindex))
                  {
                    int num = BoxOrigin[oindex].NumBoxes++;

                    if(rep == 1)
                      {
                        BoxOrigin[oindex].BoxList[num].i    = i;
                        BoxOrigin[oindex].BoxList[num].j    = j;
                        BoxOrigin[oindex].BoxList[num].k    = k;
                        BoxOrigin[oindex].BoxList[num].Rmin = Rmin;
                        BoxOrigin[oindex].BoxList[num].Rmax = Rmax;
                      }
                  }
              }
      }
    mycxxsort(BoxOrigin[oindex].BoxList, BoxOrigin[oindex].BoxList + BoxOrigin[oindex].NumBoxes, lightcone_compare_BoxList_Rmax);
  }

  lightcone_clear_boxlist(All.Time);

  mpi_printf("LIGHTCONE: Number of box replicas to check for first origin lightcone geometry settings = %d\n", BoxOrigin[0].NumBoxes);

  if(BoxOrigin[0].NumBoxes > LIGHTCONE_MAX_BOXREPLICAS)
    {
      mpi_printf(
          "\nLIGHTCONE: Your lightcone extends to such high redshift that the box needs to be replicated a huge number of "
          "times to cover it,\n"
          "more than the prescribed limit of LIGHTCONE_MAX_BOXREPLICAS=%d. We better don't do such an inefficient run, unless you "
          "override this constant.\n",
          LIGHTCONE_MAX_BOXREPLICAS);
      return 1;
    }

  return 0;
}

void lightcone::lightcone_clear_boxlist(double ascale)
{
  double time_start = log(ascale / All.TimeBegin) / All.Timebase_interval;

  double dist = Driftfac.get_comoving_distance(time_start);

#ifndef LIGHTCONE_MULTIPLE_ORIGINS
  int oindex = 0;
#else
  for(int oindex = 0; oindex < NlightconeOrigins; oindex++)
#endif
  {
    int count = 0;

    for(int i = 0; i < BoxOrigin[oindex].NumBoxes; i++)
      {
        if(dist < BoxOrigin[oindex].BoxList[i].Rmin)
          {
            BoxOrigin[oindex].BoxList[i] = BoxOrigin[oindex].BoxList[--BoxOrigin[oindex].NumBoxes];
            i--;
            count++;
          }
      }

    if(count)
      {
        mpi_printf("LIGHTCONE: Eliminated %d entries from BoxList\n", count);
        mycxxsort(BoxOrigin[oindex].BoxList, BoxOrigin[oindex].BoxList + BoxOrigin[oindex].NumBoxes, lightcone_compare_BoxList_Rmax);
      }
  }
}

bool lightcone::lightcone_box_at_corner_overlaps_at_least_with_one_cone(double *corner, double &Rmin, double &Rmax, int oindex)
{
  Rmin = 0;
  Rmax = 0;

  for(int i = 0; i < 3; i++)
    {
      double left = corner[i];
#ifdef LIGHTCONE_MULTIPLE_ORIGINS
      left -= ConeOrigins[oindex].PosOrigin[i];
#endif
      double right = left + All.BoxSize;

      double dx_min = std::min<double>(fabs(left), fabs(right));
      double dx_max = std::max<double>(fabs(left), fabs(right));

      if(left * right <= 0)
        dx_min = 0;

      Rmin += dx_min * dx_min;
      Rmax += dx_max * dx_max;
    }

  Rmin = sqrt(Rmin);
  Rmax = sqrt(Rmax);

#ifdef LIGHTCONE_PARTICLES
  for(int cone = 0; cone < Nlightcones; cone++)
    {
      if(Rmin < Cones[cone].ComDistStart && Rmax > Cones[cone].ComDistEnd)
        {
          if(Cones[cone].LcType == LC_TYPE_FULLSKY)
            {
              return true;
            }
          else if(Cones[cone].LcType == LC_TYPE_OCTANT)
            {
              return true; /* still need to make this more selective */
            }
          else if(Cones[cone].LcType == LC_TYPE_PENCIL)
            {
              return true; /* still need to make this more selective */
            }
          else if(Cones[cone].LcType == LC_TYPE_DISK)
            {
              double mindist    = MAX_DOUBLE_NUMBER;
              double first_dist = 0;

              for(int ii = 0; ii <= 1; ii++)
                for(int jj = 0; jj <= 1; jj++)
                  for(int kk = 0; kk <= 1; kk++)
                    {
                      double crn[3];
                      crn[0] = corner[0] + ii * All.BoxSize;
                      crn[1] = corner[1] + jj * All.BoxSize;
                      crn[2] = corner[2] + kk * All.BoxSize;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
                      crn[0] -= ConeOrigins[oindex].PosOrigin[0];
                      crn[1] -= ConeOrigins[oindex].PosOrigin[1];
                      crn[2] -= ConeOrigins[oindex].PosOrigin[2];
#endif
                      double dist =
                          crn[0] * Cones[cone].DiskNormal[0] + crn[1] * Cones[cone].DiskNormal[1] + crn[2] * Cones[cone].DiskNormal[2];

                      if(ii == 0 && jj == 0 && kk == 0)  // upon first iteration
                        first_dist = dist;
                      else
                        {
                          if(first_dist * dist < 0)  // points on opposite side imply overlap
                            return true;
                        }

                      if(fabs(dist) < mindist)
                        mindist = fabs(dist);

                      if(mindist < 0.5 * Cones[cone].DiskThickness)
                        return true;
                    }
            }
        }
    }
#endif

#ifdef LIGHTCONE_MASSMAPS
  if(Rmin < MassMapBoundariesComDist[0])
    return true;
#endif

  return false;
}

#ifdef LIGHTCONE_MASSMAPS

void lightcone::lightcone_init_massmaps(void)
{
  if(!All.ComovingIntegrationOn)
    Terminate("LIGHTCONE_MASSMAPS: Makes only sense for cosmological simulations with ComovingIntegrationOn enabled\n");

  Mp->Npix = nside2npix(All.LightConeMassMapsNside);

  subdivide_evenly(Mp->Npix, NTask, ThisTask, &Mp->FirstPix, &Mp->NpixLoc);

  if(ThisTask == 0)
    {
      for(int iter = 0; iter < 2; iter++)
        {
          NumMassMapBoundaries = 0;

          if(iter == 0)
            {
              double zend         = 0;
              double com_dist_end = 0;
              NumMassMapBoundaries++;

              while(zend <= All.LightConeMassMapMaxRedshift)
                {
                  com_dist_end += All.LightConeMassMapThickness;

                  double aend = Driftfac.get_scalefactor_for_comoving_distance(com_dist_end);

                  zend = 1 / aend - 1;

                  NumMassMapBoundaries++;
                }

              MassMapBoundariesAscale = (double *)Mem.mymalloc_movable(&MassMapBoundariesAscale, "MassMapBoundariesAscale",
                                                                       NumMassMapBoundaries * sizeof(double));

              mpi_printf("LIGHTCONE_MASSMAPS: %d entries\n", NumMassMapBoundaries);

              if(NumMassMapBoundaries < 2)
                Terminate("Less than two boundaries detected");
            }
          else if(iter == 1)
            {
              double zend                                   = 0;
              double com_dist_end                           = 0;
              MassMapBoundariesAscale[NumMassMapBoundaries] = 1.0;
              NumMassMapBoundaries++;

              while(zend <= All.LightConeMassMapMaxRedshift)
                {
                  com_dist_end += All.LightConeMassMapThickness;

                  double aend = Driftfac.get_scalefactor_for_comoving_distance(com_dist_end);

                  MassMapBoundariesAscale[NumMassMapBoundaries] = aend;

                  zend = 1 / aend - 1;

                  NumMassMapBoundaries++;
                }
            }
        }

      mycxxsort(MassMapBoundariesAscale, MassMapBoundariesAscale + NumMassMapBoundaries, compare_doubles);
    }

  /* now tell also all other ranks about the lightcones */

  MPI_Bcast(&NumMassMapBoundaries, 1, MPI_INT, 0, Communicator);

  if(ThisTask != 0)
    MassMapBoundariesAscale =
        (double *)Mem.mymalloc_movable(&MassMapBoundariesAscale, "MassMapBoundariesAscale", NumMassMapBoundaries * sizeof(double));

  MPI_Bcast(MassMapBoundariesAscale, NumMassMapBoundaries * sizeof(double), MPI_BYTE, 0, Communicator);

  MassMapBoundariesTime =
      (integertime *)Mem.mymalloc_movable(&MassMapBoundariesTime, "MassMapBoundariesTime", NumMassMapBoundaries * sizeof(integertime));
  MassMapBoundariesComDist =
      (double *)Mem.mymalloc_movable(&MassMapBoundariesComDist, "MassMapBoundariesComDist", NumMassMapBoundaries * sizeof(double));
}

int lightcone::lightcone_massmap_report_boundaries(void)
{
  double fac_max               = 0;
  const double allowed_fac_max = LIGHTCONE_MAX_BOXREPLICAS;

  for(int i = 0; i < NumMassMapBoundaries; i++)
    {
      MassMapBoundariesTime[i]    = log(MassMapBoundariesAscale[i] / All.TimeBegin) / All.Timebase_interval;
      MassMapBoundariesComDist[i] = Driftfac.get_comoving_distance(MassMapBoundariesTime[i]);
    }

  for(int i = 0; i < NumMassMapBoundaries - 1; i++)
    {
      MassMapBoundariesTime[i]    = log(MassMapBoundariesAscale[i] / All.TimeBegin) / All.Timebase_interval;
      MassMapBoundariesComDist[i] = Driftfac.get_comoving_distance(MassMapBoundariesTime[i]);

      double fac   = (4 * M_PI / 3.0) * pow(MassMapBoundariesComDist[i], 3) / pow(All.BoxSize, 3);
      double shell = fac - (4 * M_PI / 3.0) * pow(MassMapBoundariesComDist[i + 1], 3) / pow(All.BoxSize, 3);

      mpi_printf(
          "LIGHTCONE_MASSMAPS:   entry=%3d   scale_factor=%10.6f   redshift=%10.6f   comoving distance=%12.3f   shell-volume in "
          "units of box "
          "volume=%g\n",
          i, MassMapBoundariesAscale[i], 1 / MassMapBoundariesAscale[i] - 1, MassMapBoundariesComDist[i], shell);

      if(fac > fac_max)
        fac_max = fac;
    }

  if(fac_max > allowed_fac_max)
    {
      mpi_printf(
          "\nLIGHTCONE_MASSMAPS: Your lightcone mass maps extend to such high redshift that the box needs to be replicated a huge "
          "number of times to cover it (volume ratio %g). We better don't do such an inefficient run.\n"
          "Setting LIGHTCONE_MAX_BOXREPLICAS=%g to a higher value (at least %g) can override this, however.\n",
          fac_max, (double)LIGHTCONE_MAX_BOXREPLICAS, fac_max);
      return 1;
    }

  return 0;
}

void lightcone::lightcone_massmap_flush(int dump_allowed_flag)
{
  lightcone_massmap_binning();

  if(dump_allowed_flag)
    {
      while(All.CurrentMassMapBoundary < NumMassMapBoundaries - 1 &&
            (All.Time >= MassMapBoundariesAscale[All.CurrentMassMapBoundary + 1] || All.Ti_Current >= TIMEBASE))
        {
          lightcone_massmap_io Lcone(Mp, this, Communicator, All.SnapFormat); /* get an I/O object */
          Lcone.lightcone_massmap_save(All.CurrentMassMapBoundary++);

          lightcone_massmap_binning();
        }
    }
}

void lightcone::lightcone_massmap_binning(void)
{
  double t0 = Logs.second();

  int *Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * NTask);

  /* count how many we have of each task */
  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(int i = 0; i < Mp->NumPart; i++)
    {
      int target = Mp->P[i].Task;

      if(target < 0 || target >= NTask)
        Terminate("i=%d: target=%d target < 0 || target >= NTask", i, target);

      if(target != ThisTask)
        Send_count[target]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  Recv_offset[0] = Send_offset[0] = 0;
  int nexport = 0, nimport = 0;

  for(int i = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  lightcone_massmap_data *send_P =
      (lightcone_massmap_data *)Mem.mymalloc_movable(&send_P, "send_P", nexport * sizeof(lightcone_massmap_data));

  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(int i = 0; i < Mp->NumPart; i++)
    {
      int target = Mp->P[i].Task;

      if(target != ThisTask)
        {
          send_P[Send_offset[target] + Send_count[target]] = Mp->P[i];
          Send_count[target]++;

          Mp->P[i] = Mp->P[Mp->NumPart - 1];
          Mp->NumPart--;
          i--;
        }
    }

  if(Mp->NumPart + nimport > Mp->MaxPart)
    Mp->reallocate_memory_maxpart(Mp->NumPart + nimport);

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&send_P[Send_offset[recvTask]], Send_count[recvTask] * sizeof(lightcone_massmap_data), MPI_BYTE, recvTask,
                         TAG_DENS_A, &Mp->P[Mp->NumPart + Recv_offset[recvTask]],
                         Recv_count[recvTask] * sizeof(lightcone_massmap_data), MPI_BYTE, recvTask, TAG_DENS_A, Communicator,
                         MPI_STATUS_IGNORE);
    }

  Mp->NumPart += nimport;

  Mem.myfree_movable(send_P);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  mpi_printf("LIGHTCONE_MASSMAPS: before binning  %9d   (maxpart=%9d, memory for buffer  %g MB)\n", Mp->NumPart, Mp->MaxPart,
             (double)Mp->MaxPart * sizeof(lightcone_massmap_data) * TO_MBYTE_FAC);

  int expunge = 0;

  /* now bin on the map */
  for(int i = 0; i < Mp->NumPart; i++)
    {
      if(Mp->P[i].Task != ThisTask)
        Terminate("can't be");

      if(Mp->P[i].Ascale >= MassMapBoundariesAscale[All.CurrentMassMapBoundary] &&
         Mp->P[i].Ascale < MassMapBoundariesAscale[All.CurrentMassMapBoundary + 1])
        {
          int pix = Mp->P[i].PixIndex - Mp->FirstPix;
          if(pix < 0 || pix >= Mp->NpixLoc)
            Terminate("wrong pixel");

          MassMap[pix] += Mp->P[i].getMass();

          Mp->P[i] = Mp->P[Mp->NumPart - 1];
          Mp->NumPart--;
          i--;
        }
      else if(Mp->P[i].Ascale < MassMapBoundariesAscale[All.CurrentMassMapBoundary])
        {
          Mp->P[i] = Mp->P[Mp->NumPart - 1];
          Mp->NumPart--;
          i--;

          expunge++;
        }
    }

  if(Mp->MaxPart > LIGHTCONE_MASSMAP_ALLOC_FAC * (Sp->TotNumPart / NTask) &&
     Mp->NumPart < LIGHTCONE_MAX_FILLFACTOR * LIGHTCONE_MASSMAP_ALLOC_FAC * (Sp->TotNumPart / NTask))
    Mp->reallocate_memory_maxpart(LIGHTCONE_MASSMAP_ALLOC_FAC * (Sp->TotNumPart / NTask));

  double t1 = Logs.second();

  mpi_printf("LIGHTCONE_MASSMAPS:  after binning  %9d   (maxpart=%9d, memory for buffer  %g MB) took=%g sec expunge=%d\n", Mp->NumPart,
             Mp->MaxPart, (double)Mp->MaxPart * sizeof(lightcone_massmap_data) * TO_MBYTE_FAC, Logs.timediff(t0, t1), expunge);
}

#endif

#endif
