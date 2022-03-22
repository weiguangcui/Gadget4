/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_orphanids.cc
 *
 *  \brief determine orphan IDs so that those are flagged
 */

#include "gadgetconfig.h"

#ifdef SUBFIND_ORPHAN_TREATMENT

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/io_readsnap.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../system/system.h"

template <>
void fof<simparticles>::subfind_match_ids_of_previously_most_bound_ids(simparticles *Sp)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  mycxxsort_parallel(Sp->IdStore.ID, Sp->IdStore.ID + Sp->IdStore.NumPart, Sp->compare_IDs, Communicator);
  mycxxsort_parallel(Sp->P, Sp->P + Sp->NumPart, Sp->compare_SpP_ID, Communicator);

  MyIDType *list_min_id = (MyIDType *)Mem.mymalloc("list_min_id", NTask * sizeof(MyIDType));
  MyIDType *list_max_id = (MyIDType *)Mem.mymalloc("list_max_id", NTask * sizeof(MyIDType));

  MyIDType idmin = Sp->P[0].ID.get();
  MyIDType idmax = Sp->P[Sp->NumPart - 1].ID.get();

  MPI_Allgather(&idmin, sizeof(MyIDType), MPI_BYTE, list_min_id, sizeof(MyIDType), MPI_BYTE, Communicator);
  MPI_Allgather(&idmax, sizeof(MyIDType), MPI_BYTE, list_max_id, sizeof(MyIDType), MPI_BYTE, Communicator);

  int *num_list = (int *)Mem.mymalloc("num_list", NTask * sizeof(int));
  MPI_Allgather(&Sp->NumPart, 1, MPI_INT, num_list, 1, MPI_INT, Communicator);

  int nexport = 0, nimport = 0;
  MyIDType *import_data = NULL, *export_data = NULL;

  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < Sp->IdStore.NumPart; i++)
        {
          while(target < NTask - 1 && (num_list[target] == 0 || Sp->IdStore.ID[i] > list_max_id[target]))
            target++;

          if(num_list[target] == 0)
            Terminate("How can this be? target=%d", target);

          if(Sp->IdStore.ID[i] >= list_min_id[target] && Sp->IdStore.ID[i] <= list_max_id[target])
            {
              if(mode == 0)
                Send_count[target]++;
              else
                {
                  int off          = Send_offset[target] + Send_count[target]++;
                  export_data[off] = Sp->IdStore.ID[i];
                }
            }
        }

      if(mode == 0)
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
          Recv_offset[0] = Send_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nimport += Recv_count[j];
              nexport += Send_count[j];
              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          export_data = (MyIDType *)Mem.mymalloc("export_data", nexport * sizeof(MyIDType));
          import_data = (MyIDType *)Mem.mymalloc("import_data", nimport * sizeof(MyIDType));
        }
    }

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(MyIDType), MPI_BYTE, recvTask, TAG_DENS_B,
                         &import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(MyIDType), MPI_BYTE, recvTask, TAG_DENS_B,
                         Communicator, MPI_STATUS_IGNORE);
    }

  /* incoming data should already be sorted, so now do the match */

  int nmarked = 0;
  for(int i = 0, j = 0; i < Sp->NumPart && j < nimport;)
    {
      if(Sp->P[i].ID.get() < import_data[j])
        i++;
      else if(Sp->P[i].ID.get() > import_data[j])
        j++;
      else
        {
          if(!Sp->P[i].ID.is_previously_most_bound())
            {
              Sp->P[i].ID.mark_as_formerly_most_bound();
              nmarked++;
            }
          i++;
          j++;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);

  Mem.myfree(num_list);
  Mem.myfree(list_max_id);
  Mem.myfree(list_min_id);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  long long tot_ncheck, tot_nmarked;
  sumup_large_ints(1, &Sp->IdStore.NumPart, &tot_ncheck, Communicator);
  sumup_large_ints(1, &nmarked, &tot_nmarked, Communicator);

  mpi_printf("SUBFIND_ORPHAN_TREATMENT: Got %lld particles from previous snapshot, led to %lld additionally marked particles\n",
             tot_ncheck, tot_nmarked);
}

#endif
