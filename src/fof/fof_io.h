/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof_io.h
 *
 *  \brief declares class needed for I/O of group catalogues
 */

#ifndef FOF_IO_H
#define FOF_IO_H

#include "gadgetconfig.h"

#include "../fof/fof.h"
#include "../io/io.h"

template <typename partset>
class fof_io : public IO_Def
{
 public:
  fof<partset> *FoF;

  fof_io(fof<partset> *FoF_ptr, MPI_Comm comm, int format); /* constructor */

  void fof_subfind_save_groups(int num, const char *basename, const char *grpcat_dirbasename);
  void fof_subfind_load_groups(int num);

  /* supplied virtual functions */
  void fill_file_header(int writeTask, int lastTask, long long *nloc_part, long long *npart);
  void read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *nloc_part, long long *npart,
                        int *nstart);
  void get_datagroup_name(int grnr, char *gname);
  void write_header_fields(hid_t);
  void read_header_fields(const char *fname);
  void read_increase_numbers(int type, int n_for_this_task);
  int get_filenr_from_header(void);
  void set_filenr_in_header(int);
  void *get_base_address_of_structure(enum arrays array, int index);
  int get_type_of_element(int index);
  void set_type_of_element(int index, int type);

 public:
  struct fof_subfind_header
  {
    long long Ngroups;
    long long Nsubhalos;
    long long Nids;
    long long TotNgroups;
    long long TotNsubhalos;
    long long TotNids;
    int num_files;
    double time;
    double redshift;
    double BoxSize;
  };
  fof_subfind_header catalogue_header;

  int LegacyFormat = 0;

 private:
  struct id_list
  {
    MyIDType ID;
    MyHaloNrType GroupNr;
    int Type;
#ifdef SUBFIND
    int SubRankInGr;
    MyFloat BindingEgy;
#endif
  };
  id_list *ID_list;

  void fof_subfind_prepare_ID_list(void);

  static bool fof_subfind_compare_ID_list(const id_list &a, const id_list &b)
  {
    if(a.GroupNr < b.GroupNr)
      return true;
    if(a.GroupNr > b.GroupNr)
      return false;

#ifdef SUBFIND
    if(a.SubRankInGr < b.SubRankInGr)
      return true;
    if(a.SubRankInGr > b.SubRankInGr)
      return false;
#endif

    if(a.Type < b.Type)
      return true;
    if(a.Type > b.Type)
      return false;

#ifdef SUBFIND
    return a.BindingEgy < b.BindingEgy;
#else
    return a.ID < b.ID;
#endif
  }
};

#endif
