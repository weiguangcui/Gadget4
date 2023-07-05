/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  hdf5_util.h
 *
 *  \brief declares some helper functions for HDF5 I/O routines
 */

#ifndef HDF5_UTIL_H
#define HDF5_UTIL_H

#include "gadgetconfig.h"

#include <hdf5.h>

#define COMPRESSION_CHUNKSIZE 1000

extern hid_t Halfprec_memtype;
extern hid_t Int48_memtype;
extern hid_t Int128_memtype;

void my_create_HDF5_special_integer_types(void);
void my_create_HDF5_halfprec_handler(void);
hid_t my_H5Fcreate(const char *fname, unsigned flags, hid_t fcpl_id, hid_t fapl_id);
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint);
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id);
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
hid_t my_H5Screate_simple(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims);
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf,
                   const char *datasetname);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id);
hid_t my_H5Dopen(hid_t file_id, const char *datasetname);
hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname);
herr_t my_H5Dset_extent(hid_t dset_id, const hsize_t size[]);
herr_t my_H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buf,
                  const char *datasetname);
hid_t my_H5Gopen(hid_t loc_id, const char *groupname);
hid_t my_H5Aopen_name(hid_t loc_id, const char *attr_name);
herr_t my_H5Aread(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size);

herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname);
herr_t my_H5Gclose(hid_t group_id, const char *groupname);
herr_t my_H5Fclose(hid_t file_id, const char *fname);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);
herr_t my_H5Pclose(hid_t plist);

hid_t my_H5Tcopy(hid_t type_id);
herr_t my_H5Tclose(hid_t type_id);

herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count,
                              const hsize_t *block);
size_t my_H5Tget_size(hid_t datatype_id);
herr_t my_H5Tset_size(hid_t datatype_id, size_t size);

void write_scalar_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id);
void write_vector_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id, int length);
void write_string_attribute(hid_t handle, const char *attr_name, const char *buf);

void read_scalar_attribute(hid_t handle, const char *attr_name, void *buf, hid_t mem_type_id);
int read_scalar_attribute(hid_t handle, const char *attr_name, const char *alternative_name, void *buf, hid_t mem_type_id);
void read_vector_attribute(hid_t handle, const char *attr_name, void *buf, hid_t mem_type_id, int length);

#endif
