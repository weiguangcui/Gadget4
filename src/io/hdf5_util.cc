/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  hdf5_util.cc
 *
 *  \brief some helper functions for HDF5 I/O routines
 */

#include "gadgetconfig.h"

#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../io/hdf5_util.h"

#define HALF_ROUND_STYLE 1
#include "../half/half.hpp"
using half_float::half;

hid_t Halfprec_memtype;
hid_t Int48_memtype;
hid_t Int128_memtype;

static herr_t hdf5_conv_half(hid_t src, hid_t dst, H5T_cdata_t *cdata, size_t nelmts, size_t buf_str, size_t bkg_str, void *buf,
                             void *background, hid_t plist)
{
  size_t src_size = H5Tget_size(src);
  size_t dst_size = H5Tget_size(dst);

  char *src_buf = (char *)buf;
  char *dst_buf = (char *)buf;

  int direction;

  switch(cdata->command)
    {
      case H5T_CONV_INIT:
        /*
         * We are being queried to see if we handle this
         * conversion.
         */
        if(H5Tequal(src, Halfprec_memtype) || H5Tequal(dst, Halfprec_memtype))
          {
            cdata->need_bkg = H5T_BKG_NO;
            return 0;
          }
        else
          return -1;
        break;

      case H5T_CONV_FREE:
        break;

      case H5T_CONV_CONV:
        /*
         * Convert each element, watch out for overlap src
         * with dst on the left-most element of the buffer.
         * If the destination size is larger than the source size,
         * then we must process the elements from right to left.
         */

        if(dst_size > src_size)
          {
            direction = -1;
            src_buf += (nelmts - 1) * src_size;
            dst_buf += (nelmts - 1) * dst_size;
          }
        else
          {
            direction = 1;
          }

        for(size_t i = 0; i < nelmts; i++)
          {
            if(src_size == 2)
              {
                if(dst_size == 4)
                  {
                    *((float *)dst_buf) = (float)(*((half *)src_buf));
                  }
                else if(dst_size == 8)
                  {
                    *((double *)dst_buf) = (float)(*((half *)src_buf));
                  }
              }
            else if(dst_size == 2)
              {
                if(src_size == 4)
                  {
                    *((half *)dst_buf) = (half)(*((float *)src_buf));
                  }
                else if(src_size == 8)
                  {
                    *((half *)dst_buf) = (half)(*((double *)src_buf));
                  }
              }
            src_buf += src_size * direction;
            dst_buf += dst_size * direction;
          }

        break;

      default:
        /*
         * Unknown command.
         */
        return -1;
    }
  return 0;
}

void my_create_HDF5_halfprec_handler(void)
{
  /* define the half-precision type */

  Halfprec_memtype = H5Tcopy(H5T_NATIVE_FLOAT);

  H5Tset_fields(Halfprec_memtype, 15, 10, 5, 0, 10);
  H5Tset_ebias(Halfprec_memtype, 15);
  H5Tset_precision(Halfprec_memtype, 16);
  H5Tset_size(Halfprec_memtype, 2);

  H5Tregister(H5T_PERS_SOFT, "half-converter", Halfprec_memtype, Halfprec_memtype, hdf5_conv_half);
}

void my_create_HDF5_special_integer_types(void)
{
  Int48_memtype = H5Tcopy(H5T_NATIVE_UINT32);
  H5Tset_precision(Int48_memtype, 48);

  Int128_memtype = H5Tcopy(H5T_NATIVE_UINT64);
  H5Tset_precision(Int128_memtype, 128);
}

/*! \file hdf5_util.c
 *
 *  \brief Contains the wrapper functions to the HDF5 library functions.
 *
 *  The wrapper functions explicitly check for error conditions and terminate
 *  the run if such conditions occur.
 *
 *  The HDF5 error handler is disabled in case of termination not to repeat
 *  the error message of the handler again at the program exit.
 */

/*!
 * This routine wraps creating a file to give a nice error message
 */
hid_t my_H5Fcreate(const char *fname, unsigned int flags, hid_t fcpl_id, hid_t fapl_id)
{
  hid_t file_id = H5Fcreate(fname, flags, fcpl_id, fapl_id);

  if(file_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to create file %s\n", fname);
    }

  return file_id;
}

/*!
 * This routine wraps creating a group to give a nice error message
 */
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint)
{
  hid_t group_id = H5Gcreate(loc_id, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if(group_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to create group %s\n", groupname);
    }

  return group_id;
}

/*!
 * This routine wraps creating a dataset to give a nice error message
 */
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id)
{
  hid_t dataset_id = H5Dcreate(loc_id, datasetname, type_id, space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

  if(dataset_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to create dataset %s\n", datasetname);
    }

  return dataset_id;
}

/*!
 * This routine wraps writing a dataset to give a nice error message
 */
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf,
                   const char *datasetname)
{
  herr_t status = H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);

  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to write dataset %s\n", datasetname);
    }

  return status;
}

/*!
 * This routine wraps creating an attribute to give a nice error message
 */
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id)
{
  if(H5Aexists(loc_id, attr_name))
    H5Adelete(loc_id, attr_name);

  hid_t attribute_id = H5Acreate(loc_id, attr_name, type_id, space_id, acpl_id, H5P_DEFAULT);

  if(attribute_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to create attribute %s\n", attr_name);
    }

  return attribute_id;
}

/*!
 * This routine wraps writing an attribute to give a nice error message
 */
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name)
{
  herr_t status = H5Awrite(attr_id, mem_type_id, buf);

  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to write attribute %s\n", attr_name);
    }

  return status;
}

/*!
 * This routine wraps creating a dataspace to give a nice error message
 */
hid_t my_H5Screate(H5S_class_t type)
{
  hid_t dataspace_id = H5Screate(type);
  if(dataspace_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      switch(type)
        {
          case H5S_SCALAR:
            Terminate("Error detected in HDF5: unable to create a scalar dataspace\n");
            break;
          case H5S_SIMPLE:
            Terminate("Error detected in HDF5: unable to create a simple dataspace\n");
            break;
          default:
            Terminate("Error detected in HDF5: unknown dataspace type\n");
            break;
        }
    }
  return dataspace_id;
}

/*!
 * This routine wraps creating a simple dataspace to give a nice error message
 */
hid_t my_H5Screate_simple(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims)
{
  hid_t dataspace_id = H5Screate_simple(rank, current_dims, maximum_dims);
  if(dataspace_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to create a simple dataspace\n");
    }
  return dataspace_id;
}

/*!
 * This routine wraps opening a file to give a nice error message
 */
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id)
{
  hid_t file_id = H5Fopen(fname, flags, fapl_id);
  if(file_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to open file %s\n", fname);
    }
  return file_id;
}

/*!
 * This routine wraps opening a group to give a nice error message
 */
hid_t my_H5Gopen(hid_t loc_id, const char *groupname)
{
  hid_t group = H5Gopen(loc_id, groupname, H5P_DEFAULT);
  if(group < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to open group %s\n", groupname);
    }
  return group;
}

/*!
 * This routine wraps opening a dataset to give a nice error message
 */
hid_t my_H5Dopen(hid_t file_id, const char *datasetname)
{
  hid_t dataset = H5Dopen(file_id, datasetname, H5P_DEFAULT);
  if(dataset < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to open dataset %s\n", datasetname);
    }
  return dataset;
}

herr_t my_H5Dset_extent(hid_t dset_id, const hsize_t size[])
{
  herr_t status = H5Dset_extent(dset_id, size);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to set extent of dataset\n");
    }
  return status;
}

/*!
 * This routine wraps opening a dataset. However, in contrast to my_H5Dopen(), if the dataset
 * does not exist it does not terminate the run. This is useful while reading an ICs file
 * because in that case a non-exisitng dataset is put to zero (see also read_ic.c)
 */
hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname)
{
  /* save error handler and disable it */
  H5E_auto_t errfunc;
  void *client_data;
  H5Eget_auto(H5E_DEFAULT, &errfunc, &client_data);
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  hid_t dataset = H5Dopen(file_id, datasetname, H5P_DEFAULT);

  /* reset error handler */
  H5Eset_auto(H5E_DEFAULT, errfunc, client_data);

  return dataset;
}

/*!
 * This routine wraps opening an attribute to give a nice error message
 */
hid_t my_H5Aopen_name(hid_t loc_id, const char *attr_name)
{
  hid_t attribute_id = H5Aopen_name(loc_id, attr_name);
  if(attribute_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to open attribute %s\n", attr_name);
    }
  return attribute_id;
}

/*!
 * This routine wraps reading a dataset to give a nice error message
 */
herr_t my_H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buf,
                  const char *datasetname)
{
  herr_t status = H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to read dataset %s\n", datasetname);
    }
  return status;
}

/*!
 * This routine wraps reading an attribute to give a nice error message
 */
herr_t my_H5Aread(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size)
{
  hid_t hdf5_space   = H5Aget_space(attr_id);
  hssize_t attr_size = H5Sget_simple_extent_npoints(hdf5_space);
  H5Sclose(hdf5_space);

  if(attr_size != size)
    {
      H5E_auto_t errfunc;
      void *client_data;
      H5Eget_auto(H5E_DEFAULT, &errfunc, &client_data);
      errfunc(H5P_DEFAULT, client_data);
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: mismatch in size for attribute %s, expected size = %lld, actual attribute size = %lld\n",
                attr_name, size, attr_size);
    }

  herr_t status = H5Aread(attr_id, mem_type_id, buf);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to read attribute %s\n", attr_name);
    }
  return status;
}

/*!
 * This routine wraps closing an attribute to give a nice error message
 */
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name)
{
  herr_t status = H5Aclose(attr_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to close attribute %s\n", attr_name);
    }
  return status;
}

/*!
 * This routine wraps closing a dataset to give a nice error message
 */
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname)
{
  herr_t status = H5Dclose(dataset_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to close dataset %s\n", datasetname);
    }
  return status;
}

/*!
 * This routine wraps closing a group to give a nice error message
 */
herr_t my_H5Gclose(hid_t group_id, const char *groupname)
{
  herr_t status = H5Gclose(group_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to close group %s\n", groupname);
    }
  return status;
}

herr_t my_H5Pclose(hid_t plist)
{
  herr_t status = H5Pclose(plist);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to close property\n");
    }
  return status;
}

/*!
 * This routine wraps closing a file to give a nice error message
 */
herr_t my_H5Fclose(hid_t file_id, const char *fname)
{
  herr_t status = H5Fclose(file_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to close file %s\n", fname);
    }
  return status;
}

/*!
 * This routine wraps closing a dataspace to give a nice error message
 */
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type)
{
  herr_t status = H5Sclose(dataspace_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      switch(type)
        {
          case H5S_SCALAR:
            Terminate("Error detected in HDF5: unable to close a scalar dataspace\n");
            break;
          case H5S_SIMPLE:
            Terminate("Error detected in HDF5: unable to close a simple dataspace\n");
            break;
          default:
            Terminate("Error detected in HDF5: unknown dataspace type\n");
            break;
        }
    }
  return status;
}

/*!
 * This routine wraps copying an existing datatype to give a nice error message
 */
hid_t my_H5Tcopy(hid_t type_id)
{
  hid_t datatype_id = H5Tcopy(type_id);
  if(datatype_id < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: could not properly copy datatype\n");
    }
  return datatype_id;
}

/*!
 * This routine wraps closing a datatype to give a nice error message
 */
herr_t my_H5Tclose(hid_t type_id)
{
  herr_t status = H5Tclose(type_id);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: could not properly close datatype\n");
    }
  return status;
}

/*!
 * This routine wraps selecting a hyperslab to give a nice error message
 */
herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count,
                              const hsize_t *block)
{
  herr_t status = H5Sselect_hyperslab(space_id, op, start, stride, count, block);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: could not properly select the chosen hyperslab\n");
    }
  return status;
}

/*!
 * This routine wraps returning the size in bytes of a given datatype to give a nice error message
 */
size_t my_H5Tget_size(hid_t datatype_id)
{
  size_t size = H5Tget_size(datatype_id);
  if(size == 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: unable to determine the size of the given datatype\n");
    }
  return size;
}

/*!
 * This routine wraps setting the size in bytes of a given datatype to give a nice error message
 */
herr_t my_H5Tset_size(hid_t datatype_id, size_t size)
{
  herr_t status = H5Tset_size(datatype_id, size);
  if(status < 0)
    {
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);
      Terminate("Error detected in HDF5: could not properly set the size of the given datatype\n");
    }
  return status;
}

void read_scalar_attribute(hid_t handle, const char *attr_name, void *buf, hid_t mem_type_id)
{
  hid_t hdf5_attribute = my_H5Aopen_name(handle, attr_name);
  my_H5Aread(hdf5_attribute, mem_type_id, buf, attr_name, 1);
  my_H5Aclose(hdf5_attribute, attr_name);
}

int read_scalar_attribute(hid_t handle, const char *attr_name, const char *alternative_name, void *buf, hid_t mem_type_id)
{
  if(H5Aexists(handle, attr_name))
    {
      hid_t hdf5_attribute = my_H5Aopen_name(handle, attr_name);
      my_H5Aread(hdf5_attribute, mem_type_id, buf, attr_name, 1);
      my_H5Aclose(hdf5_attribute, attr_name);
      return 0;
    }
  else
    {
      hid_t hdf5_attribute = my_H5Aopen_name(handle, alternative_name);
      my_H5Aread(hdf5_attribute, mem_type_id, buf, alternative_name, 1);
      my_H5Aclose(hdf5_attribute, alternative_name);
      return 1;
    }
}

void read_vector_attribute(hid_t handle, const char *attr_name, void *buf, hid_t mem_type_id, int length)
{
  hid_t hdf5_attribute = my_H5Aopen_name(handle, attr_name);
  my_H5Aread(hdf5_attribute, mem_type_id, buf, attr_name, length);
  my_H5Aclose(hdf5_attribute, attr_name);
}

void write_scalar_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id)
{
  if(H5Aexists(handle, attr_name))
    H5Adelete(handle, attr_name);

  hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);

  hid_t hdf5_attribute = my_H5Acreate(handle, attr_name, mem_type_id, hdf5_dataspace, H5P_DEFAULT);

  my_H5Awrite(hdf5_attribute, mem_type_id, buf, attr_name);

  my_H5Aclose(hdf5_attribute, attr_name);
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}

void write_vector_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id, int length)
{
  if(H5Aexists(handle, attr_name))
    H5Adelete(handle, attr_name);

  hsize_t adim[1] = {(hsize_t)length};

  hid_t hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);

  hid_t hdf5_attribute = my_H5Acreate(handle, attr_name, mem_type_id, hdf5_dataspace, H5P_DEFAULT);

  my_H5Awrite(hdf5_attribute, mem_type_id, buf, attr_name);

  my_H5Aclose(hdf5_attribute, attr_name);
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);
}

void write_string_attribute(hid_t handle, const char *attr_name, const char *buf)
{
  if(H5Aexists(handle, attr_name))
    H5Adelete(handle, attr_name);

  hid_t atype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(atype, strlen(buf));

  hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hid_t hdf5_attribute = my_H5Acreate(handle, attr_name, atype, hdf5_dataspace, H5P_DEFAULT);

  my_H5Awrite(hdf5_attribute, atype, buf, attr_name);

  my_H5Aclose(hdf5_attribute, attr_name);
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}
