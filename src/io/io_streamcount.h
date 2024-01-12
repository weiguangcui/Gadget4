/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_streamcount.h
 *
 *  \brief class to measure I/O performance
 */

#ifndef IO_STREAMCOUNT_H
#define IO_STREAMCOUNT_H

#include "gadgetconfig.h"

#include <errno.h>
#include <string.h>

#include "../data/macros.h"

class io_streamcount
{
 public:
  long long byte_count = 0;

  /*! \brief  A wrapper for the fwrite() function
   *
   *  This catches I/O errors occuring for fwrite(). In this case we
   *  better stop. If stream is null, no attempt at writing is done.
   *
   *  \param ptr pointer to the beginning of data to write
   *  \param size size in bytes of a single data element
   *  \param nmemb number of elements to be written
   *  \param stream pointer to the output stream
   *  \return number of elements written to stream
   */
  size_t my_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t nwritten;

    if(!stream)
      return 0;

    if(size * nmemb > 0)
      {
        if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
          {
            Terminate("I/O error (fwrite) has occured: %s\n", strerror(errno));
          }
      }
    else
      nwritten = 0;

    byte_count += size * nmemb;

    return nwritten;
  }

  /*! \brief  A wrapper for the fread() function
   *
   *  This catches I/O errors occuring for fread(). In this case we
   *  better stop. If stream is null, no attempt at readingis done.
   *
   *  \param ptr pointer to the beginning of memory location where to store data
   *  \param size size in bytes of a single data element
   *  \param nmemb number of elements to be read
   *  \param stream pointer to the nput stream
   *  \return number of elements read from stream
   */
  size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t nread;

    if(!stream)
      return 0;

    if(size * nmemb > 0)
      {
        if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
          {
            if(feof(stream))
              {
                Terminate("I/O error (fread) has occured: end of file\n");
              }
            else
              Terminate("I/O error (fread) has occured: %s\n", strerror(errno));
          }
      }
    else
      nread = 0;

    byte_count += size * nmemb;

    return nread;
  }

  void reset_io_byte_count(void) { byte_count = 0; }

  long long get_io_byte_count(void) { return byte_count; }
};

#endif
