/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  peano.h
 *
 *  \brief declaration of function prototypes used for Peano-Hilbert keys
 */

#ifndef SORT_H
#define SORT_H

peanokey peano_hilbert_key(MyIntPosType x, MyIntPosType y, MyIntPosType z, int bits);
void peano_hilbert_key_inverse(peanokey key, int bits, MyIntPosType *x, MyIntPosType *y, MyIntPosType *z);

unsigned char peano_incremental_key(unsigned char pix, unsigned char *rotation);

#endif
