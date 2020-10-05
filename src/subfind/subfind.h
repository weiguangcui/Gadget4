/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind.h
 *
 *  \brief defines various constants used by the Subfind code
 */

#ifndef SUBFIND_H
#define SUBFIND_H

#ifdef SUBFIND

#define FIND_SMOOTHING_LENGTHS 0
#define FIND_TOTAL_DENSITIES 1

#define RECOMPUTE_ALL 0
#define UPDATE_ALL 1

#define MAX_ITER_UNBIND 500

#define TAG_POLLING_DONE 201
#define TAG_SET_ALL 202
#define TAG_GET_NGB_INDICES 204
#define TAG_GET_TAILANDLEN 205
#define TAG_GET_TAILANDLEN_DATA 206
#define TAG_SET_TAILANDLEN 207
#define TAG_SET_HEADANDNEXT 209
#define TAG_SETHEADGETNEXT_DATA 210
#define TAG_SET_NEXT 211
#define TAG_SETHEADGETNEXT 213
#define TAG_GET_NEXT 215
#define TAG_GET_NEXT_DATA 216
#define TAG_GET_HEAD 217
#define TAG_GET_HEAD_DATA 218
#define TAG_ADD_PARTICLE 219
#define TAG_ADDBOUND 220
#define TAG_NID 222
#define TAG_NID_DATA 223
#define TAG_SETRANK 224
#define TAG_SETRANK_OUT 226
#define TAG_GET_RANK 227
#define TAG_GET_RANK_DATA 228
#define TAG_MARK_PARTICLE 229
#define TAG_SET_NEWTAIL 230
#define TAG_GET_OLDTAIL 231
#define TAG_GET_TWOHEADS 232
#define TAG_GET_TWOHEADS_DATA 233

#endif
#endif
