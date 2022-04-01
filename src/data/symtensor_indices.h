/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file symtensor_indices.h
 *
 *  \brief defines some symbols for accessing the elements of (storage-optimized) symmetric tensors
 */

#ifndef SRC_DATA_SYMTENSOR_INDICES_H_
#define SRC_DATA_SYMTENSOR_INDICES_H_

/* 1-vector*/
#define vX 0
#define vY 1
#define vZ 2

/* 2-tensor element mapping, symmetric 3x3 */
#define qXX 0
#define qXY 1
#define qXZ 2
#define qYX qXY
#define qYY 3
#define qYZ 4
#define qZX qXZ
#define qZY qYZ
#define qZZ 5

/* 3-tensor element mapping, symmetric 3x3x3 */
#define dXXX 0
#define dXXY 1
#define dXXZ 2

#define dXYX dXXY
#define dXYY 3
#define dXYZ 4

#define dXZX dXXZ
#define dXZY dXYZ
#define dXZZ 5

#define dYXX dXXY
#define dYXY dXYY
#define dYXZ dXYZ

#define dYYX dXYY
#define dYYY 6
#define dYYZ 7

#define dYZX dXYZ
#define dYZY dYYZ
#define dYZZ 8

#define dZXX dXXZ
#define dZXY dXYZ
#define dZXZ dXZZ

#define dZYX dXYZ
#define dZYY dYYZ
#define dZYZ dYZZ

#define dZZX dXZZ
#define dZZY dYZZ
#define dZZZ 9

/* 4-tensor element mapping, symmetric 3x3x3x3 */
#define sXXXX 0
#define sXXXY 1
#define sXXXZ 2

#define sXXYX sXXXY
#define sXXYY 3
#define sXXYZ 4

#define sXXZX sXXXZ
#define sXXZY sXXYZ
#define sXXZZ 5

#define sXYXX sXXXY
#define sXYXY sXXYY
#define sXYXZ sXXYZ

#define sXYYX sXXYY
#define sXYYY 6
#define sXYYZ 7

#define sXYZX sXXYZ
#define sXYZY sXYYZ
#define sXYZZ 8

#define sXZXX sXXXZ
#define sXZXY sXXYZ
#define sXZXZ sXXZZ

#define sXZYX sXXYZ
#define sXZYY sXYYZ
#define sXZYZ sXYZZ

#define sXZZX sXXZZ
#define sXZZY sXYZZ
#define sXZZZ 9
//-----------------
#define sYXXX sXXXY
#define sYXXY sXXYY
#define sYXXZ sXXYZ

#define sYXYX sXXYY
#define sYXYY sXYYY
#define sYXYZ sXYYZ

#define sYXZX sXXYZ
#define sYXZY sXYYZ
#define sYXZZ sXYZZ

#define sYYXX sXXYY
#define sYYXY sXYYY
#define sYYXZ sXYYZ

#define sYYYX sXYYY
#define sYYYY 10
#define sYYYZ 11

#define sYYZX sXYYZ
#define sYYZY sYYYZ
#define sYYZZ 12

#define sYZXX sXXYZ
#define sYZXY sXYYZ
#define sYZXZ sXYZZ

#define sYZYX sXYYZ
#define sYZYY sYYYZ
#define sYZYZ sYYZZ

#define sYZZX sXYZZ
#define sYZZY sYYZZ
#define sYZZZ 13
//-----------------
#define sZXXX sXXXZ
#define sZXXY sXXYZ
#define sZXXZ sXXZZ

#define sZXYX sXXYZ
#define sZXYY sXYYZ
#define sZXYZ sXYZZ

#define sZXZX sXXZZ
#define sZXZY sXYZZ
#define sZXZZ sXZZZ

#define sZYXX sXXYZ
#define sZYXY sXYYZ
#define sZYXZ sXYZZ

#define sZYYX sXYYZ
#define sZYYY sYYYZ
#define sZYYZ sYYZZ

#define sZYZX sXYZZ
#define sZYZY sYYZZ
#define sZYZZ sYZZZ

#define sZZXX sXXZZ
#define sZZXY sXYZZ
#define sZZXZ sXZZZ

#define sZZYX sXYZZ
#define sZZYY sYYZZ
#define sZZYZ sYZZZ

#define sZZZX sXZZZ
#define sZZZY sYZZZ
#define sZZZZ 14

/* 5-tensor element mapping, symmetric 3x3x3x3x3 */
#define rXXXXX 0
#define rXXXXY 1
#define rXXXXZ 2
#define rXXXYX rXXXXY
#define rXXXYY 3
#define rXXXYZ 4
#define rXXXZX rXXXXZ
#define rXXXZY rXXXYZ
#define rXXXZZ 5
#define rXXYXX rXXXXY
#define rXXYXY rXXXYY
#define rXXYXZ rXXXYZ
#define rXXYYX rXXXYY
#define rXXYYY 6
#define rXXYYZ 7
#define rXXYZX rXXXYZ
#define rXXYZY rXXYYZ
#define rXXYZZ 8
#define rXXZXX rXXXXZ
#define rXXZXY rXXXYZ
#define rXXZXZ rXXXZZ
#define rXXZYX rXXXYZ
#define rXXZYY rXXYYZ
#define rXXZYZ rXXYZZ
#define rXXZZX rXXXZZ
#define rXXZZY rXXYZZ
#define rXXZZZ 9
#define rXYXXX rXXXXY
#define rXYXXY rXXXYY
#define rXYXXZ rXXXYZ
#define rXYXYX rXXXYY
#define rXYXYY rXXYYY
#define rXYXYZ rXXYYZ
#define rXYXZX rXXXYZ
#define rXYXZY rXXYYZ
#define rXYXZZ rXXYZZ
#define rXYYXX rXXXYY
#define rXYYXY rXXYYY
#define rXYYXZ rXXYYZ
#define rXYYYX rXXYYY
#define rXYYYY 10
#define rXYYYZ 11
#define rXYYZX rXXYYZ
#define rXYYZY rXYYYZ
#define rXYYZZ 12
#define rXYZXX rXXXYZ
#define rXYZXY rXXYYZ
#define rXYZXZ rXXYZZ
#define rXYZYX rXXYYZ
#define rXYZYY rXYYYZ
#define rXYZYZ rXYYZZ
#define rXYZZX rXXYZZ
#define rXYZZY rXYYZZ
#define rXYZZZ 13
#define rXZXXX rXXXXZ
#define rXZXXY rXXXYZ
#define rXZXXZ rXXXZZ
#define rXZXYX rXXXYZ
#define rXZXYY rXXYYZ
#define rXZXYZ rXXYZZ
#define rXZXZX rXXXZZ
#define rXZXZY rXXYZZ
#define rXZXZZ rXXZZZ
#define rXZYXX rXXXYZ
#define rXZYXY rXXYYZ
#define rXZYXZ rXXYZZ
#define rXZYYX rXXYYZ
#define rXZYYY rXYYYZ
#define rXZYYZ rXYYZZ
#define rXZYZX rXXYZZ
#define rXZYZY rXYYZZ
#define rXZYZZ rXYZZZ
#define rXZZXX rXXXZZ
#define rXZZXY rXXYZZ
#define rXZZXZ rXXZZZ
#define rXZZYX rXXYZZ
#define rXZZYY rXYYZZ
#define rXZZYZ rXYZZZ
#define rXZZZX rXXZZZ
#define rXZZZY rXYZZZ
#define rXZZZZ 14
#define rYXXXX rXXXXY
#define rYXXXY rXXXYY
#define rYXXXZ rXXXYZ
#define rYXXYX rXXXYY
#define rYXXYY rXXYYY
#define rYXXYZ rXXYYZ
#define rYXXZX rXXXYZ
#define rYXXZY rXXYYZ
#define rYXXZZ rXXYZZ
#define rYXYXX rXXXYY
#define rYXYXY rXXYYY
#define rYXYXZ rXXYYZ
#define rYXYYX rXXYYY
#define rYXYYY rXYYYY
#define rYXYYZ rXYYYZ
#define rYXYZX rXXYYZ
#define rYXYZY rXYYYZ
#define rYXYZZ rXYYZZ
#define rYXZXX rXXXYZ
#define rYXZXY rXXYYZ
#define rYXZXZ rXXYZZ
#define rYXZYX rXXYYZ
#define rYXZYY rXYYYZ
#define rYXZYZ rXYYZZ
#define rYXZZX rXXYZZ
#define rYXZZY rXYYZZ
#define rYXZZZ rXYZZZ
#define rYYXXX rXXXYY
#define rYYXXY rXXYYY
#define rYYXXZ rXXYYZ
#define rYYXYX rXXYYY
#define rYYXYY rXYYYY
#define rYYXYZ rXYYYZ
#define rYYXZX rXXYYZ
#define rYYXZY rXYYYZ
#define rYYXZZ rXYYZZ
#define rYYYXX rXXYYY
#define rYYYXY rXYYYY
#define rYYYXZ rXYYYZ
#define rYYYYX rXYYYY
#define rYYYYY 15
#define rYYYYZ 16
#define rYYYZX rXYYYZ
#define rYYYZY rYYYYZ
#define rYYYZZ 17
#define rYYZXX rXXYYZ
#define rYYZXY rXYYYZ
#define rYYZXZ rXYYZZ
#define rYYZYX rXYYYZ
#define rYYZYY rYYYYZ
#define rYYZYZ rYYYZZ
#define rYYZZX rXYYZZ
#define rYYZZY rYYYZZ
#define rYYZZZ 18
#define rYZXXX rXXXYZ
#define rYZXXY rXXYYZ
#define rYZXXZ rXXYZZ
#define rYZXYX rXXYYZ
#define rYZXYY rXYYYZ
#define rYZXYZ rXYYZZ
#define rYZXZX rXXYZZ
#define rYZXZY rXYYZZ
#define rYZXZZ rXYZZZ
#define rYZYXX rXXYYZ
#define rYZYXY rXYYYZ
#define rYZYXZ rXYYZZ
#define rYZYYX rXYYYZ
#define rYZYYY rYYYYZ
#define rYZYYZ rYYYZZ
#define rYZYZX rXYYZZ
#define rYZYZY rYYYZZ
#define rYZYZZ rYYZZZ
#define rYZZXX rXXYZZ
#define rYZZXY rXYYZZ
#define rYZZXZ rXYZZZ
#define rYZZYX rXYYZZ
#define rYZZYY rYYYZZ
#define rYZZYZ rYYZZZ
#define rYZZZX rXYZZZ
#define rYZZZY rYYZZZ
#define rYZZZZ 19
#define rZXXXX rXXXXZ
#define rZXXXY rXXXYZ
#define rZXXXZ rXXXZZ
#define rZXXYX rXXXYZ
#define rZXXYY rXXYYZ
#define rZXXYZ rXXYZZ
#define rZXXZX rXXXZZ
#define rZXXZY rXXYZZ
#define rZXXZZ rXXZZZ
#define rZXYXX rXXXYZ
#define rZXYXY rXXYYZ
#define rZXYXZ rXXYZZ
#define rZXYYX rXXYYZ
#define rZXYYY rXYYYZ
#define rZXYYZ rXYYZZ
#define rZXYZX rXXYZZ
#define rZXYZY rXYYZZ
#define rZXYZZ rXYZZZ
#define rZXZXX rXXXZZ
#define rZXZXY rXXYZZ
#define rZXZXZ rXXZZZ
#define rZXZYX rXXYZZ
#define rZXZYY rXYYZZ
#define rZXZYZ rXYZZZ
#define rZXZZX rXXZZZ
#define rZXZZY rXYZZZ
#define rZXZZZ rXZZZZ
#define rZYXXX rXXXYZ
#define rZYXXY rXXYYZ
#define rZYXXZ rXXYZZ
#define rZYXYX rXXYYZ
#define rZYXYY rXYYYZ
#define rZYXYZ rXYYZZ
#define rZYXZX rXXYZZ
#define rZYXZY rXYYZZ
#define rZYXZZ rXYZZZ
#define rZYYXX rXXYYZ
#define rZYYXY rXYYYZ
#define rZYYXZ rXYYZZ
#define rZYYYX rXYYYZ
#define rZYYYY rYYYYZ
#define rZYYYZ rYYYZZ
#define rZYYZX rXYYZZ
#define rZYYZY rYYYZZ
#define rZYYZZ rYYZZZ
#define rZYZXX rXXYZZ
#define rZYZXY rXYYZZ
#define rZYZXZ rXYZZZ
#define rZYZYX rXYYZZ
#define rZYZYY rYYYZZ
#define rZYZYZ rYYZZZ
#define rZYZZX rXYZZZ
#define rZYZZY rYYZZZ
#define rZYZZZ rYZZZZ
#define rZZXXX rXXXZZ
#define rZZXXY rXXYZZ
#define rZZXXZ rXXZZZ
#define rZZXYX rXXYZZ
#define rZZXYY rXYYZZ
#define rZZXYZ rXYZZZ
#define rZZXZX rXXZZZ
#define rZZXZY rXYZZZ
#define rZZXZZ rXZZZZ
#define rZZYXX rXXYZZ
#define rZZYXY rXYYZZ
#define rZZYXZ rXYZZZ
#define rZZYYX rXYYZZ
#define rZZYYY rYYYZZ
#define rZZYYZ rYYZZZ
#define rZZYZX rXYZZZ
#define rZZYZY rYYZZZ
#define rZZYZZ rYZZZZ
#define rZZZXX rXXZZZ
#define rZZZXY rXYZZZ
#define rZZZXZ rXZZZZ
#define rZZZYX rXYZZZ
#define rZZZYY rYYZZZ
#define rZZZYZ rYZZZZ
#define rZZZZX rXZZZZ
#define rZZZZY rYZZZZ
#define rZZZZZ 20

/* 6-tensor element mapping, symmetric 3x3x3x3x3x3 */
#define pXXXXXX 0
#define pXXXXXY 1
#define pXXXXXZ 2
#define pXXXXYY 3
#define pXXXXYZ 4
#define pXXXXZZ 5
#define pXXXYYY 6
#define pXXXYYZ 7
#define pXXXYZZ 8
#define pXXXZZZ 9
#define pXXYYYY 10
#define pXXYYYZ 11
#define pXXYYZZ 12
#define pXXYZZZ 13
#define pXXZZZZ 14
#define pXYYYYY 15
#define pXYYYYZ 16
#define pXYYYZZ 17
#define pXYYZZZ 18
#define pXYZZZZ 19
#define pXZZZZZ 20
#define pYYYYYY 21
#define pYYYYYZ 22
#define pYYYYZZ 23
#define pYYYZZZ 24
#define pYYZZZZ 25
#define pYZZZZZ 26
#define pZZZZZZ 27

/* 7-tensor element mapping, symmetric 3x3x3x3x3x3x3 */
#define tXXXXXXX 0
#define tXXXXXXY 1
#define tXXXXXXZ 2
#define tXXXXXYY 3
#define tXXXXXYZ 4
#define tXXXXXZZ 5
#define tXXXXYYY 6
#define tXXXXYYZ 7
#define tXXXXYZZ 8
#define tXXXXZZZ 9
#define tXXXYYYY 10
#define tXXXYYYZ 11
#define tXXXYYZZ 12
#define tXXXYZZZ 13
#define tXXXZZZZ 14
#define tXXYYYYY 15
#define tXXYYYYZ 16
#define tXXYYYZZ 17
#define tXXYYZZZ 18
#define tXXYZZZZ 19
#define tXXZZZZZ 20
#define tXYYYYYY 21
#define tXYYYYYZ 22
#define tXYYYYZZ 23
#define tXYYYZZZ 24
#define tXYYZZZZ 25
#define tXYZZZZZ 26
#define tXZZZZZZ 27
#define tYYYYYYY 28
#define tYYYYYYZ 29
#define tYYYYYZZ 30
#define tYYYYZZZ 31
#define tYYYZZZZ 32
#define tYYZZZZZ 33
#define tYZZZZZZ 34
#define tZZZZZZZ 35

#endif /* SRC_DATA_SYMTENSOR_INDICES_H_ */
