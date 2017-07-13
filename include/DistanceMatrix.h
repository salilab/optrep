
/**
 *  \file DistanceMatrix.h
 *  \brief Store distance matrix for a bead pair
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPOPTREP_DISTANCE_MATRIX_H
#define IMPOPTREP_DISTANCE_MATRIX_H
#include <IMP.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include "optrep_config.h"

IMPOPTREP_BEGIN_NAMESPACE

class IMPOPTREPEXPORT DistanceMatrix {
 public:
 /**
    \param[in] total_number_of_models sets the size of the distance matrix 
  */

void DistanceMatrix(std::size_t total_number_of_models); 

// total number of models
 std::size_t n;
 
 Float mindist;
 Float maxdist; 

 Floats distmat; 
 // flattened distance matrix
 
};

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_DISTANCE_MATRIX_H */
