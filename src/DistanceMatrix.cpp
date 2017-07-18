/**
*  \file DistanceMatrix.cpp
*  \brief Store distance matrix corresponding to a bead 
*
*  Copyright 2007-2017 IMP Inventors. All rights reserved.
*
*/
#include <IMP/algebra.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include <iostream>
#include <string>
#include <IMP/rmf.h>
#include <IMP/optrep/DistanceMatrix.h>

IMPOPTREP_BEGIN_NAMESPACE

DistanceMatrix::DistanceMatrix(std::size_t total_number_of_models)  : Object("distancematrix%1%") {
        n=total_number_of_models;
}


IMPOPTREP_END_NAMESPACE
