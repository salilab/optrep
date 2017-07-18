

/**
 *  \file Cluster.h
 *  \brief Store distance matrix for a bead pair
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPOPTREP_CLUSTER_H
#define IMPOPTREP_CLUSTER_H
#include <IMP.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include <IMP/Object.h>
#include "optrep_config.h"

IMPOPTREP_BEGIN_NAMESPACE

class IMPOPTREPEXPORT Cluster : public Object{
 public:
 
Cluster(Int cc, Ints cm); 

// storage
Int cluster_center;

Ints cluster_members;

IMP_OBJECT_METHODS(Cluster);

};

IMP_OBJECTS(Cluster,Clusters);

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_CLUSTER_H */
