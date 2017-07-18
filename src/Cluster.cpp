/**
*  \file Cluster.cpp
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
#include <IMP/optrep/Cluster.h>

IMPOPTREP_BEGIN_NAMESPACE

Cluster::Cluster(Int cc, Ints cm) : Object("cluster%1%"){
    
    cluster_center = cc;
    cluster_members = cm;
    
}

IMPOPTREP_END_NAMESPACE
