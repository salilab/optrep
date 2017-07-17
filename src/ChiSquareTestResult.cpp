/**
*  \file ChiSquareTestResult.cpp
*  \brief Store result of chi-square test 
*
*  Copyright 2007-2017 IMP Inventors. All rights reserved.
*
*/
#include <IMP/algebra.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include <iostream>
#include <string>
#include <IMP/rmf.h>
#include <IMP/optrep/ChiSquareTestResult.h>

IMPOPTREP_BEGIN_NAMESPACE

ChiSquareTestResult::ChiSquareTestResult(Float pv, Float cv) {
    
    pvalue = pv;
    cramersv = cv;
    
}

IMPOPTREP_END_NAMESPACE
