

/**
 *  \file ChiSquareTestResult.h
 *  \brief Store p-value and Cramer's V from chi-square test
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPOPTREP_CHISQUARE_TEST_RESULT_H
#define IMPOPTREP_CHISQUARE_TEST_RESULT_H
#include <IMP.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include <IMP/object_macros.h>
#include "optrep_config.h"

IMPOPTREP_BEGIN_NAMESPACE

class IMPOPTREPEXPORT ChiSquareTestResult {
 public:
 
ChiSquareTestResult(Float pval, Float cramersv); 

// storage
Float pvalue;
Float cramersv;


IMP_OBJECT_METHODS(ChiSquareTestResult);

};

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_CHISQUARE_TEST_RESULT_H */
