/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChainVariation.cpp
 * Author: mwittig
 * 
 * Created on July 22, 2019, 11:06 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>

#include <regex>

#include "meinetools.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CIsbtGtAllele.h"
#include "CVariantChainVariation.h"

CVariantChainVariation::CVariantChainVariation() 
{
    
}

CVariantChainVariation::CVariantChainVariation(const CVariantChainVariation& orig) 
{
    first_variant = orig.first_variant;
    second_variant= orig.second_variant;
}

CVariantChainVariation::~CVariantChainVariation() {
}

CVariantChainVariation& CVariantChainVariation::operator =(const CVariantChainVariation& orig)
{
    first_variant = orig.first_variant;
    second_variant= orig.second_variant;
    return *this;
}

bool   CVariantChainVariation::operator <(const CVariantChainVariation& orig)const
{
    if(first_variant < orig.first_variant )
        return true;
    if(first_variant > orig.first_variant)
        return false;
    if(second_variant < orig.second_variant )
        return true;
    return false;
}

bool   CVariantChainVariation::operator ==(const CVariantChainVariation& orig)const
{
    return (
            first_variant == orig.first_variant && 
            second_variant == orig.second_variant
            );
}


