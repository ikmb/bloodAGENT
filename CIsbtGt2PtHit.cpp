/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2PtHit.cpp
 * Author: mwittig
 * 
 * Created on July 25, 2019, 9:42 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>
#include <iostream>

#include "CIsbtGt2PtHit.h"

CIsbtGt2PtHit::CIsbtGt2PtHit(const CIsbtPtAllele& allele) : m_phenotype_allele(allele)
{
    m_typed_not_in_anno = 0;
    m_anno_not_in_typed = 0;
}

CIsbtGt2PtHit::CIsbtGt2PtHit(const CIsbtGt2PtHit& orig)  : m_phenotype_allele(orig.m_phenotype_allele)
{
    m_typed_not_in_anno = orig.m_typed_not_in_anno;
    m_anno_not_in_typed = orig.m_anno_not_in_typed;
}

CIsbtGt2PtHit::~CIsbtGt2PtHit() 
{
    
}

bool CIsbtGt2PtHit::sort_by_errors_asc( const CIsbtGt2PtHit& c1, const CIsbtGt2PtHit& c2 ) 
{ 
    return c1.m_anno_not_in_typed+c1.m_typed_not_in_anno < c2.m_anno_not_in_typed+c2.m_typed_not_in_anno; 
}

std::ostream& operator<<(std::ostream& os, const CIsbtGt2PtHit& me)
{
    os << me.m_phenotype_allele << " e1: " << me.m_anno_not_in_typed << " e2: " << me.m_typed_not_in_anno;
}


