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
#include <iomanip>

#include "CIsbtGt2PtHit.h"

CIsbtGt2PtHit::CIsbtGt2PtHit(const CIsbtPtAllele& allele) : m_phenotype_allele(allele)
{
    m_typed_not_in_anno = 0;
    m_anno_not_in_typed = 0;
    m_not_covered = 0;
    m_high_impact_match=0;
    m_high_impact_mismatch=0;
    m_score = 0.0f;
    m_high_impact_not_covered = 0;
    m_high_impact_anno_not_in_typed = 0;
    m_high_impact_typed_not_in_anno = 0;
    m_match = 0;
}

CIsbtGt2PtHit::CIsbtGt2PtHit(const CIsbtGt2PtHit& orig)  : m_phenotype_allele(orig.m_phenotype_allele)
{
    m_typed_not_in_anno = orig.m_typed_not_in_anno;
    m_anno_not_in_typed = orig.m_anno_not_in_typed;
    m_high_impact_match=orig.m_high_impact_match;
    m_high_impact_mismatch=orig.m_high_impact_match;
    m_not_covered = orig.m_not_covered;
    m_score = orig.m_score;
    m_high_impact_not_covered = orig.m_high_impact_not_covered;
    m_high_impact_anno_not_in_typed = orig.m_high_impact_anno_not_in_typed;
    m_high_impact_typed_not_in_anno = orig.m_high_impact_typed_not_in_anno;
    m_match = orig.m_match;
}

CIsbtGt2PtHit::~CIsbtGt2PtHit() 
{
    
}

bool CIsbtGt2PtHit::sort_by_errors_asc( const CIsbtGt2PtHit& c1, const CIsbtGt2PtHit& c2 ) 
{ 
    if(c1.m_high_impact_mismatch < c2.m_high_impact_mismatch)
        return true;
    if(c1.m_high_impact_mismatch > c2.m_high_impact_mismatch)
        return false;
    if(c1.m_high_impact_typed_not_in_anno < c2.m_high_impact_typed_not_in_anno)
        return true;
    if(c1.m_high_impact_typed_not_in_anno > c2.m_high_impact_typed_not_in_anno)
        return false;
    if(c1.m_typed_not_in_anno < c2.m_typed_not_in_anno)
        return true;
    if(c1.m_typed_not_in_anno > c2.m_typed_not_in_anno)
        return false;
    if(c1.m_high_impact_anno_not_in_typed < c2.m_high_impact_anno_not_in_typed)
        return true;
    if(c1.m_high_impact_anno_not_in_typed > c2.m_high_impact_anno_not_in_typed)
        return false;
    if(c1.m_anno_not_in_typed < c2.m_anno_not_in_typed)
        return true;
    if(c1.m_anno_not_in_typed > c2.m_anno_not_in_typed)
        return false;
    if(c1.m_high_impact_not_covered < c2.m_high_impact_not_covered)
        return true;
    if(c1.m_high_impact_not_covered > c2.m_high_impact_not_covered)
        return false;
    return c1.m_not_covered < c2.m_not_covered; 
}
bool CIsbtGt2PtHit::sort_by_score_desc( const CIsbtGt2PtHit& c1, const CIsbtGt2PtHit& c2 )
{
    return c1.score() > c2.score();
}

std::ostream& operator<<(std::ostream& os, const CIsbtGt2PtHit& me)
{
    os << "score: " << std::fixed << std::setprecision(5)  << me.m_score << " for " << me.m_phenotype_allele << " e1: " << me.m_anno_in_typed_but_not_in_current_genotype << " e2: " << me.m_anno_not_in_typed << " e3: " << me.m_typed_not_in_anno;
    return os;
}


