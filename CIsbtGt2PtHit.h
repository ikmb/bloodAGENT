/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2PtHit.h
 * Author: mwittig
 *
 * Created on July 25, 2019, 9:42 AM
 */

#ifndef CISBTGT2PTHIT_H
#define CISBTGT2PTHIT_H

#include "CIsbtPtAllele.h"


class CIsbtGt2PtHit {
public:
    CIsbtGt2PtHit(const CIsbtPtAllele& allele);
    CIsbtGt2PtHit(const CIsbtGt2PtHit& orig);
    virtual ~CIsbtGt2PtHit();
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtGt2PtHit& me);
    
    static bool sort_by_errors_asc( const CIsbtGt2PtHit& c1, const CIsbtGt2PtHit& c2 );
    static bool sort_by_score_desc( const CIsbtGt2PtHit& c1, const CIsbtGt2PtHit& c2 );
    
    int errurSum()const{return m_typed_not_in_anno+m_anno_not_in_typed;}
    /**
    * An ISBT allele relevant base change detected, but this one is not relevant for the current allele
    */
    int m_typed_not_in_anno;
    /**
    * An ISBT allele relevant base change detected, but this one is not relevant for the current allele. And it is a high impact SNP
    */
    int m_high_impact_typed_not_in_anno;
    /**
    * An ISBT base change which characterizes this allele is not present
    */
    int m_anno_not_in_typed;
    /**
    * An ISBT base change which characterizes this allele is not present. And it is a high impact SNP
    */
    int m_high_impact_anno_not_in_typed;
    int m_not_covered;  // number of SNPs that are not covered
    int m_high_impact_not_covered;  // number of SNPs that are not covered
    int m_high_impact_match;
    int m_high_impact_mismatch;
    int m_match;
    int m_null_variants;
    CIsbtPtAllele m_phenotype_allele;
    
    
    void score(const float& val){m_score = val;}
    float score()const{return m_score;}
    
private:
    float m_score;

};


#endif /* CISBTGT2PTHIT_H */

