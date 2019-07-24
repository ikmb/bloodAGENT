/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChain.h
 * Author: mwittig
 *
 * Created on July 22, 2019, 9:56 AM
 */

#ifndef CVARIANTCHAIN_H
#define CVARIANTCHAIN_H

#include "CIsbtGtAllele.h"
#include "CIsbtGt.h"
#include "ISBTAnno.h"
#include "CVcf.h"
#include "CVariantChainVariation.h"


class CVariantChain {
public:
    CVariantChain(): m_isbt_anno(NULL){}
    CVariantChain(CISBTAnno* isbtAnno): m_isbt_anno(isbtAnno){}
    CVariantChain(const CVariantChain& orig);
    CVariantChain& operator=(const CVariantChain& orig);
    virtual ~CVariantChain();
    
    friend std::ostream& operator<<(std::ostream& os, const CVariantChain& me);
    
    /// this function is usually used to add ISBT variant for genomic positions
    /// where LRG sequence is != reference sequence (e.g del261G in ABO as the del is the hg19 reference)
    /// these SNPs are always homozyguous for the alternative allele
    void addHA(const CIsbtVariant& var);
    //void addHR(const CIsbtVariant& var);
    bool add(const CVcfSnp& var);
    
    std::set<CIsbtGt> getPossibleGenotypes();
    
private:
    
    void getPossibleGenotypes(std::set<CIsbtGt>& vars, CIsbtGtAllele var, map<string,set<CVariantChainVariation>>::iterator iter, int type = 0);
    
    CISBTAnno* m_isbt_anno;
    
    /// key is the phasing chain id with two presets
    /// ha == homozygous for the alternative
    /// hr == homozygous for the reference
    /// no == non phased
    /// everything else is a unique phasing identifier
    std::map<std::string,set<CVariantChainVariation>> m_chains;

    
};

#endif /* CVARIANTCHAIN_H */

