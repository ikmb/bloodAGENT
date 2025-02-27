/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChains.h
 * Author: mwittig
 *
 * Created on July 24, 2019, 7:35 AM
 */

#ifndef CVARIANTCHAINS_H
#define CVARIANTCHAINS_H

#include <string>

#include "CIsbtGt.h"


class CVariantChains {
public:
    CVariantChains(int maxThreads=8);
    CVariantChains(CISBTAnno* isbt,int maxThreads=8);
    CVariantChains(const CVariantChains& orig);
    virtual ~CVariantChains();
    
    friend std::ostream& operator<<(std::ostream& os, const CVariantChains& me);
    
    bool init();
    std::string add(const CVcfSnp& act_snp);
    
    void setBreakPhasingVariable(bool v){m_break_phasing=v;}
    
    std::set<CIsbtGt> getPossibleGenotypes(const string& system)const;
    
    // this is required in our test data generator. Avoid calling this method!!!
    void removeReferenceSnps();
    void removeUncoveredSnps(double limit, int verbose = 0);
    
    const CISBTAnno* isbtSnps()const{return m_isbt;}
  
private:
    
    /// These are differences between LRG and genome reference
    void addReferenceSnps();

    CISBTAnno* m_isbt;
    int m_maxThreads;
    map<string,CVariantChain> m_variant_chains;
    bool m_break_phasing;
    
};

#endif /* CVARIANTCHAINS_H */

