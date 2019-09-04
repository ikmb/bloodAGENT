/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2Pt.h
 * Author: mwittig
 *
 * Created on July 25, 2019, 7:00 AM
 */

#ifndef CISBTGT2PT_H
#define CISBTGT2PT_H

#include "bigWig.h"
#include "CBigWigReader.h"
#include "CIsbtPtAllele.h"
#include "CIsbtGtAllele.h"
#include "CIsbtGt2PtHit.h"

class CIsbtGt2Pt {
public:
    CIsbtGt2Pt(const std::string& filename);
    CIsbtGt2Pt(const CIsbtGt2Pt& orig);
    virtual ~CIsbtGt2Pt();
    
    typedef std::map<CIsbtGt,std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>> typing_result;
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtGt2Pt& me);
    
    std::vector<CIsbtGt2PtHit> findMatches(const std::string& system, const CIsbtGtAllele& IsbtGt);
    typing_result type(const string& system, const CVariantChains& variants);
    
    void sort(typing_result& var);
    
    /// returns all best calls within a specific range
    /// \param system: the blood group system
    /// \param phenotype: report phenotye or allele (default is allele (false)))
    /// \param top_score_range: we multiply this value with the best score and report all calls >= this call
    /// \return one call per line
    std::string getCallAsString(const std::string& system, bool phenotype = true, float top_score_range = 0.999f)const;
    
    vector<CIsbtPtAllele> alleleVector(const string& system)const;
    
private:
    
    void init(const std::string& filename);
    float scoreHits(std::map<CIsbtGt,std::map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>&);
    float getPredictedScoreOfGenotype(const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& allele_calls)const;
    float getTopPredictedScoreOfAllGenotypes(const typing_result& genotype_calls)const;
    
    
    /// set phenotype to true to get the phenotype, otherwise you get the allele
    std::string getStringOfTypingResult(const CIsbtGt& gt,const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results, bool phenotype = true)const;
    
    std::map<std::string,vector<CIsbtPtAllele>> m_allele_vector;
    
    std::map<std::string,typing_result> m_typing_results;
    

};

#endif /* CISBTGT2PT_H */

