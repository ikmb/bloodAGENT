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
#include "ISBTAnno.h"

class CIsbtGt2Pt {
public:
    CIsbtGt2Pt(const std::string& filename,int maxThreads=8);
    CIsbtGt2Pt(const CIsbtGt2Pt& orig);
    virtual ~CIsbtGt2Pt();
    
    typedef std::map<CIsbtGt,std::multimap<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>> typing_result;
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtGt2Pt& me);
    
    std::vector<CIsbtGt2PtHit> findMatches(const std::string& system, const CIsbtGtAllele& IsbtGt, const CISBTAnno* isbt_snps, int required_coverage);
    typing_result type(const string& system, const CVariantChains& variants, int required_coverage = 10, float score_range = 1.0f);
    void doTheMatching(const std::string& system,CIsbtGt2Pt::typing_result& mRet, const CVariantChains& variants,set<CIsbtGt>::const_iterator  possible_sample_genotype, int required_coverage, float& highest_score, float score_range);
    void doCleaning(CIsbtGt2Pt::typing_result& mRet, float highest_score, float score_range);
    
    void sort(typing_result& var);
    
    /// returns all best calls within a specific range
    /// \param system: the blood group system
    /// \param phenotype: report phenotye or allele (default is allele (false)))
    /// \param top_score_range: we multiply this value with the best score and report all calls >= this call
    /// \return one call per line
    std::string getCallAsString(const CISBTAnno& isbt_anno, const std::string& system, bool phenotype = true, float top_score_range = 0.999f, const std::string& sampleId = "")const;
    nlohmann::json getCallAsJson(const CISBTAnno& isbt_anno, const CTranscriptAnno& trans_anno, const CBigWigReader& bwr, const std::string& system, bool phenotype = true, float top_score_range = 0.999f, int coverage_limit = 10)const;
    
    
    vector<CIsbtPtAllele> alleleVector(const string& system)const;
    /// this returns the CIsbtPtAllele of a given allele
    /// @param1: alleel name as string. Eg ABO*B.01
    CIsbtPtAllele alleleOf(const string& allele)const;
    string systemOf(const string& allele)const;
    
    
private:
    
    void init(const std::string& filename);
    void scoreHits(CIsbtGt2Pt::typing_result&, const string& system,const CISBTAnno* isbt_anno);
    void scoreHit(CIsbtGt2PtHit&, const string& system,const CISBTAnno* isbt_anno);
    float getPredictedScoreOfGenotype(const std::multimap<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& allele_calls)const;
    float getTopPredictedScoreOfAllGenotypes(const typing_result& genotype_calls)const;
    static bool sort_by_space_separated_entries_asc(const string& a,const string& b);
    void findAlleTaggingBaseChanges();
    
    
    /// set phenotype to true to get the phenotype, otherwise you get the allele
    //std::string getStringOfTypingResult(const CIsbtGt& gt,const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results, bool phenotype = true)const;
    nlohmann::json getJsonOfTypingResult(const CIsbtGt& gt,const std::multimap<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results)const;
    
    std::map<std::string,vector<CIsbtPtAllele>> m_allele_vector;
    std::map<std::string,vector<CIsbtPtAllele>> m_allele_vector_redundant;
    
    std::map<std::string,typing_result> m_typing_results;

    
    
    void runInThread(
        std::function<void(const std::string& ,CIsbtGt2Pt::typing_result& , const CVariantChains& , 
                           set<CIsbtGt>::const_iterator ,int , float& , float )> func,
        const std::string& system,CIsbtGt2Pt::typing_result& mRet, const CVariantChains& variants, 
        set<CIsbtGt>::const_iterator  possible_sample_genotypes,int required_coverage, float& highest_score, float score_range) const ;
    
    void waitForCompletion()const {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_condition.wait(lock, [this] { return m_activeThreads == 0; });
    }
    
    mutable int m_maxThreads;
    mutable int m_activeThreads;
    mutable std::mutex m_mutex;
    mutable std::mutex m_objectMutex;
    mutable std::condition_variable m_condition;
};

#endif /* CISBTGT2PT_H */

