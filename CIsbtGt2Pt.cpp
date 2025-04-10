/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2Pt.cpp
 * Author: mwittig
 * 
 * Created on July 25, 2019, 7:00 AM
 */

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <iomanip>
#include <libgen.h>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>

#include "meinetools.h"
#include "CBigWigReader.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtGt.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CVariantChainVariation.h"
#include "CVariantChain.h"
#include "CVariantChains.h"
#include "json/single_include/nlohmann/json.hpp"
#include "CTranscript.h"
#include "CTranscriptAnno.h"
#include "CIsbtGt2Pt.h"

#define VAL_COSINE_REF_ALLELE -1.0f
#define VAL_COSINE_ALT_ALLELE 1.0f
#define VAL_COSINE_COV_FAILED_ALLELE 0.25f
#define VAL_COSINE_STANDARD_WEIGHT 1.0f
#define VAL_COSINE_HIGH_IMPACT_WEIGHT 2.0f

using namespace std;

CIsbtGt2Pt::CIsbtGt2Pt(const string& filename,int maxThreads) 
{
    m_maxThreads = maxThreads;
    m_activeThreads=0;
    if(!CMyTools::file_exists(filename))
        throw(CMyException("File does not exist: ")+filename);
    init(filename);
    findAlleTaggingBaseChanges();
}

CIsbtGt2Pt::CIsbtGt2Pt(const CIsbtGt2Pt& orig) 
{
    m_allele_vector = orig.m_allele_vector;
    m_allele_vector_redundant = orig.m_allele_vector_redundant;
    m_typing_results=orig.m_typing_results;
    std::unique_lock<std::mutex> lock(m_mutex);
    m_maxThreads =orig.m_maxThreads;
    m_activeThreads=orig.m_activeThreads;
}

CIsbtGt2Pt::~CIsbtGt2Pt() {
}

void CIsbtGt2Pt::init(const string& filename)
{
    CParsedTextfile ptx(filename,"\t","Allele",0,true, "#");
    if(ptx.First())
        do{
            m_allele_vector[ptx["MySystemKey"]].push_back(CIsbtPtAllele(ptx["Allele"], ptx["Phenotype"], ptx["Phenotype_flat"], ptx["base_change"], ptx["acid_change"], ptx["incidence"]));
        }while(ptx.Next());
    
}

void CIsbtGt2Pt::sort(CIsbtGt2Pt::typing_result& var)
{
    for(CIsbtGt2Pt::typing_result::iterator gt_scores = var.begin(); gt_scores != var.end(); gt_scores++)
    {
        for(std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>::iterator act_alleles = gt_scores->second.begin(); act_alleles != gt_scores->second.end(); act_alleles++)
        {
            std::sort(act_alleles->second.begin(),act_alleles->second.end(),CIsbtGt2PtHit::sort_by_score_desc);
        }
    }
}

void CIsbtGt2Pt::doCleaning(CIsbtGt2Pt::typing_result& mRet, float highest_score, float score_range)
{
    waitForCompletion();
    std::unique_lock<std::mutex> lock(m_mutex);
    for(CIsbtGt2Pt::typing_result::iterator i = mRet.begin(); i != mRet.end(); )
    {
        float act_score = getPredictedScoreOfGenotype(i->second);
        if(act_score < highest_score*score_range)
            i= mRet.erase(i);
        else
            i++;
    }
}

void CIsbtGt2Pt::doTheMatching(const std::string& system,CIsbtGt2Pt::typing_result& mRet, const CVariantChains& variants,
        set<CIsbtGt>::const_iterator  possible_sample_genotype, int required_coverage, float& highest_score, float score_range)
{
    //cout << "typing " << *possible_sample_genotype << endl; // output the genotype 
    std::multiset<CIsbtGtAllele> possible_sample_alleles = possible_sample_genotype->getAlleles();
    // !!!!!!!!!!!! hier //////////////////
    //std::set<CIsbtGtAllele>::const_iterator iterSampleAlleles = possible_sample_alleles.begin();
    // evaluate each allele if it fits to the current genotype
    for(const CIsbtGtAllele& possible_sample_allele:possible_sample_alleles)
    {
        //vector<CIsbtGt2PtHit> gt2pt = findMatches(system,possible_sample_allele,variants.isbtSnps(),required_coverage);
        vector<CIsbtGt2PtHit> gt2pt= cosineSimilarityMatches(system,possible_sample_allele,variants.isbtSnps(),required_coverage);
        //cout << possible_sample_alleles.size() << " with " << possible_sample_allele << " " << gt2pt.size() << endl;
        if(gt2pt.size() == 0)
        {
            //cout << *possible_sample_genotype << " has no result " << endl;
            //act_hits.clear();
            continue;
        }
        else
        {
            //cout << *possible_sample_genotype << " has " << gt2pt.size() << endl;
            // map<CIsbtGt,multimap<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>
            std::lock_guard<std::mutex> lock(m_objectMutex);
            typing_result::iterator iterA = mRet.find(*possible_sample_genotype);
            if(iterA == mRet.end())
            {
                multimap<CIsbtGtAllele,vector<CIsbtGt2PtHit>> mTmp;
                iterA = mRet.insert(pair<CIsbtGt,multimap<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>(*possible_sample_genotype,mTmp)).first;
            }
            iterA->second.insert(pair<CIsbtGtAllele,vector<CIsbtGt2PtHit>>(possible_sample_allele,gt2pt));
            /*
            size_t a = mRet.size();
            size_t b = iterA->second.size();
            size_t c = iterA->second.find(possible_sample_allele)->second.size();
            cout << a << " * " << b << " * " << c << " = " << a*b*c << endl;
            *///cout << gt2pt[0] << " -------- " << mRet[*possible_sample_genotype][possible_sample_allele][0] << endl;
        }
    }
    
    std::lock_guard<std::mutex> lock(m_objectMutex);
    typing_result::iterator iterA = mRet.find(*possible_sample_genotype);
    if(iterA != mRet.end())
    {
        float act_score = getPredictedScoreOfGenotype(iterA->second);
        highest_score = std::max(act_score,highest_score);
        
    }
                
}



CIsbtGt2Pt::typing_result CIsbtGt2Pt::type(const string& system, const CVariantChains& variants, int required_coverage, float score_range)
{
    // map<CIsbtGt,multimap<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>  mRet;
    CIsbtGt2Pt::typing_result mRet;
    std::set<CIsbtGt> theoretical_genotypes = variants.getPossibleGenotypes(system);
    float highest_score = 0.0f;
    int cleaning_counter=0;
    //cout << "typing " << m_activeThreads << endl;
    
    // !!!!!!!!!!!!!!!!!!!!!
    // Hier alle Allele Des systems holen und dann gegen alle theoretical_genotypes abgleichen
    // Und zwar direkt unten im For loop
    // Am beste ich übergebde der Klasse ISBTAnno Objekte von CIsbtPtAllele und er wandelt diese objekte
    // in ein set oder Vector von CIsbtVariant um. Dann habe ich auch alle Infos ob high impact variant
    // dann kann ich das scoring laufen lassen
    // go through all heterozygous genetoypes
     
    for(set<CIsbtGt>::const_iterator possible_sample_genotypes = theoretical_genotypes.begin(); possible_sample_genotypes != theoretical_genotypes.end(); possible_sample_genotypes++)
    {
        
        auto func = [this](const std::string& system,CIsbtGt2Pt::typing_result& mRet, const CVariantChains& variants, 
                           set<CIsbtGt>::const_iterator  possible_sample_genotypes, 
                           int required_coverage, float& highest_score, float score_range) mutable {
                doTheMatching(system,mRet,variants,possible_sample_genotypes,
                      required_coverage,highest_score,score_range);
            };
        runInThread(func, system,mRet,variants,possible_sample_genotypes,
                      required_coverage,highest_score,score_range);
        //
        bool clean_up = false;
        {
            std::lock_guard<std::mutex> lock(m_objectMutex);
            clean_up = ++cleaning_counter%m_maxThreads == 0;
        }
        if(clean_up)
        {
            doCleaning(mRet,highest_score,score_range);
        }
        //*/
        //doTheMatching(system,mRet,variants,possible_sample_genotypes,required_coverage,highest_score,score_range);
    }
    doCleaning(mRet,highest_score,score_range);
    m_typing_results[system]=mRet;
    return mRet;
}

    
void CIsbtGt2Pt::scoreCosineSimilarity(CIsbtGt2PtHit& act_hit,const vector<float>& typedSNV, const vector<float>& insilicoSNV, const vector<float>& weights)
{
    float dotProduct = 0.0;
    float normA = 0.0;
    float normB = 0.0;
    
    for (size_t i = 0; i < typedSNV.size(); ++i) {
        dotProduct += typedSNV[i] * insilicoSNV[i] * weights[i]; // Gewichte einbeziehen
        normA += typedSNV[i] * typedSNV[i] * weights[i]; // Gewichtete Norm
        normB += insilicoSNV[i] * insilicoSNV[i] * weights[i];
    }

    // Vermeidung der Division durch Null
    if (normA == 0 || normB == 0) 
        act_hit.score(0.0f);
    else
        act_hit.score(dotProduct / (sqrt(normA) * sqrt(normB)));
}

void CIsbtGt2Pt::scoreHit(CIsbtGt2PtHit& act_hit, const string& system,const CISBTAnno* isbt_anno)
{
    // NOT In Use. We use scoreCosineSimilarity
    float system_var_count = static_cast<float>(isbt_anno->getIsbtVariantCount(system));
    float act_hit_malus = 0.0f;
    
    /**
    * An ISBT allele relevant base change detected, but this one is not relevant for the current allele
    -2 m_typed_not_in_anno;
    * An ISBT allele relevant base change detected, but this one is not relevant for the current allele. And it is a high impact SNP
    -4 m_high_impact_typed_not_in_anno;
    * An ISBT base change which characterizes this allele but is not present in the sample
    -2 m_anno_not_in_typed;
    * An ISBT base change which characterizes this allele but is not present in the sample. And it is a high impact SNP
    -4 m_high_impact_anno_not_in_typed;
    -1 m_not_covered;  // number of SNPs that are not covered
    -4 m_high_impact_not_covered;  // number of SNPs that are not covered
    +2 m_high_impact_match;
    -4 m_high_impact_mismatch;
    int m_null_variants;
    */
    
    system_var_count += act_hit.m_high_impact_match;
    act_hit_malus+= act_hit.m_typed_not_in_anno.size()*2.0f +
                    act_hit.m_high_impact_typed_not_in_anno.size()*4.0f +
                    act_hit.m_anno_not_in_typed.size()*2.0f +
                    act_hit.m_high_impact_anno_not_in_typed.size() * 4.0f +
                    act_hit.m_not_covered*1.0f +
                    act_hit.m_high_impact_not_covered*4.0f + 
                    act_hit.m_high_impact_mismatch.size()*4.0f;

    if(act_hit_malus > system_var_count)
        act_hit.score(0.0f);
    else
        act_hit.score(0.5f/system_var_count*(system_var_count-act_hit_malus));
    //cout << "score of " << act_hit << endl;
}

void CIsbtGt2Pt::scoreHits(CIsbtGt2Pt::typing_result& all_hits, const string& system,const CISBTAnno* isbt_anno)
{
    for(auto& gt_scores:all_hits)
    {
        // multimap<CIsbtGtAllele,vector<CIsbtGt2PtHit>>
        for(auto& act_alleles:gt_scores.second)
        {
            // vector<CIsbtGt2PtHit>
            for(auto& act_hit:act_alleles.second)
            {
                scoreHit(act_hit, system,isbt_anno);
             }
        }
    }
}

vector<CIsbtGt2PtHit> CIsbtGt2Pt::findMatches(const string& system, const CIsbtGtAllele& isbtGtAllele, const CISBTAnno* isbt_snps, int required_coverage)
{
    std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iterSys = m_allele_vector.find(system);
    if(iterSys == m_allele_vector.end())
        return vector<CIsbtGt2PtHit>();
    //std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iterSys = m_allele_vector_redundant.find(system);
    //if(iterSys == m_allele_vector_redundant.end())
    //    return vector<CIsbtGt2PtHit>();
    
    //cout << "find matches " << isbtGtAllele << endl;
    vector<CIsbtGt2PtHit> vRet;
    // calculate matching parameters for each annotated allele, 
    for(const CIsbtPtAllele& anno:iterSys->second)
    {
        CIsbtGt2PtHit actHit(anno);
        // for each annotated base change 
        for(const string& a:anno.baseChanges())
        {
            if(a.size() == 0)
                continue;
            CISBTAnno::variation var = isbt_snps->getIsbtVariant(system,a);
            bool isHighImpactSnp = var.isHighImpactSNP();
            double act_variant_coverage = isbt_snps->getIsbtVariant(system,a).getCoverage();
            if(static_cast<int>(act_variant_coverage) < required_coverage)
            {
                // incomplete covered
                if(isHighImpactSnp)
                    actHit.m_high_impact_not_covered++;
                else
                    actHit.m_not_covered++;
            }
            if( !isbtGtAllele.contains(a) )
            {
                if(isHighImpactSnp)
                    actHit.m_high_impact_anno_not_in_typed.insert(var);
                else
                    actHit.m_anno_not_in_typed.insert(var);
            }
            else 
            {
                if( isHighImpactSnp )
                    actHit.m_high_impact_match++;
                else
                    actHit.m_match++;
            }
        }
        // Go through ISBT annotation
        for(const CIsbtVariant& i:isbtGtAllele.variantSet())
        {
            if(!anno.containsBaseChange(i.name()))
            {
                if(i.isHighImpactSNP())
                    actHit.m_high_impact_typed_not_in_anno.insert(i);
                else
                    actHit.m_typed_not_in_anno.insert(i);
            }
        }
        //cout << actHit << endl;
        scoreHit(actHit, system,isbt_snps);
        vRet.push_back(actHit);
    }
    std::sort(vRet.begin(),vRet.end(),CIsbtGt2PtHit::sort_by_score_desc);
    return vRet;
}

void CIsbtGt2Pt::outPutCosineSim(const string& allelea,const CIsbtGtAllele& isbtGtAllele, vector<float> a, vector<float> b, vector<float> c)
{
    std::unique_lock<std::mutex> lock(m_debugMutex);
    cout << "Haplotype: " << isbtGtAllele << endl;
    cout << "Allele:    " << allelea << endl;
    for (const float &element : a) {
        std::cout << std::setprecision(1) << std::setw(3) << element;
    }
    std::cout << std::endl;
    for (const float &element : b) {
        std::cout << std::setprecision(1) << std::setw(3) << element;
    }
    std::cout << std::endl;
    for (const float &element : c) {
        std::cout << std::setprecision(1) << std::setw(3) << element;
    }
    std::cout << std::endl;
}

vector<CIsbtGt2PtHit> CIsbtGt2Pt::cosineSimilarityMatches(const string system, const CIsbtGtAllele& isbtGtAllele, const CISBTAnno* isbt_snps, int required_coverage)
{
    std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iterSys = m_allele_vector.find(system);
    if(iterSys == m_allele_vector.end())
        return vector<CIsbtGt2PtHit>();
    vector<CIsbtGt2PtHit> vRet;
    
    
     // this is a possible Genotype combination derived from the VCF file
    std::set<CIsbtVariant>  potential_haplotype = isbtGtAllele.variantSet();
    // This is the genotype_to_phenotype_annotation from the ISBT. 
    // ToDo	ABO	ABO	ABO*O.01.01	O	O	O	261delG	Thr88Profs*31	28.41%
    
    // this is a list of all variation from this system
    vector<CISBTAnno::variation> allSystemVariations= isbt_snps->getAllVariations(system);
    map<std::string,int> allSystemVariationsIdx;
    int idx = 0;
    for(CISBTAnno::variation& var:allSystemVariations)
        allSystemVariationsIdx[var.name()]=idx++;
    int systemVarCount = allSystemVariations.size();
    
    
    // calculate matching parameters for each annotated allele, 
    for(const CIsbtPtAllele& allele_specific_SNV:iterSys->second)
    {
        // set all SNVs to reference allele and all weights to 1.0f
        vector<float> typedSNV(systemVarCount, VAL_COSINE_REF_ALLELE);
        vector<float> insilicoSNV(systemVarCount, VAL_COSINE_REF_ALLELE);
        vector<float> weights(systemVarCount, VAL_COSINE_STANDARD_WEIGHT);
    
        CIsbtGt2PtHit actHit(allele_specific_SNV);
        // for each annotated base change 
        for(const string& a:allele_specific_SNV.baseChanges())
        {
            if(a.size() == 0)
                continue;
            map<std::string,int>::iterator i = allSystemVariationsIdx.find(a);
            if(i != allSystemVariationsIdx.end())
                insilicoSNV[i->second] = VAL_COSINE_ALT_ALLELE;
        }    
        // in silico build allele
        for(const CIsbtVariant& a:potential_haplotype)
        {
            map<std::string,int>::iterator i = allSystemVariationsIdx.find(a.name());
            if(i != allSystemVariationsIdx.end())
                typedSNV[i->second] = VAL_COSINE_ALT_ALLELE;
        }  
        int idx = 0;
        for(CISBTAnno::variation& var:allSystemVariations)
        {
            if(var.isHighImpactSNP())
                weights[idx]=VAL_COSINE_HIGH_IMPACT_WEIGHT;
            if(static_cast<int>(var.getCoverage()) < required_coverage)
            {
                //typedSNV[idx]=insilicoSNV[idx]=VAL_COSINE_COV_FAILED_ALLELE;
                if(var.isHighImpactSNP())
                    actHit.m_high_impact_not_covered++;
                else
                {
                    //cout << var << endl;
                    actHit.m_not_covered++;
                }
            }
            // now we count the issues
            if(typedSNV[idx]!=insilicoSNV[idx])
            {
                // an ISBT SNP which is not detected in the sample but described for this potential allele
                if(typedSNV[idx] == VAL_COSINE_REF_ALLELE)
                {
                    if(var.isHighImpactSNP())
                        actHit.m_high_impact_anno_not_in_typed.insert(var);
                    else
                        actHit.m_anno_not_in_typed.insert(var);
                }
                // an ISBT SNP which is detected in the sample but not a SNP that is needed for this potential allele
                else if(insilicoSNV[idx] == VAL_COSINE_REF_ALLELE)
                {
                    if(var.isHighImpactSNP())
                        actHit.m_high_impact_typed_not_in_anno.insert(var);
                    else
                        actHit.m_typed_not_in_anno.insert(var);
                }
            }
            else if(typedSNV[idx] !=  VAL_COSINE_COV_FAILED_ALLELE)
            {
                if(var.isHighImpactSNP())
                    actHit.m_high_impact_match++;
                else
                    actHit.m_match++;
            }
            idx++;
        }
        scoreCosineSimilarity(actHit, typedSNV,insilicoSNV,weights);
        vRet.push_back(actHit);
        //outPutCosineSim(allele_specific_SNV.name(),isbtGtAllele,typedSNV,insilicoSNV,weights);
    }
    std::sort(vRet.begin(),vRet.end(),CIsbtGt2PtHit::sort_by_score_desc);
    return vRet;
}

nlohmann::json CIsbtGt2Pt::getJsonOfTypingResult(const CIsbtGt& gt,const std::multimap<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results, bool homozygous_only)const
{
    nlohmann::json jRet;
    
    nlohmann::json haplotypes;
    for(auto act_allele : gt.getAlleles())
    {
        nlohmann::json genotypes;
        for(auto act_variant :act_allele.variantSet())
        {
            nlohmann::json genotype = act_variant.getSnpAsJson();
            genotypes.push_back(genotype);
        }
        haplotypes["genotypes"].push_back(genotypes);
        if(homozygous_only) // used for RhD for example, when I discover heterozygous RhD deletion from coverage analysis, I expect only a single second allele
            break;
    }
    jRet["haplotypes"]=haplotypes;
    nlohmann::json alleles;
    nlohmann::json phenotypes;
    nlohmann::json flat_phenotypes;
    for(const auto& act_allele:results)
    {
        float high_score = act_allele.second.front().score();
        nlohmann::json allele;
        nlohmann::json phenotype;
        nlohmann::json flat_phenotype;
        for(const auto& act_hit:act_allele.second)
        {
            if(act_hit.score() < high_score)
                break;
            allele["names"].push_back(act_hit.m_phenotype_allele.name());
            nlohmann::json metrics;
            metrics["typed_not_in_anno_count"] = act_hit.m_typed_not_in_anno.size();
            for(CISBTAnno::variation var :act_hit.m_typed_not_in_anno)
                metrics["typed_not_in_anno"].push_back(var.getSnpAsJson());
            metrics["anno_not_in_typed_count"] = act_hit.m_anno_not_in_typed.size();
            for(CISBTAnno::variation var :act_hit.m_anno_not_in_typed)
                metrics["anno_not_in_typed"].push_back(var.getSnpAsJson());
            metrics["high_impact_snp_matches"] = act_hit.m_high_impact_match;
            metrics["high_impact_mismatch_count"] = act_hit.m_high_impact_mismatch.size();
            for(CISBTAnno::variation var :act_hit.m_high_impact_mismatch)
                metrics["high_impact_mismatch"].push_back(var.getSnpAsJson());
            metrics["high_impact_typed_not_in_anno_count"] = act_hit.m_high_impact_typed_not_in_anno.size();
            for(CISBTAnno::variation var :act_hit.m_high_impact_typed_not_in_anno)
                metrics["high_impact_typed_not_in_anno"].push_back(var.getSnpAsJson());
            metrics["high_impact_anno_not_in_typed_count"] = act_hit.m_high_impact_anno_not_in_typed.size();
            for(CISBTAnno::variation var :act_hit.m_high_impact_anno_not_in_typed)
                metrics["high_impact_anno_not_in_typed"].push_back(var.getSnpAsJson());
            metrics["not_covered"] = act_hit.m_not_covered;
            metrics["high_impact_not_covered_count"] = act_hit.m_high_impact_not_covered;
            metrics["null_variants_count"] = act_hit.m_null_variants;
            allele["issues"].push_back(metrics);
            phenotype.push_back(act_hit.m_phenotype_allele.phenotype());
            flat_phenotype.push_back(act_hit.m_phenotype_allele.flatPhenotype());
        }
        alleles.push_back(allele);
        phenotypes.push_back(phenotype);
        flat_phenotypes.push_back(flat_phenotype);
        if(homozygous_only) // used for RhD for example, when I discover heterozygous RhD deletion from coverage analysis, I expect only a single second allele
            break;
    }
    jRet["alleles"]=alleles;
    jRet["phenotypes"]=phenotypes;
    jRet["flat_phenotypes"]=flat_phenotypes;
    
    double val = getPredictedScoreOfGenotype(results);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << val;
    jRet["score"]=stream.str();
    jRet["weak_score"]=jRet["score"];
    
    return jRet;
}

nlohmann::json CIsbtGt2Pt::getCallAsJson(const CISBTAnno& isbt_anno, const CTranscriptAnno& trans_anno, const CBigWigReader& bwr, const std::string& system, bool phenotype, float top_score_range, int coverage_limit)const
{
    nlohmann::json j;
    j["system"]=system;
    nlohmann::json uncovered_target_variants_list;
    std::vector<CISBTAnno::variation> vl = isbt_anno.getCoverageFailedVariants(system);
    for(auto a : vl)
        uncovered_target_variants_list.push_back(a.name());
    j["coverage_failed_variants"]=uncovered_target_variants_list;
    
    vector<double> cov = trans_anno.getCoverages(system,bwr);
    if(!cov.empty())
    {
        j["mean_coverage"]["cds"]=cov[0];
        for(size_t i = 1; i < cov.size();i++)
            j["mean_coverage"]["exons"].push_back(cov[i]);
    }
    
    std::vector<CISBTAnno::variation> all_variations = isbt_anno.getAllVariations(system);
    nlohmann::json all_target_variants_list;
    for(auto a : all_variations)
    {
        nlohmann::json act_j = a.getSnpAsJson();
        act_j["coverage_limit"] = coverage_limit;
        act_j["is_covered"] = a.isCovered(static_cast<double>(coverage_limit));
        all_target_variants_list.push_back(act_j);
    }
    j["relevant_variations"]=all_target_variants_list;
    
    std::map<std::string,typing_result>::const_iterator iRes = m_typing_results.find(system);
    if(iRes != m_typing_results.end())
    {
        bool type_by_snps = true; // for example RHD: if coverage is 0 we do not type and set this to false
        bool is_RHD_DEL_HET = false;
        // special RhD treatment
        if(system.compare("RHD") == 0 && trans_anno.hasKey("RHD") && !phenotype)
        {
            double rhd_cov  = trans_anno.getExonicCoverage("RHD",bwr);
            double rhce_cov = trans_anno.getExonicCoverage("RHCE",bwr);
            if(rhce_cov == 0.0)
            {
                nlohmann::json js;
                nlohmann::json allele;
                allele["names"].push_back("n.a.");
                js["alleles"].push_back(allele);
                js["phenotypes"].push_back("n.a.");
                js["flat_phenotypes"].push_back("n.a.");
                js["score"]=0.0f;
                js["weak_score"]=0.0f;
                j["calls"].push_back(js);
                type_by_snps = false;
            }
            else if( rhd_cov/rhce_cov <= 0.1 )
            {
                nlohmann::json js;
                nlohmann::json allele;
                allele["names"].push_back("RHD*01N.01");
                js["alleles"].push_back(allele);
                js["phenotypes"].push_back("RHD*01N.01");
                js["flat_phenotypes"].push_back("RHD*01N.01");
                js["score"]=2.0f;
                js["weak_score"]=2.0f;
                j["calls"].push_back(js);
                type_by_snps = false;
            }
            else if( rhd_cov/rhce_cov <= 0.66 )
            {
                is_RHD_DEL_HET = true;
            }
            //cout << "RHD\t- & -\tRhD-/RhD-\t2\t-" << endl;
        }
        if(type_by_snps) // we skip that if we did RhD typing by coverage nd detected dd or a complete coverage fail for RhDCE
        {
            // std::map<CIsbtGt,std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>>
            const CIsbtGt2Pt::typing_result& typing = iRes->second;
            float top_score = getTopPredictedScoreOfAllGenotypes(iRes->second);
            for(auto& act_gt:typing)
            {
                if(getPredictedScoreOfGenotype(act_gt.second) >= top_score*top_score_range)
                {
                    nlohmann::json jAct = getJsonOfTypingResult(act_gt.first,act_gt.second);
                    if(!uncovered_target_variants_list.empty())
                        jAct["score"] = 0.0;
                    //if(system.compare("RHD") == 0 && is_RHD_DEL_HET && jAct["alleles"].size() == 1)
                    if(system.compare("RHD") == 0 && is_RHD_DEL_HET)
                    {
                        nlohmann::json allele;
                        jAct = getJsonOfTypingResult(act_gt.first,act_gt.second,true);
                        allele["names"].push_back("RHD*01N.01");
                        jAct["alleles"].push_back(allele);
                    }
                    j["calls"].push_back(jAct);
                }
            }
        }
    }
    //cout << j << endl;
    return j;
}

std::string CIsbtGt2Pt::getCallAsString(const CISBTAnno& isbt_anno, const std::string& system, bool phenotype, float top_score_range, const std::string& sampleId)const
{
    ostringstream uncovered_target_variants_list("");
    /*
    if(isbt_anno.hasUncoveredVariants(system))
    {
        std::vector<CISBTAnno::variation> vl = isbt_anno.getCoverageFailedVariants(system);
        uncovered_target_variants = vl.size();
        for(auto a : vl)
        {
            if(uncovered_target_variants_list.str().empty())
                uncovered_target_variants_list << a.name();
            else
                uncovered_target_variants_list << ',' << a.name();
        }
    }
    ostringstream osr("");
    std::map<std::string,typing_result>::const_iterator iRes = m_typing_results.find(system);
    if(iRes != m_typing_results.end())
    {
        // std::map<CIsbtGt,std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>>
        const CIsbtGt2Pt::typing_result& typing = iRes->second;
        float top_score = getTopPredictedScoreOfAllGenotypes(iRes->second);
        int count = 0;
        for(auto& act_gt:typing)
        {
            if(getPredictedScoreOfGenotype(act_gt.second) >= top_score*top_score_range)
            {
                if(count++ > 0)
                    osr << endl;
                osr << (sampleId.empty() ? "" : sampleId+"\t") << system << '\t' << getStringOfTypingResult(act_gt.first,act_gt.second,phenotype) << '\t' << uncovered_target_variants << '/' << isbt_anno.getAllVariations(system).size() << '\t' << uncovered_target_variants_list.str();
            }
        }
    }
     */
    return uncovered_target_variants_list.str();
}

float CIsbtGt2Pt::getTopPredictedScoreOfAllGenotypes(const typing_result& genotype_calls)const
{
    float fRet = 0.0f;
    
    for(const auto& act_gt:genotype_calls)
    {
        float act_score = getPredictedScoreOfGenotype(act_gt.second);
        if(act_score > fRet)
            fRet= act_score;
    }
    
    return fRet;
}

void CIsbtGt2Pt::findAlleTaggingBaseChanges()
{
    // std::map<std::string,vector<CIsbtPtAllele>>
    for(auto blood_system : m_allele_vector)
    {
        map<string,set<CIsbtPtAllele>> unique_Finder;
        const string& act_system = blood_system.first;
        for(auto act_allele : blood_system.second)
        {
            std::vector<string> act_all_combinations = act_allele.getFullBaseChangeRecombinations();
            for(auto act_combination : act_all_combinations)
                unique_Finder[act_combination].insert(act_allele);
        }
        // find tagging SNPs
        for(auto act_allele : blood_system.second)
        {
            vector<string> hit_list;
            for(auto act_gt_combi : unique_Finder)
            {   // if the act_allele is the only one for this genotype set then store it
                // multiple hits are possible, so we will take the one with the fewest variations
                if(act_gt_combi.second.size() == 1 && act_gt_combi.second.begin()->operator ==(act_allele))
                    hit_list.push_back(act_gt_combi.first);
            }
            std::sort(hit_list.begin(),hit_list.end(),sort_by_space_separated_entries_asc);
            
            for(auto act_gt : hit_list)
            {
                //cerr << act_system << '\t' << act_allele.name() << '\t' << act_gt << endl;
                m_allele_vector_redundant[act_system].push_back(CIsbtPtAllele(act_allele.name(), act_allele.phenotype(), act_allele.flatPhenotype(), act_gt, "", 0.0f));
            }
        }
    }
}

bool CIsbtGt2Pt::sort_by_space_separated_entries_asc(const string& a,const string& b)
{
    size_t countA = count(a.begin(), a.end(), ' ');
    size_t countB = count(b.begin(), b.end(), ' ');
    
    if(countA < countB)
        return true;
    return false;
    
}

float CIsbtGt2Pt::getPredictedScoreOfGenotype(const std::multimap<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& allele_calls)const
{
    float fRet = 0.0f;
    for(auto& act_allele:allele_calls)
    {
        if(!act_allele.second.empty())
            fRet+=getPredictedScoreOfAllele(act_allele);
    }
    if(allele_calls.size()==1) // is homozygous?
        fRet+=fRet;
    return fRet;
}

float CIsbtGt2Pt::getPredictedScoreOfAllele(const std::pair<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& allele)const
{
    float fRet = 0.0f;
    if(!allele.second.empty()){
        fRet+=allele.second.front().score();
    }
    return fRet;
}

vector<CIsbtPtAllele> CIsbtGt2Pt::alleleVector(const string& system)const
{
    std::map<std::string,vector<CIsbtPtAllele>>::const_iterator i = m_allele_vector.find(system);
    if(i != m_allele_vector.end())
        return i->second;
    return vector<CIsbtPtAllele>();
}

CIsbtPtAllele CIsbtGt2Pt::alleleOf(const string& allele)const
{
    for(map<string,vector<CIsbtPtAllele>>::const_iterator i = m_allele_vector.begin(); i != m_allele_vector.end(); i++)
    {
        for(auto isbtptallele:i->second)
            if(isbtptallele.name().compare(allele) == 0)
                return isbtptallele;
        
    }
    return CIsbtPtAllele();
}

string CIsbtGt2Pt::systemOf(const string& allele)const
{
    for(map<string,vector<CIsbtPtAllele>>::const_iterator i = m_allele_vector.begin(); i != m_allele_vector.end(); i++)
    {
        for(auto isbtptallele:i->second)
            if(isbtptallele.name().compare(allele) == 0)
                return i->first;
        
    }
    return "";
}

void CIsbtGt2Pt::runInThread(
        std::function<void(const std::string& ,CIsbtGt2Pt::typing_result& , const CVariantChains& , 
                           set<CIsbtGt>::const_iterator ,int , float& , float )> func,
        const std::string& system,CIsbtGt2Pt::typing_result& mRet, const CVariantChains& variants, 
        set<CIsbtGt>::const_iterator  possible_sample_genotypes,int required_coverage, float& highest_score, float score_range) const
{
    std::unique_lock<std::mutex> lock(m_mutex);

    // Warten, bis die Anzahl der aktiven Threads kleiner als m_maxThreads ist
    m_condition.wait(lock, [this] { return m_activeThreads < m_maxThreads; });

    // Erhöhen Sie die Anzahl der aktiven Threads
    ++m_activeThreads;
    //cout << "threads matches " << m_activeThreads << endl;
    // Starten Sie einen neuen Thread, um die Funktion auszuführen
    std::thread([this, func, &system,&mRet, &variants,possible_sample_genotypes, 
        required_coverage,&highest_score,score_range]() {
        func(system,mRet,variants,possible_sample_genotypes, 
        required_coverage,highest_score,score_range);
            // Reduzieren Sie die Anzahl der aktiven Threads und benachrichtigen Sie andere Threads
        std::unique_lock<std::mutex> lock(m_mutex);
        --m_activeThreads;
        m_condition.notify_one();
    }).detach();
}

std::ostream& operator<<(std::ostream& os, const CIsbtGt2Pt& me)
{
    long unsigned int i = 0;
    for(std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iter = me.m_allele_vector.begin(); iter != me.m_allele_vector.end(); iter++)
    {
        for(auto x:iter->second)
        {
            os <<  iter->first << '\t' << x;
            if( ++i != iter->second.size())
                os << endl;
        }
    }
    return os;
}




