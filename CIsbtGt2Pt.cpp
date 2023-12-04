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
#include <libgen.h>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>

#include "mytools.h"
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

using namespace std;

CIsbtGt2Pt::CIsbtGt2Pt(const string& filename) 
{
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


CIsbtGt2Pt::typing_result CIsbtGt2Pt::type(const string& system, const CVariantChains& variants, int required_coverage)
{
    map<CIsbtGt,map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>  mRet;
    std::set<CIsbtGt> theoretical_genotypes = variants.getPossibleGenotypes(system);
    // !!!!!!!!!!!!!!!!!!!!!
    // Hier alle Allele Des systems holen und dann gegen alle theoretical_genotypes abgleichen
    // Und zwar direkt unten im For loop
    // Am beste ich übergebde der Klasse ISBTAnno Objekte von CIsbtPtAllele und er wandelt diese objekte
    // in ein set oder Vector von CIsbtVariant um. Dann habe ich auch alle Infos ob high impact variant
    // dann kann ich das scoring laufen lassen
    // go through all heterozygous genetoypes
    for(set<CIsbtGt>::const_iterator possible_sample_genotypes = theoretical_genotypes.begin(); possible_sample_genotypes != theoretical_genotypes.end(); possible_sample_genotypes++)
    {
        //cout << "typing " << *possible_sample_genotypes << endl; // output the genotype 
        std::set<CIsbtGtAllele> possible_sample_alleles = possible_sample_genotypes->getAlleles();
        // !!!!!!!!!!!! hier //////////////////
        //std::set<CIsbtGtAllele>::const_iterator iterSampleAlleles = possible_sample_alleles.begin();
        // evaluate each allele if it fits to the current genotype
        for(const CIsbtGtAllele& possible_sample_allele:possible_sample_alleles)
        {
            vector<CIsbtGt2PtHit> gt2pt =  findMatches(system,possible_sample_allele,variants.isbtSnps(),required_coverage);
            //for(auto i : gt2pt)
            //    cout << i << endl;
            /*if(
                gt2pt.size() != 0 && iterSampleAlleles->variantCount() < gt2pt[0].errurSum()
              )
            {
                //cout << possible_sample_genotypes << " has to many errors in " << gt2pt[0] << endl;
                //act_hits.clear();
                break;
            }
            else*/
            if(gt2pt.size() == 0)
            {
                //cout << possible_sample_genotypes << " has no result " << endl;
                //act_hits.clear();
                break;
            }
            else
            {
                mRet[*possible_sample_genotypes][possible_sample_allele]=gt2pt;
                //cout << gt2pt[0] << " -------- " << mRet[*possible_sample_genotypes][possible_sample_allele][0] << endl;
            }
        }
        
    }
    scoreHits(mRet);
    sort(mRet);
    m_typing_results[system]=mRet;
    return mRet;
}

float CIsbtGt2Pt::scoreHits(map<CIsbtGt,map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>& all_hits)
{
    float fRet = 1.0f;
    pair<int,int> range_typed_not_in_anno=pair<int,int>(0,-1); // this is for normalizing values
    pair<int,int> range_anno_not_in_typed=pair<int,int>(0,-1); // this is for normalizing values
    pair<int,int> range_high_impact_matches=pair<int,int>(0,-1); // this is for normalizing values
    pair<int,int> range_high_impact_mismatches=pair<int,int>(0,-1); // this is for normalizing values
    pair<int,int> range_anno_in_typed_but_not_in_current_genotype=pair<int,int>(0,-1); // this is for normalizing values
    pair<int,int> range_not_covered=pair<int,int>(0,-1); // this is for normalizing values
    // map<CIsbtGt,map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>
    for(auto& gt_scores:all_hits)
    {
        // map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>
        for(auto& act_alleles:gt_scores.second)
        {
            // vector<CIsbtGt2PtHit>
            for(auto& act_hit:act_alleles.second)
            {
                //range_typed_not_in_anno.first = min(range_typed_not_in_anno.first,act_hit.m_typed_not_in_anno);
                range_typed_not_in_anno.second = std::max(range_typed_not_in_anno.second,act_hit.m_typed_not_in_anno);
                //range_anno_not_in_typed.first = min(range_anno_not_in_typed.first,act_hit.m_anno_not_in_typed);
                range_anno_not_in_typed.second = std::max(range_anno_not_in_typed.second,act_hit.m_anno_not_in_typed);
                range_anno_in_typed_but_not_in_current_genotype.second = std::max(range_anno_in_typed_but_not_in_current_genotype.second,act_hit.m_anno_in_typed_but_not_in_current_genotype);
                range_not_covered.second = std::max(range_not_covered.second,act_hit.m_not_covered);
                range_high_impact_matches.second =std::max(range_high_impact_matches.second,act_hit.m_high_impact_match);
                range_high_impact_mismatches.second =std::max(range_high_impact_mismatches.second,act_hit.m_high_impact_mismatch);
            }
        }
    }
    if(range_typed_not_in_anno.first == -1)
        return fRet;
    for(auto& gt_scores:all_hits)
    {
        for(auto& act_alleles:gt_scores.second)
        {
            for(auto& act_hit:act_alleles.second)
            {
                float normed_typed_not_in_anno = 1.0f;
                if(range_typed_not_in_anno.first != range_typed_not_in_anno.second)
                    normed_typed_not_in_anno = 1.0f-static_cast<float>(act_hit.m_typed_not_in_anno-range_typed_not_in_anno.first)/static_cast<float>(range_typed_not_in_anno.second-range_typed_not_in_anno.first);
                float normed_anno_not_in_typed = 1.0f;
                if(range_anno_not_in_typed.first != range_anno_not_in_typed.second)
                    normed_anno_not_in_typed = 1.0f-static_cast<float>(act_hit.m_anno_not_in_typed-range_anno_not_in_typed.first)/static_cast<float>(range_anno_not_in_typed.second-range_anno_not_in_typed.first);
                float normed_anno_in_typed_but_not_in_current_genotype = 1.0f;
                if(range_anno_in_typed_but_not_in_current_genotype.first != range_anno_in_typed_but_not_in_current_genotype.second)
                    normed_anno_in_typed_but_not_in_current_genotype = 1.0f-static_cast<float>(act_hit.m_anno_in_typed_but_not_in_current_genotype-range_anno_in_typed_but_not_in_current_genotype.first)/static_cast<float>(range_anno_in_typed_but_not_in_current_genotype.second-range_anno_in_typed_but_not_in_current_genotype.first);
                float normed_not_covered = 1.0f;
                if(range_not_covered.first != range_not_covered.second)
                    normed_not_covered = 1.0f-static_cast<float>(act_hit.m_not_covered-range_not_covered.first)/static_cast<float>(range_not_covered.second-range_not_covered.first);
                
                [[maybe_unused]] float normed_high_impact_matches = 1.0f;
                if(range_high_impact_matches.first != range_high_impact_matches.second)
                    normed_high_impact_matches = static_cast<float>(act_hit.m_high_impact_match-range_high_impact_matches.first)/static_cast<float>(range_high_impact_matches.second-range_high_impact_matches.first);
                [[maybe_unused]] float normed_high_impact_mismatches = 1.0f;
                if(range_high_impact_mismatches.first != range_high_impact_mismatches.second)
                    normed_high_impact_mismatches = 1.0f-static_cast<float>(act_hit.m_high_impact_mismatch-range_high_impact_mismatches.first)/static_cast<float>(range_high_impact_mismatches.second-range_high_impact_mismatches.first);

                // harmonic mean
                // weigh: anno_not_in_typed as 1/3 important
                //cout << "ranges typed_not_in_anno: " << range_typed_not_in_anno.first << " - " << range_typed_not_in_anno.second << endl;
                //cout << "ranges anno_not_in_typed: " << range_anno_not_in_typed.first << " - " << range_anno_not_in_typed.second << endl;
                float score = ((2.0f+1.0f+5.0f+1.0f)*normed_anno_not_in_typed*normed_typed_not_in_anno*normed_anno_in_typed_but_not_in_current_genotype*normed_not_covered); // numerator
                float denominator = ((2.0f*normed_anno_not_in_typed)+(1.0f*normed_typed_not_in_anno)+(5.0f*normed_anno_in_typed_but_not_in_current_genotype)+(1.0f*normed_not_covered));
                if(denominator == 0.0f)
                    score = 0.0f;
                else
                    score /= denominator;
                // here comes the impact of m_high_impact_match
                // pull up the score that for those that have high impact matches
                if(range_high_impact_matches.first != range_high_impact_matches.second)
                {
                    for(int i = 0; i < act_hit.m_high_impact_match; i++)
                        score = (score + 1.0f)/2.0f;
                }
                
                act_hit.score(score);
                //cout << "score of " << act_hit << endl;
                fRet = std::max(fRet,score);
            }
        }
    }
    return fRet;
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
            bool isHighImpactSnp = isbt_snps->getIsbtVariant(system,a).isHighImpactSNP();
            double act_variant_coverage = isbt_snps->getIsbtVariant(system,a).getCoverage();
            if(static_cast<int>(act_variant_coverage+0.5) < required_coverage)
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
                    actHit.m_high_impact_anno_not_in_typed++;
                else
                    actHit.m_anno_not_in_typed++;
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
                if(isHighImpactSnp)
                    actHit.m_high_impact_typed_not_in_anno++;
                else
                    actHit.m_typed_not_in_anno++;
            }
        }
        //cout << actHit << endl;
        vRet.push_back(actHit);
    }
    std::sort(vRet.begin(),vRet.end(),CIsbtGt2PtHit::sort_by_errors_asc);
    return vRet;
}

nlohmann::json CIsbtGt2Pt::getJsonOfTypingResult(const CIsbtGt& gt,const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results)const
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
            metrics["high_impact_snp_matches"] = act_hit.m_high_impact_match;
            metrics["typed_not_in_anno"] = act_hit.m_typed_not_in_anno;
            metrics["anno_not_in_typed"] = act_hit.m_anno_not_in_typed;
            metrics["anno_in_typed_but_not_in_current_genotype"] = act_hit.m_anno_in_typed_but_not_in_current_genotype; // this is a strong indicator for a false positive, as it is well covered, typed but not in this genotype
            metrics["not_covered"] = act_hit.m_not_covered;
            allele["issues"].push_back(metrics);
            phenotype.push_back(act_hit.m_phenotype_allele.phenotype());
            flat_phenotype.push_back(act_hit.m_phenotype_allele.flatPhenotype());
        }
        alleles.push_back(allele);
        phenotypes.push_back(phenotype);
        flat_phenotypes.push_back(flat_phenotype);
    }
    jRet["alleles"]=alleles;
    jRet["phenotypes"]=phenotypes;
    jRet["flat_phenotypes"]=flat_phenotypes;
    jRet["score"]=getPredictedScoreOfGenotype(results);
    jRet["wobble_score"]=jRet["score"];
    
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
        // special RhD treatment, trans_anno only has key if parameter Trick was set
        if(system.compare("RHD") == 0 && trans_anno.hasKey("RHD"))
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
                js["wobble_score"]=0.0f;
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
                js["wobble_score"]=2.0f;
                j["calls"].push_back(js);
                type_by_snps = false;
            }
            else if( rhd_cov/rhce_cov <= 0.66 )
            {
                is_RHD_DEL_HET = true;
            }
            //cout << "RHD\t- & -\tRhD-/RhD-\t2\t-" << endl;
        }
        if(type_by_snps)
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
                    if(system.compare("RHD") == 0 && is_RHD_DEL_HET && jAct["alleles"].size() == 1)
                    {
                        nlohmann::json allele;
                        allele["names"].push_back("RHD*01N.01");
                        jAct["alleles"].push_back(allele);
                    }
                    j["calls"].push_back(jAct);
                }
            }
        }
    }
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

float CIsbtGt2Pt::getPredictedScoreOfGenotype(const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& allele_calls)const
{
    float fRet = 0.0f;
    for(auto& act_allele:allele_calls)
    {
        if(!act_allele.second.empty())
            fRet+=act_allele.second.front().score();
    }
    if(allele_calls.size()==1) // is homozygous?
        fRet+=fRet;
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




