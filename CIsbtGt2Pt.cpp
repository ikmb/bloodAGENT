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

#include "mytools.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CVariantChain.h"
#include "CVariantChains.h"
#include "CIsbtGt2Pt.h"

using namespace std;

CIsbtGt2Pt::CIsbtGt2Pt(const string& filename) 
{
    init(filename);
}

CIsbtGt2Pt::CIsbtGt2Pt(const CIsbtGt2Pt& orig) 
{
    m_allele_vector = orig.m_allele_vector;
    m_typing_results=orig.m_typing_results;
}

CIsbtGt2Pt::~CIsbtGt2Pt() {
}

void CIsbtGt2Pt::init(const string& filename)
{
    CParsedTextfile ptx(filename,"\t","Allele",0,true, "#");
    if(ptx.First())
        do{
            m_allele_vector[ptx["MySystemKey"]].push_back(CIsbtPtAllele(ptx["Allele"], ptx["Phenotype"], ptx["base_change"], ptx["acid_change"], ptx["incidence"]));
        }while(ptx.Next());
    
}

void CIsbtGt2Pt::sort(CIsbtGt2Pt::typing_result& var)
{
    for(auto& gt_scores:var)
    {
        for(auto& act_alleles:gt_scores.second)
        {
            std::sort(act_alleles.second.begin(),act_alleles.second.end(),CIsbtGt2PtHit::sort_by_score_desc);
        }
    }
}

CIsbtGt2Pt::typing_result CIsbtGt2Pt::type(const string& system, const CVariantChains& variants)
{
    map<CIsbtGt,map<CIsbtGtAllele,vector<CIsbtGt2PtHit>>>  mRet;;
    std::set<CIsbtGt> theoretical_genotypes = variants.getPossibleGenotypes(system);
    
    for(auto& possible_sample_genotypes:theoretical_genotypes)
    {
        //cout << "typing " << possible_sample_genotypes << endl; // output the genotype 
        std::set<CIsbtGtAllele> possible_sample_alleles = possible_sample_genotypes.getAlleles();
        //std::set<CIsbtGtAllele>::const_iterator iterSampleAlleles = possible_sample_alleles.begin();
        
        for(auto& possible_sample_allele:possible_sample_alleles)
        {
            vector<CIsbtGt2PtHit> gt2pt =  findMatches(system,possible_sample_allele);
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
                mRet[possible_sample_genotypes][possible_sample_allele]=gt2pt;
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
    for(auto& gt_scores:all_hits)
    {
        for(auto& act_alleles:gt_scores.second)
        {
            for(auto& act_hit:act_alleles.second)
            {
                //range_typed_not_in_anno.first = min(range_typed_not_in_anno.first,act_hit.m_typed_not_in_anno);
                range_typed_not_in_anno.second = std::max(range_typed_not_in_anno.second,act_hit.m_typed_not_in_anno);
                //range_anno_not_in_typed.first = min(range_anno_not_in_typed.first,act_hit.m_anno_not_in_typed);
                range_anno_not_in_typed.second = std::max(range_anno_not_in_typed.second,act_hit.m_anno_not_in_typed);
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

                // harmonic mean
                // weigh: anno_not_in_typed as 1/3 important
                //cout << "ranges typed_not_in_anno: " << range_typed_not_in_anno.first << " - " << range_typed_not_in_anno.second << endl;
                //cout << "ranges anno_not_in_typed: " << range_anno_not_in_typed.first << " - " << range_anno_not_in_typed.second << endl;
                float score = ((2.0f+1.0f)*normed_anno_not_in_typed*normed_typed_not_in_anno); // numerator
                float denominator = ((2.0f*normed_anno_not_in_typed)+(1.0f*normed_typed_not_in_anno));
                if(denominator == 0.0f)
                    score = 0.0f;
                else
                    score /= denominator;
                act_hit.score(score);
                //cout << "score of " << act_hit << "(3.0f+1.0f)*(" <<normed_anno_not_in_typed << '*' << normed_typed_not_in_anno << ")/((" <<normed_anno_not_in_typed << "*3.0f)+(1.0f*"<<normed_typed_not_in_anno << "))" << endl;
                fRet = std::max(fRet,score);
            }
        }
    }
    return fRet;
}

vector<CIsbtGt2PtHit> CIsbtGt2Pt::findMatches(const string& system, const CIsbtGtAllele& isbtGtAllele)
{
    std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iterSys = m_allele_vector.find(system);
    if(iterSys == m_allele_vector.end())
        return vector<CIsbtGt2PtHit>();
    
    //cout << "find matches " << isbtGtAllele << endl;
    vector<CIsbtGt2PtHit> vRet;
    // calculate matching parameters for each annotated allele, 
    for(auto anno:iterSys->second)
    {
        CIsbtGt2PtHit actHit(anno);
        // for each annotated base change 
        for(auto a:anno.baseChanges())
        {
            if(!isbtGtAllele.contains(a))
                actHit.m_anno_not_in_typed++;
        }
        for(auto i:isbtGtAllele.variantSet())
        {
            if(!anno.containsBaseChange(i.name()))
                actHit.m_typed_not_in_anno++;
        }
        vRet.push_back(actHit);
    }
    std::sort(vRet.begin(),vRet.end(),CIsbtGt2PtHit::sort_by_errors_asc);
    return vRet;
}

std::string CIsbtGt2Pt::getStringOfTypingResult(const CIsbtGt& gt,const std::map<CIsbtGtAllele,std::vector<CIsbtGt2PtHit>>& results)const
{
    ostringstream osr("");
    
    osr << gt << "\t";
    int count_outer = 0;
    for(const auto& act_allele:results)
    {
        float high_score = act_allele.second.front().score();
        int count = 0;
        if(count_outer++ != 0)
            osr << '/';
        for(const auto& act_hit:act_allele.second)
        {
            if(act_hit.score() < high_score)
                break;
            if(count++ != 0)
                osr << ';';
            osr << act_hit.m_phenotype_allele.phenotype();
        }
        
    }
    osr << "\t" << getPredictedScoreOfGenotype(results);
    return osr.str();
}


std::string CIsbtGt2Pt::getCallAsString(const std::string& system)const
{
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
            if(getPredictedScoreOfGenotype(act_gt.second) >= top_score*0.99999f)
            {
                if(count++ > 0)
                    osr << endl;
                osr << getStringOfTypingResult(act_gt.first,act_gt.second);
            }
        }
    }
    return osr.str();
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




