/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChain.cpp
 * Author: mwittig
 * 
 * Created on July 22, 2019, 9:56 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>

#include <regex>

#include "mytools.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CVariantChainVariation.h"
#include "CVariantChain.h"

using namespace std;



CVariantChain::CVariantChain(const CVariantChain& orig) : m_isbt_anno(orig.m_isbt_anno)
{
    m_chains  = orig.m_chains;
}

CVariantChain& CVariantChain::operator=(const CVariantChain& orig)
{
    m_isbt_anno = orig.m_isbt_anno;
    m_chains  = orig.m_chains;
    return *this;
}

CVariantChain::~CVariantChain() 
{
    
}
void CVariantChain::addHA(const CIsbtVariant& var)
{
    CVariantChainVariation vcv;
    vcv.first_variant = var;
    vcv.second_variant = var;
    m_chains["ha"].insert(vcv);
    
}

/*
void CVariantChain::addHR(const CIsbtVariant& var)
{
    m_chains["hr"].insert(var);
}
*/

bool CVariantChain::add(const CVcfSnp& var)
{
    if(!m_isbt_anno)
        return false;
    
    CISBTAnno::variation isbv = m_isbt_anno->getCorrespondingIsbtVariation(var);
    CVariantChainVariation clean_vcv(isbv);
    // remove variant if it is already stored somewhere
    // ususally it should be stored as "hr"
    for(map<string,set<CVariantChainVariation>>::iterator i = m_chains.begin();i != m_chains.end();)
    {
        set<CVariantChainVariation>::iterator j = i->second.find(clean_vcv);
        if( j != i->second.end())
            i->second.erase(j);
        j = i->second.find(CVariantChainVariation());
        if( j != i->second.end())
            i->second.erase(j);
        if(i->second.empty())
        {
            m_chains.erase(i);
            i = m_chains.begin();
        }
        else
            i++;
    }
    
    vector<string> alleles = var.alleles();
    if(isbv.isInDel())
        alleles = var.indelalleles();
    CVariantChainVariation vcv;
    for(int i = 0; i < 2 && i < alleles.size(); i++)
    {
        if(isbv.alternative().compare(alleles[i]) == 0)
            if(i == 0)
                vcv.first_variant = isbv;
            else
                vcv.second_variant=isbv;
    }
    
    if(vcv != CVariantChainVariation())
    {
        if(var.isPhased())
            m_chains[to_string(var.phasingID())].insert(vcv);
        else
            m_chains["no"].insert(vcv);
    }
    
    return true;
}


std::set<CIsbtGt> CVariantChain::getPossibleGenotypes()
{
    std::set<CIsbtGt> sRet;
    map<string,set<CVariantChainVariation>>::iterator i = m_chains.begin();
    
    getPossibleGenotypes(sRet, CIsbtGtAllele(),i);
    getPossibleGenotypes(sRet, CIsbtGtAllele(),i,1);
    
    return sRet;
}

void CVariantChain::getPossibleGenotypes(std::set<CIsbtGt>& vars, CIsbtGtAllele var, 
        map<string,set<CVariantChainVariation>>::iterator iter, int type)
{
    CIsbtGtAllele varCopy = var;
    if(type == 0)
    {
        for(set<CVariantChainVariation>::iterator i = iter->second.begin(); i != iter->second.end(); i++)
            if(i->first_variant != CIsbtVariant())
                var.add(i->first_variant);
    }
    else
    {
        for(set<CVariantChainVariation>::iterator i = iter->second.begin(); i != iter->second.end(); i++)
            if(i->second_variant != CIsbtVariant())
                var.add(i->second_variant);
    }
    iter++;
    if(iter == m_chains.end())
    {
        iter--;
        CIsbtGt newGt;
        newGt.add(var);
        // add second possibility
        for(set<CVariantChainVariation>::iterator i = iter->second.begin(); i != iter->second.end(); i++)
            if(i->second_variant != CIsbtVariant())
                varCopy.add(i->second_variant);
        newGt.add(varCopy);
        vars.insert(newGt);
        return;
    }
    getPossibleGenotypes(vars, var,iter);
    getPossibleGenotypes(vars, var,iter,1);
}


std::ostream& operator<<(std::ostream& os, const CVariantChain& me)
{
    /// key is the phasing chain id with two presets
    /// ha == homozygous for the alternative
    /// hr == homozygous for the reference
    /// no == non phased
    /// everything else is a unique phasing identifier
    for(std::map<std::string,set<CVariantChainVariation>>::const_iterator i = me.m_chains.begin();i != me.m_chains.end() ;i++)
    {
        for(set<CVariantChainVariation>::const_iterator j = i->second.begin(); j != i->second.end(); j++)
        {
            if(i != me.m_chains.begin() || j != i->second.begin())
                os << endl;
            os << i->first << ' ' << j->first_variant << '|' << j->second_variant;
        }
    }
    return os;
}


