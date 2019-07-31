/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChains.cpp
 * Author: mwittig
 * 
 * Created on July 24, 2019, 7:35 AM
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

#include "CVariantChains.h"

using namespace std;

CVariantChains::CVariantChains()
{
    m_isbt = NULL;
}

CVariantChains::CVariantChains(CISBTAnno* isbt) : m_isbt(isbt)
{
    init();
}

CVariantChains::CVariantChains(const CVariantChains& orig) : m_isbt(orig.m_isbt)
{
    m_variant_chains = orig.m_variant_chains;
}

CVariantChains::~CVariantChains() 
{
    
}


bool CVariantChains::init()
{
    if(!m_isbt)
        return false;
    
    set<string> loci = m_isbt->loci();
    for(set<string>::iterator locusIter = loci.begin(); locusIter != loci.end(); locusIter++)
        m_variant_chains[*locusIter]=CVariantChain(m_isbt);
    addReferenceSnps();
    return true;
}

void CVariantChains::addReferenceSnps()
{
    // add SNPs where LRG and hg19 differ
    map<string,vector<CISBTAnno::variation> > refVar = m_isbt->getReferenceVariations();
    for(map<string,vector<CISBTAnno::variation> >::const_iterator i = refVar.begin(); i != refVar.end(); i++)
    {
        for(vector<CISBTAnno::variation>::const_iterator j = i->second.begin(); j != i->second.end(); j++)
        {
            m_variant_chains[i->first].addHA(*j);
            
        }
    }
}

void CVariantChains::removeReferenceSnps()
{
    map<string,vector<CISBTAnno::variation> > refVar = m_isbt->getReferenceVariations();
    for(map<string,vector<CISBTAnno::variation> >::const_iterator i = refVar.begin(); i != refVar.end(); i++)
    {
        for(vector<CISBTAnno::variation>::const_iterator j = i->second.begin(); j != i->second.end(); j++)
        {
            m_variant_chains[i->first].removeHA(*j);
            
        }
    }
}


bool CVariantChains::add(const CVcfSnp& act_snp)
{
    if(!m_isbt)
        return false;
    
    string the_act_system = m_isbt->getSystemAt(act_snp.chrom(),act_snp.pos());
    if(!the_act_system.empty())
    {
        m_variant_chains[the_act_system].add(act_snp);
        return true;
    }
    return false;
}

std::set<CIsbtGt> CVariantChains::getPossibleGenotypes(const string& system)const
{
    std::set<CIsbtGt> sRet;
    map<string,CVariantChain>::const_iterator i = m_variant_chains.find(system);
    if(i != m_variant_chains.end())
    {
        return i->second.getPossibleGenotypes();
    }
    return sRet;
}


std::ostream& operator<<(std::ostream& os, const CVariantChains& me)
{
    /// key is the phasing chain id with two presets
    /// ha == homozygous for the alternative
    /// hr == homozygous for the reference
    /// no == non phased
    /// everything else is a unique phasing identifier
    if(me.m_isbt)
    {
        set<string> loci = me.m_isbt->loci();
        for(set<string>::iterator locusIter = loci.begin(); locusIter != loci.end(); locusIter++)
        {
            os << "----------" << endl << *locusIter << endl;
            map<string,CVariantChain>::const_iterator i = me.m_variant_chains.find(*locusIter);
            if(i != me.m_variant_chains.end())
                os << i->second << endl;

        }
    }  
    return os;
}



