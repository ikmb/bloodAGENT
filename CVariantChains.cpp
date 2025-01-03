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
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <regex>

#include "mytools.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "CIsbtGtAllele.h"
#include "CIsbtGt.h"
#include "ISBTAnno.h"
#include "CVariantChainVariation.h"
#include "CVariantChain.h"



#include "CVariantChains.h"

using namespace std;

CVariantChains::CVariantChains(int maxThreads) : m_maxThreads(maxThreads)
{
    m_isbt = NULL;
    m_break_phasing = false;
}

CVariantChains::CVariantChains(CISBTAnno* isbt,int maxThreads) : m_isbt(isbt), m_maxThreads(maxThreads)
{
    init();
    m_break_phasing = false;
}

CVariantChains::CVariantChains(const CVariantChains& orig) : m_isbt(orig.m_isbt),m_maxThreads(orig.m_maxThreads)
{
    m_variant_chains = orig.m_variant_chains;
    m_break_phasing = orig.m_break_phasing;
}

CVariantChains::~CVariantChains() 
{
    
}


bool CVariantChains::init()
{
    if(!m_isbt)
        throw CMyException("Can not init CVariantChains object without an CISBTAnno object");
    
    set<string> loci = m_isbt->loci();
    for(set<string>::iterator locusIter = loci.begin(); locusIter != loci.end(); locusIter++)
        m_variant_chains[*locusIter]=CVariantChain(m_isbt,m_maxThreads);
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

void CVariantChains::removeUncoveredSnps(double limit, int verbose)
{
    for(map<string,CVariantChain>::iterator i = m_variant_chains.begin(); i != m_variant_chains.end(); i++)
    {
        for(std::map<std::string,set<CVariantChainVariation>>::iterator iChains =  i->second.getChains().begin(); iChains !=  i->second.getChains().end(); iChains++)
        {
           for(set<CVariantChainVariation>::iterator iVchainVar = iChains->second.begin(); iVchainVar != iChains->second.end();)
           {
               if( !iVchainVar->first_variant.isCovered(limit) && !iVchainVar->second_variant.isCovered(limit) && limit > 0)
               {
                   if(verbose >= 1)
                       cerr << i->first << ": deleting too low covered variant " << iVchainVar->first_variant << '/' << iVchainVar->second_variant << endl;
                   iVchainVar = iChains->second.erase(iVchainVar);
               }
               else
                   iVchainVar++;
           }
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


string CVariantChains::add(const CVcfSnp& act_snp)
{
    if(!m_isbt)
        return "";
    
    string the_act_system = m_isbt->getSystemAt(act_snp.chrom(),act_snp.pos());
    if(!the_act_system.empty())
    {
        m_variant_chains[the_act_system].add(act_snp,m_break_phasing);
        return the_act_system;
    }
    return "";
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



