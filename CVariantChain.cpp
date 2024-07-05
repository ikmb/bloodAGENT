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


using namespace std;



CVariantChain::CVariantChain(const CVariantChain& orig) : m_isbt_anno(orig.m_isbt_anno)
{
    m_chains  = orig.m_chains;
    std::unique_lock<std::mutex> lock(m_mutex);
    m_maxThreads =orig.m_maxThreads;
    m_activeThreads=orig.m_activeThreads;
}

CVariantChain& CVariantChain::operator=(const CVariantChain& orig)
{
    m_isbt_anno = orig.m_isbt_anno;
    m_chains  = orig.m_chains;
    std::unique_lock<std::mutex> lock(m_mutex);
    m_maxThreads =orig.m_maxThreads;
    m_activeThreads=orig.m_activeThreads;
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

void CVariantChain::removeHA(const CIsbtVariant& var)
{
    CVariantChainVariation vcv;
    vcv.first_variant = var;
    vcv.second_variant = var;
    set<CVariantChainVariation>::iterator i = m_chains["ha"].find(vcv);
    if(i != m_chains["ha"].end())
        m_chains["ha"].erase(i);
}

/*
void CVariantChain::addHR(const CIsbtVariant& var)
{
    m_chains["hr"].insert(var);
}
*/

bool CVariantChain::add(const CVcfSnp& var)
{
    static int unique_unphased_id = 1;
    if(!m_isbt_anno)
        return false;
    
    CISBTAnno::variation isbv = m_isbt_anno->getCorrespondingIsbtVariation(var);
    isbv.addVcfSnp(var);
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
    for(size_t i = 0; i < 2 && i < alleles.size(); i++)
    {
        if(isbv.alternative().compare(alleles[i]) == 0)
        {
            if(i == 0)
                vcv.first_variant = isbv;
            else
                vcv.second_variant=isbv;
        }
    }
    
    if(vcv != CVariantChainVariation())
    {
        if(var.isPhased())
            m_chains[to_string(var.phasingID())].insert(vcv);
        else
        {
            m_chains[string("no_").append(std::to_string(unique_unphased_id++))].insert(vcv);
        }
    }
    
    return true;
}


std::set<CIsbtGt> CVariantChain::getPossibleGenotypes()const
{
    std::set<CIsbtGt> sRet;
    map<string,set<CVariantChainVariation>>::const_iterator i = m_chains.begin();
    
    if(i == m_chains.end())
    {
        CIsbtGt newGt;
        newGt.add(CIsbtGtAllele());
        newGt.add(CIsbtGtAllele());
        sRet.insert(newGt);
    }
    else
    {
        // This declaration and definition is necessary because getPossibleGenotypesMT is not a static function
        auto func = [this](std::set<CIsbtGt>& sRet, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B,
                     std::map<std::string, std::set<CVariantChainVariation>>::const_iterator iter,int type) mutable {
                getPossibleGenotypesMT(sRet, allele_A, allele_B, iter, type);
            };
            
        runInThread(func, sRet, CIsbtGtAllele(), CIsbtGtAllele(),i,0);
        runInThread(func, sRet, CIsbtGtAllele(), CIsbtGtAllele(),i,1);
        // Warten Sie darauf, dass alle Threads beendet sind
        waitForCompletion();

        //getPossibleGenotypes(sRet, CIsbtGtAllele(), CIsbtGtAllele(),i);
        //getPossibleGenotypes(sRet, CIsbtGtAllele(), CIsbtGtAllele(),i,1);
    }
    return sRet;
}

void CVariantChain::getPossibleGenotypesMT(std::set<CIsbtGt>& vars, CIsbtGtAllele allA, CIsbtGtAllele allB,
        map<string,set<CVariantChainVariation>>::const_iterator iter, int type)const
{
    if(type == 0)
    {
        for(set<CVariantChainVariation>::iterator i = iter->second.begin(); i != iter->second.end(); i++)
        {
            if(i->first_variant != CIsbtVariant())
                allA.add(i->first_variant);
            if(i->second_variant != CIsbtVariant())
                allB.add(i->second_variant);
        }
    }
    else
    {
        for(set<CVariantChainVariation>::iterator i = iter->second.begin(); i != iter->second.end(); i++)
        {
            if(i->first_variant != CIsbtVariant())
                allB.add(i->first_variant);
            if(i->second_variant != CIsbtVariant())
                allA.add(i->second_variant);
        }
    }
    iter++;
    if(iter == m_chains.end())
    {
        CIsbtGt newGt;
        newGt.add(allA);
        newGt.add(allB);
        std::lock_guard<std::mutex> lock(m_setMutex);
        vars.insert(newGt);
        return;
    }
    auto func = [this](std::set<CIsbtGt>& sRet, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B,
                     std::map<std::string, std::set<CVariantChainVariation>>::const_iterator iter,int type) mutable {
                getPossibleGenotypesMT(sRet, allele_A, allele_B, iter, type);
            };
    {
        m_mutex.lock();
        if(m_activeThreads < m_maxThreads)
        {
            // Erhöhen Sie die Anzahl der aktiven Threads
            ++m_activeThreads;
            m_mutex.unlock();
            runInThreadRecursive(func, vars, allA, allB,iter);
        }
        else
        {
            m_mutex.unlock();
            getPossibleGenotypesMT(vars, allA, allB,iter);
        }
    }        
    {
        m_mutex.lock();
        if(m_activeThreads < m_maxThreads)
        {
            // Erhöhen Sie die Anzahl der aktiven Threads
            ++m_activeThreads;
            m_mutex.unlock();
            runInThreadRecursive(func, vars, allA, allB,iter,1);
        }
        else
        {
            m_mutex.unlock();
            getPossibleGenotypesMT(vars, allA, allB,iter,1);
        }
    }  
    //*/      
}

void CVariantChain::runInThread(
        std::function<void(std::set<CIsbtGt>&, CIsbtGtAllele, CIsbtGtAllele,
                           std::map<std::string, std::set<CVariantChainVariation>>::const_iterator, int)> func,
        std::set<CIsbtGt>& vars, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B,
        std::map<std::string, std::set<CVariantChainVariation>>::const_iterator iter,
        int type) const 
{ // Hinzugefügtes const für runInThread
    std::unique_lock<std::mutex> lock(m_mutex);

    // Warten, bis die Anzahl der aktiven Threads kleiner als m_maxThreads ist
    m_condition.wait(lock, [this] { return m_activeThreads < m_maxThreads; });

    // Erhöhen Sie die Anzahl der aktiven Threads
    ++m_activeThreads;
    //cout << "threads genotypes " << m_activeThreads << endl;
    // Starten Sie einen neuen Thread, um die Funktion auszuführen
    std::thread([this, func, &vars, allele_A, allele_B, iter, type]() {
        func(vars, allele_A, allele_B, iter, type);

        // Reduzieren Sie die Anzahl der aktiven Threads und benachrichtigen Sie andere Threads
        std::unique_lock<std::mutex> lock(m_mutex);
        --m_activeThreads;
        m_condition.notify_one();
    }).detach();
}

void CVariantChain::runInThreadRecursive(
        std::function<void(std::set<CIsbtGt>&, CIsbtGtAllele, CIsbtGtAllele,
                           std::map<std::string, std::set<CVariantChainVariation>>::const_iterator, int)> func,
        std::set<CIsbtGt>& vars, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B,
        std::map<std::string, std::set<CVariantChainVariation>>::const_iterator iter,
        int type) const 
{ // Hinzugefügtes const für runInThread
    //cout << "threads genotypes " << m_activeThreads << endl;
    // Starten Sie einen neuen Thread, um die Funktion auszuführen
    std::thread([this, func, &vars, allele_A, allele_B, iter, type]() {
        func(vars, allele_A, allele_B, iter, type);

        // Reduzieren Sie die Anzahl der aktiven Threads und benachrichtigen Sie andere Threads
        std::unique_lock<std::mutex> lock(m_mutex);
        --m_activeThreads;
        m_condition.notify_one();
    }).detach();
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


