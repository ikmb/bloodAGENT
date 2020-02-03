/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtPtAllele.cpp
 * Author: mwittig
 * 
 * Created on July 25, 2019, 8:18 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>

#include "mytools.h"

#include "CIsbtPtAllele.h"

using namespace std;

CIsbtPtAllele::CIsbtPtAllele()
{
    init("","","","","",0.0f);
}

CIsbtPtAllele::CIsbtPtAllele(std::string name, std::string phenotype,const std::string& flat_phenotype_name, std::string base_changes, std::string acid_changes, string incidence) 
{
    try{
        init(name,phenotype,flat_phenotype_name,base_changes,acid_changes,stof(incidence));
    }catch(...)
    {
        init(name,phenotype,flat_phenotype_name,base_changes,acid_changes,-1.0f);
    }
}

CIsbtPtAllele::CIsbtPtAllele(std::string name, std::string phenotype,const std::string& flat_phenotype_name, std::string base_changes, std::string acid_changes, float incidence) 
{
    init(name,phenotype,flat_phenotype_name,base_changes,acid_changes,incidence);
}

CIsbtPtAllele::CIsbtPtAllele(const CIsbtPtAllele& orig) 
{
    m_name = orig.m_name;
    m_phenotype_name = orig.m_phenotype_name;
    m_base_changes = orig.m_base_changes; // ISBT naming like 261delG or 85A>G, ...
    m_acid_changes = orig.m_acid_changes; // ISBT naming like Ile60Leu or Ser230Ile, ...
    m_incidence = orig.m_incidence;
    m_flat_phenotype_name = orig.m_flat_phenotype_name;
}

CIsbtPtAllele::~CIsbtPtAllele() 
{
    
}
void CIsbtPtAllele::init(const std::string& name, const std::string& phenotype,const std::string& flat_phenotype_name, const std::string& base_changes, const std::string& acid_changes, float incidence)
{
    m_name = name;
    m_phenotype_name = phenotype;
    m_flat_phenotype_name = flat_phenotype_name;
    m_incidence = incidence;
    
    vector<string> parsed = CMyTools::GetParsedLine(base_changes," ");
    for(auto x:parsed)
        m_base_changes.insert(x);
    parsed = CMyTools::GetParsedLine(acid_changes," ");
    for(auto x:parsed)
        m_acid_changes.insert(x);
     
}

bool CIsbtPtAllele::containsBaseChange(const std::string& isbt_base_change)const
{
    for(auto x:m_base_changes)
        if(x.compare(isbt_base_change) == 0)
            return true;
    return false;
}

std::vector<std::string> CIsbtPtAllele::getFullBaseChangeRecombinations()const
{
    std::vector<std::string> vRet;
    for(set<string>::iterator i = m_base_changes.begin(); i != m_base_changes.end();i++)
    {
       string actSet = "";
        for(set<string>::iterator j = i; j != m_base_changes.end();j++)
        {
            actSet.append(*j).append(" ");
            vRet.push_back(actSet);
        }
    }
    return vRet;
}

bool CIsbtPtAllele::operator==(const CIsbtPtAllele& me)const
{
    return (
            me.name().compare(m_name) == 0 &&
            me.m_phenotype_name.compare(m_phenotype_name) == 0 &&
            me.m_incidence == m_incidence &&
            me.m_base_changes == m_base_changes &&
            me.m_acid_changes == m_acid_changes
            );
}

bool CIsbtPtAllele::operator<(const CIsbtPtAllele& me)const
{
    if(me.name().compare(m_name) < 0)
        return true;
    if(me.name().compare(m_name) > 0)
        return false;
    if(me.m_phenotype_name.compare(m_phenotype_name) < 0)
        return true;
    if(me.m_phenotype_name.compare(m_phenotype_name) > 0)
        return false;
    if(me.m_incidence < m_incidence)
        return true;
    if(me.m_incidence > m_incidence)
        return false;
    if(me.m_base_changes.size() < m_base_changes.size())
        return true;
    if(me.m_base_changes.size() > m_base_changes.size())
        return false;
    if(me.m_acid_changes.size() < m_acid_changes.size())
        return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const CIsbtPtAllele& me)
{
    long unsigned int i = 0;
    os << me.m_name << '\t';
    for(auto x:me.m_base_changes)
    {
        os << x << ( ++i == me.m_base_changes.size() ? "" : " ");
    }
    os << '\t';
    i = 0;
    for(auto x:me.m_acid_changes)
    {
        os << x << ( ++i == me.m_acid_changes.size() ? "" : " ");
    }
    os << '\t' << me.m_incidence;
    return os;
}





