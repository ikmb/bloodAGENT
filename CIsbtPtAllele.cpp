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

CIsbtPtAllele::CIsbtPtAllele(std::string name, std::string base_changes, std::string acid_changes, string incidence) 
{
    try{
        init(name,base_changes,acid_changes,stof(incidence));
    }catch(...)
    {
        init(name,base_changes,acid_changes,-1.0f);
    }
}

CIsbtPtAllele::CIsbtPtAllele(std::string name, std::string base_changes, std::string acid_changes, float incidence) 
{
    init(name,base_changes,acid_changes,incidence);
}

CIsbtPtAllele::CIsbtPtAllele(const CIsbtPtAllele& orig) 
{
    m_name = orig.m_name;
    m_base_changes = orig.m_base_changes; // ISBT naming like 261delG or 85A>G, ...
    m_acid_changes = orig.m_acid_changes; // ISBT naming like Ile60Leu or Ser230Ile, ...
    m_incidence = orig.m_incidence;
}

CIsbtPtAllele::~CIsbtPtAllele() 
{
    
}
void CIsbtPtAllele::init(std::string name, std::string base_changes, std::string acid_changes, float incidence)
{
    m_name = name;
    m_incidence = incidence;
    
    vector<string> parsed = CMyTools::GetParsedLine(base_changes," ");
    for(auto x:parsed)
        m_base_changes.insert(x);
    parsed = CMyTools::GetParsedLine(acid_changes," ");
    for(auto x:parsed)
        m_acid_changes.insert(x);
     
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





