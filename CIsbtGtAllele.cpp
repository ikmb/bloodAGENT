/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CPossibleGenotype.cpp
 * Author: mwittig
 * 
 * Created on July 24, 2019, 9:05 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>


#include "CIsbtGtAllele.h"

using namespace std;

CIsbtGtAllele::CIsbtGtAllele() 
{
}

CIsbtGtAllele::CIsbtGtAllele(const CIsbtGtAllele& orig) 
{
    m_gt = orig.m_gt;
}

CIsbtGtAllele& CIsbtGtAllele::operator =(const CIsbtGtAllele& orig)
{
    m_gt = orig.m_gt;
}

CIsbtGtAllele::~CIsbtGtAllele() 
{
    
}

bool CIsbtGtAllele::operator <(const CIsbtGtAllele& orig)const
{
    std::set<CIsbtVariant>::const_iterator i = m_gt.begin();
    std::set<CIsbtVariant>::const_iterator j = orig.m_gt.begin();
    for(;i != m_gt.end() && j!= orig.m_gt.end();i++,j++)
    {
        if(*i < *j)
            return true;
        if(*i > *j)
            return false;
    }
    if(m_gt.size() < orig.m_gt.size())
        return true;
    return false;
}

bool CIsbtGtAllele::operator ==(const CIsbtGtAllele& orig)const
{
    if(m_gt.size() != orig.m_gt.size())
        return false;
    std::set<CIsbtVariant>::const_iterator i = m_gt.begin();
    std::set<CIsbtVariant>::const_iterator j = orig.m_gt.begin();
    for(;i != m_gt.end() && j!= orig.m_gt.end();i++,j++)
    {
        if(*i != *j)
            return false;
    }
    return true;
}

bool CIsbtGtAllele::add(const CIsbtVariant& var)
{
    return m_gt.insert(var).second;
}


std::ostream& operator<<(std::ostream& os, const CIsbtGtAllele& me)
{
    long unsigned int i = 0;
    for(auto x:me.m_gt)
    {
        os << x << ( ++i == me.m_gt.size() ? "" : " ");
    }
    return os;
}
