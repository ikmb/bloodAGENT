/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt.cpp
 * Author: mwittig
 * 
 * Created on July 24, 2019, 9:15 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>


#include "CIsbtGtAllele.h"
#include "CIsbtGt.h"

using namespace std;

CIsbtGt::CIsbtGt() 
{
    
}

CIsbtGt::CIsbtGt(const CIsbtGt& orig) 
{
    m_gt = orig.m_gt;
}

CIsbtGt& CIsbtGt::operator =(const CIsbtGt& orig)
{
    m_gt = orig.m_gt;
    return *this;
}

CIsbtGt::~CIsbtGt() {
}

bool CIsbtGt::operator <(const CIsbtGt& orig)const
{
    std::set<CIsbtGtAllele>::const_iterator i = m_gt.begin();
    std::set<CIsbtGtAllele>::const_iterator j = orig.m_gt.begin();
    for(;i != m_gt.end() && j!= orig.m_gt.end();i++,j++)
    {
        if(*i < *j)
            return true;
    }
    if(m_gt.size() < orig.m_gt.size())
        return true;
    return false;
}

bool CIsbtGt::operator ==(const CIsbtGt& orig)const
{
    if(m_gt.size() != orig.m_gt.size())
        return false;
    std::set<CIsbtGtAllele>::const_iterator i = m_gt.begin();
    std::set<CIsbtGtAllele>::const_iterator j = orig.m_gt.begin();
    for(;i != m_gt.end() && j!= orig.m_gt.end();i++,j++)
    {
        if(*i != *j)
            return false;
    }
    return true;
}

bool CIsbtGt::add(const CIsbtGtAllele& var)
{
    return m_gt.insert(var).second;
}

std::ostream& operator<<(std::ostream& os, const CIsbtGt& me)
{
    long unsigned int i = 0;
    for(auto x:me.m_gt)
    {
        os << x << ( ++i == me.m_gt.size() ? "" : " & ");
        if(me.isHomozygous())
            os << " homozygous";
    }
    return os;
}