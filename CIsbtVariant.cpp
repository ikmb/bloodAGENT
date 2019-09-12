/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtVariant.cpp
 * Author: mwittig
 * 
 * Created on July 22, 2019, 8:08 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>
#include <limits>

#include <regex>

#include "CIsbtVariant.h"
#include "CBigWigReader.h"

using namespace std;

bool CIsbtVariant::verbose = false;

CIsbtVariant::CIsbtVariant() 
{
    m_isbt_name = "";
    m_chromosome = "";
    m_position = -1;
    m_strand = 'u';
    m_lrg_position = "";
    /// attention: m_lrg_reference will be revised in bool CIsbtVariant::parseIsbtVariant()
    m_lrg_reference = "";
    m_lrg_alternative = "";
    m_vcf_coordinate = -1;
    m_vcf_reference = "";
    m_vcf_alternative = "";
    m_coverage = std::numeric_limits<float>::quiet_NaN();
}

CIsbtVariant::CIsbtVariant(const string& lrg_anno, const string& refBase, const string& chrom, int pos, char strand, int vcfCoord, const string& vcfRef, const string& vcfAlt) 
{
    m_isbt_name = lrg_anno;
    m_chromosome = chrom;
    m_position = pos;
    m_strand = strand;
    m_lrg_position = "";
    /// attention: m_lrg_reference will be revised in bool CIsbtVariant::parseIsbtVariant()
    m_lrg_reference = refBase;
    m_lrg_alternative = "";
    m_vcf_coordinate = vcfCoord;
    m_vcf_reference = vcfRef;
    m_vcf_alternative = vcfAlt;
    m_coverage = std::numeric_limits<float>::quiet_NaN();
    parseIsbtVariant();
}

CIsbtVariant::CIsbtVariant(const CIsbtVariant& orig) 
{
    m_isbt_name = orig.m_isbt_name;
    m_lrg_position = orig.m_lrg_position;
    m_lrg_reference = orig.m_lrg_reference;
    m_lrg_alternative = orig.m_lrg_alternative;
    m_vcf_coordinate = orig.m_vcf_coordinate;
    m_vcf_reference = orig.m_vcf_reference;
    m_vcf_alternative = orig.m_vcf_alternative;
    
    m_chromosome = orig.m_chromosome;
    m_position = orig.m_position;
    m_strand = orig.m_strand;
    m_coverage = orig.m_coverage;
}

CIsbtVariant& CIsbtVariant::operator =(const CIsbtVariant& orig)
{
    m_isbt_name = orig.m_isbt_name;
    m_lrg_position = orig.m_lrg_position;
    m_lrg_reference = orig.m_lrg_reference;
    m_lrg_alternative = orig.m_lrg_alternative;
    m_vcf_coordinate = orig.m_vcf_coordinate;
    m_vcf_reference = orig.m_vcf_reference;
    m_vcf_alternative = orig.m_vcf_alternative;
    
    m_chromosome = orig.m_chromosome;
    m_position = orig.m_position;
    m_strand = orig.m_strand;
    m_coverage = orig.m_coverage;
    return *this;
}

bool   CIsbtVariant::operator <(const CIsbtVariant& orig)const
{
    int c = m_chromosome.compare(orig.m_chromosome);
    if(c < 0 )
        return true;
    if(c > 0)
        return false;
    if(m_position < orig.m_position)
        return true;
    if(m_position > orig.m_position)
        return false;
    c = m_isbt_name.compare(orig.m_isbt_name);
    if(c < 0 )
        return true;
    return false;
}

bool   CIsbtVariant::operator ==(const CIsbtVariant& orig)const
{
    return (
            m_chromosome.compare(orig.m_chromosome) == 0 && 
            m_position == orig.m_position && 
            m_strand == orig.m_strand && 
            m_isbt_name.compare(orig.m_isbt_name) == 0
            );
}

/*
 bool  CIsbtVariant::operator ==(const string& orig)const
 {
     return m_isbt_name.compare(orig) == 0;
 }
 */
 
std::ostream& operator<<(std::ostream& os, const CIsbtVariant& me)
{
    os << me.name();
    if(me.verbose)
        os << " (" << (me.m_coverage == me.m_coverage ? me.m_coverage : 0 ) << "x)";
    return os;
}

CIsbtVariant::~CIsbtVariant() 
{
    
}

bool CIsbtVariant::parseIsbtVariant()
{
    smatch m;
    regex e ("^[0-9_+-]{1,}");
    
    regex_search(m_isbt_name,m,e);
    for (auto x:m)
        m_lrg_position = x;
    
    string actBaseChange= m.suffix().str();
    size_t pos = actBaseChange.find('>');
    if(pos != string::npos)
    {
        m_lrg_reference = actBaseChange.substr(0,pos);
        m_lrg_alternative = actBaseChange.substr(pos+1);
    }
    else if(actBaseChange.substr(0,6).compare("delins")==0)
    {
        m_lrg_alternative = actBaseChange.substr(6);
        if(strand() == '-')
            m_lrg_reference = CMyTools::GetComplSequence(m_lrg_reference);
    }
    else if(actBaseChange.substr(0,3).compare("del")==0)
    {
        if(!actBaseChange.substr(3).empty())
            m_lrg_reference = actBaseChange.substr(3);
        else if(strand() == '-')
            m_lrg_reference = CMyTools::GetComplSequence(m_lrg_reference);
        m_lrg_alternative = "-";
    }
    else if(actBaseChange.substr(0,3).compare("ins")==0)
    {
        m_lrg_reference = "-";
        m_lrg_alternative = actBaseChange.substr(3);
    }
    else if(actBaseChange.substr(0,3).compare("dup")==0)
    {
        if(strand() == '-')
            m_lrg_reference = CMyTools::GetComplSequence(m_lrg_reference);
        m_lrg_alternative = m_lrg_reference + actBaseChange.substr(3);
    }
    else 
    {
        cerr << "Can not parse this variant annotation: " << m_lrg_position << " o" << actBaseChange << endl;
    }
    //cout << vRet.first << " o " << m_lrg_reference << " o " << m_lrg_alternative << endl;
    return true;
}

bool CIsbtVariant::addCoverage(const CBigWigReader& bigWig)
{
    m_coverage = bigWig.getMinCoverage(m_chromosome,m_position-1,m_position+1);
    return m_coverage == m_coverage;
}

bool CIsbtVariant::isCovered(double limit)const
{
    if( !(m_coverage == m_coverage) ) // is NaN
        return false;
    return m_coverage >= limit;
}