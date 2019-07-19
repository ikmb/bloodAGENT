/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ISBTAnno.cpp
 * Author: mwittig
 * 
 * Created on July 18, 2019, 2:05 PM
 */

#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "ISBTAnno.h"

CISBTAnno::CISBTAnno(const std::string& filename) 
{
    m_data_red = readAnnotation(filename);
}

CISBTAnno::CISBTAnno(const CISBTAnno& orig) 
{
    m_data_red = orig.m_data_red;
    m_entry_finder = orig.m_entry_finder;
    m_strand = orig.m_strand;
    m_parsed_isbt_variant = orig.m_parsed_isbt_variant;
}

CISBTAnno::~CISBTAnno() {
}

bool CISBTAnno::readAnnotation(const std::string& filename)
{
    m_vanno = CParsedTextfile(filename,"\t",-1,0,true,"#");
    m_entry_finder.clear();
    if(m_vanno.First())
    {
        do{
            ostringstream osr("");
            osr << m_vanno["Chrom (hg19)"] << '_' << m_vanno["1-based end (hg19)"];
            m_entry_finder.insert(pair<string,int>(osr.str(),m_entry_finder.size()));
            if(m_strand.find(m_vanno["system/gene"]) == m_strand.end())
                m_strand[m_vanno["system/gene"]]=m_vanno["strand (hg19)"][0];
            m_parsed_isbt_variant.push_back(parseIsbtVariant(m_vanno["Transcript annotation short"]));
        }while(m_vanno.Next());
        return true;
    }
    return false;
}

bool CISBTAnno::isVcfAlleleAnIsbtVariant(const std::string& allele, const std::string isbtVariant, const std::string system)
{
    variation var = parseIsbtVariant(isbtVariant);
    if(strand(system) == '-')
        return var.second.second.compare(CMyTools::GetComplSequence(allele)) == 0;
    return var.second.second.compare(allele) == 0;
}

CISBTAnno::variation CISBTAnno::parseIsbtVariant(string var)
{
    variation vRet("",pair<string,string>("",""));
    
    smatch m;
    regex e ("^[0-9_+-]{1,}");
    
    regex_search(var,m,e);
    for (auto x:m)
        vRet.first = x;
    
    string actBaseChange= m.suffix().str();
    size_t pos = actBaseChange.find('>');
    if(pos != string::npos)
    {
        vRet.second.first = actBaseChange.substr(0,pos);
        vRet.second.second = actBaseChange.substr(pos+1);
    }
    else if(actBaseChange.substr(0,6).compare("delins")==0)
    {
        vRet.second.first = m_vanno["Reference base (hg19)"];
        vRet.second.second = actBaseChange.substr(6);
        if(strand(m_vanno["system/gene"]) == '-')
            vRet.second.first = CMyTools::GetComplSequence(vRet.second.first);
    }
    else if(actBaseChange.substr(0,3).compare("del")==0)
    {
        vRet.second.first = actBaseChange.substr(3);
        if(vRet.second.first.size() == 0) // sequence not part of annotation, so get it from table
        {
            vRet.second.first = m_vanno["Reference base (hg19)"];
            if(strand(m_vanno["system/gene"]) == '-')
                vRet.second.first = CMyTools::GetComplSequence(vRet.second.first);
        }
        vRet.second.second = "-";
    }
    else if(actBaseChange.substr(0,3).compare("ins")==0)
    {
        vRet.second.first = "-";
        vRet.second.second = actBaseChange.substr(3);
    }
    else if(actBaseChange.substr(0,3).compare("dup")==0)
    {
        vRet.second.first = m_vanno["Reference base (hg19)"];
        if(strand(m_vanno["system/gene"]) == '-')
            vRet.second.first = CMyTools::GetComplSequence(vRet.second.first);
        vRet.second.second = vRet.second.first + actBaseChange.substr(3);
    }
    else 
    {
        cerr << "Can not parse this variant annotation: " << vRet.first << " o" << actBaseChange << endl;
    }
    //cout << vRet.first << " o " << vRet.second.first << " o " << vRet.second.second << endl;
    return vRet;
}

map<string,vector<string> > CISBTAnno::getReferenceVariations()
{
    map<string,vector<string> > mRet;
    if(m_vanno.First())
    {
        do{
            if(m_vanno["is transcript_NC == hg19_NC"].compare("FALSE")==0)
                mRet[m_vanno["system/gene"]].push_back(m_vanno["Transcript annotation short"]);
        }while(m_vanno.Next());
    }
    return mRet;
}

std::vector<std::string> CISBTAnno::getVariationsAt(std::string chrom, int pos)const
{
    vector<string> vRet;
    ostringstream osr("");
    osr << chrom << '_' << pos;
    pair<multimap<string,int>::const_iterator, multimap<string,int>::const_iterator> mRange;
    mRange = m_entry_finder.equal_range(osr.str());
    
    for (multimap<string,int>::const_iterator it=mRange.first; it!=mRange.second; ++it)
      vRet.push_back(m_vanno.cell(it->second,"Transcript annotation short"));
    
    return vRet;
}

string CISBTAnno::getCorrespondingIsbtVariation(CVcfSnp vcfsnp)const
{
    ostringstream osr("");
    osr << vcfsnp.chrom() << '_' << vcfsnp.pos();
    pair<multimap<string,int>::const_iterator, multimap<string,int>::const_iterator> mRange;
    mRange = m_entry_finder.equal_range(osr.str());
    
    for (multimap<string,int>::const_iterator it=mRange.first; it!=mRange.second; ++it)
    {  
        if(it->second >= m_parsed_isbt_variant.size())
            throw("something is out of bounds in string CISBTAnno::getCorrespondingIsbtVariations(CVcfSnp vcfsnp)const");
        variation varParsed = m_parsed_isbt_variant[it->second];
        if(strand(m_vanno.cell(it->second,"system/gene")) == '-')
        {
            varParsed.second.first  = CMyTools::GetComplSequence(varParsed.second.first);
            varParsed.second.second = CMyTools::GetComplSequence(varParsed.second.second);
        }
        int equal_counter = 0;
        int allele_counter = 0;
        for(auto s:vcfsnp.alleles())
        {
            allele_counter++;
            if(s.compare(varParsed.second.first) == 0 || s.compare(varParsed.second.second) == 0)
                equal_counter++;
        }
        // TODO: This is not robust!!!
        // If a variation and and an in/del overlap in the annotation and we have the In/del in the sample
        // The indel will be ignored. Either change the ISBT helper table to exactly match the vcf output or 
        // parse the vcf output and adapt it to the helper table. Example
        // ABO del261G in helper table: -/C
        //                      in vcf: T/TC
        if(equal_counter == allele_counter || // alleles match perfect 
           abs(mRange.second->second - mRange.first->second) == 1  ) // only one entry
            return m_vanno.cell(it->second,"Transcript annotation short");
    }
    return "";
}



string CISBTAnno::getSystemAt(std::string chrom, int pos)const
{
    ostringstream osr("");
    osr << chrom << '_' << pos;
    pair<multimap<string,int>::const_iterator, multimap<string,int>::const_iterator> mRange;
    mRange = m_entry_finder.equal_range(osr.str());
    
    if(mRange.first != mRange.second)
        return m_vanno.cell(mRange.first->second,"system/gene");
    
    
    return "";
}
