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

#include "CIsbtVariant.h"
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
    m_loci = orig.m_loci;
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
            m_isbt_variant_to_index[m_vanno["system/gene"]][m_vanno["Transcript annotation short"]]=m_entry_finder.size();
            m_entry_finder.insert(pair<string,int>(osr.str(),m_entry_finder.size()));
            if(m_strand.find(m_vanno["system/gene"]) == m_strand.end())
                m_strand[m_vanno["system/gene"]]=m_vanno["strand (hg19)"][0];
            m_loci.insert(m_vanno["system/gene"]);
            m_parsed_isbt_variant.push_back(CIsbtVariant(m_vanno["Transcript annotation short"], 
                    m_vanno["Reference base (hg19)"], m_vanno["Chrom (hg19)"], stoi(m_vanno["1-based end (hg19)"]),m_vanno["strand (hg19)"][0]));
        }while(m_vanno.Next());
        return true;
    }
    return false;
}

bool CISBTAnno::isVcfAlleleAnIsbtVariant(const std::string& allele, const std::string isbtVariant, const std::string system)
{
    CIsbtVariant var = m_parsed_isbt_variant[m_isbt_variant_to_index[system][isbtVariant]];
    return var.reference().compare(allele) == 0;
}

map<string,vector<CISBTAnno::variation> > CISBTAnno::getReferenceVariations()
{
    map<string,vector<CISBTAnno::variation> > mRet;
    int i = 0;
    if(m_vanno.First())
    {
        do{
            if(m_vanno["is transcript_NC == hg19_NC"].compare("FALSE")==0)
                mRet[m_vanno["system/gene"]].push_back(m_parsed_isbt_variant[m_isbt_variant_to_index[m_vanno["system/gene"]][m_vanno["Transcript annotation short"]]]);
            i++;
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
    {
        string act_entry = m_vanno.cell(it->second,"Transcript annotation short");
        vRet.push_back(act_entry);
    }
    
    
    return vRet;
}

CISBTAnno::variation CISBTAnno::getCorrespondingIsbtVariation(CVcfSnp vcfsnp)const
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
        int equal_counter = 0;
        int allele_counter = 0;
        for(auto s:vcfsnp.alleles())
        {
            allele_counter++;
            if(s.compare(varParsed.reference()) == 0 || s.compare(varParsed.alternative()) == 0)
                equal_counter++;
        }
        // TODO: This is not robust!!!
        // If a variation and and an in/del overlap in the annotation and we have the In/del in the sample
        // The indel will be ignored. Either change the ISBT helper table to exactly match the vcf output or 
        // parse the vcf output and adapt it to the helper table. Example
        // ABO del261G in helper table: -/C
        //                      in vcf: T/TC
        // another exception would be KEL with G for k+,A for K+ and C for Kmod at rs8176058
        if(equal_counter == allele_counter || // alleles match perfect 
           abs(mRange.second->second - mRange.first->second) == 1  ) // only one entry
            return m_parsed_isbt_variant[it->second];
    }
    return CISBTAnno::variation();
}

CISBTAnno::variation CISBTAnno::getIsbtVariant(const string& system,const string& isbt_var)const
{
    std::map<std::string,std::map<std::string,int>>::const_iterator i = m_isbt_variant_to_index.find(system);
    if(i != m_isbt_variant_to_index.end())
    {
        std::map<std::string,int>::const_iterator j = i->second.find(isbt_var);
        if(j != i->second.end())
            return m_parsed_isbt_variant[j->second];
    }
    return  variation();
}




string CISBTAnno::getSystemAt(std::string chrom, int pos)const
{
    ostringstream osr("");
    osr << chrom << '_' << pos;
    pair<multimap<string,int>::const_iterator, multimap<string,int>::const_iterator> mRange;
    mRange = m_entry_finder.equal_range(osr.str());
    
    if(mRange.first != mRange.second)
    {
        string act_entry = m_vanno.cell(mRange.first->second,"system/gene");
        return act_entry;
    }
    
    
    return "";
}
