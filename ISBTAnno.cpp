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
#include <experimental/filesystem>
#include <htslib/sam.h>

#include "CBigWigReader.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"

CISBTAnno::CISBTAnno(const std::vector<std::string> filename, const std::string build) 
{
    m_build=build;
    for(auto act : filename)
    {
        if(!CMyTools::file_exists(act))
            throw(CMyException("File does not exist: ")+act);
        m_data_red = readAnnotation(act);
    }
    generateIndex();
}

CISBTAnno::CISBTAnno(const CISBTAnno& orig) 
{
    m_data_red = orig.m_data_red;
    m_entry_finder = orig.m_entry_finder;
    m_strand = orig.m_strand;
    m_parsed_isbt_variant = orig.m_parsed_isbt_variant;
    m_loci = orig.m_loci;
    m_coverage_failed = orig.m_coverage_failed;
    m_build = orig.m_build;
}

CISBTAnno::~CISBTAnno() {
}

bool CISBTAnno::generateIndex()
{
    m_entry_finder.clear();
    if(m_vanno.First())
    {
        do{
            ostringstream osr("");
            //osr << m_vanno["Chrom ("+m_build+")"] << '_' << m_vanno["Pos ("+m_build+")"];
            // I have to skip the high impact marker "!" when I build my index. This I do here
            // The fix starts here
            int subIdx = 0;
            string var = m_vanno["Transcript annotation short"];
            if(var.size() > 0)
                while(var[subIdx]=='!')
                    subIdx++;
            var = var.substr(subIdx);
            // the fix ends here and is used in the next line
            osr << m_vanno["Chrom ("+m_build+")"] << '_' << m_vanno["Coordinate in VCF "+m_build+""];
            m_isbt_variant_to_index[m_vanno["system/gene"]][var]=m_entry_finder.size();
            m_entry_finder.insert(pair<string,int>(osr.str(),m_entry_finder.size()));
            if(m_strand.find(m_vanno["system/gene"]) == m_strand.end())
                m_strand[m_vanno["system/gene"]]=m_vanno["strand ("+m_build+")"][0];
            m_loci.insert(m_vanno["system/gene"]);
            m_parsed_isbt_variant.push_back(CIsbtVariant(m_vanno["Transcript annotation short"], 
                    m_vanno["Reference base ("+m_build+")"], m_vanno["Chrom ("+m_build+")"], stoi(m_vanno["1-based end ("+m_build+")"]),m_vanno["strand ("+m_build+")"][0],
                    stoi(m_vanno["Coordinate in VCF "+m_build+""]),m_vanno["RefAllele in VCF "+m_build+""], m_vanno["AltAllele in VCF "+m_build+""],
                    m_vanno["is transcript_NC == "+m_build+"_NC"].compare("FALSE")==0,m_vanno["TYPE"]));
        }while(m_vanno.Next());
        return true;
    }
    return false;
}

bool CISBTAnno::readAnnotation(const std::string& filename)
{
    CParsedTextfile parsed(filename,"\t",-1,0,true,"#");
    if(m_vanno.size() == 0)
    {
        m_vanno = parsed;
    }
    else if(parsed.First())
    {
        do{
            m_vanno.add(parsed.lineObject());
        }while(parsed.Next());
    }
    return true;
  }

bool CISBTAnno::isVcfAlleleAnIsbtVariant(const std::string& allele, const std::string isbtVariant, const std::string system)
{
    CIsbtVariant var = m_parsed_isbt_variant[m_isbt_variant_to_index[system][isbtVariant]];
    return var.reference().compare(allele) == 0;
}

std::set<CISBTAnno::variation> CISBTAnno::getReferenceVariations(const std::string& system)
{
    set<CISBTAnno::variation> sRet;
    int i = 0;
    if(m_vanno.First())
    {
        do{
            if(m_vanno["system/gene"].compare(system)!=0)
                continue;
            if(m_vanno["is transcript_NC == "+m_build+"_NC"].compare("FALSE")==0)
            {
                std::string snp_short = m_vanno["Transcript annotation short"];
                while(!snp_short.empty() && snp_short[0]=='!')
                    snp_short.erase(snp_short.begin());
                sRet.insert(m_parsed_isbt_variant[m_isbt_variant_to_index[m_vanno["system/gene"]][snp_short]]);
            }
            i++;
        }while(m_vanno.Next());
    }
    return sRet;
}

map<string,vector<CISBTAnno::variation> > CISBTAnno::getReferenceVariations()
{
    map<string,vector<CISBTAnno::variation> > mRet;
    int i = 0;
    if(m_vanno.First())
    {
        do{
            //cout << m_vanno.line() << endl;
            if(m_vanno["is transcript_NC == "+m_build+"_NC"].compare("FALSE")==0)
            {
                std::string snp_short = m_vanno["Transcript annotation short"];
                while(!snp_short.empty() && snp_short[0]=='!')
                    snp_short.erase(snp_short.begin());
                mRet[m_vanno["system/gene"]].push_back(m_parsed_isbt_variant[m_isbt_variant_to_index[m_vanno["system/gene"]][snp_short]]);
            }
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
        const CParsedTextfile& vanno = m_vanno;
        string act_entry = vanno.cell(it->second,"Transcript annotation short");
        vRet.push_back(act_entry);
    }
    
    
    return vRet;
}

CISBTAnno::variation CISBTAnno::getCorrespondingIsbtVariation(CVcfSnp vcfsnp)const
{
    ostringstream osr("");
    osr << vcfsnp.chrom() << '_' << vcfsnp.pos();
    pair<multimap<string,int>::const_iterator, multimap<string,int>::const_iterator> mRange;
    //cout << osr.str() << endl;
    mRange = m_entry_finder.equal_range(osr.str());
    
    for (multimap<string,int>::const_iterator it=mRange.first; it!=mRange.second; ++it)
    {  
        if(static_cast<size_t>(it->second) >= m_parsed_isbt_variant.size())
            throw("something is out of bounds in string CISBTAnno::getCorrespondingIsbtVariations(CVcfSnp vcfsnp)const");
        variation varParsed = m_parsed_isbt_variant[it->second];
        int equal_counter = 0;
        int allele_counter = 0;
        for(auto s:vcfsnp.alleles())
        {
            allele_counter++;
            if( s.compare(varParsed.vcfReference()) == 0 || s.compare(varParsed.vcfAlternative()) == 0 )
                equal_counter++;
        }
        // TODO: This is not robust!!!
        // If a variation and and an in/del overlap in the annotation and we have the In/del in the sample
        // The indel will be ignored. Either change the ISBT helper table to exactly match the vcf output or 
        // parse the vcf output and adapt it to the helper table. Example
        // ABO del261G in helper table: -/C
        //                      in vcf: T/TC
        // another exception would be KEL with G for k+,A for K+ and C for Kmod at rs8176058
        if( (equal_counter == allele_counter || // alleles match perfect 
           distance(mRange.first, mRange.second) == 1 ) &&  // only one entry
           vcfsnp.refAllele().compare(varParsed.vcfReference()) == 0   )
            return m_parsed_isbt_variant[it->second];
    }
    return CISBTAnno::variation();
}
std::vector<CISBTAnno::variation>  CISBTAnno::getIsbtVariants(const string& system)const
{
    std::vector<CISBTAnno::variation> vRet;
    std::map<std::string,std::map<std::string,int>>::const_iterator i = m_isbt_variant_to_index.find(system);
    if(i != m_isbt_variant_to_index.end())
    {
        for(auto index : i->second)
            vRet.push_back(m_parsed_isbt_variant[index.second]);
    }
    return vRet;
}

size_t  CISBTAnno::getIsbtVariantCount(const string& system)const
{
    std::map<std::string,std::map<std::string,int>>::const_iterator i = m_isbt_variant_to_index.find(system);
    if(i != m_isbt_variant_to_index.end())
    {
        return i->second.size();
    }   
    return 0;
}

int CISBTAnno::getIsbtVariantIndex(const string& system,const string& isbt_var)const
{
    std::map<std::string,std::map<std::string,int>>::const_iterator i = m_isbt_variant_to_index.find(system);
    if(i != m_isbt_variant_to_index.end())
    {
        std::map<std::string,int>::const_iterator j = i->second.find(isbt_var);
        if(j != i->second.end())
            return j->second;
    }
    return -1;
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

bool CISBTAnno::addCoverage(const CBigWigReader& bigWig, int min_req_coverage)
{
    bool bRet = true;
    for(auto& x:m_parsed_isbt_variant)
    {
        if(!x.addCoverage(bigWig))
            bRet = false;
        if(x.getCoverage() < min_req_coverage)
        {
            string system = getSystemAt(x.chrom(),x.pos());
            m_coverage_failed[system].push_back(x);
        }
        
    }
    return bRet;
}

bool   CISBTAnno::hasUncoveredVariants(const string& system)const
{
    std::map<string,std::vector<variation>>::const_iterator iter = m_coverage_failed.find(system);
    if(iter == m_coverage_failed.end() || iter->second.empty())
        return false;
    return true;
}

std::vector<CISBTAnno::variation>   CISBTAnno::getAllVariations(const string& system)const
{
    vector<variation>  vRet;
    for(auto& x:m_parsed_isbt_variant)
    {
        if(getSystemAt(x.chrom(),x.vcfCoordinate()).compare(system) == 0)
            vRet.push_back(x);
    }
    return vRet;
}

std::vector<CISBTAnno::variation>   CISBTAnno::getCoverageFailedVariants(const string& system)const
{
    std::map<string,std::vector<variation>>::const_iterator iter = m_coverage_failed.find(system);
    if(iter == m_coverage_failed.end())
        return std::vector<variation>();
    return iter->second;
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

std::ostream& operator<<(std::ostream& os, const CISBTAnno& me)
{
    bool bFirst = true;
    for(auto locus:me.m_loci)
    {
        std::map<std::string,std::map<std::string,int>>::const_iterator i = me.m_isbt_variant_to_index.find(locus);
        if(i != me.m_isbt_variant_to_index.end())
        {
            for(auto variant:i->second)
            {
                if(!bFirst)
                    os << endl;
                else
                    bFirst=false;
                os << locus << ' ' << me.m_parsed_isbt_variant[variant.second];

            }
        }
    }
    return os;
}

