/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   CMotifFinder.cpp
 * Author: mwittig
 * 
 * Created on January 5, 2022, 10:34 AM
 */
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "vcf.h"
#include "CVcfSnp.h"
#include "gzstream.h"
#include "meinetools.h"
#include "CFastqReader.h"
#include "CMotifFinder.h"

using namespace std;

CMotifFinder::CMotifFinder(const std::string& config, const std::string& filenames,int verbose) 
{
    CParsedTextfile motif_config(config,"\t",-1,0, true,"#");
    map<string,int> motifs;
    if(motif_config.First())
        do{
            vector<string> entries = CMyTools::Tokenize(motif_config["seq-A"],",");
            for(auto a : entries)          
            {
                motifs[a]=0;
                motifs[CMyTools::GetComplSequence(a)]=0;
            }
            entries = CMyTools::Tokenize(motif_config["seq-B"],",");
            for(auto a : entries)          
            {
                motifs[a]=0;
                motifs[CMyTools::GetComplSequence(a)]=0;
            }
        }while(motif_config.Next());
    vector<string> entries = CMyTools::Tokenize(filenames,",");
    for(auto a : entries)  
    {
        if(verbose >= 2)
        {
            cerr << "motif search in " << a << endl;
            findMotifs(a,motifs);
        }
    }
    storeMotifSnps(motif_config,motifs);
}

CMotifFinder::CMotifFinder(const CMotifFinder& orig)
{
    m_motifs_snps = orig.m_motifs_snps;
}

CMotifFinder::~CMotifFinder() 
{
    
}

std::map<string,int>  CMotifFinder::findMotifs(const std::string& filename,std::vector<string>& motifs)
{
    std::map<string,int> mRet;
    
    CFastqReader fastq(filename);
    string strLine;
    while(fastq.getline(strLine))
    {
        for(std::vector<string>::const_iterator i = motifs.begin(); i != motifs.end();i++)
        {
            if(strLine.find(*i) != string::npos)
                mRet[*i]++;
        }
        fastq.getline(strLine);
        fastq.getline(strLine);
        fastq.getline(strLine);
    }
    return mRet;
}

void  CMotifFinder::findMotifs(const std::string& filename,std::map<string,int>& motifs)
{
    
    CFastqReader fastq(filename);
    string strLine;
    while(fastq.getline(strLine) && fastq.getline(strLine))
    {
        for(std::map<string,int>::const_iterator i = motifs.begin(); i != motifs.end();i++)
        {
            if(strLine.find(i->first) != string::npos || strLine.find(CMyTools::GetComplSequence(i->first)) != string::npos)
                motifs[i->first]++;
        }
        ;
        fastq.getline(strLine);
        fastq.getline(strLine);
    }
}

void CMotifFinder::storeMotifSnps(CParsedTextfile& config, map<string,int> motifs)
{
    // map<std::string,CVcfSnp>    m_motifs_snps;
    if(config.First())
        do{
            
            /*
            m_chrom = 
            m_pos =       
            m_depth = 
            m_alleles.push_back(rec->d.allele[bcfA]);
            m_coverage.push_back(ad[i]);
            m_qualities.push_back(gq[i]);
            m_haplotype_qualities.push_back(hq[i]);
            m_phasing_id = ps[0];
            m_mapping_quality = mq[0];
            m_ref_allele = rec->d.allele[0];  
             */
            
            std::string                 chrom = config["Chr (hg19)"];
            long                        pos = atol(config["Coordinate in VCF hg19"].c_str());
            std::vector<std::string>    alleles;
            std::vector<int>            vcoverage;
            long                        depth=0;
            bool                        verbose = false;
            int                         phasing_id=-1;
            std::string                 ref_allele = config["RefAllele in VCF hg19"];
            int                         mapping_quality=50;

            
            int coverage = 0;
            CVcfSnp actSNP;
            vector<string> entries = CMyTools::Tokenize(config["seq-A"],",");
            for(auto a : entries)          
                coverage+=motifs[a];
            depth+=coverage;
            if(coverage != 0)
            {
                vcoverage.push_back(coverage);
                alleles.push_back(config["RefAllele in VCF hg19"]);
            }
            coverage = 0;
            entries = CMyTools::Tokenize(config["seq-B"],",");
            for(auto a : entries)          
                coverage+=motifs[a];
            depth+=coverage;
            if(coverage != 0)
            {
                vcoverage.push_back(coverage);
                alleles.push_back(config["AltAllele in VCF hg19"]);
            }
            std::vector<int>            qualities(vcoverage.size(),50);
            std::vector<int>            haplotype_qualities(vcoverage.size(),0);
            
            m_motifs_snps[config["system"]][config["short"]]=CVcfSnp(chrom, pos,alleles,vcoverage,qualities,haplotype_qualities,mapping_quality,
                                                   depth,verbose,phasing_id,ref_allele);
        }while(config.Next());
}



map<std::string,CVcfSnp> CMotifFinder::getSystemsMotifSnps(const std::string& system)const
{
    map<std::string,map<std::string,CVcfSnp> >::const_iterator i = m_motifs_snps.find(system);
    if(i != m_motifs_snps.end())
        return i->second;
    return map<std::string,CVcfSnp>();
}

std::set<std::string> CMotifFinder::getSystems()const
{
    set<std::string> sRet;
    for(map<std::string,map<std::string,CVcfSnp> >::const_iterator i = m_motifs_snps.begin(); i != m_motifs_snps.end(); i++)
        sRet.insert(i->first);
    return sRet;
}


std::ostream& operator<<(std::ostream& os, const CMotifFinder& me)
{
    bool line_feed = false;
    for(const auto& s : me.getSystems())
        for(const auto& e : me.getSystemsMotifSnps(s))
        {
            os << (line_feed ? "\n" : "") << s << '\t' << e.second;
            line_feed=true;
        }
    return os;
}


