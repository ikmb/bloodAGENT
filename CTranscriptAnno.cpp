/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CTranscriptAnno.cpp
 * Author: mwittig
 * 
 * Created on August 7, 2019, 3:00 PM
 */

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <htslib/sam.h>

#include "meinetools.h"
#include "CBigWigReader.h"
#include "CTranscript.h"
#include "CTranscriptAnno.h"

using namespace std;


CTranscriptAnno::CTranscriptAnno(const string& filename) 
{
    if(!CMyTools::file_exists(filename))
        throw(CMyException("File does not exist: ")+filename);
    CParsedTextfile anno(filename,"\t","refGene.name2",0, true,  "#");
    if(anno.First())
        do{
            CTranscript act_entry(anno["refGene.chrom"],anno["refGene.strand"],anno["refGene.txStart"],anno["refGene.txEnd"],anno["refGene.exonStarts"],anno["refGene.exonEnds"]);
            m_transcripts.insert(pair<std::string,CTranscript>(anno["System"],act_entry));
        }while(anno.Next());
}

CTranscriptAnno::CTranscriptAnno(const CTranscriptAnno& orig)
{
    m_transcripts = orig.m_transcripts;
}

CTranscriptAnno& CTranscriptAnno::operator= (const CTranscriptAnno& orig)
{
    m_transcripts = orig.m_transcripts;
    return *this;
}

CTranscriptAnno::~CTranscriptAnno() {
}

std::set<std::string> CTranscriptAnno::loci()const
{
    std::set<std::string> sRert;
    for(const auto& x:m_transcripts)
        sRert.insert(x.first);
    return sRert;
}

CTranscript CTranscriptAnno::getTranscript(const std::string& name)
{
    std::map<std::string,CTranscript>::const_iterator i =  m_transcripts.find(name);
    if(i!=m_transcripts.end())
        return i->second;
    return CTranscript("","","","","","");
    
}

double CTranscriptAnno::getExonicCoverage(const string& target, const CBigWigReader& bw)const
{
    double dRet = 0.0;
    size_t base_count = 0;
    std::map<std::string,CTranscript>::const_iterator i =  m_transcripts.find(target);
    if(i!=m_transcripts.end())
    {
         
        const CTranscript& trans = i->second;
        for(int i = 0; i < trans.exonCount(); i++)
        {
            int start = trans.exonStart(i);
            int end = trans.exonEnd(i);
            string chrom = trans.getChrom();
            
            double act_cov =  bw.getSumCoverage(chrom,start,end);
            dRet += (act_cov == act_cov ? act_cov : 0); // is not nan?
            base_count+=(end-start)+1;
        }
        dRet/=static_cast<double>(base_count);
    }
    return dRet;
}

vector<double> CTranscriptAnno::getCoverages(const string& target, const CBigWigReader& bw)const
{
    vector<double> dRet;
    dRet.push_back(0.0);
    
    std::map<std::string,CTranscript>::const_iterator i =  m_transcripts.find(target);
    size_t base_count = 0;
    if(i!=m_transcripts.end())
    {
        const CTranscript& trans = i->second;
        for(int i = 0; i < trans.exonCount(); i++)
        {
            int start = trans.exonStart(i);
            int end = trans.exonEnd(i);
            string chrom = trans.getChrom();
            
            double act_cov =  bw.getSumCoverage(chrom,start,end);
            dRet[0]+=(act_cov == act_cov ? act_cov : 0); // is not nan?
            act_cov = (act_cov == act_cov ? act_cov/static_cast<double>((end-start)+1) : 0); // is not nan?
            dRet.push_back(act_cov);
            base_count+=(end-start)+1;
        }
        dRet[0]/=static_cast<double>(base_count);
        if(i->second.getStrand().compare("-")==0 && i->second.exonStarts().size()>1)
            std::reverse(dRet.begin()+1,dRet.end());
    }
    
    return dRet;
}



