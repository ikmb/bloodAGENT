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

#include "mytools.h"
#include "CTranscript.h"
#include "CTranscriptAnno.h"

using namespace std;


CTranscriptAnno::CTranscriptAnno(const string& filename) 
{
    CParsedTextfile anno(filename,"\t","refGene.name2",0, true,  "#");
    if(anno.First())
        do{
            CTranscript act_entry(anno["refGene.chrom"],anno["refGene.txStart"],anno["refGene.txEnd"],anno["refGene.exonStarts"],anno["refGene.exonEnds"]);
            m_transcripts.insert(pair<std::string,CTranscript>(anno["refGene.name2"],act_entry));
        }while(anno.Next());
}

CTranscriptAnno::CTranscriptAnno(const CTranscriptAnno& orig)
{
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
    return CTranscript("","","","","");
    
}

double CTranscriptAnno::getExonicCoverage(const string& target, const CBigWigReader& bw)
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


