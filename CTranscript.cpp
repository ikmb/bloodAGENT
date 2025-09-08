/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CTranscript.cpp
 * Author: mwittig
 * 
 * Created on August 7, 2019, 3:01 PM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <htslib/sam.h>

#include "meinetools.h"
#include "CBigWigReader.h"
#include "CTranscript.h"

using namespace std;

CTranscript::CTranscript(const std::string& chrom, const std::string& strand, const std::string& txStart, const std::string& txEnd, const std::string& exonStarts, const std::string& exonEnds) 
{
    m_txStart = txStart;
    m_txEnd = txEnd;
    m_chrom = chrom;
    m_strand = strand;
    vector<string> parsed = CMyTools::GetParsedLine(exonStarts,",");
    for(auto x:parsed)
    {
        if(x.length() != 0)
            m_exonStarts.push_back(stoi(x));
    }
    parsed = CMyTools::GetParsedLine(exonEnds,",");
    for(auto x:parsed)
    {
        if(x.length() != 0)
            m_exonEnds.push_back(stoi(x));
    }
}

CTranscript::CTranscript(const CTranscript& orig) 
{
    m_txStart = orig.m_txStart;
    m_txEnd = orig.m_txEnd;
    m_chrom = orig.m_chrom;
    m_strand = orig.m_strand;
    m_exonStarts = orig.m_exonStarts;
    m_exonEnds = orig.m_exonEnds;
}

CTranscript::~CTranscript() 
{
    
}

