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
#include "mytools.h"
#include "CFastqReader.h"
#include "CMotifFinder.h"

using namespace std;

CMotifFinder::CMotifFinder(const std::string& config, std::string& filenames) 
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
    vector<string> entries = CMyTools::Tokenize(filenames);
    for(auto a : entries)  
        findMotifs(a,motifs);
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
            if(strLine.find(*i) != string::npos || strLine.find(CMyTools::GetComplSequence(*i)) != string::npos)
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
    while(fastq.getline(strLine))
    {
        for(std::map<string,int>::const_iterator i = motifs.begin(); i != motifs.end();i++)
        {
            if(strLine.find(i->first) != string::npos || strLine.find(CMyTools::GetComplSequence(i->first)) != string::npos)
                motifs[i->first]++;
        }
        fastq.getline(strLine);
        fastq.getline(strLine);
        fastq.getline(strLine);
    }
}