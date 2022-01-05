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

#include "gzstream.h"
#include "mytools.h"
#include "CFastqReader.h"
#include "CMotifFinder.h"

using namespace std;

CMotifFinder::CMotifFinder() {
}



CMotifFinder::~CMotifFinder() {
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