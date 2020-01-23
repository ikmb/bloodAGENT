/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CTranscript.h
 * Author: mwittig
 *
 * Created on August 7, 2019, 3:01 PM
 */

#ifndef CTRANSCRIPT_H
#define CTRANSCRIPT_H

#include <vector>

#include "CBigWigReader.h"

class CTranscript {
public:
    CTranscript(const std::string& chrom, const std::string& txStart, const std::string& txEnd, const std::string& exonStarts, const std::string& exonEnds);
    CTranscript(const CTranscript& orig);
    virtual ~CTranscript();
    
    int exonCount()const{return m_exonStarts.size();}
    std::vector<int> exonStarts()const{return m_exonStarts;}
    std::vector<int> exonEnds()const{return m_exonEnds;}
    int exonStart(int idx)const{return m_exonStarts[idx];}
    int exonEnd(int idx)const{return m_exonEnds[idx];}
    std::string getChrom()const{return m_chrom;}
    
private:
    std::string m_txStart;
    std::string m_txEnd;
    std::string m_chrom;
    std::vector<int> m_exonStarts;
    std::vector<int> m_exonEnds;
};

#endif /* CTRANSCRIPT_H */

