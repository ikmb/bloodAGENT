/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CTranscriptAnno.h
 * Author: mwittig
 *
 * Created on August 7, 2019, 3:00 PM
 */

#ifndef CTRANSCRIPTANNO_H
#define CTRANSCRIPTANNO_H

#include <map>

#include "CBigWigReader.h"
#include "CTranscript.h"


class CTranscriptAnno {
public:
    CTranscriptAnno(){}
    CTranscriptAnno(const std::string& filename);
    CTranscriptAnno(const CTranscriptAnno& orig);
    CTranscriptAnno& operator= (const CTranscriptAnno& orig);
    virtual ~CTranscriptAnno();
    
    std::set<std::string> loci()const;
    CTranscript getTranscript(const std::string& name);
    bool hasKey(const std::string& v)const{return m_transcripts.find(v)!=m_transcripts.end();}
    
    double getExonicCoverage(const string& target, const CBigWigReader& bw)const;
    
    vector<double> getCoverages(const string& target, const CBigWigReader& bw)const;
    
private:
    
    std::map<std::string,CTranscript>  m_transcripts;
};

#endif /* CTRANSCRIPTANNO_H */

