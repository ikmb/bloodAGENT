/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2Pt.h
 * Author: mwittig
 *
 * Created on July 25, 2019, 7:00 AM
 */

#ifndef CISBTGT2PT_H
#define CISBTGT2PT_H

#include "CIsbtPtAllele.h"


class CIsbtGt2Pt {
public:
    CIsbtGt2Pt(const std::string& filename);
    CIsbtGt2Pt(const CIsbtGt2Pt& orig);
    virtual ~CIsbtGt2Pt();
private:
    
    void init(const std::string& filename);
    
    std::map<std::string,vector<CIsbtPtAllele>> m_allele_vector;

};

#endif /* CISBTGT2PT_H */

