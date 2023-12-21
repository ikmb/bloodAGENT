/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   CScoreHaplotype.h
 * Author: mwittig
 *
 * Created on December 1, 2023, 3:28 PM
 */


/**************************
 
 NOT IN USE!!!
 
 **************************/

#ifndef CSCOREHAPLOTYPE_H
#define CSCOREHAPLOTYPE_H

class CScoreHaplotype {
public:
    CScoreHaplotype(const std::set<CIsbtVariant>& haplotype1, const std::set<CIsbtVariant>& haplotype2)
        : m_haplotype1(haplotype1), m_haplotype2(haplotype2) {}
    CScoreHaplotype(const CScoreHaplotype& other)
        : m_haplotype1(other.m_haplotype1), m_haplotype2(other.m_haplotype2) {}

    virtual ~CScoreHaplotype();
    
    int performAlignment();
private:
    const std::set<CIsbtVariant>& m_haplotype1; // ISBT
    const std::set<CIsbtVariant>& m_haplotype2; // SAMPLE
};

#endif /* CSCOREHAPLOTYPE_H */

