/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CPossibleGenotype.h
 * Author: mwittig
 *
 * Created on July 24, 2019, 9:05 AM
 */

#ifndef CISBTGTALLELE_H
#define CISBTGTALLELE_H

#include "CIsbtVariant.h"


class CIsbtGtAllele {
public:
    CIsbtGtAllele();
    CIsbtGtAllele(const CIsbtGtAllele& orig);
    virtual ~CIsbtGtAllele();
    
    CIsbtGtAllele& operator =(const CIsbtGtAllele& orig);
    bool          operator <(const CIsbtGtAllele& orig)const;
    bool          operator >(const CIsbtGtAllele& orig)const{return !(*this < orig || *this == orig);}
    bool          operator ==(const CIsbtGtAllele& orig)const;
    bool          operator !=(const CIsbtGtAllele& orig)const{return !(*this == orig);}
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtGtAllele& me);
    
    bool empty()const{return m_gt.empty();}
    long unsigned int size()const{return m_gt.size();}
    
    
    /// returns true if the CIsbtVariant is new
    bool add(const CIsbtVariant& var);
    
private:
    
    std::set<CIsbtVariant>  m_gt;

};

#endif /* CISBTGTALLELE_H */
