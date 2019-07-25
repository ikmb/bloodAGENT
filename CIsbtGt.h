/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt.h
 * Author: mwittig
 *
 * Created on July 24, 2019, 9:15 AM
 */

#ifndef CISBTGT_H
#define CISBTGT_H

#include "CIsbtGtAllele.h"


class CIsbtGt {
public:
    CIsbtGt();
    CIsbtGt(const CIsbtGt& orig);
    virtual ~CIsbtGt();
    
    CIsbtGt&        operator =(const CIsbtGt& orig);
    bool          operator <(const CIsbtGt& orig)const;
    bool          operator >(const CIsbtGt& orig)const{return !(*this < orig || *this == orig);}
    bool          operator ==(const CIsbtGt& orig)const;
    bool          operator !=(const CIsbtGt& orig)const{return !(*this == orig);}
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtGt& me);
    
    bool empty()const{return m_gt.empty();}
    long unsigned int size()const{return m_gt.size();}
    
    
    /// returns true if the CIsbtVariant is new
    bool add(const CIsbtGtAllele& var);
    
    std::set<CIsbtGtAllele> getAlleles()const{return m_gt;}
    
private:
    std::set<CIsbtGtAllele> m_gt;
};

#endif /* CISBTGT_H */

