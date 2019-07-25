/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtPtAllele.h
 * Author: mwittig
 *
 * Created on July 25, 2019, 8:18 AM
 */

#ifndef CISBTPTALLELE_H
#define CISBTPTALLELE_H

class CIsbtPtAllele {
public:
    CIsbtPtAllele(std::string name, std::string base_changes, std::string acid_changes, float incidence);
    CIsbtPtAllele(const CIsbtPtAllele& orig);
    virtual ~CIsbtPtAllele();
private:
    
    void init(std::string name, std::string base_changes, std::string acid_changes, float incidence);
    
    std::string m_name;
    std::set<std::string>   m_base_changes; // ISBT naming like 261delG or 85A>G, ...
    std::set<std::string>   m_acid_changes; // ISBT naming like Ile60Leu or Ser230Ile, ...
    float                   m_incidence;
    

};

#endif /* CISBTPTALLELE_H */

