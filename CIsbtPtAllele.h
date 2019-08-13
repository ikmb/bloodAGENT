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

#include <string>


class CIsbtPtAllele {
public:
    CIsbtPtAllele(std::string name, std::string phenotype, std::string base_changes, std::string acid_changes, std::string incidence);
    CIsbtPtAllele(std::string name, std::string phenotype, std::string base_changes, std::string acid_changes, float incidence);
    CIsbtPtAllele(const CIsbtPtAllele& orig);
    virtual ~CIsbtPtAllele();
    
    friend std::ostream& operator<<(std::ostream& os, const CIsbtPtAllele& me);
    
    std::set<std::string>  baseChanges()const{return m_base_changes;}
    std::set<std::string>  acidChanges()const{return m_acid_changes;}
    std::string name()const{return m_name;}
    std::string phenotype()const{if(m_phenotype_name.compare("#N/A") == 0)return m_name;return m_phenotype_name;}
    
    
    bool containsBaseChange(const std::string& isbt_base_change)const;
    
private:
    
    void init(std::string name, std::string phenotype, std::string base_changes, std::string acid_changes, float incidence);
    
    std::string m_name;
    std::string m_phenotype_name;
    std::set<std::string>   m_base_changes; // ISBT naming like 261delG or 85A>G, ...
    std::set<std::string>   m_acid_changes; // ISBT naming like Ile60Leu or Ser230Ile, ...
    float                   m_incidence;
    

};

#endif /* CISBTPTALLELE_H */

