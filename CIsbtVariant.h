/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtVariant.h
 * Author: mwittig
 *
 * Created on July 22, 2019, 8:08 AM
 */

#ifndef CISBTVARIANT_H
#define CISBTVARIANT_H

#include "mytools.h"
#include "CBigWigReader.h"



class CIsbtVariant {
public:
    CIsbtVariant();
    CIsbtVariant(std::string lrg_anno, std::string refBase, std::string chrom, int pos, char strand);
    CIsbtVariant(const CIsbtVariant& orig);
    
    CIsbtVariant& operator =(const CIsbtVariant& orig);
    bool          operator <(const CIsbtVariant& orig)const;
    bool          operator >(const CIsbtVariant& orig)const{return !(*this < orig || *this == orig);}
    bool          operator ==(const CIsbtVariant& orig)const;
    //bool          operator ==(const string& orig)const;
    bool          operator !=(const CIsbtVariant& orig)const{return !(*this == orig);};
    friend std::ostream& operator<<(std::ostream& os, const CIsbtVariant& me);
    
    virtual ~CIsbtVariant();
    
    char strand(){return m_strand;}
    
    std::string name()const{return m_isbt_name;}
    std::string posLRG()const{return m_lrg_position;}
        
    std::string chrom()const{return m_chromosome;}
    int         pos()const{return m_position;}
       
    std::string lrgReference(){return m_lrg_reference;}
    std::string lrgAlternative(){return m_lrg_alternative;}
    std::string reference()const{if(m_strand == '-')return CMyTools::GetComplSequence(m_lrg_reference); return m_lrg_reference;}
    std::string alternative()const{if(m_strand == '-')return CMyTools::GetComplSequence(m_lrg_alternative);return m_lrg_alternative;}
    
    bool isInDel()const{return (m_lrg_reference.compare("-") == 0 || m_lrg_alternative.compare("-") == 0);}
    
    bool addCoverage(const CBigWigReader& bigWig);
    double getCoverage()const{return m_coverage;}
    bool isCovered(double limit = 0.0f)const;
    
    void setVerbose(bool value = true){verbose=value;}
private:
    
    std::string m_isbt_name;
    std::string m_lrg_position;
    std::string m_lrg_reference;
    std::string m_lrg_alternative;
    
    std::string m_chromosome;
    int         m_position;
    char        m_strand;
    double      m_coverage;
    
    
    bool parseIsbtVariant();
    
    static bool verbose;
    

};

#endif /* CISBTVARIANT_H */

