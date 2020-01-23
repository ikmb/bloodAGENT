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
    CIsbtVariant(const std::string& lrg_anno, const std::string& refBase, const std::string& chrom, int pos, char strand, int vcfCoord, const std::string& vcfRef, const std::string& vcfAlt, bool are_ref_and_alt_switched_in_GRCh, const string& variation_type);
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
       
    std::string lrgReference()const{return m_lrg_reference;}
    std::string lrgAlternative()const{return m_lrg_alternative;}
    int         vcfCoordinate()const{return m_vcf_coordinate;}
    std::string vcfReference()const{return m_vcf_reference;}
    std::string vcfAlternative()const{return m_vcf_alternative;}
    std::string reference()const{if(m_strand == '-')return CMyTools::GetComplSequence(m_lrg_reference); return m_lrg_reference;}
    std::string alternative()const{if(m_strand == '-')return CMyTools::GetComplSequence(m_lrg_alternative);return m_lrg_alternative;}
    
    bool isInDel()const{return (m_variation_type.compare("del") == 0 || m_variation_type.compare("ins") == 0 || m_variation_type.compare("delins") == 0);}
    
    bool addCoverage(const CBigWigReader& bigWig);
    double getCoverage()const{if(m_coverage != m_coverage) return 0.0;return m_coverage;}
    bool isCovered(double limit = 0.0f)const;
    
    bool isRefNClikeGRChNC()const{return !m_are_ref_and_alt_switched_in_GRCh;}
    
    void setVerbose(bool value = true){verbose=value;}
private:
    
    std::string m_isbt_name;
    std::string m_lrg_position;
    std::string m_lrg_reference;
    std::string m_lrg_alternative;
    int         m_vcf_coordinate;
    std::string m_vcf_reference;
    std::string m_vcf_alternative;
    std::string m_variation_type;
    
    bool        m_are_ref_and_alt_switched_in_GRCh;
    
    std::string m_chromosome;
    int         m_position;
    char        m_strand;
    double      m_coverage;
    
    
    bool parseIsbtVariant();
    
    static bool verbose;
    

};

#endif /* CISBTVARIANT_H */

