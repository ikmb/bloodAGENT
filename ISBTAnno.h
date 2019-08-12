/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ISBTAnno.h
 * Author: mwittig
 *
 * Created on July 18, 2019, 2:05 PM
 */

#ifndef ISBTANNO_H
#define ISBTANNO_H


#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>

#include <regex>
#include <iterator>

#include "mytools.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"

class CISBTAnno {
public:
    
    ///              position,<original nuc,alternative nuc>
    /// pos has to be string due to intronic annotation (e.g. 1254-31))
    typedef CIsbtVariant variation;
    
    CISBTAnno(const std::string& filename);
    CISBTAnno(const CISBTAnno& orig);
    virtual ~CISBTAnno();
    
    friend std::ostream& operator<<(std::ostream& os, const CISBTAnno& me);
    
    /// return all variations annotated for this genomic position
    std::vector<std::string> getVariationsAt(std::string chrom, int pos)const;
    
    /// return ISBT variation for this SNP or empty string
    variation getCorrespondingIsbtVariation(CVcfSnp)const;
    
    /// return blood group system for this genomic position
    std::string getSystemAt(std::string chrom, int pos)const;
    
    /// return all variations between LRG sequence and std GRCh reference ...
    std::map<std::string,std::vector<CISBTAnno::variation> > getReferenceVariations();
    
    /// return strand +/-/u
    char strand(const std::string& system)const{if(m_strand.find(system) != m_strand.end())return m_strand.find(system)->second;return 'u';}
    
    
    bool isVcfAlleleAnIsbtVariant(const std::string& allele, const std::string isbtVariant, const std::string system);
    
    std::set<std::string>   loci()const{ return m_loci;}
    
    
    variation getIsbtVariant(const std::string& system,const std::string& isbt_var)const;
    
    std::vector<variation>              getAllVariations()const{return m_parsed_isbt_variant;}
    size_t                              size()const{return m_parsed_isbt_variant.size();}
    variation                           variant(size_t idx)const{if(idx < size())return m_parsed_isbt_variant[idx]; return variation();}
    
    bool addCoverage(const CBigWigReader& bigWig);
    
private:
    bool m_data_red;
    bool readAnnotation(const std::string& filename);
    
    
    
    CParsedTextfile m_vanno;
    std::multimap<std::string,int>                      m_entry_finder; // by chromosomal coordinate, e.g. "chr9_1234567"
    std::map<std::string,char>                          m_strand;
    std::vector<variation>                              m_parsed_isbt_variant;
    std::map<std::string,std::map<std::string,int>>     m_isbt_variant_to_index;
    std::set<std::string>                               m_loci;
    
};

#endif /* ISBTANNO_H */

