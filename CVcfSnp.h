/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVcfSnp.h
 * Author: mwittig
 *
 * Created on July 9, 2018, 9:02 AM
 */

#ifndef CVCFSNP_H
#define CVCFSNP_H

#include <string>


class CVcfSnp {
public:
    CVcfSnp();
    CVcfSnp(std::string chrom,  long pos,  std::vector<std::string> alleles,  std::vector<int> coverage,  std::vector<int> qualities,  std::vector<int> haplotype_qualities,  int mapping_quality,  long depth,  bool verbose,  int phasing_id,  std::string ref_allele);
    CVcfSnp(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& seq_names,bcf1_t *rec, bool verbose = false);
    CVcfSnp(const CVcfSnp& orig);
    CVcfSnp(const CVcfSnp* orig);
    virtual ~CVcfSnp();
    
    friend std::ostream& operator<<(std::ostream& os, const CVcfSnp& me);
    
    bool isInRegion(const std::string& chrom, long start, long end);
    bool isInRegion(const std::vector<std::string>& chrom, std::vector<long> start, std::vector<long> end);
    
    std::string getPosKey()const;
    
    
    std::string                 chrom()const{return m_chrom;}
    long                        pos()const{return m_pos;}
    std::vector<std::string>    alleles()const{return m_alleles;}
    std::vector<std::string>    indelalleles()const;
    std::string                 refAllele()const{return m_ref_allele;}
    std::string                 SNP()const;
    int                         phasingID()const{return m_phasing_id;}
    bool                        isPhased()const{return m_phasing_id != -1;}
    
    bool                        isHomozygous()const;
    bool                        isHeterozygous()const{return !isHomozygous();};
    
    std::vector<int>            genotypeQualities()const{return m_qualities;}
    std::vector<int>            haplotypeQualities()const{return m_haplotype_qualities;}
    int                         mappingQuality()const{return m_mapping_quality;}
    
    
private:
    
    std::string                 m_chrom;
    long                        m_pos;
    std::vector<std::string>    m_alleles;
    std::vector<int>            m_coverage;
    std::vector<int>            m_qualities;
    std::vector<int>            m_haplotype_qualities;
    int                         m_mapping_quality;
    long                        m_depth;
    bool                        m_verbose;
    int                         m_phasing_id;
    std::string                 m_ref_allele;
    
    void read_SNP_entry(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& m_seq_names,bcf1_t *rec);
};

#endif /* CVCFSNP_H */

