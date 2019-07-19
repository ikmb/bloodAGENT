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
    CVcfSnp(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& seq_names,bcf1_t *rec, bool verbose = false);
    CVcfSnp(const CVcfSnp& orig);
    virtual ~CVcfSnp();
    
    friend std::ostream& operator<<(std::ostream& os, const CVcfSnp& me);
    
    bool isInRegion(const std::string& chrom, long start, long end);
    bool isInRegion(const std::vector<std::string>& chrom, std::vector<long> start, std::vector<long> end);
    
    std::string getPosKey()const;
    
    
    std::string                    chrom()const{return m_chrom;}
    long                             pos()const{return m_pos;}
    std::vector<std::string>     alleles()const{return m_alleles;}
    std::vector<std::string>     indelalleles()const;
    std::string                    SNP()const;
    size_t                     phasingID()const{return m_phasing_id;}
    bool                        isPhased()const{return m_phasing_id != -1;}
    
    
private:
    
    std::string                 m_chrom;
    long                        m_pos;
    std::vector<std::string>    m_alleles;
    std::vector<int>            m_coverage;
    std::vector<int>            m_qualities;
    long                        m_depth;
    bool                        m_verbose;
    int                      m_phasing_id;
    
    void read_SNP_entry(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& m_seq_names,bcf1_t *rec);
};

#endif /* CVCFSNP_H */

