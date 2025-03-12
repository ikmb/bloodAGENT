/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVcf.h
 * Author: mwittig
 *
 * Created on July 6, 2018, 7:52 PM
 */

#ifndef CVCF_H
#define CVCF_H

#include "vcf.h"

#include "meinetools.h"
#include "CVcfSnp.h"

class CVcf {
public:
    CVcf(const std::string& filename, bool verbose = false);
    CVcf(const CVcf& orig);
    virtual ~CVcf();
    
    bool open(const std::string& filename);
    bool read_header();
    bool read_sequences();
    bool read_record();
    int  nsamples ()const;
    
    bool limit_to_sample(const std::string& samplename);
    
    CVcfSnp get_record();
    
private:
    
    htsFile *m_inf;
    bcf_hdr_t *m_hdr;
    std::vector<std::string>    m_seq_names;
    bool                        m_verbose;
    
    bcf1_t *m_rec;  // struc for storing each record
};

#endif /* CVCF_H */

