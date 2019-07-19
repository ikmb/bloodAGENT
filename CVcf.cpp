/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVcf.cpp
 * Author: mwittig
 * 
 * Created on July 6, 2018, 7:52 PM
 */

#include <cstdlib>
#include <vector>
#include <string>

#include "CVcf.h"

using namespace std;


CVcf::CVcf(const string& filename, bool verbose ) 
{
    m_inf = NULL;
    m_hdr = NULL;
    m_rec = bcf_init();
    m_verbose = verbose;
    if(open(filename) && read_header())
        read_sequences();
}

CVcf::CVcf(const CVcf& orig) 
{
    memcpy( m_inf, orig.m_inf, sizeof(orig.m_inf) );
    memcpy( m_hdr, orig.m_hdr, sizeof(orig.m_hdr) );
    m_rec = bcf_init();
    m_verbose = orig.m_verbose;
}

CVcf::~CVcf() {
    //free(m_inf);
    //free(m_hdr);
    free(m_rec);
    
    bcf_hdr_destroy(m_hdr);
    bcf_close(m_inf);
}

bool CVcf::open(const std::string& filename)
{
    m_inf = bcf_open(filename.c_str(), "r");
    return (m_inf != NULL);
}

bool CVcf::read_header()
{
    if(m_inf)
        m_hdr = bcf_hdr_read(m_inf);
    if(m_hdr && m_verbose)
        fprintf(stderr, "INFORMATION: File contains %i samples\n", bcf_hdr_nsamples(m_hdr));
    return (m_hdr != NULL);
}

bool CVcf::read_sequences()
{
    // report names of all the sequences in the VCF file
    // bcf_hdr_seqnames returns a newly allocated array of pointers to the seq names
    // caller has to deallocate the array, but not the seqnames themselves; the number
    // of sequences is stored in the int pointer passed in as the second argument.
    // The id in each record can be used to index into the array to obtain the sequence
    // name
    m_seq_names.clear();
    if(m_hdr)
    {
        const char **seqnames = NULL;
        int nseq = 0;  // number of sequences
        seqnames = bcf_hdr_seqnames(m_hdr, &nseq);
        for (int i = 0; i < nseq; i++) {
            m_seq_names.push_back(bcf_hdr_id2name(m_hdr, i));
        }
        free(seqnames);
    }
    return !m_seq_names.empty();

}

int  CVcf::nsamples ()const
{
    if(!m_hdr)
        return -1;
    return bcf_hdr_nsamples(m_hdr);
}

bool CVcf::limit_to_sample(const std::string& samplename)
{
    // limit the VCF data to the sample name passed in
    if(m_hdr)
    {
        bcf_hdr_set_samples(m_hdr, samplename.c_str(), 0);
        if (bcf_hdr_nsamples(m_hdr) != 1) {
                fprintf(stderr, "ERROR: please limit to a single sample\n");
                return false;
        }
        return true;
    }
    return false;
}

bool CVcf::read_record()
{
    if(bcf_read(m_inf, m_hdr, m_rec) == 0) 
    {
        if (!bcf_is_snp(m_rec)) 
        {
             if(m_verbose)
                fprintf(stderr, "INFORMATION: Non SNP entry read\n");
            return true;
        }
        // filter data for each call
        int nfi_arr = 0;
        int nfi     = 0;
        int *fi     = NULL;
       
        // the bcf_get_format_int32 function does not appear to reallocate
        // the array it returns for each of the samples on each call. Just
        // needs to be freed in the end. First call to bcf_get_format_*
        // takes care of calling bcf_unpack, which fills the `d` member
        // of bcf1_t
        nfi = bcf_get_format_int32(m_hdr, m_rec, "FI", &fi, &nfi_arr);
        
        
        free(fi);
        return true;
    }
    return false;
}


CVcfSnp CVcf::get_record()
{
    return CVcfSnp(m_inf, m_hdr,m_seq_names,m_rec);
}


