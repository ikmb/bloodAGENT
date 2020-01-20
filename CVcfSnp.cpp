/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVcfSnp.cpp
 * Author: mwittig
 * 
 * Created on July 9, 2018, 9:02 AM
 */
#include <string>
#include <vector>
#include <ostream>
#include <sstream>
#include <iostream>

#include "vcf.h"

#include "CVcfSnp.h"

// https://github.com/tonynaggs/htslib
// http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html

using namespace std;

CVcfSnp::CVcfSnp(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& seq_names,bcf1_t *rec, bool verbose) 
{
    m_chrom = "";
    m_pos = -1;
    m_depth = -1;
    m_verbose = verbose;
    m_phasing_id = -1;
    m_ref_allele="";
    read_SNP_entry(inf, hdr,seq_names,rec);
}

CVcfSnp::CVcfSnp(const CVcfSnp& orig) 
{
    m_chrom     = orig.m_chrom;
    m_pos       = orig.m_pos;
    m_alleles   = orig.m_alleles;
    m_coverage  = orig.m_coverage;
    m_qualities = orig.m_qualities;
    m_depth     = orig.m_depth;
    m_verbose   = orig.m_verbose;
    m_phasing_id = orig.m_phasing_id;
    m_ref_allele = orig.m_ref_allele;
}

CVcfSnp::~CVcfSnp() {
}

 bool CVcfSnp::isHomozygous()const
 {
     if(m_alleles.size() == 0)
         return false;
     for(size_t i = 1; i < m_alleles.size(); i++)
         if(m_alleles[i].compare(m_alleles[i-1]) != 0)
             return false;
     return true;
 }

std::vector<std::string>  CVcfSnp::indelalleles()const
{
    // ToDO:
    // Seltener Fall eineszum SNP benachbarten indels mit Freebyas gecallt
    // Soll sein G/A, reportiert aber mit Freebayes als:
    // 6       31106499        .       GTCCCCCCCA      ATCCCCCCA       271.543
    vector<std::string> vRet = m_alleles;
    if(vRet.size() > 1)
    {
        string hlp = m_ref_allele;
        size_t counter = 0;
        while(counter < hlp.size())
        {
            bool remove = true;
            for(size_t i = 0; i < vRet.size(); i++)
            {
                if(vRet[i].size() == 0 || vRet[i][0] != hlp[counter])
                {
                    remove = false;
                    break;
                }
            }
            if(remove)
                for(size_t i = 0; i < vRet.size(); i++)
                    vRet[i]=vRet[i].substr(1);
            counter++;
        };
    }
    for(size_t i = 0; i < vRet.size(); i++)
        if(vRet[i].size() == 0)
            vRet[i]="-";
    return vRet;
}

void CVcfSnp::read_SNP_entry(htsFile *inf, bcf_hdr_t *hdr,std::vector<std::string>& seq_names,bcf1_t *rec)
{
    // quality data for each call
    int ngq_arr = 0;
    int ngq     = 0;
    int *gq     = NULL;
    // genotype data for each call
    // genotype arrays are twice as large as
    // the other arrays as there are two values for each sample
    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    
    
    int nad_arr = 0;
    int nad     = 0;
    int *ad     = NULL;

    // depth
    int ndp_arr = 0;
    int ndp     = 0;
    int *dp     = NULL;
    
    
    // Phased
    int nps_arr = 0;
    int nps       = 0;
    int *ps     = NULL;

    ngq = bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq_arr);
    ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
    nad = bcf_get_format_int32(hdr, rec, "AD", &ad, &nad_arr);
    ndp = bcf_get_info_int32(hdr, rec, "DP", &dp, &ndp_arr);
    nps = bcf_get_format_int32(hdr, rec, "PS", &ps, &nps_arr);
    
    

    if(rec->rid < static_cast<int32_t>(seq_names.size()))
        m_chrom = seq_names[rec->rid];
    m_pos = rec->pos+1;                 // vcf LIB INTERNALLY STORES 0-BASED POSITIONS
    if(ndp > 0)
        m_depth = dp[0];
    
    for(int i = 0; i < ngt; i++)
    {
        m_alleles.push_back(rec->d.allele[bcf_gt_allele(gt[i])]);
    }
    cout << endl;
    for(int i = 0; i < nad; i++)
        m_coverage.push_back(ad[i]);
    for(int i = 0; i < ngq; i++)
        m_qualities.push_back(gq[i]);
    if(nps > 0)   
        m_phasing_id = ps[0];
    
    m_ref_allele = rec->d.allele[0];    
    
    
    free(gq);
    free(gt);
    free(ad);
    free(dp);
    free(ps);
}

std::string    CVcfSnp::SNP()const
{
    ostringstream osr("");
    for(size_t i = 0; i < m_alleles.size(); i++)
        osr << m_alleles[i] << ( i + 1 == m_alleles.size() ? '\t' : ( isPhased() ? '|' : '/'));
    return osr.str();
}

std::ostream& operator<<(std::ostream& os, const CVcfSnp& me)
{
    os << me.m_chrom << '\t' << me.m_pos << '\t' << me.m_depth << '\t';
    
    for(size_t i = 0; i < me.m_alleles.size(); i++)
        os << me.m_alleles[i] << ( i + 1 == me.m_alleles.size() ? '\t' : ( me.isPhased() ? '|' : '/'));
    
    for(size_t i = 0; i < me.m_coverage.size(); i++)
        os << me.m_coverage[i] << ( i + 1 == me.m_coverage.size() ? '\t' : '/');
    
    for(size_t i = 0; i < me.m_qualities.size(); i++)
        os << me.m_qualities[i] << ( i + 1 == me.m_qualities.size() ? '\t' : '/');
    
    if(me.isPhased())
        os << "phased\t" << me.phasingID();
    
    return os;
}

bool CVcfSnp::isInRegion(const std::string& chrom, long start, long end)
{
    if(m_chrom.compare(chrom) != 0)
        return false;
    if(m_pos >= start && m_pos <= end)
        return true;
    return false;
}

bool CVcfSnp::isInRegion(const std::vector<std::string>& chrom, std::vector<long> start, std::vector<long> end)
{
    if(chrom.size() != start.size() || start.size() != end.size())
    {
        fprintf(stderr, "ERROR: Function \"bool CVcfSnp::isInRegion(const std::vector<std::string>& chrom, std::vector<long> start, std::vector<long> end)\" requires all vectors having the same number of entries.\n");
        return false;
    }
    
    for(size_t i = 0; i < chrom.size(); i++)
        if(isInRegion(chrom[i], start[i],end[i]))
            return true;
    return false;
}


std::string CVcfSnp::getPosKey()const
{
    ostringstream osr;
    osr << m_chrom << '_' << m_pos;
    return osr.str();
}
