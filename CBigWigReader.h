/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CBigWigReader.h
 * Author: mwittig
 *
 * Created on August 7, 2019, 9:52 AM
 */

#ifndef CBIGWIGREADER_H
#define CBIGWIGREADER_H

#include "bigWig.h"
#include <limits>

class CBigWigReader {
public:
    static int                  m_instance_counter;
  
    CBigWigReader(const std::string& filename);
    CBigWigReader(const CBigWigReader& orig);
    virtual ~CBigWigReader();
    
    bool init(const std::string& filename);
    bool initBam(const std::string& filename);
    bool ready()const;
    
    double getCoverageFromBam(const std::string& chrom, int start, int end, bwStatsType type) const;
    double getPercentCoveredFromBam(const std::string& chrom, int start, int end) const;

    
    double getMinCoverage(const std::string& chrom, int start, int end)const;
    double getAverageCoverage(const std::string& chrom, int start, int end)const;
    double getMaxCoverage(const std::string& chrom, int start, int end)const;
    double getSumCoverage(const std::string& chrom, int start, int end)const;
    
    double getPercentCoveredBases(const std::string& chrom, int start, int end);
    
    bool isCovered(const std::string& chrom, int start, int end, double min_required = 0.5);
    
private:
    
    bigWigFile_t*               m_fp;
    bool       m_useBam = false;      // true wenn BAM-Modus
    samFile*   m_bam = nullptr;    // BAM handle
    bam_hdr_t* m_bamHdr = nullptr; // BAM Header
    hts_idx_t* m_bamIdx = nullptr; // BAM Index
    double getCoverage(const std::string& chrom, int start, int end, bwStatsType type)const;
};

#endif /* CBIGWIGREADER_H */

