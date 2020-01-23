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
    bool ready()const{return (m_fp != NULL);}
    
    
    double getMinCoverage(const std::string& chrom, int start, int end)const{return getCoverage(chrom, start, end, bwStatsType::min);}
    double getAverageCoverage(const std::string& chrom, int start, int end)const{return getCoverage(chrom, start, end, bwStatsType::average);}
    double getMaxCoverage(const std::string& chrom, int start, int end)const{return getCoverage(chrom, start, end, bwStatsType::max);}
    double getSumCoverage(const std::string& chrom, int start, int end)const{return getCoverage(chrom, start, end, bwStatsType::sum);}
    
    double getPercentCoveredBases(const std::string& chrom, int start, int end);
    
    bool isCovered(const std::string& chrom, int start, int end, double min_required = 0.5);
    
private:
    
    bigWigFile_t*               m_fp;
    double getCoverage(const std::string& chrom, int start, int end, bwStatsType type)const;
};

#endif /* CBIGWIGREADER_H */

