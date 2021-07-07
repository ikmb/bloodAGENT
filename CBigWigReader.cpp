/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CBigWigReader.cpp
 * Author: mwittig
 * 
 * Created on August 7, 2019, 9:52 AM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>
#include <iostream>
#include <cstring>
#include <limits>

#include "CBigWigReader.h"

using namespace std;

int CBigWigReader::m_instance_counter = 0;

CBigWigReader::CBigWigReader(const std::string& filename) 
{
    m_fp = NULL;
    init(filename);
}

CBigWigReader::CBigWigReader(const CBigWigReader& orig) 
{
    m_fp = NULL;
    memcpy( m_fp, orig.m_fp, sizeof(*orig.m_fp) );
}

CBigWigReader::~CBigWigReader() 
{
    if(m_fp)
    {
        bwClose(m_fp);
    }
    if(--m_instance_counter == 0)
        bwCleanup();
}
double CBigWigReader::getCoverage(const std::string& chrom, int start, int end, bwStatsType type)const
{
    /// enum bwStatsType: 
    /// doesNotExist = -1, /*!< This does nothing */
    /// mean = 0, /*!< The mean value */
    /// average = 0, /*!< The mean value */
    /// stdev = 1, /*!< The standard deviation of the values */
    /// dev = 1, /*!< The standard deviation of the values */
    /// max = 2, /*!< The maximum value */
    /// min = 3, /*!< The minimum value */
    /// cov = 4, /*!< The number of bases covered */
    /// coverage = 4, /*!<The number of bases covered */ 
    /// sum = 5 /*!< The sum of per-base values */
    double bRet = std::numeric_limits<float>::quiet_NaN();
    if(ready())
    {
        string str = chrom;
        char *cstr = &str[0];
        double *stats = bwStatsFromFull(m_fp, cstr,start, end, 1, type);
        if(stats) 
        {
            bRet =  stats[0];
            free(stats);
        }
    }
    return bRet;
}
bool CBigWigReader::isCovered(const std::string& chrom, int start, int end, double min_required)
{
    double bVal = getAverageCoverage( chrom, start, end);
    if(bVal != bVal || bVal < min_required)
        return false;
    return true;
}


double CBigWigReader::getPercentCoveredBases(const std::string& chrom, int start, int end)
{
    double bRet = std::numeric_limits<double>::quiet_NaN();
    if(ready())
    {
        string str = chrom;
        char *cstr = &str[0];
        bwOverlappingIntervals_t* intervals = bwGetValues(m_fp, cstr,start, end,0);
        
         typedef struct {
            uint32_t l; /**<Number of intervals held*/
            uint32_t m; /**<Maximum number of values/intervals the struct can hold*/
            uint32_t *start; /**<The start positions (0-based half open)*/
            uint32_t *end; /**<The end positions (0-based half open)*/
            float *value; /**<The value associated with each position*/
        } bwOverlappingIntervals_t;

        
        if(intervals) 
        {
            int     length = 0;
            for(uint32_t i = 0; i < intervals->l; i++)
            {
                if(intervals->value[i] == intervals->value[i])
                    length++;
            }
            bRet = 100.0*static_cast<double>(length)/static_cast<double>(end-start);
            bwDestroyOverlappingIntervals(intervals);
        }
    }
    return bRet;
}

bool CBigWigReader::init(const std::string& filename)
{
    //Initialize enough space to hold 128KiB (1<<17) of data at a time
    if(m_instance_counter++ == 0)
    {
        if(bwInit(1<<17) != 0)
            return false;
    }
    
    string str = filename;
    char *cstr = &str[0];
    
    if(!bwIsBigWig(cstr, NULL))
    {
        cerr << "An error occured while opening " << filename << ", which does not seem to be a bigWig file."<< endl;
        return false;
    }
    
    m_fp = bwOpen(cstr, NULL, "r");
    if(!m_fp) {
        cerr << "An error occured while opening " << filename << endl;
        return false;
    }
    return true;
}
 
double CBigWigReader::getMinCoverage(const std::string& chrom, int start, int end)const
{
    return getCoverage(chrom, start, end, bwStatsType::min);
}

double CBigWigReader::getAverageCoverage(const std::string& chrom, int start, int end)const
{
    return getCoverage(chrom, start, end, bwStatsType::average);
}

double CBigWigReader::getMaxCoverage(const std::string& chrom, int start, int end)const
{
    return getCoverage(chrom, start, end, bwStatsType::max);
}

double CBigWigReader::getSumCoverage(const std::string& chrom, int start, int end)const
{
    return getCoverage(chrom, start, end, bwStatsType::sum);
}

    