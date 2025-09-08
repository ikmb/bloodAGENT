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
#include <htslib/sam.h>

#include "meinetools.h"
#include "CBigWigReader.h"

using namespace std;

int CBigWigReader::m_instance_counter = 0;

CBigWigReader::CBigWigReader(const std::string& filename) 
{
    m_fp = NULL;
    if(!CMyTools::file_exists(filename))
        throw(CMyException("File does not exist: ")+filename);
    // BAM erkennen
    if (filename.size() >= 4 && filename.substr(filename.size()-4) == ".bam") {
        if (!initBam(filename)) throw(CMyException("Failed to open BAM: ")+filename);
        m_useBam = true;
    } else {
        if (!init(filename)) throw(CMyException("Failed to open BigWig: ")+filename);
        m_useBam = false;
    }
}

CBigWigReader::CBigWigReader(const CBigWigReader& orig) 
{
    m_fp = NULL;
    m_useBam = false;
    m_bam = nullptr; m_bamHdr = nullptr; m_bamIdx = nullptr;
    // memcpy( m_fp, orig.m_fp, sizeof(*orig.m_fp) );
}

CBigWigReader::~CBigWigReader() 
{
    if(m_fp) bwClose(m_fp);
    if(m_bamIdx) hts_idx_destroy(m_bamIdx);
    if(m_bamHdr) bam_hdr_destroy(m_bamHdr);
    if(m_bam) sam_close(m_bam);
    if(--m_instance_counter == 0)
        bwCleanup();
}

bool CBigWigReader::ready() const {
    if (m_useBam) return m_bam && m_bamHdr && m_bamIdx;
    return m_fp != nullptr;
}

double CBigWigReader::getCoverage(const std::string& chrom, int start, int end, bwStatsType type)const
{
    if (m_useBam) {
        return getCoverageFromBam(chrom, start, end, type);
    } 
    else {
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
        if(ready()) {
        string str = chrom;
            char *cstr = &str[0];
            double *stats = bwStatsFromFull(m_fp, cstr,start, end, 1, type);
            if(stats) {   
                bRet =  stats[0];
                free(stats);
            }
        }
        return bRet;
    }
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
    if (m_useBam) {
         return getPercentCoveredFromBam(chrom, start, end);
    } 
    else {
        double bRet = std::numeric_limits<double>::quiet_NaN();
        if(ready())
        {
            /*
             typedef struct {
            uint32_t l; //<Number of intervals held
            uint32_t m; //<Maximum number of values/intervals the struct can hold
            uint32_t *start; //<The start positions (0-based half open)
            uint32_t *end; //<The end positions (0-based half open)
            float *value; //<The value associated with each position
            } bwOverlappingIntervals_t;
             */
            string str = chrom;
            char *cstr = &str[0];
            bwOverlappingIntervals_t* intervals = bwGetValues(m_fp, cstr,start, end,0);
            if(intervals) 
            {
                int length = 0;
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

bool CBigWigReader::initBam(const std::string& filename)
{
    // BAM öffnen, Header lesen, Index laden
    m_bam = sam_open(filename.c_str(), "r");
    if(!m_bam) {
        std::cerr << "Failed to open BAM " << filename << std::endl;
        return false;
    }
    m_bamHdr = sam_hdr_read(m_bam);
    if(!m_bamHdr) {
        std::cerr << "Failed to read BAM header " << filename << std::endl;
        return false;
    }
    m_bamIdx = sam_index_load(m_bam, filename.c_str());
    if(!m_bamIdx) {
        std::cerr << "BAM index (.bai/.csi) not found for " << filename << std::endl;
        return false;
    }
    return true;
}

// -------- BAM-Coverage-Helfer --------
static inline void add_block_cov(std::vector<uint32_t>& cov, int32_t rpos, uint32_t len, int winStart, int winEnd) {
    if (len == 0) return;
    int32_t s = rpos;
    int32_t e = rpos + (int32_t)len;
    if (e <= winStart || s >= winEnd) return;
    if (s < winStart) s = winStart;
    if (e > winEnd)   e = winEnd;
    for (int32_t i = s; i < e; ++i) cov[(size_t)(i - winStart)] += 1u;
}

double CBigWigReader::getCoverageFromBam(const std::string& chrom, int start, int end, bwStatsType type) const {
    if (!ready() || start >= end) return std::numeric_limits<double>::quiet_NaN();
    int tid = bam_name2id(m_bamHdr, chrom.c_str());
    if (tid < 0) return std::numeric_limits<double>::quiet_NaN();

    const int winStart = start;
    const int winEnd   = end;
    const size_t winLen = (size_t)(winEnd - winStart);
    std::vector<uint32_t> cov(winLen, 0u);

    hts_itr_t* itr = sam_itr_queryi(m_bamIdx, tid, winStart, winEnd);
    if (!itr) return std::numeric_limits<double>::quiet_NaN();
    bam1_t* b = bam_init1();
    while (sam_itr_next(m_bam, itr, b) >= 0) {
        const bam1_core_t& c = b->core;
        if (c.tid != tid) continue;
        int32_t rpos = c.pos; // 0-based
        uint32_t* cigar = bam_get_cigar(b);
        for (uint32_t k = 0; k < c.n_cigar; ++k) {
            uint32_t op  = bam_cigar_op(cigar[k]);
            uint32_t len = bam_cigar_oplen(cigar[k]);
            switch (op) {
                case BAM_CMATCH:    // M
                case BAM_CEQUAL:    // =
                case BAM_CDIFF:     // X
                    add_block_cov(cov, rpos, len, winStart, winEnd);
                    rpos += len; break;
                case BAM_CDEL:
                case BAM_CREF_SKIP:
                    rpos += len; break;
                case BAM_CINS:
                case BAM_CSOFT_CLIP:
                case BAM_CHARD_CLIP:
                case BAM_CPAD:
                default:
                    // keine Referenzfortschreibung
                    break;
            }
        }
    }
    bam_destroy1(b);
    hts_itr_destroy(itr);

    // Statistiken
    if (type == bwStatsType::cov) {
        // Anzahl bedeckter Basen (>0)
        size_t covered = 0;
        for (auto v: cov) if (v > 0) ++covered;
        return static_cast<double>(covered);
    }
    // sum/min/max/mean über per-base coverage
    uint64_t sum = 0;
    uint32_t vmax = 0;
    uint32_t vmin = std::numeric_limits<uint32_t>::max();
    for (auto v: cov) {
        sum += v;
        if (v > vmax) vmax = v;
        if (v < vmin) vmin = v;
    }
    if (winLen == 0) return std::numeric_limits<double>::quiet_NaN();
    switch (type) {
        case bwStatsType::sum:     
            return static_cast<double>(sum);
        case bwStatsType::max:     
            return static_cast<double>(vmax);
        case bwStatsType::min:     
            return static_cast<double>(vmin == std::numeric_limits<uint32_t>::max() ? 0.0 : (double)vmin);
        case bwStatsType::average: 
            return static_cast<double>(sum) / static_cast<double>(winLen);
        default:                   
            return std::numeric_limits<double>::quiet_NaN();
    }
}

double CBigWigReader::getPercentCoveredFromBam(const std::string& chrom, int start, int end) const {
    double covBases = getCoverageFromBam(chrom, start, end, bwStatsType::cov);
    if (!(covBases == covBases)) return std::numeric_limits<double>::quiet_NaN();
    double len = static_cast<double>(end - start);
    if (len <= 0) return std::numeric_limits<double>::quiet_NaN();
    return 100.0 * covBases / len;
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

    