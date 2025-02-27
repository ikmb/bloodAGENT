/* 
 * File:   CRefGeneEntry.cpp
 * Author: mwittig
 * 
 * Created on May 22, 2017, 4:05 PM
 */

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>


#include "CMyException.h"
#include "MyTools.h"
#include "CTwoBit.h"
#include "CRefGeneEntry.h"

/*
	0 #bin
	1 name
	2 chrom
	3 strand
	4 txStart
	5 txEnd
	6 cdsStart
	7 cdsEnd
	8 exonCount
	9 exonStarts
	0 exonEnds
	1 score
	2 name2
	3 cdsStartStat
	4 cdsEndStat
	5 exonFrames
*/

CRefGeneEntry::CRefGeneEntry(vector<string>& entries) 
{
    
    
    m_bin           = stol(entries[0]);
    m_name          = entries[1];
    m_chrom         = entries[2];    
    m_strand        = entries[3][0];
    m_txStart       = stol(entries[4]);
    m_txEnd         = stol(entries[5]);
    m_cdsStart      = stol(entries[6]);
    m_cdsEnd        = stol(entries[7]);
    m_exonCount     = stol(entries[8]);
    vector<string> sub_entries = CMyTools::GetParsedLine(entries[9],",");
    for(vector<string>::iterator i = sub_entries.begin(); i != sub_entries.end();i++)
        if( !i->empty() && i->find_first_not_of( "0123456789" ) == string::npos)
            m_exonStarts.push_back(stol(*i));
    sub_entries = CMyTools::GetParsedLine(entries[10],",");
    for(vector<string>::iterator i = sub_entries.begin(); i != sub_entries.end();i++)
        if( !i->empty() && i->find_first_not_of( "0123456789" ) == string::npos)
            m_exonEnds.push_back(stol(*i));
    m_score         = stol(entries[11]);
    m_name2         = entries[12];
    m_cdsStartStat  = entries[13];
    m_cdsEndStat    = entries[14];
    sub_entries = CMyTools::GetParsedLine(entries[10],",");
    for(vector<string>::iterator i = sub_entries.begin(); i != sub_entries.end();i++)
        if( !i->empty() && i->find_first_not_of( "0123456789" ) == string::npos)
            m_exonFrames.push_back(stoi(*i));
    
    if(m_exonStarts.size() != m_exonEnds.size())
        throw(CMyException(string("Error in CRefGeneEntry(), list sizes of exonStarts and exonEnds differ.")));
}

CRefGeneEntry::CRefGeneEntry(const CRefGeneEntry& orig) 
{
    m_bin = orig.m_bin;
    m_name = orig.m_name;
    m_chrom = orig.m_chrom;
    m_strand = orig.m_strand;
    m_txStart = orig.m_txStart;
    m_txEnd = orig.m_txEnd;
    m_cdsStart = orig.m_cdsStart;
    m_cdsEnd = orig.m_cdsEnd;
    m_exonCount = orig.m_exonCount;
    m_exonStarts = orig.m_exonStarts;
    m_exonEnds = orig.m_exonEnds;
    m_score = orig.m_score;
    m_name2 = orig.m_name2;
    m_cdsStartStat = orig.m_cdsStartStat;
    m_cdsEndStat = orig.m_cdsEndStat;
    m_exonFrames = orig.m_exonFrames;

}

CRefGeneEntry::~CRefGeneEntry() 
{
    
}

CRefGeneEntry::mutation CRefGeneEntry::getMutation(long idx,char alt, CTwoBit& tb)
{
    mutation mRet;
    long aa_start = idx%3;
    long triplet_pos = aa_start-1;
    if(triplet_pos==-1)
        triplet_pos=2;
    
    long aa_nr = idx/3;
    
    if(aa_start == 0)
    {
       aa_start = idx-2;
    }
    else
    {
        aa_start = aa_nr*3+1;
        aa_nr++;
    }
    
    string triplet="";
    for(int i = 0; i < 3; i++)
        triplet+=tb.getSequence(m_chrom.c_str(),getGenomicCoordinate(aa_start+i),getGenomicCoordinate(aa_start+i));
    if(isReverseStrand())
        for(int i = 0; i < 3; i++)
            triplet[i] = CMyTools::GetComplNuc(triplet[i]);
    
    string alt_triplet = triplet;
    alt_triplet[triplet_pos] = alt;
    pair<string,char> oldAA = CMyTools::getAminoAcid(triplet);
    pair<string,char> newAA = CMyTools::getAminoAcid(alt_triplet);
    
    mRet.ref_allele = triplet[triplet_pos];
    mRet.alt_allele = alt_triplet[triplet_pos];
    mRet.cDNA_position = idx;
    mRet.ref_amino = oldAA;
    mRet.alt_amino = newAA;
    mRet.amino_position = aa_nr;
    
    return mRet;
}

pair<long,long> CRefGeneEntry::getCDNAcoordinate(long pos)const
{
    pair<long,long> pRet(0,0);
    if(isForwardStrand())
    {
        if(pos < m_cdsStart)
        {
            pRet.second = m_cdsStart-pos;
        }
        else
        {
            vector<long>::const_iterator es = m_exonStarts.begin();
            vector<long>::const_iterator ee = m_exonEnds.begin();
            long start = *es;
            // move to cds start
            while(start < m_cdsStart && es != m_exonStarts.end())
            {
                start = *es;
                while(start < m_cdsStart && start < *ee)
                    start++;
                if(start == m_cdsStart)
                    break;
                es++;ee++;
            }
            if(start == m_cdsStart && es != m_exonStarts.end())
            {
                while(start < m_cdsEnd && es != m_exonStarts.end() && start < pos)
                {
                    pRet.second = 0;
                    while(start < m_cdsEnd && start < *ee && start < pos)
                    {    
                        pRet.first++;
                        start++;
                    }
                    es++;
                    if(es == m_exonStarts.end() || pos == start)
                        break;
                    while(start < m_cdsEnd && start < *es && start < pos)
                    {    
                        pRet.second++;
                        start++;
                    }
                    ee++;
                }
            }
            if(pos > m_cdsEnd)
            {
                pRet.second=pos-m_cdsEnd;
            }
                
        }
    }
    else
    {
        if(pos > m_cdsEnd)
        {
            pRet.second = pos-m_cdsEnd;
        }
        else
        {
            vector<long>::const_reverse_iterator es = m_exonStarts.rbegin();
            vector<long>::const_reverse_iterator ee = m_exonEnds.rbegin();
            long start = *ee;
            // move to cds start
            while(start > m_cdsEnd && es != m_exonStarts.rend())
            {
                start = *ee;
                while(start > m_cdsEnd && start > *es)
                    start--;
                if(start == m_cdsEnd)
                    break;
                es++;ee++;
            }
            if(start == m_cdsEnd && es != m_exonStarts.rend())
            {
                while(start > m_cdsStart && es != m_exonStarts.rend() && start >= pos)
                {
                    pRet.second = 0;
                    while(start > m_cdsStart && start > *es && start >= pos)
                    {    
                        pRet.first++;
                        start--;
                    }
                    ee++;
                    if(ee == m_exonEnds.rend() || pos == start)
                        break;
                    while(start > m_cdsStart && start > *ee && start >= pos)
                    {    
                        pRet.second++;
                        start--;
                    }
                    es++;
                }
            }
            if(pos < m_cdsStart)
            {
                pRet.second=m_cdsStart-pos;
            }
                
        }
    }    
    
    return pRet;
}


std::string CRefGeneEntry::getSector(const string& chrom, const long& position)
{
    if(chrom.compare(m_chrom) != 0)
        return "NA";
    if(position > m_txEnd)
        return (isForwardStrand() ? "3'" : "5'");
    if(position <= m_txStart)
        return (isForwardStrand() ? "5'" : "3'");
    if(isForwardStrand())
    {
        int i = 1;
        vector<long>::iterator se = m_exonEnds.begin();
        for(vector<long>::iterator ss = m_exonStarts.begin(); ss != m_exonStarts.end(); ss++, se++, i++)
            if(position > *ss && position <= *se)
                return string("exon ")+std::to_string(i);
    }
    else
    {
        int i = 1;
        vector<long>::reverse_iterator se = m_exonEnds.rbegin();
        for(vector<long>::reverse_iterator ss = m_exonStarts.rbegin(); ss != m_exonStarts.rend(); ss++, se++, i++)
            if(position > *ss && position <= *se)
                return string("exon ")+std::to_string(i);
    }
    if(isForwardStrand() && m_exonStarts.size() > 1)
    {
        int i = 1;
        vector<long>::iterator ss = m_exonStarts.begin()+1;
        for(vector<long>::iterator se = m_exonEnds.begin(); ss != m_exonStarts.end(); ss++, se++, i++)
            if(position > *se && position <= *ss)
                return string("intron ")+std::to_string(i);
    }
    else
    {
        int i = 1;
        vector<long>::reverse_iterator se = m_exonEnds.rbegin()+1;
        for(vector<long>::reverse_iterator ss = m_exonStarts.rbegin(); se != m_exonEnds.rend(); ss++, se++, i++)
            if(position > *se && position <= *ss)
                return string("intron ")+std::to_string(i);
    }
    return "NA";
}

std::string CRefGeneEntry::getSectors()
{
    ostringstream strRet("");
    
    if(isForwardStrand())
    {
        int i = 1;
        ostringstream intron_helper("");
        vector<long>::iterator se = m_exonEnds.begin();
        for(vector<long>::iterator ss = m_exonStarts.begin(); ss != m_exonStarts.end(); ss++, se++, i++)
        {
            if(intron_helper.str().size() != 0)
            {
                strRet << intron_helper.str() << *ss-1 << endl;
                intron_helper.str("");
            }
            strRet << m_name << "\t" << m_name2 << "\t" << "exon_" << i << "\t"
                   << m_chrom << "\t" << *ss << "\t" << *se << endl;
            intron_helper << m_name << "\t" << m_name2 << "\t" << "intron_" << i << "\t"
                   << m_chrom << "\t" << *se+1 << "\t";
        }
    }
    else
    {
        int i = 1;
        ostringstream intron_helper("");
        vector<long>::reverse_iterator se = m_exonEnds.rbegin();
        for(vector<long>::reverse_iterator ss = m_exonStarts.rbegin(); ss != m_exonStarts.rend(); ss++, se++, i++)
        {
            if(intron_helper.str().size() != 0)
            {
                strRet << m_name << "\t" << m_name2 << "\t" << "intron_" << i-1 << "\t"
                   << m_chrom << "\t" << *se+1 << "\t" << intron_helper.str() << endl;
                intron_helper.str("");
            }
            strRet << m_name << "\t" << m_name2 << "\t" << "exon_" << i << "\t"
                   << m_chrom << "\t" << *ss << "\t" << *se << endl;
            intron_helper << *ss-1;
        }
    }
    return strRet.str();
}

long CRefGeneEntry::getGenomicCoordinate(long idx)const
{
    if(isForwardStrand())
    {
        int iExon = 0;
        while(m_exonStarts.size() > static_cast<std::size_t>(iExon) && m_exonStarts[static_cast<std::size_t>(iExon)] > m_cdsStart)
            iExon++;
        if( !(m_exonStarts.size() > iExon))
            throw(CMyException(string("Error in CRefGeneEntry::getCDNAat(long idx), index out of bounds.")));
        
        //long current_start = m_cdsStart;
        long cDNA_passed = m_exonEnds[iExon] - m_cdsStart;
        if(cDNA_passed > idx)
            return m_cdsStart+idx;
        else
        {
            iExon++;
            while(m_exonStarts.size() > iExon)
            {
                idx-=cDNA_passed;
                cDNA_passed = m_exonEnds[iExon] - m_exonStarts[iExon];
                if(cDNA_passed > idx)
                    return m_exonStarts[iExon]+idx;
                iExon++;
            }
        }
    }
    else
    {
        idx--;
        int iExon = m_exonStarts.size()-1;
        while(iExon >= 0 && m_exonEnds[iExon] < m_cdsEnd)
            iExon--;
        if( iExon < 0)
            throw(CMyException(string("Error in CRefGeneEntry::getCDNAat(long idx), index out of bounds.")));
        
        //long current_start = m_cdsEnd;
        long cDNA_passed = m_cdsEnd - m_exonStarts[iExon];
        if(cDNA_passed > idx)
            return m_cdsEnd-idx;
        else
        {
            iExon--;
            while(iExon >= 0)
            {
                idx-=cDNA_passed;
                cDNA_passed = m_exonEnds[iExon] - m_exonStarts[iExon];
                if(cDNA_passed > idx)
                    return m_exonEnds[iExon]-idx;
                iExon--;
            }
        }
    }
    return string::npos;
}


