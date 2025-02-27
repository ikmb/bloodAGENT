/* 
 * File:   CTwoBit.cpp
 * Author: mwittig
 * 
 * Created on February 18, 2011, 10:52 AM
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>

//#include "BGZF.h"
#include "CTwoBit.h"

using namespace std;

CTwoBit::CTwoBit() {
    initialize();
}

CTwoBit::CTwoBit(const char* filename) {
    initialize();
    m_bOpen = loadTwoBit(filename);

}

// copy constructor is incomplete (memcopy is missing!!!)
CTwoBit::CTwoBit(const CTwoBit& orig) {
    m_bOpen     = orig.m_bOpen;
    m_pFile     = orig.m_pFile;
    m_stcHeader = orig.m_stcHeader;
    m_data = orig.m_data;
}

CTwoBit::~CTwoBit() {
    for(map<string,chrom_data>::iterator iter = m_data.begin(); iter != m_data.end(); iter++)
        free(iter->second.packedDNA);
}

vector<pair<string,unsigned int> >  CTwoBit::getChromSizes()
{
    vector<pair<string,unsigned int> > vctReturn;

    for(map<string,chrom_data>::const_iterator iter = m_data.begin();iter != m_data.end();iter++)
        vctReturn.push_back( pair<string,unsigned int>(iter->first,iter->second.dnaSize) );

    return vctReturn;
}

bool CTwoBit::loadTwoBit(const char* filename)
{
    m_pFile = new ifstream(filename,ios_base::binary);
    if(!m_pFile || !m_pFile->is_open())
    {
        cerr << "failed to open 2bit file \'" << filename << '\'' << endl;
        return false;
    }
    m_pFile->seekg(ios::beg);
    m_pFile->read((char*)&m_stcHeader,16);
    if(m_stcHeader.signature != 0x1A412743)
    {
        cerr << "TwoBit file signature of \"" << filename << "\" is incorrect." << endl;
        return false;
    }
    if(!storeOffsets() || !storeData())
        return false;
    m_pFile->close();
    free(m_pFile);
    m_pFile = NULL;
    return false;
}

bool CTwoBit::storeOffsets()
{
    try
    {
        m_pFile->seekg(16,ios::beg);
        for(int idx = 0; idx < m_stcHeader.seqCount; idx++)
        {
            char name_length;
            m_pFile->read(&name_length,1);
            char* strName = (char*)malloc(name_length);
            m_pFile->read(strName,name_length);
            chrom_data initData = getChromDataStruc();
            m_pFile->read((char*)&initData.offset,4);
            m_data[string(strName).substr(0,name_length)]=initData;
            free(strName);
        }
        return true;
    }
    catch(exception ex)
    {
        cerr << ex.what() << endl;
    }
    catch(...)
    {

    }
    return false;
}

void CTwoBit::closeTwoBit()
{
    if(m_pFile && m_pFile->is_open())
        m_pFile->close();
    delete m_pFile;
    m_pFile = NULL;
}

bool CTwoBit::storeData()
{
    try
    {
        for(map<string,chrom_data>::iterator iter = m_data.begin(); iter != m_data.end(); iter++)
        {
            m_pFile->seekg(iter->second.offset,ios::beg);
            m_pFile->read((char*)&(iter->second.dnaSize),4);
            m_pFile->read((char*)&(iter->second.nBlockCount),4);
            uint32_t dummy=0;
            iter->second.nBlockStarts.reserve(iter->second.nBlockCount);
            for(uint32_t i = 0 ; i < iter->second.nBlockCount; i++)
            {
                m_pFile->read((char*)&dummy,4);
                iter->second.nBlockStarts.push_back(dummy);
            }
             iter->second.nBlockSizes.reserve(iter->second.nBlockCount);
            for(uint32_t i = 0 ; i < iter->second.nBlockCount; i++)
            {
                m_pFile->read((char*)&dummy,4);
                iter->second.nBlockSizes.push_back(dummy);
            }
            m_pFile->read((char*)&(iter->second.maskBlockCount),4);
            iter->second.maskBlockStarts.reserve(iter->second.maskBlockCount);
            for(uint32_t i = 0 ; i < iter->second.maskBlockCount; i++)
            {
                m_pFile->read((char*)&dummy,4);
                iter->second.maskBlockStarts.push_back(dummy);
            }
            iter->second.maskBlockSizes.reserve(iter->second.maskBlockCount);
            for(uint32_t i = 0 ; i < iter->second.maskBlockCount; i++)
            {
                m_pFile->read((char*)&dummy,4);
                iter->second.maskBlockSizes.push_back(dummy);
            }
            m_pFile->read((char*)&(iter->second.reserved),4);
            if(iter->second.dnaSize % 4 == 0)
            {
                iter->second.packedDNA = (char*)malloc(iter->second.dnaSize / 4);
                m_pFile->read(iter->second.packedDNA,iter->second.dnaSize / 4);
            }
            else
            {
                iter->second.packedDNA = (char*)malloc((iter->second.dnaSize/4)+1);
                m_pFile->read(iter->second.packedDNA,(iter->second.dnaSize/4)+1);
            }
            
        }
        return true;
    }
    catch(exception ex)
    {
        cerr << ex.what() << endl;
    }
    catch(...)
    {
        cerr << "unknown error loading 2bit data" << endl;
    }
    return false;
}

string CTwoBit::getFaSequence(const string& chrom,uint32_t start, uint32_t end,const char* faHeader ,uint32_t lineBreaks)
{
    string strTmp(getFaSequence(chrom,start,end,lineBreaks));
    ostringstream strReturn("");
    strReturn << faHeader << endl << strTmp;
    return strReturn.str();
}

string CTwoBit::getFaSequence(const string& chrom,uint32_t start, uint32_t end,uint32_t lineBreaks)
{
    if(lineBreaks==0)
        return getSequence(chrom,start,end);
    string strTmp(getSequence(chrom,start,end));
    ostringstream strReturn("");
    uint32_t counter = 0;
    for(string::const_iterator iter = strTmp.begin(); iter != strTmp.end(); iter++)
    {
        strReturn << *iter;
        if(++counter % lineBreaks == 0)
            strReturn << endl;
    }
    return strReturn.str();
}

string CTwoBit::getSequence(const string& chrom)
{
    map<string,chrom_data>::const_iterator iter = m_data.find(chrom);
    if(iter == m_data.end())
        return "";
    return getSequence(chrom,1,iter->second.dnaSize);
}

string CTwoBit::getSequence(const string& chrom,uint32_t start, uint32_t end, uint16_t flanks)
{
    map<string,chrom_data>::const_iterator iter = m_data.find(chrom);
    if(iter == m_data.end())
        return "";

    uint32_t start_neu  = 1;
    uint32_t end_neu    = end+flanks;
    if(start > flanks)
        start_neu = start-flanks;

    if(end_neu > iter->second.dnaSize)
        end_neu = iter->second.dnaSize;
    return getSequence(chrom,start_neu,end_neu);
}

string CTwoBit::getSequence(const string& chrom,uint32_t start, uint32_t end)
{
    ostringstream strRetrun("");
    uint32_t base_count = 0;
    map<string,chrom_data>::const_iterator iter = m_data.find(chrom);
    if(iter != m_data.end())
    {
        if(start < 1 || end < start || end > iter->second.dnaSize)
            return "";
        int8_t  startBit = 3-( (start-1) % 4);
        uint32_t startByte= (start-1) / 4;
        int8_t  endBit   = 3-(end % 4);
        uint32_t endByte  = end / 4;
        for(uint32_t idxByte = startByte; idxByte <= endByte; idxByte++ )
        {
            int8_t fromBit =3;
            int8_t toBit =0;
            uint8_t aktByte = iter->second.packedDNA[idxByte];
            if(idxByte==startByte)
                fromBit=startBit;
            if(idxByte == endByte)
                toBit=endBit+1;
            for(int8_t idxBit=fromBit; idxBit >= toBit; idxBit-- )
                strRetrun << getAllele(aktByte,idxBit);
        }
        nBlocks(strRetrun,iter,start,end);
    }
    return strRetrun.str();
}

void CTwoBit::nBlocks(ostringstream& strSequence,map<string,chrom_data>::const_iterator iter,uint32_t start,uint32_t end)
{
    vector<uint32_t>::const_iterator iterStarts = iter->second.nBlockStarts.begin();
    vector<uint32_t>::const_iterator iterSizes = iter->second.nBlockSizes.begin();
    for(;iterStarts != iter->second.nBlockStarts.end(); iterStarts++, iterSizes++)
    {
        uint32_t actStart = *iterStarts;
        uint32_t actSize  = *iterSizes;
        uint32_t actEnd  = (actStart+actSize)-1;
        if( (actStart >= start && actStart <= end) ||
            (actEnd >= start && actEnd <= end) ||
            (actStart < start && actEnd > end) )
        {
            // cout << iter->first << ':' << start << '-' << end << endl << actStart << " length " << actSize << " => " << actStart << '-' << actEnd << endl;
            if(actStart < start)
            {
                ++actSize-=(start-actStart);
                strSequence.seekp(ios::beg);
            }
            else
                strSequence.seekp(++actStart-start,ios::beg);
            if((uint32_t)strSequence.tellp()+actSize > strSequence.str().length())
                actSize=(strSequence.str().length()-strSequence.tellp());
            while(actSize-- > 0)
                strSequence.put('N');
        }
    }
}

char CTwoBit::getAllele(uint8_t theByte,uint8_t theBits)
{
    // 00   0  T
    // 01   1  C
    // 10   2  A
    // 11   3  G
    switch ( (theByte >> (theBits*2)) & 3 )
    {
        case 0:
            return 'T';
        case 1:
            return 'C';
        case 2:
            return 'A';
        case 3:
            return 'G';
    }
    return ' ';
}


