/* 
 * File:   CRefGeneTable.cpp
 * Author: mwittig
 * 
 * Created on May 22, 2017, 3:41 PM
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <set>


#include "CMyException.h"
#include "MyTools.h"
#include "CTwoBit.h"
#include "CRefGeneEntry.h"
#include "CRefGeneTable.h"

CRefGeneTable::CRefGeneTable(const char* filename) 
              : m_header_line("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames")
{
    loadData(filename);
}

CRefGeneTable::CRefGeneTable(const CRefGeneTable& orig) : m_header_line(orig.m_header_line)
{
    m_entries = orig.m_entries;
}

CRefGeneTable::~CRefGeneTable() {
}

bool CRefGeneTable::loadData(const char* filename)
{
    string strLine;
    ifstream* pFile = new ifstream(filename);
    if(!pFile || !pFile->is_open())
        throw CMyException(string("Error reading file in CRefGeneTable!")+"failed to open 2bit file \'" + filename + "\'");
    getline(*pFile,strLine);
    if(strLine.compare(m_header_line)!=0)
        throw CMyException(string("Error reading file in CRefGeneTable! Invalid header. Got:\"")+strLine+ "\", expected: \""+m_header_line+'\"');
    size_t expected_cell_count = CMyTools::GetParsedLine(strLine).size();
    int line_count = 1;
    while(getline(*pFile,strLine))
    {
        if(strLine.size() == 0) // skip empty lines
            continue;
        line_count++;
        vector<string> parsed = CMyTools::GetParsedLine(strLine);
        if(parsed.size() != expected_cell_count)
            throw CMyException(string("invalid number of cells in line ")+to_string(line_count)+ " of file \""+filename+'\"');
        m_entries.push_back(CRefGeneEntry(parsed));
        //cout << m_entries[m_entries.size()-1].name() << endl;
    }
    pFile->close();
    delete pFile;
    return true;
}

const bool CRefGeneTable::contains(const std::string& id)const
{
    try{
        this->operator [](id);
        return true;
    }
    catch(CMyException& err)
    {
        return false;
    }
}

std::string CRefGeneTable::transcriptAt(const std::string& chrom, long pos, long pos2, long flanks)
{
    set<string> setRet;
    for(vector<CRefGeneEntry>::const_iterator i = m_entries.begin(); i != m_entries.end(); i++)
    {
        long txStart = i->txStart()-flanks;
        long txEnd = i->txEnd()+flanks;
        
        if(     i->chrom().compare(chrom) == 0 && 
            (
                (txStart <= pos && txEnd >= pos) ||
                (txStart <= pos2 && txEnd >= pos2) ||
                (pos < txStart && pos2 > txEnd)
            )
          )
          setRet.insert(i->name2());
    }
    ostringstream strRet("");
    std::copy(setRet.begin(), setRet.end(), std::ostream_iterator<std::string>(strRet, ","));
    return strRet.str().substr(0,strRet.str().size()-1);
}

const CRefGeneEntry& CRefGeneTable::operator[](const std::string& id)const
{
    for(vector<CRefGeneEntry>::const_iterator i = m_entries.begin(); i != m_entries.end(); i++)
    {
        if(i->name().compare(id) == 0 || i->name2().compare(id) == 0)
            return *i;
    }
    throw CMyException("Id not found");
}



