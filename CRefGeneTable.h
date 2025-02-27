/* 
 * File:   CRefGeneTable.h
 * Author: mwittig
 *
 * Created on May 22, 2017, 3:41 PM
 */

#ifndef CREFGENETABLE_H
#define	CREFGENETABLE_H

#include "CRefGeneEntry.h"


using namespace std;

class CRefGeneTable {
public:
    CRefGeneTable(const char *);
    CRefGeneTable(const CRefGeneTable& orig);
    virtual ~CRefGeneTable();
    
    const std::vector<CRefGeneEntry>& entries()const{return m_entries;}
    const CRefGeneEntry& operator[](const std::string& id)const;
    const bool contains(const std::string& id)const;
    
    std::string transcriptAt(const std::string& chrom, long pos, long pos2 = -1,long flanks = 0);
    
private:
    
    
    const string m_header_line;
    std::vector<CRefGeneEntry> m_entries;
    
    bool loadData(const char*);
    
    
};

#endif	/* CREFGENETABLE_H */

