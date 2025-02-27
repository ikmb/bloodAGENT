/* 
 * File:   CTwoBit.h
 * Author: mwittig
 *
 * Created on February 18, 2011, 10:52 AM
 */

#ifndef _CTWOBIT_H
#define	_CTWOBIT_H

using namespace std;


class CTwoBit {

    typedef struct{
        uint32_t signature;
        uint32_t version;
        uint32_t seqCount;
        uint32_t reserved;
    }header;

    typedef struct{
        uint32_t            offset;
        uint32_t            dnaSize;
        uint32_t            nBlockCount;
        vector<uint32_t>    nBlockStarts;
        vector<uint32_t>    nBlockSizes;
        uint32_t            maskBlockCount;
        vector<uint32_t>    maskBlockStarts;
        vector<uint32_t>    maskBlockSizes;
        uint32_t            reserved;
        char*               packedDNA;
    }chrom_data;


public:
    CTwoBit();
    CTwoBit(const char*);
    CTwoBit(const CTwoBit& orig);
    virtual ~CTwoBit();
    
    /// coordinates are 1-based index
    string getSequence(const std::string& chrom,uint32_t start, uint32_t end, uint16_t flanks);
    /// coordinates are 1-based index
    string getSequence(const std::string& chrom,uint32_t start, uint32_t end);
    /// coordinates are 1-based index
    string getSequence(const std::string& chrom);
    /// coordinates are 1-based index
    string getFaSequence(const std::string& chrom,uint32_t start, uint32_t end,uint32_t lineBreaks = 50);
    /// coordinates are 1-based index
    string getFaSequence(const std::string& chrom,uint32_t start, uint32_t end,const char* faHeader ,uint32_t lineBreaks = 50);

    vector<pair<string,unsigned int> >  getChromSizes();

private:
    bool        m_bOpen;
    ifstream*   m_pFile;
    header      m_stcHeader;

    map<string,chrom_data> m_data;

    bool storeOffsets();
    bool storeData();
    char getAllele(uint8_t,uint8_t);
    void nBlocks(ostringstream& strSequence,map<string,chrom_data>::const_iterator iter,uint32_t start,uint32_t end);
    
    bool loadTwoBit(const char*);
    void closeTwoBit();

    void initialize(){
        m_bOpen = false;
        m_pFile = NULL;
    }
    
    chrom_data getChromDataStruc(){
        chrom_data cdReturn;
        cdReturn.dnaSize = 0;
        cdReturn.maskBlockCount = 0;
        cdReturn.nBlockCount = 0;
        cdReturn.offset = 0;
        cdReturn.packedDNA = NULL;
        cdReturn.reserved = 0;
        return cdReturn;
    }

};

#endif	/* _CTWOBIT_H */

