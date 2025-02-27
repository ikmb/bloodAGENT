/* 
 * File:   CRefGeneEntry.h
 * Author: mwittig
 *
 * Created on May 22, 2017, 4:05 PM
 */

#ifndef CREFGENEENTRY_H
#define	CREFGENEENTRY_H



class CRefGeneEntry {
  
    public:  
    
    class mutation{
    public:
        mutation(){ref_allele ='V';alt_allele='N';ref_amino=std::pair<std::string,char>("Xaa",'X');alt_amino=std::pair<std::string,char>("Xaa",'X');cDNA_position=0;amino_position=0;}
        char ref_allele;
        char alt_allele;
        long cDNA_position;
        std::pair<std::string,char> ref_amino;
        std::pair<std::string,char> alt_amino;
        long amino_position;
        
        
        std::string aaChangeAsString()const{return ref_amino.first+std::to_string(amino_position)+alt_amino.first;}
        
    };
    
   CRefGeneEntry(vector<std::string>&);
    CRefGeneEntry(const CRefGeneEntry& orig);
    virtual ~CRefGeneEntry();
    
    
    std::string name()const{return m_name;}
    std::string name2()const{return m_name2;}
    std::string chrom()const{return m_chrom;}
    char strand()const{return m_strand;}
    long  txStart()const{return m_txStart;}
    long  txEnd()const{return m_txEnd;}
    
    bool isForwardStrand()const{return m_strand == '+';}
    bool isReverseStrand()const{return m_strand == '-';}
    
    
    /// get the genomic coordinate if a given cDNA position
    /// cDNA positions are 0 based!!!
    long getGenomicCoordinate(long idx)const;
    /// get the cDNA coordinate if a given genomic position
    /// genomic positions are 1 based!!!
    /// @return (pair<long,long>): first=CDS position; second=offset e.g. 337+5 if it is the fifth intronic NT of the intron after CDS-pos 337
    std::pair<long,long> getCDNAcoordinate(long pos)const;
    /// Annotate an mutation
    /// InDels not supported currently
    /// @param (int): cDNA coordinate of single point mutation
    /// @param (char): the alternative nucleotide
    /// @param (CTwoBit): the two bit reference (must be the same build the coordinates are taken from!)
    CRefGeneEntry::mutation getMutation(long idx,char alt, CTwoBit& tb);
    
    /// Return the sector of the genomic position
    /// Where in the transcript is this position? 5', exon, intron, 3'
    /// wrong chromosome returns NA
    /// @param (string): chromosome
    /// @param (long): position
    /// @return (string): NA, 5', exon [1..N], intron [1..N], 3'
    std::string getSector(const string& chrom, const long& position);
    
    /// Return the sector of a transcript
    /// @return (string): list with transcript compartment chrom start end
    std::string getSectors();
    
private:

    size_t          m_bin;
    string          m_name;
    string          m_chrom;
    char            m_strand;
    long            m_txStart;
    long            m_txEnd;
    long            m_cdsStart;
    long            m_cdsEnd;
    long            m_exonCount;
    vector<long>    m_exonStarts;
    vector<long>    m_exonEnds;
    long            m_score;
    string          m_name2;
    string          m_cdsStartStat;
    string          m_cdsEndStat;
    vector<int>     m_exonFrames;

};

#endif	/* CREFGENEENTRY_H */

