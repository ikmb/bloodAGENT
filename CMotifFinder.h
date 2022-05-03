/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   CMotifFinder.h
 * Author: mwittig
 *
 * Created on January 5, 2022, 10:34 AM
 */

#ifndef CMOTIFFINDER_H
#define CMOTIFFINDER_H

class CMotifFinder {
public:
    CMotifFinder(const std::string& config, const std::string& filenames,int verbose = 0);
    CMotifFinder(const CMotifFinder& orig);
    virtual ~CMotifFinder();
    
    friend std::ostream& operator<<(std::ostream& os, const CMotifFinder& me);
    
    static std::map<string,int>  findMotifs(const std::string& filename,std::vector<string>& motifs);
    static void  findMotifs(const std::string& filename,std::map<string,int>& motifs);
    
    map<std::string,CVcfSnp> getSystemsMotifSnps(const std::string& system)const;
    std::set<std::string> getSystems()const;
    
private:

    void storeMotifSnps(CParsedTextfile& config, map<string,int> motifs);
    // first string is system
    // second string is SNP name in short (see variation annotation)
    map<std::string,map<std::string,CVcfSnp> >    m_motifs_snps;
    
};

#endif /* CMOTIFFINDER_H */

