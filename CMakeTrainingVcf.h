/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CMakeTrainingVcf.h
 * Author: mwittig
 *
 * Created on July 30, 2019, 12:42 PM
 */

#ifndef CMAKETRAININGVCF_H
#define CMAKETRAININGVCF_H

class CMakeTrainingVcf {
public:
    CMakeTrainingVcf();
    CMakeTrainingVcf(const CMakeTrainingVcf& orig);
    virtual ~CMakeTrainingVcf();
    
    
    static std::string getHomEntries(const std::string& system, const CIsbtPtAllele& allele, const CISBTAnno& anno);
    static std::string getHetEntries(const std::string& system, const CIsbtPtAllele& alleleA, const CIsbtPtAllele& alleleB, CISBTAnno& anno, bool phased, int dropout_prob, int haplotype_crack);
    static int         getRandomInteger(const int start, const int end);
    static void        removeEmptyStringsFromSet(std::set<std::string>& variations);
    
private:

};

#endif /* CMAKETRAININGVCF_H */

