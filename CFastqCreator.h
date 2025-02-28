/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CFastqCreator.h
 * Author: mwittig
 *
 * Created on August 2, 2019, 12:31 PM
 */

#ifndef CFASTQCREATOR_H
#define CFASTQCREATOR_H

#include "meinetools.h"

class CFastqCreator {
public:
    CFastqCreator(const std::string& filename);
    CFastqCreator(const CFastqCreator& orig);
    virtual ~CFastqCreator();
    
    /**
    * generate paired end fast q files
    * @param minInsert min length of the paired end fragment
    * @param maxInsert max  length of the paired end fragment
    * @param readLength the leangth of each read
    * @param coverage we aim for this coverage
    * @param basename the basename extended by .R1.fastq and .R2.fastq
    */
    void makeIlluminaPairedEnd(int minInsert, int maxInsert, int readLength, int coverage, std::string basename);
    
    /**
    * generate PacBio ccs reads in sam format
    * @param minSize min length of the read
    * @param maxSize max  length of the read
    * @param coverage we aim for this coverage
    * @param basename output filename (name it *.sam to avois confusion)
    */
    void makePacBioRead(int minSize, int maxSize, int coverage, std::string filename);
    
private:

    CMultiFasta m_fa;
};

#endif /* CFASTQCREATOR_H */

