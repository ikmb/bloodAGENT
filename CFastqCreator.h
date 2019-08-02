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

#include "mytools.h"

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
    * @param basename the basename extended by _R1.fastq and _R2.fastq
    */
    void makeIlluminaPairedEnd(int minInsert, int maxInsert, int readLength, int coverage, std::string basename);
    
private:

    CMultiFasta m_fa;
};

#endif /* CFASTQCREATOR_H */

