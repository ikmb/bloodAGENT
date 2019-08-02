/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CFastqCreator.cpp
 * Author: mwittig
 * 
 * Created on August 2, 2019, 12:31 PM
 */
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>

#include <regex>
#include <iterator>

#include "mytools.h"
#include "CFastqCreator.h"

using namespace std;

CFastqCreator::CFastqCreator(const std::string& filename) 
{
    m_fa = CMultiFasta(filename);
}

CFastqCreator::CFastqCreator(const CFastqCreator& orig) 
{
    
}

CFastqCreator::~CFastqCreator() {
}

void CFastqCreator::makeIlluminaPairedEnd(int minInsert, int maxInsert, int readLength, int coverage, std::string basename)
{
    
    ofstream R1(string(basename).append("_R1.fastq").c_str());
    ofstream R2(string(basename).append("_R2.fastq").c_str());
    
    
    if(maxInsert < minInsert)
        swap(minInsert,maxInsert);
    
    if(minInsert < readLength)
    {
        cerr << "Min insert size must be >= readLength. Exit." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(!R1 || !R2)
    {
        cerr << "could not create output files" << endl;
        exit(EXIT_FAILURE);
    }
    
    int ref_counter = 1;
    string qualities = string(readLength,'J');
    for(CMultiFasta::const_iterator i = m_fa.begin(); i != m_fa.end(); i++)
    {
        string name = i->first;
        string seq = i->second;

        int tilings = ((seq.size()/readLength)+1)*coverage;
        int inner_counter=0;
        
        while(tilings-- != 0)
        {
            int insert_size = minInsert + rand() % (( maxInsert + 1 ) - minInsert);
            int ref_start = rand() % ( seq.length()-insert_size + 1 );
            
            string fragment = seq.substr(ref_start,insert_size);
            
            
           R1 << "@HWI-ST778:164:D2C3NACXX:7:1101:"  << std::setfill('0') << std::setw(4) << ref_counter << ':' << std::setfill('0') << std::setw(4) << ++inner_counter << " 1:N:0:TCTTCACA" << endl;
           R1 << fragment.substr(0,readLength) << endl << '+' << endl;
           R1 << qualities << endl;
           R2 << "@HWI-ST778:164:D2C3NACXX:7:1101:"  << std::setfill('0') << std::setw(4) << ref_counter << ':' << std::setfill('0') << std::setw(4) << inner_counter << " 2:N:0:TCTTCACA" << endl;
           R2 << fragment.substr(fragment.length()-readLength) << endl << '+' << endl;
           R2 << qualities << endl;
        }
        ref_counter++;
        
    }
    R1.close();
    R2.close();
}
