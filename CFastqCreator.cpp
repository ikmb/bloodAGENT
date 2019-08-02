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

void CFastqCreator::makePacBioRead(int minSize, int maxSize, int coverage, std::string filename)
{
    ofstream R1(filename.c_str());
    
    if(!R1 )
    {
        cerr << "could not create output files" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Write header
    R1 << "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1" << endl;
    R1 << "@RG\tID:a92d5a0b\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=101-500-400;SEQUENCINGKIT=101-427-800;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000;BarcodeFile=/opt/pacbio/smrtlink/userdata/jobs_root/000/000065/pacbio-barcodes/c51caf4d_0a32_4cbd_a508_3e92ccf4589b_index_pacbio_univ/barcodeset.xml;BarcodeHash=fef9cfb8d4d08f346f0dbfec2e5e3ae4;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\tPU:m54349_190510_125704\tPM:SEQUEL" << endl;
    R1 << "@PG\tID:ccs-3.4.1\tPN:ccs\tVN:3.4.1\tDS:Generate circular consensus sequences (ccs) from subreads.\tCL:ccs ccs --minLength=1500 --maxLength=18000 --numThreads=15 --force --minPasses=5 /inPath/bc1001.bam in_silico_generated.bam" << endl;

    if(minSize < maxSize)
        swap(minSize,maxSize);
    
    int counter = 1;
    for(CMultiFasta::const_iterator i = m_fa.begin(); i != m_fa.end(); i++)
    {
        string name = i->first;
        string seq = i->second;

        int tilings = ((seq.size()/ (minSize+(maxSize-minSize)/2) )+1)*coverage;
        while(tilings-- != 0)
        {
            int fragment_size = minSize + rand() % (( maxSize + 1 ) - minSize);
            int ref_start = rand() % ( seq.length()-fragment_size + 1 );
            
            string fragment = seq.substr(ref_start,fragment_size);
            
           R1 << "m54349_190510_125704/1"   << std::setfill('0') << std::setw(7)  << counter << "/ccs\t4\t*\t0\t255\t*\t*\t0\t0\t";
           R1 << seq.substr(ref_start,fragment_size) << '\t' << string(fragment_size,'~') << '\t';
           R1 << "RG:Z:a92d5a0b\tbc:B:S,0,0\tbq:i:85\tnp:i:18\trq:f:0.999881\tsn:B:f,4.86057,9.51716,4.88941,8.44968\tt1:f:1.109\tt2:f:0.041\tt3:f:6.429\tzm:i:1" << std::setfill('0') << std::setw(7)  << counter << endl;
        }
        counter++;
        
    }
    R1.close();
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
    string qualities = string(readLength,'~');
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
