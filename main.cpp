/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: mwittig
 *
 * Created on July 18, 2019, 12:29 PM
 */

#include <cstdlib>
#include <vector>
#include <string>
#include <libgen.h>

#include <regex>
#include <iterator>

#include "api/BamIndex.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
//#include "shared/bamtools_global.h"


// command line parsing
//#include <tclap/CmdLine.h>
#include <complex>

#include "mytools.h"
#include "CBigWigReader.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CVariantChain.h"
#include "CVariantChains.h"
#include "CTranscript.h"
#include "CTranscriptAnno.h"

#include "CIsbtGt2Pt.h"
#include "CMakeTrainingVcf.h"
#include "CFastqCreator.h"

using namespace std;
using namespace BamTools;


/*
 * 
 */
int main(int argc, char** argv) 
{
    //cout << argc << endl;
    
    //CFastqCreator fq("/home/mwittig/data/Genotypisierung/Haemocarta/Erythrogene_Tables/IndelDups/target.hg38.purplevariation.fasta");
    //fq.makePacBioRead(1500,9000,50,"/home/mwittig/ramDisk/purple.ccs.5passes.sam");
    //return 0;
    
    /*
    CVcf vcf("/home/mwittig/data/Genotypisierung/Haemocarta/Erythrogene_Tables/IndelDups/purple.hg19.bwa.variants.ISBT.vcf.gz");
    while(vcf.read_record())
    {
        CVcfSnp act_snp = vcf.get_record();
        cerr << act_snp << endl;
    };
    return 0;
    */
    /* ********************
    // ON MWMOB
    CTranscriptAnno trans_anno("/home/mwittig/coding/cpp/deepBlood/data/config/exonic_annotation.hg19.abotarget.txt");
    CBigWigReader bwr("/home/mwittig/coding/cpp/deepBlood/data/example/bc1002.asm20.hg19.ccs.5passes.abotarget.bw");
    CISBTAnno  isbt("/home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation.dat");
    CIsbtGt2Pt isbTyper("/home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation.dat");
    CVcf vcf_file("/home/mwittig/coding/cpp/deepBlood/data/example/bc1002.asm20.hg19.ccs.5passes.phased.phenotype.SNPs.vcf.gz");
     */
    //* ********************
    // ON Cluster
    if(argc != 6)
    {
        cerr << "Please provide 5 parameters in the right order:" << endl
             << basename(argv[0]) << " variation_annotation.dat \\" << endl
             << "   genotype_to_phenotype_annotation.dat\\" << endl
             << "   exonic_annotation.hg19.abotarget.txt\\" << endl
             << "   sample.bw\\" << endl
             << "   sample.vcf.gz\\" << endl;
                
        exit(EXIT_FAILURE);
    }
    CISBTAnno  isbt(argv[1]);
    CIsbtGt2Pt isbTyper(argv[2]);
    CTranscriptAnno trans_anno(argv[3]);
    CBigWigReader bwr(argv[4]);
    CVcf vcf_file(argv[5]);
    //*/
    isbt.addCoverage(bwr);
    std::set<string> loci = isbt.loci();
    
    // generate VCF Test Data
    
    /*// make test data
    //    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mwittig/coding/cpp/MyTools/dist/Debug/GNU-Linux/:/home/mwittig/coding/cpp/bamtools/build/src/api/:/home/mwittig/coding/cpp/htslib
    for(auto locus:loci)
    {
        for(auto known_allele:isbTyper.alleleVector(locus))
        {
             CMyTools::cmd("cat /home/mwittig/coding/cpp/deepBlood/data/testdata/VCFheaderTestData.vcf > /home/mwittig/ramDisk/simulator.vcf",false);
            ofstream outfile;
            outfile.open("/home/mwittig/ramDisk/simulator.vcf", std::ios_base::app);
            outfile << CMakeTrainingVcf::getHomEntries(locus,known_allele,isbt) << endl;
            outfile.close();
            /// if you analyze simulated VCF do the following on your CVariantChains object
            // vcs.removeReferenceSnps();
        }
    }//*/
    
    CVariantChains vcs(&isbt);
    //CVcf vcf_file(argv[3]);
    while(vcf_file.read_record())
    {
        CVcfSnp act_snp = vcf_file.get_record();
        //cerr << act_snp << endl;
        bool adding_successfull = vcs.add(act_snp);
        if(!adding_successfull)
            cerr << act_snp << "\tadding SNP failed" << endl;
    };
    
    for(auto locus:loci)
    {
        isbTyper.type(locus,vcs);
        cout << locus << '\t' << isbTyper.getCallAsString(locus,false) << endl;
    }

    
    // parse results
    // awk 'BEGIN{FS="\t"}{if(NR%2==0)print $2;else printf "%s\t%s\t",$1,$2}' simulator.txt > simulation.csv
    return 0;
    
   
    
    
    
    return 0;
}

