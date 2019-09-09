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
#include <tclap/CmdLine.h>
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


void phenotype(const string& arg_target_anno,const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_vcf_file,const string& arg_bigWig,int arg_coverage, int arg_verbose, float arg_top_hits = 1.0, const string& arg_locus = "");

/*
 * 
 */
int main(int argc, char** argv) 
{
    try
    {
        ostringstream tclapMessage;
        tclapMessage << basename(argv[0]) << ':'  << endl 
                     << "the NGS ISBT typer of the Institute of Clinical Molecular Biology in Kiel, Germany"  << endl
                     << "" << endl 
                     << "";
        TCLAP::CmdLine cmd   (tclapMessage.str().c_str(), ' ', "0.9");
        TCLAP::CmdLine cmdjob(tclapMessage.str().c_str(), ' ', "0.9");

        tclapMessage.clear();
        tclapMessage << "define which kind of job should be performed:"  << endl
                     << "  phenotype: report phenotypes of all systems" << endl 
                     << "  pht: report phenotypes of one system";
        TCLAP::ValueArg<string> tc_jobType     ("j","job","define which kind of job should be performed",true,"","string");
        TCLAP::ValueArg<string> tc_jobTypeInner("j","job","define which kind of job should be performed",true,"","string");
        cmd.add(tc_jobType);
        cmdjob.add(tc_jobTypeInner);

        cmd.ignoreUnmatched(true);
        cmd.parse(argc,argv);
        cerr << "job type " << tc_jobType.getValue() << ". looking for required parameters ..." << endl;
        TCLAP::ValueArg<int> tc_verbose("d","verbose","Set verbose level. 0-1. 0 == no verbose, 1 == full verbose.",false,0,"int");
        cmdjob.add(tc_verbose);
        if(tc_jobType.getValue().compare("phenotype") == 0 || tc_jobType.getValue().compare("pht") == 0)
        {
            /*
             "${OUTPUT_PATH}" -j anchor -a "/home/mwittig/data/Genotypisierung/Haemocarta/NGS/Paralogs/anchors.bed" -b "/home/mwittig/data/Genotypisierung/Haemocarta/NGS/Paralogs/homolgous_parts.bed" -i "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.bam" -p "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.meets.bam" -f "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.fails.bam"
             */
            TCLAP::ValueArg<string> tc_abo_target_annotation("t","target","A text file containing the transcript annotation of all blood group targets (UCSC table export). This file comes with the package an can be usually found in the subfolder data.",true,"data/config/exonic_annotation.hg19.abotarget.txt","string");
            cmdjob.add(tc_abo_target_annotation);
            TCLAP::ValueArg<string> tc_variants("s","variants","A text file containing the variant annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/genotype_to_phenotype_annotation.dat","string");
            cmdjob.add(tc_variants);
            TCLAP::ValueArg<string> tc_gt2pt("g","gt2pt","A text file containing the genotype to phenotype annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/genotype_to_phenotype_annotation.dat","string");
            cmdjob.add(tc_gt2pt);
            TCLAP::ValueArg<string> tc_vcf("v","vcf","A vcf file with the variants of the sample. Please be aware of different genome build. This file should fit to the config files from parameters target, variants and gt2pt.",true,"","string");
            cmdjob.add(tc_vcf);
            TCLAP::ValueArg<string> tc_bigwig("b","bigwig","The output bam file to which all alignments go that failed the filtering criteria.",true,"","string");
            cmdjob.add(tc_bigwig);
            TCLAP::ValueArg<int> tc_coverage("c","coverage","The minimum required coverage for a solid call.",false,10,"int");
            cmdjob.add(tc_coverage);
            TCLAP::ValueArg<float> tc_tophits("r","scoreRange","Value between 0.0 and 1.0. After genotyping multiply the best score with this value and report all genotypes better than the resulting value.",false,1.0,"int");
            cmdjob.add(tc_tophits);

            
            // ln -s ../mnts/sftp\:host\=ikmbhead.rz.uni-kiel.de\,user\=sukko545/ifs/data/nfs_share/sukko545/haemo/DZHK/190233/190233.hg19.bwa.variants.ISBT.vcf.gz SNPs.vcf.gz
            // ln -s ../mnts/sftp\:host\=ikmbhead.rz.uni-kiel.de\,user\=sukko545/ifs/data/nfs_share/sukko545/haemo/DZHK/190233/190233.hg19.bwa.bw coverage.bw

            if(tc_jobType.getValue().compare("phenotype") == 0)
            {
                cmdjob.parse(argc,argv);
                cerr << "Everything found. starting run ..." << endl;
                phenotype(tc_abo_target_annotation.getValue(),
                        tc_variants.getValue(),
                        tc_gt2pt.getValue(),
                        tc_vcf.getValue(),
                        tc_bigwig.getValue(),
                        tc_coverage.getValue(), 
                        tc_verbose.getValue(),
                        tc_tophits.getValue());
                exit(EXIT_SUCCESS);
            }
            else if(tc_jobType.getValue().compare("pht") == 0)
            {
                TCLAP::ValueArg<string> tc_system("t","target","The system that should be analyzed",true,"ABO","string");
                cmdjob.add(tc_system);
                cmdjob.parse(argc,argv);
                cerr << "Everything found. starting run ..." << endl;
                phenotype(tc_abo_target_annotation.getValue(),
                        tc_variants.getValue(),
                        tc_gt2pt.getValue(),
                        tc_vcf.getValue(),
                        tc_bigwig.getValue(),
                        tc_coverage.getValue(), 
                        tc_verbose.getValue(),
                        tc_tophits.getValue(),
                        tc_system.getValue());
                exit(EXIT_SUCCESS);
            }
        }
        else
        {
            cerr << "job " << tc_jobType.getValue() << " not implemented." << endl;
            exit(EXIT_SUCCESS);
        }
    }
    catch(const TCLAP::ArgException& err)
    {
        cerr << "Error: " << err.what() << endl;
        exit(EXIT_FAILURE);
    }
    catch(const CMyException& err)
    {
        cerr << "Error parsing command line: " << err.what() << endl;
        exit(EXIT_FAILURE);
    }
    catch(...)
    {
        cerr << "Unexpected error." << endl;
        exit(EXIT_FAILURE);
    }
    
    
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
    //*/
    //* ********************
    // ON MWMOB
    return 0;
}
        

void phenotype(const string& arg_target_anno,const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_vcf_file,const string& arg_bigWig,int arg_coverage, int arg_verbose, float arg_top_hits, const string& arg_locus)
{
    try
    {
        CTranscriptAnno trans_anno(arg_target_anno);
        if(arg_verbose >= 2)
            cerr << "transcript annotation loaded from:"  << arg_target_anno << endl;
        CISBTAnno  isbt(arg_isbt_SNPs);
        if(arg_verbose >= 2)
            cerr << "ISBT variations loaded from:"  << arg_isbt_SNPs << endl;
        CIsbtGt2Pt isbTyper(arg_genotype_to_phenotype);
        if(arg_verbose >= 2)
            cerr << "ISBT genotype to phenotype translation loaded from:"  << arg_genotype_to_phenotype << endl;
        CVcf vcf_file(arg_vcf_file);
        if(arg_verbose >= 2)
            cerr << "VCF file loaded from:"  << arg_vcf_file << endl;
        CBigWigReader bwr(arg_bigWig);
        if(arg_verbose >= 2)
            cerr << "BigWig file loaded from:"  << arg_bigWig << endl;

        isbt.addCoverage(bwr);
        std::set<string> loci = isbt.loci();


        if(arg_verbose >= 2)
        {
            cerr << "ISBT annotation:" << endl << isbt << endl << endl;
            cerr << "ISBT genotype to phenotype translation:" << endl << isbTyper << endl << endl;

         }

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
            if(!adding_successfull && arg_verbose >= 1)
                cerr << act_snp << "\tSNP not added as it is not ISBT relevant" << endl;
            else if(adding_successfull && arg_verbose == 2)
                cerr << act_snp << "\tISBT relevant SNP added" << endl;
        };

        if(arg_verbose >= 1)
            cerr << "Variant chains of current sample: " << vcs << endl;    
        for(auto locus:loci)
        {
            if( arg_locus.length() == 0 ||  locus.compare(arg_locus) == 0 )
            {
                isbTyper.type(locus,vcs);
                cout << locus << '\t' << isbTyper.getCallAsString(locus,false,arg_top_hits) << endl;
            }
        }
        if(arg_verbose >= 2)
            cerr << isbTyper << endl;
    }
    catch(const CMyException& err)
    {
        throw(err);
    }
    catch(...)
    {
        throw(CMyException("Unexpected Error in phenotype function"));
    }
}

