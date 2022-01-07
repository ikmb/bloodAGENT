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
#include <experimental/filesystem>

#include "api/BamIndex.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
//#include "shared/bamtools_global.h"


// command line parsing
#include <tclap/CmdLine.h>
#include "json/single_include/nlohmann/json.hpp"
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

#include "CMotifFinder.h"

#include "CIsbtGt2Pt.h"
#include "CMakeTrainingVcf.h"
#include "CFastqCreator.h"

using namespace std;
using namespace BamTools;

#define APP_VERSION_DEEPBLOOD "v0.0.1"

void phenotype(const string& arg_target_anno,const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_vcf_file,const string& arg_bigWig,
        const string& arg_fastqgz, const string& arg_motifs,int arg_coverage, int arg_verbose, float arg_top_hits = 1.0, const string& arg_locus = "", 
        bool arg_is_in_silico = false, const string& sampleId = "", const string& outfile = "");
void inSilicoVCF(const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_allele_A,const string& arg_allele_B, int arg_verbose);
string getArgumentList(TCLAP::CmdLine& args);

/*
 * export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/mwittig/coding/cpp/MyTools/dist/Debug/GNU-Linux/:/home/mwittig/coding/fremd/htslib:/home/mwittig/coding/fremd/libBigWig
 * ln -s ../mnts/gvfs/sftp\:host\=medcluster.medfdm.uni-kiel.de\,user\=sukko545/work_ifs/sukko545/haemo/PacBio PacBio
 * 
 * StatusPlus
 * http://134.245.63.197:4200/
 * --locus ABO,RHD,LU,KEL,FY,JK,DI
 * test, test
 * 
  "${OUTPUT_PATH}" --job phenotype --target /home/mwittig/coding/cpp/deepBlood/data/config/exonic_annotation.hg19.abotarget.txt 
                   --variants /home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation.dat 
                   --gt2pt /home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation.dat --vcf /home/mwittig/ramDisk/SNPs.vcf 
                   --bigwig /home/mwittig/coding/cpp/deepBlood/data/config/fake_50x_coverage_at_target.bw 
                   --coverage 10 --verbose 1 --insilicovcf --locus RHD --scoreRange 0.0
 
  "${OUTPUT_PATH}"  --job vcf --variants /home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation.dat 
                    --gt2pt /home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation.dat -a "RHD*01EL.24" -b "RHD*01W.72"
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
        TCLAP::ValueArg<int> tc_verbose("d","verbose","Set verbose level. 0-2. 0 == no verbose, 1 == warnings, 2 == status, 3 == full",false,0,"int");
        cmdjob.add(tc_verbose);
        
        cmd.ignoreUnmatched(true);
        cmd.parse(argc,argv);
        if(tc_verbose.getValue() >= 2)
            cerr << "job type " << tc_jobType.getValue() << ". looking for required parameters ..." << endl;
         if(tc_jobType.getValue().compare("phenotype") == 0)
        {
            /*
             "${OUTPUT_PATH}" -j anchor -a "/home/mwittig/data/Genotypisierung/Haemocarta/NGS/Paralogs/anchors.bed" -b "/home/mwittig/data/Genotypisierung/Haemocarta/NGS/Paralogs/homolgous_parts.bed" -i "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.bam" -p "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.meets.bam" -f "/home/mwittig/data/tmp/Blood/G06322.hg38.PEonly.fails.bam"
             */
            TCLAP::ValueArg<string> tc_abo_target_annotation("t","target","A text file containing the transcript annotation of all blood group targets (UCSC table export). This file comes with the package an can be usually found in the subfolder data.",true,"data/config/exonic_annotation.hg19.abotarget.txt","string");
            cmdjob.add(tc_abo_target_annotation);
            TCLAP::ValueArg<string> tc_variants("s","variants","A text file containing the variant annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/variation_annotation.dat","string");
            cmdjob.add(tc_variants);
            TCLAP::ValueArg<string> tc_gt2pt("g","gt2pt","A text file containing the genotype to phenotype annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/genotype_to_phenotype_annotation.dat","string");
            cmdjob.add(tc_gt2pt);
            TCLAP::ValueArg<string> tc_vcf("v","vcf","A vcf file with the variants of the sample. Please be aware of different genome build. This file should fit to the config files from parameters target, variants and gt2pt.",true,"","string");
            cmdjob.add(tc_vcf);
            TCLAP::ValueArg<string> tc_bigwig("b","bigwig","The big wig or wig file that contains the coverage data..",true,"","string");
            cmdjob.add(tc_bigwig);
            TCLAP::ValueArg<string> tc_fastqgz("z","fastq","The fastq.gz files comma separated",true,"","string");
            cmdjob.add(tc_fastqgz);
            TCLAP::ValueArg<string> tc_motifs("m","motif","Configuration file that lists specific sequence motifs. These motifs identify SNPs that usually are not present in vcf files",true,"","string");
            cmdjob.add(tc_motifs);
            TCLAP::ValueArg<int> tc_coverage("c","coverage","The minimum required coverage for a solid call.",false,10,"int");
            cmdjob.add(tc_coverage);
            TCLAP::ValueArg<float> tc_tophits("r","scoreRange","Value between 0.0 and 1.0. After genotyping multiply the best score with this value and report all genotypes better than the resulting value.",false,1.0,"int");
            cmdjob.add(tc_tophits);
            TCLAP::SwitchArg tc_isInSilico("i","insilicovcf","If the input vcf is in silico generated VCF file, this trigger must be set! If not, the differences between LRG and GRCh are always analyzed resulting in wrong calls.",false);
            cmdjob.add(tc_isInSilico);
            TCLAP::ValueArg<string> tc_locus("l","locus","Get typing of specific loci only. Provide a comma separated list without spaces E.g. ABO,RHD, ...",false,"","string");
            cmdjob.add(tc_locus);
            TCLAP::ValueArg<string> tc_Id("f","id","provide a sample identifier that will be used for result output",false,"","string");
            cmdjob.add(tc_Id);
            TCLAP::ValueArg<string> tc_output("o","out","provide output file",false,"","string");
            cmdjob.add(tc_output);

            // ln -s ../mnts/sftp\:host\=ikmbhead.rz.uni-kiel.de\,user\=sukko545/ifs/data/nfs_share/sukko545/haemo/DZHK/190233/190233.hg19.bwa.variants.ISBT.vcf.gz SNPs.vcf.gz
            // ln -s ../mnts/sftp\:host\=ikmbhead.rz.uni-kiel.de\,user\=sukko545/ifs/data/nfs_share/sukko545/haemo/DZHK/190233/190233.hg19.bwa.bw coverage.bw

            cmdjob.parse(argc,argv);
            if(tc_verbose.getValue() >= 2)
            {
                cerr << "Parameter validation passed. Starting run ..." << endl;
            }
            
            phenotype(tc_abo_target_annotation.getValue(),
                    tc_variants.getValue(),
                    tc_gt2pt.getValue(),
                    tc_vcf.getValue(),
                    tc_bigwig.getValue(),
                    tc_fastqgz.getValue(),
                    tc_motifs.getValue(),
                    tc_coverage.getValue(), 
                    tc_verbose.getValue(),
                    tc_tophits.getValue(),
                    tc_locus.getValue(),
                    tc_isInSilico.getValue(),
                    tc_Id.getValue(),
                    tc_output.getValue());
            exit(EXIT_SUCCESS);
        }
        if(tc_jobType.getValue().compare("vcf") == 0)
        {
            TCLAP::ValueArg<string> tc_variants("s","variants","A text file containing the variant annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/genotype_to_phenotype_annotation.dat","string");
            cmdjob.add(tc_variants);
            TCLAP::ValueArg<string> tc_gt2pt("g","gt2pt","A text file containing the genotype to phenotype annotation of the ISBT. This file comes with the package an can be usually found in the subfolder data.",true,"data/config/genotype_to_phenotype_annotation.dat","string");
            cmdjob.add(tc_gt2pt);
            TCLAP::ValueArg<string> tc_alleleA("a","alleleA","First allele of the requested in silico data",true,"ABO*O.01.01","string");
            cmdjob.add(tc_alleleA);
            TCLAP::ValueArg<string> tc_alleleB("b","alleleB","First allele of the requested in silico data",true,"ABO*A2.01","string");
            cmdjob.add(tc_alleleB);
            cmdjob.parse(argc,argv);
            if(tc_verbose.getValue() >= 2)
                cerr << "Everything found. starting run ..." << endl;
            inSilicoVCF(tc_variants.getValue(),
                    tc_gt2pt.getValue(),
                    tc_alleleA.getValue(),
                    tc_alleleB.getValue(), 
                    tc_verbose.getValue());
            exit(EXIT_SUCCESS);
            
        }
        else
        {
            cerr << "job " << tc_jobType.getValue() << " not implemented." << endl;
            exit(EXIT_SUCCESS);
        }
    }
    catch(const TCLAP::ArgException& err)
    {
        cerr << "Error parsing command line: " << err.what() << endl;
        exit(EXIT_FAILURE);
    }
    catch(const CMyException& err)
    {
        cerr << err.what() << endl;
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
        
// ln -s ~/coding/cpp/deepBlood/data/example/bc1001.asm20.hg19.ccs.5passes.abotarget.bw coverage.bw
// ln -s ~/coding/cpp/deepBlood/data/example/bc1001.asm20.hg19.ccs.5passes.phased.phenotype.SNPs.vcf.gz SNPs.vcf.gz
void phenotype(const string& arg_target_anno,const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_vcf_file,
               const string& arg_bigWig,const string& arg_fastqgz, const string& arg_motifs,int arg_coverage, int arg_verbose, float arg_top_hits, const string& arg_locus, 
               bool arg_is_in_silico,const string& sampleId, const string& outfile)
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
        CMotifFinder mf(arg_motifs,arg_fastqgz,arg_verbose);
        cerr << mf << endl;
        if(arg_verbose >= 2)
            cerr << "Motif SNPs defined in " << arg_motifs << " looked up in "  << arg_fastqgz << endl;
        isbt.addCoverage(bwr,arg_coverage);
        std::set<string> loci = isbt.loci();
        
        
        vector<std::string> v = CMyTools::Tokenize(arg_locus,",");
        set<std::string> arg_loci_set(v.begin(), v.end());

        if(arg_verbose >= 3)
        {
            cerr << "ISBT annotation:" << endl << isbt << endl << endl;
            cerr << "ISBT genotype to phenotype translation:" << endl << isbTyper << endl << endl;
            CIsbtVariant::setVerbose();
         }
        // generate VCF Test Data
        CVariantChains vcs(&isbt);
        /// if you analyze simulated VCF do the following on your CVariantChains object
        if(arg_is_in_silico)
            vcs.removeReferenceSnps();
        //CVcf vcf_file(argv[3]);
        while(vcf_file.read_record())
        {
            CVcfSnp act_snp = vcf_file.get_record();
            //cerr << act_snp << endl;
            bool adding_successfull = vcs.add(act_snp);
            if(!adding_successfull && arg_verbose >= 3)
                cerr << act_snp << "\tSNP not added as it is not ISBT relevant" << endl;
            else if(adding_successfull && arg_verbose >= 2)
                cerr << act_snp << "\tISBT relevant SNP added" << endl;
        };
        vcs.removeUncoveredSnps(static_cast<double>(arg_coverage),arg_verbose);
        
        
        if(arg_verbose >= 3)
            cerr << "Variant chains of current sample: " << vcs << endl;   
        ofstream out_file;
        if(!outfile.empty())
        {
            out_file.open (outfile.c_str(), ios::out);
            if(!out_file.is_open())
                cerr << "failed to open " << outfile << endl;
        }
        
        nlohmann::json j;
        for(auto locus:loci)
        {
            if( arg_locus.length() == 0 || arg_loci_set.find(locus) != arg_loci_set.end() )
            {
                isbTyper.type(locus,vcs,arg_coverage);
                nlohmann::json jCall = isbTyper.getCallAsJson(isbt,trans_anno,bwr,locus,false,arg_top_hits,arg_coverage);
                j["loci"][locus]=jCall;
            }
        }
        j["sample_id"]=sampleId;
        j["version"]=APP_VERSION_DEEPBLOOD;
        j["genome"]="hg19";
        j["parameters"]["--target"]=arg_target_anno;
        j["parameters"]["--variants"]=arg_isbt_SNPs;
        j["parameters"]["--gt2pt"]=arg_genotype_to_phenotype;
        j["parameters"]["--vcf"]=arg_vcf_file;
        j["parameters"]["--bigwig"]=arg_bigWig;
        j["parameters"]["--coverage"]=arg_coverage;
        j["parameters"]["--verbose"]=arg_verbose;
        j["parameters"]["--scoreRange"]=arg_top_hits;
        j["parameters"]["--locus"]=arg_locus;
        j["parameters"]["--insilicovcf"]=arg_is_in_silico;
        j["parameters"]["--id"]=sampleId;
        j["parameters"]["--out"]=outfile;
        
        if(arg_verbose >= 3)
            cerr << isbTyper << endl;
        if(out_file.is_open())
        {
            out_file << j.dump();
            if(arg_verbose >= 2)
                cerr << "results written to " << outfile << endl;
        }
        else
            cout  << j.dump();
        if(out_file.is_open())
            out_file.close();
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

// # extract all alleles and mix
// cut -f 4 data/config/genotype_to_phenotype_annotation.dat | grep ABO | shuf > ~/ramDisk/mixed_ABO
// generate pairs of allelses from that shuffeled output (in excel)
// Example file mixed_ABO:
// ABO*O.01.26     ABO*cisAB.01
// ABO*A2.12       ABO*O.01.58
// ABO*AW.22       ABO*BW.04
// ...
// export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/mwittig/coding/cpp/MyTools/dist/Debug/GNU-Linux/:/home/mwittig/coding/fremd/htslib:/home/mwittig/coding/fremd/libBigWig
/* 
   for i in `cat mixed_ABO`
   do
    A1=`echo $i | cut -f 1`
    A2=`echo $i | cut -f 2`
    /home/mwittig/coding/cpp/deepBlood/dist/Debug/GNU-Linux/deepblood --job vcf --variants /home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation.dat --gt2pt /home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation.dat -a $A1 -b $A2 1> ${A1}_${A2}.vcf
   done
*/


void inSilicoVCF(const string& arg_isbt_SNPs,const string& arg_genotype_to_phenotype,const string& arg_allele_A,const string& arg_allele_B, int arg_verbose)
{
    try
    {
        CISBTAnno  isbt(arg_isbt_SNPs);
        if(arg_verbose >= 2)
            cerr << "ISBT variations loaded from:"  << arg_isbt_SNPs << endl;
        CIsbtGt2Pt isbTyper(arg_genotype_to_phenotype);
        if(arg_verbose >= 2)
            cerr << "ISBT genotype to phenotype translation loaded from:"  << arg_genotype_to_phenotype << endl;
        
            
        string systemA = isbTyper.systemOf(arg_allele_A);
        string systemB = isbTyper.systemOf(arg_allele_B);
        
        if(systemA.compare(systemB) != 0)
            throw CMyException(string("The system of both alleles must be the same. You have ")+arg_allele_A+"/"+arg_allele_B+" with '"+systemA+"'/'"+systemB+"'");
        
        CIsbtPtAllele alleleA = isbTyper.alleleOf(arg_allele_A);
        CIsbtPtAllele alleleB = isbTyper.alleleOf(arg_allele_B);
        
        cout << CMakeTrainingVcf::getHetEntries(systemA,alleleA,alleleB,isbt);
        
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


string getArgumentList(TCLAP::CmdLine& args)
{
    ostringstream osr("");
    
    const std::list<TCLAP::Arg*>& argL = args.getArgList();
    for(const auto& j:argL)
        osr << j->getName() << endl;
    
    return osr.str();
}


