/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CMakeTrainingVcf.cpp
 * Author: mwittig
 * 
 * Created on July 30, 2019, 12:42 PM
 */

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>

#include "mytools.h"
#include "vcf.h"
#include "CVcf.h"
#include "CVcfSnp.h"
#include "CIsbtVariant.h"
#include "ISBTAnno.h"
#include "CVariantChain.h"
#include "CVariantChains.h"
#include "CIsbtPtAllele.h"
#include "CIsbtGtAllele.h"
#include "CIsbtGt2PtHit.h"

#include "CMakeTrainingVcf.h"

using namespace std;
/*
 * 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	UnnamedSample
chr1	25747230	.	G	C	457.03	.	AC=2;AF=1.0;AN=2;DP=16;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941	GT:AD:DP:GQ:PL	1/1:0,16:16:48:471,48,0

 chr4	145041720	.	A	G	251.60	.	AC=1;AF=0.5;AN=2;BaseQRankSum=-2.915;DP=17;ExcessHet=3.0103;FS=0.0;MLEAC=1;MLEAF=0.5;MQ=60.0;MQRankSum=0.0;QD=14.8;ReadPosRankSum=-0.293;SOR=0.813	GT:AD:DP:GQ:PL:PS	0|1:7,10:17:99:259,0,193:14502904
  
 */

CMakeTrainingVcf::CMakeTrainingVcf() 
{
    
}

CMakeTrainingVcf::CMakeTrainingVcf(const CMakeTrainingVcf& orig) 
{
    
}

CMakeTrainingVcf::~CMakeTrainingVcf() 
{
    
}

std::string CMakeTrainingVcf::getHomEntries(const std::string& system, const CIsbtPtAllele& allele, const CISBTAnno& anno)
{
    ostringstream osr("");
    int count = 0;
    for(const auto& variation:allele.baseChanges())
    {
        CISBTAnno::variation act_variant = anno.getIsbtVariant(system,variation);
        if(act_variant == CISBTAnno::variation())
        {
            cerr << "skipping unassigned " << variation << " of " << allele << endl;
            continue;
        }
        if(count++ != 0)
            osr << endl;
        osr << act_variant.chrom() << '\t'
            << act_variant.vcfCoordinate() << "\t.\t"
            << act_variant.vcfReference() << '\t'
            << act_variant.vcfAlternative() << '\t'
            << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t1/1:0,20:20:48:471,48,0";
    }
    return osr.str();
}

std::string CMakeTrainingVcf::getHetEntries(const std::string& system, const CIsbtPtAllele& alleleA, const CIsbtPtAllele& alleleB, const CISBTAnno& anno, bool phased)
{
    ostringstream osr("");
    std::set<std::string> variationsA = alleleA.baseChanges();
    std::set<std::string> variationsB = alleleB.baseChanges();
 
    set<CISBTAnno::variation> varSet;
    
    for(auto& var:variationsA)
        varSet.insert(anno.getIsbtVariant(system,var));
    for(auto& var:variationsB)
        varSet.insert(anno.getIsbtVariant(system,var));
    
    int count = 0;
    int phase_id = 0;
    for(auto& actVar:varSet)
    {
        bool isInA = variationsA.find(actVar.name()) != variationsA.end();
        bool isInB = variationsB.find(actVar.name()) != variationsB.end();
        bool hetA =  isInA && !isInB;
        bool hetB =  !isInA && isInB;
        bool homo =  isInA && isInB;
        
        if(!hetA && !hetB && !homo)
        {
            cerr << "skipping unassigned " << actVar << " of " << alleleA << '/' << alleleB << endl;
            continue;
        }
        if(count++ != 0)
            osr << endl;
        else
            phase_id=actVar.pos();
        osr << actVar.chrom() << '\t'
            << actVar.vcfCoordinate() << "\t.\t"
            << actVar.vcfReference() << '\t'
            << actVar.vcfAlternative() << '\t';
        if(phased)
        {
            if(hetA)
            {
                if(actVar.isRefNClikeGRChNC())
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL:PS\t1|0:0,20:20:48:471,48,0:"<<phase_id;
                else
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL:PS\t0|1:0,20:20:48:471,48,0:"<<phase_id;
            }
            if(hetB)
            {
                if(actVar.isRefNClikeGRChNC())
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL:PS\t0|1:0,20:20:48:471,48,0:"<<phase_id;
                else
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL:PS\t1|0:0,20:20:48:471,48,0:"<<phase_id;
            }
           if(homo)
           {
               if(actVar.isRefNClikeGRChNC())
                   osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t1|1:0,20:20:48:471,48,0";
               else
                   osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t0|0:0,20:20:48:471,48,0";
           }
        }
        else
        {
            // UNPHASED
            if(hetA)
            {
                osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t0/1:0,20:20:48:471,48,0";
            }
            if(hetB)
            {
                osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t0/1:0,20:20:48:471,48,0";
            }
            if(homo)
            {
                if(actVar.isRefNClikeGRChNC())
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t1/1:0,20:20:48:471,48,0";
                else
                    osr << "450.0\t.\tAC=2;AF=1.0;AN=2;DP=20;ExcessHet=3.0103;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=59.69;QD=28.56;SOR=0.941\tGT:AD:DP:GQ:PL\t0/0:0,20:20:48:471,48,0";
            }
        }
        
       
        
    }
     return osr.str();
}




