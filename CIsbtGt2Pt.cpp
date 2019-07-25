/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CIsbtGt2Pt.cpp
 * Author: mwittig
 * 
 * Created on July 25, 2019, 7:00 AM
 */

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <libgen.h>

#include "mytools.h"
#include "CIsbtGt2Pt.h"

using namespace std;

CIsbtGt2Pt::CIsbtGt2Pt(const string& filename) 
{
    init(filename);
}

CIsbtGt2Pt::CIsbtGt2Pt(const CIsbtGt2Pt& orig) 
{
    m_allele_vector = orig.m_allele_vector;
}

CIsbtGt2Pt::~CIsbtGt2Pt() {
}

void CIsbtGt2Pt::init(const string& filename)
{
    CParsedTextfile ptx(filename,"\t","Phenotype",0,true, "#");
    if(ptx.First())
        do{
            m_allele_vector[ptx["System"]].push_back(CIsbtPtAllele(ptx["Phenotype"], ptx["base_change"], ptx["acid_change"], ptx["incidence"]));
        }while(ptx.Next());
    
}

std::ostream& operator<<(std::ostream& os, const CIsbtGt2Pt& me)
{
    long unsigned int i = 0;
    for(std::map<std::string,vector<CIsbtPtAllele>>::const_iterator iter = me.m_allele_vector.begin(); iter != me.m_allele_vector.end(); iter++)
    {
        for(auto x:iter->second)
        {
            os <<  iter->first << '\t' << x;
            if( ++i != iter->second.size())
                os << endl;
        }
    }
    return os;
}



