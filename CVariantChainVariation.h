/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChainVariation.h
 * Author: mwittig
 *
 * Created on July 22, 2019, 11:06 AM
 */

#ifndef CVARIANTCHAINVARIATION_H
#define CVARIANTCHAINVARIATION_H

#include "CIsbtVariant.h"


class CVariantChainVariation
{
public:
    CVariantChainVariation();
    CVariantChainVariation(const CIsbtVariant& var){first_variant=second_variant=var;}
    CVariantChainVariation(const CVariantChainVariation& orig);
    virtual ~CVariantChainVariation();
    
    CVariantChainVariation& operator =(const CVariantChainVariation& orig);
    bool                    operator <(const CVariantChainVariation& orig)const;
    bool                    operator ==(const CVariantChainVariation& orig)const;
    bool                    operator !=(const CVariantChainVariation& orig)const{return !(*this == orig);}
    
    
    
    CIsbtVariant first_variant;
    CIsbtVariant second_variant;
    
    
    
private:

    
    
};

#endif /* CVARIANTCHAINVARIATION_H */

