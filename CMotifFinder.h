/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   CMotifFinder.h
 * Author: mwittig
 *
 * Created on January 5, 2022, 10:34 AM
 */

#ifndef CMOTIFFINDER_H
#define CMOTIFFINDER_H

class CMotifFinder {
public:
    CMotifFinder();
    virtual ~CMotifFinder();
    
    static std::map<string,int>  findMotifs(const std::string& filename,std::vector<string>& motifs);
    static void  findMotifs(const std::string& filename,std::map<string,int>& motifs);
    
private:

};

#endif /* CMOTIFFINDER_H */

