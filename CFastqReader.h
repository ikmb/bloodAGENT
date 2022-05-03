/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   CFastqReader.h
 * Author: mwittig
 *
 * Created on January 5, 2022, 10:02 AM
 */

#ifndef CFASTQREADER_H
#define CFASTQREADER_H

class CFastqReader {
public:
    CFastqReader(const std::string& filename);
    CFastqReader(const CFastqReader& orig);
    virtual ~CFastqReader();
    
    bool getline(std::string& strLine);
    
private:
    
    igzstream* m_in;
};

#endif /* CFASTQREADER_H */

