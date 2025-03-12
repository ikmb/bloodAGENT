/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   CFastqReader.cpp
 * Author: mwittig
 * 
 * Created on January 5, 2022, 10:02 AM
 */

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "gzstream.h"
#include "meinetools.h"
#include "CFastqReader.h"

using namespace std;

CFastqReader::CFastqReader(const std::string& filename) {
    m_in = NULL;
    if(!CMyTools::file_exists(filename))
        throw(CMyException("File does not exist: ")+filename);
    m_in = new igzstream(filename.c_str());
}

CFastqReader::CFastqReader(const CFastqReader& orig)
{
    m_in = orig.m_in;
}

CFastqReader::~CFastqReader() {
    if(m_in)
        delete m_in;
}

bool CFastqReader::getline(string& strLine)
{
    static string dummy;
    if(std::getline(*m_in,strLine))
    {
        return true;
    }
    return false;
}

