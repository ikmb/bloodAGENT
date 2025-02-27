/* 
 * File:   CMyException.cpp
 * Author: mwittig
 * 
 * Created on June 18, 2010, 11:13 AM
 */

#include <string>

using namespace std;

#include "CMyException.h"

CMyException::CMyException()throw() 
{
    
}
CMyException::CMyException(const char* orig) throw()
{
    m_strErr = orig;
}
CMyException::CMyException(const std::string& orig) throw() {
    m_strErr = orig;
}
CMyException::CMyException(const CMyException& orig) throw(){
    m_strErr = orig.m_strErr;
}

CMyException& CMyException::operator= (const CMyException& orig) throw()
{
    m_strErr = orig.m_strErr;
    return *this;
}

CMyException CMyException::operator+(const std::string& value)
{
    m_strErr.append(value);
    return *this;
}

CMyException CMyException::operator+(const int& v)
{
    m_strErr.append(std::to_string(v));
    return *this;
}

CMyException CMyException::operator+(const size_t& v)
{
    m_strErr.append(std::to_string(v));
    return *this;
}

string CMyException::what() const throw()
{
    return m_strErr;
}