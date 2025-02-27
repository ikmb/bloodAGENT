/* 
 * File:   CMyException.h
 * Author: mwittig
 *
 * Created on June 18, 2010, 11:13 AM
 */

#ifndef _CMYEXCEPTION_H
#define	_CMYEXCEPTION_H


class CMyException
{
public:
    CMyException () throw();
    CMyException (const std::string&) throw();
    CMyException (const char*) throw();
    CMyException (const CMyException&) throw();
    CMyException& operator= (const CMyException&) throw();
    CMyException operator+(const std::string&);
    CMyException operator+(const int&);
    CMyException operator+(const size_t&);
    friend std::ostream& operator<<(std::ostream& os, const CMyException& me){os << me.what(); return os;}
    virtual std::string what() const throw();

private:
    std::string  m_strErr;
};



#endif	/* _CMYEXCEPTION_H */

