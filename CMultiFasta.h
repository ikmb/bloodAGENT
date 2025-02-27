/* 
 * File:   CMultiFasta.h
 * Author: mwittig
 *
 * Created on July 20, 2012, 1:00 PM
 */

#ifndef _CMULTIFASTA_H
#define	_CMULTIFASTA_H

using namespace std;

class CMultiFasta {
public:
    CMultiFasta();
    CMultiFasta(const string& filename);
    CMultiFasta(const CMultiFasta& orig);
    virtual ~CMultiFasta();

    size_t size(){return m_mapSequences.size();}

    map<string,string>::iterator begin(){return m_mapSequences.begin(); }
    map<string,string>::iterator end(){return m_mapSequences.end();}

    typedef map<string,string>::iterator iterator;
    typedef map<string,string>::const_iterator const_iterator;

    string operator[](const string&);
    
     
private:

    void readSequence(const string& filename);

    map<string,string> m_mapSequences;
    

};

#endif	/* _CMULTIFASTA_H */

