/*************************************************
 * 
    CMyTools , is a collection of frequently use functions
    Copyright (C) [2006]  [Michael Wittig] (m.wittig@mucosa.de  "www.ikmb.uni-kiel.de")

    This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
    as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program; 
    if not, write to the Free Software Foundation, 
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110, USA

***************************************************/

#include <string>
#include <string.h>
#include <assert.h>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <dirent.h>
#include <stdlib.h>
#include <set>
#include "MyTools.h"

using namespace std;

CMyTools::CMyTools()
{
}

CMyTools::~CMyTools()
{
}

bool CMyTools::file_exists(const string& fileName)
{
    ifstream infile(fileName.c_str());
    return infile.good();
}

vector<string> CMyTools::GetParsedLine(string str,const string& sep)
{
    vector<string> tokens;
    // Skip delimiters at beginning.
    //string::size_type lastPos = str.find_first_not_of(sep, 0);
    string::size_type lastPos = 0;
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(sep, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        if(str.length() > pos && sep.find(str[pos]) != string::npos )
            lastPos=pos+1;
        else
            lastPos = str.find_first_not_of(sep, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(sep, lastPos);
    }
    return tokens;
}


// Bitte bei Gelegenheit optimieren !!!!!!!!!!!!
string	CMyTools::ReplaceSubstr(std::string& result, const std::string& replaceWhat,  const std::string& replaceWithWhat)
{
  while(1)
  {
    const int pos = result.find(replaceWhat);
    if (pos==-1) break;
    result.replace(pos,replaceWhat.size(),replaceWithWhat);
  }
  return result;
}



vector<string> CMyTools::GetDirectoryListing(const string& path)
{
    DIR *dirp;
    struct dirent *entry;
    vector<string> strReturn;
    
    if(dirp = opendir(path.c_str()))
    {
        while(entry = readdir(dirp))
            strReturn.push_back(entry->d_name);
        closedir(dirp);
    }   
    return strReturn;
}

vector<string> CMyTools::GetDirectoryListing(const string& path,const char* contains, bool full_path)
{
    vector<string> strReturn;
    vector<string> strhelp = GetDirectoryListing(path);
    for(vector<string>::iterator iter = strhelp.begin();iter != strhelp.end();iter++)
        if(iter->find(contains)!=string::npos)
        {
            if(full_path != true)
                strReturn.push_back(*iter);
            else
                strReturn.push_back(string(path)+"/"+(*iter));
        }
    
    return strReturn;
}
	

bool CMyTools::isNumeric(const char *val)
{
	try
	{
		ostringstream o;
		o << atof(val);
//		cout << "(" << val << ")" << std::atof(val) << " - (" << o.str() << ")" << std::atof(o.str().c_str()) << " = " << std::atof(o.str().c_str())-std::atof(val) << endl;
		if(atof(o.str().c_str())-atof(val) != 0)
			throw o.str();
		return true;
	}
	catch(string v)
	{
		throw string("converting const char* to float failed for ")+val+" -> "+v;	
	}
	catch(...)
	{
		throw string("converting const char* to float failed for ")+val;	
	}
}


string CMyTools::GetComplSequence(const string& strSequence)
{
	try
    {
        int intPos = 0;
    	size_t intSize = strSequence.length();
    	char* buffer=(char*)malloc(intSize+1);
        if(!buffer)
            throw(string("can't allocate enough memory."));
        buffer[intSize]=0;
        while(intSize > 0)
        {
            buffer[intPos++]=GetComplNuc(strSequence.at(--intSize));
        }        
       	string strReturn(buffer);
        free(buffer);
        return strReturn;
    }
    catch(string strError)
    {
        throw(string("error in CMyTools::GetSeqFromFile, ")+strError); 
    }  
    catch(...)  
    {
        throw("unknown error in CMyTools::GetSeqFromFile\n");                 
    }
	
}

vector<unsigned int> CMyTools::FindApproxSubs(unsigned int k, string & P, string & T)
{
  // Aufgabe: Finde alle Positionen in T, an denen ein Substring von T mit P 
  // ann�hernd �bereinstimmt (mit Abstand <= k), und gebe diese in Out aus.
	
  unsigned int n = (unsigned int )T.length();
  unsigned int m = (unsigned int )P.length();
  unsigned int *C = new unsigned int[m];
  unsigned int last = k + 1;
  unsigned int pC,nC, i, pos;
	
  vector<unsigned int>  Out;

  for(i = 0; i < m; i++) { C[i] = i + 1; }

  // Searching
  for(pos = 0; pos < n; pos++) 
  {
		pC = 0;
		nC = 0;
		//for(i = 0; i < last; i++) 
		for(i = 0; i <  last; i++) 
		{
			if(P[i] == T[pos]) 
			{
				nC = pC;
			} 
			else 
			{
				if(pC < nC) 
					nC = pC;
				if(C[i] < nC) 
					nC = C[i];
				nC++;
			}
			pC = C[i];
			C[i] = nC;
		}
		while(C[last - 1] > k) 
			last--;
		if(last == m) 
		  Out.push_back(pos + 1);
		else
			last++; 
  }
  delete C;
  return Out;
}


char CMyTools::GetComplNuc(const char& source)
{
	switch(source)
	{
	case 'A':
		return 'T';
	case 'a':
		return 't';
	case 'C':
		return 'G';
	case 'c':
		return 'g';
	case 'G':
		return 'C';
	case 'g':
		return 'c';
	case 'T':
	case 'U':
		return 'A';
	case 't':
	case 'u':
		return 'a';
	case 'M':
		return 'K';
	case 'm':
		return 'k';
	case 'R':
		return 'Y';
	case 'r':
		return 'y';
	case 'Y':
		return 'R';
	case 'y':
		return 'r';
	case 'K':
		return 'M';
	case 'k':
		return 'm';
	case 'V':
		return 'B';
	case 'v':
		return 'b';
	case 'H':
		return 'D';
	case 'h':
		return 'd';
	case 'D':
		return 'H';
	case 'd':
		return 'h';
	case 'B':
		return 'V';
	case 'b':
		return 'v';
	case '(':
         return ')';
    case ')':
         return '(';
	case '[':
         return ']';
    case ']':
         return '[';
	}
	return source;
}

std::pair<std::string, char> CMyTools::getAminoAcid(std::string triplet)
{
    for(int i = 0; i < triplet.size(); i++)
        triplet[i] = toupper(triplet[i]);
    if(triplet.size() != 3)
        return pair<string,char>("Xaa",'X');
    switch(triplet[0])
    {
        case 'A':
            switch(triplet[1])
            {
                case 'A':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Lys",'K');
                        case 'C':
                            return pair<string,char>("Asn",'N');
                        case 'G':
                            return pair<string,char>("Lys",'K');
                        case 'T':
                            return pair<string,char>("Asn",'N');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'C':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Thr",'T');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'G':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Arg",'R');
                        case 'C':
                            return pair<string,char>("Ser",'S');
                        case 'G':
                            return pair<string,char>("Arg",'R');
                        case 'T':
                            return pair<string,char>("Ser",'S');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'T':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                            return pair<string,char>("Ile",'I');
                        case 'G':
                            return pair<string,char>("Met",'M');
                        case 'T':
                            return pair<string,char>("Ile",'I');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                default:
                    return pair<string,char>("Xaa",'X');
            }
            break;
        case 'C':
            switch(triplet[1])
            {
                case 'A':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Gln",'Q');
                        case 'C':
                            return pair<string,char>("His",'H');
                        case 'G':
                            return pair<string,char>("Gln",'Q');
                        case 'T':
                            return pair<string,char>("His",'H');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'C':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Pro",'P');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'G':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Arg",'R');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'T':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Leu",'L');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                default:
                    return pair<string,char>("Xaa",'X');
            }
            break;
        case 'G':
            switch(triplet[1])
            {
                case 'A':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Glu",'E');
                        case 'C':
                            return pair<string,char>("Asp",'D');
                        case 'G':
                            return pair<string,char>("Glu",'E');
                        case 'T':
                            return pair<string,char>("Asp",'D');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'C':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Ala",'A');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'G':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Gly",'G');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'T':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Val",'V');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                default:
                    return pair<string,char>("Xaa",'X');
            }
            break;
        case 'T':
            switch(triplet[1])
            {
                case 'A':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Stop",'.');
                        case 'C':
                            return pair<string,char>("Tyr",'Y');
                        case 'G':
                            return pair<string,char>("Stop",'.');
                        case 'T':
                            return pair<string,char>("Tyr",'Y');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'C':
                    switch(triplet[2])
                    {
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                            return pair<string,char>("Ser",'S');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'G':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Stop",'.');
                        case 'C':
                            return pair<string,char>("Cys",'C');
                        case 'G':
                            return pair<string,char>("Trp",'W');
                        case 'T':
                            return pair<string,char>("Cys",'C');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                case 'T':
                    switch(triplet[2])
                    {
                        case 'A':
                            return pair<string,char>("Leu",'L');
                        case 'C':
                            return pair<string,char>("Phe",'F');
                        case 'G':
                            return pair<string,char>("Leu",'L');
                        case 'T':
                            return pair<string,char>("Phe",'F');
                        default:
                            return pair<string,char>("Xaa",'X');
                    } 
                default:
                    return pair<string,char>("Xaa",'X');
            }
        default:
            return pair<string,char>("Xaa",'X');
        
    }
}

char CMyTools::decodeBaseChange(char base, char cc)
{
	try
	{
		switch(base)
		{
			case 'a':
			case 'A':
				switch(cc)
				{
					case '0':return 'A';break;
					case '1':return 'C';break;
					case '2':return 'G';break;
					case '3':return 'T';break;
					default:return 'N';
				}
				break;
			case 'c':
			case 'C':
				switch(cc)
				{
					case '0':return 'C';break;
					case '1':return 'A';break;
					case '2':return 'T';break;
					case '3':return 'G';break;
					default:return 'N';
				}
				break;
			case 'g':
			case 'G':
				switch(cc)
				{
					case '0':return 'G';break;
					case '1':return 'T';break;
					case '2':return 'A';break;
					case '3':return 'C';break;
					default:return 'N';
				}
				break;
			case 't':
			case 'T':
				switch(cc)
				{
					case '0':return 'T';break;
					case '1':return 'G';break;
					case '2':return 'C';break;
					case '3':return 'A';break;
					default:return 'N';
				}
				break;
			default:
				return 'N';
		}		
	}
	catch(const char* err)
	{
		throw(err);
	}
	catch(string& err)
	{
		throw(err);
	}
	catch (exception& err)
        {
      		throw(err);
    	}
    	catch(...)
	{
		throw(1);
	}	
}


string CMyTools::getBasespace(string cs,bool keep_lead_base)
{
	try
	{
		string	bs;
		char	lead_base;
		string::iterator iter = cs.begin();
                if(keep_lead_base)bs=*iter;
		if(iter!=cs.end())
		{
			lead_base=*iter;
			for(iter++;iter!=cs.end();iter++)
			{
				lead_base=decodeBaseChange(lead_base,*iter);
				bs+=lead_base;
			}
		}
		return bs;		
	}
	catch(const char* err)
	{
		throw(err);
	}
	catch(string& err)
	{
		throw(err);
	}
	catch (exception& err)
    {
  		throw(err);
	}
	catch(...)
	{
		throw(1);
	}			
}

string CMyTools::getBasespaceRev(string cs,bool keep_lead_base)
{
	try
	{
		string	bs;
                char	lead_base;
		string::reverse_iterator iter = cs.rbegin();
                if(keep_lead_base)bs=*iter;
		if(iter!=cs.rend())
		{
			lead_base=*iter;
			for(iter++;iter!=cs.rend();iter++)
			{
				lead_base=decodeBaseChange(lead_base,*iter);
				bs+=lead_base;
			}
		}
		return bs;		
	}
	catch(const char* err)
	{
		throw(err);
	}
	catch(string& err)
	{
		throw(err);
	}
	catch (exception& err)
    {
  		throw(err);
	}
	catch(...)
	{
		throw(1);
	}			
}

template<class T> void CMyTools::GetMem(T*& OldMem, int Elems)
{
    typedef int cntr; // Type of element cntr
    const int CountSize = sizeof(cntr); // And size
    const int TypeSize = sizeof(T);
    if( Elems == 0) 
    {
        free( &(((cntr*)OldMem)[-1]) );
        return;
    }

    T* p = OldMem;
    cntr OldCount = 0;
    if(p) 
    { // Previously allocated memory
        cntr* tmp = reinterpret_cast<cntr*>(p);
        p = reinterpret_cast<T*>(--tmp);
        OldCount = *(cntr*)p; // Previous # Elems
    }

    T* m = (T*)realloc(p, Elems * TypeSize + CountSize);
    assert(m != 0);
    *((cntr*)m) = Elems; // Keep track of count
    const cntr Increment = Elems - OldCount;
    if( Increment > 0) 
    {
        // Starting address of data:
        long StartAdr = (long)&(m[OldCount]);
        StartAdr += CountSize;
        // Zero the additional new memory:
        memset((void*)StartAdr, 0, Increment * TypeSize);
    }
    // Return the address beyond the count:
    OldMem = (T*)&(((cntr*)m)[1]);
}

string CMyTools::dtoa(long value, int base) 
{
    string strReturn;
    strReturn.reserve( MYTOOLS_H_MAX_ITOA_DIGITS_SIZE ); 
    const   char  translator[] = "0123456789abcdef";
    
    bool    negative = value < 0;
    value   = abs(value);
    if (base < 2 || base > 16) 
        return strReturn;

     while (value)
     {
        strReturn += translator[value % base];
        value /= base;
     }

    if (negative) 
        strReturn += '-';
    
    reverse( strReturn.begin(), strReturn.end() );
    return strReturn;
}

bool CMyTools::doIntersect(int a_start,int a_end,int b_start,int b_end)
{
    if(a_start > a_end)
        swap(a_start,a_end);
    if(b_start > b_end)
        swap(b_start,b_end);
    
    if(a_end < b_start || a_start > b_end)
        return false;
    return true;
}

size_t CMyTools::intersectCount(int a_start,int a_end,int b_start,int b_end)
{
    if(a_start > a_end)
        swap(a_start,a_end);
    if(b_start > b_end)
        swap(b_start,b_end);
    
    if(a_end < b_start || a_start > b_end)
        return 0;
    else
    {
        size_t r_start = std::max(a_start,b_start);
        size_t r_end = std::min(a_end,b_end);
        return (r_end-r_start)+1;
    }
}

size_t CMyTools::insertSize(int a_start,int a_end,int b_start,int b_end)
{
    size_t r_start = std::min(std::min(a_start,b_start),std::min(a_end,b_end));
    size_t r_end = std::max(std::max(a_start,b_start),std::max(a_end,b_end));
    return (r_end-r_start)+1;
}

float CMyTools::bytesToFloat(unsigned char* byte_array) 
{ 
    float result;
    std::copy(reinterpret_cast<const char*>(&byte_array[0]),
              reinterpret_cast<const char*>(&byte_array[4]),
              reinterpret_cast<char*>(&result));
    return result;
} 

void CMyTools::erase_windows_line_ending(string& s)
{
    if (!s.empty() && s[s.size() - 1] == '\r')
        s.erase(s.size() - 1);
}

int CMyTools::cmd(const string& commando, bool write_output)
{
    FILE *in;
    char buff[512];

    if(!(in = popen(commando.c_str(), "r")))
    {
        return 1;
    }

    while(fgets(buff, sizeof(buff), in)!=NULL)
        if(write_output)
            cout << buff;
    pclose(in);
    if(write_output)
        cout << endl;
    return 0;
}

int CMyTools::cmd(const string& commando, std::string& get_output)
{
    FILE *in;
    char buff[512];

    if(!(in = popen(commando.c_str(), "r")))
    {
        return 1;
    }

    while(fgets(buff, sizeof(buff), in)!=NULL)
        get_output.append(buff);
    pclose(in);
    if (!get_output.empty() && get_output[get_output.length()-1] == '\n') 
    {
        get_output.erase(get_output.length()-1);
    }
    return 0;
}

string CMyTools::toString(std::set<std::string> v, const char sep)
{
    ostringstream s;
    
    for(std::set<std::string>::iterator i = v.begin(); i != v.end(); i++)
    {
        if(i!=v.begin())
            s << sep <<  *i;
        else 
            s << *i;
    }
    
    return s.str();
}

std::string CMyTools::trim(const std::string &s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && isspace(*it))
        it++;

    std::string::const_reverse_iterator rit = s.rbegin();
    while (rit.base() != it && isspace(*rit))
        rit++;

    return std::string(it, rit.base());
}

