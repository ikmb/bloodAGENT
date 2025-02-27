#include <string>
#include <map>
#include <vector>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <set>

#include "CMyException.h"
#include "MyTools.h"
#include "ParsedTextfile.h"


using namespace std;

CParsedTextfile::CParsedTextfile()  : m_comment_char("")
{
    cursor=vctLines.begin();
    m_strValidSeparator="\t";
    errors.str("");
    m_skip = 0;
    m_row_key = "";
    m_row_key_column = -1;
}


CParsedTextfile::CParsedTextfile(const CParsedTextfile& orig) : m_comment_char(orig.m_comment_char)
{
	*this=orig;
}

CParsedTextfile::CParsedTextfile(const string& filename, const string& pchrSeperator,const string& row_key, int skip, bool read_header, const string& comment) : m_comment_char(comment)
{
    ifstream in(filename);
    if(!in)
        throw CMyException(string("failed to open file: ")+filename);
    m_skip = skip;
    errors.str("");
    m_strValidSeparator=pchrSeperator;
    m_row_key = row_key;
    m_row_key_column = -1;
    if(!ReadFile(in,read_header))
        throw CMyException("Error reading flat file!");
    First();

    in.close();
}

CParsedTextfile::CParsedTextfile(const string& filename, const string& pchrSeperator,const string& row_key, const string& read_until, bool read_header, const string& comment) : m_comment_char(comment)
{
    ifstream in(filename);
    if(!in)
        throw CMyException(string("failed to open file: ")+filename);
    m_skip = 0;
    errors.str("");
    m_strValidSeparator=pchrSeperator;
    m_row_key = row_key;
    m_row_key_column = -1;
    if(!ReadFile(in,read_header,read_until))
        throw CMyException("Error reading flat file!");
    First();

    in.close();
}

CParsedTextfile::CParsedTextfile(const string& filename,const string& sep,int row_key, int skip, bool read_header, const string& comment) : m_comment_char(comment)
{
    ifstream in(filename);
    if(!in)
        throw CMyException(string("failed to open file: ")+filename);
    m_skip = skip;
    errors.str("");
    m_strValidSeparator=sep;
    m_row_key = "";
    m_row_key_column = row_key;
    if(!ReadFile(in,read_header))
        throw CMyException("Error reading flat file!");
    First();

    in.close();
}


CParsedTextfile::CParsedTextfile(ifstream& f, const string& pchrSeperator,int row_key, int skip, bool read_header, const string& comment) : m_comment_char(comment)
{
   try
    {
        m_skip = skip;
        errors.str("");
        m_strValidSeparator=pchrSeperator;
        m_row_key_column = row_key;
        m_row_key = "";
        if(!ReadFile(f,read_header))
            throw CMyException("Error reading flat file!");
        First();
    }
    catch(const char* err)
    {
        cerr << "error reading file." << endl;
    } 
}
	
CParsedTextfile::CParsedTextfile(ifstream& f, const string& pchrSeperator,const string& row_key, int skip, bool read_header, const string& comment) : m_comment_char(comment)
{
    try
    {
        m_skip = skip;
        errors.str("");
        m_strValidSeparator=pchrSeperator;
        m_row_key = row_key;
        m_row_key_column = -1;
        if(!ReadFile(f,read_header))
            throw CMyException("Error reading flat file!");
        First();
    }
    catch(const char* err)
    {
        cerr << "error reading file." << endl;
    }
}

CParsedTextfile::~CParsedTextfile()
{
}

CParsedTextfile& CParsedTextfile::operator=(const CParsedTextfile& orig)
{
    mapHeader=orig.mapHeader;
    vctLines=orig.vctLines;
    errors << orig.errors.str().c_str();
    cursor=vctLines.begin();
    m_strValidSeparator=orig.m_strValidSeparator;
    m_skip = orig.m_skip;
    m_row_key = orig.m_row_key;
    m_row_key_column = orig.m_row_key_column;
    mapRow = orig.mapRow;
    return *this;
}


bool CParsedTextfile::CreateHeader(const string& chrLine)
{
    try
    {
        bool bReturn = true;
        int i=0;

        vector<string>  entries(CMyTools::Tokenize(chrLine,m_strValidSeparator));
        if(mapHeader.size()==0)
            for(vector<string>::iterator iter=entries.begin();iter != entries.end();iter++, i++)
                mapHeader[*iter]=i;
        else
            bReturn = false;		
        return bReturn;
    }
    catch(...)
    {
            return false;	
    }
}

void CParsedTextfile::CoutHeader()
{
    for(map<string,int>::iterator i = mapHeader.begin(); i != mapHeader.end(); i++)
        cout << i->first << " - " << i->second << endl;
}

string CParsedTextfile::getHeaderAsString(const string& sep)
{
    ostringstream strRet("");
    map<int,string> sorted;
    for(map<string,int>::iterator i = mapHeader.begin(); i != mapHeader.end(); i++)
        sorted[i->second]=i->first;
    for(map<int,string>::iterator i = sorted.begin(); i != sorted.end(); i++)
        strRet << (i == sorted.begin() ? "" : sep) << i->second;
    return strRet.str();
}

vector<string> CParsedTextfile::getHeader()
{
    vector<string> vRet;
    
    vRet.resize(mapHeader.size());
    for(map<string,int>::iterator i = mapHeader.begin(); i != mapHeader.end(); i++)
        vRet[i->second]=i->first;
    return vRet;
}

bool CParsedTextfile::ReadFile(ifstream& f, bool read_header, string read_until)
{
    try
    {
        size_t lineCOunt=0;
        bool bReturn = true;
        string	strLine;
        int line_entries = -1;
        for( int i = 0; i < m_skip; i++)
        {
            if(!getline(f,strLine))
                break;
            lineCOunt++;
        }
        if(!read_until.empty())
        {
            while(getline(f,strLine))
            {
                lineCOunt++;
                CMyTools::erase_windows_line_ending(strLine);
                if(strLine.compare(read_until) == 0)
                {
                    cerr << lineCOunt << " lines until " << read_until << endl;
                    break;
                }
            }
        }
        if(read_header)
        {
            while(getline(f,strLine))
            {
                lineCOunt++;
                if(strLine.size() == 0 || (!m_comment_char.empty() && strLine.find_first_of(m_comment_char) == 0 ) )
                    continue;
                CMyTools::erase_windows_line_ending(strLine);
                CreateHeader(strLine);
                //if( CreateHeader(strLine) )
                //    cerr << "header red." << endl;
                break;
            }
        }
        
        bool make_row_map = m_row_key != "" && mapHeader.find(m_row_key) != mapHeader.end();
        int rowKeyColumnNr = 0;
        if(make_row_map)
            rowKeyColumnNr = mapHeader.find(m_row_key)->second;
        else if(m_row_key_column != -1)
        {
            rowKeyColumnNr = m_row_key_column;
            make_row_map = true;
        }
        lineCOunt++;
        while(getline(f,strLine))	
        {
            lineCOunt++;
            if(strLine.size() == 0 || strLine[0] == '#')
                continue;
            CMyTools::erase_windows_line_ending(strLine);
            ParsedTextfile_OneLine	s=CMyTools::Tokenize(strLine,m_strValidSeparator.c_str());
            
            if( (read_header && s.size() == mapHeader.size()) || s.size() == line_entries || (line_entries == -1 && read_header == false) )
            {
                vctLines.push_back(s);
                if(make_row_map)
                {
                    if(mapRow.find(s[rowKeyColumnNr]) != mapRow.end())
                    {
                        cerr << "Error: row key " << s[rowKeyColumnNr] << " was already used in row " << mapRow[s[rowKeyColumnNr]] << endl;
                        return false;
                    }
                    mapRow[s[rowKeyColumnNr]]=vctLines.size()-1;
                }
            }
            else
            {
                if( line_entries != -1 && s.size() != line_entries)
                    cerr << "invalid entry at line " << lineCOunt << ": last line had " << line_entries << " columns, current line has " 
                         << s.size() << ". (" << strLine << ")" << endl;
                else
                    cerr << "invalid entry at line " << lineCOunt << ": header has " << mapHeader.size() << " columns, entry has " 
                         << s.size() << ". (" << strLine << ")" << endl;
                //bReturn = false;
            }
            line_entries = s.size();
        }
        return bReturn;
    }
    catch(...)
    {
        cerr << errors.str() << endl;
        return false;	
    }	
}

bool CParsedTextfile::CreateRowMap()
{
    bool make_row_map = (m_row_key != "" && mapHeader.find(m_row_key) != mapHeader.end()) || m_row_key_column != -1;
    if(make_row_map)
    {
        int rowKeyColumnNr = m_row_key_column;
        if(m_row_key != "" && mapHeader.find(m_row_key) != mapHeader.end())
            rowKeyColumnNr = mapHeader.find(m_row_key)->second;
        int counter = 0;
        mapRow.clear();
        for(vector<ParsedTextfile_OneLine>::iterator i = vctLines.begin(); i != vctLines.end(); i++)
        {
            mapRow[(*i)[rowKeyColumnNr]]=counter++;
        }
    }  
    return true;
}

bool CParsedTextfile::Next()
{
    if(!isEOF())		
    {
            cursor++;
            if(!isEOF())
                    return true;
            cursor--;
    }
    return false;		
}

bool  CParsedTextfile::Prev()
{
    if(cursor != vctLines.begin() && !isEOF())		
    {
        cursor--;
        return true;
    }
    return false;
}

bool  CParsedTextfile::First()
{
    cursor=vctLines.begin();
    if(!isEOF())	
        return true;
    return false;
}

bool  CParsedTextfile::Last()
{
    if(vctLines.size()>0)
    {
            cursor=vctLines.end();
            cursor--;
            return true;
    }
    return false;		
}

bool CParsedTextfile::erase(CParsedTextfile::iterator del_this)
{
    if(del_this != vctLines.end())
    {
        cursor = vctLines.erase(del_this);
        if(cursor != vctLines.begin())
            cursor--;
        CreateRowMap();
        return true;
    }
    return false;
}


bool CParsedTextfile::isEOF() 
{
    if(cursor==vctLines.end())	
            return true;
    return false;	
}

string CParsedTextfile::operator[](const string& strCol)
{
    if(isEOF())
        throw CMyException(string("End of file reached"));
    if(mapHeader.find(strCol) == mapHeader.end())
    {
        throw CMyException(string("Invalid header column name: ")+strCol);
    }
    return ((*cursor)[mapHeader[string(strCol)]]);
}

string CParsedTextfile::operator[](const int idx)
{
    if(isEOF())
        throw CMyException(string("End of file reached"));
    if( !(idx >= 0 && idx < cursor->size()) )
    {
        throw CMyException(string("Column index out of bounds"));
    }
    return ((*cursor)[idx]);
}

string CParsedTextfile::cell(const string& rowKey,const string& colKey)const
{
    map<string,int>::const_iterator     imapHeader = mapHeader.find(colKey);
    map<string,size_t>::const_iterator  imapRow    = mapRow.find(rowKey);
	
    if(imapRow == mapRow.end())
    {
        throw CMyException(string("Row key not found: ")+rowKey);
    }
    
    if(imapHeader == mapHeader.end())
    {
        throw CMyException(string("Column key not found: ")+colKey);
    }
    
    return cell(imapRow->second,imapHeader->second);
}

CParsedTextfile::iterator CParsedTextfile::find(const string& rowKey)
{
    map<string,size_t>::iterator  imapRow    = mapRow.find(rowKey);
	
    if(imapRow == mapRow.end())
    {
        cursor = vctLines.end();;
    }
    else
        cursor = vctLines.begin()+imapRow->second;
    return cursor;
}

string CParsedTextfile::cell(const string& rowKey,int colIdx)const
{
    map<string,size_t>::const_iterator  imapRow    = mapRow.find(rowKey);
	
    if(imapRow == mapRow.end())
    {
        throw CMyException(string("Row key not found: ")+rowKey);
    }
    return cell(imapRow->second,colIdx);
}

string CParsedTextfile::cell(int row,const string& colKey)const
{
    map<string,int>::const_iterator     imapHeader = mapHeader.find(colKey);
    if(imapHeader == mapHeader.end())
    {
        throw CMyException(string("Column key not found: ")+colKey);
    }
    if(row < vctLines.size())
    {
        return cell(row,imapHeader->second);
    }
    else
        throw CMyException(string("Row index out of bounds: "));
}

string CParsedTextfile::cell(int row,int col)const
{
    if(row < vctLines.size())
    {
        const ParsedTextfile_OneLine& line = vctLines[row];
        if(col < line.size())
            return line[col];
        throw CMyException(string("Column index out of bounds"));
    }
    else
        throw CMyException(string("Row index out of bounds: "));
}

vector<string> CParsedTextfile::cellTok(const string& rowKey,const string& colKey, const string& sep)
{
    return CMyTools::Tokenize(cell(rowKey,colKey),sep);
}

vector<string> CParsedTextfile::cellTok(const string& rowKey,int col, const string& sep)
{
    return CMyTools::Tokenize(cell(rowKey,col),sep);
}

vector<string> CParsedTextfile::cellTok(int row,const string& colKey, const string& sep)
{
    return CMyTools::Tokenize(cell(row,colKey),sep);
}

vector<string> CParsedTextfile::cellTok(int row,int col, const string& sep)
{
    return CMyTools::Tokenize(cell(row,col),sep);
}    

int CParsedTextfile::GetColId(const string& strCol)
{
    if(mapHeader.find(strCol)==mapHeader.end())
            return -1;
    return mapHeader[string(strCol)];	
}

string	CParsedTextfile::line()
{
    return line(m_strValidSeparator.c_str());	
}

string	CParsedTextfile::line(int idx, const string& separator)
{
    ostringstream ostrRet("");
    int count = 0;
    if(idx < vctLines.size() )
    {
        for(ParsedTextfile_OneLine::iterator iter=vctLines[idx].begin();iter!=vctLines[idx].end();iter++)
        {
            if(count != 0)
                ostrRet << separator;
            ostrRet << *iter;
            count++;
        }
    }
    return ostrRet.str();
}
string	CParsedTextfile::line(const ParsedTextfile_OneLine& given_iter, const string& separator)
{
    ostringstream ostrReturn;
    for(ParsedTextfile_OneLine::const_iterator iter=given_iter.begin();iter!=given_iter.end();iter++)
        ostrReturn << (iter == given_iter.begin() ? *iter : separator+(*iter));
	
    return ostrReturn.str();
}

string	CParsedTextfile::line(const string& separator)
{
    ostringstream ostrReturn;
    for(ParsedTextfile_OneLine::iterator iter=(*cursor).begin();iter!=(*cursor).end();iter++)
        ostrReturn << *iter << separator;
	
    return ostrReturn.str().substr(0,ostrReturn.str().length()-1);
}

vector<string>	CParsedTextfile::lineVector()
{
    vector<string> vReturn;
    for(ParsedTextfile_OneLine::iterator iter=(*cursor).begin();iter!=(*cursor).end();iter++)
        vReturn.push_back(*iter);
	
    return vReturn;
}

string	CParsedTextfile::GetErrorlog()
{
	string strReturn(errors.str());
	errors.str("");
	return strReturn;	
}

void    CParsedTextfile::add(ParsedTextfile_OneLine val)
{
	if(mapHeader.size()==val.size())
		vctLines.push_back(val);
}
	
	
bool CParsedTextfile::SetHeader(const string& val)
{
	if(mapHeader.size()==0)
            return CreateHeader(val);
        return false;
}



bool CParsedTextfile::GoMilestone()
{
    try
    {
        for(vector<ParsedTextfile_OneLine>::iterator iter=vctLines.begin();iter!=vctLines.end();iter++)
            if(milestone==iter)
            {
                cursor=milestone;
                return true;
            }
        return false;
    }
    catch(...)
    {
            return false;	
	}	
}

std::map<std::string,size_t>    CParsedTextfile::getCustomRowKey(std::vector<int> columns, std::string sep)
{
    map<std::string,size_t> mRet;
    size_t row_count = 0;
    for(vector<ParsedTextfile_OneLine>::iterator i = vctLines.begin(); i != vctLines.end(); i++, row_count++)
    {
        ostringstream key_maker;
        for(int j = 0; j < columns.size(); j++)
        {
            if(i->size() <= columns[j])
                throw CMyException(string("Failed to generate Custom row key because column is out of bounds"));
            key_maker << (*i)[columns[j]];
            if(j+1 < columns.size())
                key_maker <<  sep;
        }
        if(mRet.find(key_maker.str()) != mRet.end())
            throw CMyException(string("Failed to generate Custom row key because the following key is already in use (duplicate unique key?) ")+key_maker.str());
        mRet[key_maker.str()]=row_count;
    }
    return mRet;
}


