#ifndef ParsedTextfile_H_
#define ParsedTextfile_H_

#include <map>
#include <vector>


typedef std::vector<std::string> ParsedTextfile_OneLine;


class CParsedTextfile : private CMyTools
{
    
        
public:
    typedef vector<ParsedTextfile_OneLine>::const_iterator const_iterator;
    typedef vector<ParsedTextfile_OneLine>::iterator             iterator;
    
	CParsedTextfile(); 
	CParsedTextfile(const CParsedTextfile&); 
	CParsedTextfile(const string&,const string& sep = "\t",const string& row_key = "", int skip = 0, bool read_header = true, const string& comment = "#");
	CParsedTextfile(const string&,const string& sep = "\t",int row_key = -1          , int skip = 0, bool read_header = true, const string& comment = "#");
	CParsedTextfile(const string&,const string& sep = "\t",const string& row_key = "", const string& read_until = "", bool read_header = true, const string& comment = "#");
	CParsedTextfile(ifstream&,const string& sep = "\t",const string& row_key = "", int skip = 0, bool read_header = true, const string& comment = "#");
	CParsedTextfile(ifstream&,const string& sep = "\t",const string& row_key = "", const string& read_until = "", bool read_header = true, const string& comment = "#");
	CParsedTextfile(ifstream&,const string& sep = "\t",int row_key = -1, int skip = 0, bool read_header = true, const string& comment = "#");
	
	virtual ~CParsedTextfile();
		
	string operator[](const string&);
        string operator[](const int idx);
        string cell(const string& rowKey,const string& colKey)const;
        string cell(const string& rowKey,int col)const;
        string cell(int row,const string& colKey)const;
        string cell(int row,int col)const;
        /// get cell tokenized
        vector<string> cellTok(const string& rowKey,const string& colKey, const string& sep = ";");
        /// get cell tokenized
        vector<string> cellTok(const string& rowKey,int col, const string& sep = ";");
        /// get cell tokenized
        vector<string> cellTok(int row,const string& colKey, const string& sep = ";");
        /// get cell tokenized
        vector<string> cellTok(int row,int col, const string& sep = ";");
        iterator find(const string& rowKey);
	CParsedTextfile& operator=(const CParsedTextfile&);
	
	void    add(ParsedTextfile_OneLine val);
        vector<string> getHeader();
        string  getHeaderAsString(const string& sep = "\t");
	bool 	SetHeader(const string& val);
        void    CoutHeader();
	
	string	line();
	string	line(int idx, const string& separator = "\t");
	string	line(const string& );
        std::vector<std::string> lineVector();
	string	line(const ParsedTextfile_OneLine& given_iter, const string& separator = "\t");
	ParsedTextfile_OneLine  lineObject(){return *cursor;}
	int GetColId(const string& );
	size_t size(){return vctLines.size();}
	size_t lineSize(){if(cursor != vctLines.end())return cursor->size();return string::npos;}
	
	string	GetErrorlog();
	bool 	hasErrors(){if(errors.str().size()!=0)return true;return false;}
        
        bool hasRowKey(const std::string& key)const{return (mapRow.find(key) != mapRow.end());}
	bool hasColKey(const std::string& key)const{return (mapHeader.find(key) != mapHeader.end());}
	
	bool Next();
	bool Prev();
	bool First();
	bool Last();
	bool isEOF();
        bool erase(){return erase(cursor);};
        bool erase(CParsedTextfile::iterator del_this);
        vector<ParsedTextfile_OneLine>::iterator getCursor(){return cursor;}
        
        /// 
        /// \param columns, select columns that should be combined to create a unique key
        /// \param sep, the separator to separate the column entries for key generation
        /// \return the map with the key first and the corresponding row second
        std::map<std::string,size_t>    getCustomRowKey(std::vector<int> columns, std::string sep="_");
        
        
	vector<ParsedTextfile_OneLine>::iterator begin(){return vctLines.begin();}
	vector<ParsedTextfile_OneLine>::iterator end(){return vctLines.end();}
	
	void SetMilestone(){milestone=cursor;}
	bool GoMilestone();
	
	void clear(){vctLines.clear();cursor=milestone=vctLines.end();}
        
        size_t rowCount()const{return vctLines.size();}
        size_t colCount()const{return cursor->size();}
        
        string getActSeparator(){return m_strValidSeparator;}
        
	
private:
	
	map<string,int>                                 mapHeader;
	map<string,size_t>                              mapRow;
	vector<ParsedTextfile_OneLine>		 	vctLines;
	vector<ParsedTextfile_OneLine>::iterator	cursor;
	vector<ParsedTextfile_OneLine>::iterator	milestone;
	ostringstream					errors;
        
        const string& m_comment_char;
	
	string	m_strValidSeparator;
        string  m_row_key;
        int     m_row_key_column;
        int     m_skip;
        
	bool CreateHeader(const string&);
	bool CreateRowMap();
	bool ReadFile(ifstream& file, bool read_header = true, std::string read_until="");
	
	
};

#endif /*ParsedTextfile_H_*/


