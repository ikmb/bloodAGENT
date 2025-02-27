#ifndef MYTOOLS_H_
#define MYTOOLS_H_

using namespace std;

#define MYTOOLS_H_MAX_ITOA_DIGITS_SIZE 35

class CMyTools
{
public:
	CMyTools();
	virtual ~CMyTools();

	// String manipulation	
	static vector<string> GetParsedLine(string,const string& sep="\t");
	static vector<string> Tokenize(string a,const string& sep="\t"){return GetParsedLine(a,sep.c_str());}
	static string	      ReplaceSubstr(std::string& result, const std::string& replaceWhat,  const std::string& replaceWithWhat);
        static string         dtoa(long value, int base);
        static void           erase_windows_line_ending(string&);
	// File System
	static vector<string> GetDirectoryListing(const string& path);
        static vector<string> GetDirectoryListing(const string& path,const char* contains, bool full_path = false);
	static int cmd(const string& commando, bool write_output);
        static int cmd(const string& commando, std::string& get_output);
        static bool file_exists(const std::string& fileName);
	static bool isNumeric(const char *);
	
	// genetic code
	static string GetComplSequence(const string& strSequence);
	static char   GetComplNuc(const char& source);
        static std::pair<std::string, char> getAminoAcid(std::string triplet);
	static std::vector<unsigned int> FindApproxSubs(unsigned int k, string & P, string & T);

        template<class T> static void GetMem(T*& OldMem, int Elems);
        
	// AB color code
	static string getBasespace(string cs,bool keep_lead_base=true);
	static string getBasespaceRev(string cs,bool keep_lead_base=true);
	static char   decodeBaseChange(char base, char cc);
        
        // diverse
        static bool doIntersect(int a_start,int a_end,int b_start,int b_end);
        static size_t intersectCount(int a_start,int a_end,int b_start,int b_end);
        static size_t insertSize(int a_start,int a_end,int b_start,int b_end);
        static float bytesToFloat(unsigned char* byte_array);
        
        static string toString(std::set<std::string> v, const char sep = ',');
        
        static std::string trim(const std::string &s);
        
private:

};

#endif /*MYTOOLS_H_*/


