// SimpleTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include "expatpplib.h"
#include <string.h>

class myParser : public expatpp {
public:
	myParser() : mDepth(0) {};
	
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
	
private:
	void WriteIndent();
	int	mDepth;	
};


void 
myParser::WriteIndent()
{
  for (int i = 0; i < mDepth; i++)
    putchar('\t');
}


void 
myParser::startElement(const char* name, const char** atts)
{
  WriteIndent();
  puts(name);
  if (atts) { /* write list of attributes indented below element */
    int i;
    for (i=0; atts[i]; i++) {
      WriteIndent();
      putchar('-'); putchar(' ');
      puts(atts[i]);
    }
  }
  mDepth++;
}


void 
myParser::endElement(const char*)
{
  mDepth--;
}


void 
myParser::charData(const XML_Char *s, int len)
{
  const int leadingSpace = skipWhiteSpace(s);
  if (len==0 || len==leadingSpace)
  	return;  // called with whitespace between elements
  	
  WriteIndent();

/* write out the user data bracketed by ()*/
  putchar('(');
  fwrite(s, len, 1, stdout);
  puts(")");
}


int main(int argc, char* argv[])
{	
	myParser parser;
	char filename[80];
	FILE* xmlFile;
	for (;;) {	  
		int depth = 0;
		
		puts("\n\nXML test: enter filename");
		gets(filename);
		if (strlen(filename)==0)
			break;
		
		xmlFile = fopen(filename, "r");
		if (!xmlFile)
			break;
		
		if (!parser.parseFile(xmlFile)) {
			fprintf(stderr,
				"%s at line %d\n",
				XML_ErrorString(parser.XML_GetErrorCode()),
				parser.XML_GetCurrentLineNumber()
			);
			return 1;
		}
	}  // loop asking for and parsing files
	puts("\nfinished!");
	return 0;
}
