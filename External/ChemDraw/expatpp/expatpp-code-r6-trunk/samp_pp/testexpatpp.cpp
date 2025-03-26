/* This is simple demonstration of how to use expat. This program
reads an XML document from standard input and writes a line with the
name of each element to standard output indenting child elements by
one tab stop more than their parent element. 

Copied from elements.c 98/08/31 AD to open specified file rather than stdin
and added a CharacterDataHandler
*/

#include <stdio.h>
#include "expatpp.h"
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
myParser::startElement(const XML_Char* name, const XML_Char** atts)
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
myParser::endElement(const XML_Char*)
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


int main()
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
}
