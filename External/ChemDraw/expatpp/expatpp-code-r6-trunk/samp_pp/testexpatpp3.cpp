/* 
	\file testexpatpp3.cpp
	Demonstrates use of nesting parsers and parsing attributes
	Structure of topParser is similar to myParser in testexpatpp.cpp
	
	\invoking nested parsers with attributes
	An important point when writing a nested parser is that the element
	used to invoke the parser, eg: <settings>, is "consumed" by the parent parser.
	That's fine if you have an element-oriented architecture like oofRep.
	
	If you expect to have important attributes on the starting element for your
	sub-parser then use a handleAtts() approach as shown below.
	
	If the ONLY thing in your sub-parser is handling attributes and it will never
	be invoked with nested elements then consider not bothering making it a sub-parser
	at all and just make it a method of the current parent parser. 
	
	\par expected output
\verbatim
BEEPX.MSG
- channel
- 1
- msgno
- 1
    ignored
settingsParser::ctor
   ComPort='1'
   DataSource='device'
   RecordFile='off'
   WindowingStyle='hanning'
   Points='256'


settingsParser::ctor


  settingsParser::startElement
     anAtt='Andy'
     anotherAtt='woz here'



finished!
\endverbatim
	
*/

#include <stdio.h>
#include "expatpp.h"
#include <string.h>

/// sample with a nested element with many attributes
const char kSample[] = 
	"<BEEPX.MSG channel='1' msgno='1'>"
		"<ignored />"
		"<settings \n"  // note this is a top-level settings object with attributes
		"ComPort='1'\n"
		"DataSource='device'\n"
		"RecordFile='off'\n"
		"WindowingStyle='hanning'\n"
		"Points='256'\n"
		" />"
		"<settings>"
			"<anElementWithinSettings anAtt='Andy' anotherAtt='woz here' />"
		"</settings>"	
	"</BEEPX.MSG>"
;

class topParser : public expatppNesting {
public:
	topParser() : mDepth(0) {};
	
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
	
private:
	void WriteIndent();
	int	mDepth;	
};


class settingsParser : public expatppNesting {
public:
	settingsParser(expatppNesting* parent, const XML_Char **atts);
	virtual void startElement(const XML_Char* name, const XML_Char** atts);
	
private:
	void handleAtts(const XML_Char** atts);
};
	

void 
topParser::WriteIndent()
{
  for (int i = 0; i < mDepth; i++)
    putchar('\t');
}


void 
topParser::startElement(const char* name, const char** atts)
{
	if(strcmp(name,"settings")==0)
		new settingsParser(this, atts);  // transfer control to settings parser
	else {
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
}


void 
topParser::endElement(const char*)
{
  mDepth--;
}


void 
topParser::charData(const XML_Char *s, int len)
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


settingsParser::settingsParser(expatppNesting* parent, const XML_Char **atts) :
	expatppNesting(parent)
{
	printf("settingsParser::ctor\n");
	handleAtts(atts);
	printf("\n\n");
}


void 
settingsParser::startElement(const XML_Char* name, const XML_Char** atts)
{
	printf("  settingsParser::startElement\n");
	handleAtts(atts);
	printf("\n\n");
}


void 
settingsParser::handleAtts(const XML_Char** atts)
{
	for (int i=0; atts[i]!=0; i++) {
		printf("   %s='%s'\n", atts[i], atts[++i]);
	}
}



int main()
{	
	topParser parser;
	char filename[80];
	int depth = 0;
	
	if (!parser.parseString(kSample)) {
		fprintf(stderr,
			"%s at line %d\n",
			XML_ErrorString(parser.XML_GetErrorCode()),
			parser.XML_GetCurrentLineNumber()
		);
		return 1;
	}
	puts("\nfinished!");
}
