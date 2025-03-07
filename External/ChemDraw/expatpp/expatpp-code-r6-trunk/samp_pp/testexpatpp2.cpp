/**
\file testexpatpp2.cpp
Uses OOFILE XML classes to help write an XML file to the name you specify
then parses it back into a data structure. This shows the simple persistence
you would use to store an XML-based preferences file.
*/

#include <stdio.h>
#include "expatpp.h"
#include <string.h>
#ifndef H_OOFXML
	#include "oofxml.h"
#endif 

/**
Simple data class with public members so it can be updated
by parser.

writes itself as XML file like
<?xml version="1.0" standalone="no" ?>
<simpleData>
	<name>Andy Dent</name>
	<isProgrammer/>
	<age>37</age>
</simpleData>
*/


// -------------------------------------------------------
//      s i m p l e D a t a
// -------------------------------------------------------
class simpleData {
public:
	simpleData(const char* inName=0, int inAge=0, bool inProg=false);
	
	static void writeSampleXML(FILE* xmlFile);
	void writeXML(FILE* xmlFile);
	void dumpData();

	oofString	mName;
	int	mAge;
	bool mIsProgrammer;
};


simpleData::simpleData(const char* inName, int inAge, bool inProg) :
	mName(inName),
	mAge(inAge),
	mIsProgrammer(inProg)
{
}


void
simpleData::writeXML(FILE* xmlFile)
{
	oofXMLwriter theXML;
	theXML.startElement("simpleData");
		theXML.addSimpleElement("name", mName);
		theXML.addSimpleElement("age", mAge);
		if (mIsProgrammer)
			theXML.addEmptyElement("isProgrammer");
	theXML.endElement();
	
	assert(theXML.topLevelClosed());
	
	fputs(theXML.generatedXML(), xmlFile);
	fseek(xmlFile, 0, SEEK_SET); // reset for reading
}


void
simpleData::dumpData()
{
	printf("simpleData\n");
	printf("   name=%s\n", mName.chars());
	printf("   age=%d\n", mAge);
	if (mIsProgrammer)
		printf("   is Programmer\n");
	else
		printf("   is NOT Programmer\n");
}


void
simpleData::writeSampleXML(FILE* xmlFile)
{
	simpleData theSample("Andy Dent", 37, true);
	theSample.writeXML(xmlFile);
}


// -------------------------------------------------------
//      s d P a r s e r
// -------------------------------------------------------
class sdParser : public expatpp {
public:
	sdParser(simpleData&);
	
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
	
private:
	enum elementStateT {eNone, eInName, eInAge};
	elementStateT mState;
	simpleData& mData;	
};


sdParser::sdParser(simpleData& inData) : 
	mState(eNone),
	mData(inData) 
{}


/**
Either change state or set bool indicating hit empty element
*/
void 
sdParser::startElement(const XML_Char* name, const XML_Char** /*atts*/)
{
	if(strcmp(name,"name")==0){
		// ignore the attributes
		assert(mState==eNone);  // enforce the structure we expect
		mState = eInName;
	}
	else if(strcmp(name,"age")==0){
		assert(mState==eNone);  // enforce the structure we expect
		mState = eInAge;
	}
	else if(strcmp(name,"isProgrammer")==0) {
		mData.mIsProgrammer = true;
		mState = eNone;
	}
}


void 
sdParser::endElement(const XML_Char*)
{
	mState = eNone;
}


/**
Based on state - read name or age.
If you are developing to a rigid XML format it is a really good idea to start
with a parser that has assertions for unknowns like this one because it will quickly
help you identify problems with the data being written. Flexible XML formats make it hard
to decide if bugs are in the writer or reader.
*/
void 
sdParser::charData(const XML_Char *s, int len)
{
	if (emptyCharData(s, len))
		return;
	
	switch(mState) {
		case eInName :
			mData.mName.setChars(s, len);
			break;
			
		case eInAge :
		{
			int i;
			std::sscanf(s,"%d", &i);
			mData.mAge = i;
		}
			break;
		
		default : {
			assert(mState==eNone);  // more flexible parser would ignore unknown states
		}	
	};
}


// -------------------------------------------------------
//      m a i n
// -------------------------------------------------------
/**
sample showing parsing of file in loop reading minimal sized buffer
to cut down on overhead.
*/
int main()
{
	simpleData theData;	
	sdParser parser(theData);
	char filename[80];
	FILE* xmlFile;
	char buf[BUFSIZ];
	
	int done;
	int depth = 0;
	
	puts("\n\nXML test: enter filename to write");
	gets(filename);
	if (strlen(filename)==0)
		return 0;
	
	xmlFile = fopen(filename, "wb+");
	if (!xmlFile)
		return 0;
	
	puts("sample data in original empty state!");
	theData.dumpData();
	
	simpleData::writeSampleXML(xmlFile);
  
// read back thefile
	if (!parser.parseFile(xmlFile)) {
		fprintf(stderr,
			"%s at line %d\n",
			XML_ErrorString(parser.XML_GetErrorCode()),
			parser.XML_GetCurrentLineNumber()
		);
		return 1;
	}

	puts("\n\nsample data after reading back!\n");
	theData.dumpData();
	
	puts("\nfinished!");
}

