// expatpp
#ifdef UNDER_CE
	#include <string.h>
	#include <windows.h>
	#include <dbgapi.h>
	#define assert ASSERT
#else
	#include <string>
        #include <string.h>
	using namespace std;
	#include <assert.h>
#endif
#include "expatpp.h"

	
// may be defined in xmltchar.h or elsewhere
#ifndef tcscmp
	#ifdef XML_UNICODE
		#define tcscmp wcscmp
	#else
		#define tcscmp strcmp
	#endif  // XML_UNICODE
#endif  // tcscmp


#ifndef BUFSIZ
	#define BUFSIZ 4096
#endif

expatpp::expatpp(bool createParser) :
	mParser(0),  // in case of exception below
	mHaveParsed(false)
{
  if (createParser) {
  // subclasses may call this ctor after parser created!
		mParser = XML_ParserCreate(0);
		SetupHandlers();
	}
}


void
expatpp::SetupHandlers()
{
	::XML_SetUserData(mParser, this);
	::XML_SetElementHandler(mParser, startElementCallback, endElementCallback);
	::XML_SetCharacterDataHandler(mParser, charDataCallback);
	::XML_SetProcessingInstructionHandler(mParser, processingInstructionCallback);
	::XML_SetDefaultHandler(mParser, defaultHandlerCallback);
	::XML_SetUnparsedEntityDeclHandler(mParser, unParsedEntityDeclCallback);
	::XML_SetNotationDeclHandler(mParser, notationDeclCallback);
	::XML_SetNotStandaloneHandler(mParser, notStandaloneHandlerCallback);
	::XML_SetNamespaceDeclHandler(mParser, startNamespaceCallback, endNamespaceCallback);
#ifndef EXPATPP_COMPATIBLE_EXPAT12
	::XML_SetAttlistDeclHandler(mParser, attlistDeclCallback);
	::XML_SetCdataSectionHandler(mParser, startCdataSectionCallback, endCdataSectionCallback);
	::XML_SetCommentHandler(mParser, commentCallback);
	::XML_SetDoctypeDeclHandler(mParser, startDoctypeDeclCallback, endDoctypeDeclCallback);
	::XML_SetElementDeclHandler(mParser, elementDeclCallback);
	::XML_SetEntityDeclHandler(mParser, entityDeclCallback);
	::XML_SetSkippedEntityHandler(mParser, skippedEntityCallback);
	::XML_SetXmlDeclHandler(mParser, xmlDeclCallback);		  
#endif
}


expatpp::~expatpp()
{
	if (mParser)  // allows subclasses to avoid finishing parsing
	  ReleaseParser();
}


/**
	Provide single point that will call XML_ParserFree.
	Nothing else in this code should call XML_ParserFree!
*/
void 
expatpp::ReleaseParser()
{
	::XML_ParserFree(mParser);
	mParser = 0;
}



/**
	Provide single point that will call XML_ParserReset.
	Guarded against trivial reset before use in case that breaks
	expat or creates overhead.
	
	\todo pass in encoding to XML_ParserReset when we support encodings
*/
void 
expatpp::ResetParser()
{
#ifdef EXPATPP_COMPATIBLE_EXPAT12
	assert(!"Reset not available in earlier than expat 1.95.3");s
#else
	if (mHaveParsed) {
		::XML_ParserReset(mParser, NULL);
		SetupHandlers();
		mHaveParsed = false;
	}
#endif
}


/**
	Parse entire file, basically copy of the loop from the elements.c example.
*/
XML_Status
expatpp::parseFile(FILE* inFile)
{	
	ResetParser();
	
	char buf[BUFSIZ];
	int done;
	if (!inFile)
	  	return XML_STATUS_ERROR;
	fseek(inFile, 0, SEEK_SET); // reset for reading
	do {
		size_t len = fread(buf, 1, sizeof(buf), inFile);
		done = len < sizeof(buf);
		enum XML_Status parseStatus;
		if ((parseStatus = XML_Parse(buf, len, done))!=XML_STATUS_OK) {
			return parseStatus;
		}
	} while (!done);
	return XML_STATUS_OK;
}


XML_Status
expatpp::XML_Parse(const char *s, int len, int isFinal)
{
	mHaveParsed = true;
	const XML_Status retStatus = ::XML_Parse(mParser, s, len, isFinal);
	if (isFinal)
		CheckFinalStatus(retStatus);
	return retStatus;
}


XML_Error
expatpp::XML_GetErrorCode()
{
	return ::XML_GetErrorCode(mParser);
}


int
expatpp::XML_GetCurrentLineNumber()
{
	return ::XML_GetCurrentLineNumber(mParser);
}


int
expatpp::XML_GetCurrentColumnNumber()
{
	return ::XML_GetCurrentColumnNumber(mParser);
}




/**
	Parse string which is assumed to be entire XML document.
	Written to stop stupid errors of being off by one in the string length causing 
	wasted debugging time, such as:
\verbatim	
	const char[] kSampleSettings = "<settings/>";
	const int sampleSize = sizeof(kSampleSettings)-1;  // unless you remember to subtract one here will get invalid token error
	if (!parser.XML_Parse(kSampleSettings, sampleSize, 1)) {
\endverbatim	
*/
XML_Status
expatpp::parseString(const char* inString)
{
	ResetParser();
	const int inLen = strlen(inString);
	return XML_Parse(inString, inLen, 1);	
}

void 
expatpp::startElementCallback(void *userData, const XML_Char* name, const XML_Char** atts)
{
	((expatpp*)userData)->startElement(name, atts);
}


void 
expatpp::endElementCallback(void *userData, const XML_Char* name)
{
	((expatpp*)userData)->endElement(name);
}


void 
expatpp::startNamespaceCallback(void *userData, const XML_Char* prefix, const XML_Char* uri)
{
	((expatpp*)userData)->startNamespace(prefix, uri);
}


void 
expatpp::endNamespaceCallback(void *userData, const XML_Char* prefix)
{
	((expatpp*)userData)->endNamespace(prefix);
}


void 
expatpp::charDataCallback(void *userData, const XML_Char* s, int len)
{
	((expatpp*)userData)->charData(s, len);
}


void
expatpp:: processingInstructionCallback(void *userData, const XML_Char* target, const XML_Char* data)
{
	((expatpp*)userData)->processingInstruction(target, data);
}


void
expatpp::defaultHandlerCallback(void* userData, const XML_Char* s, int len)
{
	((expatpp*)userData)->defaultHandler(s, len);
}


int
expatpp::notStandaloneHandlerCallback(void* userData)
{
	return ((expatpp*)userData)->notStandaloneHandler();
}


void
expatpp::unParsedEntityDeclCallback(void* userData, const XML_Char* entityName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId, const XML_Char* notationName)
{
	((expatpp*)userData)->unparsedEntityDecl(entityName, base, systemId, publicId, notationName);
}


void
expatpp::notationDeclCallback(void *userData, const XML_Char* notationName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId)
{
	((expatpp*)userData)->notationDecl(notationName, base, systemId, publicId);
}


void 
expatpp::startElement(const XML_Char*, const XML_Char**)
{}


void 
expatpp::endElement(const XML_Char*)
{}


void 
expatpp::startNamespace(const XML_Char* /* prefix */, const XML_Char* /* uri */)
{}


void 
expatpp::endNamespace(const XML_Char*)
{}


void 
expatpp::charData(const XML_Char*, int )
{
}


void
expatpp::processingInstruction(const XML_Char*, const XML_Char*)
{
}


void
expatpp::defaultHandler(const XML_Char*, int)
{
}


int
expatpp::notStandaloneHandler()
{
	return 0;
}


void
expatpp::unparsedEntityDecl(const XML_Char*, const XML_Char*, const XML_Char*, const XML_Char*, const XML_Char*)
{
}


void
expatpp::notationDecl(const XML_Char*, const XML_Char*, const XML_Char*, const XML_Char*)
{
}


int 
expatpp::skipWhiteSpace(const XML_Char* startFrom)
{
	// use our own XML definition of white space
	// TO DO - confirm this is correct!
	const XML_Char* s = startFrom;
	XML_Char c = *s;
	while ((c==' ') || (c=='\t') || (c=='\n') || (c=='\r')) {
		s++;
		c = *s;
	}
	const int numSkipped = s - startFrom;
	return numSkipped;
}


/**
	Iterate the paired attribute name/value until find a pair with matching name.
	\return pointer to the value or null if not found.
*/
const XML_Char* 
expatpp::getAttribute(const XML_Char* matchingName, const XML_Char** atts)
{
	for (int i=0; atts[i]; i++) {
		const XML_Char* attributeName = atts[i++];
		assert(attributeName);  // shouldn't fail this because of loop test above
		if(tcscmp(attributeName, matchingName)==0) {  
			return atts[i];  // if 2nd item was missing, this returns 0 safely indicating failure
		}
	}
	return 0;
}


/**
\bug will always return 0 for PPC
*/
bool 
expatpp::getIntegerAttribute(const XML_Char *matchingName, const XML_Char **atts, int& outAtt)
{
	const XML_Char* attStr = getAttribute(matchingName, atts);
	if (!attStr)
		return false;
	int i=0;
#ifdef XML_UNICODE
fail to compile because need this now
#else
	sscanf(attStr, "%d", &i);
#endif
	outAtt = i;
	return true;
}


/**
\bug will always return 0 for PPC
*/
bool 
expatpp::getDoubleAttribute(const XML_Char *matchingName, const XML_Char **atts, double& outAtt)
{
	const XML_Char* attStr = getAttribute(matchingName, atts);
	if (!attStr)
		return false;
	float f = 0.0;  // sscanf doesn't allow point to double
#ifdef XML_UNICODE
fail to compile because need this now
#else
	sscanf(attStr, "%f", &f);
#endif
	outAtt = f;
	return true;
}


bool 
expatpp::emptyCharData(const XML_Char *s, int len)
{
// usually call from top of overriden charData methods
	if (len==0)
		return true;  //*** early exit - empty string, may never occur??
		
// skip newline and empty whitespace
	if (
		((len==1) && ( (s[0]=='\n') || (s[0]=='\r')) ) ||  // just CR or just LF
		((len==2) && (s[0]=='\r') && (s[1]=='\n'))  // DOS-style CRLF
	)
		return true;  //*** early exit - newline
		
	const int lastCharAt = len-1;
	if (s[lastCharAt]==' ') {  // maybe all whitespace
		int i;
		for (i=0; i<lastCharAt; i++) {
			if (s[i]!=' ')
				break;
		}
		if (i==lastCharAt)
			return true;	  //*** early exit - all spaces
	}
	return false;
}


//-------- Added for expat 1.95.5---------------
void
expatpp::attlistDeclCallback(void *userData, 
	const XML_Char *elname,
	const XML_Char *attname,
	const XML_Char *att_type,
	const XML_Char *dflt,
	int             isrequired)
{
	((expatpp*)userData)->attlistDecl(elname, attname, att_type, dflt, isrequired);
}


void
expatpp::commentCallback(void *userData, const XML_Char *data)
{
	((expatpp*)userData)->comment(data);
}


void
expatpp::elementDeclCallback(void *userData, const XML_Char *name, XML_Content *model)
{
	((expatpp*)userData)->elementDecl(name, model);
}


void
expatpp::endCdataSectionCallback(void *userData)
{
	((expatpp*)userData)->endCdataSection();
}


void
expatpp::endDoctypeDeclCallback(void *userData)
{
	((expatpp*)userData)->endDoctypeDecl();
}


void
expatpp::entityDeclCallback(void *userData,
	const XML_Char *entityName,
	int is_parameter_entity,
	const XML_Char *value,
	int value_length,
	const XML_Char *base,
	const XML_Char *systemId,
	const XML_Char *publicId,
	const XML_Char *notationName)
{
	((expatpp*)userData)->entityDecl(entityName, is_parameter_entity, value, value_length, base, systemId, publicId, notationName);
}


void
expatpp::skippedEntityCallback(void *userData, const XML_Char *entityName, int is_parameter_entity)
{
	((expatpp*)userData)->skippedEntity(entityName, is_parameter_entity);
}


void
expatpp::startCdataSectionCallback(void *userData)
{
	((expatpp*)userData)->startCdataSection();
}


void
expatpp::startDoctypeDeclCallback(void *userData, 
		const XML_Char *doctypeName,
        const XML_Char *sysid,
        const XML_Char *pubid,
        int has_internal_subset)
{
	((expatpp*)userData)->startDoctypeDecl(doctypeName, sysid, pubid, has_internal_subset);
}


void
expatpp::xmlDeclCallback(void *userData, const XML_Char      *version,
                                    const XML_Char      *encoding,
                                    int                  standalone)
{
	((expatpp*)userData)->xmlDecl(version, encoding, standalone);
}


void
expatpp::attlistDecl( 
	const XML_Char *elname,
	const XML_Char *attname,
	const XML_Char *att_type,
	const XML_Char *dflt,
	int             isrequired)
{
}


void
expatpp::comment( const XML_Char *data)
{
}


void
expatpp::elementDecl( const XML_Char *name, XML_Content *model)
{
}


void 
expatpp::endCdataSection()
{
}


void
expatpp::endDoctypeDecl()
{
}


void
expatpp::entityDecl(
	const XML_Char *entityName,
	int is_parameter_entity,
	const XML_Char *value,
	int value_length,
	const XML_Char *base,
	const XML_Char *systemId,
	const XML_Char *publicId,
	const XML_Char *notationName)
{
}


void
expatpp::skippedEntity( const XML_Char *entityName, int is_parameter_entity)
{
}


void 
expatpp::startCdataSection()
{
}


void
expatpp::startDoctypeDecl(const XML_Char *doctypeName,
        const XML_Char *sysid,
        const XML_Char *pubid,
        int has_internal_subset)
{
}


void
expatpp::xmlDecl( const XML_Char      *version,
                                    const XML_Char      *encoding,
                                    int                  standalone)
{
}




// -------------------------------------------------------
//      e x p a t p p N e s t i n g
// -------------------------------------------------------
/**
	\param parent can be null in which case this is root parser
	
	\note The handlers set in here MUST be also set in SetupHandlers
	which is a virtual method invoked by expatpp::ResetParser. Otherwise
	you can have subtle bugs with a nested parser not properly returning
	after reusing a parser (nasty and found rapidly only via extensive unit
	tests and plentiful assertions!).
	
	\WARNING 
	The assumption that is not obvious here is that if you want to use 
	nested parsers, then your topmost parser must also be an expatppNesting
	subclass, NOT an expatpp subclass, because we need the 
	nestedStartElementCallback and nestedEndElementCallback
	callbacks to override those in the expatpp ctor.
	
	
	
	\todo go back over code in detail and confirm above warning still valid
	I think if we used expat's functions to invoke the registered callback
	might be safer - the explicit function call we have in nestedEndElementCallback
	certainly assumes the parent type.
*/
expatppNesting::expatppNesting(expatppNesting* parent) :
	expatpp(parent==0),  // don't create parser - we're taking over from parent if given
	mDepth(0),
	mParent(parent),
	mOwnedChild(0),
	mSelfDeleting(true)
{
	if ( parent )
	{
		RegisterWithParentXMLParser();
		parent->AdoptChild(this);
	}
	else
	{
		// No parent - the expatpp constructor will have created a new mParser (expat parser)
		::XML_SetElementHandler(mParser, nestedStartElementCallback, nestedEndElementCallback);
	}
	assert(mParser);  // either we created above or expatpp 
}


expatppNesting::~expatppNesting()
{
	assert(!mParent);  // if we are a sub-parser, should not delete without calling returnToParent
	DeleteChild();
}


/**
	Call parent version then override same as in our ctor.
*/
void
expatppNesting::SetupHandlers()
{
	expatpp::SetupHandlers();
	::XML_SetElementHandler(mParser, nestedStartElementCallback, nestedEndElementCallback);
}

/**
	Must use if you have adopted a child parser and want to dispose of it early.
*/
void
expatppNesting::DeleteChild()
{
	delete mOwnedChild;
	mOwnedChild = 0;
}


/**
	Invoked as a callback from a child ctor when we pass in a parent pointer.
	OR used from switchToNewSubParser, in which case it may be the 2nd time
	we're called for a given child (see scenarios in expatppNesting class comment).
*/
void
expatppNesting::AdoptChild(expatppNesting* adoptingChild)
{
	if ( mOwnedChild != adoptingChild )
	{
		delete mOwnedChild;
		mOwnedChild = adoptingChild;
	}
}


/**
	 to use parent's underlying expat parser
*/
void
expatppNesting::RegisterWithParentXMLParser()
{
	mParser = mParent->mParser;
	::XML_SetUserData(mParser, this);
}


/**
	User code (typically the startElement handler of user parsers derived from expatppNesting) 
	may call 
		switchToNewSubParser( new UserChildParser() );
	to hand off the current document to a child parser that understands the next segment of XML.
	Control will be returned to the original (parent) parser when the end of the child element 
	is reached.
	In its lifetime a 'parent' parser may switch control to several child parsers (one at a time 
	of course) as it moves through the document encoutering various types of child element.
	
	A child to which older code (eg: OOFILE) has just switched control by
	new childParser(this) will be self-deleting and will clear our mOwnedChild in its dtor. 
*/
void expatppNesting::switchToNewSubParser( expatppNesting* pAdoptedChild )
{
	assert(pAdoptedChild);
	AdoptChild(pAdoptedChild);
	pAdoptedChild->BeAdopted(this);
}


/**
	If this is root parser, nestedEndElementCallback won't call returnToParent.
	Therefore it is safe to put parsers on the stack.
*/
expatppNesting* 
expatppNesting::returnToParent()
{
	expatppNesting* ret = mParent;
	::XML_SetUserData(mParser, mParent);
	mParent=0;
	mParser=0;  // prevent parser shutdown by expatpp::~expatpp!!
	if (mSelfDeleting) {
		ret->OwnedChildOrphansItself(this);
		delete this;  // MUST BE LAST THING CALLED IN NON-VIRTUAL FUNCTION, NO MEMBER ACCESS
	}
	return ret;
}


void 
expatppNesting::nestedStartElementCallback(void *userData, const XML_Char* name, const XML_Char** atts)
{
	assert(userData);
	expatppNesting* nestedParser = (expatppNesting*)userData;
	nestedParser->mDepth++;
	nestedParser->startElement(name, atts);  // probably user override
}


/**
	If this is root parser, will never hit nestedEndElementCallback after closing element,
	except for when we call it.
	\param userData should be non-nil except for specific case of ending root
*/
void 
expatppNesting::nestedEndElementCallback(void *userData, const XML_Char* name)
{
	if (!userData)
		return;  //  end tag for root
		
	expatppNesting* nestedParser = (expatppNesting*)userData;
// we don't know until we hit a closing tag 'outside' us that our run is done 	
	if (nestedParser->mDepth==0) {
		expatppNesting* parentParser = nestedParser->returnToParent();
		nestedEndElementCallback(parentParser, name);   // callbacks for expatppNesting stay registered, so safe 
		//if we don't invoke their callback, they will not balance their mDepth		
	}
	else {
	// end of an element this parser has started - normal case
		nestedParser->endElement(name);  // probably user override
		nestedParser->mDepth--;
	}
}


/**
	Called by switchToNewSubParser to indicate a newly created child parser
	is now the currently active child for adoptingParent and the child
	isn't expected to be self deleting.
	
	Normal code to create an owned child would be either
		switchToNewSubParser( new UserChildParser(this) );
	where this is the currently active parser and you want to be deleting it, or
		new UserChildParser(this);
	to have a child parser self-delete
	
	\par Important Safety Note
	Copes with the situation of people forgetting to pass 
	in the parent parser (and hence creating a new one by default)
	if invoked by switchToNewSubParser( new UserChildParser() )
	by somewhat wastefully deleting the parser created in expatpp::expatpp
	by us being a root parser.	
*/
void
expatppNesting::BeAdopted(expatppNesting* adoptingParent)
{
	if (mParent) {
		assert(mParent==adoptingParent);
	}
	else {  // root parser being adopted, cleanup!
		ReleaseParser();
		mParent = adoptingParent;
		RegisterWithParentXMLParser();
	}
	mSelfDeleting = false;
}




