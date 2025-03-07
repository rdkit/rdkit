// expatpp
#ifndef H_EXPATPP
#define H_EXPATPP

#ifdef EXPATPP_COMPATIBLE_EXPAT12 // earlier versions of expat up to v1.2 
	#include "xmlparse.h"
#else
	#include "expat.h"  // since some version of expat moved to SourceForge
#endif
#include <stdio.h>
#include <assert.h>


/**
\file expatpp.h
Latest version 29-Dec-2002 compatible with expat 1.95.6
*/

/**
expatpp follows a simple pattern for converting the semi-OOP callback design of 
expat into a true class which allows you to override virtual methods to supply
callbacks.

\par USING expatpp
see testexpatpp.cpp for a detailed example

1) decide which callbacks you wish to use, eg: just startElement

2) declare a subclass of expatpp, eg:
class myExpat : public expatpp {
	virtual void startElement(const XML_Char* name, const XML_Char** atts);
};

3) create an instance of your object and pass in a buffer to parse
myExpat parser;
parser.XML_Parse(buf, len, done)


\par HOW IT WORKS
The User Data which expat maintains is simply a pointer to an instance of your object.

Inline static functions are specified as the callbacks to expat.
These static functions take the user data parameter returned from expat and cast it
to a pointer to an expatpp object.

Using that typed pointer they then call the appropriate virtual method.

If you have overriden a given virtual method then your version will be called, otherwise
the (empty) method in the base expatpp class is called.

\par Possible Efficiency Tactic
For efficiency, you could provide your own constructor and set some of the callbacks
to 0, so expat doesn't call the static functions. (untested idea).

\par Naming Conventions
The virtual functions violate the usual AD Software convention of lowercase first letter 
for public methods but this was a late change to protected and too much user code out there.


\todo Possibly implement some handling for XML_SetExternalEntityRefHandler which does NOT
receive user data, just the parser, so can't use normal pattern for invoking virtual methods

\todo Possibly implement handling for XML_UnknownEncodingHandler.

\todo review design for nested calls - not happy that it is the right thing that they don't see
their start and ending elements - makes it harder to unit test them in isolation.

\todo unit tests

\todo especially test abort mechanism

\todo reinstate copy constrution and assignment with child parser cleanup

\todo allow specification of encoding
*/
class expatpp {
public:
	expatpp(bool createParser=true);
	virtual ~expatpp();

	operator XML_Parser() const;
	
protected:  // callback virtuals should only be invoked through our Callback static functions
	bool emptyCharData(const XML_Char* s, int len);  // utility often used in overridden charData

// overrideable callbacks
	virtual void startElement(const XML_Char* name, const XML_Char** atts);
	virtual void endElement(const XML_Char*);
	virtual void charData(const XML_Char*, int len);
	virtual void processingInstruction(const XML_Char* target, const XML_Char* data);
	virtual void defaultHandler(const XML_Char*, int len);
	virtual int notStandaloneHandler();
	virtual void unparsedEntityDecl(const XML_Char* entityName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId, const XML_Char* notationName);
	virtual void notationDecl(const XML_Char* notationName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId);
	virtual void startNamespace(const XML_Char* prefix, const XML_Char* uri);
	virtual void endNamespace(const XML_Char*);
/// \name Callbacks added to support expat 1.95.5
//@{
	virtual void attlistDecl( 
		const XML_Char *elname,
		const XML_Char *attname,
		const XML_Char *att_type,
		const XML_Char *dflt,
		int             isrequired);	
	virtual void endCdataSection();	
	virtual void endDoctypeDecl();	
	virtual void comment( const XML_Char *data);	
	virtual void elementDecl( const XML_Char *name, XML_Content *model);	
	virtual void entityDecl(
		const XML_Char *entityName,
		int is_parameter_entity,
		const XML_Char *value,
		int value_length,
		const XML_Char *base,
		const XML_Char *systemId,
		const XML_Char *publicId,
		const XML_Char *notationName);	
	virtual void skippedEntity(const XML_Char *entityName, int is_parameter_entity);	
	virtual void startCdataSection();
	virtual void startDoctypeDecl(const XML_Char *doctypeName,
	    const XML_Char *sysid,
	    const XML_Char *pubid,
	    int has_internal_subset);	
	virtual void xmlDecl( const XML_Char      *version,
	                                    const XML_Char      *encoding,
	                                    int                  standalone);
//@}

public:	
/// \name XML interfaces
//@{
	XML_Status XML_Parse(const char* buffer, int len, int isFinal);
	virtual XML_Status  parseFile(FILE* inFile);
	virtual XML_Status  parseString(const char*);
	XML_Error XML_GetErrorCode();
	int XML_GetCurrentLineNumber();
	int XML_GetCurrentColumnNumber();
//@}
	
protected:
	XML_Parser mParser;
	bool mHaveParsed;
	
/// \name overrideables to customise behaviour, must call parent
//@{
	virtual void ReleaseParser();
	virtual void ResetParser();
	virtual void SetupHandlers();
//@}

/**
	Override so subclass can react to an error causing exit from parse.
	rather than leave it for application code to check status.
	Useful point to insert logging to silently grab failed parses
*/
	virtual void CheckFinalStatus(XML_Status) {};
		
// static interface functions for callbacks
public:
	static void startElementCallback(void *userData, const XML_Char* name, const XML_Char** atts);
	static void endElementCallback(void *userData, const XML_Char* name);
	static void startNamespaceCallback(void *userData, const XML_Char* prefix, const XML_Char* uri);
	static void endNamespaceCallback(void *userData, const XML_Char* prefix);
	static void charDataCallback(void *userData, const XML_Char* s, int len);
	static void processingInstructionCallback(void *userData, const XML_Char* target, const XML_Char* data);
	static void defaultHandlerCallback(void* userData, const XML_Char* s, int len);
	static int notStandaloneHandlerCallback(void* userData);	
	static void unParsedEntityDeclCallback(void* userData, const XML_Char* entityName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId, const XML_Char* notationName);
	static void notationDeclCallback(void *userData, const XML_Char* notationName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId);
/// \name Callback interfacess added to support expat 1.95.5
//@{
	static void attlistDeclCallback(void *userData,  
		const XML_Char *elname,
		const XML_Char *attname,
		const XML_Char *att_type,
		const XML_Char *dflt,
		int             isrequired);	
	static void commentCallback(void *userData,  const XML_Char *data);	
	static void elementDeclCallback(void *userData,  const XML_Char *name, XML_Content *model);	
	static void endCdataSectionCallback(void *userData);
	static void endDoctypeDeclCallback(void *userData);	
	static void entityDeclCallback(void *userData, 
		const XML_Char *entityName,
		int is_parameter_entity,
		const XML_Char *value,
		int value_length,
		const XML_Char *base,
		const XML_Char *systemId,
		const XML_Char *publicId,
		const XML_Char *notationName);	
	static void skippedEntityCallback(void *userData,  const XML_Char *entityName, int is_parameter_entity);	
	static void startCdataSectionCallback(void *userData);
	static void startDoctypeDeclCallback(void *userData, 
		const XML_Char *doctypeName,
        const XML_Char *sysid,
        const XML_Char *pubid,
        int has_internal_subset);	
	static void xmlDeclCallback(void *userData,  const XML_Char      *version,
	                                    const XML_Char      *encoding,
	                                    int                  standalone);
//@}

	
// utilities
	static int skipWhiteSpace(const XML_Char*);
	static const XML_Char* getAttribute(const XML_Char *matchingName, const XML_Char **atts);
	static bool getIntegerAttribute(const XML_Char *matchingName, const XML_Char **atts, int& outAtt);
	static bool getDoubleAttribute(const XML_Char *matchingName, const XML_Char **atts, double& outAtt);
};


/**
	subclass to support a hierarchy of parsers, in a sort of recursion or
	'nesting' approach, where a top-level parser might create sub-parsers 
	for part of a file.
	
	The currently active child parser is owned (mOwnedChild) and is deleted
	by DeleteChild (invoked from the dtor) so error handling can propagate 
	up the tree, closing parsers, without leaks.

	\par Switching to sub-parsers
	You can transfer to a sub-parser with
	- new UserChildParser(this)  // carries on using our parser, is self-deleting
	- switchToNewSubParser( someVar = new UserChildParser(this) )  // if want to get values back after end parsing

	\warning You can accidentally invoke a new parser without it doing anything
	- new UserChildParser()  // will be new top-level parser, nothing to do with our XML
	
	\par Self-deletion
	If you transfer control to a sub-parser with just new UserChildParser(this) then
	it will be automatically self-deleting in its returnToParent method and
	will invoke OwnedChildOrphansItself to clear our mOwnedChild.
	
	The reason for self-deletion being governed by a somewhat complex chain of
	calls rather than simply a boolean flag is because expatpp has been in use
	worldwide for many years and it was deemed too unfriendly to break code in
	a manner which could cause unwanted side effects - the current approach safely
	preserves self-deletion but also allows for expatpp to have parent parsers
	own and delete children, without compiling with different options.
	
	\note 
	If you invoke a sub-parser with switchToNewSubParser( new UserChildParser() );
	then the user child parser will start with a new XML parser instance
	created by the expatpp ctor. This is safe but slightly wasteful of processing 
	as the new parser will be discarded by BeAdopted().

	\par Switching to child	and explicitly deleting
	switchToNewSubParser( somevar = new UserChildParser(this) ) allows you to get values
	back out of the child parser, in the context of the parent, eg:
	
\verbatim

void MultiFilterParser::startElement(const XML_Char* name, const XML_Char **atts)
{
	if(strcmp(name,"FilterRequest")==0) {
		switchToNewSubParser( 
			mCurrentFilterParser = new FilterRequestParser(this, atts) 
		);  // we own and will have to explicitly delete 
...
}
		
void MultiFilterParser::endElement(const XML_Char *name)
{
	if(strcmp(name,"FilterRequest")==0) {
		assert(mCurrentFilterParser);
		FilterClause* newClause = mCurrentFilterParser->orphanBuiltClause();  // retrieve data built by sub-parser
...
		mCurrentFilterParser = 0;
		DeleteChild();
	}
}
\endverbatim
*/
class expatppNesting : public expatpp {

public:
	expatppNesting(expatppNesting* parent=0);  ///< NOT a copy ctor!! this is a recursive situation
	virtual ~expatppNesting();
	
	void switchToNewSubParser( expatppNesting* pAdoptedChild );
	expatppNesting* returnToParent();

protected:
	void BeAdopted(expatppNesting* adoptingParent);
	void OwnedChildOrphansItself(expatppNesting* callingChild);
	void RegisterWithParentXMLParser();
	virtual void AdoptChild(expatppNesting* adoptingChild);
	virtual void DeleteChild();
	
	int	mDepth;	
	bool mSelfDeleting;   ///< only valid if mParent not null
	expatppNesting* mParent; ///< may be null the parent owns this object
	expatppNesting* mOwnedChild;	///< owned, optional currently active child (auto_ptr not used to avoid STL dependency)

public:
/// \name interface functions for callbacks
//@{
	static void nestedStartElementCallback(void* userData, const XML_Char* name, const XML_Char** atts);
	static void nestedEndElementCallback(void* userData, const XML_Char* name);
//@}


/// \name overrideables to customise behaviour, must call parent
//@{
	virtual void SetupHandlers();
//@}

private:
	// Forbid copy-construction and assignment, to prevent double-deletion of mOwnedChild
						expatppNesting( const expatppNesting & );
	expatppNesting &	operator=( const expatppNesting & );
};


// inlines

// -------------------------------------------------------
//      e x p a t p p
// -------------------------------------------------------
inline 	
expatpp::operator XML_Parser() const
{
	return mParser;
}


// -------------------------------------------------------
//      e x p a t p p N e s t i n g
// -------------------------------------------------------
inline void
expatppNesting::OwnedChildOrphansItself(expatppNesting* callingChild)
{
	assert(callingChild==mOwnedChild);
	mOwnedChild = 0;
}


	
#endif   // H_EXPATPP
