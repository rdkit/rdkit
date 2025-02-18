// CommonCS/LibCommon/Hdr/XMLParser.h
// Contains: Definitions for objects that hold the content of CML files
// Copyright © 1999-2004, CambridgeSoft Corp., All Rights Reserved

#pragma once

// --> https://github.com/gittiver/libexpatpp/blob/main/src/expatpp.hpp
#include "expatpp.h"

#include "CoreChemistryAPI.h"
#include "CDXObject.h"

#include <stack>
#include <memory>

// Defined in this header file:
class XMLParser;

// *************************
// class XMLElementHandler
//
// XMLElementHandler is an abstract base class for classes which handle elements
// that do not correspond to CDX Objects, such as string, fonttable, etc.
// *************************

class CORE_CHEMISTRY_API XMLElementHandler
{
	// Keep track of the parser so we can access things like the current CDX object
	XMLParser * const m_parser;

protected:
	XMLParser *Parser() const { return m_parser; }

public:
	XMLElementHandler(XMLParser *parser) : m_parser(parser) {}
	virtual ~XMLElementHandler() {};

	// The endElement method is called when the element end tag is read.
	// This is used to store the element in parser->CurrentObject().
	// It need not be overridden if the data is stored as the chardata and subelements are read.
	virtual void endElement(const XML_Char* name);

	// The charData method is called when some character data is contained in this element.
	// It need not be overriden if no character data is expected, but if it's not, and some
	// data shows up, the parse is aborted.
	virtual void charData(const XML_Char *s, int len);

	// The startElement method is called when an element is found inside of
	// the element which this handler is handling.
	// The default is to push an unknown element handler on the patser's
	// element stack, so that all elements and attributes will be ignored.
	virtual void startElement(const XML_Char* name, const XML_Char** atts);
};

class XMLUnknownHandler
	: public XMLElementHandler
{
public:
	~XMLUnknownHandler() {}
	XMLUnknownHandler(XMLParser *parser, const XML_Char **atts);
};

template<class T>
struct XMLCharPtrLess
{
	inline bool operator()(const T& t1, const T& t2) const
		{return (strcmp(t1, t2) < 0); }
};

// *************************
// class XMLObjectIDTranslator
// *************************
class XMLObjectIDTranslator
{
protected:
	typedef std::map<const char *, CDXObjectID, XMLCharPtrLess<const char *> > IDMap;
	IDMap m_idMap;

public:
	XMLObjectIDTranslator() {};
	CDXObjectID XMLObjectID(const char *s);
};

class CORE_CHEMISTRY_API XMLParser : public expatpp
{
protected:
	CDXObjectFactory m_CDXObjectFactory;
	std::stack<CDXObject *> m_ancestors;
	std::stack<XMLElementHandler *> m_elementHandler;
	CDXObject *m_root;

	XMLObjectIDTranslator* m_ObjectIDTranslator;
//	static CDXObjectID XMLObjectID(const std::string &s);

	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
	virtual int notStandaloneHandler();

	virtual void startCDXObject(CDXTag objTag, const XML_Char **atts);
	
	// Helper method implemented by subclasses
	virtual CDXTag getAttributeID(const char *s) const = 0;

public:
	XMLParser(XMLObjectIDTranslator* objectIDTranslator);

	virtual ~XMLParser();

	CDXObject *GetRoot() { return m_root; }
	const CDXObject *GetRoot() const { return m_root; }
	CDXObject *CurrentObject() { return m_ancestors.empty() ? NULL : m_ancestors.top(); }
	const CDXObject *CurrentObject() const { return m_ancestors.empty() ? NULL : m_ancestors.top(); }
	void PushHandler(XMLElementHandler *eh) { m_elementHandler.push(eh); }

	virtual bool IgnoreTag(CDXTag tag, const CDXObject *parent) const { return false; }
};
