// CommonCS/LibCommon/Hdr/XMLParser.h
// Contains: Definitions for objects that hold the content of CML files
// Copyright © 1999-2004, CambridgeSoft Corp., All Rights Reserved

// BSD 3-Clause License
// 
// Copyright (c) 1986-2025, CambridgeSoft Corp, Revvity Inc and others.
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
