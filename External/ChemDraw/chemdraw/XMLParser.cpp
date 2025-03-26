// CommonCS/LibCommon/Src/XMLParser.cpp
// Copyright © 1999-2007, CambridgeSoft Corp., All Rights Reserved

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


#include "XMLParser.h"
#include "CDXUnicode.h"
#include "CDXMLNames.h"		// For kCDXML_id

#include <fstream>
#include <sstream>

using namespace std;


// *************************
// class XMLElementHandler
//
// XMLElementHandler is an abstract base class for classes which handle elements
// that do not correspond to CDX Objects, such as string, fonttable, etc.
// *************************

void XMLElementHandler::endElement(const XML_Char* /* name */)
{
	// The default endElement handler does nothing.
	// This is OK for elements which store their info as the chardata and subelements come along.
}

void XMLElementHandler::charData(const XML_Char* /* s */, int /* len */)
{
	// The default chardata handler ignores the characters.
	// This is OK for elements which aren't supposed to have any chardata, but
	// we throw an exception if some comes in, because we really don't want to be doing this.
	// A little whitespace is allowed, though.
#if 0 // default is to ignore this data
	for (int i = 0;  i < len;  ++i)
		if (s[i] != ' ' && s[i] != 0x0A && s[i] != 0x0D && s[i] != '\t')
			throw runtime_error("unexpected character data");
#endif
}

void XMLElementHandler::startElement(const XML_Char* /* name */, const XML_Char** atts)
{
	// The startElement method is called when an element is found inside of
	// the element which this handler is handling.
	// The default is to push an unknown element handler on the patser's
	// element stack, so that all elements and attributes will be ignored.

	m_parser->PushHandler(new XMLUnknownHandler(m_parser, atts));
}


#pragma mark
XMLUnknownHandler::XMLUnknownHandler(XMLParser *parser, const XML_Char **)
	: XMLElementHandler(parser)
{
	// Ignore all contained attributes and elements
}


CDXObjectID
XMLObjectIDTranslator::XMLObjectID(const char *s)
{
	IDMap::const_iterator i = m_idMap.find(s);
	return (i == m_idMap.end()) ? kCDXObj_UnknownObject : i->second;
}


// *************************
// class XMLParser
//
// wrapper for the expat parser for XML
// *************************
XMLParser::XMLParser(XMLObjectIDTranslator* objectIDTranslator)
	: expatpp(true) // tell base class constructor to create a parser
	, m_root(NULL)
	, m_ObjectIDTranslator(NULL) // set in subclass
{
	m_ObjectIDTranslator = objectIDTranslator;
}

XMLParser::~XMLParser()
{
	// This should only be non-empty if an exception was thrown while we were handling an element.
	while (!m_elementHandler.empty())
	{
		delete m_elementHandler.top();
		m_elementHandler.pop();
	}
	delete m_root;
	m_root = NULL;
}

void XMLParser::startElement(const XML_Char *name, const XML_Char **atts)
{
	CDXTag objTag = m_ObjectIDTranslator->XMLObjectID(name);

	if (objTag == 0 || objTag == kCDXObj_UnknownObject)
		PushHandler(new XMLUnknownHandler(this, atts));
	else
	{
		// Everything else is a normal CDX object
		startCDXObject(objTag, atts);
	}
}

void XMLParser::startCDXObject(CDXTag objTag, const XML_Char **atts)
{
	CDXObject *obj = NULL;

	// Figure out the CDXTag for each attribute
	const XML_Char **p;
	for (p = atts;  *p != 0;  p += 2)
	{
		if (strcmp(p[0], kCDXML_id) == 0)
		{
			CDXObjectID objID = atoi(p[1]);

			// Allocate the object, and store the attributes in it
			obj = m_CDXObjectFactory.AllocateObject(objID, objTag);
			break;
		}
	}

	if (obj == NULL)
		obj = m_CDXObjectFactory.AllocateObject(0, objTag);

	for (p = atts;  *p != 0;  p += 2)
	{
		if (strcmp(p[0], kCDXML_id) != 0)
		{
			CDXTag attributeID = getAttributeID(p[0]);

			//ASSERT(attributeID != 0);
			if (attributeID != 0	// unrecognized
				&& !IgnoreTag(attributeID, obj))
			{
				obj->XMLStoreAttribute(attributeID, p[1]);
			}
		}
	}
	
	// Remember this object in the ancestors stack
	if (m_ancestors.empty())
		m_root = obj;
	else
		CurrentObject()->AddChild(obj);
	m_ancestors.push(obj);
}

void XMLParser::endElement(const XML_Char* name)
{
	// If there is an element handler (e.g. for string), delegate to it.
	if (!m_elementHandler.empty())
	{
		m_elementHandler.top()->endElement(name);
		delete m_elementHandler.top();
		m_elementHandler.pop();
	}
	else
	{
		// Finish the object at the top of the ancestors stack, and pop the stack
		assert(m_ObjectIDTranslator->XMLObjectID(name) == CurrentObject()->GetTag());
		CurrentObject()->FinishReading();
		m_ancestors.pop();
	}
}

void XMLParser::charData(const XML_Char *s, int len)
{
	// If there is an element handler (e.g. for string, delegate to it.
	if (!m_elementHandler.empty())
		m_elementHandler.top()->charData(s, len);
	else
		CurrentObject()->XMLStoreCharacterData(string(s, len));
}

int XMLParser::notStandaloneHandler()
{
	return 1;
}
