// CommonCS/LibCommon/Src/CDXMLParser.cpp
// Copyright 1999-2004, CambridgeSoft Corp., All Rights Reserved

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

#include "CDXStdObjects.h"
#include "CDXMLParser.h"
#include "CDXMLNames.h"
#include "CDXUnicode.h"
#include "CDXFontTable.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "cs_lock.h"

static cs::CriticalSection s_lock;

using namespace std;

CDXMLObjectIDTranslator::CDXMLObjectIDTranslator()
{	
	PROTECT_GLOBAL_AND_STATIC_DATA(s_lock);
	m_idMap.insert(IDMap::value_type(kCDXML_CDXML,				kCDXObj_Document));
	m_idMap.insert(IDMap::value_type(kCDXML_page,				kCDXObj_Page));
	m_idMap.insert(IDMap::value_type(kCDXML_group,				kCDXObj_Group));
	m_idMap.insert(IDMap::value_type(kCDXML_fragment,			kCDXObj_Fragment));
	m_idMap.insert(IDMap::value_type(kCDXML_node,				kCDXObj_Node));
	m_idMap.insert(IDMap::value_type(kCDXML_bond,				kCDXObj_Bond));
	m_idMap.insert(IDMap::value_type(kCDXML_text,				kCDXObj_Text));
	m_idMap.insert(IDMap::value_type(kCDXML_graphic,			kCDXObj_Graphic));
	m_idMap.insert(IDMap::value_type(kCDXML_curve,				kCDXObj_Curve));
	m_idMap.insert(IDMap::value_type(kCDXML_embeddedobject,		kCDXObj_EmbeddedObject));
	m_idMap.insert(IDMap::value_type(kCDXML_altgroup,			kCDXObj_NamedAlternativeGroup));
	m_idMap.insert(IDMap::value_type(kCDXML_templategrid,		kCDXObj_TemplateGrid));
	m_idMap.insert(IDMap::value_type(kCDXML_regnum,				kCDXObj_RegistryNumber));
	m_idMap.insert(IDMap::value_type(kCDXML_scheme,				kCDXObj_ReactionScheme));
	m_idMap.insert(IDMap::value_type(kCDXML_step,				kCDXObj_ReactionStep));
	m_idMap.insert(IDMap::value_type(kCDXML_objectdefinition,	kCDXObj_ObjectDefinition));
	m_idMap.insert(IDMap::value_type(kCDXML_spectrum,			kCDXObj_Spectrum));
	m_idMap.insert(IDMap::value_type(kCDXML_objecttag,			kCDXObj_ObjectTag));
	m_idMap.insert(IDMap::value_type(kCDXML_sequence,			kCDXObj_Sequence));
	m_idMap.insert(IDMap::value_type(kCDXML_crossreference,		kCDXObj_CrossReference));
	m_idMap.insert(IDMap::value_type(kCDXML_splitter,			kCDXObj_Splitter));
	m_idMap.insert(IDMap::value_type(kCDXML_table,				kCDXObj_Table));
	m_idMap.insert(IDMap::value_type(kCDXML_bracketedgroup,		kCDXObj_BracketedGroup));
	m_idMap.insert(IDMap::value_type(kCDXML_bracketattachment,	kCDXObj_BracketAttachment));
	m_idMap.insert(IDMap::value_type(kCDXML_crossingbond,		kCDXObj_CrossingBond));
	m_idMap.insert(IDMap::value_type(kCDXML_border,				kCDXObj_Border));
	m_idMap.insert(IDMap::value_type(kCDXML_geometry,			kCDXObj_Geometry));
	m_idMap.insert(IDMap::value_type(kCDXML_constraint,			kCDXObj_Constraint));
	m_idMap.insert(IDMap::value_type(kCDXML_tlcplate,			kCDXObj_TLCPlate));
	m_idMap.insert(IDMap::value_type(kCDXML_gepplate,			kCDXObj_GEPPlate));
	m_idMap.insert(IDMap::value_type(kCDXML_tlclane,			kCDXObj_TLCLane));
	m_idMap.insert(IDMap::value_type(kCDXML_geplane,			kCDXObj_GEPLane));
	m_idMap.insert(IDMap::value_type(kCDXML_tlcspot,			kCDXObj_TLCSpot));
	m_idMap.insert(IDMap::value_type(kCDXML_gepband,			kCDXObj_GEPBand));
	m_idMap.insert(IDMap::value_type(kCDXML_chemicalproperty,	kCDXObj_ChemicalProperty));
	m_idMap.insert(IDMap::value_type(kCDXML_arrow,				kCDXObj_Arrow));
	m_idMap.insert(IDMap::value_type(kCDXML_stoichiometrygrid,	kCDXObj_StoichiometryGrid));
	m_idMap.insert(IDMap::value_type(kCDXML_sgcomponent,		kCDXObj_SGComponent));
	m_idMap.insert(IDMap::value_type(kCDXML_sgdatum,			kCDXObj_SGDatum));
	m_idMap.insert(IDMap::value_type(kCDXML_bioshape,			kCDXObj_BioShape));
	m_idMap.insert(IDMap::value_type(kCDXML_plasmidmap,			kCDXObj_PlasmidMap));
	m_idMap.insert(IDMap::value_type(kCDXML_plasmidmarker,		kCDXObj_PlasmidMarker));
	m_idMap.insert(IDMap::value_type(kCDXML_marker,				kCDXObj_Marker));
	m_idMap.insert(IDMap::value_type(kCDXML_plasmidregion,		kCDXObj_PlasmidRegion));
	m_idMap.insert(IDMap::value_type(kCDXML_rlogic,				kCDXObj_RLogic));
	m_idMap.insert(IDMap::value_type(kCDXML_rlogicitem,			kCDXObj_RLogicItem));
	m_idMap.insert(IDMap::value_type(kCDXML_annotation,			kCDXObj_Annotation));
    m_idMap.insert(IDMap::value_type(kCDXML_documentproperties, kCDXObj_DocumentProperties));
    m_idMap.insert(IDMap::value_type(kCDXML_property,           kCDXObj_Property));
    m_idMap.insert(IDMap::value_type(kCDXML_coloredmoleculararea, kCDXObj_ColoredMolecularArea));
}

CDXMLObjectIDTranslator *
CDXMLObjectIDTranslator::GetTranslator()
{
	// This functione exists just to hold the static object, to avoid constructing it at startup.
	// This static object will be constructed the first time we're called.
	static CDXMLObjectIDTranslator theTranslator;
	return &theTranslator;
}

// *************************
// class CDXMLStringHandler
//
// CDXMLStringHandler is the handler for string elements
// A string element represents a chunk of text in a single font, face, size, and color.
// *************************
#pragma mark
class CDXMLStringHandler : public XMLElementHandler
{
	CDXString m_curString;
public:
	~CDXMLStringHandler() {}
	CDXMLStringHandler(CDXMLParser *parser, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
};

CDXMLStringHandler::CDXMLStringHandler(CDXMLParser *parser, const XML_Char **atts)
	:	XMLElementHandler(parser)
	,	m_curString("")
{
	// Store this string's attributes if any were specified
	int stringFont = UINT16(-1);
	int stringSize = 0;
	int stringFace = 0;
	int stringColor = 3;
	int stringAlpha = 0xFFFF;
	bool has_alpha(false);

	for (const XML_Char **p = atts;  *p != 0;  p += 2)
		if (strcmp(p[0], kCDXML_font) == 0)
			stringFont = atoi(p[1]);
		else if (strcmp(p[0], kCDXML_size) == 0)
			stringSize = atof(p[1]) * 20.0;
		else if (strcmp(p[0], kCDXML_face) == 0)
			stringFace = atoi(p[1]);
		else if (strcmp(p[0], kCDXML_color) == 0)
			stringColor = atoi(p[1]);
		else if (strcmp(p[0], kCDXML_alpha) == 0)
		{
			// Read alpha value
			stringAlpha = atoi(p[1]);
			has_alpha = true;
		}

	if (stringFont != UINT16(-1))
		m_curString.AddStyle(CDXStyle(0, stringFont, stringFace, stringSize, stringColor, stringAlpha, has_alpha));
}

void CDXMLStringHandler::endElement(const XML_Char* name)
{
	assert(strcmp(name, kCDXML_string) == 0);
	if (!m_curString.styles().empty())
	{
		try {
			CDXMLParser* parser = dynamic_cast<CDXMLParser*>(Parser());
			if (parser != NULL)
			{
				CDXFontTableEntry f = parser->GetDocument()->LookupFont(m_curString.style(0).family);
				CDXCharSet actualCharSet = f.m_type;
				m_curString.SetText(UnicodeToString(f.m_family, &actualCharSet, m_curString.str()));
			}
		}
		catch (...)
		{
			m_curString.SetText(m_curString.str());
		}
	}
	else
	{
		m_curString.SetText(m_curString.str());
	}
	Parser()->CurrentObject()->AppendText(m_curString);
}

void CDXMLStringHandler::charData(const XML_Char *s, int len)
{
	// The XML parser normalizes line ending chars to newline ('\n').
	// ChemDraw, due to its Macintosh origins, expects carriage returns
	// ('\r') to represent newlines in a string. Translate them here.
	//
	string str(s, len);
	replace(str.begin(), str.end(), '\n', '\r');
	m_curString += str;
}

#pragma mark
class CDXMLColorHandler : public XMLElementHandler
{
public:
	~CDXMLColorHandler() {}
	CDXMLColorHandler(CDXMLParser *parser, const XML_Char **atts);
};

CDXMLColorHandler::CDXMLColorHandler(CDXMLParser *parser, const XML_Char **atts)
	: XMLElementHandler(parser)
{
	CDXColor c(0,0,0);
	for (const XML_Char **p = atts;  *p != 0;  p += 2)
	{
		if (strcmp(p[0], kCDXML_red) == 0)
			c.m_red = atof(p[1]) * 255;
		else if (strcmp(p[0], kCDXML_green) == 0)
			c.m_green = atof(p[1]) * 255;
		else if (strcmp(p[0], kCDXML_blue) == 0)
			c.m_blue = atof(p[1]) * 255;
	}
	unsigned short comp;
	comp = c.m_red;		c.m_red = (comp << 8) + comp;
	comp = c.m_green;	c.m_green = (comp << 8) + comp;
	comp = c.m_blue;	c.m_blue = (comp << 8) + comp;

	parser->GetDocument()->m_colorTable.AddColor(c);
}

// *************************
// class CDXMLColorTableHandler
//
// element handler for colortable elements
// *************************
#pragma mark
class CDXMLColorTableHandler : public XMLElementHandler
{
	CDXDocument *m_document;
public:
	~CDXMLColorTableHandler() {}
	CDXMLColorTableHandler(CDXMLParser *parser, const XML_Char **atts);
	void startElement(const XML_Char* name, const XML_Char** atts);
};

CDXMLColorTableHandler::CDXMLColorTableHandler(CDXMLParser *parser, const XML_Char **atts)
	: XMLElementHandler(parser)
	, m_document(dynamic_cast<CDXDocument *>(Parser()->CurrentObject()))
{
	if (m_document == 0)
		throw runtime_error("colortable not contained in a document");
	m_document->m_colorTable.m_colors.resize(2);	// Eliminate all but the basic black and white
}

void CDXMLColorTableHandler::startElement(const XML_Char* name, const XML_Char** atts)
{
	XMLElementHandler* handler;
	CDXMLParser* parser;
	if (strcmp(name, kCDXML_colorElement) == 0 && ((parser = dynamic_cast<CDXMLParser*>(Parser())) != NULL) )
		handler = new CDXMLColorHandler(parser, atts);
	else
		handler = new XMLUnknownHandler(Parser(), atts);

	Parser()->PushHandler(handler);
}

#pragma mark
class CDXMLPropRepHandler : public XMLElementHandler
{
public:
	~CDXMLPropRepHandler() {}
	CDXMLPropRepHandler(CDXMLParser *parser, const XML_Char **atts);
};

CDXMLPropRepHandler::CDXMLPropRepHandler(CDXMLParser *parser, const XML_Char **atts)
	: XMLElementHandler(parser)
{
	CDXObjectID objID = kCDXUndefinedId;
	CDXTag propTag = 0;
	for (const XML_Char **p = atts;  *p != 0;  p += 2)
	{
		if (strcmp(p[0], kCDXML_attribute) == 0)
			propTag = CDXMLAttributeID(p[1]);
		else if (strcmp(p[0], kCDXML_object) == 0)
			objID = atoi(p[1]);
	}
	if (objID != kCDXUndefinedId && propTag != 0)
	{
		CDXGraphicObject *obj = dynamic_cast<CDXGraphicObject *>(parser->CurrentObject());
		// We only handle PropRep objects in graphic objects, ignore them otherwise
		if (obj != 0)
			obj->PropertyRepresented(CDXPropRep(objID, propTag));
	}
}

// *************************
// class CDXMLParser
//
// wrapper for the expat parser for CDXML
// *************************
CDXMLParser::CDXMLParser()
	: XMLParser(CDXMLObjectIDTranslator::GetTranslator())
{
	CDXAddStandardFactories(m_CDXObjectFactory);
}

CDXMLParser::~CDXMLParser()
{
}

void CDXMLParser::startElement(const XML_Char *name, const XML_Char **atts)
{
	if (!m_elementHandler.empty())
		m_elementHandler.top()->startElement(name, atts);
	else if (strcmp(name, kCDXML_string) == 0)
		PushHandler(new CDXMLStringHandler(this, atts));
	else if (strcmp(name, kCDXML_fonttable) == 0)
		PushHandler(new CDXMLFontTableHandler(this, atts));
	else if (strcmp(name, kCDXML_colortable) == 0)
		PushHandler(new CDXMLColorTableHandler(this, atts));
	else if (strcmp(name, kCDXML_representElement) == 0)
		PushHandler(new CDXMLPropRepHandler(this, atts));
	else
	{
		CDXTag objTag = m_ObjectIDTranslator->XMLObjectID(name);

		if (objTag == 0 || objTag == kCDXObj_UnknownObject)
			PushHandler(new XMLUnknownHandler(this, atts));
		else if (IgnoreObject(objTag, CurrentObject()))
			PushHandler(new XMLUnknownHandler(this, atts));
		else
			startCDXObject(objTag, atts);
	}
}

CDXDocument *CDXMLParser::GetDocument()
{
	return dynamic_cast<CDXDocument*>(m_root);
}

std::unique_ptr<CDXDocument> CDXMLParser::ReleaseDocument()
{
	CDXDocument *doc = dynamic_cast<CDXDocument*>(m_root);
	m_root = NULL;
	std::unique_ptr<CDXDocument> w(doc);
	return w;
}

CDXTag CDXMLParser::getAttributeID(const char *s) const
{
	CDXTag retVal = CDXMLAttributeID(s);
	//ASSERT(retVal != 0);	// unknown attribute name
	return retVal;
}


#ifdef FULLASSERT
void TestCDXParser(char *s, unsigned long siz);

void TestCDXParser(char *s, unsigned long siz)
{
	CDXMLParser parser;
	if (!parser.XML_Parse(s, siz, true))
	{
		ostringstream os;
		os << "Error reading CDXML document at line " << parser.XML_GetCurrentLineNumber() << ": " << XML_ErrorString(parser.XML_GetErrorCode());
		throw runtime_error( os.str() );
	}
	else
	{
		ofstream dest("parseout.xml");
		dest << kCDXML_HeaderString;
		XMLDataSink ds(dest);
		parser.GetDocument()->XMLWrite(ds);
		delete parser.GetDocument();
	}
}
#endif // FULLASSERT
