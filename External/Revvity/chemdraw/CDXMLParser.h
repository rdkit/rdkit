// CommonCS/LibCommon/Hdr/CDXMLParser.h
// Copyright (c) 1999-2004, CambridgeSoft Corp., All Rights Reserved

#pragma once

#include "CoreChemistryAPI.h"
#include "XMLParser.h"

// Defined in this header file:
class CDXMLElementHandler;
class CDXMLParser;

// Used in this file
class CDXDocument;

// *************************
// class XMLObjectIDTranslator
//
// this is a singleton, used to hold the translation table for cdxml element names
// *************************
class CDXMLObjectIDTranslator : public XMLObjectIDTranslator
{
public:
	CDXMLObjectIDTranslator();
	static CDXMLObjectIDTranslator* GetTranslator();
};

class CORE_CHEMISTRY_API CDXMLParser : public XMLParser
{
protected:
	virtual void startElement(const XML_Char *name, const XML_Char **atts);

	virtual CDXTag getAttributeID(const char *s) const;

	virtual bool IgnoreObject(CDXTag objTag, const CDXObject *parent) const { return false; }

public:
	CDXMLParser();

	virtual ~CDXMLParser();

	// Note: GetDocument does not relinquish ownership of the document, while ReleaseDocument does.
	CDXDocument *GetDocument();
	std::unique_ptr<CDXDocument> ReleaseDocument();
};
