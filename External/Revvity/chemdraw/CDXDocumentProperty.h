// CommonCS/LibCommon/Hdr/CDXDocumentProperty.h
// Copyright (c) 2016, PerkinElmer, Inc., All Rights Reserved

#pragma once

#include "CoreChemistryAPI.h"
#include "CDXStdObjects.h"

// ***************************************
// ***** class CDXDocumentProperty  ******
// ***************************************
//
// Specialization of CDXObject for CDXDocumentProperty objects

#pragma mark
class CORE_CHEMISTRY_API CDXDocumentProperty : public CDXObject
{
public:
    CDXDocumentProperty(CDXObjectID id = 0);

    virtual CDXObject*	Clone()	const;

	const CDXString& GetName() const { return m_name; }
	void SetName(const CDXString& name) { m_name = name; }

	CDXPropertyRule GetRule() const { return m_rule; }
	void SetRule(CDXPropertyRule rule) { m_rule = rule; }

	CDXPropertyDataType GetType() const { return m_type; }
	void SetType(CDXPropertyDataType type) { m_type = type; }

	const CDXString& GetValue() const { return m_value; }
	CDXString& GetValue() { return m_value; }
	void SetValue(const CDXString& val) { m_value = val; }

	static const std::string& RuleToString(CDXPropertyRule rule);
	static CDXPropertyRule StringToRule(const std::string& str);
	static const std::string& DataTypeToString(CDXPropertyDataType rule);
	static CDXPropertyDataType StringToDataType(const std::string& str);

protected:
    virtual CDXTag	GetTag() const;
    virtual std::string XMLObjectName() const;

    // Serialize and Deserialize method for CDX 
    virtual void StoreAttribute(CDXDataSource& ds, CDXTag tag, size_t size);
    virtual void WriteAttributesTo(CDXDataSink& ds) const;

    // Serialize and Deserialize method for CDXML 
    virtual void XMLStoreAttribute(CDXTag tag, const std::string& val);
    virtual void XMLStoreCharacterData(const std::string &);
    virtual void XMLWriteAttributes(XMLDataSink& ds) const;
    virtual bool XMLNeedToWriteContent() const;
    virtual void XMLWriteContent(XMLDataSink &) const;

private:
    CDXString	    m_name;
    CDXPropertyRule m_rule;
    CDXPropertyDataType m_type;
    CDXString       m_value;
};

