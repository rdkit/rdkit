// CommonCS/LibCommon/Src/CDXDocumentProperty.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 2016, PerkinElmer, Inc., All Rights Reserved

#include "CDXStdObjects.h"
#include "CDXDocumentProperty.h"
#include "CDXMLNames.h"
#include "cs_stringUtils.h"

#include <memory>
#include <vector>
#include <stdexcept>

using std::string;
using std::invalid_argument;

const string DocumentPropertyTypeIntegerString = "Integer";
const string DocumentPropertyTypeStringString = "String";
const string DocumentPropertyTypeFloatString = "Float";

const string DocumentPropertyRuleOptionalString = "Optional";
const string DocumentPropertyRuleRecommendedString = "Recommended";
const string DocumentPropertyRuleRequiredString = "Required";

const string DocumentPropertyCommonErrorString = "There was a problem while reading document properties. ";
const string DocumentPropertyRuleErrorString = DocumentPropertyCommonErrorString + "Unknown document property rule found.";
const string DocumentPropertyTypeErrorString = DocumentPropertyCommonErrorString + "Unknown document property type found.";

CDXDocumentProperty::CDXDocumentProperty(CDXObjectID id) :
    CDXObject(CDXDatumID::kCDXObj_Property, id),
    m_name(),
    m_rule(CDXPropertyRule::kCDXPropertyRule_Optional),
    m_type(CDXPropertyDataType::kCDXPropertyDataType_String)
{
}

CDXTag CDXDocumentProperty::GetTag() const
{
    return kCDXObj_Property;
}

std::string CDXDocumentProperty::XMLObjectName() const
{
    return kCDXML_property;
}

CDXObject* CDXDocumentProperty::Clone() const
{
    return new CDXDocumentProperty(*this);
}

void CDXDocumentProperty::StoreAttribute(CDXDataSource& ds, CDXTag tag, size_t size)
{
    switch (tag)
    {
    case kCDXProp_Name:
        m_name.Read(ds, size);
        break;
    case kCDXProp_Property_Rule:
        m_rule = (CDXPropertyRule)ds.GetINT(size);
        break;
    case kCDXProp_Property_DataType:
        m_type = (CDXPropertyDataType)ds.GetINT(size);
        break;
    case kCDXProp_Property_Value:
        m_value.Read(ds, size);
        break;
    default:
        CDXObject::StoreAttribute(ds, tag, size);
        break;
    }
}

void CDXDocumentProperty::WriteAttributesTo(CDXDataSink& ds) const
{
    CDXObject::WriteAttributesTo(ds);

    ds.PutAttributeCDXString(kCDXProp_Name, m_name);
    ds.PutAttribute(kCDXProp_Property_Rule, INT8(m_rule));
    ds.PutAttribute(kCDXProp_Property_DataType, INT8(m_type));
    ds.PutAttributeCDXString(kCDXProp_Property_Value, m_value);
}

void CDXDocumentProperty::XMLStoreAttribute(CDXTag tag, const std::string& str)
{
    switch (tag)
    {
    case kCDXProp_Name:
        m_name.SetText(str);
        break;
    case kCDXProp_Property_Rule:
        m_rule = StringToRule(str);
        break;
    case kCDXProp_Property_DataType:
        m_type = StringToDataType(str);
        break;
    default:
        CDXObject::XMLStoreAttribute(tag, str);
        break;
    }
}

void CDXDocumentProperty::XMLStoreCharacterData(const std::string& str)
{
    m_value.SetText(m_value.str() + str);
}

void CDXDocumentProperty::XMLWriteAttributes(XMLDataSink& ds) const
{
    CDXObject::XMLWriteAttributes(ds);

    CDXMLPutAttribute(ds, kCDXProp_Name, m_name.str());
    CDXMLPutAttribute(ds, kCDXProp_Property_Rule, RuleToString(m_rule));
    CDXMLPutAttribute(ds, kCDXProp_Property_DataType, DataTypeToString(m_type));
}

bool CDXDocumentProperty::XMLNeedToWriteContent() const
{
    return true;
}

void CDXDocumentProperty::XMLWriteContent(XMLDataSink& ds) const
{
    XMLPut(ds, m_value.str());
}

const std::string& CDXDocumentProperty::RuleToString(CDXPropertyRule rule)
{
    switch (rule)
    {
    case kCDXPropertyRule_Optional:
        return DocumentPropertyRuleOptionalString;
    case kCDXPropertyRule_Recommended:
        return DocumentPropertyRuleRecommendedString;
    case kCDXPropertyRule_Required:
        return DocumentPropertyRuleRequiredString;
    }

    throw invalid_argument(DocumentPropertyRuleErrorString);
}

CDXPropertyRule CDXDocumentProperty::StringToRule(const std::string& str)
{
    if (cs::str_compareNoCase(str, DocumentPropertyRuleOptionalString) == 0)
    {
        return kCDXPropertyRule_Optional;
    }
    else if (cs::str_compareNoCase(str, DocumentPropertyRuleRecommendedString) == 0)
    {
        return kCDXPropertyRule_Recommended;
    }
    else if (cs::str_compareNoCase(str, DocumentPropertyRuleRequiredString) == 0)
    {
        return kCDXPropertyRule_Required;
    }

    throw invalid_argument(DocumentPropertyRuleErrorString);        
}


const std::string& CDXDocumentProperty::DataTypeToString(CDXPropertyDataType rule)
{
    switch (rule)
    {
    case kCDXPropertyDataType_String:
        return DocumentPropertyTypeStringString;
    case kCDXPropertyDataType_Int:
        return DocumentPropertyTypeIntegerString;
    case kCDXPropertyDataType_Float:
        return DocumentPropertyTypeFloatString;
    }

    throw invalid_argument(DocumentPropertyTypeErrorString);
}

CDXPropertyDataType CDXDocumentProperty::StringToDataType(const std::string& str)
{
    if (cs::str_compareNoCase(str, DocumentPropertyTypeStringString) == 0)
    {
        return kCDXPropertyDataType_String;
    }
    else if (cs::str_compareNoCase(str, DocumentPropertyTypeIntegerString) == 0)
    {
        return kCDXPropertyDataType_Int;
    }
    else if (cs::str_compareNoCase(str, DocumentPropertyTypeFloatString) == 0)
    {
        return kCDXPropertyDataType_Float;
    }

    throw invalid_argument(DocumentPropertyTypeErrorString);
}
