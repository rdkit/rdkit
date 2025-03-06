// CommonCS/LibCommon/Src/CDXDocumentProperty.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 2016, PerkinElmer, Inc., All Rights Reserved

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
