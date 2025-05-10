// CommonCS/LibCommon/Hdr/CDXDocumentProperty.h
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

