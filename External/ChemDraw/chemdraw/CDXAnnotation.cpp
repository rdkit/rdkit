// CommonCS/LibCommon/Src/CDXAnnotation.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 2001-2007, CambridgeSoft Corp., All Rights Reserved

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
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <sstream>

// *************************
// ** class CDXAnnotation **
// *************************
//
// Specialization of CDXObject for CDXAnnotation objects

CDXAnnotation::CDXAnnotation(CDXObjectID id_arg, const string &keyword, const string &content)
	:	CDXObject(kCDXObj_Annotation, id_arg)
	,	m_keyword(keyword)
	,	m_content(content)
{
}

CDXAnnotation::CDXAnnotation (const CDXAnnotation &src)
	:	CDXObject(src)
	,	CDX_INHERIT_SRC (m_keyword)
	,	CDX_INHERIT_SRC (m_content)
{
}

CDXObject*	CDXAnnotation::Clone() const
{
	return new CDXAnnotation(*this);
}

std::string CDXAnnotation::XMLObjectName() const
{
	return kCDXML_annotation;
}

void CDXAnnotation::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Annotation_Keyword:
		m_keyword.Read(src_arg, size_arg);
		break;
	case kCDXProp_Annotation_Content:
		m_content.Read(src_arg, size_arg);
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXAnnotation::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Annotation_Keyword:
		m_keyword = CDXString(UnicodeToString(attribValue_arg));
		break;
	case kCDXProp_Annotation_Content:
		m_content = CDXString(UnicodeToString(attribValue_arg));
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXAnnotation::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);

	if (!m_keyword.empty())
		sink_arg.PutAttributeCDXString( kCDXProp_Annotation_Keyword, m_keyword );

	if (!m_content.empty())
		sink_arg.PutAttributeCDXString( kCDXProp_Annotation_Content, m_content );
}

void CDXAnnotation::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);

	if (!m_keyword.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Annotation_Keyword, StringToUnicode("", m_keyword.str()) );

	if (!m_content.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Annotation_Content, StringToUnicode("", m_content.str()) );

}



