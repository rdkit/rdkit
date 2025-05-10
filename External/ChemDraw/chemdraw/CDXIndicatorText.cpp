// CommonCS/LibCommon/Src/CDXIndicatorText.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2004, CambridgeSoft Corp., All Rights Reserved

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

// ***********************
// ** class CDXSequence **
// ***********************
//
// Specialization of CDXObject for CDXSequence objects

CDXSequence::CDXSequence(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_Sequence, id_arg)
{
}

CDXSequence::CDXSequence (const CDXSequence &src)
	:	CDXGraphicObject(src)
	,	CDX_INHERIT_SRC (m_identifier)
{
}

CDXObject*	CDXSequence::Clone() const
{
	return new CDXSequence(*this);
}

std::string CDXSequence::XMLObjectName() const
{
	return kCDXML_sequence;
}

void CDXSequence::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Sequence_Identifier:
		m_identifier.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_identifier[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_identifier [0];
			int	numSkipBytes = (int)m_identifier.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_identifier.erase(0, numSkipBytes);
		}
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXSequence::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Sequence_Identifier:
		m_identifier = attribValue_arg;
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXSequence::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (!m_identifier.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_Sequence_Identifier );
        CDXPut(sink_arg, m_identifier);
	}
}

void CDXSequence::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (!m_identifier.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_Sequence_Identifier, m_identifier );
}



// *****************************
// ** class CDXCrossReference **
// *****************************
//
// Specialization of CDXObject for CDXCrossReference objects

CDXCrossReference::CDXCrossReference(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_CrossReference, id_arg)
{
}

CDXCrossReference::CDXCrossReference (const CDXCrossReference &src)
	:	CDXGraphicObject(src)
	,	CDX_INHERIT_SRC (m_container)
	,	CDX_INHERIT_SRC (m_document)
	,	CDX_INHERIT_SRC (m_identifier)
	,	CDX_INHERIT_SRC (m_sequence)
{
}

CDXObject*	CDXCrossReference::Clone() const
{
	return new CDXCrossReference(*this);
}

std::string CDXCrossReference::XMLObjectName() const
{
	return kCDXML_crossreference;
}

void CDXCrossReference::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_CrossReference_Container:
		m_container.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_container[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_container [0];
			int	numSkipBytes = (int)m_container.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_container.erase(0, numSkipBytes);
		}
		break;
	case kCDXProp_CrossReference_Document:
		m_document.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_document[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_document [0];
			int	numSkipBytes = (int)m_document.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_document.erase(0, numSkipBytes);
		}
		break;
	case kCDXProp_CrossReference_Identifier:
		m_identifier.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_identifier[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_identifier [0];
			int	numSkipBytes = (int)m_identifier.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_identifier.erase(0, numSkipBytes);
		}
		break;
	case kCDXProp_CrossReference_Sequence:
		m_sequence.assign(src_arg.GetString(size_arg));
		// The count of styles is normally zero, but it might not be.
		if (size_arg > 1  &&  m_sequence[1] == 0)	// They included the styles header (we assume there couldn't be > 16 styles)
		{
			INT16 numStyles = m_sequence [0];
			int	numSkipBytes = (int)m_sequence.length();
			if (numSkipBytes > 2 + numStyles * 10)
				numSkipBytes = 2 + numStyles * 10;
			m_sequence.erase(0, numSkipBytes);
		}
		break;
	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXCrossReference::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_CrossReference_Container:
		m_identifier = attribValue_arg;
		break;
	case kCDXProp_CrossReference_Document:
		m_identifier = attribValue_arg;
		break;
	case kCDXProp_CrossReference_Identifier:
		m_identifier = attribValue_arg;
		break;
	case kCDXProp_CrossReference_Sequence:
		m_identifier = attribValue_arg;
		break;
	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXCrossReference::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (!m_container.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_CrossReference_Container );
		sink_arg.Put( UINT16( m_container.size() + 2 ) );
		sink_arg.Put( UINT16 (0) );	// the styles
		sink_arg.Put( (INT8 *)m_container.data(), m_container.size() );
	}
	if (!m_document.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_CrossReference_Document );
		sink_arg.Put( UINT16( m_document.size() + 2 ) );
		sink_arg.Put( UINT16 (0) );	// the styles
		sink_arg.Put( (INT8 *)m_document.data(), m_document.size() );
	}
	if (!m_identifier.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_CrossReference_Identifier );
		sink_arg.Put( UINT16( m_identifier.size() + 2 ) );
		sink_arg.Put( UINT16 (0) );	// the styles
		sink_arg.Put( (INT8 *)m_identifier.data(), m_identifier.size() );
	}
	if (!m_sequence.empty())
	{	// Expects a CDXString, even though only the text is of interest.
		sink_arg.PutTag( kCDXProp_CrossReference_Sequence );
		sink_arg.Put( UINT16( m_sequence.size() + 2 ) );
		sink_arg.Put( UINT16 (0) );	// the styles
		sink_arg.Put( (INT8 *)m_sequence.data(), m_sequence.size() );
	}
}

void CDXCrossReference::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (!m_container.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_CrossReference_Container, m_container );

	if (!m_document.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_CrossReference_Document, m_document );

	if (!m_identifier.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_CrossReference_Identifier, m_identifier );

	if (!m_sequence.empty())
		CDXMLPutAttribute(sink_arg, kCDXProp_CrossReference_Sequence, m_sequence );
}
