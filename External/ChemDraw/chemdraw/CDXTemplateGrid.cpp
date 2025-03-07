// CommonCS/LibCommon/Src/CDXTemplateGrid.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2004, CambridgeSoft Corp., All Rights Reserved

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
#include "CDXMLNames.h"
#include <stdlib.h>	// for atoi and friends
#include <ostream>

// ***************************
// ** class CDXTemplateGrid **
// ***************************
//
// Specialization of CDXObject for CDXTemplateGrid objects

CDXTemplateGrid::CDXTemplateGrid(CDXObjectID id_arg)
	: CDXObject(kCDXObj_TemplateGrid, id_arg)
	, m_paneHeight(0)
	, m_numRows(0)
	, m_numColumns(0)
	, m_extent(0, 0)
{
}

CDXTemplateGrid::CDXTemplateGrid (const CDXTemplateGrid &src)
	: CDXObject		(src)
	, m_paneHeight	(src.m_paneHeight)
	, m_numRows		(src.m_numRows)
	, m_numColumns	(src.m_numColumns)
	, m_extent		(src.m_extent)
{
}

CDXTemplateGrid::~CDXTemplateGrid()
{
}

CDXObject*	CDXTemplateGrid::Clone() const
{
	return new CDXTemplateGrid (*this);
}

std::string CDXTemplateGrid::XMLObjectName() const
{
	return kCDXML_templategrid;
}

void CDXTemplateGrid::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Template_PaneHeight:
		m_paneHeight = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Template_NumRows:
		m_numRows = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_Template_NumColumns:
		m_numColumns = src_arg.GetUINT(size_arg);
		break;
	case kCDXProp_2DExtent:
	{	if (size_arg != 8)
			throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
		CDXCoordinate y = src_arg.GetINT32();
		CDXCoordinate x = src_arg.GetINT32();
		m_extent = CDXPoint2D(x, y);
		break;
	}
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXTemplateGrid::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Template_PaneHeight:
		m_paneHeight = atoi(value_arg.c_str());
		break;
	case kCDXProp_Template_NumRows:
		m_numRows = atoi(value_arg.c_str());
		break;
	case kCDXProp_Template_NumColumns:
		m_numColumns = atoi(value_arg.c_str());
		break;
	case kCDXProp_2DExtent:
		m_extent = StringToCDXPoint2D(value_arg);
		break;
	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, value_arg);
		break;
	}
}

void CDXTemplateGrid::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);
	sink_arg.PutAttribute( kCDXProp_Template_PaneHeight, m_paneHeight );
	sink_arg.PutAttribute( kCDXProp_Template_NumRows, m_numRows );
	sink_arg.PutAttribute( kCDXProp_Template_NumColumns, m_numColumns );
	sink_arg.PutAttribute( kCDXProp_2DExtent, m_extent );
}

void CDXTemplateGrid::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXObject::XMLWriteAttributes(sink_arg);
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_PaneHeight, m_paneHeight );
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_NumRows, m_numRows );
	CDXMLPutAttribute( sink_arg, kCDXProp_Template_NumColumns, m_numColumns );
	CDXMLPutAttribute( sink_arg, kCDXProp_2DExtent, m_extent );
}
