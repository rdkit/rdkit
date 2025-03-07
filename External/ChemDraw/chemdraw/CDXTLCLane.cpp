// CommonCS/LibCommon/Src/CDXTLCLane.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright 2001-2004, CambridgeSoft Corp., All Rights Reserved

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
// ** class CDXPlateLane  **
// ***********************
//
// Specialization of CDXObject for CDXPlateLane objects

CDXPlateLane::CDXPlateLane(CDXTag tag, CDXObjectID id)
	:	CDXGraphicObject(tag, id)
{
}

CDXPlateLane::CDXPlateLane (const CDXPlateLane &src)
	:	CDXGraphicObject(src)
{
}

CDXTLCLane::CDXTLCLane(CDXObjectID id)
	: CDXPlateLane(kCDXObj_GEPLane, id)
{
}

CDXTLCLane::CDXTLCLane(const CDXTLCLane &a)
	: CDXPlateLane(a)
{
}

CDXObject*	CDXTLCLane::Clone() const
{
	return new CDXTLCLane(*this);
}

std::string CDXTLCLane::XMLObjectName() const
{
	return kCDXML_tlclane;
}

// CDXGEPLane class

CDXGEPLane::CDXGEPLane(CDXObjectID id)
	: CDXPlateLane(kCDXObj_GEPLane, id)
{
}

CDXGEPLane::CDXGEPLane(const CDXGEPLane &a)
	: CDXPlateLane(a)
{
}

CDXObject*	CDXGEPLane::Clone() const
{
	return new CDXGEPLane(*this);
}

std::string CDXGEPLane::XMLObjectName() const
{
	return kCDXML_geplane;
}

void CDXGEPLane::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	CDXPlateLane::XMLStoreAttribute(attribTag_arg, attribValue_arg);
}

void CDXGEPLane::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXPlateLane::XMLWriteAttributes(sink_arg);
}

void CDXGEPLane::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXPlateLane::WriteAttributesTo(sink_arg);
}

void CDXGEPLane::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	CDXPlateLane::StoreAttribute(src_arg, attribTag_arg, size_arg);
}