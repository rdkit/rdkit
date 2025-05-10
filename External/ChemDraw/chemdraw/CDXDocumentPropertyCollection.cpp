// CommonCS/LibCommon/Src/CDXDocumentPropertyCollection.cpp
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
#include "CDXDocumentPropertyCollection.h"
#include "CDXMLNames.h"

#include <memory>
#include <vector>
#include <cctype>
#include <stdexcept>

CDXDocumentPropertyCollection::CDXDocumentPropertyCollection(CDXObjectID id)
    : CDXObject(kCDXObj_DocumentProperties, id)
{
}

CDXTag CDXDocumentPropertyCollection::GetTag() const
{
    return kCDXObj_DocumentProperties;
}

std::string CDXDocumentPropertyCollection::XMLObjectName() const
{
    return kCDXML_documentproperties;
}

CDXObject*	CDXDocumentPropertyCollection::Clone()	const
{
    return new CDXDocumentPropertyCollection(*this);
}

void CDXDocumentPropertyCollection::AddProperty(CDXDocumentProperty *prop)
{
    AddChild(prop);
}

void CDXDocumentPropertyCollection::RemoveProperty(CDXDocumentProperty *prop)
{
    RemoveChild(prop, true);
}
