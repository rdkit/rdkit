//
// CDXColoredMolecularArea.cpp
//
// Copyright (c) 2018 PerkinElmer, Inc. All rights reserved.

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
//

#include "CDXStdObjects.h"
#include "CDXColoredMolecularArea.h"
#include "CDXMLNames.h"

CDXColoredMolecularArea::CDXColoredMolecularArea(CDXObjectID id) :
    CDXObject(CDXDatumID::kCDXObj_ColoredMolecularArea, id),
    backgroundColor(0)
{
}

CDXTag CDXColoredMolecularArea::GetTag() const
{
    return kCDXObj_ColoredMolecularArea;
}

std::string CDXColoredMolecularArea::XMLObjectName() const
{
    return kCDXML_coloredmoleculararea;
}

CDXObject* CDXColoredMolecularArea::Clone() const
{
    return new CDXColoredMolecularArea(*this);
}

void CDXColoredMolecularArea::StoreAttribute(CDXDataSource& ds, CDXTag tag, size_t size)
{
    switch (tag)
    {
    case kCDXProp_BackgroundColor:
        backgroundColor = UINT16(ds.GetUINT(size));
        break;

    case kCDXProp_BasisObjects:
        basisObjects = ReadObjectIDList(ds, size);
        break;

    default:
        CDXObject::StoreAttribute(ds, tag, size);
        break;
    }
}

void CDXColoredMolecularArea::WriteAttributesTo(CDXDataSink& ds) const
{
    CDXObject::WriteAttributesTo(ds);

    ds.PutAttribute(kCDXProp_BackgroundColor, backgroundColor);
    ds.PutAttributeForObjectIDList(kCDXProp_BasisObjects, basisObjects);
}

void CDXColoredMolecularArea::XMLStoreAttribute(CDXTag tag, const std::string& str)
{
    switch (tag)
    {
    case kCDXProp_BackgroundColor:
        backgroundColor = atoi(str.c_str());
        break;

    case kCDXProp_BasisObjects:
        {
            basisObjects.clear();
            std::istringstream is(str);

            std::copy(std::istream_iterator<CDXObjectID>(is),
                std::istream_iterator<CDXObjectID>(),
                std::back_inserter(basisObjects));
        }
        break;

    default:
        CDXObject::XMLStoreAttribute(tag, str);
        break;
    }
}

void CDXColoredMolecularArea::XMLWriteAttributes(XMLDataSink& ds) const
{
    CDXObject::XMLWriteAttributes(ds);

    CDXMLPutAttribute(ds, kCDXProp_BackgroundColor, backgroundColor);
    CDXMLPutAttributeForObjectIDList(ds, kCDXProp_BasisObjects, basisObjects);
}
