//
// CDXColoredMolecularArea.cpp
//
// Copyright (c) 2018 PerkinElmer, Inc. All rights reserved.
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
