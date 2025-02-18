// CommonCS/LibCommon/Src/CDXDocumentPropertyCollection.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 2016, PerkinElmer, Inc., All Rights Reserved

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
