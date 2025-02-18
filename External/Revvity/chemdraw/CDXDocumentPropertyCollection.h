// CommonCS/LibCommon/Hdr/CDXDocumentPropertyCollection.h
// Copyright (c) 2016, PerkinElmer, Inc., All Rights Reserved

#pragma once

#include "CoreChemistryAPI.h"
#include "cs_univDefs.h"
#include "CDXStdObjects.h"
#include "CDXDocumentProperty.h"

// *************************************************
// ***** class CDXDocumentPropertyCollection  ******
// *************************************************
//
// Specialization of CDXObject for CDXDocumentPropertyCollection objects

class CORE_CHEMISTRY_API CDXDocumentPropertyCollection : public CDXObject
{
public:
    CDXDocumentPropertyCollection(CDXObjectID id = 0);

    virtual CDXObject*	Clone()	const;

    void AddProperty(CDXDocumentProperty *prop);
    void RemoveProperty(CDXDocumentProperty *prop);

protected:
    virtual CDXTag	GetTag() const;
    virtual std::string XMLObjectName() const;
};

