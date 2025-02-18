//
// CDXColoredMolecularArea.h
//
// Copyright (c) 2018 PerkinElmer, Inc. All rights reserved.
//

#pragma once

#include "CoreChemistryAPI.h"
#include "CDXStdObjects.h"

class CORE_CHEMISTRY_API CDXColoredMolecularArea : public CDXObject
{
public:

    CDXColoredMolecularArea(CDXObjectID id = 0);

    CDXObject* Clone() const override;

    UINT16 GetBackgroundColor() const { return backgroundColor; }
    void SetBackgroundColor(UINT16 bgColor) { backgroundColor = bgColor; }

    std::vector<CDXObjectID> GetBasisObjects() const { return basisObjects; }
    void SetBasisObjects(const std::vector<CDXObjectID>& inBasisObjects) { basisObjects = inBasisObjects; }

protected:

    CDXTag GetTag() const override;
    std::string XMLObjectName() const override;

    // Serialize and Deserialize method for CDX 
    void StoreAttribute(CDXDataSource& ds, CDXTag tag, size_t size) override;
    void WriteAttributesTo(CDXDataSink& ds) const override;

    // Serialize and Deserialize method for CDXML 
    void XMLStoreAttribute(CDXTag tag, const std::string& val) override;
    void XMLWriteAttributes(XMLDataSink& ds) const override;

private:

    UINT16 backgroundColor;
    std::vector<CDXObjectID> basisObjects;
};
