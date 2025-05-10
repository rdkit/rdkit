//
// CDXColoredMolecularArea.h
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
