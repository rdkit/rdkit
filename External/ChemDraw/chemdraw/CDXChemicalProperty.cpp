// CommonCS/LibCommon/Src/CDXChemicalProperty.cpp
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
#include "CDXUnicode.h"
#include <stdlib.h>	// for atoi and friends
#include <sstream>
#include <stdexcept>
extern XMLPutEnum<CDXPositioningType> sXMLPositioningType;

namespace
{
    std::istream& operator>>(std::istream& is, CDXChemicalProperty::ExternalBond& externalBond)
    {
        return is >> externalBond.first >> externalBond.second;
    }
}

CDXChemicalProperty::~CDXChemicalProperty()
{
}

CDXChemicalProperty::CDXChemicalProperty(CDXObjectID id)
    : CDXObject(kCDXObj_ChemicalProperty, id)
    , m_isActive(false)
    , m_isChemicallySignificant(false)
    , m_displayID(0)
    , m_positioning(kCDXPositioningType_Auto)
    , m_positioningAngle(0.0)
    , m_positioningOffset(CDXPoint2D(0, 0))
{
}

CDXChemicalProperty::CDXChemicalProperty(const CDXChemicalProperty& src)
    : CDXObject(src)
    , CDX_INHERIT_SRC(m_name)
    , CDX_INHERIT_SRC(m_isActive)
    , CDX_INHERIT_SRC(m_isChemicallySignificant)
    , CDX_INHERIT_SRC(m_displayID)
    , CDX_INHERIT_SRC(m_basisObjects)
    , CDX_INHERIT_SRC(m_propertyTypes)
    , CDX_INHERIT_SRC(m_positioning)
    , CDX_INHERIT_SRC(m_positioningAngle)
    , CDX_INHERIT_SRC(m_positioningOffset)
    , CDX_INHERIT_SRC(m_externalBonds)
{
}

CDXObject* CDXChemicalProperty::Clone() const
{
    return new CDXChemicalProperty(*this);
}

string CDXChemicalProperty::XMLObjectName() const
{
    return kCDXML_chemicalproperty;
}

void CDXChemicalProperty::StoreAttribute(CDXDataSource& src_arg, CDXTag attribTag_arg, size_t size_arg)
{
    switch (attribTag_arg)
    {
    case kCDXProp_Name:
        m_name = src_arg.GetString(size_arg);
        break;

    case kCDXProp_ChemicalPropertyIsActive:
        m_isActive = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
        break;

    case kCDXProp_ChemicalPropertyIsChemicallySignificant:
        m_isChemicallySignificant = (size_arg == 0) || (src_arg.GetUINT(size_arg) != 0);
        break;

    case kCDXProp_ChemicalPropertyDisplayID:
        m_displayID = src_arg.GetINT(size_arg);
        break;

    case kCDXProp_ChemicalPropertyType:
    {
        if ((size_arg % 4) != 0)
        {
            throw std::runtime_error("Bad property types attribute");
        }

        size_t num = size_arg / 4;
        m_propertyTypes.clear();

        while (num--)
        {
            const CDXChemicalPropertyType iObj = (CDXChemicalPropertyType)src_arg.GetUINT32();
            if (iObj < 0)
            {
                continue;
            }

            m_propertyTypes.push_back(iObj);
        }

        break;
    }
    
    case kCDXProp_Positioning:
        m_positioning = (CDXPositioningType)src_arg.GetINT(size_arg);
        break;

    case kCDXProp_PositioningAngle:
        m_positioningAngle = src_arg.GetINT32();
        break;

    case kCDXProp_PositioningOffset:
    {
        if (size_arg != 8)
        {
            throw invalid_cdx_error(GetObjectID(), GetTag(), attribTag_arg);
        }

        CDXCoordinate y = src_arg.GetINT32();
        CDXCoordinate x = src_arg.GetINT32();
        m_positioningOffset = CDXPoint2D(x, y);
        break;
    }

    case kCDXProp_BasisObjects:
        m_basisObjects = ReadObjectIDList(src_arg, size_arg);
        break;

    case kCDXProp_ChemicalPropertyExternalBonds:
    {
        if ((size_arg % sizeof(ExternalBond)) != 0)
        {
            throw std::runtime_error("Bad external bonds attribute");
        }

        m_externalBonds.clear();
        size_t numPairs = size_arg / sizeof(ExternalBond);

        while (numPairs--)
        {
            const auto first = src_arg.GetINT32();
            const auto second = src_arg.GetINT32();
            m_externalBonds.push_back(std::make_pair(first, second));
        }

        break;
    }

    default:
        CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
        break;
    }
}

void CDXChemicalProperty::XMLStoreAttribute(CDXTag tag, const string& attribValue_arg)
{
    switch (tag)
    {
    case kCDXProp_Name:
        m_name = attribValue_arg;
        break;

    case kCDXProp_ChemicalPropertyIsActive:
        m_isActive = (attribValue_arg == "yes") || (attribValue_arg == "true");
        break;

    case kCDXProp_ChemicalPropertyIsChemicallySignificant:
        m_isChemicallySignificant = (attribValue_arg == "yes") || (attribValue_arg == "true");
        break;

    case kCDXProp_ChemicalPropertyDisplayID:
        m_displayID = atoi(attribValue_arg.c_str());
        break;

    case kCDXProp_ChemicalPropertyType:
    {
        vector< int > myTypes;
        std::istringstream is(attribValue_arg);

        std::copy(std::istream_iterator< int >(is),
            std::istream_iterator< int >(),
            std::back_inserter(myTypes));

        m_propertyTypes.clear();
        for (vector<int>::const_iterator i = myTypes.begin(); i != myTypes.end(); ++i)
        {
            m_propertyTypes.push_back((CDXChemicalPropertyType)*i);
        }

        break;
    }
    
    case kCDXProp_Positioning:
        m_positioning = sXMLPositioningType.lookup(attribValue_arg);
        break;

    case kCDXProp_PositioningAngle:
        m_positioningAngle = DegreesToCdxAngle(atof(attribValue_arg.c_str()));
        break;

    case kCDXProp_PositioningOffset:
        m_positioningOffset = StringToCDXPoint2D(attribValue_arg);
        break;

    case kCDXProp_BasisObjects:
    {
        m_basisObjects.clear();
        std::istringstream is(attribValue_arg);

        std::copy(std::istream_iterator< CDXObjectID >(is),
            std::istream_iterator< CDXObjectID >(),
            std::back_inserter(m_basisObjects));
        break;
    }

    case kCDXProp_ChemicalPropertyExternalBonds:
    {
        m_externalBonds.clear();

        std::istringstream is(attribValue_arg);
        while (!is.eof())
        {
            ExternalBond externalBond;
            is >> externalBond;
            m_externalBonds.push_back(externalBond);
        }

        break;
    }

    
    default:
        CDXObject::XMLStoreAttribute(tag, attribValue_arg);
        break;
    }
}

void CDXChemicalProperty::WriteAttributesTo(CDXDataSink& sink_arg) const
{
    CDXObject::WriteAttributesTo(sink_arg);
    sink_arg.PutAttribute(kCDXProp_ChemicalPropertyDisplayID, (UINT32)m_displayID);

    if (!m_name.empty())
    {
        sink_arg.PutAttributeString(kCDXProp_Name, m_name);
    }

    if (m_isActive)
    {
        sink_arg.PutAttribute(kCDXProp_ChemicalPropertyIsActive);
    }

    if (m_isChemicallySignificant)
    {
        sink_arg.PutAttribute(kCDXProp_ChemicalPropertyIsChemicallySignificant);
    }

    if (!m_propertyTypes.empty())
    {
        const size_t numBytes = m_propertyTypes.size() * 4;
        sink_arg.PutTag(kCDXProp_ChemicalPropertyType);
        sink_arg.Put(UINT16(numBytes));
        for (vector<CDXChemicalPropertyType>::const_iterator it = m_propertyTypes.begin(); it != m_propertyTypes.end(); ++it)
        {
            sink_arg.Put((UINT32)*it);
        }
    }

    if (m_positioning != kCDXPositioningType_Auto)
    {
        sink_arg.PutAttribute(kCDXProp_Positioning, (UINT8)m_positioning);
    }

    if (m_positioningAngle != 0)
    {
        sink_arg.PutAttribute(kCDXProp_PositioningAngle, m_positioningAngle);
    }

    if ((m_positioningOffset.x != 0) || (m_positioningOffset.y != 0))
    {
        sink_arg.PutAttribute(kCDXProp_PositioningOffset, m_positioningOffset);
    }

    sink_arg.PutAttributeForObjectIDList(kCDXProp_BasisObjects, m_basisObjects);
    
    if (!m_externalBonds.empty())
    {
        const size_t numBytes = m_externalBonds.size() * sizeof(ExternalBond);
        sink_arg.PutTag(kCDXProp_ChemicalPropertyExternalBonds);
        sink_arg.Put(UINT16(numBytes));
        for (vector<ExternalBond>::const_iterator it = m_externalBonds.begin(); it != m_externalBonds.end(); ++it)
        {
            sink_arg.Put(it->first);
            sink_arg.Put(it->second);
        }
    }
}

void CDXChemicalProperty::XMLWriteAttributes(XMLDataSink& sink_arg) const
{
    CDXObject::XMLWriteAttributes(sink_arg);
    CDXMLPutAttribute(sink_arg, kCDXProp_ChemicalPropertyDisplayID, m_displayID);

    if (!m_name.empty())
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Name, StringToUnicode("", m_name));
    }

    if (m_isActive)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_ChemicalPropertyIsActive);
    }

    if (m_isChemicallySignificant)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_ChemicalPropertyIsChemicallySignificant, true);
    }

    if (!m_propertyTypes.empty())
    {
        sink_arg.os << string(" ") << CDXMLAttributeName(kCDXProp_ChemicalPropertyType) << "=\"";

        for (vector< CDXChemicalPropertyType >::const_iterator it = m_propertyTypes.begin(); it != m_propertyTypes.end(); ++it)
        {
            if (it != m_propertyTypes.begin())
            {
                sink_arg.os << " ";
            }

            sink_arg.os << *it;
        }

        sink_arg.os << "\"" << GetTextEOL();
    }

    if (m_positioning != kCDXPositioningType_Auto)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_Positioning, m_positioning);
    }

    if (m_positioningAngle != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_PositioningAngle, CdxAngleToDegrees(m_positioningAngle));
    }

    if ((m_positioningOffset.x != 0) || (m_positioningOffset.y != 0))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_PositioningOffset, m_positioningOffset);
    }

    CDXMLPutAttributeForObjectIDList(sink_arg, kCDXProp_BasisObjects, m_basisObjects);

    if (!m_externalBonds.empty())
    {
        sink_arg.os << string(" ") << CDXMLAttributeName(kCDXProp_ChemicalPropertyExternalBonds) << "=\"";

        for (vector<ExternalBond>::const_iterator it = m_externalBonds.begin(); it != m_externalBonds.end(); ++it)
        {
            if (it != m_externalBonds.begin())
            {
                sink_arg.os << " ";
            }

            sink_arg.os << it->first << " " << it->second;
        }

        sink_arg.os << "\"" << GetTextEOL();
    }
}

bool CDXChemicalProperty::ShouldWriteToMDL() const
{
    if ((m_propertyTypes.size() == 1) && (m_propertyTypes[0] == kCDXChemicalPropertyTypeIUPACAtomNumber))
    {
        // Suppress storage of isolated IUPAC atom number property
        return false;
    }
    return true;
}

bool CDXChemicalProperty::HasAttachedDataType() const
{    
    for (const auto& type : m_propertyTypes)
    {
        if (type == kCDXChemicalPropertyTypeAttachedData)
        {
            return true;
        }
    }

    return false;
}