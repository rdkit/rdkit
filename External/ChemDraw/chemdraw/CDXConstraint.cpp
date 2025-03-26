// CommonCS/LibCommon/Src/CDXConstraint.cpp
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

#include <ostream>

#include "CDXStdObjects.h"
#include "XMLPutEnum.h"
#include "CDXMLNames.h"
#include <sstream>
#include <stdexcept>

#if TARGET_OS_WIN32
#pragma warning( disable : 4786 ) // identifier was truncated to 'number' characters in the debug information
#endif

XMLPutEnum<CDXConstraintType>::Value sXMLConstraintTypeValues[] = {
	{kCDXConstraintType_Undefined,			"Undefined"},
	{kCDXConstraintType_Distance,			"Distance"},
	{kCDXConstraintType_Angle,				"Angle"},
	{kCDXConstraintType_ExclusionSphere,	"ExclusionSphere"}
};

XMLPutEnum<CDXConstraintType> sXMLConstraintType(sXMLConstraintTypeValues, sizeof sXMLConstraintTypeValues, kCDXConstraintType_Undefined);
void XMLPut(XMLDataSink &sink, CDXConstraintType  v) {	sink.os << sXMLConstraintType.lookup(v); }

// ***************************
// **  class CDXConstraint  **
// ***************************
//
// Specialization of CDXObject for CDXConstraint objects

CDXConstraint::CDXConstraint(CDXObjectID id) 
: 	CDXGraphicObject(kCDXObj_Constraint, id), 
	m_constraintType(kCDXConstraintType_Undefined), 
	m_min(0), m_max(0), 
	m_ignoreUnconnectedAtoms(false),
	m_dihedralIsChiral(false)
{
}

CDXConstraint::CDXConstraint (const CDXConstraint &src) 
:	CDXGraphicObject(src), 
	m_constraintType(src.m_constraintType), 
	m_min(src.m_min), 
	m_max(src.m_max),
	m_ignoreUnconnectedAtoms(src.m_ignoreUnconnectedAtoms),
	m_dihedralIsChiral(src.m_dihedralIsChiral),
	m_basisObjects(src.m_basisObjects)
{
}


CDXObject* CDXConstraint::Clone() const
{
	return new CDXConstraint(*this);
}


std::string CDXConstraint::XMLObjectName() const
{
	return kCDXML_constraint;
}

void CDXConstraint::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch (attribTag_arg)
	{
	case kCDXProp_ConstraintType:
		m_constraintType = (CDXConstraintType) src_arg.GetUINT(size_arg);
		break;

	case kCDXProp_ConstraintMin:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_min = SwapBytes(f64);
		break;

	case kCDXProp_ConstraintMax:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_max = SwapBytes(f64);
		break;

	case kCDXProp_IgnoreUnconnectedAtoms:
		m_ignoreUnconnectedAtoms = size_arg == 0 || src_arg.GetUINT(size_arg) != 0;
		break;

	case kCDXProp_DihedralIsChiral:
		m_dihedralIsChiral = size_arg == 0 || src_arg.GetUINT(size_arg) != 0;
		break;

	case kCDXProp_BasisObjects:
        m_basisObjects = ReadObjectIDList(src_arg, size_arg);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}

}

void CDXConstraint::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_ConstraintType:
		m_constraintType = sXMLConstraintType.lookup(attribValue_arg);
		break;

	case kCDXProp_ConstraintMin:
		m_min = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_ConstraintMax:
		m_max = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_IgnoreUnconnectedAtoms:
		m_ignoreUnconnectedAtoms = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_DihedralIsChiral:
		m_dihedralIsChiral = attribValue_arg == "yes" || attribValue_arg == "true";
		break;

	case kCDXProp_BasisObjects:
		{
			std::istringstream is(attribValue_arg);
			m_basisObjects.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_basisObjects));
		}
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}


void CDXConstraint::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	if (m_constraintType != kCDXConstraintType_Undefined)
		sink_arg.PutAttribute( kCDXProp_ConstraintType, (INT8) m_constraintType );
	
	sink_arg.PutAttribute( kCDXProp_ConstraintMin, FLOAT64(m_min));
	sink_arg.PutAttribute( kCDXProp_ConstraintMax, FLOAT64(m_max));

	if (m_ignoreUnconnectedAtoms)
		sink_arg.PutAttribute( kCDXProp_IgnoreUnconnectedAtoms );

	if (m_dihedralIsChiral)
		sink_arg.PutAttribute( kCDXProp_DihedralIsChiral );

    sink_arg.PutAttributeForObjectIDList(kCDXProp_BasisObjects, m_basisObjects);
}


void CDXConstraint::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	if (m_constraintType != kCDXConstraintType_Undefined)
		CDXMLPutAttribute(sink_arg, kCDXProp_ConstraintType, m_constraintType );

	CDXMLPutAttribute(sink_arg, kCDXProp_ConstraintMin, m_min);
	CDXMLPutAttribute(sink_arg, kCDXProp_ConstraintMax, m_max);

	if (m_ignoreUnconnectedAtoms)
		CDXMLPutAttribute(sink_arg, kCDXProp_IgnoreUnconnectedAtoms );

	if (m_dihedralIsChiral)
		CDXMLPutAttribute(sink_arg, kCDXProp_DihedralIsChiral );

    CDXMLPutAttributeForObjectIDList(sink_arg, kCDXProp_BasisObjects, m_basisObjects);
}



