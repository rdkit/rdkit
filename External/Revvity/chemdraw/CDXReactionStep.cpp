// CommonCS/LibCommon/Src/CDXReactionStep.cpp
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
#include "CDXMLNames.h"
#include <sstream>
#include <stdexcept>

// *****************************
// ** class CDXReactionScheme **
// *****************************
//
CDXReactionScheme::CDXReactionScheme(CDXObjectID id_arg)
	:	CDXObject(kCDXObj_ReactionScheme, id_arg)
{
}

CDXReactionScheme::CDXReactionScheme (const CDXReactionScheme &src)
	:	CDXObject				(src)
{
}

CDXReactionScheme::~CDXReactionScheme()
{
}

CDXObject*	CDXReactionScheme::Clone() const
{
	return new CDXReactionScheme (*this);
}

std::string CDXReactionScheme::XMLObjectName() const
{
	return kCDXML_scheme;
}



// ***************************
// ** class CDXReactionStep **
// ***************************
//
CDXReactionStep::CDXReactionStep(CDXObjectID id_arg)
	:	CDXObject(kCDXObj_ReactionStep, id_arg)
{
}

CDXReactionStep::CDXReactionStep (const CDXReactionStep &src)
	:	CDXObject				(src)
	,	m_reactants				(src.m_reactants)
	,	m_products				(src.m_products)
	,	m_plusses				(src.m_plusses)
	,	m_arrows				(src.m_arrows)
	,	m_aamap					(src.m_aamap)
	,	m_aamapManual			(src.m_aamapManual)
	,	m_aamapAuto				(src.m_aamapAuto)
	,	m_objectsAboveArrow		(src.m_objectsAboveArrow)
	,	m_objectsBelowArrow		(src.m_objectsBelowArrow)
{
}

CDXReactionStep::~CDXReactionStep()
{
}

CDXObject*	CDXReactionStep::Clone() const
{
	return new CDXReactionStep (*this);
}

std::string CDXReactionStep::XMLObjectName() const
{
	return kCDXML_step;
}

void CDXReactionStep::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_ReactionStep_Atom_Map:
		{
			if ((size_arg % 8) != 0)
				throw std::runtime_error ("Bad atom-atom map attribute");
			size_t	numPairs = size_arg / 8;
			m_aamap.clear();
			while (numPairs--)
			{
				const CDXObjectID	iObj_1 = src_arg.GetUINT32(),
									iObj_2 = src_arg.GetUINT32();
				if (iObj_1 < 0  ||  iObj_2 < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in atom-atom map");
				m_aamap.push_back (std::make_pair (iObj_1, iObj_2));
			}
		}
		break;
	case kCDXProp_ReactionStep_Atom_Map_Manual:
		{
			if ((size_arg % 8) != 0)
				throw std::runtime_error ("Bad atom-atom map attribute");
			size_t	numPairs = size_arg / 8;
			m_aamapManual.clear();
			while (numPairs--)
			{
				const CDXObjectID	iObj_1 = src_arg.GetUINT32(),
									iObj_2 = src_arg.GetUINT32();
				if (iObj_1 < 0  ||  iObj_2 < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in atom-atom map");
				m_aamapManual.push_back (std::make_pair (iObj_1, iObj_2));
			}
		}
		break;
	case kCDXProp_ReactionStep_Atom_Map_Auto:
		{
			if ((size_arg % 8) != 0)
				throw std::runtime_error ("Bad atom-atom map attribute");
			size_t	numPairs = size_arg / 8;
			m_aamapAuto.clear();
			while (numPairs--)
			{
				const CDXObjectID	iObj_1 = src_arg.GetUINT32(),
									iObj_2 = src_arg.GetUINT32();
				if (iObj_1 < 0  ||  iObj_2 < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in atom-atom map");
				m_aamapAuto.push_back (std::make_pair (iObj_1, iObj_2));
			}
		}
		break;
	case kCDXProp_ReactionStep_Reactants:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad reactants attribute");
			size_t	num = size_arg / 4;
			m_reactants.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in reactants map");
				m_reactants.push_back(iObj);
			}
		}
		break;
	case kCDXProp_ReactionStep_Products:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad products attribute");
			size_t	num = size_arg / 4;
			m_products.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in products map");
				m_products.push_back(iObj);
			}
		}
		break;
	case kCDXProp_ReactionStep_Plusses:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad plusses attribute");
			size_t	num = size_arg / 4;
			m_plusses.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in plusses map");
				m_plusses.push_back(iObj);
			}
		}
		break;
	case kCDXProp_ReactionStep_Arrows:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad arrows attribute");
			size_t	num = size_arg / 4;
			m_arrows.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in arrows map");
				m_arrows.push_back(iObj);
			}
		}
		break;
	case kCDXProp_ReactionStep_ObjectsAboveArrow:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad objects above arrow attribute");
			size_t	num = size_arg / 4;
			m_objectsAboveArrow.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in objects above arrow map");
				m_objectsAboveArrow.push_back(iObj);
			}
		}
		break;
	case kCDXProp_ReactionStep_ObjectsBelowArrow:
		{
			if ((size_arg % 4) != 0)
				throw std::runtime_error ("Bad objects below arrow attribute");
			size_t	num = size_arg / 4;
			m_objectsBelowArrow.resize (0);
			while (num--)
			{
				const CDXObjectID	iObj = src_arg.GetUINT32();
				if (iObj < 0)	// todo: check upper limit
					throw std::runtime_error ("Illegal obj# in objects below arrow map");
				m_objectsBelowArrow.push_back(iObj);
			}
		}
		break;
	default:
		CDXObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

std::istream& operator>>(std::istream &is, std::pair<CDXObjectID, CDXObjectID> &p);
std::istream& operator>>(std::istream &is, std::pair<CDXObjectID, CDXObjectID> &p)
{
	return is >> p.first >> p.second;
}

void CDXReactionStep::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_ReactionStep_Atom_Map:
		{
			std::istringstream is(attribValue_arg);
			m_aamap.clear();
			while (!is.eof())
			{
				std::pair<CDXObjectID, CDXObjectID> p;
				is >> p;
				m_aamap.push_back(p);
			}
		}
		break;
	case kCDXProp_ReactionStep_Atom_Map_Manual:
		{
			std::istringstream is(attribValue_arg);
			m_aamapManual.clear();
			while (!is.eof())
			{
				std::pair<CDXObjectID, CDXObjectID> p;
				is >> p;
				m_aamapManual.push_back(p);
			}
		}
		break;
	case kCDXProp_ReactionStep_Atom_Map_Auto:
		{
			std::istringstream is(attribValue_arg);
			m_aamapAuto.clear();
			while (!is.eof())
			{
				std::pair<CDXObjectID, CDXObjectID> p;
				is >> p;
				m_aamapAuto.push_back(p);
			}
		}
		break;
	case kCDXProp_ReactionStep_Reactants:
		{
			std::istringstream is(attribValue_arg);
			m_reactants.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_reactants));
		}
		break;
	case kCDXProp_ReactionStep_Products:
		{
			std::istringstream is(attribValue_arg);
			m_products.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_products));
		}
		break;
	case kCDXProp_ReactionStep_Plusses:
		{
			std::istringstream is(attribValue_arg);
			m_plusses.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					 std::istream_iterator< CDXObjectID >(),
					 std::back_inserter(m_plusses));
		}
		break;
	case kCDXProp_ReactionStep_Arrows:
		{
			std::istringstream is(attribValue_arg);
			m_arrows.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					 std::istream_iterator< CDXObjectID >(),
					 std::back_inserter(m_arrows));
		}
		break;
	case kCDXProp_ReactionStep_ObjectsAboveArrow:
		{
			std::istringstream is(attribValue_arg);
			m_objectsAboveArrow.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_objectsAboveArrow));
		}
		break;
	case kCDXProp_ReactionStep_ObjectsBelowArrow:
		{
			std::istringstream is(attribValue_arg);
			m_objectsBelowArrow.clear();
			std::copy(std::istream_iterator< CDXObjectID >(is),
					  std::istream_iterator< CDXObjectID >(),
					  std::back_inserter(m_objectsBelowArrow));
		}
		break;

	default:
		CDXObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXReactionStep::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXObject::WriteAttributesTo(sink_arg);
	if (!m_reactants.empty())
	{
		const size_t	numBytes = m_reactants.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Reactants );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_reactants.begin();  it != m_reactants.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_products.empty())
	{
		const size_t	numBytes = m_products.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Products );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_products.begin();  it != m_products.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_plusses.empty())
	{
		const size_t	numBytes = m_plusses.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Plusses );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_plusses.begin();  it != m_plusses.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_arrows.empty())
	{
		const size_t	numBytes = m_arrows.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Arrows );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_arrows.begin();  it != m_arrows.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_objectsAboveArrow.empty())
	{
		const size_t	numBytes = m_objectsAboveArrow.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_ObjectsAboveArrow );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_objectsAboveArrow.begin();  it != m_objectsAboveArrow.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_objectsBelowArrow.empty())
	{
		const size_t	numBytes = m_objectsBelowArrow.size() * 4;
		sink_arg.PutTag ( kCDXProp_ReactionStep_ObjectsBelowArrow );
		sink_arg.Put( UINT16( numBytes ) );
		for (std::vector<CDXObjectID>::const_iterator it = m_objectsBelowArrow.begin();  it != m_objectsBelowArrow.end();  it++)
			sink_arg.Put ((UINT32)*it);
	}
	if (!m_aamap.empty())
	{
		const size_t	numBytes = m_aamap.size() * 8;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Atom_Map );
		sink_arg.Put( UINT16( numBytes ) );
		for (Map::const_iterator it = m_aamap.begin();  it != m_aamap.end();  it++)
		{
			sink_arg.Put ((UINT32) it->first);
			sink_arg.Put ((UINT32) it->second);
		}
	}
	if (!m_aamapManual.empty())
	{
		const size_t	numBytes = m_aamapManual.size() * 8;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Atom_Map_Manual );
		sink_arg.Put( UINT16( numBytes ) );
		for (Map::const_iterator it = m_aamapManual.begin();  it != m_aamapManual.end();  it++)
		{
			sink_arg.Put ((UINT32) it->first);
			sink_arg.Put ((UINT32) it->second);
		}
	}
	if (!m_aamapAuto.empty())
	{
		const size_t	numBytes = m_aamapAuto.size() * 8;
		sink_arg.PutTag ( kCDXProp_ReactionStep_Atom_Map_Auto );
		sink_arg.Put( UINT16( numBytes ) );
		for (Map::const_iterator it = m_aamapAuto.begin();  it != m_aamapAuto.end();  it++)
		{
			sink_arg.Put ((UINT32) it->first);
			sink_arg.Put ((UINT32) it->second);
		}
	}
}

void CDXReactionStep::XMLWriteAttributes(XMLDataSink &sink) const
{
	CDXObject::XMLWriteAttributes(sink);
	if (!m_reactants.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Reactants) << "=\"";
		sink.os << m_reactants;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_products.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Products) << "=\"";
		sink.os << m_products;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_plusses.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Plusses) << "=\"";
		sink.os << m_plusses;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_arrows.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Arrows) << "=\"";
		sink.os << m_arrows;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_objectsAboveArrow.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_ObjectsAboveArrow) << "=\"";
		sink.os << m_objectsAboveArrow;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_objectsBelowArrow.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_ObjectsBelowArrow) << "=\"";
		sink.os << m_objectsBelowArrow;
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_aamap.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Atom_Map) << "=\"";
		for (Map::const_iterator it = m_aamap.begin();  it != m_aamap.end();  it++)
		{
			if (it != m_aamap.begin())
				sink.os << " ";
			sink.os << it->first << " " << it->second;
		}
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_aamapManual.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Atom_Map_Manual) << "=\"";
		for (Map::const_iterator it = m_aamapManual.begin();  it != m_aamapManual.end();  it++)
		{
			if (it != m_aamapManual.begin())
				sink.os << " ";
			sink.os << it->first << " " << it->second;
		}
		sink.os << "\"" << GetTextEOL();
	}
	if (!m_aamapAuto.empty())
	{
		sink.os << std::string(" ") << CDXMLAttributeName(kCDXProp_ReactionStep_Atom_Map_Auto) << "=\"";
		for (Map::const_iterator it = m_aamapAuto.begin();  it != m_aamapAuto.end();  it++)
		{
			if (it != m_aamapAuto.begin())
				sink.os << " ";
			sink.os << it->first << " " << it->second;
		}
		sink.os << "\"" << GetTextEOL();
	}
}
