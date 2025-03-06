// CommonCS/LibCommon/Src/CDXColorTable.cpp
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

// *******************************
// ***** class CDXColorTable *****
// *******************************
//

namespace
{
	static const CDXColorIndex kUndefinedColorIndex = (CDXColorIndex)-1;
}

// Return the current size of the color table
CDXColorIndex CDXColorTable::NumColors() const
{
    return (CDXColorIndex)m_colors.size();
}

const CDXColor&	CDXColorTable::GetColor(CDXColorIndex index) const
{
    assert(index < NumColors());
    return m_colors[index];
}

CDXColorIndex CDXColorTable::FindOrAddColor(const CDXColor& color)
{
    const int val = GetColorsIndex(color);
    if (val != kUndefinedColorIndex)
    {
        return val;
    }

    return AddColor(color);
}

CDXColorIndex CDXColorTable::AddColor(const CDXColor& color)
{
    assert(m_colors.size() < UINT16_MAX);

    if (m_colors.size() < UINT16_MAX)
    {
		m_colors.push_back(color);
		return NumColors() - 1;
	}
	else
	{
		return kUndefinedColorIndex;
	}
}

CDXColorTable::CDXColorTable()
{
	m_colors.resize(4);
	m_colors [0] = CDXColor::Black();
	m_colors [1] = CDXColor::White();
	m_colors [2] = CDXColor::White();
	m_colors [3] = CDXColor::Black();
}

CDXColorIndex CDXColorTable::GetColorsIndex (const CDXColor& color, bool skipHardcoded /* = false */, bool getNearest /* = false */) const
{
	int index = -1;
	std::vector<CDXColor>::const_iterator it,
										  beg = m_colors.begin();
	if (skipHardcoded)
		beg += 2;
	it = std::find (beg, m_colors.end(), color);
	if (it == m_colors.end())
		index = -1;
	else
		index = int(it - m_colors.begin());

	// Try to find the nearest color by comparing "distance" between colors
	if (index == -1 && getNearest && m_colors.size() > 0)
	{
		UINT16 r1 = color.m_red;
		UINT16 g1 = color.m_green;
		UINT16 b1 = color.m_blue;

		long minDistSq = 65535;
		int i = 0;
		for (it = m_colors.begin(); it != m_colors.end(); ++it, i++)
		{
			UINT16 r2 = it->m_red;
			UINT16 g2 = it->m_green;
			UINT16 b2 = it->m_blue;

			long distSq = (r1 - r2) * (r1 - r2) + (g1 - g2) * (g1 - g2) + (b1 - b2) * (b1 - b2);
			if (distSq < minDistSq)
			{
				minDistSq = distSq;
				index = i;
			}
		}
	}

	return index;
}

void CDXColorTable::Read (CDXDataSource &src_arg, size_t size_arg)
{
	if (size_arg < 2)
		throw std::runtime_error ("Bad color table");
	UINT16 numEntries = src_arg.GetUINT16();
	if (size_arg - 2 != numEntries * 2 * 3)	// two bytes per R/G/B, three of those per color.
		throw std::runtime_error ("Bad color table");
	m_colors.resize (2+numEntries);	// keep the default first two values
	for (int i = 0;  i < numEntries;  i++)
	{
		m_colors[i+2].m_red   = src_arg.GetUINT16();
		m_colors[i+2].m_green = src_arg.GetUINT16();
		m_colors[i+2].m_blue  = src_arg.GetUINT16();

		unsigned short comp;
		comp = m_colors[i+2].m_red / 257;
		m_colors[i+2].m_red = (comp << 8) + comp;
		comp = m_colors[i+2].m_green / 257;
		m_colors[i+2].m_green = (comp << 8) + comp;
		comp = m_colors[i+2].m_blue / 257;
		m_colors[i+2].m_blue = (comp << 8) + comp;
	}
}

void CDXColorTable::Write(CDXDataSink &sink_arg) const	// write the color table to a binary CDX stream
{
	// A binary color table is the number of colors (UINT16)
	// followed by the r,g,b for each color (3 * UINT16)
	// We don't write the first two entries, which are always black and white.
	const int	numBytes = (NumColors() - 2) * 3 * 2  + 2;
	sink_arg.PutTag (kCDXProp_ColorTable);
	sink_arg.Put (UINT16 (numBytes));

	sink_arg.Put (UINT16(NumColors() - 2));
	for (std::vector<CDXColor>::const_iterator i = m_colors.begin() + 2;  i != m_colors.end();  i++)
	{
		sink_arg.Put (i->m_red);
		sink_arg.Put (i->m_green);
		sink_arg.Put (i->m_blue);
	}
}

std::istream & operator>>(std::istream &is, CDXColor& c)
{
	return is >> c.m_red >> c.m_green >> c.m_blue;
}

void CDXColorTable::XMLRead(const std::string &attribValue_arg)
{
	// Read from XML attribute string

	std::istringstream is(attribValue_arg);
	m_colors.resize(2);
	m_colors [0] = CDXColor::Black();
	m_colors [1] = CDXColor::White();
	std::copy(std::istream_iterator< CDXColor >(is),
			  std::istream_iterator< CDXColor >(),
			  std::back_inserter(m_colors));
}

std::ostream & operator<<(std::ostream &os, const CDXColor &c)
{
	return os << GetTextEOL() << "<" << kCDXML_colorElement
		<< " " << kCDXML_red   << "=\"" << CDXValueToString(c.m_red  /65535.0, 4) << "\""
		<< " " << kCDXML_green << "=\"" << CDXValueToString(c.m_green/65535.0, 4) << "\""
		<< " " << kCDXML_blue  << "=\"" << CDXValueToString(c.m_blue /65535.0, 4) << "\""
		<< "/>";
}

void CDXColorTable::XMLWrite(XMLDataSink &sink) const
{
	sink.os << "<" << kCDXML_colortable << ">";
	// We don't write the first two entries, which are always black and white.
	std::copy(m_colors.begin() + 2, m_colors.end(), std::ostream_iterator<CDXColor>(sink.os));
	sink.os << GetTextEOL() << "</" << kCDXML_colortable << ">";
}
