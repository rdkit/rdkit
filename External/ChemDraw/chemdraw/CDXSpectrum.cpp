// CommonCS/LibCommon/Src/CDXSpectrum.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright (c) 1986-2009, CambridgeSoft Corp., All Rights Reserved

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
#include "cs_auto_buffer.h"
#include "XMLPutEnum.h"
#include <stdlib.h>	// for atoi and friends
#include <math.h>	// for floor
#include <stdio.h>	// for sprintf
#include <sstream>
#include <cmath>

static std::string lCharDataBuffer;

XMLPutEnum<CDXSpectrumYType>::Value sXMLSpectrumYTypeValues[] = {
	{kCDXSpectrumYType_Unknown,				"Unknown"},
	{kCDXSpectrumYType_Absorbance,			"Absorbance"},
	{kCDXSpectrumYType_Transmittance,		"Transmittance"},
	{kCDXSpectrumYType_PercentTransmittance,"PercentTransmittance"},
	{kCDXSpectrumYType_Other,				"Other"},
	{kCDXSpectrumYType_ArbitraryUnits,		"ArbitraryUnits"}
};

XMLPutEnum<CDXSpectrumYType> sXMLSpectrumYType(sXMLSpectrumYTypeValues, sizeof sXMLSpectrumYTypeValues, kCDXSpectrumYType_Unknown);
void XMLPut(XMLDataSink &sink, CDXSpectrumYType  v) {	sink.os << sXMLSpectrumYType.lookup(v); }

XMLPutEnum<CDXSpectrumXType>::Value sXMLSpectrumXTypeValues[] = {
	{kCDXSpectrumXType_Unknown,			"Unknown"},
	{kCDXSpectrumXType_Wavenumbers,		"Wavenumbers"},
	{kCDXSpectrumXType_Microns,			"Microns"},
	{kCDXSpectrumXType_Hertz,			"Hertz"},
	{kCDXSpectrumXType_MassUnits,		"MassUnits"},
	{kCDXSpectrumXType_PartsPerMillion,	"PartsPerMillion"},
	{kCDXSpectrumXType_Other,			"Other"}
};

XMLPutEnum<CDXSpectrumXType> sXMLSpectrumXType(sXMLSpectrumXTypeValues, sizeof sXMLSpectrumXTypeValues, kCDXSpectrumXType_Unknown);
void XMLPut(XMLDataSink &sink, CDXSpectrumXType  v) {	sink.os << sXMLSpectrumXType.lookup(v); }

XMLPutEnum<CDXSpectrumClass>::Value sXMLSpectrumClassValues[] = {
	{kCDXSpectrumClass_Unknown,			"Unknown"},
	{kCDXSpectrumClass_Chromatogram,	"Chromatogram"},
	{kCDXSpectrumClass_Infrared,		"Infrared"},
	{kCDXSpectrumClass_UVVis,			"UVVis"},
	{kCDXSpectrumClass_XRayDiffraction,	"XRayDiffraction"},
	{kCDXSpectrumClass_MassSpectrum,	"MassSpectrum"},
	{kCDXSpectrumClass_NMR,				"NMR"},
	{kCDXSpectrumClass_Raman,			"Raman"},
	{kCDXSpectrumClass_Fluorescence,	"Fluorescence"},
	{kCDXSpectrumClass_Atomic,			"Atomic"}
};

XMLPutEnum<CDXSpectrumClass> sXMLSpectrumClass(sXMLSpectrumClassValues, sizeof sXMLSpectrumClassValues, kCDXSpectrumClass_Unknown);
void XMLPut(XMLDataSink &sink, CDXSpectrumClass  v) {	sink.os << sXMLSpectrumClass.lookup(v); }

// ***********************
// ** class CDXSpectrum **
// ***********************
//
// Specialization of CDXObject for CDXSpectrum objects

CDXSpectrum::CDXSpectrum(CDXObjectID id_arg)
	:	CDXGraphicObject(kCDXObj_Spectrum, id_arg)
	,	m_xSpacing (0)
	,	m_xLow (0)
	,	m_xType (kCDXSpectrumXType_Unknown)
	,	m_yType (kCDXSpectrumYType_Unknown)
	,	m_class (kCDXSpectrumClass_Unknown)
	,	m_yLow (0.0)
	,	m_yScale (0.0)
{
}

CDXSpectrum::CDXSpectrum (const CDXSpectrum &src)
	:	CDXGraphicObject(src)
	,	m_dataPoints (src.m_dataPoints)
	,	m_xSpacing (src.m_xSpacing)
	,	m_xLow (src.m_xLow)
	,	m_xType (src.m_xType)
	,	m_yType (src.m_yType)
	,	m_class (src.m_class)
	,	m_yLow (src.m_yLow)
	,	m_yScale (src.m_yScale)
{
}

CDXSpectrum::~CDXSpectrum()
{
}

CDXObject*	CDXSpectrum::Clone() const
{
	return new CDXSpectrum (*this);
}

std::string CDXSpectrum::XMLObjectName() const
{
	return kCDXML_spectrum;
}


// Take a CDX input stream and create a CDX structure
void CDXSpectrum::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	FLOAT64 f64;
	switch (attribTag_arg)
	{
	case kCDXProp_Spectrum_XSpacing:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_xSpacing = SwapBytes(f64);
		break;

	case kCDXProp_Spectrum_XLow:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_xLow = SwapBytes(f64);
		break;

	case kCDXProp_Spectrum_XType:
		m_xType=(CDXSpectrumXType)src_arg.GetUINT16();
		break;

	case kCDXProp_Spectrum_YType:
		m_yType=(CDXSpectrumYType)src_arg.GetUINT16();
		break;

	case kCDXProp_Spectrum_Class:
		m_class=(CDXSpectrumClass)src_arg.GetUINT16();
		break;

	case kCDXProp_Spectrum_XAxisLabel:
		m_xAxisLabel.Read(src_arg, size_arg);
		break;

	case kCDXProp_Spectrum_YAxisLabel:
		m_yAxisLabel.Read(src_arg, size_arg);
		break;

	case kCDXProp_Spectrum_YLow:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_yLow = SwapBytes(f64);
		break;

	case kCDXProp_Spectrum_YScale:
		src_arg.GetBytes((char *)&f64, sizeof(f64));
		m_yScale = SwapBytes(f64);
		break;

	case kCDXProp_Spectrum_DataPoint:
	{
		size_t numDataPoints = size_arg / sizeof(FLOAT64);
		m_dataPoints.resize(numDataPoints);
		for (int i=0; i<numDataPoints; i++)
		{
			src_arg.GetBytes((char *)&f64, sizeof(f64));
			m_dataPoints[i] = SwapBytes(f64);
		}
		
		if (m_yScale == 0)
		{
			// For writing XML, we scale the output data to a 10000:1 dynamic range.  This is plenty for imaging purposes (though
			// not good enough for spectral storage for later use as spectral data.)

			// Find the dynamic range of the spectrum.
			m_yLow=m_dataPoints[0];
			double yHigh=m_dataPoints[0];
			for (int i=1; i<m_dataPoints.size(); i++)
			{
				if (m_dataPoints[i] < m_yLow)
					m_yLow=m_dataPoints[i];
				if (m_dataPoints[i] > yHigh)
					yHigh=m_dataPoints[i];
			}

			// Scale the data as we put it out.
			m_yScale=(yHigh-m_yLow)/10000.0;
		}

		break;
	}

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
	lCharDataBuffer.erase();
}


// Read CDXML and create a CDX structure
void CDXSpectrum::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_Spectrum_XSpacing:
		m_xSpacing = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_Spectrum_XLow:
		m_xLow = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_Spectrum_XType:
		m_xType = sXMLSpectrumXType.lookup(attribValue_arg);
		break;

	case kCDXProp_Spectrum_YType:
		m_yType = sXMLSpectrumYType.lookup(attribValue_arg);
		break;

	case kCDXProp_Spectrum_Class:
		m_class = sXMLSpectrumClass.lookup(attribValue_arg);
		break;

	case kCDXProp_Spectrum_XAxisLabel:
		m_xAxisLabel=CDXString(UnicodeToString(attribValue_arg));
		break;

	case kCDXProp_Spectrum_YAxisLabel:
		m_yAxisLabel=CDXString(UnicodeToString(attribValue_arg));
		break;

	case kCDXProp_Spectrum_YLow:
		m_yLow = atof(attribValue_arg.c_str());
		break;

	case kCDXProp_Spectrum_YScale:
		m_yScale = atof(attribValue_arg.c_str());
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
	lCharDataBuffer.erase();
}

// Take a CDX structure and create a CDX output stream.
void CDXSpectrum::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

	sink_arg.PutAttribute( kCDXProp_Spectrum_XSpacing, FLOAT64(m_xSpacing));

	sink_arg.PutAttribute( kCDXProp_Spectrum_XLow, FLOAT64(m_xLow));

	if (m_xType != kCDXSpectrumXType_Unknown)
		sink_arg.PutAttribute(kCDXProp_Spectrum_XType, UINT16(m_xType));

	if (m_yType != kCDXSpectrumYType_Unknown)
		sink_arg.PutAttribute(kCDXProp_Spectrum_YType, UINT16(m_yType));

	if (m_class != kCDXSpectrumClass_Unknown)
		sink_arg.PutAttribute(kCDXProp_Spectrum_Class, UINT16(m_class));

	if (!m_xAxisLabel.empty())
		sink_arg.PutAttributeCDXString( kCDXProp_Spectrum_XAxisLabel, m_xAxisLabel);

	if (!m_yAxisLabel.empty())
		sink_arg.PutAttributeCDXString( kCDXProp_Spectrum_YAxisLabel, m_yAxisLabel);

	if (m_yLow != 0)
		sink_arg.PutAttribute( kCDXProp_Spectrum_YLow, FLOAT64(m_yLow));

	if (m_yScale != 0)
		sink_arg.PutAttribute( kCDXProp_Spectrum_YScale, FLOAT64(m_yScale));

	if (m_dataPoints.size() != 0)
	{
		int floatSize = sizeof(FLOAT64);
		auto_buffer<char> ab(m_dataPoints.size() * floatSize);
		char *p = ab.ptr();
		for (std::vector<double>::const_iterator i = m_dataPoints.begin(); i != m_dataPoints.end(); ++i, p += floatSize)
			SwapBytes(&(*i), p);
		sink_arg.PutAttribute(kCDXProp_Spectrum_DataPoint, (const INT8 *)ab.ptr(), ab.length());
	}
}


// Take a CDX structure and write a CDXML output stream
void CDXSpectrum::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

	streamsize oldPrecision=sink_arg.os.precision(8);
	CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_XSpacing, m_xSpacing);

	CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_XLow, m_xLow);
	sink_arg.os.precision(oldPrecision);

	CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_XType, m_xType);

	CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_YType, m_yType);

	CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_Class, m_class);

	if (m_xAxisLabel.length() != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_XAxisLabel, StringToUnicode("", m_xAxisLabel.str()));

	if (m_yAxisLabel.length() != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_YAxisLabel, StringToUnicode("", m_yAxisLabel.str()));

	if (m_yLow != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_YLow, m_yLow);

	if (m_yScale != 0)
		CDXMLPutAttribute(sink_arg, kCDXProp_Spectrum_YScale, m_yScale);

}

// CDXSpectrum::XMLNeedToWriteContent
//
// Override this to indicate that we need to write content between the start-tag and end-tag.
bool CDXSpectrum::XMLNeedToWriteContent() const
{
	return true;
}

// CDXSpectrum::XMLWriteContent
//
// Override this to write content between the start-tag and end-tag.
void CDXSpectrum::XMLWriteContent(XMLDataSink &sink_arg) const
{
	// If yScale is 0 or 1, write out unscaled floating point data
	if (m_yScale == 0 || m_yScale == 1)
	{
		for (int i=0; i<m_dataPoints.size(); i++)
			sink_arg.os << (i%8 == 0 ? GetTextEOL() : " ") << m_dataPoints[i]-m_yLow;
	}
	else
	{
		// This will convert [0,0.5) to 0, [0.5,1.5) to 1, ..., [9999.5,10000] to 10000.
		for (int i=0; i<m_dataPoints.size(); i+=16)
		{
			char buf[1024];
			if (m_dataPoints.size()-i > 15)
			{
				int v[16];
				for (int j=0; j<16; j++)
					v[j]=(int) floor(0.5 + (m_dataPoints[i+j]-m_yLow)/m_yScale);
				snprintf(buf, 1024, "%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i",
							v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15]);
				sink_arg.os << buf << GetTextEOL();
			}
			else
			{
				for (int j = 0; j < m_dataPoints.size() - i; j++)
				{
					int theVal = (int) floor(0.5 + (m_dataPoints[i+j]-m_yLow)/m_yScale);
					sink_arg.os << theVal;
				}
				sink_arg.os << GetTextEOL();
			}
		}
	}

/*
	// If yScale is 0 or 1, write out unscaled floating point data
	if (m_yScale == 0 || m_yScale == 1)
	{
		for (int i=0; i<m_dataPoints.size(); i++)
			sink_arg.os << (i%8 == 0 ? GetTextEOL() : " ") << m_dataPoints[i]-m_yLow;
	}
	else
	{
		// This will convert [0,0.5) to 0, [0.5,1.5) to 1, ..., [9999.5,10000] to 10000.
		for (int i=0; i<m_dataPoints.size(); i++)
			sink_arg.os << (i%16 == 0 ? GetTextEOL() : " ") << (int) floor(0.5 + (m_dataPoints[i]-m_yLow)/m_yScale);
	}
*/
}

void CDXSpectrum::XMLProcessTokens(const std::string &charData_arg)
{
	size_t origSize = m_dataPoints.size();
	std::istringstream is(charData_arg);
	std::copy(std::istream_iterator<double>(is), std::istream_iterator<double>(), std::back_inserter(m_dataPoints));
	if (m_yScale != 0 && m_yScale != 1)
		for (size_t i=origSize; i<m_dataPoints.size(); i++)
			m_dataPoints[i]*=m_yScale;
	if (m_yLow)
		for (size_t i=origSize; i<m_dataPoints.size(); i++)
			m_dataPoints[i]+=m_yLow;
}

void CDXSpectrum::XMLStoreCharacterData(const std::string &charData_arg)
{
	lCharDataBuffer += charData_arg;
	size_t lastSpace = lCharDataBuffer.find_last_of(' ');
	if (lastSpace != std::string::npos)
	{
		XMLProcessTokens(lCharDataBuffer.substr(0,lastSpace));
		lCharDataBuffer.erase(0,lastSpace);
	}
}

void CDXSpectrum::FinishReading()
{
	if (lCharDataBuffer.size() != 0)
	{
		XMLProcessTokens(lCharDataBuffer);
		lCharDataBuffer.erase();
	}
}
