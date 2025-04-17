// CommonCS/LibCommon/Hdr/CDXStdObjects.h
// Contains: Program-independent class library for managing CDX objects
//           This file describes the interface for the standard subclasses of CDXObject.
// Copyright 1986-2009, CambridgeSoft Corp., All Rights Reserved

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

#pragma once

#if TARGET_OS_WIN32
// Don't warn about truncation of long symbols
#pragma warning (disable:4786)
#endif

#include "CoreChemistryAPI.h"
#include "CDXObject.h"
#include "CDXFontTable.h"
#include <set>
#include <string>
#include <vector>
#include <iterator>
#include <stdint.h>

using std::set;
using std::vector;
using std::pair;
// Classes defined in this file:
class CDXGraphicObject;			// abstract class for visible objects which can have position, bounds, ordering
class CDXDocument;				// subclasses for each type of object found in CDX files
class CDXGroup;
class CDXFragment;
class CDXNode;
class CDXBond;
class CDXText;
class CDXGraphic;
class CDXCurve;
class CDXEmbeddedObject;
class CDXNamedAlternativeGroup;
class CDXTemplateGrid;
class CDXObjectDefinition;
class CDXSpectrum;
class CDXObjectTagDictionaryEntry;
class CDXObjectTag;
class CDXDocumentPropertyCollection;
class CDXDocumentProperty;

// Add factory functions for all the objects defined in this file
CORE_CHEMISTRY_API void CDXAddStandardFactories(CDXObjectFactory &);

// Don't warn about unknown pragmas
#if TARGET_OS_WIN32
#pragma warning (disable:4068)
#endif

typedef UINT16 CDXColorIndex;

// ConnectBonds
// This defines the m_beginNode (etc) derived members from the m_beginNodeID
// values that are stored in the file.  It is normally called by a FinishReading() method.
CORE_CHEMISTRY_API void ConnectBonds(CDXObject *obj);

CORE_CHEMISTRY_API std::string MakeStringSafe(const std::string& inInternalUTF8String);

// *******************************
// ***** class CDXColorTable *****
// *******************************
//
// Helper class for manipulating CDX color tables

struct CORE_CHEMISTRY_API CDXColor
{
	UINT16	m_red,
			m_green,
			m_blue;
			CDXColor() : m_red (0), m_green (0), m_blue (0)	{}
			CDXColor (UINT16 red, UINT16 green, UINT16 blue) : m_red (red), m_green (green), m_blue (blue)	{}
	bool	operator == (const CDXColor &op2) const
			{ return m_red == op2.m_red  &&  m_green == op2.m_green  &&  m_blue == op2.m_blue; }

	static CDXColor	White()	{ return CDXColor ((UINT16)-1,(UINT16)-1,(UINT16)-1); }
	static CDXColor	Black()	{ return CDXColor ((UINT16) 0,(UINT16) 0,(UINT16) 0); }
	static CDXColor	Red()	{ return CDXColor ((UINT16)-1,(UINT16) 0,(UINT16) 0); }
	static CDXColor	Green()	{ return CDXColor ((UINT16) 0,(UINT16)-1,(UINT16) 0); }
	static CDXColor	Blue()	{ return CDXColor ((UINT16) 0,(UINT16) 0,(UINT16)-1); }
	static CDXColor	Yellow(){ return CDXColor ((UINT16)-1,(UINT16)-1,(UINT16) 0); }
	static CDXColor	Cyan()	{ return CDXColor ((UINT16) 0,(UINT16)-1,(UINT16)-1); }
	static CDXColor	Magenta(){return CDXColor ((UINT16)-1,(UINT16) 0,(UINT16)-1); }
};


class CORE_CHEMISTRY_API CDXColorTable
{
public:
    std::vector<CDXColor> m_colors;
        // The first two entries are always black and white.
        // The next two are the background and foreground colors.

    // Constructor
    CDXColorTable();

    // Read from binary data stream
    void Read(CDXDataSource& src_arg, size_t size_arg);

    // Write to binary data stream
    void Write(CDXDataSink& sink_arg) const;

    // Read from XML attribute string
    void XMLRead(const std::string&);

    // Write as XML attribute string
    void XMLWrite(XMLDataSink&) const;

    // Return the current size of the color table
	CDXColorIndex NumColors() const;

    // Return the resulting size of the color table
	CDXColorIndex AddColor(const CDXColor& color);

    // Add a color to the color table, and return its index
	CDXColorIndex FindOrAddColor(const CDXColor& color);

    // Get a color from the table
    const CDXColor&	GetColor(CDXColorIndex index) const;

    // Return the index of the first instance of the given color
	CDXColorIndex GetColorsIndex(const CDXColor& color, bool skipHardcoded = false, bool getNearest = false) const;
};


// ***********************
// ** class CDXDocument **
// ***********************
//
// Specialization of CDXObject for CDXDocument objects

class CORE_CHEMISTRY_API CDXDocument : public CDXObject
{
public:

	CDXString		m_creationProgram;
	CDXString		m_name;
	CDXString		m_comment;
	std::string		m_cartridgeData;
	CDXRectangle	m_boundingBox;
	Point32			m_windowPosition;
	Point32			m_windowSize;
	bool			m_windowIsZoomed;
	bool			m_fractionalWidths;
	bool			m_interpretChemically;
	bool			m_showAtomQuery;
	bool			m_showAtomStereo;
	bool			m_showAtomEnhancedStereo;
	bool			m_showAtomNumber;
	bool			m_showResidueID;
	bool			m_showBondQuery;
	bool			m_showBondRxn;
	bool			m_showBondStereo;
	bool			m_showTerminalCarbons;
	bool			m_showNonTerminalCarbons;
	bool			m_hideImplicitHydrogens;
	UINT16			m_magnification;
	CDXFontStyle	m_labelStyle;
	CDXFontStyle	m_captionStyle;
	CDXCoordinate	m_hashSpacing;
	CDXCoordinate	m_marginWidth;
	CDXCoordinate	m_lineWidth;
	CDXCoordinate	m_boldWidth;
	CDXCoordinate	m_bondLength;
	UINT16			m_bondSpacing;
	CDXCoordinate	m_bondSpacingAbs;
	CDXBondSpacingType		m_bondSpacingType;
	UINT32			m_chainAngle;
	CDXCoordinate	m_labelLineHeight;
	CDXCoordinate	m_captionLineHeight;
	CDXTextJustification	m_labelJustification;
	CDXTextJustification	m_captionJustification;
	CDXAminoAcidTermini		m_aminoAcidTermini;
	bool			m_showSequenceTermini;
	bool			m_showSequenceBonds;
	bool			m_showSequenceUnlinkedBranches;
	UINT16			m_residueWrapCount;
	UINT16			m_residueBlockCount;
	CDXRectangle	m_printMargins;
	CDXPoint2D		m_fixInPlaceExtent;
	CDXPoint2D		m_fixInPlaceGap;
	CDXColorTable	m_colorTable;
    CDXColorIndex	m_fgColor;
    CDXColorIndex	m_bgColor;
	UINT16			m_fgAlpha;
	UINT16			m_bgAlpha;
	CDXFontTable	m_fontTable;
	char			m_macPrintRecord[120];
	UINT32			m_flags;
	UINT32			m_flags2;
	vector<CDXChemProp> m_chemicalProperties;
    CDXLineType     m_lineType;
    CDXFillType     m_fillType;
    CDXRxnAutonumberStyle m_rxnAutonumberStyle;
    UINT16          m_rxnAutonumberStart;
    bool            m_rxnAutonumberConditions;
    CDXString       m_rxnAutonumberFormat;
    CDXString       m_monomerRenderingStyle;

public:

	// This is a set of bits set for attributes which have been given values.
	// Note that some attributes have defaults; they are initialized with
	// the defaults and always have values.
	// These are the ones which don't have defaults:
	enum CORE_CHEMISTRY_API CDXDocumentProperty1 
    {
		has_boundingBox			= 0x00000002,
		has_windowPosition		= 0x00000004,
		has_windowSize			= 0x00000008,
		has_windowIsZoomed		= 0x00000010,
		has_fractionalWidths	= 0x00000020,
		has_magnification		= 0x00000040,
		has_labelStyle			= 0x00000080,
		has_captionStyle		= 0x00000100,
		has_hashSpacing			= 0x00000200,
		has_marginWidth			= 0x00000400,
		has_lineWidth			= 0x00000800,
		has_boldWidth			= 0x00001000,
		has_bondLength			= 0x00002000,
		has_bondSpacing			= 0x00004000,
		has_bondSpacingAbs		= 0x00008000,
		has_chainAngle			= 0x00010000,
		has_labelLineHeight		= 0x00020000,
		has_captionLineHeight	= 0x00040000,
		has_labelJustification  = 0x00080000,
		has_captionJustification= 0x00100000,
		has_printMargins		= 0x00200000,
		has_printRecord			= 0x00400000,
		has_showAtomQuery		= 0x00800000,
		has_showAtomStereo		= 0x01000000,
		has_showAtomNumber		= 0x02000000,
		has_showBondQuery		= 0x04000000,
		has_showBondStereo		= 0x08000000,
		has_fixInPlaceExtent	= 0x10000000,
		has_showBondRxn			= 0x20000000,
		has_interpretChemically = 0x40000000,
		has_fixInPlaceGap		= 0x80000000
	};

	enum CORE_CHEMISTRY_API CDXDocumentProperty2 {
		has_showTerminalCarbons	= 0x00000002,
		has_showNonTerminalCarbons	= 0x00000004,
		has_hideImplicitHydrogens	= 0x00000008,
		has_showAtomEnhancedStereo	= 0x00000010,
		has_fgColor				= 0x00000020,
		has_bgColor				= 0x00000040,
		has_fgAlpha				= 0x00000080,
		has_bgAlpha				= 0x00000100,
		has_aminoAcidTermini	= 0x00000200,
		has_showSequenceTermini	= 0x00000400,
		has_showSequenceBonds	= 0x00000800,
		has_residueWrapCount	= 0x00001000,
		has_residueBlockCount	= 0x00002000,
		has_showResidueID		= 0x00010000,
        has_lineStyle           = 0x00020000,
        has_fillStyle           = 0x00040000,
        has_rxnAutonumberStyle  = 0x00080000,
        has_rxnAutonumberStart  = 0x00100000,
        has_rxnAutonumberConditions = 0x00200000,
        has_rxnAutonumberFormat = 0x00400000,
		has_showSequenceUnlinkedBranches = 0x00800000,
        
	};

protected:

	CDXTag	GetTag() const override	{ return kCDXObj_Document; }
	void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg) override;
	void XMLStoreAttribute(CDXTag, const std::string &) override;
	void WriteAttributesTo(CDXDataSink &sink) const override;
	std::string XMLObjectName() const override;
	void XMLWriteAttributes(XMLDataSink &) const override;
	void XMLWriteContent(XMLDataSink &) const override;
	void XMLWriteContentDone(XMLDataSink &) const override;
	bool XMLNeedToWriteContent() const override;
    void FinishReading() override;

private:
	
    CDXDocument& operator=(const CDXDocument &a); /* assignment is not implemented */

public:

	CDXDocument(CDXObjectID id);
	CDXDocument(const CDXDocument &a);
	virtual	~CDXDocument();
	CDXObject*	Clone()	const override;

	bool Known(CDXDocumentProperty1 f)	const { return (m_flags & f) != 0; }
	void Known(CDXDocumentProperty1 f, bool val){ if (val) m_flags |= f;  else m_flags &= ~f; }
	bool Known(CDXDocumentProperty2 f)	const { return (m_flags2 & f) != 0; }
	void Known(CDXDocumentProperty2 f, bool val){ if (val) m_flags2 |= f;  else m_flags2 &= ~f; }
	
	CDXString CreationProgram() const { return m_creationProgram; }
	CDXString Name()			const { return m_name; }
	CDXString Comment()			const { return m_comment; }
	CDXRectangle BoundingBox()	const { return m_boundingBox; }
	
	void CreationProgram(const CDXString &s)	{ m_creationProgram = s; }
	void Name			(const CDXString &s)	{ m_name = s; }
	void Comment		(const CDXString &s)	{ m_comment = s; }
	void BoundingBox	(const CDXRectangle &r)	{ m_flags |= has_boundingBox; m_boundingBox = r; }

	CDXFontTableEntry LookupFont(int family) const	{ return m_fontTable.LookupFont(family); }

	size_t NumChemicalProperties() const { return m_chemicalProperties.size(); }
	CDXChemProp ChemicalProperty(size_t i) const { return m_chemicalProperties[i]; }
	void ChemicalProperty(const CDXChemProp &p) { m_chemicalProperties.push_back(p); }

    const CDXDocumentPropertyCollection* GetDocumentProperties() const;
    CDXDocumentPropertyCollection* GetDocumentProperties();
    
    void NormalizeTextEncoding();
};

/*
+===========================================================================================+
|								CDXDocument HELPER FUNCTIONS								|
+===========================================================================================+
*/
CORE_CHEMISTRY_API CDXDocument*	CDXReadDocFromStorage(CDXDataSource &src, bool doThrow = true);
CORE_CHEMISTRY_API CDXDocument*	CDXReadDocFromStorage(CDXDataSource &src, std::string& errMsg);
CORE_CHEMISTRY_API void			CDXWriteDocToStorage (const CDXDocument* pDoc, CDXDataSink &sink);
CORE_CHEMISTRY_API CDXDocument*	ReadCdxDiskFile		 (const char *filename, bool doThrow = true);
CORE_CHEMISTRY_API bool			WriteCdxDiskFile	 (const CDXDocument *pDoc, const char *filename, std::string* pErrorMsg = NULL);
CORE_CHEMISTRY_API bool			WriteCdxmlDiskFile	 (const CDXDocument *pDoc, const char *filename, std::string* pErrorMsg = NULL);



// ****************************
// ** class CDXGraphicObject **
// ****************************
//
// Intermediate abstract class for visible objects.  These objects typically have
// position, bounds, and front-to-back ordering.

class CORE_CHEMISTRY_API CDXGraphicObject : public CDXObject
{
public:
	CDXObjectID 	m_supersededBy;		// This ID of an object that replaces this one
	CDXPoint2D		m_2dPosition;
	CDXPoint3D		m_3dPosition;
	CDXRectangle	m_boundingBox;
	UINT16			m_zOrder;
	CDXColorIndex	m_fgColor;
	CDXColorIndex	m_bgColor;
	CDXColorIndex   highlightColor;
	UINT16			m_fgAlpha;
	UINT16			m_bgAlpha;
	bool			m_visible;
	std::vector<CDXPropRep> *m_propertiesRepresented;
	std::vector<CDXChemProp> *m_chemicalProperties;
	bool			m_ignoreErrors;
	CDXString		*m_errorDescription;
	bool			m_fractionalWidths;
	bool			m_interpretChemically;
	bool			m_showAtomQuery;
	bool			m_showAtomStereo;
	bool			m_showAtomEnhancedStereo;
	bool			m_showAtomNumber;
	bool			m_showResidueID;
	bool			m_showBondQuery;
	bool			m_showBondRxn;
	bool			m_showBondStereo;
	bool			m_showTerminalCarbons;
	bool			m_showNonTerminalCarbons;
	bool			m_hideImplicitHydrogens;
	CDXFontStyle	m_labelStyle;
	CDXFontStyle	m_captionStyle;
	CDXCoordinate	m_hashSpacing;
	CDXCoordinate	m_marginWidth;
	CDXCoordinate	m_lineWidth;
	CDXCoordinate	m_boldWidth;
	CDXCoordinate	m_bondLength;
	UINT16			m_bondSpacing;
	CDXCoordinate	m_bondSpacingAbs;
	CDXBondSpacingType		m_bondSpacingType;
	UINT32			m_chainAngle;
	CDXCoordinate	m_labelLineHeight;
	CDXCoordinate	m_captionLineHeight;
	CDXTextJustification	m_labelJustification;
	CDXTextJustification	m_captionJustification;
	CDXAminoAcidTermini		m_aminoAcidTermini;
	bool			m_showSequenceTermini;
	bool			m_showSequenceBonds;
	bool			m_showSequenceUnlinkedBranches;
	UINT16			m_residueWrapCount;
	UINT16			m_residueBlockCount;
	CDXLineType		m_lineType;
	CDXFillType		m_fillType;
	UINT32			m_flags;
	UINT32			m_flags2;
    CDXRxnAutonumberStyle m_rxnAutonumberStyle;
    UINT16          m_rxnAutonumberStart;
    bool            m_rxnAutonumberConditions;
    CDXString       m_rxnAutonumberFormat;
    bool            showPerspective;
    
		// This is a set of bits set for attributes which have been given values.
		// Note that some attributes have defaults; they are initialized with
		// the defaults and always have values.
		// These are the ones which don't have defaults:
public:
	enum CDXGraphicObjectProperty1 {
		has_2dPosition			= 0x00000001,
		has_3dPosition			= 0x00000002,
		has_boundingBox			= 0x00000004,
		has_zOrder				= 0x00000008,
		has_fgColor				= 0x00000010,
		has_bgColor				= 0x00000020,
		has_fractionalWidths	= 0x00000040,
		has_labelStyle			= 0x00000080,
		has_captionStyle		= 0x00000100,
		has_hashSpacing			= 0x00000200,
		has_marginWidth			= 0x00000400,
		has_lineWidth			= 0x00000800,
		has_boldWidth			= 0x00001000,
		has_bondLength			= 0x00002000,
		has_bondSpacing			= 0x00004000,
		has_bondSpacingAbs		= 0x00008000,
		has_chainAngle			= 0x00010000,
		has_labelLineHeight		= 0x00020000,
		has_captionLineHeight	= 0x00040000,
		has_labelJustification  = 0x00080000,
		has_captionJustification= 0x00100000,
		has_showAtomQuery		= 0x00200000,
		has_showAtomStereo		= 0x00400000,
		has_showAtomNumber		= 0x00800000,
		has_showBondQuery		= 0x01000000,
		has_showBondStereo		= 0x02000000,
		has_showBondRxn			= 0x04000000,
		has_interpretChemically = 0x08000000,
		has_showTerminalCarbons = 0x10000000,
		has_showNonTerminalCarbons = 0x20000000,
		has_hideImplicitHydrogens = 0x40000000,
		has_showAtomEnhancedStereo = 0x80000000
	};
	enum CDXGraphicObjectProperty2 {
		has_fgAlpha				= 0x00000001,
		has_bgAlpha				= 0x00000002,
		has_aminoAcidTermini	= 0x00000004,
		has_showSequenceTermini	= 0x00000008,
		has_showSequenceBonds	= 0x00000010,
		has_residueWrapCount	= 0x00000020,
		has_residueBlockCount	= 0x00000040,
		has_showResidueID		= 0x00000200,
        has_lineStyle           = 0x00000400,
        has_fillStyle           = 0x00000800,
        has_rxnAutonumberConditions = 0x00001000,
        has_rxnAutonumberStart  = 0x00002000,
        has_rxnAutonumberStyle  = 0x00004000,
        has_rxnAutonumberFormat = 0x00008000,
		has_showSequenceUnlinkedBranches = 0x000010000,
		hasHighlightColor       = 0x000020000
	};

private:
	CDXGraphicObject& operator=(const CDXGraphicObject &a); /* assignment is not implemented */

public:
	CDXGraphicObject(CDXTag tag, CDXObjectID id);
	CDXGraphicObject(const CDXGraphicObject &a);
	virtual	~CDXGraphicObject() = 0;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual bool XMLNeedToWriteContent() const;
	virtual void XMLWriteContent(XMLDataSink &sink_arg) const;

	bool Known(CDXGraphicObjectProperty1 f)				const	{ return (m_flags & f) != 0; }
	void Known(CDXGraphicObjectProperty1 f, bool val)			{ if (val) m_flags |= f;  else m_flags &= ~f; }
	bool Known(CDXGraphicObjectProperty2 f)				const	{ return (m_flags2 & f) != 0; }
	void Known(CDXGraphicObjectProperty2 f, bool val)			{ if (val) m_flags2 |= f;  else m_flags2 &= ~f; }
	bool KnownPosition()				const	{ return (m_flags & has_2dPosition ) != 0; }
	bool KnownPosition3D()				const	{ return (m_flags & has_3dPosition ) != 0; }
	bool KnownBoundingBox()				const	{ return (m_flags & has_boundingBox) != 0; }
	bool KnownZOrder()					const	{ return (m_flags & has_zOrder     ) != 0; }
	bool KnownFGColor()					const	{ return (m_flags & has_fgColor    ) != 0; }
	bool KnownBGColor()					const	{ return (m_flags & has_bgColor    ) != 0; }
	bool KnownFGAlpha()					const	{ return (m_flags2 & has_fgAlpha    ) != 0; }
	bool KnownBGAlpha()					const	{ return (m_flags2 & has_bgAlpha    ) != 0; }
	bool KnownShowAtomStereo()			const	{ return (m_flags & has_showAtomStereo) != 0; }
	bool KnownShowBondStereo()			const	{ return (m_flags & has_showBondStereo) != 0; }
	bool HasHighlightColor()			const;

	const CDXPoint2D&	Position()		const	{ assert (KnownPosition());		return m_2dPosition;  }
	const CDXPoint3D&	Position3D()	const	{ assert (KnownPosition3D());	return m_3dPosition;  }
	const CDXRectangle&	BoundingBox()	const	{ assert (KnownBoundingBox());	return m_boundingBox; }
	UINT32				ZOrder()		const	{ assert (KnownZOrder());		return m_zOrder; }
	CDXColorIndex		FGColor()		const	{ assert (KnownFGColor());		return m_fgColor; }
	CDXColorIndex		BGColor()		const	{ assert (KnownBGColor());		return m_bgColor; }
	bool				ShowAtomStereo()const	{ assert (KnownShowAtomStereo());	return m_showAtomStereo; }
	bool				ShowBondStereo()const	{ assert (KnownShowBondStereo());	return m_showAtomStereo; }

	void Position	(CDXPoint2D   p)			{ m_flags |= has_2dPosition;  m_2dPosition  = p; }
	void Position3D	(CDXPoint3D   p)			{ m_flags |= has_3dPosition;  m_3dPosition  = p; }
	void BoundingBox(CDXRectangle r)			{ m_flags |= has_boundingBox; m_boundingBox = r; }
	void ZOrder		(UINT32       z)			{ m_flags |= has_zOrder;      m_zOrder      = z; }
	void FGColor	(CDXColorIndex  colorIndex)	{ m_flags |= has_fgColor;     m_fgColor     = colorIndex; }
	void BGColor	(CDXColorIndex  colorIndex)	{ m_flags |= has_bgColor;     m_bgColor     = colorIndex; }
	void ShowAtomStereo	(bool	  s)			{ m_flags |= has_showAtomStereo; m_showAtomStereo = s; }
	void ShowBondStereo	(bool	  s)			{ m_flags |= has_showBondStereo; m_showAtomStereo = s; }
    
	void SetHighlightColor(CDXColorIndex  colorIndex);
	CDXColorIndex GetHighlightColor() const;
    
    void SetShowPerspective(bool inShowPerspective);
    bool GetShowPerspective() const;

	size_t NumPropertiesRepresented() const { return ((m_propertiesRepresented == NULL) ? 0 : m_propertiesRepresented->size()); }
	const CDXPropRep &PropertyRepresented(size_t i) const { if (m_propertiesRepresented != NULL && m_propertiesRepresented->size() > i) return (*m_propertiesRepresented)[i]; static const CDXPropRep emptyProp; return emptyProp; }
	void PropertyRepresented(const CDXPropRep &p) { if (m_propertiesRepresented == NULL) m_propertiesRepresented = new std::vector<CDXPropRep>; m_propertiesRepresented->push_back(p); }

	size_t NumChemicalProperties() const { return ((m_chemicalProperties == NULL) ? 0 : m_chemicalProperties->size()); }
	const CDXChemProp &ChemicalProperty(size_t i) const { if (m_chemicalProperties != NULL && m_chemicalProperties->size() > i) return (*m_chemicalProperties)[i]; static const CDXChemProp emptyProp; return emptyProp;  }
	void ChemicalProperty(const CDXChemProp &p) { if (m_chemicalProperties == NULL) m_chemicalProperties = new std::vector<CDXChemProp>; m_chemicalProperties->push_back(p); }

	// Figure out the bounds of this object and all its children
	virtual CDXRectangle *BoundsIncludingChildren() const;
};


// *******************
// ** class CDXPage **
// *******************
//
// Specialization of CDXObject for CDXPage objects

class CORE_CHEMISTRY_API CDXPage : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Page; }
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
    virtual void FinishReading();

public:
	CDXRectangle			m_boundingBox;
	UINT16					m_bgColor;
	CDXCoordinate			m_pageWidth;
	CDXCoordinate			m_pageHeight;
	CDXCoordinate			m_headerPosition;
	CDXCoordinate			m_footerPosition;
	CDXCoordinate			m_pageOverlap;
	CDXString				m_header;
	CDXString				m_footer;
	bool					m_printTrimMarks;
	CDXDrawingSpaceType		m_drawingSpaceType;
	UINT16					m_heightPages;
	UINT16					m_widthPages;
	std::vector<CDXCoordinate>	m_splitterPositions;
	CDXPageDefinition		m_pageDefinition;
	CDXRectangle			m_boundsInParent;
	UINT32					m_flags;
		// This is a set of bits set for attributes which have been given values.
		// Note that some attributes have defaults; they are initialized with
		// the defaults and always have values.
		// These are the ones which don't have defaults:
		static const UINT32 has_boundingBox;
		static const UINT32 has_pageWidth;
		static const UINT32 has_pageHeight;
		static const UINT32 has_headerPosition;
		static const UINT32 has_footerPosition;
		static const UINT32 has_pageOverlap;
		static const UINT32 has_printTrimMarks;
		static const UINT32 has_heightPages;
		static const UINT32 has_widthPages;
		static const UINT32 has_bgColor;
		static const UINT32 has_boundsInParent;

	// default page width/height in points
	static const UINT16 DefaultPageWidth = 540;
	static const UINT16 DefaultPageHeight = 720;

private:
	CDXPage& operator=(const CDXPage &a); /* assignment is not implemented */


public:
	CDXPage(CDXObjectID id);
	CDXPage(const CDXPage &a);
	virtual	~CDXPage();
	virtual CDXObject*	Clone()	const;

	bool			Known(UINT32 f)	const	{ return (m_flags & f) != 0; }
	void			Known(UINT32 f, bool val) { if (val) m_flags |= f;  else m_flags &= ~f; }
	bool			KnownBGColor()	const	{ return (m_flags & has_bgColor) != 0; }

	CDXRectangle	BoundingBox()	const	{ return m_boundingBox; }
	UINT32			BGColor()		const	{ assert (KnownBGColor());		return m_bgColor; }
	CDXRectangle	BoundsInParent()const	{ return m_boundsInParent; }

	void			BoundingBox		(const CDXRectangle &r)	{ m_flags |= has_boundingBox;		m_boundingBox = r; }
	void			BGColor			(UINT32       c)		{ m_flags |= has_bgColor;			m_bgColor     = c; }
	void			BoundsInParent	(const CDXRectangle &r)	{ m_flags |= has_boundsInParent;	m_boundsInParent = r; }
};

// ******************************
// ** class CDXGroupOrFragment **
// ******************************
//
// There are some operations for which Fragments or Groups of Fragments should be permitted,
// so this abstract base class is used to encompass them both.

class CORE_CHEMISTRY_API CDXGroupOrFragment : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag     () const = 0;
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &sink_arg) const;

private:
	CDXGroupOrFragment& operator=(const CDXGroupOrFragment &a); /* assignment is not implemented */

public:
	std::vector<CDXObjectID>		m_connectionOrdering;
	CDXSeqType						m_sequenceType;
	bool							m_integral;

	CDXGroupOrFragment(CDXTag tag, CDXObjectID id);
	CDXGroupOrFragment(const CDXGroupOrFragment &a);
	virtual	~CDXGroupOrFragment() = 0;
	virtual void FinishReading() { ConnectBonds(this); }

	CDXNode *		GetAttachmentNode(size_t attachmentIndex);
	const CDXNode *	GetAttachmentNode(size_t attachmentIndex) const	{ return const_cast<CDXGroupOrFragment*>(this)->GetAttachmentNode (attachmentIndex); }
	void AddAttachmentNode (const CDXNode *pNode);
};

// ********************
// ** class CDXGroup **
// ********************
//
// Specialization of CDXObject for CDXGroup objects

class CORE_CHEMISTRY_API CDXGroup : public CDXGroupOrFragment
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Group; }
	virtual std::string XMLObjectName() const;

private:
	CDXGroup& operator=(const CDXGroup &a); /* assignment is not implemented */

public:
	CDXGroup(CDXObjectID id);
	CDXGroup(const CDXGroup &a);
	CDXGroup(const CDXFragment &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXGroup();
};

// ***********************
// ** class CDXFragment **
// ***********************
//
// Specialization of CDXObject for CDXFragment objects

class CORE_CHEMISTRY_API CDXFragment : public CDXGroupOrFragment
{
protected:

	CDXTag	GetTag() const override { return kCDXObj_Fragment; }
	std::string XMLObjectName() const  override;
    void StoreAttribute(CDXDataSource& src, CDXTag tag, size_t size)  override;
    void XMLStoreAttribute(CDXTag, const std::string&) override;
    void WriteAttributesTo(CDXDataSink& sink) const override;
    void XMLWriteAttributes(XMLDataSink&) const override;
    
private:

	CDXFragment& operator=(const CDXFragment &a); /* assignment is not implemented */

public:

	// constructors (normal, copy, virtual) and destructor
	CDXFragment(CDXObjectID id);
	CDXFragment(const CDXFragment &a);
	CDXObject*	Clone()	const override;
	~CDXFragment() override;

	void CountElements(std::map<INT16, INT16, std::less<INT16> > &elementList);

public:

    bool isFromGuidedStereo = false;
    bool isComplement = false;
};

// *******************
// ** class CDXNode **
// *******************
//
// Specialization of CDXObject for CDXNode objects

#define kNumHydrogenUnspecified (UINT16)-1
#define kNumHydrogenHidden		(UINT16) 0

class CORE_CHEMISTRY_API CDXNode : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Node; }
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXNode& operator=(const CDXNode &a); // assignment is not implemented

public:
	// constructors (normal, copy, virtual) and destructor
	CDXNode(CDXObjectID id);
	CDXNode(const CDXNode &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXNode();

	INT32				m_charge;
	CDXNodeType			m_nodeType;
	UINT16				m_elementNum;
	UINT16				m_numHydrogens;		// #implicit H. Default=kNumHydrogenUnspecified. kNumHydrogenHidden=none present in label.
	INT16				m_isotope;
	CDXRadical			m_radical;			// Formerly (< 3/23/99) was INT8
	CDXTag				m_substituentCode;
	INT8				m_substituentCount;	// only defined if m_substituentCode defined
	CDXRingBondCount	m_ringBondCount;
	CDXUnsaturation		m_unsaturation;		// Formerly (< 3/23/99) was INT8
	CDXReactionStereo	m_rxnStereo;
	CDXTranslation		m_translation;
	INT8				m_implicitHydrogens;// kCDXProp_Atom_RestrictImplicitHydrogens
	INT8				m_abnormalValence;	// a boolean value
	INT8				m_restrictRxnChange;
	CDXAbundance		m_isotopicAbundance;
	CDXExternalConnectionType	m_externalConnectionType;
    INT8                m_externalConnectionNum;
	CDXAtomGeometry		m_geometry;
	CDXAtomCIPType		m_CIP;
	CDXObjectID			m_altGroupID;
	CDXLabelDisplay		m_labelDisplay;
	CDXTag				m_hStereo;			// zero, kCDXProp_Atom_HDot, or kCDXProp_Atom_HDash
	std::string			m_genericNickname;	// kCDXProp_Atom_GenericNickname
	std::string			m_atNum;
	std::string			m_residueID;
   	std::vector<CDXBond *>		m_alphaBonds;
	std::vector<CDXNode *>		m_alphaAtoms;
	std::vector<CDXObjectID> *	m_bondOrdering;
	std::vector<CDXObjectID> *	m_attachments;	// IDs of nodes to which variable or multiple attachment node is attached
    std::vector<CDXObjectID>    m_hydrogenBonds;
    std::vector<CDXObjectID>    m_hydrogenBondAttachmentAtoms;
    std::vector<INT16> *		m_elementList;	// atomic numbers of permitted or disallowed elements
	std::vector<std::string> *	m_genericList;	// names of permitted or disallowed generic nicknames
	bool				m_negativeList;			// true if m_elementList and m_genericlist are disallowed items, false if permitted items
	INT16				m_linkCountLow;
	INT16				m_linkCountHigh;
	CDXEnhancedStereoType m_enhancedStereoType;
	INT16				m_enhancedStereoGroupNum;
	bool				m_needsClean;
    bool			    showAtomID;
    UINT32              atomID;

	const CDXNode* GetAttachmentNodeForBond(CDXObjectID bondID) const;

	bool KnownAtNum()			const	{ return (m_flags & has_showAtomNumber) != 0; }
	bool KnownResidueID()		const	{ return (m_flags2 & has_showResidueID) != 0; }
	void ShowAtNum	(bool b)			{ m_flags |= has_showAtomNumber;  m_showAtomNumber = b; }
	void ShowResidueID	(bool b)		{ m_flags2 |= has_showResidueID;  m_showResidueID = b; }
	bool GetAtNumShown()		const	{ return m_showAtomNumber; }
	bool GetResidueIDShown()	const	{ return m_showResidueID; }
};

// *******************
// ** class CDXBond **
// *******************
//
// Specialization of CDXObject for CDXBond objects

class CORE_CHEMISTRY_API CDXBond : public CDXGraphicObject
{
public:
    
    // constructors (normal, copy, virtual) and destructor
    CDXBond(CDXObjectID id);
    CDXBond(const CDXBond &a);
    virtual CDXObject*    Clone()    const;
    virtual    ~CDXBond();
    
    void        Connects(CDXNode *, CDXNode *);
    CDXNode*    OtherNode(CDXNode* thisNode);
    
    void        SetIsRouted(bool toRoute) { isRouted = toRoute; }
    bool        IsRouted() const { return isRouted; }
    bool        MustBothAtomsResideInSameFragment() const;
    
protected:
    
	virtual CDXTag	GetTag() const	{ return kCDXObj_Bond; }
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	void		 WriteAttributesTo(CDXDataSink &sink) const;
	void		 FinishReading();
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual void XMLStoreAttribute(CDXTag attribTag_arg, const std::string &value_arg);

private:
    
	CDXBond& operator=(const CDXBond &a); /* assignment is not implemented */
    
public:

	CDXObjectID 	m_beginNodeID;
	CDXObjectID 	m_endNodeID;
	CDXBondOrder	m_bondOrder;	// actually, a bitmask of CDXBondOrder's
	CDXBondDisplay	m_display;
	CDXBondDisplay	m_display2;
	CDXConnectivity m_connectivity;
	CDXBondDoublePosition m_doublePosition;
	CDXBondTopology	m_topology;
	CDXBondReactionParticipation m_rxnParticipation;
	CDXBondCIPType	m_CIP;
	CDXNode *		m_beginNode;
	CDXNode *		m_endNode;
	std::vector<CDXObjectID> *	m_bondOrdering;
	std::vector<CDXObjectID> *	m_crossingBonds;
	INT16			m_beginAttach;
	INT16			m_endAttach;
    INT8            m_beginExternalNum;
    INT8            m_endExternalNum;

private:
    
    bool            isRouted = false;
};

// *******************
// ** class CDXText **
// *******************
//
// Specialization of CDXObject for CDXText objects

class CORE_CHEMISTRY_API CDXText : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Text; }
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual bool XMLNeedToWriteContent() const;
	virtual void XMLWriteContent(XMLDataSink &) const;

private:
	CDXText& operator=(const CDXText &a); // assignment is not implemented

public:
	// constructors (normal, copy, virtual) and destructor
	CDXText(CDXObjectID id, const std::string &text = "");
	CDXText(const CDXText &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXText();

	std::vector<UINT16>	*	m_lineStarts;
	UINT32					m_rotationAngle;
	CDXCoordinate			m_lineHeight;
	CDXCoordinate			m_wrapWidth;
	CDXTextJustification	m_justification;
	CDXLabelDisplay			m_labelAlignment;
    
    const CDXString&        GetText() const;
	void                    SetText(const CDXString& inString);
	virtual void            AppendText(const CDXString& inString);
    
    virtual void FinalizeTextAfterRead(const CDXFontTable& inFontTable);
    virtual void PrepareTextForWrite(const CDXFontTable& inFontTable);

private:
    
    // We have two representations of our text - a legacy one for backwards-compatibility
    // and a UTF-8 version
    CDXString				m_legacyText;
    CDXString               m_utf8Text;
};

// ***********************
// ** class CDXSequence **
// ***********************
//
// Specialization of CDXObject for CDXSequence objects

class CORE_CHEMISTRY_API CDXSequence : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Sequence; }
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXSequence& operator=(const CDXSequence &a); // assignment is not implemented

public:
	// constructors (normal, copy, virtual) and destructor
	CDXSequence(CDXObjectID id);
	CDXSequence(const CDXSequence &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXSequence() {}

	std::string				m_identifier;
};

// *****************************
// ** class CDXCrossReference **
// *****************************
//
// Specialization of CDXObject for CDXCrossReference objects

class CORE_CHEMISTRY_API CDXCrossReference : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_CrossReference; }
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXCrossReference& operator=(const CDXCrossReference &a); // assignment is not implemented

public:
	// constructors (normal, copy, virtual) and destructor
	CDXCrossReference(CDXObjectID id);
	CDXCrossReference(const CDXCrossReference &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXCrossReference() {}

	std::string				m_container;
	std::string				m_document;
	std::string				m_identifier;
	std::string				m_sequence;
};

// **********************
// ** class CDXGraphic **
// **********************
//
// Specialization of CDXObject for CDXGraphic objects

class CDXBracketedGroup;

class CORE_CHEMISTRY_API CDXGraphic : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Graphic; }
	virtual std::string XMLObjectName() const;
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXGraphic& operator=(const CDXGraphic &a); // assignment is not implemented

public:
	CDXGraphic(CDXObjectID id);
	CDXGraphic(const CDXGraphic &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXGraphic();

	std::vector<CDXObjectID> GetOtherBracketsInGroup() const;
	const CDXBracketedGroup* GetBracketedGroupForBracket() const;

	CDXGraphicType			m_graphicType;
	CDXArrowType			m_arrowType;
	CDXBracketType			m_bracketType;
	CDXRectangleType		m_rectangleType;
	CDXOvalType				m_ovalType;
	CDXSymbolType			m_symbolType;
	CDXOrbitalType			m_orbitalType;
	CDXFrameType			m_frameType;
	INT16					m_headSize;
	INT16					m_angularSize;
	INT16					m_bracketLipSize;
	CDXBracketUsage			m_bracketUsage;
	CDXPolymerRepeatPattern	m_polymerRepeatPattern;
	CDXPolymerFlipType		m_polymerFlipType;
	CDXPoint3D				m_3dHead;
	bool					m_has3dHead;
	CDXPoint3D				m_3dTail;
	bool					m_has3dTail;
	INT16					m_shadowSize;
	CDXPoint3D				m_center;
	bool					m_hasCenter;
	CDXPoint3D				m_majorAxisEnd;
	bool					m_hasMajorAxisEnd;
	CDXPoint3D				m_minorAxisEnd;
	bool					m_hasMinorAxisEnd;
	INT16					m_cornerRadius;
	INT16					m_fadePercent;
};

// ********************
// ** class CDXArrow **
// ********************
//
// Specialization of CDXObject for CDXArrow objects

class CORE_CHEMISTRY_API CDXArrow : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Arrow; }
	virtual std::string XMLObjectName() const;
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXArrow& operator=(const CDXArrow &a); // assignment is not implemented

public:
	CDXArrow(CDXObjectID id);
	CDXArrow(const CDXArrow &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXArrow();

	CDXArrowHeadType		m_arrowHeadType;
	CDXArrowHeadPosition	m_arrowHeadHead;
	CDXArrowHeadPosition	m_arrowHeadTail;
	INT16					m_headSize;
	INT16					m_headCenterSize;
	INT16					m_headWidth;
	INT16					m_angularSize;
	INT16					m_shaftSpacing;
	INT16					m_equilibriumRatio;
	INT16					m_fadePercent;
	CDXNoGoType				m_nogo;
	bool					m_dipole;
	CDXPoint3D				m_3dHead;
	bool					m_has3dHead;
	CDXPoint3D				m_3dTail;
	bool					m_has3dTail;
	CDXPoint3D				m_center;
	bool					m_hasCenter;
	CDXPoint3D				m_majorAxisEnd;
	bool					m_hasMajorAxisEnd;
	CDXPoint3D				m_minorAxisEnd;
	bool					m_hasMinorAxisEnd;
	INT16					m_sourceID;
	INT16					m_targetID;
};


// ********************
// ** class CDXCurve **
// ********************
//
// Specialization of CDXObject for CDXCurve objects

class CORE_CHEMISTRY_API CDXCurve : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Curve; }
	virtual std::string XMLObjectName() const;
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXCurve& operator=(const CDXCurve &a); // assignment is not implemented

public:
	CDXCurve(CDXObjectID id);
	CDXCurve(const CDXCurve &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXCurve();

	INT16 m_curveType;
	CDXArrowHeadType		m_arrowHeadType;
	CDXArrowHeadPosition	m_arrowHeadAtEnd;
	CDXArrowHeadPosition	m_arrowHeadAtStart;
	INT16					m_headSize;
	INT16					m_headCenterSize;
	INT16					m_headWidth;
	INT16					m_spacing;
	INT16					m_fadePercent;
	bool					m_closed;
	std::vector<CDXPoint2D> m_curvePoints;
	std::vector<CDXPoint3D> m_curvePoints3D;
};

// *****************************
// ** class CDXEmbeddedObject **
// *****************************
//
// Specialization of CDXObject for CDXEmbeddedObject objects

class CORE_CHEMISTRY_API CDXEmbeddedObject : public CDXGraphicObject
{
public:

    CDXEmbeddedObject(CDXObjectID id);
    CDXEmbeddedObject(const CDXEmbeddedObject& a);
    CDXObject* Clone() const override;

    bool HasFormatData(const std::string& mimeType) const;
    std::string GetFormatData(const std::string& mimeType) const;
    uint32_t GetFormatUncompressedSize(const std::string& mimeType) const;
    void SetFormatData(const std::string& mimeType, const std::string& data);
    void SetFormatUncompressedSize(const std::string& mimeType, const uint32_t& size);
    std::vector<std::string> GetAvailableFormats() const;

public:

    UINT32 m_rotationAngle;

protected:
	
    CDXTag GetTag() const override { return kCDXObj_EmbeddedObject; }
	std::string XMLObjectName() const override;
	void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz) override;
	void XMLStoreAttribute(CDXTag, const std::string &) override;
	void WriteAttributesTo(CDXDataSink &) const override;
	void XMLWriteAttributes(XMLDataSink &) const override;

private:
	
    // Assignment is not implemented
    CDXEmbeddedObject& operator=(const CDXEmbeddedObject &a) = delete;

    void PutTagData(CDXDataSink& sink, const std::string& dataMimeType) const;

private:

    struct FormatInfo
    {
        std::string data;
        uint32_t uncompressedSize = 0;
    };

    std::map<std::string, FormatInfo> mimeTypeToFormatInfoMap;
};

// ************************************
// ** class CDXNamedAlternativeGroup **
// ************************************
//
// Specialization of CDXObject for CDXNamedAlternativeGroup objects

class CORE_CHEMISTRY_API CDXNamedAlternativeGroup : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_NamedAlternativeGroup; }
	virtual std::string XMLObjectName() const;
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t siz);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXNamedAlternativeGroup& operator=(const CDXNamedAlternativeGroup &a); /* assignment is not implemented */

	CDXRectangle m_textFrame;
	CDXRectangle m_groupFrame;
	INT16 m_valence;
	UINT16 m_altGroupFlags;

	static const UINT16 has_textFrame;
	static const UINT16 has_groupFrame;

public:
	CDXNamedAlternativeGroup(CDXObjectID id);
	CDXNamedAlternativeGroup(const CDXNamedAlternativeGroup &a);
	virtual	~CDXNamedAlternativeGroup();
	virtual CDXObject*	Clone()	const;

	bool KnownTextFrame()					const	{ return (m_altGroupFlags & has_textFrame) != 0; }
	const CDXRectangle&	TextFrame()			const	{ assert (KnownTextFrame()); return m_textFrame; }
	void TextFrame(const CDXRectangle &r)			{ m_altGroupFlags |= has_textFrame; m_textFrame = r; }

	bool KnownGroupFrame()					const	{ return (m_altGroupFlags & has_groupFrame) != 0; }
	const CDXRectangle&	GroupFrame()		const	{ assert (KnownGroupFrame()); return m_groupFrame; }
	void GroupFrame(const CDXRectangle &r)			{ m_altGroupFlags |= has_groupFrame; m_groupFrame = r; }

	INT16 Valence() const { return m_valence; }
	void Valence(INT16 v) { m_valence = v; }
};

// ***************************
// ** class CDXTemplateGrid **
// ***************************
//
// Specialization of CDXObject for CDXTemplateGrid objects

class CORE_CHEMISTRY_API CDXTemplateGrid : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_TemplateGrid; }
	virtual std::string XMLObjectName() const;
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXTemplateGrid& operator=(const CDXTemplateGrid &a); /* assignment is not implemented */

public:
	CDXCoordinate m_paneHeight;
	INT16 m_numRows;
	INT16 m_numColumns;
	CDXPoint2D m_extent;
	
	CDXTemplateGrid(CDXObjectID id);
	CDXTemplateGrid(const CDXTemplateGrid &a);
	virtual	~CDXTemplateGrid();
	virtual CDXObject*	Clone()	const;
};

// *****************************
// ** class CDXReactionScheme **
// *****************************
//
class CORE_CHEMISTRY_API CDXReactionScheme : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_ReactionScheme; }
//	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
//	virtual void XMLStoreAttribute(CDXTag, const std::string &);
//	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual std::string XMLObjectName() const;
//	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXReactionScheme& operator=(const CDXReactionScheme &a); /* assignment is not implemented */

public:
	CDXReactionScheme(CDXObjectID id);
	CDXReactionScheme(const CDXReactionScheme &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXReactionScheme();
};

// ***************************
// ** class CDXReactionStep **
// ***************************
//
class CORE_CHEMISTRY_API CDXReactionStep : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_ReactionStep; }
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual std::string XMLObjectName() const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

private:
	CDXReactionStep& operator=(const CDXReactionStep &a); /* assignment is not implemented */

public:
	CDXReactionStep(CDXObjectID id);
	CDXReactionStep(const CDXReactionStep &a);
	virtual CDXObject*	Clone()	const;
	virtual	~CDXReactionStep();

	typedef std::vector<std::pair<CDXObjectID,CDXObjectID> >	Map;
	std::vector<CDXObjectID>	m_reactants;
	std::vector<CDXObjectID>	m_products;
	std::vector<CDXObjectID>	m_plusses;
	std::vector<CDXObjectID>	m_arrows;
	std::vector<CDXObjectID>	m_objectsAboveArrow;
	std::vector<CDXObjectID>	m_objectsBelowArrow;
	Map							m_aamap;
	Map							m_aamapManual;
	Map							m_aamapAuto;
};

// *******************************
// ** class CDXObjectDefinition **
// *******************************
//
// Specialization of CDXObject for CDXObjectDefinition objects

class CORE_CHEMISTRY_API CDXObjectDefinition : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_ObjectDefinition; }
	virtual std::string XMLObjectName() const;

private:
	CDXObjectDefinition& operator=(const CDXObjectDefinition &a); /* assignment is not implemented */

public:
	CDXObjectDefinition(CDXObjectID id);
	CDXObjectDefinition(const CDXObjectDefinition &a);
	virtual	~CDXObjectDefinition();
	virtual CDXObject*	Clone()	const;
};

// ***********************
// ** class CDXSpectrum **	Specialization of CDXObject for CDXSpectrum objects
// ***********************
class CORE_CHEMISTRY_API CDXSpectrum  : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Spectrum; }
	virtual std::string XMLObjectName() const;

private:
	CDXSpectrum& operator=(const CDXSpectrum &a); /* assignment is not implemented */

	void XMLProcessTokens(const std::string &);

public:
				CDXSpectrum(CDXObjectID id);
				CDXSpectrum(const CDXSpectrum &a);
	virtual		CDXObject*	Clone()	const;
	virtual		~CDXSpectrum();
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual bool XMLNeedToWriteContent() const;
	virtual void XMLWriteContent(XMLDataSink &) const;
	virtual void XMLStoreCharacterData(const std::string &);
	virtual void FinishReading();

	double				m_xSpacing;
	double				m_xLow;
	CDXSpectrumXType	m_xType;
	CDXSpectrumYType	m_yType;
	CDXSpectrumClass	m_class;
	// Todo: m_xAxisLabel and m_yAxisLabel should probably become CDXStrings
	CDXString			m_xAxisLabel;
	CDXString			m_yAxisLabel;
	double				m_yLow;
	double				m_yScale;		// If yScale is non-zero, write out data in scaled integer form
	std::vector<double>	m_dataPoints;
};

// *******************************
// ***** class CDXObjectTag ******
// *******************************
//
// Specialization of CDXObject for CDXObjectTag objects

class CORE_CHEMISTRY_API CDXObjectTag : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_ObjectTag; }
	virtual std::string XMLObjectName() const;

private:
	CDXObjectTag& operator=(const CDXObjectTag &a); /* assignment is not implemented */

public:
	CDXObjectTag(CDXObjectID id);
	CDXObjectTag(const CDXObjectTag &a);
	
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXObjectTagType	m_Type;
	std::string			m_Name;
	bool				m_Tracking;
	bool				m_Persistent;
    int32_t             m_Int32Val;
    double              m_DoubleVal;
    std::string         m_StringVal;
	CDXPositioningType	m_positioning;
	UINT32				m_positioningAngle;	// in degrees
	CDXPoint2D			m_positioningOffset;
};

std::istream & operator>>(std::istream &is, CDXColor& c);
CORE_CHEMISTRY_API std::ostream & operator<<(std::ostream &os, const CDXColor &c);
std::ostream & operator<<(std::ostream &os, const CDXFontTableEntry &f);
std::ostream & operator<<(std::ostream &os, const CDXPropRep &p);

// *******************************
// ***** class CDXSplitter  ******
// *******************************
//
// Specialization of CDXObject for CDXSplitter objects

class CORE_CHEMISTRY_API CDXSplitter : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Splitter; }
	virtual std::string XMLObjectName() const;

private:
	CDXSplitter& operator=(const CDXSplitter &a); /* assignment is not implemented */

public:
	CDXSplitter(CDXObjectID id);
	CDXSplitter(const CDXSplitter &a);
	virtual	~CDXSplitter() {}
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXPoint2D			m_pos;
	CDXPageDefinition	m_def;
};

// *******************************
// *****   class CDXTable   ******
// *******************************
//
// Specialization of CDXObject for CDXTable objects

class CORE_CHEMISTRY_API CDXTable : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Table; }
	virtual std::string XMLObjectName() const;

private:
	CDXTable& operator=(const CDXTable &t); /* assignment is not implemented */

public:
	CDXTable(CDXObjectID id);
	CDXTable(const CDXTable &a);
	virtual	~CDXTable() {}
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

};

// *******************************************
// *****   class CDXPlasmidMap   ******
// *******************************************
//
// Specialization of CDXObject for CDXPlasmidMap objects

class CORE_CHEMISTRY_API CDXPlasmidMap : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_PlasmidMap; }
	virtual std::string XMLObjectName() const;

private:
	CDXPlasmidMap& operator=(const CDXPlasmidMap &t); /* assignment is not implemented */

public:
	CDXPlasmidMap(CDXObjectID id);
	CDXPlasmidMap(const CDXPlasmidMap &a);
	virtual	~CDXPlasmidMap() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	INT32			m_numberBasePairs;
	CDXCoordinate	m_ringRadius;
};



// *******************************************
// *****   class CDXPlasmidMarker   ******
// *******************************************
//
// Specialization of CDXObject for CDXPlasmidMarker objects

class CORE_CHEMISTRY_API CDXPlasmidMarker : public CDXObjectTag
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_PlasmidMarker; }
	virtual std::string XMLObjectName() const;

private:
	CDXPlasmidMarker& operator=(const CDXPlasmidMarker &t); /* assignment is not implemented */

public:
	CDXPlasmidMarker(CDXObjectID id);
	CDXPlasmidMarker(const CDXPlasmidMarker &a);
	virtual	~CDXPlasmidMarker() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXCoordinate	m_markerOffset;
	CDXCoordinate	m_markerAngle;

};

// *******************************************
// *****   class CDXPlasmidRegion   ******
// *******************************************
//
// Specialization of CDXObject for CDXPlasmidRegion objects

class CORE_CHEMISTRY_API CDXPlasmidRegion : public CDXArrow
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_PlasmidRegion; }
	virtual std::string XMLObjectName() const;

private:
	CDXPlasmidRegion& operator=(const CDXPlasmidRegion &t); /* assignment is not implemented */

public:
	CDXPlasmidRegion(CDXObjectID id);
	CDXPlasmidRegion(const CDXPlasmidRegion &a);
	virtual	~CDXPlasmidRegion() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	INT32			m_regionStart;
	INT32			m_regionEnd;
	CDXCoordinate	m_regionOffset;


};

// *******************************************
// *****   class CDXBandMarker   ******
// *******************************************
//
// Specialization of CDXObject for CDXPlasmidMarker objects

class CORE_CHEMISTRY_API CDXBandMarker : public CDXObjectTag
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Marker; }
	virtual std::string XMLObjectName() const;

private:
	CDXBandMarker& operator=(const CDXBandMarker &t); /* assignment is not implemented */

public:
	CDXBandMarker(CDXObjectID id);
	CDXBandMarker(const CDXBandMarker &a);
	virtual	~CDXBandMarker() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
};

// *******************************************
// *****   class CDXRLogic   ******
// *******************************************
//
// Specialization of CDXObject for CDXRLogic objects

class CORE_CHEMISTRY_API CDXRLogic : public CDXText
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_RLogic; }
	virtual std::string XMLObjectName() const;

private:
	CDXRLogic& operator=(const CDXRLogic &t); /* assignment is not implemented */

public:
	CDXRLogic(CDXObjectID id);
	CDXRLogic(const CDXRLogic &a);
	virtual	~CDXRLogic() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
};


// *******************************************
// *****   class CDXRLogicItem   ******
// *******************************************
//
// Specialization of CDXObject for CDXRLogicItem objects

class CORE_CHEMISTRY_API CDXRLogicItem : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_RLogicItem; }
	virtual std::string XMLObjectName() const;

private:
	CDXRLogicItem& operator=(const CDXRLogicItem &t); /* assignment is not implemented */

public:
	CDXRLogicItem(CDXObjectID id);
	CDXRLogicItem(const CDXRLogicItem &a);
	virtual	~CDXRLogicItem() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	string		m_occurrence;
	string		m_group;
	string		m_itgroup;
	bool		m_restH;
};


// *******************************************
// *****   class CDXStoichiometryGrid   ******
// *******************************************
//
// Specialization of CDXObject for CDXStoichiometryGrid objects

class CORE_CHEMISTRY_API CDXStoichiometryGrid : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_StoichiometryGrid; }
	virtual std::string XMLObjectName() const;

private:
	CDXStoichiometryGrid& operator=(const CDXStoichiometryGrid &t); /* assignment is not implemented */

public:
	CDXStoichiometryGrid(CDXObjectID id);
	CDXStoichiometryGrid(const CDXStoichiometryGrid &a);
	virtual	~CDXStoichiometryGrid() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
};

// ************************************
// *****   class CDXSGComponent   ******
// ************************************
//
// Specialization of CDXObject for CDXSGComponent objects

class CORE_CHEMISTRY_API CDXSGComponent : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_SGComponent; }
	virtual std::string XMLObjectName() const;

private:
	CDXSGComponent& operator=(const CDXSGComponent &t); /* assignment is not implemented */

public:
	CDXSGComponent(CDXObjectID id);
	CDXSGComponent(const CDXSGComponent &a);
	virtual	~CDXSGComponent() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXCoordinate		m_width;
	bool				m_isReactant;
	bool				m_isHeader;
	INT16				m_referenceID;
};


// *********************************
// *****   class CDXSGDatum   ******
// *********************************
//
// Specialization of CDXObject for CDXSGDatum objects

class CORE_CHEMISTRY_API CDXSGDatum : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_SGDatum; }
	virtual std::string XMLObjectName() const;

private:
	CDXSGDatum& operator=(const CDXSGDatum &t); /* assignment is not implemented */

public:
	CDXSGDatum(CDXObjectID id);
	CDXSGDatum(const CDXSGDatum &a);
	virtual	~CDXSGDatum() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXSGPropertyType	m_type;
	CDXSGDataType		m_datatype;
	string				m_propertystring;
	double				m_propertyvalue;
	bool				m_isedited;
	bool				m_isreadonly;
	bool				m_ishidden;
};

// *******************************
// *****  class CDXPlateBase ******
// *******************************
//
// Specialization of CDXObject for CDXPlateBase objects

class CORE_CHEMISTRY_API CDXPlateBase : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const = 0;
	virtual std::string XMLObjectName() const = 0;

private:
	CDXPlateBase& operator=(const CDXPlateBase &t); /* assignment is not implemented */

public:
	CDXPlateBase(CDXObjectID id);
	CDXPlateBase(const CDXPlateBase &a);
	virtual	~CDXPlateBase() {};
	virtual CDXObject*	Clone()	const = 0;
	virtual CDXDatumID GetItemTypeID() const = 0;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXPoint2D		m_topleft;
	CDXPoint2D		m_topright;
	CDXPoint2D		m_bottomright;
	CDXPoint2D		m_bottomleft;
	bool			m_showBorders;
	bool			m_transparent;
	bool			m_hasTransparent;
};

class CORE_CHEMISTRY_API CDXTLCPlate : public CDXPlateBase
{
public:
	CDXTLCPlate(CDXObjectID id);
	CDXTLCPlate(const CDXTLCPlate &a);
	CDXTag	GetTag() const { return kCDXObj_TLCPlate; };

	void XMLWriteAttributes(XMLDataSink &sink_arg) const;
	void XMLStoreAttribute(CDXTag, const std::string &);
	void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	void WriteAttributesTo(CDXDataSink &sink_arg) const;

	virtual std::string XMLObjectName() const;
	CDXObject*	Clone()	const;
	CDXDatumID GetItemTypeID() const { return kCDXObj_TLCSpot; }
	
	double			m_originFraction;		// distance from the bottom to the origin, as a fraction of the total height
	double			m_solventFrontFraction;	// distance from the solvent front to the top, as a fraction of the total height
	bool			m_showOrigin;
	bool			m_showSolventFront;
	bool			m_showSideTicks;
};

class CORE_CHEMISTRY_API CDXGEPPlate : public CDXPlateBase
{
public:	
	CDXGEPPlate(CDXObjectID id);
	CDXGEPPlate(const CDXGEPPlate &a);
	CDXTag	GetTag() const { return kCDXObj_GEPPlate; };

	void XMLWriteAttributes(XMLDataSink &sink_arg) const;
	void XMLStoreAttribute(CDXTag, const std::string &);
	void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	void WriteAttributesTo(CDXDataSink &sink_arg) const;

	virtual std::string XMLObjectName() const;
	CDXObject*	Clone()	const;
	CDXDatumID GetItemTypeID() const { return kCDXObj_GEPBand; }

	INT32	m_scaleUnitID;
	bool	m_showScale;
	INT32	m_minRange;
	INT32	m_maxRange;
	INT32	m_labelsAngle;
	double	m_axisWidth;
};

// *******************************
// *****  class CDXPlateLane  ******
// *******************************
//
// Specialization of CDXObject for CDXPlateLane objects

class CORE_CHEMISTRY_API CDXPlateLane : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const = 0;
	virtual std::string XMLObjectName() const = 0;

private:
	CDXPlateLane& operator=(const CDXPlateLane &t); /* assignment is not implemented */

public:
	CDXPlateLane(CDXTag tag, CDXObjectID id);
	CDXPlateLane(const CDXPlateLane &a);
	virtual	~CDXPlateLane() {};
	virtual CDXObject*	Clone()	const = 0;
};

class CORE_CHEMISTRY_API CDXTLCLane : public CDXPlateLane
{
protected:
	CDXTag	GetTag() const	{ return kCDXObj_TLCLane; }
	std::string XMLObjectName() const;

public:
	CDXTLCLane(CDXObjectID id);//(CDXTag tag, CDXObjectID id)
	CDXTLCLane(const CDXTLCLane &a);
	CDXObject*	Clone()	const;
};

class CORE_CHEMISTRY_API CDXGEPLane : public CDXPlateLane
{
protected:
	CDXTag	GetTag() const	{ return kCDXObj_GEPLane; }
	std::string XMLObjectName() const;

public:
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);

	CDXGEPLane(CDXObjectID id);//(CDXTag tag, CDXObjectID id)
	CDXGEPLane(const CDXGEPLane &a);
	CDXObject*	Clone()	const;
};

// *******************************
// *****  class CDXPlateItemBase  ******
// *******************************
//
// Specialization of CDXObject for CDXPlateItemBase objects

class CORE_CHEMISTRY_API CDXPlateItemBase : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const = 0;
	virtual std::string XMLObjectName() const = 0;

private:
	CDXPlateItemBase& operator=(const CDXPlateItemBase &t); /* assignment is not implemented */

public:
	CDXPlateItemBase(CDXObjectID id);
	CDXPlateItemBase(const CDXPlateItemBase &a);
	virtual	~CDXPlateItemBase() {};
	virtual CDXObject*	Clone()	const = 0;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	double				m_value;
	CDXCoordinate		m_width;
	CDXCoordinate		m_height;
	INT16				m_displayType;
	bool				m_showValue;
};

class CORE_CHEMISTRY_API CDXGEPBand: public CDXPlateItemBase
{
public:
	CDXGEPBand(CDXObjectID id);
	CDXGEPBand(const CDXGEPBand &a);

	void XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg);
	void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	void XMLWriteAttributes(XMLDataSink &sink_arg) const;
	void WriteAttributesTo(CDXDataSink &sink_arg) const;

	CDXTag	GetTag() const	{ return kCDXObj_GEPBand; }
	std::string XMLObjectName() const;
	CDXObject*	Clone() const;

	INT32	m_mass;
};

class CORE_CHEMISTRY_API CDXTLCSpot: public CDXPlateItemBase
{
public:
	CDXTLCSpot(CDXObjectID id);
	CDXTLCSpot(const CDXTLCSpot &a);

	void XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg);
	void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	void XMLWriteAttributes(XMLDataSink &sink_arg) const;
	void WriteAttributesTo(CDXDataSink &sink_arg) const;

	CDXTag	GetTag() const	{ return kCDXObj_TLCSpot; }
	std::string XMLObjectName() const;
	CDXObject*	Clone() const;
	
	CDXCoordinate		m_tail;
};

// ****************************************
// *****   class CDXBracketedGroup   ******
// ****************************************
//
// Specialization of CDXObject for CDXBracketedGroup objects

class CORE_CHEMISTRY_API CDXBracketedGroup : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_BracketedGroup; }
	virtual std::string XMLObjectName() const;

private:
	CDXBracketedGroup& operator=(const CDXBracketedGroup &t); /* assignment is not implemented */

public:
	CDXBracketedGroup(CDXObjectID id);
	CDXBracketedGroup(const CDXBracketedGroup &a);
	virtual	~CDXBracketedGroup();
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	std::vector<CDXObjectID>	m_bracketedObjects;
	CDXBracketUsage				m_usage;
	CDXPolymerRepeatPattern		m_repeatPattern;
	CDXPolymerFlipType			m_flipType;
	double						m_repeatCount;
	int							m_componentOrder;
	std::string					m_SRULabel;
};

// *******************************************
// *****   class CDXBracketAttachment   ******
// *******************************************
//
// Specialization of CDXObject for CDXBracketAttachment objects

class CORE_CHEMISTRY_API CDXBracketAttachment : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_BracketAttachment; }
	virtual std::string XMLObjectName() const;

private:
	CDXBracketAttachment& operator=(const CDXBracketAttachment &t); /* assignment is not implemented */

public:
	CDXBracketAttachment(CDXObjectID id);
	CDXBracketAttachment(const CDXBracketAttachment &a);
	virtual	~CDXBracketAttachment();
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXObjectID		m_graphicID;
};

// **************************************
// *****   class CDXCrossingBond   ******
// **************************************
//
// Specialization of CDXObject for CDXCrossingBond objects

class CORE_CHEMISTRY_API CDXCrossingBond : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_CrossingBond; }
	virtual std::string XMLObjectName() const;

private:
	CDXCrossingBond& operator=(const CDXCrossingBond &t); /* assignment is not implemented */

public:
	CDXCrossingBond(CDXObjectID id);
	CDXCrossingBond(const CDXCrossingBond &a);
	virtual	~CDXCrossingBond();
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXObjectID		m_bondID;
	CDXObjectID		m_innerAtomID;
};


// *******************************
// *****   class CDXBorder   *****
// *******************************
//
// An object to hold the border of a page

class CORE_CHEMISTRY_API CDXBorder : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Border; }
	virtual std::string XMLObjectName() const;

private:
	CDXBorder& operator=(const CDXBorder &t); /* assignment is not implemented */

public:
	CDXBorder(CDXObjectID id);
	CDXBorder(const CDXBorder &a);
	virtual	~CDXBorder() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	bool			has_lineWidth;
	bool			has_color;
	bool			has_alpha;

	CDXLineType		m_lineType;
	CDXCoordinate	m_lineWidth;
	int				m_color;
	UINT16			m_alpha;
	CDXSideType		m_side;
};


// *********************************
// *****   class CDXGeometry   *****
// *********************************
//
// An object to hold a geometrical relation between other objects

class CORE_CHEMISTRY_API CDXGeometry : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Geometry; }
	virtual std::string XMLObjectName() const;

private:
	CDXGeometry& operator=(const CDXGeometry &t); /* assignment is not implemented */

public:
	CDXGeometry(CDXObjectID id);
	CDXGeometry(const CDXGeometry &a);
	virtual	~CDXGeometry() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXGeometricFeature			m_geometricFeature;
	double						m_relationValue;
	bool						m_pointIsDirected;
	std::vector<CDXObjectID>	m_basisObjects;
};


// ***********************************
// *****   class CDXConstraint   *****
// ***********************************
//
// An object to hold a distance or angle constraint between other objects

class CORE_CHEMISTRY_API CDXConstraint : public CDXGraphicObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Constraint; }
	virtual std::string XMLObjectName() const;

private:
	CDXConstraint& operator=(const CDXConstraint &t); /* assignment is not implemented */

public:
	CDXConstraint(CDXObjectID id);
	CDXConstraint(const CDXConstraint &a);
	virtual	~CDXConstraint() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXConstraintType			m_constraintType;
	double						m_min;
	double						m_max;
	bool						m_ignoreUnconnectedAtoms;
	bool						m_dihedralIsChiral;
	std::vector<CDXObjectID>	m_basisObjects;
};

// *******************************
// ** class CDXChemicalProperty **
// *******************************
//
// Specialization of CDXText for CDXChemicalProperty objects

class CORE_CHEMISTRY_API CDXChemicalProperty : public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_ChemicalProperty; }
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void StoreAttribute(CDXDataSource &src, CDXTag tag, size_t size);
	virtual void WriteAttributesTo(CDXDataSink &) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual std::string XMLObjectName() const;

private:
	CDXText& operator=(const CDXChemicalProperty &a); // assignment is not implemented

public:

    typedef std::pair<CDXObjectID, CDXObjectID> ExternalBond;

	// constructors (normal, copy, virtual) and destructor
	CDXChemicalProperty(const CDXChemicalProperty &src);
	CDXChemicalProperty(CDXObjectID id);
	virtual CDXObject* Clone() const;
	virtual	~CDXChemicalProperty();
    bool ShouldWriteToMDL() const;
    bool HasAttachedDataType() const;

	std::string							m_name;
	bool								m_isActive;
    bool                                m_isChemicallySignificant;
	CDXObjectID							m_displayID;
	CDXPositioningType					m_positioning;
	UINT32								m_positioningAngle;	// in degrees
	CDXPoint2D							m_positioningOffset;

	vector< CDXObjectID >				m_basisObjects;
	vector< CDXChemicalPropertyType >	m_propertyTypes;
    std::vector<ExternalBond>           m_externalBonds;
};


// *******************************
// *****  class CDXBioShape  *****
// *******************************
class CORE_CHEMISTRY_API CDXBioShape : public CDXGraphicObject
{
public:
	CDXBioShape(CDXObjectID id_arg);
	CDXBioShape(const CDXBioShape &src);
	virtual ~CDXBioShape();
	virtual CDXTag	GetTag() const	{ return kCDXObj_BioShape; }
	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;
	virtual std::string XMLObjectName() const;
	virtual CDXObject* Clone() const;

	enum CDXBioShapeProperty {
		has_majorAxisEnd					= 0x00000001,
		has_minorAxisEnd					= 0x00000002,
		has_1SubstrateEnzyme_ReceptorSize	= 0x00000004,
		has_receptor_NeckWidth				= 0x00000008,
		has_immunoglobin_Height				= 0x00000010,
		has_immunoglobin_Width				= 0x00000020,
		has_gprotein_UpperHeight			= 0x00000040,
		has_membrane_ElementSize			= 0x00000200,
		has_membrane_MajorAxisSize			= 0x00000400,
		has_membrane_MinorAxisSize			= 0x00000800,
		has_membrane_StartAngle				= 0x00001000,
		has_membrane_EndAngle				= 0x00002000,
		has_dna_WaveHeight					= 0x00004000,
		has_dna_WaveLength					= 0x00008000,
		has_dna_WaveWidth					= 0x00010000,
		has_dna_Offset						= 0x00020000,
		has_helixProtein_CylinderWidth		= 0x00040000,
		has_helixProtein_CylinderHeight		= 0x00080000,
		has_helixProtein_CylinderDistance	= 0x00100000,
		has_helixProtein_PipeWidth			= 0x00200000,
		has_helixProtein_Extra				= 0x00400000,
		has_fadePercent                     = 0x00800000
	};

	bool Known(CDXBioShapeProperty f)				const	{ return (m_flags & f) != 0; }
	void Known(CDXBioShapeProperty f, bool val)				{ if (val) m_flags |= f;  else m_flags &= ~f; }

	CDXBioShapeType			m_type;
	CDXPoint3D				m_majorAxisEnd;
	CDXPoint3D				m_minorAxisEnd;
	INT16					m_fadePercent;
	CDXCoordinate			m_1SubstrateEnzyme_ReceptorSize;
	CDXCoordinate			m_receptor_NeckWidth;
	CDXCoordinate			m_gprotein_UpperHeight;
	CDXCoordinate			m_membrane_ElementSize;
	INT32					m_membrane_StartAngle;
	INT32					m_membrane_EndAngle;
	CDXCoordinate			m_dna_WaveHeight;
	CDXCoordinate			m_dna_WaveLength;
	CDXCoordinate			m_dna_WaveWidth;
	CDXCoordinate			m_dna_Offset;
	CDXCoordinate			m_helixProtein_CylinderWidth;
	CDXCoordinate			m_helixProtein_CylinderHeight;
	CDXCoordinate			m_helixProtein_CylinderDistance;
	CDXCoordinate			m_helixProtein_PipeWidth;
	CDXCoordinate			m_helixProtein_Extra;

private:
	UINT32					m_flags;
};

// *********************************
// ***** class CDXAnnotation  ******
// *********************************
//
// Specialization of CDXObject for CDXAnnotation objects

class CORE_CHEMISTRY_API CDXAnnotation	: public CDXObject
{
protected:
	virtual CDXTag	GetTag() const	{ return kCDXObj_Annotation; }
	virtual std::string XMLObjectName() const;

private:
	CDXAnnotation& operator=(const CDXAnnotation &a); /* assignment is not implemented */

public:
	CDXAnnotation(CDXObjectID id, const string &keyword = string(), const string &content = string());
	CDXAnnotation(const CDXAnnotation &a);
	virtual	~CDXAnnotation() {};
	virtual CDXObject*	Clone()	const;

	virtual void StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg);
	virtual void XMLStoreAttribute(CDXTag, const std::string &);
	virtual void WriteAttributesTo(CDXDataSink &sink_arg) const;
	virtual void XMLWriteAttributes(XMLDataSink &) const;

	CDXString		m_keyword;
	CDXString		m_content;
};

