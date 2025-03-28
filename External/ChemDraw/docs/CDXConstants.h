/*
	$Workfile: CDXConstants.h $
	$Revision: 82 $
	$Date: 7/02/02 6:39p $
	Copyright:	© 1986-2002 CambridgeSoft Corp., all rights reserved.
	
	Description:	Constants defined by the CDX file format Specification
*/

#ifndef _H_CDXConstants
#define _H_CDXConstants

// ------------------------------------------------------------
// An ANSI replacement for four-byte character constants.
// ------------------------------------------------------------
#ifndef QUADCONST
#if defined( _X86_ ) && !defined( _CS3D )		// All 3D projects uses the second QUADCONST definition for both MAC and Windows
// littleEndian
#define QUADCONST(a, b, c, d)  \
	 (((long) ((d) & 0xff) << 24)	\
	| ((long) ((c) & 0xff) << 16)	\
	| ((long) ((b) & 0xff) << 8)	\
	| ((long) ((a) & 0xff)))

#else // defined( _X86_ ) && !defined( _CS3D )
#define QUADCONST(a, b, c, d)		\
	 (((long) ((a) & 0xff) << 24)	\
	| ((long) ((b) & 0xff) << 16)	\
	| ((long) ((c) & 0xff) << 8)	\
	| ((long) ((d) & 0xff)))

#endif //  defined( _X86_ ) && !defined( _CS3D )
#endif // QUADCONST

typedef UINT16 CDXTag;
typedef INT32  CDXObjectID; // signed for now, due to mac compiler bug?
const CDXObjectID kCDXUndefinedId = (CDXObjectID)-1;

const int kCDX_HeaderStringLen = 8;
#define kCDX_HeaderString "VjCD0100"
#define kCDX_Signature	   QUADCONST('V','j','C','D')
#define kCDX_HeaderLength 28

#define kCDXML_HeaderString \
	"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << GetTextEOL() << \
	"<!DOCTYPE CDXML SYSTEM \"http://www.cambridgesoft.com/xml/cdxml.dtd\" >" << GetTextEOL()

const UINT16 kCDXTag_Object			= 0x8000;
const UINT16 kCDXTag_UserDefined	= 0x4000;

enum CDXDatumID {
	// General properties.
	kCDXProp_EndObject,						// 0x0000 Marks end of object.
	kCDXProp_CreationUserName,				// 0x0001 The name of the creator (program user's name) of the document. (CDXString)
	kCDXProp_CreationDate,					// 0x0002 The time of object creation. (CDXDate)
	kCDXProp_CreationProgram,				// 0x0003 The name of the program, including version and platform, that created the associated CDX object. ChemDraw 4.0 uses ChemDraw 4.0 as the value of CreationProgram. (CDXString)
	kCDXProp_ModificationUserName,			// 0x0004 The name of the last modifier (program user's name) of the document. (CDXString)
	kCDXProp_ModificationDate,				// 0x0005 Time of the last modification. (CDXDate)
	kCDXProp_ModificationProgram,			// 0x0006 The name of the program, including version and platform, of the last program to perform a modification. ChemDraw 4.0 uses ChemDraw 4.0 as the value of CreationProgram. (CDXString)
	kCDXProp_Unused1,						// 0x0007 Table of contents. (obsolete)
	kCDXProp_Name,							// 0x0008 Name of an object. (CDXString)
	kCDXProp_Comment,						// 0x0009 An arbitrary string intended to be meaningful to a user. (CDXString)
	kCDXProp_ZOrder,						// 0x000A Back-to-front ordering index in 2D drawing. (INT16)
	kCDXProp_RegistryNumber,				// 0x000B A registry or catalog number of a molecule object. (CDXString)
	kCDXProp_RegistryAuthority,				// 0x000C A string that specifies the authority which issued a registry or catalog number. Some examples of registry authorities are CAS, Beilstein, Aldrich, and Merck. (CDXString)
	kCDXProp_Unused2,						// 0x000D Indicates that this object (the reference object) is an alias to an object elsewhere in the document (the target object). The attributes and contained objects should be taken from the target object. (obsolete)
	kCDXProp_RepresentsProperty,			// 0x000E Indicates that this object represents some property in some other object. (CDXRepresentsProperty)
	kCDXProp_IgnoreWarnings,				// 0x000F Signifies whether chemical warnings should be suppressed on this object. (CDXBooleanImplied)
	kCDXProp_ChemicalWarning,				// 0x0010 A warning concerning possible chemical problems with this object. (CDXString)
	kCDXProp_Visible,						// 0x0011 The object is visible if non-zero. (CDXBoolean)

	// Fonts.
	kCDXProp_FontTable = 0x0100,			// 0x0100 A list of fonts used in the document. (CDXFontTable)

	// Coordinates.
	kCDXProp_2DPosition = 0x0200,			// 0x0200 The 2D location (in the order of vertical and horizontal locations) of an object. (CDXPoint2D)
	kCDXProp_3DPosition,					// 0x0201 The 3D location (in the order of X-, Y-, and Z-locations in right-handed coordinate system) of an object in CDX coordinate units. The precise meaning of this attribute varies depending on the type of object. (CDXPoint3D)
	kCDXProp_2DExtent,						// 0x0202 The width and height of an object in CDX coordinate units. The precise meaning of this attribute varies depending on the type of object. (CDXPoint2D)
	kCDXProp_3DExtent,						// 0x0203 The width, height, and depth of an object in CDX coordinate units (right-handed coordinate system). The precise meaning of this attribute varies depending on the type of object. (CDXPoint3D)
	kCDXProp_BoundingBox,					// 0x0204 The smallest rectangle that encloses the graphical representation of the object. (CDXRectangle)
	kCDXProp_RotationAngle,					// 0x0205 The angular orientation of an object in degrees * 65536. (INT32)
	kCDXProp_BoundsInParent,				// 0x0206 The bounds of this object in the coordinate system of its parent (used for pages within tables). (CDXRectangle)
	kCDXProp_3DHead,						// 0x0207 The 3D location (in the order of X-, Y-, and Z-locations in right-handed coordinate system) of the head of an object in CDX coordinate units. The precise meaning of this attribute varies depending on the type of object. (CDXPoint3D)
	kCDXProp_3DTail,						// 0x0208 The 3D location (in the order of X-, Y-, and Z-locations in right-handed coordinate system) of the tail of an object in CDX coordinate units. The precise meaning of this attribute varies depending on the type of object. (CDXPoint3D)
	kCDXProp_TopLeft,						// 0x0209 The location of the top-left corner of a quadrilateral object, possibly in a rotated or skewed frame. (CDXPoint2D)
	kCDXProp_TopRight,						// 0x020A The location of the top-right corner of a quadrilateral object, possibly in a rotated or skewed frame. (CDXPoint2D)
	kCDXProp_BottomRight,					// 0x020B The location of the bottom-right corner of a quadrilateral object, possibly in a rotated or skewed frame. (CDXPoint2D)
	kCDXProp_BottomLeft,					// 0x020C The location of the bottom-left corner of a quadrilateral object, possibly in a rotated or skewed frame. (CDXPoint2D)

	// Colors.
	kCDXProp_ColorTable = 0x0300,			// 0x0300 The color palette used throughout the document. (CDXColorTable)
	kCDXProp_ForegroundColor,				// 0x0301 The foreground color of an object represented as the two-based index into the object's color table. (UINT16)
	kCDXProp_BackgroundColor,				// 0x0302 The background color of an object represented as the two-based index into the object's color table. (INT16)
	
	// Atom properties.
	kCDXProp_Node_Type = 0x0400,			// 0x0400 The type of a node object. (INT16)
	kCDXProp_Node_LabelDisplay,				// 0x0401 The characteristics of node label display. (INT8)
	kCDXProp_Node_Element,					// 0x0402 The atomic number of the atom representing this node. (INT16)
	kCDXProp_Atom_ElementList,				// 0x0403 A list of atomic numbers. (CDXElementList)
	kCDXProp_Atom_Formula,					// 0x0404 The composition of a node representing a fragment whose composition is known, but whose connectivity is not. For example, C<sub>4</sub>H<sub>9</sub> represents a mixture of the 4 butyl isomers. (CDXFormula)
	kCDXProp_Atom_Isotope = 0x0420,			// 0x0420 The absolute isotopic mass of an atom (2 for deuterium, 14 for carbon-14). (INT16)
	kCDXProp_Atom_Charge,					// 0x0421 The atomic charge of an atom. (INT8)
	kCDXProp_Atom_Radical,					// 0x0422 The atomic radical attribute of an atom. (UINT8)
	kCDXProp_Atom_RestrictFreeSites,		// 0x0423 Indicates that up to the specified number of additional substituents are permitted on this atom. (UINT8)
	kCDXProp_Atom_RestrictImplicitHydrogens,// 0x0424 Signifies that implicit hydrogens are not allowed on this atom. (CDXBooleanImplied)
	kCDXProp_Atom_RestrictRingBondCount,	// 0x0425 The number of ring bonds attached to an atom. (INT8)
	kCDXProp_Atom_RestrictUnsaturatedBonds,	// 0x0426 Indicates whether unsaturation should be present or absent. (INT8)
	kCDXProp_Atom_RestrictRxnChange,		// 0x0427 If present, signifies that the reaction change of an atom must be as specified. (CDXBooleanImplied)
	kCDXProp_Atom_RestrictRxnStereo,		// 0x0428 The change of stereochemistry of an atom during a reaction. (INT8)
	kCDXProp_Atom_AbnormalValence,			// 0x0429 Signifies that an abnormal valence for an atom is permitted. (CDXBooleanImplied)
	kCDXProp_Unused3,						// 0x042A 
	kCDXProp_Atom_NumHydrogens,				// 0x042B The number of (explicit) hydrogens in a labeled atom consisting of one heavy atom and (optionally) the symbol H (e.g., CH<sub>3</sub>). (UINT16)
	kCDXProp_Unused4,						// 0x042C 
	kCDXProp_Unused5,						// 0x042D 
	kCDXProp_Atom_HDot,						// 0x042E Signifies the presence of an implicit hydrogen with stereochemistry specified equivalent to an explicit H atom with a wedged bond. (CDXBooleanImplied)
	kCDXProp_Atom_HDash,					// 0x042F Signifies the presence of an implicit hydrogen with stereochemistry specified equivalent to an explicit H atom with a hashed bond. (CDXBooleanImplied)
	kCDXProp_Atom_Geometry,					// 0x0430 The geometry of the bonds about this atom. (INT8)
	kCDXProp_Atom_BondOrdering,				// 0x0431 An ordering of the bonds to this node, used for stereocenters, fragments, and named alternative groups with more than one attachment. (CDXObjectIDArray)
	kCDXProp_Node_Attachments,				// 0x0432 For multicenter attachment nodes or variable attachment nodes, a list of IDs of the nodes which are multiply or variably attached to this node. (CDXObjectIDArrayWithCounts)
	kCDXProp_Atom_GenericNickname,			// 0x0433 The name of the generic nickname. (CDXString)
	kCDXProp_Atom_AltGroupID,				// 0x0434 The ID of the alternative group object that describes this node. (CDXObjectID)
	kCDXProp_Atom_RestrictSubstituentsUpTo,	// 0x0435 Indicates that substitution is restricted to no more than the specified value. (UINT8)
	kCDXProp_Atom_RestrictSubstituentsExactly,	// 0x0436 Indicates that exactly the specified number of substituents must be present. (UINT8)
	kCDXProp_Atom_CIPStereochemistry,		// 0x0437 The node's absolute stereochemistry according to the Cahn-Ingold-Prelog system. (INT8)
	kCDXProp_Atom_Translation,				// 0x0438 Provides for restrictions on whether a given node may match other more- or less-general nodes. (INT8)
	kCDXProp_Atom_AtomNumber,				// 0x0439 Atom number, as text. (CDXString)
	kCDXProp_Atom_ShowQuery,				// 0x043A Show the query indicator if non-zero. (CDXBoolean)
	kCDXProp_Atom_ShowStereo,				// 0x043B Show the stereochemistry indicator if non-zero. (CDXBoolean)
	kCDXProp_Atom_ShowAtomNumber,			// 0x043C Show the atom number if non-zero. (CDXBoolean)
	kCDXProp_Atom_LinkCountLow,				// 0x043D Low end of repeat count for link nodes. (INT16)
	kCDXProp_Atom_LinkCountHigh,			// 0x043E High end of repeat count for link nodes. (INT16)
	kCDXProp_Atom_IsotopicAbundance,		// 0x043F Isotopic abundance of this atom's isotope. (INT8)
	kCDXProp_Atom_ExternalConnectionType,	// 0x0440 Type of external connection, for atoms of type kCDXNodeType_ExternalConnectionPoint. (INT8)

	// Molecule properties.
	kCDXProp_Mole_Racemic = 0x0500,			// 0x0500 Indicates that the molecule is a racemic mixture. (CDXBoolean)
	kCDXProp_Mole_Absolute,					// 0x0501 Indicates that the molecule has known absolute configuration. (CDXBoolean)
	kCDXProp_Mole_Relative,					// 0x0502 Indicates that the molecule has known relative stereochemistry, but unknown absolute configuration. (CDXBoolean)
	kCDXProp_Mole_Formula,					// 0x0503 The molecular formula representation of a molecule object. (CDXFormula)
	kCDXProp_Mole_Weight,					// 0x0504 The average molecular weight of a molecule object. (FLOAT64)
	kCDXProp_Frag_ConnectionOrder,			// 0x0505 An ordered list of attachment points within a fragment. (CDXObjectIDArray)

	// Bond properties.
	kCDXProp_Bond_Order = 0x0600,			// 0x0600 The order of a bond object. (INT16)
	kCDXProp_Bond_Display,					// 0x0601 The display type of a bond object. (INT16)
	kCDXProp_Bond_Display2,					// 0x0602 The display type for the second line of a double bond. (INT16)
	kCDXProp_Bond_DoublePosition,			// 0x0603 The position of the second line of a double bond. (INT16)
	kCDXProp_Bond_Begin,					// 0x0604 The ID of the CDX node object at the first end of a bond. (CDXObjectID)
	kCDXProp_Bond_End,						// 0x0605 The ID of the CDX node object at the second end of a bond. (CDXObjectID)
	kCDXProp_Bond_RestrictTopology,			// 0x0606 Indicates the desired topology of a bond in a query. (INT8)
	kCDXProp_Bond_RestrictRxnParticipation,	// 0x0607 Specifies that a bond is affected by a reaction. (INT8)
	kCDXProp_Bond_BeginAttach,				// 0x0608 Indicates where within the Bond_Begin node a bond is attached. (UINT8)
	kCDXProp_Bond_EndAttach,				// 0x0609 Indicates where within the Bond_End node a bond is attached. (UINT8)
	kCDXProp_Bond_CIPStereochemistry,		// 0x060A The bond's absolute stereochemistry according to the Cahn-Ingold-Prelog system. (INT8)
	kCDXProp_Bond_BondOrdering,				// 0x060B Ordered list of attached bond IDs. (CDXObjectIDArray)
	kCDXProp_Bond_ShowQuery,				// 0x060C Show the query indicator if non-zero. (CDXBoolean)
	kCDXProp_Bond_ShowStereo,				// 0x060D Show the stereochemistry indicator if non-zero. (CDXBoolean)
	kCDXProp_Bond_CrossingBonds,			// 0x060E Unordered list of IDs of bonds that cross this one (either above or below). (CDXObjectIDArray)
	kCDXProp_Bond_ShowRxn,					// 0x060F Show the reaction-change indicator if non-zero. (CDXBoolean)

	// Text properties.
	kCDXProp_Text = 0x0700,					// 0x0700 The text of a text object. (CDXString)
	kCDXProp_Justification,					// 0x0701 The horizontal justification of a text object. (INT8)
	kCDXProp_LineHeight,					// 0x0702 The line height of a text object. (UINT16)
	kCDXProp_WordWrapWidth,					// 0x0703 The word-wrap width of a text object. (INT16)
	kCDXProp_LineStarts,					// 0x0704 The number of lines of a text object followed by that many values indicating the zero-based text position of each line start. (INT16ListWithCounts)
	kCDXProp_LabelAlignment,				// 0x0705 The alignment of the text with respect to the node position. (INT8)
	kCDXProp_LabelLineHeight,				// 0x0706 Text line height for atom labels (INT16)
	kCDXProp_CaptionLineHeight,				// 0x0707 Text line height for non-atomlabel text objects (INT16)
	kCDXProp_InterpretChemically,			// 0x0708 Signifies whether to the text label should be interpreted chemically (if possible). (CDXBooleanImplied)

	// Document properties.
	kCDXProp_MacPrintInfo = 0x0800,			// 0x0800 The 120 byte Macintosh TPrint data associated with the CDX document object. Refer to Macintosh Toolbox manual for detailed description. (Unformatted)
	kCDXProp_WinPrintInfo,					// 0x0801 The Windows DEVMODE structure associated with the CDX document object. (Unformatted)
	kCDXProp_PrintMargins,					// 0x0802 The outer margins of the Document. (CDXRectangle)
	kCDXProp_ChainAngle,					// 0x0803 The default chain angle setting in degrees * 65536. (INT32)
	kCDXProp_BondSpacing,					// 0x0804 The spacing between segments of a multiple bond, measured relative to bond length. (INT16)
	kCDXProp_BondLength,					// 0x0805 The default bond length. (CDXCoordinate)
	kCDXProp_BoldWidth,						// 0x0806 The default bold bond width. (CDXCoordinate)
	kCDXProp_LineWidth,						// 0x0807 The default line width. (CDXCoordinate)
	kCDXProp_MarginWidth,					// 0x0808 The default amount of space surrounding atom labels. (CDXCoordinate)
	kCDXProp_HashSpacing,					// 0x0809 The default spacing between hashed lines used in wedged hashed bonds. (CDXCoordinate)
	kCDXProp_LabelStyle,					// 0x080A The default style for atom labels. (CDXFontStyle)
	kCDXProp_CaptionStyle,					// 0x080B The default style for non-atomlabel text objects. (CDXFontStyle)
	kCDXProp_CaptionJustification,			// 0x080C The horizontal justification of a caption (non-atomlabel text object) (INT8)
	kCDXProp_FractionalWidths,				// 0x080D Signifies whether to use fractional width information when drawing text. (CDXBooleanImplied)
	kCDXProp_Magnification,					// 0x080E The view magnification factor (INT16)
	kCDXProp_WidthPages,					// 0x080F The width of the document in pages. (INT16)
	kCDXProp_HeightPages,					// 0x0810 The height of the document in pages. (INT16)
	kCDXProp_DrawingSpaceType,				// 0x0811 The type of drawing space used for this document. (INT8)
	kCDXProp_Width,							// 0x0812 The width of an object in CDX coordinate units, possibly in a rotated or skewed frame. (CDXCoordinate)
	kCDXProp_Height,						// 0x0813 The height of an object in CDX coordinate units, possibly in a rotated or skewed frame. (CDXCoordinate)
	kCDXProp_PageOverlap,					// 0x0814 The amount of overlap of pages when a poster is tiled. (CDXCoordinate)
	kCDXProp_Header,						// 0x0815 The text of the header. (CDXString)
	kCDXProp_HeaderPosition,				// 0x0816 The vertical offset of the header baseline from the top of the page. (CDXCoordinate)
	kCDXProp_Footer,						// 0x0817 The text of the footer. (CDXString)
	kCDXProp_FooterPosition,				// 0x0818 The vertical offset of the footer baseline from the bottom of the page. (CDXCoordinate)
	kCDXProp_PrintTrimMarks,				// 0x0819 If present, trim marks are to printed in the margins. (CDXBooleanImplied)
	kCDXProp_LabelStyleFont,				// 0x081A The default font family for atom labels. (INT16)
	kCDXProp_CaptionStyleFont,				// 0x081B The default font style for captions (non-atom-label text objects). (INT16)
	kCDXProp_LabelStyleSize,				// 0x081C The default font size for atom labels. (INT16)
	kCDXProp_CaptionStyleSize,				// 0x081D The default font size for captions (non-atom-label text objects). (INT16)
	kCDXProp_LabelStyleFace,				// 0x081E The default font style for atom labels. (INT16)
	kCDXProp_CaptionStyleFace,				// 0x081F The default font face for captions (non-atom-label text objects). (INT16)
	kCDXProp_LabelStyleColor,				// 0x0820 The default color for atom labels (INT16)
	kCDXProp_CaptionStyleColor,				// 0x0821 The default color for captions (non-atom-label text objects). (INT16)
	kCDXProp_BondSpacingAbs,				// 0x0822 The absolute distance between segments of a multiple bond. (CDXCoordinate)
	kCDXProp_LabelJustification,			// 0x0823 The default justification for atom labels. (INT8)
	kCDXProp_FixInplaceExtent,				// 0x0824 Defines a size for OLE In-Place editing. (CDXPoint2D)
	kCDXProp_Side,							// 0x0825 A specific side of an object (rectangle). (INT16)
	kCDXProp_FixInplaceGap,					// 0x0826 Defines a padding for OLE In-Place editing. (CDXPoint2D)

	// Window properties.
	kCDXProp_Window_IsZoomed = 0x0900,		// 0x0900 Signifies whether the main viewing window is zoomed (maximized). (CDXBooleanImplied)
	kCDXProp_Window_Position,				// 0x0901 The top-left position of the main viewing window. (CDXPoint2D)
	kCDXProp_Window_Size,					// 0x0902 Height and width of the document window. (CDXPoint2D)
	
	// Graphic object properties.
	kCDXProp_Graphic_Type = 0x0A00,			// 0x0A00 The type of graphical object. (INT16)
	kCDXProp_Line_Type,						// 0x0A01 The type of a line object. (INT16)
	kCDXProp_Arrow_Type,					// 0x0A02 The type of arrow object, which represents line, arrow, arc, rectangle, or orbital. (INT16)
	kCDXProp_Rectangle_Type,				// 0x0A03 The type of a rectangle object. (INT16)
	kCDXProp_Oval_Type,						// 0x0A04 The type of an arrow object that represents a circle or ellipse. (INT16)
	kCDXProp_Orbital_Type,					// 0x0A05 The type of orbital object. (INT16)
	kCDXProp_Bracket_Type,					// 0x0A06 The type of symbol object. (INT16)
	kCDXProp_Symbol_Type,					// 0x0A07 The type of symbol object. (INT16)
	kCDXProp_Curve_Type,					// 0x0A08 The type of curve object. (INT16)
	kCDXProp_Arrow_HeadSize = 0x0A20,		// 0x0A20 The size of the arrow's head. (INT16)
	kCDXProp_Arc_AngularSize,				// 0x0A21 The size of an arc (in degrees * 10, so 90 degrees = 900). (INT16)
	kCDXProp_Bracket_LipSize,				// 0x0A22 The size of a bracket. (INT16)
	kCDXProp_Curve_Points,					// 0x0A23 The B&eacute;zier curve's control point locations. (CDXCurvePoints)
	kCDXProp_Bracket_Usage,					// 0x0A24 The syntactical chemical meaning of the bracket (SRU, mer, mon, xlink, etc). (INT8)
	kCDXProp_Polymer_RepeatPattern,			// 0x0A25 The head-to-tail connectivity of objects contained within the bracket. (INT8)
	kCDXProp_Polymer_FlipType,				// 0x0A26 The flip state of objects contained within the bracket. (INT8)
	kCDXProp_BracketedObjects,				// 0x0A27 The set of objects contained in a BracketedGroup. (CDXObjectIDArray)
	kCDXProp_Bracket_RepeatCount,			// 0x0A28 The number of times a multiple-group BracketedGroup is repeated. (INT16)
	kCDXProp_Bracket_ComponentOrder,		// 0x0A29 The component order associated with a BracketedGroup. (INT16)
	kCDXProp_Bracket_SRULabel,				// 0x0A2A The label associated with a BracketedGroup that represents an SRU. (CDXString)
	kCDXProp_Bracket_GraphicID,				// 0x0A2B The ID of a graphical object (bracket, brace, or parenthesis) associated with a Bracket Attachment. (CDXObjectID)
	kCDXProp_Bracket_BondID,				// 0x0A2C The ID of a bond that crosses a Bracket Attachment. (CDXObjectID)
	kCDXProp_Bracket_InnerAtomID,			// 0x0A2D The ID of the node located within the Bracketed Group and attached to a bond that crosses a Bracket Attachment. (CDXObjectID)
	kCDXProp_Curve_Points3D,				// 0x0A2E The B&eacute;zier curve's control point locations. (CDXCurvePoints3D)

	// Embedded pictures.
	kCDXProp_Picture_Edition = 0x0A60,		// 0x0A60 The section information (SectionHandle) of the Macintosh Publish & Subscribe edition embedded in the CDX picture object. (Unformatted)
	kCDXProp_Picture_EditionAlias,			// 0x0A61 The alias information of the Macintosh Publish & Subscribe edition embedded in the CDX picture object. (Unformatted)
	kCDXProp_MacPICT,						// 0x0A62 A Macintosh PICT data object. (Unformatted)
	kCDXProp_WindowsMetafile,				// 0x0A63 A Microsoft Windows Metafile object. (Unformatted)
	kCDXProp_OLEObject,						// 0x0A64 An OLE object. (Unformatted)
	kCDXProp_EnhancedMetafile,				// 0x0A65 A Microsoft Windows Enhanced Metafile object. (Unformatted)

	// Spectrum properties
	kCDXProp_Spectrum_XSpacing = 0x0A80,	// 0x0A80 The spacing in logical units (ppm, Hz, wavenumbers) between points along the X-axis of an evenly-spaced grid. (FLOAT64)
	kCDXProp_Spectrum_XLow,					// 0x0A81 The first data point for the X-axis of an evenly-spaced grid. (FLOAT64)
	kCDXProp_Spectrum_XType,				// 0x0A82 The type of units the X-axis represents. (INT16)
	kCDXProp_Spectrum_YType,				// 0x0A83 The type of units the Y-axis represents. (INT16)
	kCDXProp_Spectrum_XAxisLabel,			// 0x0A84 A label for the X-axis. (CDXString)
	kCDXProp_Spectrum_YAxisLabel,			// 0x0A85 A label for the Y-axis. (CDXString)
	kCDXProp_Spectrum_DataPoint,			// 0x0A86 The Y-axis values for the spectrum. It is an array of double values corresponding to X-axis values. (FLOAT64)
	kCDXProp_Spectrum_Class,				// 0x0A87 The type of spectrum represented. (INT16)
	kCDXProp_Spectrum_YLow,					// 0x0A88 Y value to be used to offset data when storing XML. (FLOAT64)
	kCDXProp_Spectrum_YScale,				// 0x0A89 Y scaling used to scale data when storing XML. (FLOAT64)

	// TLC properties
	kCDXProp_TLC_OriginFraction = 0x0AA0,	// 0x0AA0 The distance of the origin line from the bottom of a TLC Plate, as a fraction of the total height of the plate. (FLOAT64)
	kCDXProp_TLC_SolventFrontFraction,		// 0x0AA1 The distance of the solvent front from the top of a TLC Plate, as a fraction of the total height of the plate. (FLOAT64)
	kCDXProp_TLC_ShowOrigin,				// 0x0AA2 Show the origin line near the base of the TLC Plate if non-zero. (CDXBoolean)
	kCDXProp_TLC_ShowSolventFront,			// 0x0AA3 Show the solvent front line near the top of the TLC Plate if non-zero. (CDXBoolean)
	kCDXProp_TLC_ShowBorders,				// 0x0AA4 Show borders around the edges of the TLC Plate if non-zero. (CDXBoolean)
	kCDXProp_TLC_ShowSideTicks,				// 0x0AA5 Show tickmarks up the side of the TLC Plate if non-zero. (CDXBoolean)
	kCDXProp_TLC_Rf = 0x0AB0,				// 0x0AB0 The Retention Factor of an individual spot. (FLOAT64)
	kCDXProp_TLC_Tail,						// 0x0AB1 The length of the "tail" of an individual spot. (CDXCoordinate)
	kCDXProp_TLC_ShowRf,					// 0x0AB2 Show the spot's Retention Fraction (Rf) value if non-zero. (CDXBoolean)

	// Alternate Group properties
	kCDXProp_NamedAlternativeGroup_TextFrame = 0x0B00,	// 0x0B00 The bounding box of upper portion of the Named Alternative Group, containing the name of the group. (CDXRectangle)
	kCDXProp_NamedAlternativeGroup_GroupFrame,			// 0x0B01 The bounding box of the lower portion of the Named Alternative Group, containing the definition of the group. (CDXRectangle)
	kCDXProp_NamedAlternativeGroup_Valence,				// 0x0B02 The number of attachment points in each alternative in a named alternative group. (INT16)

	// Geometry and Constraint properties
	kCDXProp_GeometricFeature = 0x0B80,		// 0x0B80 The type of the geometrical feature (point, line, plane, etc.). (INT8)
	kCDXProp_RelationValue,					// 0x0B81 The numeric relationship (if any) among the basis objects used to define this object. (INT8)
	kCDXProp_BasisObjects,					// 0x0B82 An ordered list of objects used to define this object. (CDXObjectIDArray)
	kCDXProp_ConstraintType,				// 0x0B83 The constraint type (distance or angle). (INT8)
	kCDXProp_ConstraintMin,					// 0x0B84 The minimum value of the constraint (FLOAT64)
	kCDXProp_ConstraintMax,					// 0x0B85 The maximum value of the constraint (FLOAT64)
	kCDXProp_IgnoreUnconnectedAtoms,		// 0x0B86 Signifies whether unconnected atoms should be ignored within the exclusion sphere. (CDXBooleanImplied)
	kCDXProp_DihedralIsChiral,				// 0x0B87 Signifies whether a dihedral is signed or unsigned. (CDXBooleanImplied)
	kCDXProp_PointIsDirected,				// 0x0B88 For a point based on a normal, signifies whether it is in a specific direction relative to the reference point. (CDXBooleanImplied)

	// Reaction properties
	kCDXProp_ReactionStep_Atom_Map = 0x0C00,// 0x0C00 Represents pairs of mapped atom IDs; each pair is a reactant atom mapped to to a product atom. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Reactants,		// 0x0C01 An order list of reactants present in the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Products,			// 0x0C02 An order list of products present in the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Plusses,			// 0x0C03 An ordered list of pluses used to separate components of the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Arrows,			// 0x0C04 An ordered list of arrows used to separate components of the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_ObjectsAboveArrow,// 0x0C05 An order list of objects above the arrow in the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_ObjectsBelowArrow,// 0x0C06 An order list of objects below the arrow in the Reaction Step. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Atom_Map_Manual,	// 0x0C07 Represents pairs of mapped atom IDs; each pair is a reactant atom mapped to to a product atom. (CDXObjectIDArray)
	kCDXProp_ReactionStep_Atom_Map_Auto,	// 0x0C08 Represents pairs of mapped atom IDs; each pair is a reactant atom mapped to to a product atom. (CDXObjectIDArray)

	// CDObjectTag properties
	kCDXProp_ObjectTag_Type = 0x0D00,		// 0x0D00 The tag's data type. (INT16)
	kCDXProp_Unused6,						// 0x0D01 obsolete (obsolete)
	kCDXProp_Unused7,						// 0x0D02 obsolete (obsolete)
	kCDXProp_ObjectTag_Tracking,			// 0x0D03 The tag will participate in tracking if non-zero. (CDXBoolean)
	kCDXProp_ObjectTag_Persistent,			// 0x0D04 The tag will be resaved to a CDX file if non-zero. (CDXBoolean)
	kCDXProp_ObjectTag_Value,				// 0x0D05 The value is a INT32, FLOAT64 or unformatted string depending on the value of ObjectTag_Type. (varies)
	kCDXProp_Positioning,					// 0x0D06 How the indicator should be positioned with respect to its containing object. (INT8)
	kCDXProp_PositioningAngle,				// 0x0D07 Angular positioning, in radians * 65536. (INT32)
	kCDXProp_PositioningOffset,				// 0x0D08 Offset positioning. (CDXPoint2D)

	// CDSequence properties
	kCDXProp_Sequence_Identifier = 0x0E00,	// 0x0E00 A unique (but otherwise random) identifier for a given Sequence object. (CDXString)

	// CDCrossReference properties
	kCDXProp_CrossReference_Container = 0x0F00,	// 0x0F00 An external object containing (as an embedded object) the document containing the Sequence object being referenced. (CDXString)
	kCDXProp_CrossReference_Document,		// 0x0F01 An external document containing the Sequence object being referenced. (CDXString)
	kCDXProp_CrossReference_Identifier,		// 0x0F02 A unique (but otherwise random) identifier for a given Cross-Reference object. (CDXString)
	kCDXProp_CrossReference_Sequence,		// 0x0F03 A value matching the SequenceIdentifier of the Sequence object to be referenced. (CDXString)

	// Miscellaneous properties.
	kCDXProp_Template_PaneHeight = 0x1000,	// 0x1000 The height of the viewing window of a template grid. (CDXCoordinate)
	kCDXProp_Template_NumRows,				// 0x1001 The number of rows of the CDX TemplateGrid object. (INT16)
	kCDXProp_Template_NumColumns,			// 0x1002 The number of columns of the CDX TemplateGrid object. (INT16)

	kCDXProp_Group_Integral = 0x1100,		// 0x1100 The group is considered to be integral (non-subdivisible) if non-zero. (CDXBoolean)

	kCDXProp_SplitterPositions = 0x1ff0,	// 0x1FF0 An array of vertical positions that subdivide a page into regions. (CDXObjectIDArray)
	kCDXProp_PageDefinition,				// 0x1FF1 An array of vertical positions that subdivide a page into regions. (CDXObjectIDArray)

	// User defined properties
	// First 1024 tags are reserved for temporary tags used only during the runtime.
	kCDXUser_TemporaryBegin = kCDXTag_UserDefined,
	kCDXUser_TemporaryEnd = kCDXTag_UserDefined + 0x0400,

	// Objects.
	kCDXObj_Document = kCDXTag_Object,	// 0x8000
	kCDXObj_Page,						// 0x8001
	kCDXObj_Group,						// 0x8002
	kCDXObj_Fragment,					// 0x8003
	kCDXObj_Node,						// 0x8004
	kCDXObj_Bond,						// 0x8005
	kCDXObj_Text,						// 0x8006
	kCDXObj_Graphic,					// 0x8007
	kCDXObj_Curve,						// 0x8008
	kCDXObj_EmbeddedObject,				// 0x8009
	kCDXObj_NamedAlternativeGroup,		// 0x800a
	kCDXObj_TemplateGrid,				// 0x800b
	kCDXObj_RegistryNumber,				// 0x800c
	kCDXObj_ReactionScheme,				// 0x800d
	kCDXObj_ReactionStep,				// 0x800e
	kCDXObj_ObjectDefinition,			// 0x800f
	kCDXObj_Spectrum,					// 0x8010
	kCDXObj_ObjectTag,					// 0x8011
	kCDXObj_OleClientItem,				// 0x8012	// obsolete
	kCDXObj_Sequence,                   // 0x8013
	kCDXObj_CrossReference,             // 0x8014
	kCDXObj_Splitter,				    // 0x8015
	kCDXObj_Table,					    // 0x8016
	kCDXObj_BracketedGroup,				// 0x8017
	kCDXObj_BracketAttachment,			// 0x8018
	kCDXObj_CrossingBond,				// 0x8019
	kCDXObj_Border,						// 0x8020
	kCDXObj_Geometry,					// 0x8021
	kCDXObj_Constraint,					// 0x8022
	kCDXObj_TLCPlate,					// 0x8023
	kCDXObj_TLCLane,					// 0x8024
	kCDXObj_TLCSpot,					// 0x8025
	// Add new objects here
	kCDXObj_UnknownObject = 0x8FFF
};

enum CDXNodeType {
	kCDXNodeType_Unspecified,
	kCDXNodeType_Element,
	kCDXNodeType_ElementList,
	kCDXNodeType_ElementListNickname,
	kCDXNodeType_Nickname,
	kCDXNodeType_Fragment,
	kCDXNodeType_Formula,
	kCDXNodeType_GenericNickname,
	kCDXNodeType_AnonymousAlternativeGroup,
	kCDXNodeType_NamedAlternativeGroup,
	kCDXNodeType_MultiAttachment,
	kCDXNodeType_VariableAttachment,
	kCDXNodeType_ExternalConnectionPoint,
	kCDXNodeType_LinkNode
};

enum CDXLabelDisplay {
	kCDXLabelDisplay_Auto,
	kCDXLabelDisplay_Left,
	kCDXLabelDisplay_Center,
	kCDXLabelDisplay_Right,
	kCDXLabelDisplay_Above,
	kCDXLabelDisplay_Below,
	kCDXLabelDisplay_BestInitial
};

enum CDXRadical {	// Same as MDL codes
	kCDXRadical_None				= 0,
	kCDXRadical_Singlet				= 1,	// diradical singlet  (two dots)
	kCDXRadical_Doublet				= 2,	// monoradical		  (one dot)
	kCDXRadical_Triplet				= 3		// diradical triplet  (two dots)
};

enum CDXIsotope {
	kCDXIsotope_Natural				= 0
};

enum CDXRingBondCount {
	kCDXRingBondCount_Unspecified	= -1,
	kCDXRingBondCount_NoRingBonds	= 0,
	kCDXRingBondCount_AsDrawn		= 1,
	kCDXRingBondCount_SimpleRing	= 2,
	kCDXRingBondCount_Fusion		= 3,
	kCDXRingBondCount_SpiroOrHigher	= 4
};

enum CDXUnsaturation {
	kCDXUnsaturation_Unspecified	= 0,
	kCDXUnsaturation_MustBeAbsent	= 1,
	kCDXUnsaturation_MustBePresent	= 2,
	kCDXUnsaturationLastEnum
};

enum CDXReactionStereo {
	kCDXReactionStereo_Unspecified	= 0,
	kCDXReactionStereo_Inversion	= 1,
	kCDXReactionStereo_Retention	= 2
};

enum CDXTranslation {
	kCDXTranslation_Equal	= 0,
	kCDXTranslation_Broad	= 1,
	kCDXTranslation_Narrow	= 2,
	kCDXTranslation_Any		= 3
};

enum CDXAbundance {
	kCDXAbundance_Unspecified	= 0,
	kCDXAbundance_Any			= 1,
	kCDXAbundance_Natural		= 2,
	kCDXAbundance_Enriched		= 3,
	kCDXAbundance_Deficient		= 4,
	kCDXAbundance_Nonnatural	= 5
};

enum CDXExternalConnectionType {
	kCDXExternalConnection_Unspecified	= 0,
	kCDXExternalConnection_Diamond		= 1,
	kCDXExternalConnection_Star			= 2,
	kCDXExternalConnection_PolymerBead	= 3,
	kCDXExternalConnection_Wavy			= 4
};

enum CDXAtomGeometry {
	kCDXAtomGeometry_Unknown				=  0,
	kCDXAtomGeometry_1Ligand				=  1,
	kCDXAtomGeometry_Linear					=  2,
	kCDXAtomGeometry_Bent					=  3,
	kCDXAtomGeometry_TrigonalPlanar			=  4,
	kCDXAtomGeometry_TrigonalPyramidal		=  5,
	kCDXAtomGeometry_SquarePlanar			=  6,
	kCDXAtomGeometry_Tetrahedral			=  7,
	kCDXAtomGeometry_TrigonalBipyramidal	=  8,
	kCDXAtomGeometry_SquarePyramidal		=  9,
	kCDXAtomGeometry_5Ligand				= 10,
	kCDXAtomGeometry_Octahedral				= 11,
	kCDXAtomGeometry_6Ligand				= 12,
	kCDXAtomGeometry_7Ligand				= 13,
	kCDXAtomGeometry_8Ligand				= 14,
	kCDXAtomGeometry_9Ligand				= 15,
	kCDXAtomGeometry_10Ligand				= 16
};

enum CDXBondOrder {
	kCDXBondOrder_Single		= 0x0001,
	kCDXBondOrder_Double		= 0x0002,
	kCDXBondOrder_Triple		= 0x0004,
	kCDXBondOrder_Quadruple		= 0x0008,
	kCDXBondOrder_Quintuple		= 0x0010,
	kCDXBondOrder_Sextuple		= 0x0020,
	kCDXBondOrder_Half			= 0x0040,
	kCDXBondOrder_OneHalf		= 0x0080,
	kCDXBondOrder_TwoHalf		= 0x0100,
	kCDXBondOrder_ThreeHalf		= 0x0200,
	kCDXBondOrder_FourHalf		= 0x0400,
	kCDXBondOrder_FiveHalf		= 0x0800,
	kCDXBondOrder_Dative		= 0x1000,
	kCDXBondOrder_Ionic			= 0x2000,
	kCDXBondOrder_Hydrogen		= 0x4000,
	kCDXBondOrder_ThreeCenter	= 0x8000,
	kCDXBondOrder_SingleOrDouble = kCDXBondOrder_Single | kCDXBondOrder_Double,
	kCDXBondOrder_SingleOrAromatic = kCDXBondOrder_Single | kCDXBondOrder_OneHalf,
	kCDXBondOrder_DoubleOrAromatic = kCDXBondOrder_Double | kCDXBondOrder_OneHalf,
	kCDXBondOrder_Any = -1
};
// Permit combination of CDXBondOrder values
inline CDXBondOrder &operator |= (CDXBondOrder &lhs, const CDXBondOrder &rhs)
{
	return lhs = CDXBondOrder(UINT32(lhs) | UINT32(rhs));
}

enum CDXBondDisplay {
	kCDXBondDisplay_Solid				=  0,
	kCDXBondDisplay_Dash				=  1,
	kCDXBondDisplay_Hash				=  2,
	kCDXBondDisplay_WedgedHashBegin		=  3,
	kCDXBondDisplay_WedgedHashEnd		=  4,
	kCDXBondDisplay_Bold				=  5,
	kCDXBondDisplay_WedgeBegin			=  6,
	kCDXBondDisplay_WedgeEnd			=  7,
	kCDXBondDisplay_Wavy				=  8,
	kCDXBondDisplay_HollowWedgeBegin	=  9,
	kCDXBondDisplay_HollowWedgeEnd		= 10,
	kCDXBondDisplay_WavyWedgeBegin		= 11,
	kCDXBondDisplay_WavyWedgeEnd		= 12,
	kCDXBondDisplay_Dot					= 13,
	kCDXBondDisplay_DashDot				= 14
};

enum CDXBondDoublePosition {
	kCDXBondDoublePosition_AutoCenter	= 0x0000,
	kCDXBondDoublePosition_AutoRight	= 0x0001,
	kCDXBondDoublePosition_AutoLeft		= 0x0002,
	kCDXBondDoublePosition_UserCenter	= 0x0100,
	kCDXBondDoublePosition_UserRight	= 0x0101,
	kCDXBondDoublePosition_UserLeft		= 0x0102
};	
		
enum CDXBondTopology {
	kCDXBondTopology_Unspecified	= 0,
	kCDXBondTopology_Ring			= 1,
	kCDXBondTopology_Chain			= 2,
	kCDXBondTopology_RingOrChain	= 3
};

enum CDXBondReactionParticipation {
	kCDXBondReactionParticipation_Unspecified		= 0,
	kCDXBondReactionParticipation_ReactionCenter	= 1,
	kCDXBondReactionParticipation_MakeOrBreak		= 2,
	kCDXBondReactionParticipation_ChangeType		= 3,
	kCDXBondReactionParticipation_MakeAndChange		= 4,
	kCDXBondReactionParticipation_NotReactionCenter	= 5,
	kCDXBondReactionParticipation_NoChange			= 6,
	kCDXBondReactionParticipation_Unmapped			= 7
};

enum CDXTextJustification {
	kCDXTextJustification_Right = -1,
	kCDXTextJustification_Left,
	kCDXTextJustification_Center,
	kCDXTextJustification_Full,
	kCDXTextJustification_Above,
	kCDXTextJustification_Below,
	kCDXTextJustification_Auto,
	kCDXTextJustification_BestInitial
};

#define kCDXTagType_Unknown			"unknown"
#define kCDXTagType_Query			"query"
#define kCDXTagType_Rxn				"reaction"
#define kCDXTagType_Stereo			"stereo"
#define kCDXTagType_Number			"number"
#define kCDXTagType_Heading			"heading"
#define kCDXTagType_IDTerm			"idterm"
#define kCDXTagType_BracketUsage	"bracketusage"
#define kCDXTagType_PolymerRepeat	"polymerrepeat"
#define kCDXTagType_PolymerFlip		"polymerflip"
#define kCDXTagType_Deviation		"deviation"
#define kCDXTagType_Distance		"distance"
#define kCDXTagType_Angle			"angle"
#define kCDXTagType_Rf				"rf"

enum CDXPositioningType {
	kCDXPositioningType_Auto = 0,
	kCDXPositioningType_Angle,
	kCDXPositioningType_Offset,
	kCDXPositioningType_Absolute
};

enum CDXPageDefinition {
	kCDXPageDefinition_Undefined = 0,
	kCDXPageDefinition_Center,
	kCDXPageDefinition_TL4,
	kCDXPageDefinition_IDTerm,
	kCDXPageDefinition_FlushLeft,
	kCDXPageDefinition_FlushRight,
	kCDXPageDefinition_Reaction1,
	kCDXPageDefinition_Reaction2,
	kCDXPageDefinition_MulticolumnTL4,
	kCDXPageDefinition_MulticolumnNonTL4,
	kCDXPageDefinition_UserDefined
};

#define kCDXLineHeight_Variable  0
#define kCDXLineHeight_Automatic 1

enum CDXGraphicType {
	kCDXGraphicType_Undefined = 0,
	kCDXGraphicType_Line,
	kCDXGraphicType_Arc,
	kCDXGraphicType_Rectangle,
	kCDXGraphicType_Oval,
	kCDXGraphicType_Orbital,
	kCDXGraphicType_Bracket,
	kCDXGraphicType_Symbol
};

enum CDXBracketType
{
	kCDXBracketType_RoundPair,
	kCDXBracketType_SquarePair,
	kCDXBracketType_CurlyPair,
	kCDXBracketType_Square,
	kCDXBracketType_Curly,
	kCDXBracketType_Round
};

enum CDXRectangleType
{
	kCDXRectangleType_Plain = 0x0000,
	kCDXRectangleType_RoundEdge = 0x0001,
	kCDXRectangleType_Shadow = 0x0002,
	kCDXRectangleType_Shaded = 0x0004,
	kCDXRectangleType_Filled = 0x0008,
	kCDXRectangleType_Dashed = 0x0010,
	kCDXRectangleType_Bold = 0x0020
};

enum CDXOvalType
{
	kCDXOvalType_Circle = 0x0001,
	kCDXOvalType_Shaded = 0x0002,
	kCDXOvalType_Filled = 0x0004,
	kCDXOvalType_Dashed = 0x0008,
	kCDXOvalType_Bold   = 0x0010,
	kCDXOvalType_Shadowed   = 0x0020
};

enum CDXSymbolType
{
	kCDXSymbolType_LonePair,
	kCDXSymbolType_Electron,
	kCDXSymbolType_RadicalCation,
	kCDXSymbolType_RadicalAnion,
	kCDXSymbolType_CirclePlus,
	kCDXSymbolType_CircleMinus,
	kCDXSymbolType_Dagger,
	kCDXSymbolType_DoubleDagger,
	kCDXSymbolType_Plus,
	kCDXSymbolType_Minus,
	kCDXSymbolType_Racemic,
	kCDXSymbolType_Absolute,
	kCDXSymbolType_Relative
};

enum CDXLineType
{
	kCDXLineType_Solid  = 0x0000,
	kCDXLineType_Dashed = 0x0001,
	kCDXLineType_Bold	= 0x0002,
	kCDXLineType_Wavy	= 0x0004
};

enum CDXArrowType
{
	kCDXArrowType_NoHead			=  0,
	kCDXArrowType_HalfHead			=  1,
	kCDXArrowType_FullHead			=  2,
	kCDXArrowType_Resonance			=  4,
	kCDXArrowType_Equilibrium		=  8,
	kCDXArrowType_Hollow			= 16,
	kCDXArrowType_RetroSynthetic	= 32
};

enum CDXOrbitalType
{
	kCDXOrbitalType_s,					// s orbital
	kCDXOrbitalType_oval,				// Oval-shaped sigma or pi orbital
	kCDXOrbitalType_lobe,				// One lobe of a p orbital
	kCDXOrbitalType_p,					// Complete p orbital
	kCDXOrbitalType_hybridPlus,			// hydrid orbital
	kCDXOrbitalType_hybridMinus,		// hydrid orbital (opposite shading)
	kCDXOrbitalType_dz2Plus,			// dz2 orbital
	kCDXOrbitalType_dz2Minus,			// dz2 orbital (opposite shading)
	kCDXOrbitalType_dxy,				// dxy orbital

	kCDXOrbitalType_sShaded = 0x0100,	// shaded s orbital
	kCDXOrbitalType_ovalShaded,			// shaded Oval-shaped sigma or pi orbital
	kCDXOrbitalType_lobeShaded,			// shaded single lobe of a p orbital
	kCDXOrbitalType_pShaded,			// shaded Complete p orbital
	
	kCDXOrbitalType_sFilled = 0x0200,	// filled s orbital
	kCDXOrbitalType_ovalFilled,			// filled Oval-shaped sigma or pi orbital
	kCDXOrbitalType_lobeFilled,			// filled single lobe of a p orbital
	kCDXOrbitalType_pFilled,			// filled Complete p orbital
	kCDXOrbitalType_hybridPlusFilled,	// filled hydrid orbital
	kCDXOrbitalType_hybridMinusFilled,	// filled hydrid orbital (opposite shading)
	kCDXOrbitalType_dz2PlusFilled,		// filled dz2 orbital
	kCDXOrbitalType_dz2MinusFilled,		// filled dz2 orbital (opposite shading)
	kCDXOrbitalType_dxyFilled			// filled dxy orbital

};

enum CDXBracketUsage
{
	kCDXBracketUsage_Unspecified = 0,
	kCDXBracketUsage_Anypolymer = 18,
	kCDXBracketUsage_Component = 13,
	kCDXBracketUsage_Copolymer = 6,
	kCDXBracketUsage_CopolymerAlternating = 7,
	kCDXBracketUsage_CopolymerBlock = 9,
	kCDXBracketUsage_CopolymerRandom = 8,
	kCDXBracketUsage_Crosslink = 10,
	kCDXBracketUsage_Generic = 17,
	kCDXBracketUsage_Graft = 11,
	kCDXBracketUsage_Mer = 5,
	kCDXBracketUsage_MixtureOrdered = 15,
	kCDXBracketUsage_MixtureUnordered = 14,
	kCDXBracketUsage_Modification = 12,
	kCDXBracketUsage_Monomer = 4,
	kCDXBracketUsage_MultipleGroup = 16,
	kCDXBracketUsage_SRU = 3,
	kCDXBracketUsage_Unused1 = 1,
	kCDXBracketUsage_Unused2 = 2
};

enum CDXPolymerRepeatPattern
{
	kCDXPolymerRepeatPattern_HeadToTail = 0,
	kCDXPolymerRepeatPattern_HeadToHead,
	kCDXPolymerRepeatPattern_EitherUnknown
};

enum CDXPolymerFlipType
{
	kCDXPolymerFlipType_Unspecified = 0,
	kCDXPolymerFlipType_NoFlip,
	kCDXPolymerFlipType_Flip
};

enum CDXSpectrumYType
{
	kCDXSpectrumYType_Unknown,
	kCDXSpectrumYType_Absorbance, 
	kCDXSpectrumYType_Transmittance, 
	kCDXSpectrumYType_PercentTransmittance, 
	kCDXSpectrumYType_Other, 
	kCDXSpectrumYType_ArbitraryUnits
};

enum CDXSpectrumXType
{
	kCDXSpectrumXType_Unknown,
	kCDXSpectrumXType_Wavenumbers,
	kCDXSpectrumXType_Microns,
	kCDXSpectrumXType_Hertz,
	kCDXSpectrumXType_MassUnits,
	kCDXSpectrumXType_PartsPerMillion,
	kCDXSpectrumXType_Other
};

enum CDXSpectrumClass
{
	kCDXSpectrumClass_Unknown,
	kCDXSpectrumClass_Chromatogram,
	kCDXSpectrumClass_Infrared,
	kCDXSpectrumClass_UVVis,
	kCDXSpectrumClass_XRayDiffraction,
	kCDXSpectrumClass_MassSpectrum,
	kCDXSpectrumClass_NMR,
	kCDXSpectrumClass_Raman,
	kCDXSpectrumClass_Fluorescence,
	kCDXSpectrumClass_Atomic
};

enum CDXDrawingSpaceType
{
	kCDXDrawingSpace_Pages,
	kCDXDrawingSpace_Poster
};

enum CDXAtomCIPType
{
	kCDXCIPAtom_Undetermined			= 0,
	kCDXCIPAtom_None,
	kCDXCIPAtom_R,
	kCDXCIPAtom_S,
	kCDXCIPAtom_r,
	kCDXCIPAtom_s,
	kCDXCIPAtom_Unspecified						// No hash/wedge, but if there were one, it would have stereochemistry.
};

enum CDXBondCIPType
{
	kCDXCIPBond_Undetermined			= 0,
	kCDXCIPBond_None,
	kCDXCIPBond_E,
	kCDXCIPBond_Z
};

enum CDXObjectTagType
{
	kCDXObjectTagType_Undefined			= 0,
	kCDXObjectTagType_Double,
	kCDXObjectTagType_Long,
	kCDXObjectTagType_String
};

enum CDXSideType
{
	kCDXSideType_Undefined				= 0,
	kCDXSideType_Top,
	kCDXSideType_Left,
	kCDXSideType_Bottom,
	kCDXSideType_Right
};

enum CDXGeometricFeature
{
	kCDXGeometricFeature_Undefined				= 0,
	kCDXGeometricFeature_PointFromPointPointDistance,
	kCDXGeometricFeature_PointFromPointPointPercentage,
	kCDXGeometricFeature_PointFromPointNormalDistance,
	kCDXGeometricFeature_LineFromPoints,
	kCDXGeometricFeature_PlaneFromPoints,
	kCDXGeometricFeature_PlaneFromPointLine,
	kCDXGeometricFeature_CentroidFromPoints,
	kCDXGeometricFeature_NormalFromPointPlane
};

enum CDXConstraintType
{
	kCDXConstraintType_Undefined			= 0,
	kCDXConstraintType_Distance,
	kCDXConstraintType_Angle,
	kCDXConstraintType_ExclusionSphere
};

enum CDXCharSet
{
	kCDXCharSetUnknown = 0,
	kCDXCharSetEBCDICOEM = 37,
	kCDXCharSetMSDOSUS = 437,
	kCDXCharSetEBCDIC500V1 = 500,
	kCDXCharSetArabicASMO708 = 708,
	kCDXCharSetArabicASMO449P,
	kCDXCharSetArabicTransparent,
	kCDXCharSetArabicTransparentASMO = 720,
	kCDXCharSetGreek437G = 737,
	kCDXCharSetBalticOEM = 775,
	kCDXCharSetMSDOSLatin1 = 850,
	kCDXCharSetMSDOSLatin2 = 852,
	kCDXCharSetIBMCyrillic = 855,
	kCDXCharSetIBMTurkish = 857,
	kCDXCharSetMSDOSPortuguese = 860,
	kCDXCharSetMSDOSIcelandic,
	kCDXCharSetHebrewOEM,
	kCDXCharSetMSDOSCanadianFrench,
	kCDXCharSetArabicOEM,
	kCDXCharSetMSDOSNordic,
	kCDXCharSetMSDOSRussian,
	kCDXCharSetIBMModernGreek = 869,
	kCDXCharSetThai = 874,
	kCDXCharSetEBCDIC,
	kCDXCharSetJapanese = 932,
	kCDXCharSetChineseSimplified = 936, // PRC, Singapore
	kCDXCharSetKorean = 949,
	kCDXCharSetChineseTraditional = 950, // Taiwan, Hong Kong
	kCDXCharSetUnicodeISO10646 = 1200,
	kCDXCharSetWin31EasternEuropean = 1250,
	kCDXCharSetWin31Cyrillic,
	kCDXCharSetWin31Latin1,
	kCDXCharSetWin31Greek,
	kCDXCharSetWin31Turkish,
	kCDXCharSetHebrew,
	kCDXCharSetArabic,
	kCDXCharSetBaltic,
	kCDXCharSetVietnamese,
	kCDXCharSetKoreanJohab = 1361,
	kCDXCharSetMacRoman = 10000,
	kCDXCharSetMacJapanese,
	kCDXCharSetMacTradChinese,
	kCDXCharSetMacKorean,
	kCDXCharSetMacArabic,
	kCDXCharSetMacHebrew,
	kCDXCharSetMacGreek,
	kCDXCharSetMacCyrillic,
	kCDXCharSetMacReserved,
	kCDXCharSetMacDevanagari,
	kCDXCharSetMacGurmukhi,
	kCDXCharSetMacGujarati,
	kCDXCharSetMacOriya,
	kCDXCharSetMacBengali,
	kCDXCharSetMacTamil,
	kCDXCharSetMacTelugu,
	kCDXCharSetMacKannada,
	kCDXCharSetMacMalayalam,
	kCDXCharSetMacSinhalese,
	kCDXCharSetMacBurmese,
	kCDXCharSetMacKhmer,
	kCDXCharSetMacThai,
	kCDXCharSetMacLao,
	kCDXCharSetMacGeorgian,
	kCDXCharSetMacArmenian,
	kCDXCharSetMacSimpChinese,
	kCDXCharSetMacTibetan,
	kCDXCharSetMacMongolian,
	kCDXCharSetMacEthiopic,
	kCDXCharSetMacCentralEuroRoman,
	kCDXCharSetMacVietnamese,
	kCDXCharSetMacExtArabic,
	kCDXCharSetMacUninterpreted,
	kCDXCharSetMacIcelandic = 10079,
	kCDXCharSetMacTurkish = 10081
};

#endif // _H_CDXConstants
