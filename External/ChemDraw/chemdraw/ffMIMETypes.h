// CommonCS/LibCommon/Hdr/ffMIMETypes.h
// Contains: Definitions for utility functions relating to the reading of foreign file formats
// Copyright © 1985-2004, CambridgeSoft Corp., All Rights Reserved

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
// The following are the preferred types that are supported by the common code and that can be freely exposed to customers
// See discussion of Chemical MIME types on the web at various places, including http://www.ch.ic.ac.uk/chemime/
//
// MIME types are part of international standards; this list cannot be extended arbitrarily.
// At the very least, an attempt should be made to contact the owner/author of any new file type
// and ask *them* if *they* have a preference for the correct MIME type used to describe their file type.
//
// The "x-" infix indicates provisional types that have not completed the full approval process
// to be recognized as "official" MIME types.  As of now, *all* MIME types starting with
// "chemical/" are provisional, and hence all of them start with "chemical/x-"
#define kffMimeTypeCDX          "chemical/x-cdx"
#define kffMimeTypeXML          "text/xml"
#define kffMimeTypeMDLMolfile   "chemical/x-mdl-molfile"
#define kffMimeTypeMDLMolV3000  "chemical/x-mdl-molfile-v3000"	// Used for output only; use the unversioned kffMimeTypeMDLMolfile for input, since the shared code will perform automatic version recognition. 
#define kffMimeTypeMDLSDfile    "chemical/x-mdl-sdfile"
#define kffMimeTypeMDLSDV3000   "chemical/x-mdl-sdfile-v3000"	// Used for output only; use the unversioned kffMimeTypeMDLSDfile for input, since the shared code will perform automatic version recognition.
#define kffMimeTypeMDLTGF       "chemical/x-mdl-tgf"
#define kffMimeTypeMDLRXN       "chemical/x-mdl-rxn"
#define kffMimeTypeMDLRXNV3000  "chemical/x-mdl-rxn-v3000"		// Used for output only; use the unversioned kffMimeTypeMDLRXN for input, since the shared code will perform automatic version recognition. 
#define kffMimeTypeSMILES       "chemical/x-daylight-smiles"
#define kffMimeTypeCANSMILES	"chemical/x-daylight-cansmiles"
#define kffMimeTypeSMILES2      "chemical/x-smiles"
#define kffMimeTypeMDLSketch    "chemical/x-mdl-isis"
#define kffMimeTypeMSIMolfile   "chemical/x-msi-molfile"
#define kffMimeTypeSMD          "chemical/x-smd"
#define kffMimeTypeConnTab      "chemical/x-ct"
#define kffMimeTypeCML		    "chemical/x-cml"
#define kffMimeTypeInChI	    "chemical/x-inchi"
#define kffMimeTypeInChIKey	    "chemical/x-inchikey"

// The following are future-looking definitions for MIME types that may become valid some day
// It would be nice if we support them silently in case they do become valid,
// but as of now they MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeCDX1         "chemical/cdx"					// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLMolfile1  "chemical/mdl-molfile"			// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLMolV3000a "chemical/mdl-molfile-v3000"	// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLSDfile1   "chemical/mdl-sdfile"			// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLSDV3000a  "chemical/mdl-sdfile-v3000"		// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLTGF1      "chemical/mdl-tgf"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLRXN1      "chemical/mdl-rxn"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLRXNV3000a "chemical/mdl-rxn-v3000"		// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeSMILES1      "chemical/daylight-smiles"		// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeSMILES3      "chemical/smiles"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMDLSketch1   "chemical/mdl-isis"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeMSIMolfile1  "chemical/msi-molfile"			// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeSMD1         "chemical/smd"					// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeConnTab1     "chemical/ct"					// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeCML1	        "chemical/cml"					// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeInChI1       "chemical/inchi"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS
#define kffMimeTypeInChIKey1    "chemical/inchikey"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS


// When Stew created CDXML, he said that "text/xml" was appropriate and that it wasn't necessary
// to create a new CDXML-specific MIME type.  That hasn't stopped people outside of CambridgeSoft
// from creating this MIME type on their own and using it, so I guess we have to recognize it also.
// The Official CambridgeSoft Policy remains that "text/xml is the preferred MIME type for CDXML files.
#define kffMimeTypeCDXML        "chemical/x-cdxml"				// Should not be exposed to customers				
#define kffMimeTypeCDXML1       "chemical/cdxml"				// *** Future-looking MIME type definition.  MUST NOT BE EXPOSED TO CUSTOMERS				

// This is the MIME type for the old CHM file format used by versions of ChemDraw earlier than ChemDraw 4.0.
// The CHM format is primarily a graphical format and does not include enough chemical information to be reliably
// interpretable by shared code.  A partial conversion is provided in shared code, but the full conversion is
// available only by using the ChemDraw application (or plugin or ActiveX control)
#define kffMimeTypeChemDraw     "chemical/x-chemdraw"			// Can only reliably be read by ChemDraw directly.  Only partial support offered in shared code


// The following MIME types are technically valid,
// but are not currently supported for reading or writing in any CambridgeSoft application
#define kffMimeTypeF1           "chemical/x-questel-f1"			// Not supported in any CambridgeSoft application
#define kffMimeTypeF1Q          "chemical/x-questel-f1-query"	// Not supported in any CambridgeSoft application

#define kffMimeTypePict             "image/x-pict"
#define kffMimeTypePictCompressed   "image/x-pict-compressed"
#define kffMimeTypeEMF              "image/x-emf"
#define kffMimeTypeEMF1             "image/emf"
#define kffMimeTypeEMF1Compressed   "image/emf-compressed"
#define kffMimeTypeWMF              "image/x-wmf"
#define kffMimeTypeWMF1             "image/wmf"
#define kffMimeTypeWMF1Compressed   "image/wmf-compressed"
#define kffMimeTypeWMF2             "image/x-wmfs"
#define kffMimeTypeGIF              "image/gif"
#define kffMimeTypeBMP              "image/bmp"
#define kffMimeTypeTIFF             "image/tiff"
#define kffMimeTypePNG              "image/x-png"
#define kffMimeTypePNG1             "image/png"
#define kffMimeTypeJPEG             "image/jpeg"
#define kffMimeTypePDF              "application/pdf"
#define kffMimeTypeOLE              "application/ole"
#define kffMimeTypeOLECompressed    "application/ole-compressed"

// The CDAX Const has a frozen interface, therefore to add new functionality to the interface we have to use existing API calls.
// The most flexible way is to send a special string through the CDAX set_Data() call and catch it in the CDAX's put_Data().
// This header is a convenient location for defining any such commands as they will be used with a mimetype.
#define kCDAXSetDataDisableMonomerEditorTag  "DISABLEMONOMEREDITOR"
