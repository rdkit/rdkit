//
// CDXEmbeddedObject.cpp
//
// Copyright � 2022 PerkinElmer, Inc. All rights reserved.

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

#include <ostream>

#include "CDXStdObjects.h"
#include "CDXMLNames.h"
#include "ffMIMETypes.h"

using std::stoul;

CDXEmbeddedObject::CDXEmbeddedObject(CDXObjectID id_arg)
 : CDXGraphicObject(kCDXObj_EmbeddedObject, id_arg)
 , m_rotationAngle(0)
{
}

CDXEmbeddedObject::CDXEmbeddedObject(const CDXEmbeddedObject &src)
	:	CDXGraphicObject(src)
    ,   CDX_INHERIT_SRC(mimeTypeToFormatInfoMap)
	,	CDX_INHERIT_SRC (m_rotationAngle)
{
}

void CDXEmbeddedObject::StoreAttribute(CDXDataSource &src_arg, CDXTag attribTag_arg, size_t size_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PDF:
        SetFormatData(kffMimeTypePDF, src_arg.GetString(size_arg));
		break;

	case kCDXProp_MacPICT:
        SetFormatData(kffMimeTypePict, src_arg.GetString(size_arg));
		break;

	case kCDXProp_WindowsMetafile:
        SetFormatData(kffMimeTypeWMF1, src_arg.GetString(size_arg));
		break;

	case kCDXProp_EnhancedMetafile:
        SetFormatData(kffMimeTypeEMF1, src_arg.GetString(size_arg));
		break;

	case kCDXProp_OLEObject:
        SetFormatData(kffMimeTypeOLE, src_arg.GetString(size_arg));
		break;

	case kCDXProp_Compressed_MacPICT:
        SetFormatData(kffMimeTypePictCompressed, src_arg.GetString(size_arg));
		break;

	case kCDXProp_Compressed_WindowsMetafile:
        SetFormatData(kffMimeTypeWMF1Compressed, src_arg.GetString(size_arg));
		break;

	case kCDXProp_Compressed_EnhancedMetafile:
        SetFormatData(kffMimeTypeEMF1Compressed, src_arg.GetString(size_arg));
		break;

	case kCDXProp_Compressed_OLEObject:
        SetFormatData(kffMimeTypeOLECompressed, src_arg.GetString(size_arg));
		break;

	case kCDXProp_Uncompressed_MacPICT_Size:
        SetFormatUncompressedSize(kffMimeTypePictCompressed, src_arg.GetINT(size_arg));
		break;

	case kCDXProp_Uncompressed_WindowsMetafile_Size:
        SetFormatUncompressedSize(kffMimeTypeWMF1Compressed, src_arg.GetINT(size_arg));
		break;

	case kCDXProp_Uncompressed_EnhancedMetafile_Size:
        SetFormatUncompressedSize(kffMimeTypeEMF1Compressed, src_arg.GetINT(size_arg));
		break;

	case kCDXProp_Uncompressed_OLEObject_Size:
        SetFormatUncompressedSize(kffMimeTypeOLECompressed, src_arg.GetINT(size_arg));
		break;

	case kCDXProp_GIF:
        SetFormatData(kffMimeTypeGIF, src_arg.GetString(size_arg));
		break;

	case kCDXProp_TIFF:
        SetFormatData(kffMimeTypeTIFF, src_arg.GetString(size_arg));
		break;

	case kCDXProp_PNG:
        SetFormatData(kffMimeTypePNG1, src_arg.GetString(size_arg));
		break;

	case kCDXProp_JPEG:
        SetFormatData(kffMimeTypeJPEG, src_arg.GetString(size_arg));
		break;

	case kCDXProp_BMP:
        SetFormatData(kffMimeTypeBMP, src_arg.GetString(size_arg));
		break;

	case kCDXProp_RotationAngle:
		m_rotationAngle = src_arg.GetUINT(size_arg);
		break;

	default:
		CDXGraphicObject::StoreAttribute(src_arg, attribTag_arg, size_arg);
		break;
	}
}

void CDXEmbeddedObject::XMLStoreAttribute(CDXTag attribTag_arg, const std::string &attribValue_arg)
{
	switch (attribTag_arg)
	{
	case kCDXProp_PDF:
        SetFormatData(kffMimeTypePDF, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_MacPICT:
        SetFormatData(kffMimeTypePict, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_WindowsMetafile:
        SetFormatData(kffMimeTypeWMF1, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_EnhancedMetafile:
        SetFormatData(kffMimeTypeEMF1, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_OLEObject:
        SetFormatData(kffMimeTypeOLE, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_Compressed_MacPICT:
        SetFormatData(kffMimeTypePictCompressed, csBase64ToBinary(attribValue_arg));
		break;

	case kCDXProp_Compressed_WindowsMetafile:
        SetFormatData(kffMimeTypeWMF1Compressed, csBase64ToBinary(attribValue_arg));
		break;

	case kCDXProp_Compressed_EnhancedMetafile:
        SetFormatData(kffMimeTypeEMF1Compressed, csBase64ToBinary(attribValue_arg));
		break;

	case kCDXProp_Compressed_OLEObject:
        SetFormatData(kffMimeTypeOLECompressed, csBase64ToBinary(attribValue_arg));
		break;

	case kCDXProp_Uncompressed_MacPICT_Size:
        SetFormatUncompressedSize(kffMimeTypePictCompressed, static_cast<uint32_t>(stoul(attribValue_arg.c_str())));
		break;

	case kCDXProp_Uncompressed_WindowsMetafile_Size:
        SetFormatUncompressedSize(kffMimeTypeWMF1Compressed, static_cast<uint32_t>(stoul(attribValue_arg.c_str())));
		break;

	case kCDXProp_Uncompressed_EnhancedMetafile_Size:
        SetFormatUncompressedSize(kffMimeTypeEMF1Compressed, static_cast<uint32_t>(stoul(attribValue_arg.c_str())));
		break;

	case kCDXProp_Uncompressed_OLEObject_Size:
        SetFormatUncompressedSize(kffMimeTypeOLECompressed, static_cast<uint32_t>(stoul(attribValue_arg.c_str())));
		break;

	case kCDXProp_GIF:
        SetFormatData(kffMimeTypeGIF, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_TIFF:
        SetFormatData(kffMimeTypeTIFF, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_PNG:
        SetFormatData(kffMimeTypePNG1, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_JPEG:
        SetFormatData(kffMimeTypeJPEG, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_BMP:
        SetFormatData(kffMimeTypeBMP, HexToBinary(attribValue_arg));
		break;

	case kCDXProp_RotationAngle:
        m_rotationAngle = (UINT32) std::strtoul(attribValue_arg.c_str(),NULL,10);
		break;

	default:
		CDXGraphicObject::XMLStoreAttribute(attribTag_arg, attribValue_arg);
		break;
	}
}

void CDXEmbeddedObject::WriteAttributesTo(CDXDataSink &sink_arg) const
{
	CDXGraphicObject::WriteAttributesTo(sink_arg);

    if (HasFormatData(kffMimeTypePict))
	{
		sink_arg.PutTag(kCDXProp_MacPICT);
        PutTagData(sink_arg, kffMimeTypePict);
	}
	
    if (HasFormatData(kffMimeTypePDF))
    {
		sink_arg.PutTag(kCDXProp_PDF);
        PutTagData(sink_arg, kffMimeTypePDF);
	}

    if (HasFormatData(kffMimeTypeWMF1))
    {
		sink_arg.PutTag(kCDXProp_WindowsMetafile);
        PutTagData(sink_arg, kffMimeTypeWMF1);
	}

    if (HasFormatData(kffMimeTypeEMF1))
    {
		sink_arg.PutTag(kCDXProp_EnhancedMetafile);
        PutTagData(sink_arg, kffMimeTypeEMF1);
	}

    if (HasFormatData(kffMimeTypeOLE))
    {
		sink_arg.PutTag(kCDXProp_OLEObject);
        PutTagData(sink_arg, kffMimeTypeOLE);
	}

    if (HasFormatData(kffMimeTypePictCompressed))
    {
		sink_arg.PutAttribute( kCDXProp_Uncompressed_MacPICT_Size, INT32(GetFormatUncompressedSize(kffMimeTypePictCompressed)));
		
        sink_arg.PutTag(kCDXProp_Compressed_MacPICT);
        PutTagData(sink_arg, kffMimeTypePictCompressed);
	}

    if (HasFormatData(kffMimeTypeWMF1Compressed))
    {
		sink_arg.PutAttribute( kCDXProp_Uncompressed_WindowsMetafile_Size, INT32(GetFormatUncompressedSize(kffMimeTypeWMF1Compressed)));
		
        sink_arg.PutTag(kCDXProp_Compressed_WindowsMetafile);
        PutTagData(sink_arg, kffMimeTypeWMF1Compressed);
	}

    if (HasFormatData(kffMimeTypeEMF1Compressed))
    {
		sink_arg.PutAttribute( kCDXProp_Uncompressed_EnhancedMetafile_Size, INT32(GetFormatUncompressedSize(kffMimeTypeEMF1Compressed)));
		
        sink_arg.PutTag(kCDXProp_Compressed_EnhancedMetafile);
        PutTagData(sink_arg, kffMimeTypeEMF1Compressed);
	}

    if (HasFormatData(kffMimeTypeOLECompressed))
    {
		sink_arg.PutAttribute( kCDXProp_Uncompressed_OLEObject_Size, INT32(GetFormatUncompressedSize(kffMimeTypeOLECompressed)) );
		
        sink_arg.PutTag(kCDXProp_Compressed_OLEObject);
        PutTagData(sink_arg, kffMimeTypeOLECompressed);
    }

    if (HasFormatData(kffMimeTypeGIF))
    {
		sink_arg.PutTag(kCDXProp_GIF);
        PutTagData(sink_arg, kffMimeTypeGIF);
	}

    if (HasFormatData(kffMimeTypeTIFF))
	{
		sink_arg.PutTag(kCDXProp_TIFF);
        PutTagData(sink_arg, kffMimeTypeTIFF);
	}

    if (HasFormatData(kffMimeTypePNG1))
    {
		sink_arg.PutTag(kCDXProp_PNG);
        PutTagData(sink_arg, kffMimeTypePNG1);
	}

    if (HasFormatData(kffMimeTypeJPEG))
    {
		sink_arg.PutTag(kCDXProp_JPEG);
        PutTagData(sink_arg, kffMimeTypeJPEG);
	}

    if (HasFormatData(kffMimeTypeBMP))
    {
		sink_arg.PutTag(kCDXProp_BMP);
        PutTagData(sink_arg, kffMimeTypeBMP);
	}

    if (m_rotationAngle != 0)
    {
        sink_arg.PutAttribute(kCDXProp_RotationAngle, m_rotationAngle);
    }
}

void CDXEmbeddedObject::XMLWriteAttributes(XMLDataSink &sink_arg) const
{
	CDXGraphicObject::XMLWriteAttributes(sink_arg);

    if (HasFormatData(kffMimeTypePict))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_MacPICT, BinaryToHex(GetFormatData(kffMimeTypePict)));
    }

    if (HasFormatData(kffMimeTypePDF))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_PDF, BinaryToHex(GetFormatData(kffMimeTypePDF)));
    }

    if (HasFormatData(kffMimeTypeWMF1))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_WindowsMetafile, BinaryToHex(GetFormatData(kffMimeTypeWMF1)));
    }

    if (HasFormatData(kffMimeTypeEMF1))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_EnhancedMetafile, BinaryToHex(GetFormatData(kffMimeTypeEMF1)));
    }

    if (HasFormatData(kffMimeTypeOLE))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_OLEObject, BinaryToHex(GetFormatData(kffMimeTypeOLE)));
    }

    if (HasFormatData(kffMimeTypePictCompressed))
    {
		CDXMLPutAttribute(sink_arg, kCDXProp_Compressed_MacPICT, (string)csBinaryToBase64(GetFormatData(kffMimeTypePictCompressed)) );
		CDXMLPutAttribute(sink_arg, kCDXProp_Uncompressed_MacPICT_Size, (INT32)GetFormatUncompressedSize(kffMimeTypePictCompressed));
	}

    if (HasFormatData(kffMimeTypeWMF1Compressed))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Compressed_WindowsMetafile, (string)csBinaryToBase64(GetFormatData(kffMimeTypeWMF1Compressed)) );
		CDXMLPutAttribute(sink_arg, kCDXProp_Uncompressed_WindowsMetafile_Size, (INT32)GetFormatUncompressedSize(kffMimeTypeWMF1Compressed));
	}

    if (HasFormatData(kffMimeTypeEMF1Compressed))
    {
		CDXMLPutAttribute(sink_arg, kCDXProp_Compressed_EnhancedMetafile, (string)csBinaryToBase64(GetFormatData(kffMimeTypeEMF1Compressed)) );
		CDXMLPutAttribute(sink_arg, kCDXProp_Uncompressed_EnhancedMetafile_Size, (INT32)GetFormatUncompressedSize(kffMimeTypeEMF1Compressed));
	}

    if (HasFormatData(kffMimeTypeOLECompressed))
	{
		CDXMLPutAttribute(sink_arg, kCDXProp_Compressed_OLEObject, (string)csBinaryToBase64(GetFormatData(kffMimeTypeOLECompressed)) );
		CDXMLPutAttribute(sink_arg, kCDXProp_Uncompressed_OLEObject_Size, (INT32)GetFormatUncompressedSize(kffMimeTypeOLECompressed));
	}

    if (HasFormatData(kffMimeTypeGIF))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_GIF, BinaryToHex(GetFormatData(kffMimeTypeGIF)));
    }

    if (HasFormatData(kffMimeTypeTIFF))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_TIFF, BinaryToHex(GetFormatData(kffMimeTypeTIFF)));
    }

    if (HasFormatData(kffMimeTypePNG1))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_PNG, BinaryToHex(GetFormatData(kffMimeTypePNG1)));
    }

    if (HasFormatData(kffMimeTypeJPEG))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_JPEG, BinaryToHex(GetFormatData(kffMimeTypeJPEG)));
    }

    if (HasFormatData(kffMimeTypeBMP))
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_BMP, BinaryToHex(GetFormatData(kffMimeTypeBMP)));
    }

    if (m_rotationAngle != 0)
    {
        CDXMLPutAttribute(sink_arg, kCDXProp_RotationAngle, m_rotationAngle);
    }
}

bool CDXEmbeddedObject::HasFormatData(const string& mimeType) const
{
    return mimeTypeToFormatInfoMap.count(mimeType) != 0;
}

string CDXEmbeddedObject::GetFormatData(const string& mimeType) const
{
    return mimeTypeToFormatInfoMap.at(mimeType).data;
}

uint32_t CDXEmbeddedObject::GetFormatUncompressedSize(const string& mimeType) const
{
    return mimeTypeToFormatInfoMap.at(mimeType).uncompressedSize;
}

void CDXEmbeddedObject::SetFormatData(const string& mimeType, const string& data)
{
    if (!data.empty())
    {
        mimeTypeToFormatInfoMap[mimeType].data = data;
    }
}

void CDXEmbeddedObject::SetFormatUncompressedSize(const string& mimeType, const uint32_t& size)
{
    mimeTypeToFormatInfoMap[mimeType].uncompressedSize = size;
}

std::vector<std::string> CDXEmbeddedObject::GetAvailableFormats() const
{
    vector<string> keys;
    for (auto& it : mimeTypeToFormatInfoMap)
    {
        keys.push_back(it.first);
    }

    return keys;
}

CDXObject* CDXEmbeddedObject::Clone() const
{
	return new CDXEmbeddedObject(*this);
}

std::string CDXEmbeddedObject::XMLObjectName() const
{
	return kCDXML_embeddedobject;
}

void CDXEmbeddedObject::PutTagData(CDXDataSink& sink, const string& dataMimeType) const
{
    auto formatData = GetFormatData(dataMimeType);
    size_t size = formatData.size();
    if (size <= 0xFFFE)
    {
        sink.Put(UINT16(size));
    }
    else
    {
        sink.Put(UINT16(0xFFFF));
        sink.Put(UINT32(size));
    }

    sink.Put((INT8*)formatData.data(), size);
}
