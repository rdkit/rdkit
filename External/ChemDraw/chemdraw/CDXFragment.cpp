// CommonCS/LibCommon/Src/CDXFragment.cpp
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

#include <stdexcept>
// ***********************
// ** class CDXFragment **
// ***********************
//
// Specialization of CDXObject for CDXFragment objects

CDXFragment::CDXFragment(CDXObjectID id_arg)
 : CDXGroupOrFragment(kCDXObj_Fragment, id_arg)
{
}

CDXFragment::CDXFragment(const CDXFragment& src) :
    CDXGroupOrFragment(src),	// <- Cloning of fragment's children occurs here.
    isFromGuidedStereo(src.isFromGuidedStereo)
{
    FinishReading();
}

CDXFragment::~CDXFragment()
{
}

CDXObject*	CDXFragment::Clone() const
{
	return new CDXFragment (*this);
}

std::string CDXFragment::XMLObjectName() const
{
	return kCDXML_fragment;
}

void CDXFragment::CountElements(std::map<INT16, INT16> &elementList)
{
	// Loop over all the atoms, and count the various elements
	CDXObjectsRange atoms = ContainedObjects(kCDXObj_Node);
	for (CDXObjectsByTag::const_iterator i = atoms.begin();  i != atoms.end();  ++i)
	{
		CDXNode *p = dynamic_cast<CDXNode *>(GetObject(i));
		if (p == 0) throw std::logic_error("Object with node tag is not a node");
		if (p->m_nodeType == kCDXNodeType_Element)
			++ elementList[p->m_elementNum];

		// Count the elements in the fragment contained in this atom, if any
		CDXObjectsRange frags = p->ContainedObjects(kCDXObj_Fragment);
		for (CDXObjectsByTag::const_iterator i = frags.begin();  i != frags.end();  ++i)
		{
			CDXFragment *f = dynamic_cast<CDXFragment *>(GetObject(i));
			if (f == 0) throw std::logic_error("Object with fragment tag is not a fragment");
			f->CountElements(elementList);
		}

	}
}

void CDXFragment::StoreAttribute(CDXDataSource& src, CDXTag tag, size_t size)
{
    if (tag == kCDXProp_Frag_IsFromGuidedStereo)
    {
        isFromGuidedStereo = (src.GetUINT8() != 0);
    }
    else
    {
        CDXGroupOrFragment::StoreAttribute(src, tag, size);
    }
}

void CDXFragment::XMLStoreAttribute(CDXTag tag, const std::string& value)
{
    if (tag == kCDXProp_Frag_IsFromGuidedStereo)
    {
        isFromGuidedStereo = (value == "yes");
    }
    else
    {
        CDXGroupOrFragment::XMLStoreAttribute(tag, value);
    }
}

void CDXFragment::WriteAttributesTo(CDXDataSink& sink) const
{
    CDXGroupOrFragment::WriteAttributesTo(sink);

    if (isFromGuidedStereo)
    {
        sink.PutAttribute(kCDXProp_Frag_IsFromGuidedStereo, UINT8(isFromGuidedStereo));
    }
}

void CDXFragment::XMLWriteAttributes(XMLDataSink& sink) const
{
    CDXGroupOrFragment::XMLWriteAttributes(sink);

    if (isFromGuidedStereo)
    {
        CDXMLPutAttribute(sink, kCDXProp_Frag_IsFromGuidedStereo, isFromGuidedStereo);
    }    
}

