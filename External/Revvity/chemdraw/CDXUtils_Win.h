// CommonCS/LibCommon/Hdr/CDXUtils_Win.h
// Copyright © 1997-2010, CambridgeSoft Corp., All Rights Reserved

#pragma once

// The following limitation could be weakened to "_WINDOWS" if our CDXSharedFileSource did
// not use a CSharedFile::.  -heh 5/28/09
#ifdef _AFX
using std::runtime_error;
class CDXDocument;

yours HGLOBAL		CdxDoc_To_CdxHandle	(const CDXDocument *pCdxDoc);
bool				CdxHandle_To_CdxDoc	(yours HANDLE hdl, yours CDXDocument*& pDoc, string *pErrMsg);
#endif // _AFX

