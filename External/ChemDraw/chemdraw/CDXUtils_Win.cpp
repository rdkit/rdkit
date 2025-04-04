/*
File:		CDXUtils_Win.cpp
Purpose:	Windows-only utility functions for manipulating CDX.
Copyright:	(c) 2002-2010, CambridgeSoft Corp., All Rights Reserved
*/

#include "CDXUtils.h"
#include "CDXUtils_Win.h"

/*
+===========================================================================================+
|										WINDOWS-ONLY										|
+===========================================================================================+
*/
#ifdef _AFX
/*
+===========================================================================================+
| CDXSharedFileSink		Helper class to write CDX to a CSharedFile.							|
+===========================================================================================+
*/
class CDXSharedFileSink : public CDXDataSink
{
public:
	CDXSharedFileSink() : m_file (GMEM_DDESHARE | GMEM_MOVEABLE, 1000)	{}
	virtual void Put (const INT8 *data, size_t n )	{ m_file.Write (data, n); }

	HGLOBAL		ReleaseMemory()	{ return m_file.Detach(); }

protected:
	CSharedFile	m_file;
};

/*
+===========================================================================================+
| CDXSharedFileSource	Helper class to read CDX from a CSharedFile.						|
+===========================================================================================+
*/
class CDXSharedFileSource : public CDXDataSource
{
public:
	CDXSharedFileSource (yours HANDLE hdl)	{ 	m_file.SetHandle (hdl);  m_fAttached = true;}
	~CDXSharedFileSource()					{ 	Detach(); }
	void		Detach()	{	if (!m_fAttached) return;
								m_fAttached = false;
								HANDLE	hand = m_file.Detach();	// returns same handle as initially provided.
								::GlobalUnlock (hand);	// CSharedFile::Detach() doesn't do this
							}

	virtual void GetBytes (char *p, size_t n)
	{
		 try
		 {
			 if (m_file.Read (p, n) == n)
				return;
		 }
		 catch (CFileException*) {}
		 throw runtime_error ("Error reading CDX");
	}
	virtual void SkipBytes (size_t n)
	{
		 try
		 {
			 m_file.Seek (n, CFile::current);
		 }
		 catch (CFileException*) {}
		 throw runtime_error ("Error reading CDX");
	}

protected:
	CSharedFile			m_file;
	bool				m_fAttached;
};
//---------------------------------------------------------------------

/*
+===========================================================================================+
| CdxDoc_To_CdxHandle	Given a CDX doc, return a blob of CDX binary.						|
+===========================================================================================+
*/
yours HGLOBAL CdxDoc_To_CdxHandle (const CDXDocument *pCdxDoc)
{
	if (!pCdxDoc)
		return NULL;

	CDXSharedFileSink	sink;
	try
	{
		CDXWriteDocToStorage (pCdxDoc, sink);
	}
	catch (exception)
	{
		return NULL;
	}
	HGLOBAL		hdl = sink.ReleaseMemory();
	return hdl;
}

/*
+-------------------------------------------------------------------------------------------+
| HandleToStream_local	Move data between a global object and a std stream.					|
+-------------------------------------------------------------------------------------------+
*/
inline void HandleToStream_local (HANDLE handle, ostream& os, bool includeLen = false)
{
	if (!handle)
		return;
	UINT8*	pObj = (LPBYTE) GlobalLock (handle);
	DWORD	len = GlobalSize (handle);
	if (includeLen)
	{
		len = *(DWORD*)pObj;
		pObj += sizeof (DWORD);
	}
	os.write ((const char*)pObj, len);
	GlobalUnlock (handle);
}

/*
+===========================================================================================+
| CdxHandle_To_CdxDoc	Given a blob of CDX binary, return CDX doc.							|
|						The fre is only defined if there is an error, but may return NULL	|
|						even if no error occurred.											|
+===========================================================================================+
*/
bool CdxHandle_To_CdxDoc (yours HANDLE hdl, yours CDXDocument*& pDoc, string *pErrMsg)
{
	pDoc = NULL;
	if (!hdl)
		return true;
	const int	numBytes = GlobalSize (hdl);
	if (!numBytes)
		return true;

//	HDBG ( _CrtCheckMemory(); )
	{
		CDXSharedFileSource	source (hdl);
		string	errMsg;
		pDoc = CDXReadDocFromStorage (source, errMsg);
		if (!errMsg.empty())
		{
			ASSERT (pDoc == NULL);
			if (pErrMsg != NULL)
				*pErrMsg = errMsg;
			return false;
		}
	}

	if (0)	// This passage allows you to store the blob as a disk file
	{
//		ofstream	ofs ("c:/Cdx.cdx", ios::out | ios::binary);
//		HandleToStream_local (hdl, ofs);
	}

//	HDBG ( _CrtCheckMemory(); )
//	DBG ( if (0) )
//	DBG (	if (pDoc) { ccOut << "CdxHandle_To_CdxDoc:\n"; CDXMLDump (ccOut.GetStream(), pDoc); } )
//	DBG ( if (0) )
//	DBG (	if (pDoc) { ccOut << "CdxHandle_To_CdxDoc:\n"; ccDivertOutput divo (ccOut, "c:/cdx.xml", ccOutput::kOverwriteFileIfExists);  pDoc->XMLWrite (XMLDataSink (ccOut.GetStream())); } )
	return true;
}

#endif // _AFX
