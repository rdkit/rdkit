// CommonCS/LibCommon/Hdr/XmlDoc.h
// Copyright © 2003-2004, CambridgeSoft Corp., All Rights Reserved

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
// turn off warnings of size_t to int.
#pragma warning (disable : 4267)
#endif

//#include "cs_assert.h"
#include <algorithm>
#include <istream>
#include <functional>
#include <map>
#include <ostream>
#include <stack>
#include <string>
#include <vector>

using std::ostringstream;
using std::istringstream;
using std::stringstream;
using std::ostream;
using std::ofstream;
using std::streamsize;
using std::exception;
using std::nothrow_t;
using std::bad_alloc;
using std::swap;
using std::make_pair;
using std::fill_n;
using std::ios;
using std::multimap;
using std::streambuf;
using std::endl;
using std::hex;
using std::vector;
using std::ostream;
using std::min;
using std::max;
using std::string;
using std::map;
using std::istream;

class XMLDataSink;

class CXmlNode
{
public:
	CXmlNode(const string& name);
	CXmlNode(const CXmlNode& node) { *this = node; }
	virtual ~CXmlNode() { for (int i = 0; i < (int)m_Nodes.size(); ++i) delete m_Nodes[i]; }

	CXmlNode& operator =(const CXmlNode& node);

	const string&	GetNamespace() const			{ return m_Namespace; }
	const string&	GetName() const					{ return m_Name; }

	const string&	GetData() const					{ return m_Data; }
	void			SetData(const string& data)		{ m_Data = data; }
	void			AppendData(const string& data)	{ m_Data += data; }

	void			AddComment(const string &s)		{ m_Comments.push_back(s); }
	size_t			GetCommentsCount() const		{ return m_Comments.size(); }
	string			GetNthComment(int i) const		{ return (i >= 0 && i < (int)m_Comments.size()) ? m_Comments[i] : ""; }

	size_t			GetAttributesCount() const								{ return m_Attributes.size(); }
#ifdef _AFX
#if defined UNICODE || defined _UNICODE
	void			AddAttribute(LPCTSTR name, LPCTSTR value)				{ AddAttribute(string(CT2A(name)), string(CT2A(value))); }
#endif // defined UNICODE || defined _UNICODE
	void			AddAttribute(LPCSTR name, LPCTSTR value)				{ AddAttribute(string(name), string(CT2A(value))); }
#endif // _AFX
	void			AddAttribute(const string& name, const string& value)	{ m_Attributes.insert(make_pair(name, value)); }
	const string &	GetAttribute(const string& name)	const;
	bool			GetNthAttribute(int i, string& name, string& value) const;

	size_t				GetChildCount() const		{ return m_Nodes.size(); }
	const CXmlNode *	GetNthChild(int i) const	{ return (i >= 0 && i < (int)m_Nodes.size()) ? m_Nodes[i] : NULL; }
	size_t				SelectChildNodes(const string& path, vector<const CXmlNode *>& nodes) const;
	size_t				GetChildNodes(const string& tag, vector<const CXmlNode *>& nodes) const;
	const CXmlNode *	SelectFirstChild(const string& path) const;
	const CXmlNode *	GetFirstChild(const string& tagname, const char* attname = NULL, const char* attvalue = NULL) const;
	CXmlNode *			CreateChild(const string& name);
	CXmlNode *			CreateChild(const string& name, const string& data);

	bool Store(ostream& os) const;
	bool GetXml(string& xml) const;

	static void		XMLPut1(XMLDataSink &sink, const string &s, const string& which = "");
	static string	Encode(const string& s, bool onlyMust);
	static string	Decode(const string& s);
	static size_t	ParsePath(const string& path, vector<string>& heritance);
	static string	MakePath(const vector<string>& heritance);
	static bool		IsValidTag(const string& tag);

private:
	string m_Namespace;
	string m_Name;
	string m_Data;

	vector<string>		m_Comments;
	vector<CXmlNode *>	m_Nodes;
	map<string, string> m_Attributes;
};

class CXmlDoc
{
public:
	CXmlDoc() : m_Root(NULL) {}
	CXmlDoc(const CXmlDoc& doc) { *this = doc; }
	virtual ~CXmlDoc() { Clear(); }

	CXmlDoc& operator =(const CXmlDoc& doc);

	bool	Open(istream& is);
	bool	Store(ostream& os) const;
	string	GetErrorStr() const { return m_ErrorStr; }

	void			Clear();
	CXmlNode *		CreateRootNode(const string& name)	{ delete m_Root; return (m_Root = new CXmlNode(name)); }
	const CXmlNode *GetRootNode() const					{ return m_Root; }

	string	GetInstruction(const string& target)		{ return (!m_Instructions.empty() && m_Instructions.find(target) != m_Instructions.end()) ? m_Instructions[target] : ""; }
	size_t	GetInstructionsCount() const				{ return m_Instructions.size(); }
	bool	GetNthInstruction(int i, string& target, string& data) const;
	void	AddInstruction(const string& target, const string& value)	{ m_Instructions.insert(std::make_pair(target, value)); }

	size_t	GetCommentsCount() const		{ return m_Comments.size(); }
	void	AddComment(const string &s)		{ m_Comments.push_back(s); }
	string	GetNthComment(int i) const		{ return (i >= 0 && i < (int)m_Comments.size()) ? m_Comments[i] : ""; }

	size_t	GetNonCommsCount() const		{ return m_NonComms.size(); }
	void	AddNonComm(const string &s)		{ m_NonComms.push_back(s); }
	string	GetNthNonComm(int i) const		{ return (i >= 0 && i < (int)m_NonComms.size()) ? m_NonComms[i] : ""; }

protected:
	CXmlNode *			m_Root;
	map<string, string> m_Instructions;
	vector<string>		m_Comments;
	vector<string>		m_NonComms;
	string				m_ErrorStr;
};
