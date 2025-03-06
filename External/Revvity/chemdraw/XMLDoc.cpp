// CommonCS/LibCommon/Src/XMLDoc.cpp
// Copyright � 2003-2004, CambridgeSoft Corp., All Rights Reserved

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

// BK Not included - #include "cs_stringUtils.h"
#include <sstream>
// BK looks like this is the sourceforge
//  https://sourceforge.net/p/expatpp/code/HEAD/tree/trunk/src_pp/expatpp.h
#include "../expatpp/expatpp-code-r6-trunk/src_pp/expatpp.h"
#include "XMLDoc.h"
#include "CDXObject.h"
#include "CDXUnicode.h"

#define ASSERT assert

#define BUFFSIZE	8192

using namespace std;
class CXmlParser : public expatpp 
{
public:
	CXmlParser(CXmlDoc *doc) : m_Step(0), m_XmlDoc(doc) {}

	// overrideable callbacks
	virtual void startElement(const XML_Char* name, const XML_Char** atts);
	virtual void endElement(const XML_Char*);
	virtual void charData(const XML_Char*, int len);
	virtual void processingInstruction(const XML_Char* target, const XML_Char* data);
	virtual void defaultHandler(const XML_Char*, int len);
	virtual int notStandaloneHandler() { return 1; } // override this to return zero if you want errors for non-standalone documents
	//virtual void unparsedEntityDecl(const XML_Char* entityName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId, const XML_Char* notationName);
	//virtual void notationDecl(const XML_Char* notationName, const XML_Char* base, const XML_Char* systemId, const XML_Char* publicId);
	//virtual void startNamespace(const XML_Char* prefix, const XML_Char* uri);
	//virtual void endNamespace(const XML_Char*);

	inline bool IsFinial() const { return m_Step == 0 && m_Parents.empty(); }

private:
	int					m_Step;
	CXmlDoc *			m_XmlDoc;
	stack<CXmlNode *>	m_Parents;
};

//////////////////////////////////////////////////////////////////////////////////////
// class CXmlNode

CXmlNode::CXmlNode(const string& name)
{
	string::size_type colonPos = name.find(':');
	if (colonPos == string::npos)
	{
		m_Name = name;
	}
	else
	{
		m_Namespace = name.substr(0, colonPos);
		m_Name = name.substr(colonPos + 1);
	}
}

CXmlNode& CXmlNode::operator =(const CXmlNode& node)
{
	m_Name			= node.m_Name;
	m_Data			= node.m_Data;
	m_Comments		= node.m_Comments;
	m_Attributes	= node.m_Attributes;

	for (int i = 0; i < node.GetChildCount(); ++i)
	{
		const CXmlNode *child = node.GetNthChild(i);
		*CreateChild(string()) = *child;
	}

	return *this;
}

// copied from XMLPut()
// not escape \r and \n
void CXmlNode::XMLPut1(XMLDataSink &sink, const string &s, const string& which1)
{
	string which = which1;
	if (which.empty())
		which = "&<>\'\"\n\r";

	string::size_type	curPos = 0,
						nextPos;
	while ((nextPos = s.find_first_of(which, curPos)) != string::npos)
	{
		if (nextPos > curPos)
			sink.os.write(s.data() + curPos, nextPos - curPos);
		switch (s[nextPos])
		{
		case '&':
			sink.os << "&amp;";
			break;
		case '<':
			sink.os << "&lt;";
			break;
		case '>':
			sink.os << "&gt;";
			break;
		case '\'':
			sink.os << "&apos;";
			break;
		case '\"':
			sink.os << "&quot;";
			break;
		case '\n':
			sink.os << "&#010;";
			break;
		case '\r':
			sink.os << "&#013;";
			break;
		default:
			ASSERT(false); // should never get here as the list of chars should all be in the switch
			break;
		}
		curPos = nextPos + 1;
	}

	if (s.size() > curPos)
		sink.os.write(s.data() + curPos, s.size() - curPos);
}

string CXmlNode::Encode(const string& s, bool onlyMust)
{
	std::ostringstream os;
	XMLDataSink sink(os);

	XMLPut1(sink, s, onlyMust ? "&<>" : "&<>\'\"\n\r");
	return os.str();
}

string CXmlNode::Decode(const string& s)
{
	string out;
	XMLTranslateEscapedString(UnicodeToString(s), out);
	return out;
}

size_t CXmlNode::SelectChildNodes(const string& path, vector<const CXmlNode *>& nodes) const
{
	vector<string> h;
	size_t n = ParsePath(path, h);
	if (n == 0)
		return 0;
	
	string tag = h.back();
	h.pop_back();

	const CXmlNode *directParent = h.empty() ? this : SelectFirstChild(MakePath(h));
	if (directParent == NULL)
		return 0;
	return directParent->GetChildNodes(tag, nodes);
}

size_t CXmlNode::GetChildNodes(const string& tag, vector<const CXmlNode *>& nodes) const
{
	nodes.clear();

	for (int i = 0; i < m_Nodes.size(); ++i)
	{
		if (m_Nodes[i]->GetName() == tag)
			nodes.push_back(m_Nodes[i]);
	}

	return nodes.size();
}

string CXmlNode::MakePath(const vector<string>& heritance)
{
	string ret;

	if (heritance.empty())
		return ret;

	if (heritance.size() == 1)
		return heritance[0];

	ret = heritance[0];
	for (int i = 1; i < heritance.size(); ++i)
	{
		ret += '/';
		ret += heritance[i];
	}

	return ret;
}

size_t CXmlNode::ParsePath(const string& path, vector<string>& heritance)
{
	heritance.clear();

	string s;
	for (int i = 0; i < path.length(); ++i)
	{
		if (path[i] == '\\' || path[i] == '/')
		{
			if (s.empty() && i > 0)
				return 0;
			
			heritance.push_back(s);
			s.clear();
		}
		else
		{
			s += path[i];
		}
	}

	if (!s.empty())
		heritance.push_back(s);

	return heritance.size();
}

const CXmlNode *CXmlNode::GetFirstChild(const string& tagname, const char* attname, const char* attvalue) const
{
	string att = attname == NULL ? "" : attname;
	for (int i = 0; i < m_Nodes.size(); ++i)
	{
		if (m_Nodes[i]->GetName() == tagname &&
			(attname == NULL || m_Nodes[i]->GetAttribute(att) == attvalue))
			return m_Nodes[i];
	}

	return NULL;
}

const CXmlNode *CXmlNode::SelectFirstChild(const string& path) const
{
	const CXmlNode *ret = NULL;

	vector<string> h;
	if (ParsePath(path, h) > 0)
	{
		ret = GetFirstChild(h[0]);
		if (ret != NULL)
		{
			for (int i = 1; i < h.size(); ++i)
			{
				ret = ret->GetFirstChild(h[i]);
				if (ret == NULL)
					break;
			}
		}
	}

	return ret;
}

const string & CXmlNode::GetAttribute(const string& name)	const
{
	if (!m_Attributes.empty())
	{
		map<string, string>::const_iterator p = m_Attributes.find(name);
		if (p != m_Attributes.end())
			return p->second;
	}

	static string emptyString;
	return emptyString;
}

bool CXmlNode::GetNthAttribute(int i, string& name, string& value) const
{
	if (i >= 0 && i < m_Attributes.size())
	{
		int N = 0;
		for (map<string, string>::const_iterator p = m_Attributes.begin(); p != m_Attributes.end(); ++p)
		{
			if (N == i)
			{
				name = p->first;
				value = p->second;
				break;
			}

			++N;
		}
		return true;
	}

	return false;
}

bool CXmlNode::Store(ostream& os) const
{
	if (GetName().empty())
		return false;

	// tag name
	os << '<';
	if (!GetNamespace().empty())
		os << GetNamespace() << ':';
	os << GetName();

	// attributes
	for (int i = 0; i < GetAttributesCount(); ++i)
	{
		string name, value;
		if (!GetNthAttribute(i, name, value))
			return false;

		os << ' ' << name << "=\"" << StringToUnicode("", Encode(value, false)) << "\"";
	}

	bool nochild = GetChildCount() == 0 && GetCommentsCount() == 0;
	if (nochild && GetData().empty())
	{
		// no data, no children, no comments
		os << "/>" << GetTextEOL();
	}
	else
	{
		os << '>';

		if (!nochild)	// no child nodes, no comments
			os << GetTextEOL();

		// comments
		for (int i = 0; i < GetCommentsCount(); ++i)
			os << "<!-- " << GetNthComment(i) << " -->" << GetTextEOL();

		// data
		if (!GetData().empty())
		{
			os << StringToUnicode("", Encode(GetData(), true));
			if (!nochild)	// no child nodes, no comments
				os << GetTextEOL();
		}
		
		// child nodes
		for (int i = 0; i < GetChildCount(); ++i)
			GetNthChild(i)->Store(os);

		os << "</";
		if (!GetNamespace().empty())
			os << GetNamespace() << ":";
		os << GetName() << '>' << GetTextEOL();
	}

	return true;
}

bool CXmlNode::GetXml(string& xml) const
{
	std::ostringstream os;

	bool ret = Store(os);
	if (ret)
	{
	//	os.freeze();
	//	xml = string(os.str(), os.pcount());
		xml = os.str();
	}

	return ret;
}

CXmlNode *CXmlNode::CreateChild(const string& name)
{
	CXmlNode *node = NULL;
	if (IsValidTag(name))
	{
		node = new CXmlNode(name);
		m_Nodes.push_back(node);
	}

	return node;
}

bool CXmlNode::IsValidTag(const string& tag)
{
	return true;
}

CXmlNode *CXmlNode::CreateChild(const string& name, const string& data)
{
	CXmlNode *node = CreateChild(name);
	if (node != NULL)
		node->SetData(data);

	return node;
}

//////////////////////////////////////////////////////////////////////////////////////
// class CXmlParser

void CXmlParser::startElement(const XML_Char* name, const XML_Char** attr)
{
	CXmlNode *parent = m_Parents.empty() ? NULL : m_Parents.top();
	CXmlNode *node = NULL;
	if (parent != NULL)
		node = parent->CreateChild(name);
	else
		node = m_XmlDoc->CreateRootNode(name);

	for (int i = 0; attr[i]; i += 2) 
		node->AddAttribute(attr[i], UnicodeToString(attr[i + 1]));

	m_Parents.push(node);
	++m_Step;
}

void CXmlParser::endElement(const XML_Char*)
{
	m_Parents.pop();
	--m_Step;
}

void CXmlParser::charData(const XML_Char *data, int len)
{
	// tag data
	CXmlNode *parent = m_Parents.empty() ? NULL : m_Parents.top();
//	if (parent != NULL && !(len == 1 && (data[0] == '\r' || data[0] =='\n')))
	if (parent != NULL)
		parent->AppendData(UnicodeToString(string(data, len)));
}

void CXmlParser::processingInstruction(const XML_Char* target, const XML_Char* data)
{
	// document instruction
	if (m_XmlDoc != NULL)
		m_XmlDoc->AddInstruction(target, data);
}

void CXmlParser::defaultHandler(const XML_Char *data, int len)
{
	if (isprint(data[0]))
	{
		string s(data, len);
		if (s.length() > 4 && s.substr(0, 2) == "<?" && s.substr(s.length() - 2, 2) == "?>")
		{
			// none comments
			if (m_XmlDoc != NULL)
				m_XmlDoc->AddNonComm(s);
		}
		else if (s.length() > 7 && s.substr(0, 4) == "<!--" && s.substr(s.length() - 3, 3) == "-->")
		{
			s = s.substr(4, s.length() - 7);

			CXmlNode *parent = m_Parents.empty() ? NULL : m_Parents.top();
			if (parent != NULL)
				parent->AddComment(s);		// node comments
			else if (m_XmlDoc != NULL)
				m_XmlDoc->AddComment(s);	// document comments
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////
// class CXmlDoc

CXmlDoc& CXmlDoc::operator =(const CXmlDoc& doc)
{
	m_Instructions	= doc.m_Instructions;
	m_Comments		= doc.m_Comments;
	m_NonComms		= doc.m_NonComms;
	m_ErrorStr		= doc.m_ErrorStr;

	*CreateRootNode(string()) = *doc.GetRootNode();

	return *this;
}

bool CXmlDoc::Open(istream& is)
{
	CXmlParser parser(this);

	char buf[BUFFSIZE + 1];
	bool succeed = true;
	while (!is.eof() && is.good())
	{
		is.read(buf, BUFFSIZE);
		streamsize len = is.gcount();
		if (len == 0)
			break;
		buf[len] = 0;

		if (!parser.XML_Parse(buf, len, is.eof()))
		{
			m_ErrorStr = "Error: [line:" + cs::NumToStr(parser.XML_GetCurrentLineNumber()) + "] ";
			m_ErrorStr += XML_ErrorString(parser.XML_GetErrorCode());

			succeed = false;
			break;
		}
	}

	succeed = succeed && (parser.IsFinial() && GetRootNode() != NULL);
	if (!succeed)
		Clear();

	return succeed;
}

bool CXmlDoc::Store(ostream& os) const
{
	// instructions must go first
	for (int i = 0; i < GetInstructionsCount(); ++i)
	{
		string target, data;
		if (!GetNthInstruction(i, target, data))
			return false;

		os << "<?" << target << ' ' << data << "?>" << GetTextEOL();
	}

	// none comments
	for (int i = 0; i < GetNonCommsCount(); ++i)
		os << GetNthNonComm(i) << GetTextEOL();

	// comments
	for (int i = 0; i < GetCommentsCount(); ++i)
		os << "<!-- " << GetNthComment(i) << " -->" << GetTextEOL();

	// nodes
	bool ret = false;
	const CXmlNode *root = GetRootNode();
	if (root != NULL)
		ret = root->Store(os);

	return ret;
}

bool CXmlDoc::GetNthInstruction(int i, string& target, string& data) const
{
	if (i >= 0 && i < m_Instructions.size())
	{
		int N = 0;
		for (map<string, string>::const_iterator p = m_Instructions.begin(); p != m_Instructions.end(); ++p)
		{
			if (N == i)
			{
				target = p->first;
				data = p->second;
				break;
			}

			++N;
		}
		return true;
	}

	return false;
}

void CXmlDoc::Clear()
{ 
	delete m_Root; 
	m_Root = NULL; 
	
	m_Instructions.clear(); 
	m_Comments.clear();
}
