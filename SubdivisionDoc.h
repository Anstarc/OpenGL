// SubdivisionDoc.h : interface of the CSubdivisionDoc class
//


#pragma once

#include "Mesh3D.h"
#include "string.h"
#include "StopWatch.h"

using namespace std;

class CSubdivisionDoc : public CDocument
{
protected: // create from serialization only
	CSubdivisionDoc();
	DECLARE_DYNCREATE(CSubdivisionDoc)

// Attributes
public:
	Mesh3D* m_pmesh;
	Mesh3D* m_pmesh_origin;
// Operations
public:

// Overrides
	public:
	//virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CSubdivisionDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	virtual BOOL OnSaveDocument(LPCTSTR lpszPathName);
	afx_msg void OnImport();
	afx_msg void OnExport();
	afx_msg void OnSubdivisionDoo();
	afx_msg void OnSubdivisionCatmull();
	afx_msg void OnUpdateSubdivisionDoo(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSubdivisionCatmull(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSea(CCmdUI *pCmdUI);////
	afx_msg void OnOperationNew();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	//inline const char * CString2constchar(CString str)
	//{
	//	const size_t strsize = (str.GetLength() + 1) * 2; // 宽字符的长度;
	//	char * pstr = new char[strsize]; //分配空间;
	//	size_t sz = 0;
	//	wcstombs_s(&sz, pstr, strsize, str, _TRUNCATE);
	//	return pstr;
	//	//int n = atoi(const char*)pstr); // 字符串已经由原来的CString 转换成了 const char*
	//}
	//	
	//char * CString2constchar(CString str) {
	//	int x = 0;
	//	string s = "";
	//	while (x < str.GetLength()) {
	//		char c = str.GetAt(x++);
	//		s += c;
	//	}
	//	char * output = (char *)calloc(str.GetLength() + 1, sizeof(char));
	//	memcpy(output, s.c_str(), str.GetLength() + 1);
	//	return output;
	//}



private:
	CStopWatch timer;
	float x=0;
	LPDWORD ThreadId;
	HANDLE hThread;

public:
	static bool loop;
	void LoopDraw();
	static void ThreadFunc(CSubdivisionDoc* p);
};


