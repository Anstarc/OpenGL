// SubdivisionView.h : interface of the CSubdivisionView class
//


#pragma once




#include "Lib/arcball/Arcball.h"
#include "Lib/arcball/Camera.h"
#include "Lib/arcball/Viewport.h"
#define MAX_CHARS 256

class CSubdivisionView : public CView
{
protected: // create from serialization only
	CSubdivisionView();
	DECLARE_DYNCREATE(CSubdivisionView)

// Attributes
public:
	CSubdivisionDoc* GetDocument() const;

private:
	// mouse
	bool m_LeftButtonDown;
	bool m_RightButtonDown;
	CPoint m_LeftDownPos;
	CPoint m_RightDownPos;
	bool m_Moving;

	// arcball
	CCamera   m_Camera;
	CViewport m_Viewport;
	CArcball  m_Arcball;

	// OpenGL
	HGLRC m_hGLContext;
	int m_GLPixelIndex;
	int m_PolygonMode;
	bool m_SmoothShading;
	bool m_FirstView;
	bool m_showedge;

	float m_BackColor[3];
	float m_MeshColor[3];
	float m_EdgeColor[3];
	float m_VertexColor[3];


	UINT g_FontListID;
	HFONT hOldFont;
	int FONT_HEIGHT;

// Operations
public:
	// arcball
	CArcball *arcball() { return &m_Arcball; }

	// OpenGL specific
	BOOL SetWindowPixelFormat(HDC hDC);
	BOOL CreateViewGLContext(HDC hDC);

	// camera
	void InitCamera();
	void ViewAll(bool check_first = true);

	// mouse 
	void HandleMouseButton(int x, int y);


	//For Font Display
	UINT CreateOpenGLFont(LPWSTR strFontName, int height);
	void PositionText( int x, int y );
	virtual void glDrawText(int x, int y, const char *strString, ...);
	void DestroyFont();
	int GetFontHeight() { return FONT_HEIGHT; }
// Overrides
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	void ChangeMaterial(CString &string, bool update);
// Implementation
public:
	virtual ~CSubdivisionView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	virtual void OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnDestroy();
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnShademodeSmooth();
	afx_msg void OnUpdateShademodeSmooth(CCmdUI *pCmdUI);
	afx_msg void OnViewViewall();
	afx_msg void OnShadeWireframe();
	afx_msg void OnUpdateShadeWireframe(CCmdUI *pCmdUI);
	afx_msg void OnViewShowedge();
	afx_msg void OnUpdateViewShowedge(CCmdUI *pCmdUI);
	afx_msg void OnMaterialSilver();
	afx_msg void OnMaterialGold();
	afx_msg void OnMaterialPearl();
	afx_msg void OnMaterialBrass();
	afx_msg void OnMaterialBlackplastic();
	afx_msg void OnMaterialRuby();
	afx_msg void OnMaterialJade();
	afx_msg void OnMaterialLightblue();
	afx_msg void OnMaterialEmerald();
	afx_msg void OnMaterialPolishedsilver();
	afx_msg void OnMaterialChrome();
	afx_msg void OnMaterialCopper();
	afx_msg void OnMaterialPolishedgold();
	afx_msg void OnMaterialPewter();
	afx_msg void OnMaterialObsidian();
	afx_msg void OnMaterialPolishedbronze();
	afx_msg void OnMaterialPolishedcopper();
	afx_msg void OnMaterialTurquoise();
	afx_msg void OnSea();////
	bool LoadTexture(LPTSTR szFileName, GLuint & texid);
	void CopyTo(Mesh3D *m_pmesh1, Mesh3D *m_pmesh2);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	Mesh3D* OnCreateSea(Mesh3D* m_pmesh, float a);

	GLuint textureID = 0;
	float time = 0;
};

#ifndef _DEBUG  // debug version in SubdivisionView.cpp
inline CSubdivisionDoc* CSubdivisionView::GetDocument() const
   { return reinterpret_cast<CSubdivisionDoc*>(m_pDocument); }
#endif

