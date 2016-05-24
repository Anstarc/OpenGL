// SubdivisionView.cpp : implementation of the CSubdivisionView class
//


#include "stdafx.h"
#include "Subdivision.h"

#include "SubdivisionDoc.h"
#include "SubdivisionView.h"
#include ".\subdivisionview.h"

#include "CreateSea.h"
//#include <gl/glut.h>
//#include <gl/gl.h>
//#include "glew.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


#define GL_SEPARATE_SPECULAR_COLOR 0x81FA
#define GL_LIGHT_MODEL_COLOR_CONTROL 0x81F8


#define PI 3.1415926536
#define zDepth -10.0
#define yDepth -0.2
#define xWith 10.0
#define STEP 0.5

// CSubdivisionView

IMPLEMENT_DYNCREATE(CSubdivisionView, CView)

BEGIN_MESSAGE_MAP(CSubdivisionView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ACTIVATE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_COMMAND(ID_SHADEMODE_SMOOTH, OnShademodeSmooth)
	ON_UPDATE_COMMAND_UI(ID_SHADEMODE_SMOOTH, OnUpdateShademodeSmooth)
	ON_COMMAND(ID_VIEW_VIEWALL, OnViewViewall)
	ON_COMMAND(ID_SHADE_WIREFRAME, OnShadeWireframe)
	ON_UPDATE_COMMAND_UI(ID_SHADE_WIREFRAME, OnUpdateShadeWireframe)
	ON_COMMAND(ID_VIEW_SHOWEDGE, OnViewShowedge)
	ON_UPDATE_COMMAND_UI(ID_VIEW_SHOWEDGE, OnUpdateViewShowedge)
	ON_COMMAND(ID_MATERIAL_SILVER, OnMaterialSilver)
	ON_COMMAND(ID_MATERIAL_GOLD, OnMaterialGold)
	ON_COMMAND(ID_MATERIAL_PEARL, OnMaterialPearl)
	ON_COMMAND(ID_MATERIAL_BRASS, OnMaterialBrass)
	ON_COMMAND(ID_MATERIAL_BLACKPLASTIC, OnMaterialBlackplastic)
	ON_COMMAND(ID_MATERIAL_RUBY, OnMaterialRuby)
	ON_COMMAND(ID_MATERIAL_JADE, OnMaterialJade)
	ON_COMMAND(ID_MATERIAL_LIGHTBLUE, OnMaterialLightblue)
	ON_COMMAND(ID_MATERIAL_EMERALD, OnMaterialEmerald)
	ON_COMMAND(ID_MATERIAL_POLISHEDSILVER, OnMaterialPolishedsilver)
	ON_COMMAND(ID_MATERIAL_CHROME, OnMaterialChrome)
	ON_COMMAND(ID_MATERIAL_COPPER, OnMaterialCopper)
	ON_COMMAND(ID_MATERIAL_POLISHEDGOLD, OnMaterialPolishedgold)
	ON_COMMAND(ID_MATERIAL_PEWTER, OnMaterialPewter)
	ON_COMMAND(ID_MATERIAL_OBSIDIAN, OnMaterialObsidian)
	ON_COMMAND(ID_MATERIAL_POLISHEDBRONZE, OnMaterialPolishedbronze)
	ON_COMMAND(ID_MATERIAL_POLISHEDCOPPER, OnMaterialPolishedcopper)
	ON_COMMAND(ID_MATERIAL_TURQUOISE, OnMaterialTurquoise)
	ON_COMMAND(ID_SEA, OnSea)////
	ON_WM_TIMER()
END_MESSAGE_MAP()

// CSubdivisionView construction/destruction

CSubdivisionView::CSubdivisionView()
{
	// TODO: add construction code here
	// OpenGL
	m_hGLContext = NULL;
	m_GLPixelIndex = 0;

	// rendering options
	m_FirstView = true;
	m_SmoothShading = false;
	m_PolygonMode = GL_FILL;
	m_showedge = true;
	// other options
	m_BackColor[0] = m_BackColor[1] = m_BackColor[2] = 0.0f;
	m_MeshColor[0] = m_MeshColor[1] = m_MeshColor[2] = 1.0f;
	m_EdgeColor[0] = m_EdgeColor[1] = m_EdgeColor[2] = 0.0f;
	m_VertexColor[0] = m_VertexColor[1] = m_VertexColor[2] = 0.0;

	// mouse
	m_LeftButtonDown = false;
	m_RightButtonDown = false;
	m_Moving = false;

	FONT_HEIGHT = 20;

}

CSubdivisionView::~CSubdivisionView()
{
}

BOOL CSubdivisionView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CSubdivisionView drawing

void CSubdivisionView::OnDraw(CDC* /*pDC*/)
{
	//CSubdivisionDoc* pDoc = GetDocument();
	//ASSERT_VALID(pDoc);
	//if (!pDoc)
	//	return;

	// TODO: add draw code for native data here
}

//////////////////////////////////////////////
// OPENGL
//////////////////////////////////////////////

//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CSubdivisionView::SetWindowPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pixelDesc;

	pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion = 1;

	pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
		PFD_DOUBLEBUFFER | PFD_STEREO_DONTCARE;

	pixelDesc.iPixelType = PFD_TYPE_RGBA;
	pixelDesc.cColorBits = 32;
	pixelDesc.cRedBits = 8;
	pixelDesc.cRedShift = 16;
	pixelDesc.cGreenBits = 8;
	pixelDesc.cGreenShift = 8;
	pixelDesc.cBlueBits = 8;
	pixelDesc.cBlueShift = 0;
	pixelDesc.cAlphaBits = 0;
	pixelDesc.cAlphaShift = 0;
	pixelDesc.cAccumBits = 64;
	pixelDesc.cAccumRedBits = 16;
	pixelDesc.cAccumGreenBits = 16;
	pixelDesc.cAccumBlueBits = 16;
	pixelDesc.cAccumAlphaBits = 0;
	pixelDesc.cDepthBits = 32;
	pixelDesc.cStencilBits = 8;
	pixelDesc.cAuxBuffers = 0;
	pixelDesc.iLayerType = PFD_MAIN_PLANE;
	pixelDesc.bReserved = 0;
	pixelDesc.dwLayerMask = 0;
	pixelDesc.dwVisibleMask = 0;
	pixelDesc.dwDamageMask = 0;

	m_GLPixelIndex = ChoosePixelFormat(hDC,&pixelDesc);
	if(m_GLPixelIndex == 0) // Choose default
	{
		m_GLPixelIndex = 1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,
			sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;
	}

	if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
		return FALSE;

	return TRUE;
}

//********************************************
// CreateViewGLContext
// Create an OpenGL rendering context
//********************************************
BOOL CSubdivisionView::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);

	if(m_hGLContext==NULL)
		return FALSE;

	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;

	g_FontListID = CreateOpenGLFont(_T("Arial"), FONT_HEIGHT);

	return TRUE;	
}

//********************************************
// Create an OpenGL framework
//********************************************
int CSubdivisionView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if(CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);

	if(SetWindowPixelFormat(hDC)==FALSE)
		return -1;

	if(CreateViewGLContext(hDC)==FALSE)
		return -1;

	// Default mode
	glEnable(GL_DEPTH_TEST);
	glClearColor(m_BackColor[0],
		m_BackColor[1],
		m_BackColor[2],
		1.0f);

	// lighting
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	GLfloat ambientLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };//白色主光源环境光  
	GLfloat diffuseLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };//白色主光源漫反射  
	//GLfloat specularLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };//白色主光源镜面光  
	float lpos[4] = { 0.0f, -20.0f, -20.0f, 1.0f };
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);   //设置环境光源  
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);   //设置漫反射光源  
	//glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight); //设置镜面光源  
	glLightfv(GL_LIGHT1,GL_POSITION,lpos);
	//glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);
	glEnable(GL_LIGHT1);
	glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

	// back material
	float MatAmbientBack[]  = {0.0f, 1.0f, 0.0f, 1.0f};
	//glMaterialfv(GL_BACK,GL_AMBIENT,MatAmbientBack);
	ChangeMaterial(CString("Gold"),false);

	// init camera
	InitCamera();
	return 0;
}


// resize client
void CSubdivisionView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set OpenGL viewport and perspective
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	m_Viewport.SetSize(cx,cy);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);
	::ReleaseDC(hWnd,hDC);
}

void CSubdivisionView::OnActivateView(BOOL bActivate, 
									  CView* pActivateView, 
									  CView* pDeactiveView) 
{
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);
	::ReleaseDC(hWnd,hDC);
	CView::OnActivateView(bActivate, 
		pActivateView, 
		pDeactiveView);
}

// avoid flickering
BOOL CSubdivisionView::OnEraseBkgnd(CDC*)
{
	return TRUE;
}

// destroy rendering context
void CSubdivisionView::OnDestroy()
{
	if(wglGetCurrentContext() != NULL)
		wglMakeCurrent(NULL,NULL);

	if(m_hGLContext != NULL)
	{
		wglDeleteContext(m_hGLContext);
		m_hGLContext = NULL;
	}
}

// handle mouse button
void CSubdivisionView::HandleMouseButton(int x, 
										 int y)
{
	CVector3d vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,
		m_Camera,m_Viewport);
	m_Arcball.EndDrag(vec);
	m_Arcball.SetMode(m_LeftButtonDown+2*m_RightButtonDown);
	vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,m_Camera,m_Viewport);
	m_Arcball.BeginDrag(vec);
}

// mouse move
void CSubdivisionView::OnMouseMove(UINT nFlags, CPoint point)
{
	if(m_LeftButtonDown || m_RightButtonDown)
	{
		CVector3d vec = m_Arcball.Intersect(point.x,m_Viewport.yRes()-point.y,
			m_Camera,m_Viewport);
		m_Arcball.Motion(vec);
		m_Moving = true;
		InvalidateRect(NULL,FALSE);
	}
	else
		m_Moving = false;
}

void CSubdivisionView::OnLButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_LeftButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonDown(nFlags, point);
}

void CSubdivisionView::OnLButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_LeftButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnRButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_RightButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonDown(nFlags, point);
}

void CSubdivisionView::OnRButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_RightButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}

//////////////////////////////////////////////
// CAMERA
//////////////////////////////////////////////

// view all the current object
void CSubdivisionView::ViewAll(bool check_first)
{
	// make it only once
	if(!m_FirstView && check_first)
		return;

	CSubdivisionDoc *pDoc = (CSubdivisionDoc *)GetDocument();
	ASSERT(pDoc != NULL);
	m_FirstView = false;

	// set up the camera to visualize the whole object
	TRACE("setup camera.");
	//pDoc->m_pMesh->compute_bounding_box();
	CMatrix44 ArcballMatrix = m_Arcball.GetMatrix();

	CVector3d minBound, maxBound;
	minBound.Set(pDoc->m_pmesh->xmin,
		pDoc->m_pmesh->ymin,
		pDoc->m_pmesh->zmin);
	maxBound.Set(pDoc->m_pmesh->xmax,
		pDoc->m_pmesh->ymax,
		pDoc->m_pmesh->zmax);
	minBound = ArcballMatrix * minBound;
	maxBound = ArcballMatrix * maxBound;
	m_Camera.ViewAll(minBound[0],
		maxBound[0],
		minBound[1],
		maxBound[1],
		minBound[2],
		maxBound[2],
		m_Viewport);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	m_Camera.glDraw(m_Viewport);
	glMatrixMode(GL_MODELVIEW);
	TRACE("..ok\n");
}

// init camera parameters
void CSubdivisionView::InitCamera()
{
	// perspective
	CRect rect;
	GetClientRect(&rect);

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set viewport and camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	m_Viewport.SetOrigin(0,0);
	m_Viewport.SetSize(rect.Width(),rect.Height());
	m_Camera.SetHeightAngle(45.);
	m_Camera.SetPosition(0.,0.,5.);
	m_Camera.SetOrientation(0.,1.,0.,0);
	m_Camera.SetNearDistance(.1);
	m_Camera.SetFarDistance(1000.);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);

	::ReleaseDC(hWnd,hDC);
}



void CSubdivisionView::OnPaint()
{
	// device context for painting
	CPaintDC dc(this); 

	// model is stored in Document
	CSubdivisionDoc *pDoc = (CSubdivisionDoc *)GetDocument();
	if(pDoc->m_pmesh == NULL)
		return;
	
	// setup camera (only once)
	ViewAll();

	// useful in multidoc templates
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	
	//wglMakeCurrent(hDC,m_hGLContext);
	wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// shading option
	if(m_SmoothShading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	// culling option
	//glEnable(GL_CULL_FACE);

	// polygon mode (point, line or fill)
	glPolygonMode(GL_FRONT_AND_BACK,m_PolygonMode);

	// set mesh color
	glColor3f(m_MeshColor[0],m_MeshColor[1],m_MeshColor[2]);


	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	glLineWidth(1.5f);

	// drawing
	glPushMatrix();

	// setup viewpoint from current arcball
	m_Arcball.glDraw();

	// enable polygon offset
// 	glEnable(GL_POLYGON_OFFSET_FILL);
// 	glPolygonOffset(3.0f,1.0f);


	// draw the mesh 

	// lighting option
	glEnable(GL_LIGHTING);
	//glDisable(GL_BLEND);

	pDoc->m_pmesh->gl_draw(m_SmoothShading);


	glDisable(GL_LIGHTING);

	glColor3f(1.0f,1.0f,0.0f);

	CString str;

	str.Format(_T("vertices: %d\t faces:  %d\t edges:%d"), 
		pDoc->m_pmesh->get_num_of_vertices_list(),
		pDoc->m_pmesh->get_num_of_faces_list(),
		pDoc->m_pmesh->get_num_of_edges_list()/2);
	
	CRect mr;
	GetClientRect(mr);
	int nHeight = mr.Height();
	glDrawText(10, nHeight - 20, (LPSTR)(LPCTSTR)str);


	// set line mode
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	
	//// edge color
	glColor3f(m_EdgeColor[0],m_EdgeColor[1],m_EdgeColor[2]);

	if(m_showedge)
		pDoc->m_pmesh->gl_draw();

	glDisable(GL_POLYGON_OFFSET_FILL);

	glPopMatrix();

	// Double buffer
	SwapBuffers(dc.m_ps.hdc);
	glFlush();

	// release
	::ReleaseDC(hWnd,hDC);
}

// CSubdivisionView printing

BOOL CSubdivisionView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CSubdivisionView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CSubdivisionView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}


// CSubdivisionView diagnostics

#ifdef _DEBUG
void CSubdivisionView::AssertValid() const
{
	CView::AssertValid();
}

void CSubdivisionView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CSubdivisionDoc* CSubdivisionView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CSubdivisionDoc)));
	return (CSubdivisionDoc*)m_pDocument;
}
#endif //_DEBUG


// CSubdivisionView message handlers




void CSubdivisionView::OnShademodeSmooth()
{
	// TODO: Add your command handler code here
	m_SmoothShading = !m_SmoothShading;
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnUpdateShademodeSmooth(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	pCmdUI->SetCheck(m_SmoothShading==TRUE);
}

void CSubdivisionView::OnViewViewall()
{
	// TODO: Add your command handler code here
	ViewAll(false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnShadeWireframe()
{
	// TODO: Add your command handler code here
	if(	m_PolygonMode == GL_FILL)
		m_PolygonMode = GL_LINE;
	else
		m_PolygonMode =GL_FILL;

	InvalidateRect(NULL,FALSE);                                                                                                                                                    
}

void CSubdivisionView::OnUpdateShadeWireframe(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	pCmdUI->SetCheck(m_PolygonMode==GL_LINE);
}

void CSubdivisionView::OnViewShowedge()
{
	// TODO: Add your command handler code here
	m_showedge = !m_showedge;
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnUpdateViewShowedge(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	pCmdUI->SetCheck(m_showedge==true);
}


// change material
void CSubdivisionView::ChangeMaterial(CString &string,
									  bool update)
{
	float	ambient[]  = {0.0f,0.0f,0.0f,1.0f};
	float	diffuse[]  = {0.0f,0.0f,0.0f,1.0f};
	float	specular[]  = {0.0f,0.0f,0.0f,1.0f};
	float	emission[]  = {0.3f,0.3f,0.3f,1.0f};
	float shininess[] = {0.0f};

	if(string == "Silver")
	{
		// Ambient
		ambient[0] = 0.19225f;
		ambient[1] = 0.19225f;
		ambient[2] = 0.19225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.50754f;
		diffuse[1] = 0.50754f;
		diffuse[2] = 0.50754f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.508273f;
		specular[1] = 0.508273f;
		specular[2] = 0.508273f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(string == "Gold")
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.1995f;
		ambient[2] = 0.0745f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.75164f;
		diffuse[1] = 0.60648f;
		diffuse[2] = 0.22648f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.928281f;
		specular[1] = 0.855802f;
		specular[2] = 0.666065f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(string == "Jade")
	{
		// Ambient
		ambient[0] = 0.135f;
		ambient[1] = 0.2225f;
		ambient[2] = 0.1575f;
		ambient[3] = 0.95f;
		// Diffuse
		diffuse[0] = 0.54f;
		diffuse[1] = 0.89f;
		diffuse[2] = 0.63f;
		diffuse[3] = 0.95f;
		// Specular
		specular[0] = 0.316228f;
		specular[1] = 0.316228f;
		specular[2] = 0.316228f;
		specular[3] = 0.95f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(string == "Light blue")
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.5f;
		ambient[2] = 0.75f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.0f;
		diffuse[1] = 0.5f;
		diffuse[2] = 1.0f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.75f;
		specular[1] = 0.75f;
		specular[2] = 0.75f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 64.0f;
	}

	else if(string == "Emerald")
	{
		// Ambient
		ambient[0] = 0.0215f;
		ambient[1] = 0.1745f;
		ambient[2] = 0.0215f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.07568f;
		diffuse[1] = 0.61424f;
		diffuse[2] = 0.07568f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.633f;
		specular[1] = 0.727811f;
		specular[2] = 0.633f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(string == "Polished silver")
	{
		// Ambient
		ambient[0] = 0.23125f;
		ambient[1] = 0.23125f;
		ambient[2] = 0.23125f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.2775f;
		diffuse[1] = 0.2775f;
		diffuse[2] = 0.2775f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.773911f;
		specular[1] = 0.773911f;
		specular[2] = 0.773911f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 89.6f;
	}

	else if(string == "Chrome")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.25f;
		ambient[2] = 0.25f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.4f;
		diffuse[2] = 0.4f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.774597f;
		specular[2] = 0.774597f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(string == "Copper")
	{
		// Ambient
		ambient[0] = 0.19125f;
		ambient[1] = 0.0735f;
		ambient[2] = 0.0225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.7038f;
		diffuse[1] = 0.27048f;
		diffuse[2] = 0.0828f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.256777f;
		specular[1] = 0.137622f;
		specular[2] = 0.086014f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(string == "Polished gold")
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.2245f;
		ambient[2] = 0.0645f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.34615f;
		diffuse[1] = 0.3143f;
		diffuse[2] = 0.0903f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.797357f;
		specular[1] = 0.723991f;
		specular[2] = 0.208006f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 83.2f;
	}

	else if(string == "Pewter")
	{
		// Ambient
		ambient[0] = 0.105882f;
		ambient[1] = 0.058824f;
		ambient[2] = 0.113725f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.427451f;
		diffuse[1] = 0.470588f;
		diffuse[2] = 0.541176f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.333333f;
		specular[1] = 0.333333f;
		specular[2] = 0.521569f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 9.84615f;
	}

	else if(string == "Obsidian")
	{
		// Ambient
		ambient[0] = 0.05375f;
		ambient[1] = 0.05f;
		ambient[2] = 0.06625f;
		ambient[3] = 0.82f;
		// Diffuse
		diffuse[0] = 0.18275f;
		diffuse[1] = 0.17f;
		diffuse[2] = 0.22525f;
		diffuse[3] = 0.82f;
		// Specular
		specular[0] = 0.332741f;
		specular[1] = 0.328634f;
		specular[2] = 0.346435f;
		specular[3] = 0.82f;
		// Shininess
		shininess[0] = 38.4f;
	}

	else if(string == "Black plastic")
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.0f;
		ambient[2] = 0.0f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.01f;
		diffuse[1] = 0.01f;
		diffuse[2] = 0.01f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.5f;
		specular[1] = 0.5f;
		specular[2] = 0.5f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 32.0f;
	}

	else if(string == "Polished bronze")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.148f;
		ambient[2] = 0.006475f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.2368f;
		diffuse[2] = 0.1036f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.458561f;
		specular[2] = 0.200621f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}


	else if(string == "Polished copper")
	{
		// Ambient
		ambient[0] = 0.2295f;
		ambient[1] = 0.08825f;
		ambient[2] = 0.0275f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.5508f;
		diffuse[1] = 0.2118f;
		diffuse[2] = 0.066f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.580594f;
		specular[1] = 0.223257f;
		specular[2] = 0.0695701f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(string == "Pearl")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.20725f;
		ambient[2] = 0.20725f;
		ambient[3] = 0.922f;
		// Diffuse
		diffuse[0] = 1.0f;
		diffuse[1] = 0.829f;
		diffuse[2] = 0.829f;
		diffuse[3] = 0.922f;
		// Specular
		specular[0] = 0.296648f;
		specular[1] = 0.296648f;
		specular[2] = 0.296648f;
		specular[3] = 0.922f;
		// Shininess
		shininess[0] = 11.264f;
	}

	else if(string == "Ruby")
	{
		// Ambient
		ambient[0] = 0.1745f;
		ambient[1] = 0.01175f;
		ambient[2] = 0.01175f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.61424f;
		diffuse[1] = 0.04136f;
		diffuse[2] = 0.04136f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.727811f;
		specular[1] = 0.626959f;
		specular[2] = 0.626959f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(string == "Turquoise")
	{
		// Ambient
		ambient[0] = 0.1f;
		ambient[1] = 0.18725f;
		ambient[2] = 0.1745f;
		ambient[3] = 0.8f;
		// Diffuse
		diffuse[0] = 0.396f;
		diffuse[1] = 0.74151f;
		diffuse[2] = 0.69102f;
		diffuse[3] = 0.8f;
		// Specular
		specular[0] = 0.297254f;
		specular[1] = 0.30829f;
		specular[2] = 0.306678f;
		specular[3] = 0.8f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(string == "Brass")
	{
		// Ambient
		ambient[0] = 0.329412f;
		ambient[1] = 0.223529f;
		ambient[2] = 0.027451f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.780392f;
		diffuse[1] = 0.268627f;
		diffuse[2] = 0.113725f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.992157f;
		specular[1] = 0.741176f;
		specular[2] = 0.807843f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 27.8974f;
	}

	else
		;
	// apply
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambient);
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diffuse);
	glMaterialfv( GL_FRONT, GL_SPECULAR,  specular);
	glMaterialfv( GL_FRONT, GL_SHININESS, shininess);
	glMaterialfv( GL_FRONT, GL_EMISSION,  emission);

	if(update)
		InvalidateRect(NULL,FALSE);
}



void CSubdivisionView::OnMaterialSilver()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Silver"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialGold()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Gold"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialPearl()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Pearl"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialBrass()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Brass"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialBlackplastic()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Black plastic"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialRuby()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Ruby"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialJade()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Jade"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialLightblue()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Light blue"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialEmerald()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Emerald"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialPolishedsilver()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Polished silver"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialChrome()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Chrome"),false);
	InvalidateRect(NULL,FALSE);

}

void CSubdivisionView::OnMaterialCopper()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Copper"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialPolishedgold()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Polished gold"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialPewter()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Pewter"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialObsidian()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Obsidian"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialPolishedbronze()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Polished bronze"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialPolishedcopper()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Polished copper"),false);
	InvalidateRect(NULL,FALSE);
}

void CSubdivisionView::OnMaterialTurquoise()
{
	// TODO: Add your command handler code here
	ChangeMaterial(CString("Turquoise"),false);
	InvalidateRect(NULL,FALSE);

}

///////////////////////////////// CREATE OPENGL FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This function creates a windows font like "Arial" and returns a display list ID
/////
///////////////////////////////// CREATE OPENGL FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
UINT CSubdivisionView::CreateOpenGLFont(LPWSTR strFontName, int height)	// Build Our Bitmap Font
{
	UINT	fontListID = 0;								// This will hold the base ID for our display list
	HFONT	hFont;										// This will store the handle to our font

	// Here we generate the lists for each character we want to use.
	// This function then returns the base pointer, which will be 1 because
	// we haven't created any other lists.  If we generated another list after
	// this, the base pointer would be at 257 since the last one used was 256 (which is MAX_CHARS)
	fontListID = glGenLists(MAX_CHARS);					// Generate the list for the font

	// Now we actually need to create the font.  We use a windows function called:
	// CreateFont() that returns a handle to a font (HFONT).  Our CreateOpenGLFont()
	// function allows us to pass in a name and height.  For simplistic reasons, I left
	// other options out, but feel free to add them to your function (like bold, italic, width..)

	hFont = CreateFont(	height,							// Our desired HEIGHT of the font
		0,								// The WIDTH (If we leave this zero it will pick the best width depending on the height)
		0,								// The angle of escapement
		0,								// The angle of orientation
		FW_BOLD,						// The font's weight (We want it bold)
		FALSE,							// Italic - We don't want italic
		FALSE,							// Underline - We don't want it underlined
		FALSE,							// Strikeout - We don't want it strikethrough
		ANSI_CHARSET,					// This is the type of character set
		OUT_TT_PRECIS,					// The Output Precision
		CLIP_DEFAULT_PRECIS,			// The Clipping Precision
		ANTIALIASED_QUALITY,			// The quality of the font - We want anitaliased fonts
		FF_DONTCARE|DEFAULT_PITCH,		// The family and pitch of the font.  We don't care.
		strFontName);					// The font name (Like "Arial", "Courier", etc...)

	CClientDC dc(this);
	// Now that we have created a new font, we need to select that font into our global HDC.
	// We store the old font so we can select it back in when we are done to avoid memory leaks.
	hOldFont = (HFONT)SelectObject(dc.GetSafeHdc(), hFont);

	// This function does the magic.  It takes the current font selected in
	// the hdc and makes bitmaps out of each character.  These are called glyphs.
	// The first parameter is the HDC that holds the font to be used.
	// The second parameters is the ASCII value to start from, which is zero in our case.
	// The third parameters is the ASCII value to end on (255 is the last of the ASCII values so we minus 1 from MAX_CHARS)
	// The last parameter is the base pointer for the display lists being used.  This should be 1.

	wglUseFontBitmaps(dc.GetSafeHdc(), 0, MAX_CHARS - 1, fontListID);	// Builds 255 bitmap characters

	return fontListID;									// Return the ID to the display list to use later
}

///////////////////////////////// POSITION TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This function sets the drawing position for the text
/////
///////////////////////////////// POSITION TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
void CSubdivisionView::PositionText( int x, int y )
{
	// If you are to use this font code for your applications,
	// you must be aware that you cannot position the font in 3D,
	// which means you can't rotate and scale it.  That will be covered in
	// the next font tutorial.  BUT, though that might be a drag, this code
	// is useful because when you display the text, it will always be on top
	// of everything else.  This is good if the camera is moving around, and you
	// don't want the text to move.  If the text was positioned in 3D you would have
	// to come up with a tricky way of making it always render in front of the camera.
	// To do this, we need to set the Raster Position.  That is the position that OpenGL
	// starts drawing at.  Since it's in floating point, it's not very intuitive, so what
	// we do is create a new view port, and then always draw the text at (0, 0, 0) in that
	// view port.  The weird part is that the Y is flipped, so (0, 0) is the bottom left corner.
	// Below we do some simple math to flip it back to normal.

	// Before we create a new view port, we need to save the current one we have.
	// This saves our transform (matrix) information and our current viewport information.
	// At the end of this function we POP it back.
	glPushAttrib( GL_TRANSFORM_BIT | GL_VIEWPORT_BIT );

	// Here we use a new projection and modelview matrix to work with.
	glMatrixMode( GL_PROJECTION );						// Set our matrix to our projection matrix
	glPushMatrix();										// Push on a new matrix to work with
	glLoadIdentity();									// reset the matrix
	glMatrixMode( GL_MODELVIEW );						// Set our matrix to our model view matrix
	glPushMatrix();										// Push on a new matrix to work with
	glLoadIdentity();									// Reset that matrix

	// Because the Y is flipped, we want 0 to be at the top, not bottom.
	// If we subtract the font height from the screen height, that should display the
	// font at the top of the screen (if they passed in 0 for Y), but then we subtract
	// the Y from that to get the desired position.  Since the font's drawing point is
	// at the base line of the font, we needed to subtract the font height to make sure
	// if they passed in (0, 0) it wouldn't be off screen.  If you view this in window mode,
	// the top of the window will cut off part of the font, but in full screen it works fine.
	// You just need to add about 25 to the Y to fix that for window mode.

	// Now we create another view port (that is why we saved the old one above).
	// Since glViewPort takes the lower LEFT corner, we needed to change the Y
	// to make it more intuitive when using PositionText().  We minus 1 from the X and Y
	// because 0 is taken into account with the position.  The next 2 parameters are set
	// to 0 for the width and height so it will always draw in the middle of that position.
	// glRasterPos4f() takes (0, 0, 0) as the middle of the viewport, so if we give it a small
	// width/height it will draw at the X and Y given.  Sounds strange, to test this, try
	// using glRasterPos4f(0, 0, 0, 1) instead of PositionText() and you will see, everything
	// will be drawn from the middle.

	glViewport( x, y, 0, 0 );					// Create a new viewport to draw into

	// This is the most important function in here.  This actually positions the text.
	// The parameters are (x, y, z, w).  w should always be 1 , it's a clip coordinate.
	// don't worry about that though.  Because we set the projection and modelview matrix
	// back to the beginning (through LoadIdentity()), the view port is looking at (0, 0, 0).
	// This is the middle, so if we set the drawing position to the middle, it will draw at our
	// X and Y because the width/height of the viewport is 0, starting at X and Y.
	// You can actually call this function (or glRasterPos2f(0, 0)) instead of PositionText(),
	// but it is in floating point and doesn't work as nicely.  You will see why if you try.

	glRasterPos4f( 0, 0, 0, 1 );						// Set the drawing position

	// Now that we positioned the raster position, any text we draw afterwards will start
	// from that position.  Now we just have to put everything else back to normal.

	glPopMatrix();										// Pop the current modelview matrix off the stack
	glMatrixMode( GL_PROJECTION );						// Go back into projection mode
	glPopMatrix();										// Pop the projection matrix off the stack

	glPopAttrib();										// This restores our TRANSFORM and VIEWPORT attributes
}

///////////////////////////////// GL DRAW TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This draws 2D text onto the screen using OpenGL, given an X and Y position 
/////
///////////////////////////////// GL DRAW TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
void CSubdivisionView::glDrawText(int x, int y, const char *strString, ...)
{
	char		strText[256];							// This will hold our text to display
	va_list		argumentPtr;							// This will hold the pointer to the argument list

	// If you have never used a va_list, listen up.  Remember printf()?
	// or sprintf()?  Well, you can add unlimited arguments into the text like:
	// printf("My name is %s and I am %d years old!", strName, age);
	// Well, that is what va_list's do.  

	// First we need to check if there was even a string given
	if (strString == NULL)								// Check if a string was given
		return;											// Don't render anything then

	// First we need to parse the string for arguments given
	// To do this we pass in a va_list variable that is a pointer to the list of arguments.
	// Then we pass in the string that holds all of those arguments.
	va_start(argumentPtr, strString);					// Parse the arguments out of the string

	// Then we use a special version of sprintf() that takes a pointer to the argument list.
	// This then does the normal sprintf() functionality.
	vsprintf(strText, strString, argumentPtr);			// Now add the arguments into the full string

	va_end(argumentPtr);								// This resets and frees the pointer to the argument list.

	// Before we draw the text, we need to position it with our own function.
	PositionText(x, y);									// Call our own function to position the text on screen

	// Now, before we set the list base, we need to save off the current one.
	glPushAttrib(GL_LIST_BIT);							// This saves the list base information

	// Then we want to set the list base to the font's list base, which should be 1 in our case.
	// That way when we call our display list it will start from the font's lists'.
	glListBase(g_FontListID);							// This sets the lists base

	// Now comes the actually rendering.  We pass in the length of the string,
	// then the data types (which are characters so its a UINT), then the actually char array.
	// This will then take the ASCII value of each character and associate it with a bitmap.
	glCallLists((GLsizei)strlen(strText), GL_UNSIGNED_BYTE, strText);
	glPopAttrib();										// Return the display list back to it's previous state
}

void CSubdivisionView::DestroyFont()										
{
	glDeleteLists(g_FontListID, MAX_CHARS);				// Free the display list
	//SelectObject(g_hDC, hOldFont);						// Select the old font back in so we don't have memory leaks
}




bool CSubdivisionView::LoadTexture(LPTSTR szFileName, GLuint &texid)
{
	HBITMAP hBMP;              // Handle Of The Bitmap
	BITMAP BMP;              // Bitmap Structure

	glGenTextures(1, &texid);           // Create The Texture
	hBMP = (HBITMAP)LoadImage(GetModuleHandle(NULL), szFileName, IMAGE_BITMAP, 0, 0, LR_CREATEDIBSECTION | LR_LOADFROMFILE);

	if (!hBMP)               // Does The Bitmap Exist?
		return FALSE;             // If Not Return False

	GetObject(hBMP, sizeof(BMP), &BMP);         // Get The Object
	// hBMP:        Handle To Graphics Object
	// sizeof(BMP): Size Of Buffer For Object Information
	// &BMP:        Buffer For Object Information

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);        // Pixel Storage Mode (Word Alignment / 4 Bytes)

	// Typical Texture Generation Using Data From The Bitmap
	glBindTexture(GL_TEXTURE_2D, texid);         // Bind To The Texture ID
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // Linear Min Filter
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Linear Mag Filter
	glTexImage2D(GL_TEXTURE_2D, 0, 3, BMP.bmWidth, BMP.bmHeight, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, BMP.bmBits);

	DeleteObject(hBMP);             // Delete The Object

	return TRUE;              // Loading Was Successful
}


void CSubdivisionView::OnSea(){
	LoadTexture(_T("Wave.bmp"), textureID);
	SetTimer(1, 10, NULL);
}

// void CSubdivisionView::CopyTo(Mesh3D *m_pmesh1, Mesh3D *m_pmesh2)
// {
// 	PTR_FACE_LIST faces_list = m_pmesh1->get_faces_list();
// 	FACE_ITER fiter = faces_list->begin();
// 
// 	PTR_VERTEX_LIST vertices_list = m_pmesh1->get_vertices_list();
// 	VERTEX_ITER viter = vertices_list->begin();
// 
// 	HE_edge *ee, *tempEdge;
// 
// 	for (; viter != vertices_list->end(); viter++){
// 		m_pmesh2->insert_vertex((*viter)->x, (*viter)->y, (*viter)->z);
// 	}
// 
// 	for (; fiter != faces_list->end(); fiter++){
// 		VERTEX_LIST vl;
// 		ee = (*fiter)->edge;
// 		tempEdge = ee;
// 		do
// 		{
// 			vl.push_back(tempEdge->vert);
// 			tempEdge = tempEdge->next;
// 		} while (tempEdge != ee);
// 
// 		m_pmesh2->insert_face(vl);
// 	}
// 
// 
// }

void CSubdivisionView::OnTimer(UINT_PTR nIDEvent)
{
	if (nIDEvent == 1){

		static float time = 0;

		//boat swing
		
			// device context for painting
			CPaintDC dc(this);

			// model is stored in Document
			CSubdivisionDoc *pDoc = (CSubdivisionDoc *)GetDocument();
			if (pDoc->m_pmesh == NULL)
				return;

			// setup camera (only once)
			ViewAll();

			// useful in multidoc templates
			HWND hWnd = GetSafeHwnd();
			HDC hDC = ::GetDC(hWnd);

			//wglMakeCurrent(hDC,m_hGLContext);
			wglMakeCurrent(dc.m_ps.hdc, m_hGLContext);

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			// shading option
			if (m_SmoothShading)
				glShadeModel(GL_SMOOTH);
			else
				glShadeModel(GL_FLAT);

			// culling option
			//glEnable(GL_CULL_FACE);

			// polygon mode (point, line or fill)
			glPolygonMode(GL_FRONT_AND_BACK, m_PolygonMode);

			// set mesh color
			glColor3f(m_MeshColor[0], m_MeshColor[1], m_MeshColor[2]);


			glEnable(GL_LINE_SMOOTH);
			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			glLineWidth(1.5f);

			// drawing
			glPushMatrix();

			// setup viewpoint from current arcball
			m_Arcball.glDraw();

			// enable polygon offset
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0f,1.0f);


			// draw the mesh 

			// lighting option
			glEnable(GL_LIGHTING);

			glRotatef(10 * sin(time + 0.5), -1, 0, 0);
			pDoc->m_pmesh->gl_draw(m_SmoothShading);

			glDisable(GL_LIGHTING);


			// set line mode
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

			//// edge color
			glColor3f(m_EdgeColor[0], m_EdgeColor[1], m_EdgeColor[2]);

			if (m_showedge)
				pDoc->m_pmesh->gl_draw();

			glDisable(GL_POLYGON_OFFSET_FILL);

			glPopMatrix();

/////////////////////////////////////////////////////////////////////////////
		

		//wave

		Mesh3D* m_pwave=new Mesh3D;

		//Create Sea
		OnCreateSea(m_pwave, time);
		time += 0.1;

		//wglMakeCurrent(hDC,m_hGLContext);
		wglMakeCurrent(dc.m_ps.hdc, m_hGLContext);

		//load texture
		glBindTexture(GL_TEXTURE_2D, textureID);


		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// shading option
		if (m_SmoothShading)
			glShadeModel(GL_SMOOTH);
		else
			glShadeModel(GL_FLAT);

		// culling option
		//glEnable(GL_CULL_FACE);

		// polygon mode (point, line or fill)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// set mesh color
		//glColor3f(m_MeshColor[0], m_MeshColor[1], m_MeshColor[2]);

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(0);


		//texture
		glEnable(GL_TEXTURE_2D);
		glShadeModel(GL_FLAT);
		glBindTexture(GL_TEXTURE_2D, textureID);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		// drawing
		glPushMatrix();

		// setup viewpoint from current arcball
		m_Arcball.glDraw();

		// enable polygon offset
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(-1.0f,-1.0f);


		// draw the mesh 

		// lighting option
		glEnable(GL_LIGHTING);

		glDepthMask(GL_FALSE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4f(1.0f, 1.0f, 1.0f, 0.7f);

		//计算贴图
		float xLength = 2 * fabs(xWith) + 1;
		float zLength = 2 * fabs(zDepth) + 1;
		m_pwave->gl_draw_wave(xLength,zLength,m_SmoothShading);

		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_TEXTURE_2D);
		glDepthMask(GL_TRUE);

		glPopMatrix();

		///////////////////////////////////////////////////////////////////

		delete m_pwave;

		// Double buffer
		SwapBuffers(dc.m_ps.hdc);
		glFlush();

		// release
		::ReleaseDC(hWnd, hDC);
	}
	CView::OnTimer(nIDEvent);
}



Mesh3D* CSubdivisionView::OnCreateSea(Mesh3D* m_pmesh, float a){

	float x, y, z;
	HE_vert* vv;
	HE_face* ff;
	std::vector<int>  vert_id_list;
	VERTEX_LIST  v_list;

	for (x = xWith; x >- xWith; x -= STEP){
		vv = m_pmesh->insert_vertex(x, yDepth + sin(a)*0.3, STEP - zDepth);
		vert_id_list.push_back(vv->id);
	}


	for (z = -zDepth / 2; z > zDepth; z -= STEP){

		for (float x = xWith; x >- xWith; x -= STEP){

			vv = m_pmesh->get_vertex(vert_id_list[0]);
			vert_id_list.erase(vert_id_list.begin());

			x = vv->x; y = vv->y; z = vv->z;
			m_pmesh->insert_vertex(x, y, z);
			v_list.push_back(vv);
			vv = m_pmesh->insert_vertex(x - STEP, y, z);
			v_list.push_back(vv);

			z -= STEP; y = sin(z + a)*0.3 + yDepth;
			vv = m_pmesh->insert_vertex(x - STEP, y, z);
			v_list.push_back(vv);
			vv = m_pmesh->insert_vertex(x, y, z);
			v_list.push_back(vv);

			vert_id_list.push_back(vv->id);

			ff = m_pmesh->insert_face(v_list);
			ff->tag = true;

			v_list.clear();
		}
	}
	vert_id_list.clear();
	m_pmesh->update_mesh();

	return m_pmesh;
}