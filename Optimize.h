#pragma once

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include "HLBFGS/HLBFGS.h"
#include "HLBFGS/Lite_Sparse_Matrix.h"

#include <iostream>

#include "Mesh3D.h"
#include "SubdivisionDoc.h"


Lite_Sparse_Matrix* m_sparse_matrix;
Mesh3D *m_pmesh;


void Init();
void evalfunc(int N, double* x, double *prev_x, double* f, double* g);
void newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm);
void evalfunc_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian);
void Optimize_by_HLBFGS(int N, double *init_x, int num_iter, int M, int T, bool with_hessian);
void UpdateMesh(double* x);



void Init()
{
	m_sparse_matrix=0;

#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	std::cout.precision(16);
	std::cout << std::scientific;

	int N = m_pmesh->get_num_of_vertices_list()*3;
	std::vector<double> x(N);

	for (int i = 0; i < m_pmesh->get_num_of_vertices_list(); i++)
	{
		x[i * 3] = m_pmesh->get_vertex(i)->x;
		x[i * 3 + 1] = m_pmesh->get_vertex(i)->y;
		x[i * 3 + 2] = m_pmesh->get_vertex(i)->z;
	}

	int M = 7;
	int T = 0;

	//use Hessian
	// if M = 0, T = 0, it is Newton
	//Optimize_by_HLBFGS(N, &x[0], 1000, M, T, true);

	//without Hessian
	Optimize_by_HLBFGS(N, &x[0], 1000, M, T, false);  // it is LBFGS(M) actually, T is not used

	if (m_sparse_matrix)
		delete m_sparse_matrix;
}

