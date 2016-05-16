#pragma once

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include "HLBFGS/HLBFGS.h"
#include "HLBFGS/Lite_Sparse_Matrix.h"

#include <iostream>


class MeshOptimization
{

private:
	Lite_Sparse_Matrix* m_sparse_matrix = 0;


public:
	MeshOptimization();
	~MeshOptimization();

	void evalfunc(int N, double* x, double *prev_x, double* f, double* g);
	void newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm);
	void evalfunc_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian);
	void Optimize_by_HLBFGS(int N, double *init_x, int num_iter, int M, int T, bool with_hessian);



};

MeshOptimization::MeshOptimization()
{
#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	std::cout.precision(16);
	std::cout << std::scientific;

	int N = 1000;
	std::vector<double> x(N);

	for (int i = 0; i < N / 2; i++)
	{
		x[2 * i] = -1.2;
		x[2 * i + 1] = 1.0;
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

MeshOptimization::~MeshOptimization()
{
}
