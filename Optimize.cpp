#include "StdAfx.h"
#include "Optimize.h"


//////////////////////////////////////////////////////////////////////////
void MeshOptimization::evalfunc(int N, double* x, double *prev_x, double* f, double* g)
{
	*f = 0;
	for (int i = 0; i < N; i+=2)
	{
		double T1 = 1 - x[i];
		double T2 = 10*(x[i+1]-x[i]*x[i]);
		*f += T1*T1+T2*T2;
		g[i+1]   = 20*T2;
		g[i] = -2*(x[i]*g[i+1]+T1);
	}
}
//////////////////////////////////////////////////////////////////////////
void MeshOptimization::newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm)
{
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
}
//////////////////////////////////////////////////////////////////////////
void MeshOptimization::evalfunc_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)
{
	//the following code is not optimal if the pattern of hessian matrix is fixed.
	if (m_sparse_matrix)
	{
		delete m_sparse_matrix;
	}

	m_sparse_matrix = new Lite_Sparse_Matrix(N, N, SYM_LOWER, CCS, FORTRAN_TYPE, true);

	m_sparse_matrix->begin_fill_entry();

	static bool first = true;
	double *diag = m_sparse_matrix->get_diag();

	if (first)
	{
		// you need to update f and g
		*f = 0;
		double tmp;
		for (int i = 0; i < N; i+=2)
		{
			tmp = x[i]*x[i];
			double T1 = 1 - x[i];
			double T2 = 10*(x[i+1]-tmp);
			*f += T1*T1+T2*T2;
			g[i+1]   = 20*T2;
			g[i] = -2*(x[i]*g[i+1]+T1);
			diag[i] = 2+1200*tmp-400*x[i+1];
			diag[i+1] = 200;
			m_sparse_matrix->fill_entry(i, i+1, -400*x[i]);
		}
	}
	else
	{
		for (int i = 0; i < N; i+=2)
		{
			diag[i] = 2+1200*x[i]*x[i]-400*x[i+1];
			diag[i+1] = 200;
			m_sparse_matrix->fill_entry(i, i+1, -400*x[i]);
		}
	}

	m_sparse_matrix->end_fill_entry();

	hessian.set_diag(m_sparse_matrix->get_diag());
	hessian.set_values(m_sparse_matrix->get_values());
	hessian.set_rowind(m_sparse_matrix->get_rowind());
	hessian.set_colptr(m_sparse_matrix->get_colptr());
	hessian.set_nonzeros(m_sparse_matrix->get_nonzero());
	first = false;
}
//////////////////////////////////////////////////////////////////////////
void MeshOptimization::Optimize_by_HLBFGS(int N, double *init_x, int num_iter, int M, int T, bool with_hessian)
{
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	info[4] = num_iter;
	info[6] = T;
	info[7] = with_hessian?1:0;
	info[10] = 0;
	info[11] = 1;

	if (with_hessian)
	{
		HLBFGS(N, M, init_x, evalfunc, evalfunc_h, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}
	else
	{
		HLBFGS(N, M, init_x, evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}

}
