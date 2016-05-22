#include "StdAfx.h"
#include "Optimize.h"


void MeshOptimization::Init()// :m_pmesh(m), m_sparse_matrix(0)
{

#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	std::cout.precision(16);
	std::cout << std::scientific;

	int N = opp->m_pmesh->get_num_of_vertices_list() * 3;
	std::vector<double> x(N);

	for (int i = 0; i < opp->m_pmesh->get_num_of_vertices_list(); i++)
	{
		x[i * 3] = opp->m_pmesh->get_vertex(i)->x;
		x[i * 3 + 1] = opp->m_pmesh->get_vertex(i)->y;
		x[i * 3 + 2] = opp->m_pmesh->get_vertex(i)->z;
	}

	int M = 7;
	int T = 0;

	//use Hessian
	// if M = 0, T = 0, it is Newton
	//Optimize_by_HLBFGS(N, &x[0], 1000, M, T, true);

	//without Hessian
	Optimize_by_HLBFGS(N, &x[0], 1000, M, T, false);  // it is LBFGS(M) actually, T is not used

	if (opp->m_sparse_matrix)
		delete opp->m_sparse_matrix;
}

/////////////////////////////////////////////////////////////////////////
void UpdateMesh(double* x)
{
	PTR_VERTEX_LIST vertex_list = MeshOptimization::opp->m_pmesh->get_vertices_list();
	VERTEX_ITER     viter = vertex_list->begin();

	for (int i = 0; viter != vertex_list->end(); viter++, i++)
	{
		(*viter)->x = x[i];
		(*viter)->y = x[i + 1];
		(*viter)->z = x[i + 2];
	}

	MeshOptimization::opp->m_pmesh->compute_all_normal();
}
//////////////////////////////////////////////////////////////////////////
void evalfunc(int N, double* x, double *prev_x, double* f, double* g)
{
	UpdateMesh(x);

	*f = 0;
	double product;
	double edge_x, edge_y, edge_z;

	PTR_VERTEX_LIST vertex_list = MeshOptimization::opp->m_pmesh->get_vertices_list();
	VERTEX_ITER     viter = vertex_list->begin();

	for (int i = 0; viter != vertex_list->end(); viter++, i++)
	{
		HE_edge* t_edge = (*viter)->edge;
		do 
		{
			//对于当前结点指出去的边
			edge_x = t_edge->vert->x - (*viter)->x;
			edge_y = t_edge->vert->y - (*viter)->y;
			edge_z = t_edge->vert->z - (*viter)->z;

			product = edge_x*(*viter)->normal[0] + edge_y*(*viter)->normal[1] + edge_z*(*viter)->normal[2];
			*f += product*product;

			g[i * 3] = -2 * product*(*viter)->normal[0];
			g[i * 3 + 1] = -2 * product*(*viter)->normal[1];
			g[i * 3 + 2] = -2 * product*(*viter)->normal[2];


			//对于指向当前结点的边，在此顶点中不做处理
			t_edge = t_edge->next;

		} while (t_edge!=(*viter)->edge);
		
	}
}
//////////////////////////////////////////////////////////////////////////
void newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm)
{
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
}
//////////////////////////////////////////////////////////////////////////
void evalfunc_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)
{
	//the following code is not optimal if the pattern of hessian matrix is fixed.
	if (MeshOptimization::opp->m_sparse_matrix)
	{
		delete MeshOptimization::opp->m_sparse_matrix;
	}

	MeshOptimization::opp->m_sparse_matrix = new Lite_Sparse_Matrix(N, N, SYM_LOWER, CCS, FORTRAN_TYPE, true);

	MeshOptimization::opp->m_sparse_matrix->begin_fill_entry();

	static bool first = true;
	double *diag = MeshOptimization::opp->m_sparse_matrix->get_diag();

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
			MeshOptimization::opp->m_sparse_matrix->fill_entry(i, i + 1, -400 * x[i]);
		}
	}
	else
	{
		for (int i = 0; i < N; i+=2)
		{
			diag[i] = 2+1200*x[i]*x[i]-400*x[i+1];
			diag[i+1] = 200;
			MeshOptimization::opp->m_sparse_matrix->fill_entry(i, i + 1, -400 * x[i]);
		}
	}

	MeshOptimization::opp->m_sparse_matrix->end_fill_entry();

	hessian.set_diag(MeshOptimization::opp->m_sparse_matrix->get_diag());
	hessian.set_values(MeshOptimization::opp->m_sparse_matrix->get_values());
	hessian.set_rowind(MeshOptimization::opp->m_sparse_matrix->get_rowind());
	hessian.set_colptr(MeshOptimization::opp->m_sparse_matrix->get_colptr());
	hessian.set_nonzeros(MeshOptimization::opp->m_sparse_matrix->get_nonzero());
	first = false;
}
//////////////////////////////////////////////////////////////////////////
void Optimize_by_HLBFGS(int N, double *init_x, int num_iter, int M, int T, bool with_hessian)
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

	//typedef void(*evalfunc)(int N, double* x, double *prev_x, double* f, double* g);

	if (with_hessian)
	{
		HLBFGS(N, M, init_x, evalfunc, evalfunc_h, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}
	else
	{
		HLBFGS(N, M, init_x, evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}

}
//////////////////////////////////////////////////////////////////////////
