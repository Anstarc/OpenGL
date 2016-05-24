#include "StdAfx.h"
#include "Optimize.h"


OptimizeParameter* MeshOptimization::opp;
/////////////////////////////////////////////////////////////////////////
void MeshOptimization::UpdateMesh(double* x)
{
	PTR_VERTEX_LIST vertex_list = opp->mm_pmesh->get_vertices_list();
	VERTEX_ITER     viter = vertex_list->begin();

	for (int i = 0; viter != vertex_list->end(); viter++, i++)
	{
		opp->mm_pmesh->get_vertex(i)->x = x[i * 3];
		opp->mm_pmesh->get_vertex(i)->y = x[i * 3+1];
		opp->mm_pmesh->get_vertex(i)->z = x[i * 3+2];
// 		(*viter)->x = x[(*viter)->id * 3];
// 		(*viter)->y = x[(*viter)->id * 3 + 1];
// 		(*viter)->z = x[(*viter)->id * 3 + 2];
	}

	opp->mm_pmesh->compute_all_normal();
}

//////////////////////////////////////////////////////////////////////////
void MeshOptimization::evalfunc(int N, double* x, double *prev_x, double* f, double* g)
{
	UpdateMesh(x);

	*f = 0;
	double product;
	double edge_x, edge_y, edge_z;

	PTR_VERTEX_LIST vertex_list = opp->mm_pmesh->get_vertices_list();
	VERTEX_ITER     viter = vertex_list->begin();

	for (; viter != vertex_list->end(); viter++)
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

			g[(*viter)->id * 3] += -2 * product*(*viter)->normal[0];
			g[(*viter)->id * 3 + 1] += -2 * product*(*viter)->normal[1];
			g[(*viter)->id * 3 + 2] += -2 * product*(*viter)->normal[2];

			//TRACE("导数%d :%f\n", (*viter)->id * 3, g[(*viter)->id * 3]);

			//对于指向当前结点的边，在此顶点中不做处理
			t_edge = t_edge->next;

		} while (t_edge!=(*viter)->edge);
		
	}
}
//////////////////////////////////////////////////////////////////////////
void MeshOptimization::newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm)
{
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
	TRACE("########################\n");
	TRACE("迭代次数:%d\n", iter);
	TRACE("函数调用次数:%d\n", call_iter);
	TRACE("F(x):%f\n", *f);
	TRACE("dF(x):%f\n", *g);
}
//////////////////////////////////////////////////////////////////////////
void MeshOptimization::evalfunc_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)
{
	//the following code is not optimal if the pattern of hessian matrix is fixed.
	if (opp->m_sparse_matrix)
	{
		delete opp->m_sparse_matrix;
	}

	opp->m_sparse_matrix = new Lite_Sparse_Matrix(N, N, SYM_LOWER, CCS, FORTRAN_TYPE, true);

	opp->m_sparse_matrix->begin_fill_entry();

	static bool first = true;
	double *diag = opp->m_sparse_matrix->get_diag();

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
			opp->m_sparse_matrix->fill_entry(i, i + 1, -400 * x[i]);
		}
	}
	else
	{
		for (int i = 0; i < N; i+=2)
		{
			diag[i] = 2+1200*x[i]*x[i]-400*x[i+1];
			diag[i+1] = 200;
			opp->m_sparse_matrix->fill_entry(i, i + 1, -400 * x[i]);
		}
	}

	opp->m_sparse_matrix->end_fill_entry();

	hessian.set_diag(opp->m_sparse_matrix->get_diag());
	hessian.set_values(opp->m_sparse_matrix->get_values());
	hessian.set_rowind(opp->m_sparse_matrix->get_rowind());
	hessian.set_colptr(opp->m_sparse_matrix->get_colptr());
	hessian.set_nonzeros(opp->m_sparse_matrix->get_nonzero());
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
	info[5] = 1;
	info[6] = T;
	info[7] = with_hessian ? 1 : 0;
	info[10] = 0;
	//info[11] = 1;

	if (with_hessian)
	{
		HLBFGS(N, M, init_x, evalfunc, evalfunc_h, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}
	else
	{
		HLBFGS(N, M, init_x, evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}

}
/////////////////////////////////////////////////////////////////////////
void MeshOptimization::Init()
{

#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	std::cout.precision(16);
	std::cout << std::scientific;

	int N = opp->mm_pmesh->get_num_of_vertices_list() * 3;
	std::vector<double> x(N);


	PTR_VERTEX_LIST vertex_list = opp->mm_pmesh->get_vertices_list();
	VERTEX_ITER     viter = vertex_list->begin();

	for (int i=0; viter != vertex_list->end(); viter++,i++)
	{
		x[i * 3] = opp->mm_pmesh->get_vertex(i)->x;
		x[i * 3 + 1] = opp->mm_pmesh->get_vertex(i)->y;
		x[i * 3+2] = opp->mm_pmesh->get_vertex(i)->z;
// 		x[(*viter)->id * 3] = (*viter)->x;
// 		x[(*viter)->id * 3 + 1] = (*viter)->y;
// 		x[(*viter)->id * 3 + 2] = (*viter)->z;
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
