#pragma once

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include "HLBFGS/Lite_Sparse_Matrix.h"

#include "Mesh3D.h"

class OptimizeParameter
{
public:
	OptimizeParameter(Mesh3D* m);
	~OptimizeParameter();

	Lite_Sparse_Matrix* m_sparse_matrix;
	Mesh3D *mm_pmesh;

};
