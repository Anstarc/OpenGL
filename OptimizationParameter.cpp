#include "stdafx.h"

#include "OptimizationParameter.h"

OptimizeParameter::OptimizeParameter(Mesh3D* m) :mm_pmesh(m)
{
	m_sparse_matrix = 0;
}

OptimizeParameter::~OptimizeParameter()
{
}