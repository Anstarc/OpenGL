#include "stdafx.h"
#include "CreateSea.h"
#include <vector>


#define PI 3.1415926536
#define zDepth -10.0
#define yDepth -0.4
#define xWith 10.0
#define STEP 0.5


Mesh3D* CreateSea::OnCreateSea(Mesh3D* m_pmesh,float a){
	
	float x, y, z;
	HE_vert* vv;
	HE_face* ff;
	std::vector<int>  vert_id_list;

	for (x = -xWith; x <= xWith; x += STEP){
		vv=m_pmesh->insert_vertex(x, yDepth+sin(a), STEP-zDepth);
		vert_id_list.push_back(vv->id);
	}


		

	for (z = -zDepth/2; z > zDepth; z -= STEP){

		for (float x = -xWith; x <= xWith; x += STEP){

			VERTEX_LIST  v_list;

			vv = m_pmesh->get_vertex(vert_id_list[0]);
			vert_id_list.erase(vert_id_list.begin());

			x = vv->x; y = vv->y; z = vv->z;
			m_pmesh->insert_vertex(x, y, z);
			v_list.push_back(vv);
			vv = m_pmesh->insert_vertex(x + STEP, y, z);
			v_list.push_back(vv);

			z -= STEP; y = sin(z+a)*0.3+yDepth;
			vv = m_pmesh->insert_vertex(x + STEP, y, z);
			v_list.push_back(vv);
			vv = m_pmesh->insert_vertex(x, y, z);
			v_list.push_back(vv);

			vert_id_list.push_back(vv->id);

			ff=m_pmesh->insert_face(v_list);
			ff->tag = true;
		}
	}

	//m_pmesh->update_mesh();
	return m_pmesh;
}
