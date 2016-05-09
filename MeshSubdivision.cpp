#include "StdAfx.h"
#include ".\meshsubdivision.h"


//--------------------------------------------------------------
#define PI 3.1415926536
#include <vector>

Mesh3D* MeshSubdivision::Doo_Sabin()
{
	if (m_pmesh == NULL) {
		return NULL;
	}

	Mesh3D* New_mesh = new Mesh3D;


	////////////////////////////////////////////////

	//your implmentation
	//New_mesh:  subdivisioned mesh
	//m_pmesh:  original mesh

	// insert all vertex

	//==================
	PTR_FACE_LIST face_list = m_pmesh->get_faces_list();
	FACE_ITER     fiter = face_list->begin();

	std::vector<int>  vert_id_list;
	vert_id_list.assign(m_pmesh->get_num_of_edges_list(), 0);//将所有点id赋值为0
	// count of vertex

	//=================


	// for each face, compute new vertex, insert vertex
	// and insert F-faces 
	// save the id of these newly inserted vertex , so they can
	// be acquired using the face and vertex

	//=================
	//int  dd = New_mesh->get_num_of_vertices_list();
	for (; fiter != face_list->end(); fiter++)
	{
		//对每一个面
		int  n = (*fiter)->valence;//边数
		HE_edge * anchor_edge = (*fiter)->edge;// as anchor edge
		HE_edge * temp_edge;
		VERTEX_LIST  v_list;

		float new_x, new_y, new_z;

		do{
			temp_edge = anchor_edge;
			// curretnt vertex is anchor_edge->vert

			//计算一个新顶点的坐标 new_x, new_y, new_z
// 			if (!m_pmesh->isborder(temp_edge->vert)){
				for (int i = 0; i < n; i++, temp_edge = temp_edge->next){
					if (temp_edge == anchor_edge){
						new_x = temp_edge->vert->x*(n + 5) / 4 / n;
						new_y = temp_edge->vert->y*(n + 5) / 4 / n;
						new_z = temp_edge->vert->z*(n + 5) / 4 / n;
					}
					else{
						new_x += temp_edge->vert->x*(3 + 2 * cos(2 * PI*i / n)) / 4 / n;
						new_y += temp_edge->vert->y*(3 + 2 * cos(2 * PI*i / n)) / 4 / n;
						new_z += temp_edge->vert->z*(3 + 2 * cos(2 * PI*i / n)) / 4 / n;
					}
				}
// 			}
// 			else
// 			{
// 				new_x = temp_edge->vert->x;
// 				new_y = temp_edge->vert->y;
// 				new_z = temp_edge->vert->z;
// 			}


			HE_vert *vv = New_mesh->insert_vertex(new_x, new_y, new_z);
			// insert a new vertex
			//for every half_edge, save the id of the corresponding vertex
			vert_id_list[anchor_edge->id] = vv->id;

			v_list.push_back(vv);

			anchor_edge = anchor_edge->next;

		} while (anchor_edge != (*fiter)->edge);


		// one face is ok, 
		New_mesh->insert_face(v_list);// F-faces

	}

	//=================


	//
	//search all vertex, construct V-faces 
	//

	//===============
	PTR_VERTEX_LIST  vert_list = m_pmesh->get_vertices_list();
	VERTEX_ITER       viter = vert_list->begin();
	for (; viter != vert_list->end(); viter++)
	{// for one vertex

		
		VERTEX_LIST   v_list;
		HE_edge * edge = (*viter)->edge;
		HE_edge * t_edge = NULL;  HE_vert * vert = NULL;

		if (m_pmesh->isborder(*viter))
			continue;

		do
		{
			t_edge = edge;
			do { t_edge = t_edge->next; } while (t_edge->next != edge);
			//找到上一条边
			vert = New_mesh->get_vertex(vert_id_list[t_edge->id]);
			v_list.push_back(vert);
			edge = t_edge->pair;
		} while (edge != (*viter)->edge);

		New_mesh->insert_face(v_list);

	}

	//===============


	//
	// search all edges, to build E-faces
	//

	//======================
	PTR_EDGE_LIST edge_list = m_pmesh->get_edges_list();
	EDGE_ITER     eiter = edge_list->begin();
	for (; eiter != edge_list->end(); eiter++){
		// note if the pair of this edge has been visited
		// the E-face is generated already
		if (m_pmesh->isborder(*eiter)||(*eiter)->pair->tag) continue;   else (*eiter)->tag = true;

		VERTEX_LIST   v_list;
		HE_vert * vert;   HE_edge * edge = (*eiter);

		vert = New_mesh->get_vertex(vert_id_list[edge->id]);
		v_list.push_back(vert);

		do{ edge = edge->next; } while (edge->next != (*eiter));
		vert = New_mesh->get_vertex(vert_id_list[edge->id]);
		v_list.push_back(vert);
		
		edge = edge->next->pair;
		vert = New_mesh->get_vertex(vert_id_list[edge->id]);
		v_list.push_back(vert);
		
		do{ edge = edge->next; } while (edge->next != (*eiter)->pair);
		vert = New_mesh->get_vertex(vert_id_list[edge->id]);
		v_list.push_back(vert);
		
		New_mesh->insert_face(v_list);
	}

	//======================


	/////////////////////////////////////////////////

	New_mesh->update_mesh();

	delete m_pmesh;
	m_pmesh = NULL;
	return New_mesh;


}

//--------------------------------------------------------------

Mesh3D* MeshSubdivision::Catmull_Clark()
{

	if (m_pmesh == NULL) {
		return NULL;
	}


}
