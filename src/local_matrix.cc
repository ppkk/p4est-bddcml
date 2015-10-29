#include "local_matrix.h"
#include "p4est/my_p4est_interface.h"

using namespace std;

LocalMatrixComponent::LocalMatrixComponent(int ndofs) : ndofs(ndofs){
   mat.resize(ndofs, vector<real>(ndofs, 0.0));
}

LocalVectorComponent::LocalVectorComponent(int ndofs) : ndofs(ndofs){
   vec.resize(ndofs, 0.0);
}

LocalMatrix::LocalMatrix(int ncomponents, int ndofs) : ncomponents(ncomponents), ndofs(ndofs) {
   comps.resize(ncomponents, vector<LocalMatrixComponent>(ncomponents, LocalMatrixComponent(ndofs)));
}

LocalVector::LocalVector(int ncomponents, int ndofs) : ncomponents(ncomponents), ndofs(ndofs) {
   comps.resize(ncomponents, LocalVectorComponent(ndofs));
}

void LocalMatrixComponent::copy(const LocalMatrixComponent &copy_from)
{
   for(int i = 0; i < ndofs; i++) {
      for(int j = 0; j < ndofs; j++) {
         this->mat[i][j] = copy_from.mat[i][j];
      }
   }
}

void LocalVectorComponent::copy(const LocalVectorComponent &copy_from)
{
   for(int i = 0; i < ndofs; i++) {
      this->vec[i] = copy_from.vec[i];
   }
}

void LocalMatrixComponent::clear() {
   for(int i_dof = 0; i_dof < ndofs; i_dof++) {
      for(int j_dof = 0; j_dof < ndofs; j_dof++) {
         mat[i_dof][j_dof] = 0.0;
      }
   }
}

void LocalVectorComponent::clear() {
   for(int i_dof = 0; i_dof < ndofs; i_dof++) {
      vec[i_dof] = 0.0;
   }
}

void LocalMatrix::clear()
{
   for(int i_comp = 0; i_comp < ncomponents; i_comp++) {
      for(int j_comp = 0; j_comp < ncomponents; j_comp++) {
         comps[i_comp][j_comp].clear();
      }
   }
}

void LocalVector::clear()
{
   for(int i_comp = 0; i_comp < ncomponents; i_comp++) {
      comps[i_comp].clear();
   }
}


HangingInfo::HangingInfo(const P4estClass &p4est) : p4est(p4est){
   active_elem_idx = -1;
   faces.resize(Def::d()->num_faces, 0.0);
   edges.resize(Def::d()->num_edges, 0.0);

   coefs.resize(Def::d()->num_element_nodes, vector<double>(Def::d()->num_element_nodes, 0.0));
}


// coefs[local_dof][independent_dof]
// independent_dof is actually \sum_ld coefs[id][ld] * shape_ld


void HangingInfo::init_coefs(int elem_idx)
{
   if(elem_idx == active_elem_idx)
      return;

   for(int i = 0; i < Def::d()->num_element_nodes; i++) {
      for(int j = 0; j < Def::d()->num_element_nodes; j++) {
         coefs[i][j] = (i == j) ? 1.0 : 0.0;
      }
   }

   std::fill(faces.begin(), faces.end(), 0.0);
   std::fill(edges.begin(), edges.end(), 0.0);
   active_elem_anyhang = p4est.get_hanging_info(elem_idx, this);

   if(!active_elem_anyhang)
      return;

   int indep_node, hanging_node;
   vector<int> indep_nodes;

   if(p4est.get_num_dim() == 2) {
      for(int face_idx = 0; face_idx < Def::d()->num_faces; face_idx++) {
         if(faces[face_idx] == -1)
            continue;

         if(faces[face_idx] == 0) {
            // it is a first part of a larger face
            indep_node = Def::d()->face_corners[face_idx][0];
            hanging_node = Def::d()->face_corners[face_idx][1];
         }
         else if (faces[face_idx] == 1) {
            // it is a first part of a larger face
            indep_node = Def::d()->face_corners[face_idx][1];
            hanging_node = Def::d()->face_corners[face_idx][0];
         }
         else {
            assert(0);
         }

         coefs[hanging_node][hanging_node] = 0.5;
         coefs[hanging_node][indep_node] = 0.5;
      }

   }
   else if(p4est.get_num_dim() == 3) {
      for(int edge_idx = 0; edge_idx < Def::d()->num_edges; edge_idx++) {
         if(edges[edge_idx] == -1)
            continue;

         if((edges[edge_idx] == 0) || (edges[edge_idx] == 2)) {
            indep_node = Def::d()->edge_corners[edge_idx][0];
            hanging_node = Def::d()->edge_corners[edge_idx][1];
         }
         else if((edges[edge_idx] == 1) || (edges[edge_idx] == 3)) {
            indep_node = Def::d()->edge_corners[edge_idx][1];
            hanging_node = Def::d()->edge_corners[edge_idx][0];
         }
         else {
            assert(edges[edge_idx] == 4); // inner edge, will be dealt with using faces
            continue;
         }

         coefs[hanging_node][hanging_node] = 0.5;
         coefs[hanging_node][indep_node] = 0.5;
      }
      for(int face_idx = 0; face_idx < Def::d()->num_faces; face_idx++) {
         if(faces[face_idx] == -1)
            continue;

         indep_nodes.clear();
         for(int loc_face_corner = 0; loc_face_corner < Def::d()->num_face_corners; loc_face_corner++) {
            // finding the hanging for this situation
            if(loc_face_corner == 3 - faces[face_idx]) {
               hanging_node = Def::d()->face_corners[face_idx][loc_face_corner];
            }
            else {
               indep_nodes.push_back(Def::d()->face_corners[face_idx][loc_face_corner]);
            }
         }

         coefs[hanging_node][hanging_node] = 0.25;
         for(int indep_node : indep_nodes) {
            coefs[hanging_node][indep_node] = 0.25;
         }
      }
   }
   else {
      assert(0);
   }
   active_elem_idx = elem_idx;
}

void HangingInfo::apply_constraints(int elem_idx, const LocalMatrixComponent &in, LocalMatrixComponent *out)
{
   init_coefs(elem_idx);

   if(!active_elem_anyhang) {
      out->copy(in);
      return;
   }

   out->clear();
   for(int i_loc = 0; i_loc < in.ndofs; i_loc++) {
      for(int i_glob = 0; i_glob < in.ndofs; i_glob++) {
         if(coefs[i_loc][i_glob] == 0.0)
            continue;

         for(int j_loc = 0; j_loc < in.ndofs; j_loc++) {
            for(int j_glob = 0; j_glob < in.ndofs; j_glob++) {
               if(coefs[j_loc][j_glob] == 0.0)
                  continue;

               out->mat[i_glob][j_glob] += in.mat[i_loc][j_loc] * coefs[i_loc][i_glob] * coefs[j_loc][j_glob];
            }
         }
      }
   }
}

void HangingInfo::apply_constraints(int elem_idx, const LocalVectorComponent &in, LocalVectorComponent *out)
{
   init_coefs(elem_idx);

   if(!active_elem_anyhang) {
      out->copy(in);
      return;
   }

   out->clear();
   for(int i_loc = 0; i_loc < in.ndofs; i_loc++) {
      for(int i_glob = 0; i_glob < in.ndofs; i_glob++) {
         if(coefs[i_loc][i_glob] == 0.0)
            continue;

         out->vec[i_glob] += in.vec[i_loc] * coefs[i_loc][i_glob];
      }
   }
}

void HangingInfo::apply_constraints(int elem_idx, const LocalMatrix &in, LocalMatrix *out)
{
   assert(in.ncomponents == out->ncomponents);
   assert(in.ndofs == out->ndofs);
   for(int i_comp = 0; i_comp < in.ncomponents; i_comp++) {
      for(int j_comp = 0; j_comp < in.ncomponents; j_comp++) {
         apply_constraints(elem_idx, in.comps[i_comp][j_comp], &out->comps[i_comp][j_comp]);
      }
   }
}

void HangingInfo::apply_constraints(int elem_idx, const LocalVector &in, LocalVector *out)
{
   assert(in.ncomponents == out->ncomponents);
   assert(in.ndofs == out->ndofs);
   for(int i_comp = 0; i_comp < in.ncomponents; i_comp++) {
      apply_constraints(elem_idx, in.comps[i_comp], &out->comps[i_comp]);
   }
}









