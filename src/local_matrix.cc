#include <set>

#include "local_matrix.h"
#include "p4est/my_p4est_interface.h"
#include "shapefun.h"
#include "integration_cell.h"

using namespace std;

LocalMatrixComponent::LocalMatrixComponent(int ndofs) : ndofs(ndofs){
   mat.resize(ndofs, vector<real>(ndofs, 0.0));
}

LocalVectorComponent::LocalVectorComponent(int ndofs) : ndofs(ndofs){
   vec.resize(ndofs, 0.0);
}

LocalMatrix::LocalMatrix(int ncomponents, int nnodes) : ncomponents(ncomponents), nnodes(nnodes) {
   comps.resize(ncomponents, vector<LocalMatrixComponent>(ncomponents, LocalMatrixComponent(nnodes)));
}

LocalVector::LocalVector(int ncomponents, int nnodes) : ncomponents(ncomponents), nnodes(nnodes) {
   comps.resize(ncomponents, LocalVectorComponent(nnodes));
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

void LocalMatrixComponent::print() const {
   cout << PrintVec2D<double>(mat);
}

void LocalVectorComponent::print() const {
   cout << PrintVec<double>(vec) << endl;
}

void LocalMatrix::print() const {
   for(int i = 0; i < ncomponents; i++) {
      for(int j = 0; j < ncomponents; j++) {
         cout << "component (" << i << ", " << j << ")" << endl;
         comps[i][j].print();
      }
   }
   cout << endl;
}

void LocalVector::print() const {
   for(int i = 0; i < ncomponents; i++) {
      cout << "component (" << i << ")" << endl;
      comps[i].print();
   }
   cout << endl;
}


HangingInfo::HangingInfo(const P4estClass &p4est) : p4est(p4est){
   active_elem_idx = -1;
   faces.resize(Def::d()->num_faces, 0.0);
   edges.resize(Def::d()->num_edges, 0.0);

   coefs.resize(Def::d()->num_element_nodes, vector<double>(Def::d()->num_element_nodes, 0.0));
}


// coefs[local_dof][independent_dof]
// independent_dof is actually \sum_ld coefs[id][ld] * shape_ld

void HangingInfo::init_coefs(int elem_idx, const IntegrationCell &cell) {
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
   active_elem_idx = elem_idx;

   if(!active_elem_anyhang)
      return;

   ReferenceElement ref_elem(Def::d()->num_dim, Def::d()->order);
   std::set<int> hanging_nodes;

   for(int face_idx = 0; face_idx < Def::d()->num_faces; face_idx++) {
      if(faces[face_idx] != -1) {
         // it is hanging. find all its nodes
         for(int face_node : ref_elem.face_nodes[face_idx]) {
            hanging_nodes.insert(face_node);
         }
      }
   }

   for(int edge_idx = 0; edge_idx < Def::d()->num_edges; edge_idx++) {
      if(edges[edge_idx] != -1) {
         // it is hanging. find all its nodes
         for(int edge_node : ref_elem.edge_nodes[edge_idx]) {
            hanging_nodes.insert(edge_node);
         }
      }
   }

   for(int hanging_node : hanging_nodes) {
//      cout << PrintVec<double>(ref_elem.children_nodes_parent_basis_values[cell.child_position][hanging_node]) << endl;
      coefs[hanging_node] = ref_elem.children_nodes_parent_basis_values[cell.child_position][hanging_node];
   }

}

void HangingInfo::apply_constraints(int elem_idx, const IntegrationCell &cell,
                                    const LocalMatrixComponent &in, LocalMatrixComponent *out)
{
   init_coefs(elem_idx, cell);

   if(!active_elem_anyhang) {
      out->copy(in);
      return;
   }

   num_hanging_elements++;

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

void HangingInfo::apply_constraints(int elem_idx, const IntegrationCell &cell, bool inverse,
                                    const LocalVectorComponent &in, LocalVectorComponent *out)
{
   init_coefs(elem_idx, cell);

   if(!active_elem_anyhang) {
      out->copy(in);
      return;
   }

   out->clear();
   for(int i_loc = 0; i_loc < in.ndofs; i_loc++) {
      for(int i_glob = 0; i_glob < in.ndofs; i_glob++) {

         if(inverse)
            out->vec[i_glob] += in.vec[i_loc] * coefs[i_glob][i_loc];
         else
            out->vec[i_glob] += in.vec[i_loc] * coefs[i_loc][i_glob];
      }
   }
}

void HangingInfo::apply_constraints(int elem_idx, const IntegrationCell &cell,
                                    const LocalMatrix &in, LocalMatrix *out)
{
   assert(in.ncomponents == out->ncomponents);
   assert(in.nnodes == out->nnodes);
   for(int i_comp = 0; i_comp < in.ncomponents; i_comp++) {
      for(int j_comp = 0; j_comp < in.ncomponents; j_comp++) {
         apply_constraints(elem_idx, cell, in.comps[i_comp][j_comp], &out->comps[i_comp][j_comp]);
      }
   }
}

void HangingInfo::apply_constraints(int elem_idx, const IntegrationCell &cell,
                                    const LocalVector &in, LocalVector *out)
{
   assert(in.ncomponents == out->ncomponents);
   assert(in.nnodes == out->nnodes);
   for(int i_comp = 0; i_comp < in.ncomponents; i_comp++) {
      apply_constraints(elem_idx, cell, false, in.comps[i_comp], &out->comps[i_comp]);
   }
}

void HangingInfo::apply_constraints_inverse(int elem_idx, const IntegrationCell &cell,
                                    const LocalVector &in, LocalVector *out)
{
   assert(in.ncomponents == out->ncomponents);
   assert(in.nnodes == out->nnodes);
   for(int i_comp = 0; i_comp < in.ncomponents; i_comp++) {
      apply_constraints(elem_idx, cell, true, in.comps[i_comp], &out->comps[i_comp]);
   }
}


int HangingInfo::num_hanging_elements;







