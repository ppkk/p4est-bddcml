#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "vtk_output.h"
#include "element.h"
#include "local_solution.h"
#include "integration_cell.h"
#include "shapefun.h"
#include "bddcml/arrays.h"

using namespace std;

VtkOutput::VtkOutput(const P4estClass &p4est, const NodalElementMesh &mesh, const std::vector<double> &solutions) :
            p4est(p4est), mesh(mesh), solutions(solutions) {
}

void VtkOutput::output_in_corners(const string &filename) {
   clear_arrays();
   init_arrays(1);
   prepare_arrays_corners();
   output(filename);
}

void VtkOutput::output_in_nodes(const string &filename) {
   clear_arrays();
   int displayed_in_elem = Def::d()->order * Def::d()->order;
   if(Def::d()->num_dim == 3)
      displayed_in_elem *= Def::d()->order;

   init_arrays(displayed_in_elem);
   prepare_arrays_nodes(nullptr);
   output(filename);
}

void VtkOutput::output_exact_sol_in_nodes(const string &filename, exact_fn exact_sol) {
   clear_arrays();
   int displayed_in_elem = Def::d()->order * Def::d()->order;
   if(Def::d()->num_dim == 3)
      displayed_in_elem *= Def::d()->order;

   init_arrays(displayed_in_elem);
   prepare_arrays_nodes(exact_sol);
   output(filename);

}

void VtkOutput::output_pvtu(const string &filename) {
   stringstream str;
   str << filename << ".pvtu";
   ofstream f(str.str());

   f << "<?xml version=\"1.0\"?>" << endl;
   f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
   f << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl;
   f << "    <PPoints>" << endl;
   f << "      <PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"/>" << endl;
   f << "    </PPoints>" << endl;
   f << "    <PCellData Scalars=\"treeid,level,mpirank\">" << endl;
   f << "      <PDataArray type=\"Int32\" Name=\"treeid\" format=\"ascii\"/>" << endl;
   f << "      <PDataArray type=\"UInt8\" Name=\"level\" format=\"ascii\"/>" << endl;
   f << "      <PDataArray type=\"Int32\" Name=\"mpirank\" format=\"ascii\"/>" << endl;
   f << "    </PCellData>" << endl;
   f << "    <PPointData>" << endl;
   for(auto solname : solutions_names)
      f << "      <PDataArray type=\"Float32\" Name=\"" << solname <<"\" format=\"ascii\"/>" << endl;
   f << "    </PPointData>" << endl;

   for(int rank = 0; rank < mpi_size; rank++) {
      f << "    <Piece Source=\"" << filename << "_" << rank << ".vtu" << "\"/>" << endl;
   }
   f << "  </PUnstructuredGrid>" << endl;
   f << "</VTKFile>" << endl;

   f.close();

   str.str("");
   str.clear();
   str << filename << ".visit";
   f.open(str.str());
   f << "!NBLOCKS " << mpi_size << endl;
   for(int rank = 0; rank < mpi_size; rank++) {
      f << filename << "_" << rank << ".vtu" << endl;
   }

   f.close();
}

void VtkOutput::output(const string &filename) {
   if(mpi_rank == 0)
      output_pvtu(filename);

   stringstream str;
   str << filename << "_" << mpi_rank << ".vtu";
   ofstream f(str.str());

   f << "<?xml version=\"1.0\"?>" << endl;
   f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
   f << "  <UnstructuredGrid>" << endl;
   f << "    <Piece NumberOfPoints=\"" << num_displayed_elements * Def::d()->num_corners
     << "\" NumberOfCells=\"" << num_displayed_elements <<"\">" << endl;
   f << "      <Points>" << endl;
   f << "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
   f << PrintVec2D<double>(coords, " ", "          ");
   f << "        </DataArray>" << endl;
   f << "      </Points>" << endl;
   f << "      <Cells>" << endl;
   f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
   f << PrintVec2D<int>(conectivity, " ", "          ");
   f << "        </DataArray>" << endl;
   f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
   f << "          " << PrintVec<int>(offsets, " ") << endl;
   f << "        </DataArray>" << endl;
   f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
   f << "          " << PrintVec<int>(types, " ") << endl;
   f << "        </DataArray>" << endl;
   f << "      </Cells>" << endl;
   f << "      <CellData Scalars=\"treeid,level,mpirank\">" << endl;
   f << "        <DataArray type=\"Int32\" Name=\"treeid\" format=\"ascii\">" << endl;
   f << "          " << PrintVec<int>(treeid, " ") << endl;
   f << "        </DataArray>" << endl;
   f << "        <DataArray type=\"UInt8\" Name=\"level\" format=\"ascii\">" << endl;
   f << "          " << PrintVec<int>(level, " ") << endl;
   f << "        </DataArray>" << endl;
   f << "        <DataArray type=\"Int32\" Name=\"mpirank\" format=\"ascii\">" << endl;
   f << "          " << PrintVec<int>(mpirank, " ") << endl;
   f << "        </DataArray>" << endl;
   f << "      </CellData>" << endl;
   f << "      <PointData Scalars=\"";
   for(auto solname : solutions_names)
      f << solname << ((solname == solutions_names.back()) ? "" : ",");
   f << "\">" << endl;
   for(unsigned component = 0; component < solutions_names.size(); component++) {
      f << "        <DataArray type=\"Float32\" Name=\"" << solutions_names[component] << "\" format=\"ascii\">" << endl;
      f << PrintVec<double>(solutions_values[component], "\n", "          ");
      f << "        </DataArray>" << endl;
   }

   f << "      </PointData>" << endl;
   f << "    </Piece>" << endl;
   f << "  </UnstructuredGrid>" << endl;
   f << "</VTKFile>" << endl;

   f.close();
}

void VtkOutput::prepare_arrays_corners() {
   ReferenceElement ref_elem(Def::d()->num_dim, Def::d()->order);
   vector<vector<double> > elem_corners_coords;

   for(const NodalElement& elem : mesh.elements) {
      // coordinates
      elem.cell.corners_coords(&elem_corners_coords);
      for(auto corner_coords : elem_corners_coords) {
         if(corner_coords.size() == 2)
            corner_coords.push_back(0.0);
         coords.push_back(corner_coords);
      }

      // todo: pass vector<>
      LocalSolution loc_sol(p4est, elem, ref_elem, solutions);
      // solution values
      for(int comp = 0; comp < Def::d()->num_components; comp++) {
         for(int corner_node : ref_elem.corner_nodes) {
            solutions_values[comp].push_back(loc_sol.loc_vec.comps[comp].vec[corner_node]);
         }
      }
   }
}

void VtkOutput::prepare_arrays_nodes(exact_fn exact) {
   ReferenceElement ref_elem(Def::d()->num_dim, Def::d()->order);
   vector<vector<double> > elem_nodes_coords;

   for(const NodalElement& elem : mesh.elements) {
      elem.cell.nodes_coords(Def::d()->order + 1, &elem_nodes_coords);
      LocalSolution loc_sol(p4est, elem, ref_elem, solutions);

      for(auto subelement_start : Def::d()->cartesian_ids_plot_subelements) {
         for(auto corner_difs : Def::d()->cartesian_ids_corners) {
            // transfer node id from pair/triplet to node_idx
            int node_idx = 0;
            for(int dim = Def::d()->num_dim - 1; dim >= 0; dim--) {
               node_idx *= (Def::d()->order + 1);
               node_idx += subelement_start[dim] + corner_difs[dim];
            }
            // coordinates
            vector<double> node_coords = elem_nodes_coords[node_idx];
            if(node_coords.size() == 2)
               node_coords.push_back(0.0);
            coords.push_back(node_coords);

            // solution values
            for(int comp = 0; comp < Def::d()->num_components; comp++) {
               if(exact == nullptr) {
                  solutions_values[comp].push_back(loc_sol.loc_vec.comps[comp].vec[node_idx]);
               }
               else {
                  solutions_values[comp].push_back(exact(node_coords)[comp]);
               }
            }

         }
      }
   }
}

void VtkOutput::init_arrays(int num_displayed_elems_in_each_elem) {
   num_displayed_elements = mesh.num_elements() * num_displayed_elems_in_each_elem;
   solutions_values.resize(Def::d()->num_components);
   for(int comp = 0; comp < Def::d()->num_components; comp++) {
      stringstream str;
      str << "solution_" << comp+1;
      solutions_names.push_back(str.str());
      solutions_values[comp].reserve(Def::d()->num_corners * num_displayed_elements);
   }

   coords.reserve(Def::d()->num_corners * num_displayed_elements);

   treeid.resize(num_displayed_elements, 0);
   types.resize(num_displayed_elements, (Def::d()->num_dim == 2) ? 8 : 11);
   mpirank.resize(num_displayed_elements, mpi_rank);

   conectivity.reserve(num_displayed_elements);
   offsets.reserve(num_displayed_elements);
   level.reserve(num_displayed_elements);

   int corner_id = 0;
   for(int elem = 0; elem < mesh.num_elements(); elem++) {
      for(int subelem_num = 0; subelem_num < num_displayed_elems_in_each_elem; subelem_num++) {
         vector<int> elem_conectivity;
         for(int corner = 0; corner < Def::d()->num_corners; corner++) {
            elem_conectivity.push_back(corner_id++);
         }
         conectivity.push_back(elem_conectivity);
         if(offsets.empty())
            offsets.push_back(Def::d()->num_corners);
         else
            offsets.push_back(offsets.back() + Def::d()->num_corners);
         level.push_back(mesh.elements[elem].cell.refinement_level);
      }
   }
}


void VtkOutput::clear_arrays() {
   coords.clear();
   conectivity.clear();
   offsets.clear();
   types.clear();
   treeid.clear();
   level.clear();
   mpirank.clear();
   solutions_values.clear();
   solutions_names.clear();
}
