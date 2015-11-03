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

VtkOutput::VtkOutput(const P4estClass &p4est, const NodalElementMesh &mesh, const RealArray &solutions) :
            p4est(p4est), mesh(mesh), solutions(solutions) {
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
   clear_arrays();
   prepeare_arrays();

   if(mpi_rank == 0)
      output_pvtu(filename);

   stringstream str;
   str << filename << "_" << mpi_rank << ".vtu";
   ofstream f(str.str());

   f << "<?xml version=\"1.0\"?>" << endl;
   f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
   f << "  <UnstructuredGrid>" << endl;
   f << "    <Piece NumberOfPoints=\"" << mesh.num_elements() * Def::d()->num_corners
     << "\" NumberOfCells=\"" << mesh.num_elements() <<"\">" << endl;
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

void VtkOutput::prepeare_arrays() {
   init_arrays();

   int corner_id = 0;

   ReferenceElement ref_elem(Def::d()->num_dim, Def::d()->order);

   for(int comp = 0; comp < Def::d()->num_components; comp++) {
      stringstream str;
      str << "solution_" << comp+1;
      solutions_names.push_back(str.str());
   }


   for(const NodalElement& elem : mesh.elements) {
      vector<int> elem_conectivity;
      for(auto corner_coords : elem.cell.corners_coords()) {
         if(corner_coords.size() == 2)
            corner_coords.push_back(0.0);
         coords.push_back(corner_coords);
         elem_conectivity.push_back(corner_id++);
      }
      conectivity.push_back(elem_conectivity);
      if(offsets.empty())
         offsets.push_back(Def::d()->num_corners);
      else
         offsets.push_back(offsets.back() + Def::d()->num_corners);
      treeid.push_back(0);
      types.push_back((Def::d()->num_dim == 2) ? 8 : 11);
      level.push_back(elem.cell.refinement_level);
      mpirank.push_back(mpi_rank);

      LocalSolution loc_sol(p4est, elem, ref_elem, solutions.val);

      for(int comp = 0; comp < Def::d()->num_components; comp++) {
         for(int corner_node : ref_elem.corner_nodes) {
            solutions_values[comp].push_back(loc_sol.loc_vec.comps[comp].vec[corner_node]);
         }

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

void VtkOutput::init_arrays() {
   solutions_values = vector<vector<double> >(Def::d()->num_components);
}
