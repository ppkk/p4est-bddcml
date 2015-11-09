#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <string>

#include "definitions.h"

class NodalElementMesh;
class RealArray;

class VtkOutput
{
public:
   VtkOutput(const P4estClass &p4est, const NodalElementMesh &mesh, const std::vector<double> &solutions);

   void output_in_corners(const std::string &filename);
   void output_in_nodes(const std::string &filename);
   void output_exact_sol_in_nodes(const std::string &filename, exact_fn exact_sol);

private:
   void clear_arrays();
   void init_arrays(int num_displayed_elems_in_each_elem);
   void prepare_arrays_corners();
   void prepare_arrays_nodes(exact_fn exact);
   void output_pvtu(const std::string &filename);
   void output(const std::string &filename);

   int num_displayed_elements;

   std::vector<std::vector<double> > coords;
   std::vector<std::vector<int> > conectivity;
   std::vector<int> offsets;
   std::vector<int> types;
   std::vector<int> treeid;
   std::vector<int> level;
   std::vector<int> mpirank;
   std::vector<std::vector<double> > solutions_values;
   std::vector<std::string> solutions_names;

   const P4estClass &p4est;
   const NodalElementMesh &mesh;
   const std::vector<double> &solutions;
};

#endif // VTK_OUTPUT_H
