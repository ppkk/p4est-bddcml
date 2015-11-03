#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <string>

#include "definitions.h"

class NodalElementMesh;
class RealArray;

class VtkOutput
{
public:
   VtkOutput(const P4estClass &p4est, const NodalElementMesh &mesh, const RealArray &solutions);

   void output(const std::string &filename);

private:
   void clear_arrays();
   void init_arrays();
   void prepeare_arrays();
   void output_pvtu(const std::string &filename);

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
   const RealArray &solutions;
};

#endif // VTK_OUTPUT_H
