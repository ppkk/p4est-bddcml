#ifndef INTEGRATION_CELL_H
#define INTEGRATION_CELL_H

#include "definitions.h"

class IntegrationCell
{
public:
   // creates empty element
   IntegrationCell() {}
   IntegrationCell(std::vector<double> coords, double size) : position(coords), size(size) {}
   // creates square
   IntegrationCell(double x, double y, double size);
   // creates cube
   IntegrationCell(double x, double y, double z, double size);

   void clear();

   void pokus();

   int n_dimensions() const  {assert(! position.empty()); return position.size();}
   std::vector<std::vector<double> > corners_coords() const;
   std::vector<std::vector<double> > nodes_coords(int num_nodes_1d) const;

   void fill_parent_cell(IntegrationCell *parent) const;

public:
   std::vector<double> position;
   double size;
   int refinement_level;
   int child_position;
};


class IntegrationMesh
{
public:
   void clear() {cells.clear(); }
   int num_elements() const {return cells.size(); }
public:
   std::vector<IntegrationCell> cells;
};


#endif // INTEGRATION_CELL_H
