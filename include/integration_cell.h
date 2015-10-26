#ifndef INTEGRATION_CELL_H
#define INTEGRATION_CELL_H

#include <assert.h>
#include <vector>

//struct p4est;
//struct p4est_lnodes;

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

   int n_dimensions() const  {assert(! position.empty()); return position.size();}
   std::vector<std::vector<double> > corners() const;

public:
   std::vector<double> position;
   double size;
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
