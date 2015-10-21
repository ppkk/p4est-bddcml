#ifndef GEOMETRY_MESH_H
#define GEOMETRY_MESH_H

#include "element.h"

struct p4est;
struct p4est_lnodes;

class GeometryMesh
{
public:
   void clear() {elements.clear(); }
   int num_elements() const {return elements.size(); }
public:
   std::vector<Element> elements;
};


#endif // GEOMETRY_MESH_H
