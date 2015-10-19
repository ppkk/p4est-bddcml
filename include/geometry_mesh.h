#ifndef GEOMETRY_MESH_H
#define GEOMETRY_MESH_H

#include "element.h"

struct p4est;
struct p4est_lnodes;

class GeometryMesh
{
public:
   void prepare_subdomain_mesh(p4est *p4est, p4est_lnodes *lnodes);

public:
   std::vector<Element> elements;
};


#endif // GEOMETRY_MESH_H
