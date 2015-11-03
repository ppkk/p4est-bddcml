#ifndef ELEMENT_H
#define ELEMENT_H

#include <assert.h>

#include "definitions.h"

class IntegrationCell;
class ReferenceElement;

// todo: this is not definitive. The idea is to separate IntegrationCell, which is only geometric entity
// used for integration.
// NodalElement should represent something, which than can be equipped with dofs, but not yet. It, however,
// contains nodes, which are GLOBALY numbered by p4est::lnodes and this numbering is used in BDDCML
// NodalElementComponent is 1 component of field (velocity, pressure, etc)
// The relationship of NodalElementComponent and ReferenceElement, however, has to be refined...

class NodalElementComponent
{
public:
   NodalElementComponent(const IntegrationCell &cell, const ReferenceElement &reference_element);
   int order() const;
   int ndofs() const {return dofs.size(); }
public:
   // nodal element corresponds to 1 integration cell
   const IntegrationCell &cell;
   const ReferenceElement &reference_element;

   std::vector<int> dofs;
};

class NodalElement
{
public:
   // constructor which creates Nodal element with ncomponents SAME components
   // there might be different constructors later (think of velocity + pressure for flow)
   NodalElement(int elem_idx, int ncomponents, const IntegrationCell &cell, const ReferenceElement &reference_element);
   int ncomponents() const {return components.size(); }

   void add_node_and_dofs(int subdomain_node_idx);

public:
   const IntegrationCell &cell;
   std::vector<NodalElementComponent> components;

   int elem_idx;
   std::vector<int> nodes;
};

class NodalElementMesh
{
public:
   //NodalElementMesh();
   void clear() {elements.clear(); }
   int num_elements() const {return elements.size(); }

public:
   std::vector<NodalElement> elements;
};

#endif // ELEMENT_H
