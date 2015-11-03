#include "element.h"
#include "integration_cell.h"
#include "shapefun.h"

NodalElementComponent::NodalElementComponent(const IntegrationCell &cell,
                                             const ReferenceElement &reference_element) :
         cell(cell), reference_element(reference_element) {

}

int NodalElementComponent::order() const {
   return reference_element.order;
}


NodalElement::NodalElement(int elem_idx, int ncomponents, const IntegrationCell &cell,
                           const ReferenceElement &reference_element) :
         cell(cell), elem_idx(elem_idx) {
   for(int i = 0; i < ncomponents; i++) {
      components.push_back(NodalElementComponent(cell, reference_element));
   }
}

void NodalElement::add_node_and_dofs(int subdomain_node_idx) {
   nodes.push_back(subdomain_node_idx);
   for(int comp = 0; comp < Def::d()->num_components; comp++) {
      components[comp].dofs.push_back(Def::d()->num_components * subdomain_node_idx + comp);
   }
}
