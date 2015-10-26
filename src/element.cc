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


NodalElement::NodalElement(int ncomponents, const IntegrationCell &cell,
                           const ReferenceElement &reference_element) :
         cell(cell) {
   for(int i = 0; i < ncomponents; i++) {
      components.push_back(NodalElementComponent(cell, reference_element));
   }
}
