#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

class Element
{
public:
   // creates empty element
   Element() {}
   Element(std::vector<double> coords, double size) : position(coords), size(size) {}
   // creates square
   Element(double x, double y, double size);
   // creates cube
   Element(double x, double y, double z, double size);
   int n_dimensions() const  {return position.size();}
   std::vector<std::vector<double> > nodes() const;

public:
   std::vector<double> position;
   double size;
};

#endif // ELEMENT_H
