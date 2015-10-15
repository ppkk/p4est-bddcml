#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

class Element
{
public:
   int n_dimensions() const  {return position.size();}
   std::vector<std::vector<double> > nodes() const;

public:
   std::vector<double> position;
   double size;
};

#endif // ELEMENT_H
