#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include <vector>

class IntegrationCell;

typedef double (*LevelSetFunctionPtr)(std::vector<double> point);

enum class LevelSetValue
{
   Inside,
   Outside,
   Border,
   Undefined
};

class LevelSet
{
public:
   LevelSet(LevelSetFunctionPtr function) : fn(function) {}
   inline double operator() (std::vector<double> point) const { return fn(point); }

   LevelSetValue apply(const std::vector<double> point) const;

   // it is done in element nodes only!
   LevelSetValue apply(const IntegrationCell* cell) const;

private:
   LevelSetFunctionPtr fn;
};

#endif // LEVEL_SET_H
