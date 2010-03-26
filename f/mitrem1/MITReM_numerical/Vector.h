//---------------------------------------------------------------------------

#ifndef VectorH
#define VectorH

//---------------------------------------------------------------------------

#include "TypeDefs.h"

//---------------------------------------------------------------------------

class Vector
{
public :
  Vector(unsigned nComponents_);
  ~Vector();

  void      setComponents(unsigned c, double value);

  unsigned    getNComponents() const;
  DoubleVector  getComponents() const;
  double      getComponents(unsigned c) const;

private :
  unsigned      nComponents;
  EmptyDoubleVector  components;        // array of components of the vector
};

//---------------------------------------------------------------------------

inline unsigned Vector::getNComponents() const
{
  return nComponents;
}

//---------------------------------------------------------------------------

inline DoubleVector Vector::getComponents() const
{
  return components;
}

//---------------------------------------------------------------------------

inline double Vector::getComponents(unsigned c) const
{
  return components[c];
}

//---------------------------------------------------------------------------

#endif
