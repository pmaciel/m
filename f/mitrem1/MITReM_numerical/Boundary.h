//---------------------------------------------------------------------------

#ifndef BoundaryMitremH
#define BoundaryMitremH

//---------------------------------------------------------------------------

#include <string>
#include "TypeDefs.h"

//---------------------------------------------------------------------------

class Boundary
{
public :
  Boundary(unsigned nBoundaryElements_);
  ~Boundary();

  void    setBoundaryElements(unsigned b, unsigned boundaryElement);
  void    setType(const std::string &type);

  unsigned  getNBoundaryElements() const;
  IndexList  getBoundaryElements() const;
  unsigned  getBoundaryElements(unsigned b) const;
  std::string  getType() const;

private :
  unsigned    nBoundaryElements;
  EmptyIndexList  boundaryElements;    // array of boundary element indices
  std::string    type;          // type of boundary
};

//---------------------------------------------------------------------------

inline void Boundary::setType(const std::string &type)
{
  this->type = type;
}

//---------------------------------------------------------------------------

inline unsigned Boundary::getNBoundaryElements() const
{
  return nBoundaryElements;
}

//---------------------------------------------------------------------------

inline IndexList Boundary::getBoundaryElements() const
{
  return boundaryElements;
}

//---------------------------------------------------------------------------

inline unsigned Boundary::getBoundaryElements(unsigned b) const
{
  return boundaryElements[b];
}

//---------------------------------------------------------------------------

inline std::string Boundary::getType() const
{
  return type;
}

//---------------------------------------------------------------------------

#endif
