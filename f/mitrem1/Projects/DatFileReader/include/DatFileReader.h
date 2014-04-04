//---------------------------------------------------------------------------

#ifndef DatFileReaderH
#define DatFileReaderH

//---------------------------------------------------------------------------

#include <cstdlib>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//---------------------------------------------------------------------------

class DatFileReader
{
public :
  DatFileReader(const std::string &filename);

  std::string getFileName();
  std::string getVersion();
  void printData();

  void readScalar(const std::string &parameter, int& returnval);
  void readScalar(const std::string &parameter, unsigned& returnval);
  void readScalar(const std::string &parameter, double& returnval);
  void readScalar(const std::string &parameter, bool& returnval);
  void readScalar(const std::string &parameter, std::string &returnval);

  void readvector(); // not yet implemented

  unsigned readMultipleVector_nVectors(const std::string &name);
  void readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, int& returnval);
  void readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, unsigned& returnval);
  void readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, double& returnval);
  void readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, bool& returnval);
  void readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, std::string &returnval);

  void checkusage();

private :
  void errorFileDoesNotExist(const std::string &file);
  void errorParameterDoesNotExist(const std::string &file, const std::string &parameter);
  void errorVectorcompDoesNotExist(const std::string &file, const std::string &vectorcomp);
  void errorOutOfBounds(const std::string &file, const std::string &multivector);
  void errorWrongType(const std::string &file, const std::string &type);

  std::string filename_;

  typedef struct {
    std::string content;
    bool used;
  } DataElement;

  std::vector<DataElement> Data;
};

//---------------------------------------------------------------------------
#endif
