//---------------------------------------------------------------------------
#include "DatFileReader.h"
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DatFileReader::DatFileReader(const std::string &filename) //constructor
{
  filename_ = filename;

  // - data inlezen en voor verdere verwerking in array van std::string steken -

  std::ifstream istrm(filename_.c_str());
  if (istrm == NULL) errorFileDoesNotExist(filename_);
  istrm.setf(std::ios::scientific);

  std::string temp;
  getline(istrm,temp);  // eerste regel (versie-info) overslaan

  DataElement DataElement_;
  char c;
  while( !istrm.eof() ) {
    c = istrm.get();

    // whitespace skippen
    while ((c == int('\n'))||(c == int(' '))||(c == int('\t')))
      c = istrm.get();

    // zien of er commentaar in de dat-file zit
    if (c == int('/')) {
      c = istrm.get();
      if (c == int('/')) { // is dit een lijn commentaar?
        getline(istrm,temp); //commentaarlijn binnenhalen
        //std::cout << "Commentaar tegengekomen in de dat-file : " << "\"" << temp << "\"" << std::endl;
      }
      if (c == int('*')) { // is dit een blok commentaar?
        //std::cout << "Blok commentaar tegengekomen in de dat-file : " << "\"" ;
        c = istrm.get();
        while (!(( c == int('*'))&&(istrm.peek()==int('/')))) {
        //std::cout << c;
        c = istrm.get();
        }
        //std::cout << "\"" << std::endl;
      }
    }

    // zien of we er een std::string in zijn geheel moeten uitpakken
    else if (c == int('"')) {
      getline(istrm,temp,'"'); // pointer wordt ineens achter de laatste " gezet (handig he)
      //std::cout << "STRING aangetroffen! : " << temp << std::endl;
      DataElement_.content = temp;
      DataElement_.used = 0;
      Data.push_back(DataElement_);
      }
      // "=" negeren

    else if (c == int('=')) {
      // we doen niks, en laten de = dus gewoon verdwijnen
    }

    // data in array steken
    else {
      istrm.unget();
      istrm >> temp;

      DataElement_.content = temp;
      DataElement_.used = 0;
      Data.push_back(DataElement_);
    }
  }
  istrm.close();
}
//---------------------------------------------------------------------------
std::string DatFileReader::getFileName()
{
  return filename_;
}
//---------------------------------------------------------------------------
std::string DatFileReader::getVersion()
{
  std::ifstream istrm(filename_.c_str());
  istrm.setf(std::ios::scientific);
  std::string FileString;
  getline(istrm,FileString);
  istrm.close();

  return FileString;
}
//---------------------------------------------------------------------------
void DatFileReader::readScalar(const std::string &parameter, std::string &returnval)
{
  bool found = false;
  for (std::vector<DataElement>::iterator it = Data.begin(); it != Data.end(); it++) {
    if ((*it).content == parameter) {
      (*it).used = 1;
      it++;
      (*it).used = 1;
      returnval = (*it).content;
      found = true;
    }
  }

  if (!found) {
    errorParameterDoesNotExist(filename_,parameter);
  }
}
//---------------------------------------------------------------------------
void DatFileReader::readScalar(const std::string &parameter, int& returnval)
{
  std::string string_;
  readScalar(parameter, string_);

  returnval = atoi(string_.c_str());
}
//---------------------------------------------------------------------------
void DatFileReader::readScalar(const std::string &parameter, unsigned& returnval)
{
  std::string string_;
  readScalar(parameter, string_);

  returnval = atoi(string_.c_str());
  if (atoi(string_.c_str()) < 0) errorWrongType(filename_,"unsigned");
}
//---------------------------------------------------------------------------
void DatFileReader::readScalar(const std::string &parameter, double& returnval)
{
  std::string string_;
  readScalar(parameter,string_);

  char * pEnd;
  returnval = strtod(string_.c_str(),&pEnd);
}
//---------------------------------------------------------------------------
void DatFileReader::readScalar(const std::string &parameter, bool& returnval)
{
  std::string string_;
  readScalar(parameter,string_);

  if (string_ == "false") returnval = false;
  else if (string_ == "true") returnval = true;
  else errorWrongType(filename_,"bool");
}
//---------------------------------------------------------------------------
void DatFileReader::readvector()
{
  std::cout << "WARNING: readvector not implemented yet!";
}
//---------------------------------------------------------------------------
void DatFileReader::readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, std::string &returnval)
{
  // ______ zou wel wa overzichtelijker gecodeerd kunnen worden in de toekomst __________

  bool foundname = false;
  unsigned foundparameters = 0;

  // eerst op zoek naar de naam van de vectoren, bvb [points]
  std::vector<DataElement>::iterator it;
  for (it = Data.begin(); it != Data.end(); it++) {
    if ((*it).content == name) {
      (*it).used = 1;
      foundname = true;
      break;
    }
  }

  if (!foundname) {
    errorParameterDoesNotExist(filename_,name);
  }

  // daarna op zoek naar de naam van de parameter, bvb <x>
  std::vector<DataElement>::iterator last_it = it + 1;
  for (unsigned n=0 ; n<=number ; n++) {
    for (it = last_it; it != Data.end(); it++) {
      if ((*it).content == parameter) {
        last_it = it + 1;
        foundparameters++;
        break;
      }
    }
  }

  if (foundparameters == 0) {
    errorVectorcompDoesNotExist(filename_,parameter);
  }
  else if (foundparameters < number+1){
    errorOutOfBounds(filename_,name);
  }
  else if (foundparameters == number+1){
    (*--last_it).used = 1;
    last_it++;
    (*last_it).used = 1;
    returnval = (*last_it).content;
    return;
  }
  else exit(1); // impossible error
}
//---------------------------------------------------------------------------
void DatFileReader::readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, int& returnval)
{
  std::string string_;
  readMultipleVector(name,number,parameter,string_);

  returnval = atoi(string_.c_str());
}
//---------------------------------------------------------------------------
void DatFileReader::readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, unsigned& returnval)
{
  std::string string_;
  readMultipleVector(name,number,parameter,string_);

  returnval = atoi(string_.c_str());
  if (atoi(string_.c_str()) < 0) errorWrongType(filename_,"unsigned");
}
//---------------------------------------------------------------------------
void DatFileReader::readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, double& returnval)
{
  std::string string_;
  readMultipleVector(name,number,parameter,string_);

  char * pEnd;
  returnval = strtod(string_.c_str(),&pEnd);
}
//---------------------------------------------------------------------------
void DatFileReader::readMultipleVector(const std::string &name, unsigned number, const std::string &parameter, bool& returnval)
{
  std::string string_;
  readMultipleVector(name,number,parameter,string_);

  if (string_ == "false") returnval = false;
  else if (string_ == "true") returnval = true;
  else errorWrongType(filename_,"bool");
}
//---------------------------------------------------------------------------
unsigned DatFileReader::readMultipleVector_nVectors(const std::string &name)
{
  for (std::vector<DataElement>::iterator it = Data.begin(); it != Data.end(); it++) {
    if ((*it).content == name) {
      (*it).used = 1;   // naam gevonden
      it++;
      (*it).used = 1;   // inhoud
      return atoi((*it).content.c_str()); // de std::string typecasten naar int
    }
  }

  errorParameterDoesNotExist(filename_,name);
  return 0;
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void DatFileReader::errorFileDoesNotExist(const std::string &file)
{
  std::cout << "ERROR IN DatFileReader.cpp.\
      \nTHE FILE " << file <<" DOES NOT EXIST." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------
void DatFileReader::errorParameterDoesNotExist(const std::string &file, const std::string &parameter)
{
  std::cout << "ERROR IN " << file << ".\
      \nTHE PARAMETER " << parameter <<" DOES NOT EXIST.\
      \nMAKE SURE YOU SPELLED THE NAME CORRECTLY AND PUT IT BETWEEN []." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------
void DatFileReader::errorVectorcompDoesNotExist(const std::string &file, const std::string &vectorcomp)
{
  std::cout << "ERROR IN " << file << ".\
      \nTHE VECTOR COMPONENT " << vectorcomp <<" DOES NOT EXIST.\
      \nMAKE SURE YOU SPELLED THE NAME CORRECTLY AND PUT IT BETWEEN <>." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------
void DatFileReader::errorOutOfBounds(const std::string &file, const std::string &multivector)
{
  std::cout << "ERROR IN " << file << ".\
      \nYOU ARE TRYING TO ACCESS A VECTOR OF " << multivector <<" THAT DOES NOT EXIST.\
      \n\nPOSSIBLE CAUSES:\
      \n- YOU DEFINED LESS VECTORS FOR IT THAN YOU DECLARED,\
      \n- YOU MISSPELLED ONE OF THE VECTOR COMPONENTS." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------
void DatFileReader::errorWrongType(const std::string &file, const std::string &type)
{
  std::cout << "ERROR IN " << file << ".\
      \nTHE VALUE IS NOT OF TYPE " << type << "." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------
