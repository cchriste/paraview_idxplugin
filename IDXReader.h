#ifndef __IDXReader_h
#define __IDXReader_h
 
//#include "vtkPolyDataAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGrid.h"
#include "vtkRectilinearGridAlgorithm.h"
#include "vtkRectilinearGrid.h"
#include "vtkImageAlgorithm.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include <visuscpp/idx/visus_idx.h>

using std::string;
//VISUS_USE_NAMESPACE

class IDXReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  vtkTypeMacro(IDXReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  int DataTypeIndex(char* dname);
  int ReadFile(char* field_name, int* offset, int* count, int* r_offset, int* r_count, vtkRectilinearGrid* output);
  void SetTransformMode(int mode);
  void SetPercentToRemove(double percent);
  //AM void SetFieldName(char* field_name);
 
  static IDXReader *New();
 
  // Description:
  // Specify file name of the .idx file.
  vtkSetStringMacro(FileName);

  // Specify field name from .idx file.
  //AM vtkGetStringMacro(FieldName);

  //int GetTransformMode();
  //int SetPercentToRemove(int setValue);
protected:
  IDXReader();
  ~IDXReader(){}
 
  //
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
private:
  IDXReader(const IDXReader&);       // Not implemented.
  void operator=(const IDXReader&);  // Not implemented.
 
  double PercentToRemove;
  char* FileName;
  char* FieldName;
  int TransformMode;

  //ScopedPtr<Dataset> dataset;
  //ScopedPtr<Application> app;
  //ScopedPtr<Access> accessdata;

};
 
#endif
