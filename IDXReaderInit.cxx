#include "vtkClientServerInterpreter.h"

#ifndef PARAVIEW_BUILD_SHARED_LIBS
#define PARAVIEW_BUILD_SHARED_LIBS
#endif
#if defined(PARAVIEW_BUILD_SHARED_LIBS) && defined(_WIN32)
# define VTK_WRAP_CS_EXPORT __declspec(dllexport)
#else
# define VTK_WRAP_CS_EXPORT
#endif

extern void IDXReader_Init(vtkClientServerInterpreter* csi);


extern "C" void VTK_WRAP_CS_EXPORT IDXReader_Initialize(
  vtkClientServerInterpreter *csi)
{
	// printf("[A] : Inititalization\n");
	std::cout << "IDXReaderInit: IDXReaqder_Initialize" << endl;
  IDXReader_Init(csi);

}
