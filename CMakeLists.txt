cmake_minimum_required(VERSION 2.8)
 
FIND_PACKAGE(ParaView REQUIRED)

SET(IDX_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
SET(IDX_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib)
SET(IDX_LIBRARIES visuskernel visusidx)

# tinyxml dependency
#ADD_SUBDIRECTORY(../nvisusio/libs/tinyxml        /usr/local/lib/tinyxml)


#ADD_DEPENDENCIES(visusviewer visusappkit)
#TARGET_LINK_LIBRARIES(visusviewer visusappkit)


INCLUDE(${PARAVIEW_USE_FILE})
ADD_PARAVIEW_PLUGIN(IDXReader "1.0"
  SERVER_MANAGER_XML IDXReader.xml
  SERVER_MANAGER_SOURCES IDXReader.cxx
  GUI_RESOURCE_FILES IDXReaderGUI.xml)

ADD_DEPENDENCIES(IDXReader visusidx visusdb visusdataflow visuskernel tinyxml curl crypto ssl FreeImage)
TARGET_LINK_LIBRARIES(IDXReader visusidx visusdb visusdataflow visuskernel tinyxml curl crypto ssl FreeImage)
