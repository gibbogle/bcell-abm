#CONFIG      += uitools
CONFIG += release

#INCLUDEPATH +=  C:\VTK-build-4.4.0-static C:\VTK-build-4.4.0-static\rendering \
#                C:\VTK-src\GUISupport\Qt C:\VTK-src\common C:\VTK-src\rendering \
#                C:\VTK-src\graphics C:\VTK-src\filtering C:\VTK-src\IO \
#                C:\VTK-src\imaging \
INCLUDEPATH  += C:\users\mbog002\VTK-build4.4.0-static\include\vtk-5.8
INCLUDEPATH  += c:\users\mbog002\qwt-5.2.0\src

FORMS         = Bcell_GUI.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h myvtk.h misc.h \
                        libpara32.h result_set.h graphs.h
RESOURCES     += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp \
                myvtk.cpp misc.cpp lognormal.cpp graphs.cpp

# See cmake_link_command.txt for the full list of libraries that CMake links
LIBS += -L"C:\users\mbog002\VTK-build4.4.0-static\lib\vtk-5.8" -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32  \
-lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lpthread -lvtkalglib \
-lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32
# -lVPIC -lCosmo

#LIBS += -LC:\users\gib\abm\build32msvs\release -lpara32
LIBS += -Lc:\users\mbog002\qwt-5.2.0\lib -lqwt5

DEPENDPATH   += c:\users\mbog002\qwt-5.2.0\lib

QT           += network

QMAKE_CXXFLAGS += -Wno-deprecated -Wno-write-strings
#NOTE: fix for link with Qt-4.7
QMAKE_LFLAGS    += -Wl,-enable-auto-import

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS icons # ABM_GUI.pro
sources.path = .
INSTALLS += target sources

#DEFINES += __GFORTRAN_DLL__
DEFINES += __MSVS_DLL__
LIBS += C:\bin\bcell32_ms.dll

DEFINES += __COMPILETIME_LOADING__
# Note: On Windows compile-time DLL loading apparently does not permit OpenMP to use multiple threads

