// myvtk.h
#ifndef MYVTK_H
#define MYVTK_H

#include <QtGui>
#include <QtCore>
#include <QIODevice>
#include <QVTKWidget.h>
#include <vtkRenderer.h> 
#include <vtkRenderWindow.h>
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include <vtkMPEG2Writer.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>

#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

//#include <vtkConfigure.h>

using namespace std;

struct cell_pos {
	int tag;
	int x, y, z;
	double diameter;
//	double state;
	int state;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int BCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

class MyVTK
{
public:
	MyVTK(QWidget *);
	~MyVTK();

	void createMappers();
	void read_cell_positions(QString, QString, bool);
	void get_cell_positions(bool fast);
	void init();
	void cleanup();
	void unpack(int x, double *, double *, double *);
	void renderCells(bool,bool);
	void process_Bcells();
    void process_Dcells(bool);
    void process_bonds();
	bool startPlayer(QString, QTimer *, bool);
	bool nextFrame();
	void pause();
	void playon();
	void saveSnapshot(QString, QString);
	void stop();

	QList<CELL_POS > BCpos_list;
	QList<CELL_POS > DCpos_list;
	QList<BOND_POS > bondpos_list;
	QList<vtkActor *> B_Actor_list;
	QList<vtkActor *> D_Actor_list;
	QList<vtkActor *> Bnd_Actor_list;

	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
	vtkPolyDataMapper *BcellMapper;
	vtkPolyDataMapper *DcellMapper;
	vtkPolyDataMapper *bondMapper;
	vtkPolyDataMapper *FDcellMapper;
	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
	vtkSmartPointer<vtkJPEGWriter> writer;
//	vtkSmartPointer<vtkTIFFWriter> writer;
	vtkSmartPointer<vtkImageCast> castFilter;
	vtkWindowToImageFilter *w2i;
	vtkWindowToImageFilter *w2img;
//	vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkJPEGWriter> jpgwriter;

	char msg[2048];
	double zoomlevel;
	double Pi;
	bool DCmotion;
	bool DCfade;
	bool first_VTK;
	bool playing;
	bool paused;
	bool save_image;
	QString casename;
	int framenum;
	QTimer *timer;
	QString infile;
	QFile *playerData;
	QTextStream *playerStream;

};

#endif
