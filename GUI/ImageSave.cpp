#include "ImageSave.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>

#include "vtkSmartPointer.h"

#include "log.h"
#include <QFileDialog>

LOG_USE();


// Constructor
ImageSave::ImageSave(vtkSmartPointer<vtkRenderWindow> renWin)
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".", 
		tr("Image files (*.png *.jpg *.tif *.bmp)"));    
    if (fileName.compare("") == 0) {
        return;
    }
    QFileInfo fi(fileName);
    QString imgType = fi.suffix();
    vtkSmartPointer<vtkWindowToImageFilter> w2img = vtkWindowToImageFilter::New();
    w2img->SetInput(renWin);
    if (imgType.compare("png") == 0) {
        vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
        pngwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        pngwriter->SetFileName((fileName.toStdString()).c_str());
        pngwriter->Write();
//		pngwriter->Delete();	// Note: using vtkSmartPointer, delete is not necessary.
    } else if (imgType.compare("jpg") == 0) {
        vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
        jpgwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        jpgwriter->SetFileName((fileName.toStdString()).c_str());
        jpgwriter->Write();
//		jpgwriter->Delete();
    } else if (imgType.compare("tif") == 0) {
        vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
        tifwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        tifwriter->SetFileName((fileName.toStdString()).c_str());
        tifwriter->Write();
//		tifwriter->Delete();
    } else if (imgType.compare("bmp") == 0) {
        vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
        bmpwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        bmpwriter->SetFileName((fileName.toStdString()).c_str());
        bmpwriter->Write();
    }

};


