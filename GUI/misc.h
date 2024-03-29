#ifndef MISC_H
#define MISC_H

#include <QtGui>

#include <QThread>
#include <QTcpServer>

#include "libBcell32.h"
#include "SimpleView2DUI.h"

class SocketHandler : public QThread
 {
    Q_OBJECT

public:
    SocketHandler(int newport, QObject *parent = 0);
    ~SocketHandler();
    void run();
	void stop();
	void end();

	QTcpServer *tcpServer;
	int port;
	bool exiting;
	bool stopped;
	quint16 blockSize;
	QTcpSocket *socket;
	char msg[2048];
	static const int CPORT0 = 5000;
	static const int CPORT1 = 5001;

private slots:
	 void processor();
signals:
	 void sh_connected();
	 void sh_disconnected();
	 void sh_output(QString);
	 void sh_error(QString);
};

class ExecThread: public QThread
{
	Q_OBJECT
	QString inputFile;
public:
	ExecThread(QString);
	void run();
	void snapshot();
    void saveGradient2D(int k);
	void pause();
	void unpause();
	void stop();
	int ncpu;
	int nsteps;
	bool paused;
	bool stopped;
    QMutex mutex1, mutex2;
signals:
	void display();
	void summary();
};

bool quitMessage(QString);

#ifdef _WIN32
#include "windows.h"
#define sleep(n) Sleep(n)   // milliseconds
#else
#define sleep(n) usleep(1000*n)
#endif

#endif
