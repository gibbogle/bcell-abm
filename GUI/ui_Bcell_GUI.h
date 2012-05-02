/********************************************************************************
** Form generated from reading UI file 'Bcell_GUI.ui'
**
** Created: Tue 1. May 11:28:18 2012
**      by: Qt User Interface Compiler version 4.6.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BCELL_GUI_H
#define UI_BCELL_GUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMdiArea>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QRadioButton>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStackedWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTextBrowser>
#include <QtGui/QTextEdit>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qmylabel.h"
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action_saveAs;
    QAction *action_save;
    QAction *action_open_input;
    QAction *action_stop;
    QAction *action_run;
    QAction *action_inputs;
    QAction *action_outputs;
    QAction *action_VTK;
    QAction *action_pause;
    QAction *action_load_results;
    QAction *action_add_graph;
    QAction *action_remove_graph;
    QAction *action_remove_all;
    QAction *action_play_VTK;
    QAction *action_set_speed;
    QAction *action_save_snapshot;
    QWidget *centralwidget;
    QGridLayout *gridLayout_5;
    QStackedWidget *stackedWidget;
    QWidget *page_input;
    QVBoxLayout *verticalLayout;
    QTabWidget *tabs;
    QWidget *tab_B;
    QMyLabel *label_BC_AVIDITY_MEDIAN;
    QSlider *slider_TC_AVIDITY_MEDIAN;
    QLineEdit *line_BC_AVIDITY_MEDIAN;
    QMyLabel *label_BC_AVIDITY_SHAPE;
    QSlider *slider_TC_AVIDITY_SHAPE;
    QLineEdit *line_BC_AVIDITY_SHAPE;
    QWidget *layoutWidget;
    QGridLayout *gridLayout;
    QMyLabel *label_BC_COGNATE_FRACTION;
    QSpacerItem *horizontalSpacer_3;
    QSlider *slider_TC_COGNATE_FRACTION;
    QLineEdit *line_BC_COGNATE_FRACTION;
    QMyLabel *label_BC_STIM_RATE_CONSTANT;
    QSpacerItem *horizontalSpacer_4;
    QSlider *slider_TC_STIM_RATE_CONSTANT;
    QLineEdit *line_BC_STIM_RATE_CONSTANT;
    QMyLabel *label_BC_STIM_HALFLIFE;
    QSpacerItem *horizontalSpacer_5;
    QLineEdit *line_BC_STIM_HALFLIFE;
    QLabel *units_TC_STIM_HALFLIFE;
    QMyLabel *label_MOTILITY_BETA;
    QSpacerItem *horizontalSpacer_10;
    QLineEdit *line_MOTILITY_BETA;
    QMyLabel *label_MOTILITY_RHO;
    QSpacerItem *horizontalSpacer_11;
    QLineEdit *line_MOTILITY_RHO;
    QSlider *slider_TC_STIM_HALFLIFE;
    QSlider *slider_MOTILITY_BETA;
    QSlider *slider_MOTILITY_RHO;
    QwtPlot *qwtPlot_TC_AVIDITY;
    QwtPlot *qwtPlot_DIVIDE1;
    QMyLabel *label_DIVIDE1_MEDIAN;
    QSlider *slider_DIVIDE1_MEDIAN;
    QLineEdit *line_DIVIDE1_MEDIAN;
    QMyLabel *label_DIVIDE1_SHAPE;
    QSlider *slider_DIVIDE1_SHAPE;
    QLineEdit *line_DIVIDE1_SHAPE;
    QLabel *alabel_dist;
    QLineEdit *line_DIVIDE2_SHAPE;
    QSlider *slider_DIVIDE2_MEDIAN;
    QLineEdit *line_DIVIDE2_MEDIAN;
    QMyLabel *label_DIVIDE2_MEDIAN;
    QMyLabel *label_DIVIDE2_SHAPE;
    QSlider *slider_DIVIDE2_SHAPE;
    QwtPlot *qwtPlot_DIVIDE2;
    QWidget *tab_DC;
    QWidget *layoutWidget1;
    QGridLayout *gridLayout_DC;
    QMyLabel *label_DC_BIND_DELAY;
    QLabel *units_DC_BIND_DELAY;
    QMyLabel *label_DC_DENS_HALFLIFE;
    QLineEdit *line_DC_DENS_HALFLIFE;
    QMyLabel *label_MAX_TC_BIND;
    QMyLabel *label_MAX_COG_BIND;
    QLabel *units_DC_DENS_HALFLIFE;
    QSpinBox *spin_MAX_TC_BIND;
    QSpinBox *spin_MAX_COG_BIND;
    QLineEdit *line_DC_BIND_DELAY;
    QMyLabel *label_DC_ANTIGEN_SHAPE;
    QwtPlot *qwtPlot_DC_ANTIGEN;
    QLineEdit *line_DC_ANTIGEN_SHAPE;
    QMyLabel *label_DC_ANTIGEN_MEDIAN;
    QLineEdit *line_DC_ANTIGEN_MEDIAN;
    QSlider *slider_DC_ANTIGEN_SHAPE;
    QSlider *slider_DC_ANTIGEN_MEDIAN;
    QMyLabel *label_DC_LIFETIME_SHAPE;
    QwtPlot *qwtPlot_DC_LIFETIME;
    QLineEdit *line_DC_LIFETIME_SHAPE;
    QMyLabel *label_DC_LIFETIME_MEDIAN;
    QLineEdit *line_DC_LIFETIME_MEDIAN;
    QSlider *slider_DC_LIFETIME_SHAPE;
    QSlider *slider_DC_LIFETIME_MEDIAN;
    QLabel *alabel_dist_2;
    QWidget *tab_chemo;
    QCheckBox *cbox_USE_S1PR1;
    QWidget *gridLayoutWidget_2;
    QGridLayout *gridLayout_S1P_2;
    QMyLabel *label_CCL21_BDRY_CONC;
    QLineEdit *line_CCL21_BDRY_CONC;
    QMyLabel *label_CCL21_DIFF_COEFF;
    QLineEdit *line_CCL21_DIFF_COEFF;
    QMyLabel *label_CCL21_HALFLIFE;
    QLineEdit *line_CCL21_HALFLIFE;
    QLabel *units_CCL21conc;
    QLabel *units_CCL21diff;
    QLabel *units_CCL21life;
    QLabel *label_CCL21_STRENGTH;
    QLineEdit *line_CCL21_STRENGTH;
    QLabel *label_CCL21_BDRY_RATE;
    QLineEdit *line_CCL21_BDRY_RATE;
    QLabel *units_CCL21rate;
    QCheckBox *cbox_USE_CCR7;
    QCheckBox *cbox_USE_EBI2;
    QCheckBox *cbox_USE_CXCR5;
    QCheckBox *cbox_USE_S1PR2;
    QRadioButton *rbut_CCL21_BDRY_0;
    QWidget *gridLayoutWidget_5;
    QGridLayout *gridLayout_S1P_5;
    QMyLabel *label_CXCL13_BDRY_CONC;
    QLineEdit *line_CXCL13_BDRY_CONC;
    QMyLabel *label_CXCL13_DIFF_COEFF;
    QLineEdit *line_CXCL13_DIFF_COEFF;
    QMyLabel *label_CXCL13_HALFLIFE;
    QLineEdit *line_CXCL13_HALFLIFE;
    QLabel *units_CXCL13conc;
    QLabel *units_CXCL13diff;
    QLabel *units_CXCL13life;
    QLabel *label_CXCL13_STRENGTH;
    QLineEdit *line_CXCL13_STRENGTH;
    QLabel *label_CXCL13_BDRY_RATE;
    QLineEdit *line_CXCL13_BDRY_RATE;
    QLabel *units_CXCL13rate;
    QRadioButton *rbut_CXCL13_BDRY_0;
    QWidget *gridLayoutWidget_4;
    QGridLayout *gridLayout_S1P_4;
    QMyLabel *label_OXY_BDRY_CONC;
    QLineEdit *line_OXY_BDRY_CONC;
    QMyLabel *label_OXY_DIFF_COEFF;
    QLineEdit *line_OXY_DIFF_COEFF;
    QMyLabel *label_OXY_HALFLIFE;
    QLineEdit *line_OXY_HALFLIFE;
    QLabel *units_OXYconc;
    QLabel *units_OXYdiff;
    QLabel *units_OXYlife;
    QLabel *label_OXY_STRENGTH;
    QLineEdit *line_OXY_STRENGTH;
    QLabel *label_OXY_BDRY_RATE;
    QLineEdit *line_OXY_BDRY_RATE;
    QLabel *units_OXYrate;
    QRadioButton *rbut_OXY_BDRY_0;
    QWidget *gridLayoutWidget_3;
    QGridLayout *gridLayout_S1P_3;
    QMyLabel *label_S1P_BDRY_CONC;
    QLineEdit *line_S1P_BDRY_CONC;
    QMyLabel *label_S1P_DIFF_COEFF;
    QLineEdit *line_S1P_DIFF_COEFF;
    QMyLabel *label_S1P_HALFLIFE;
    QLineEdit *line_S1P_HALFLIFE;
    QLabel *units_S1Pconc;
    QLabel *units_S1Pdiff;
    QLabel *units_S1Plife;
    QLabel *label_S1P_STRENGTH_POS;
    QLineEdit *line_S1P_STRENGTH_POS;
    QLabel *label_S1P_STRENGTH_NEG;
    QLineEdit *line_S1P_STRENGTH_NEG;
    QLabel *label_S1P_BDRY_RATE;
    QLineEdit *line_S1P_BDRY_RATE;
    QLabel *units_S1Prate;
    QRadioButton *rbut_S1P_BDRY_0;
    QLabel *label_chemokine;
    QRadioButton *rbut_S1P_BDRY_1;
    QRadioButton *rbut_CCL21_BDRY_1;
    QRadioButton *rbut_OXY_BDRY_1;
    QRadioButton *rbut_CXCL13_BDRY_1;
    QWidget *tab_TCR;
    QWidget *layoutWidget_2;
    QGridLayout *gridLayout_3;
    QMyLabel *label_IL2_THRESHOLD;
    QMyLabel *label_ACTIVATION_THRESHOLD;
    QMyLabel *label_FIRST_DIVISION_THRESHOLD;
    QMyLabel *label_DIVISION_THRESHOLD;
    QMyLabel *label_EXIT_THRESHOLD;
    QLineEdit *line_IL2_THRESHOLD;
    QLineEdit *line_ACTIVATION_THRESHOLD;
    QLineEdit *line_FIRST_DIVISION_THRESHOLD;
    QLineEdit *line_DIVISION_THRESHOLD;
    QLineEdit *line_EXIT_THRESHOLD;
    QMyLabel *label_STIMULATION_LIMIT;
    QLineEdit *line_STIMULATION_LIMIT;
    QLabel *label_40;
    QWidget *tab_run;
    QWidget *layoutWidget2;
    QGridLayout *gridLayout_4;
    QMyLabel *label_NX;
    QSpacerItem *horizontalSpacer_12;
    QSpinBox *spin_NX;
    QMyLabel *label_BLOB_RADIUS;
    QSpacerItem *horizontalSpacer_13;
    QLineEdit *line_BLOB_RADIUS;
    QMyLabel *label_BC_FRACTION;
    QSpacerItem *horizontalSpacer_14;
    QLineEdit *line_BC_FRACTION;
    QSpacerItem *horizontalSpacer_15;
    QLineEdit *line_FLUID_FRACTION;
    QMyLabel *label_RESIDENCE_TIME;
    QSpacerItem *horizontalSpacer_21;
    QLineEdit *line_RESIDENCE_TIME;
    QLabel *units_RESIDENCE_TIME;
    QMyLabel *label_INFLAMM_DAYS1;
    QSpacerItem *horizontalSpacer_22;
    QLineEdit *line_INFLAMM_DAYS1;
    QMyLabel *label_INFLAMM_DAYS2;
    QSpacerItem *horizontalSpacer_23;
    QLineEdit *line_INFLAMM_DAYS2;
    QMyLabel *label_INFLAMM_LEVEL;
    QSpacerItem *horizontalSpacer_24;
    QLineEdit *line_INFLAMM_LEVEL;
    QMyLabel *label_CHEMO_RADIUS;
    QSpacerItem *horizontalSpacer_27;
    QLineEdit *line_CHEMO_RADIUS;
    QMyLabel *label_BASE_EXIT_PROB;
    QSpacerItem *horizontalSpacer_28;
    QLineEdit *line_BASE_EXIT_PROB;
    QMyLabel *label_NDAYS;
    QSpacerItem *horizontalSpacer_29;
    QMyLabel *label_SEED1;
    QSpacerItem *horizontalSpacer_30;
    QSpinBox *spin_SEED1;
    QMyLabel *label_SEED2;
    QSpacerItem *horizontalSpacer_31;
    QSpinBox *spin_SEED2;
    QLineEdit *line_NDAYS;
    QMyLabel *label_FLUID_FRACTION;
    QMyLabel *label_NT_ANIMATION;
    QSpacerItem *horizontalSpacer_32;
    QSpinBox *spin_NT_ANIMATION;
    QLabel *units_CHEMO_RADIUS;
    QLabel *units_NDAYS;
    QLabel *units_INFLAMM1;
    QLabel *units_INFLAMM2;
    QLabel *units_BLOB_RADIUS;
    QMyLabel *label_NCPU;
    QSpacerItem *horizontalSpacer_34;
    QSpinBox *spin_NCPU;
    QCheckBox *cbox_savepos;
    QCheckBox *cbox_IV_SHOW_NONCOGNATE;
    QCheckBox *cbox_USE_TRAFFIC;
    QRadioButton *rbut_SPECIES_1;
    QRadioButton *rbut_SPECIES_0;
    QCheckBox *cbox_USE_EXIT_CHEMOTAXIS;
    QCheckBox *cbox_COMPUTED_OUTFLOW;
    QMyLabel *label_INPUT_FILE;
    QLineEdit *text_INPUT_FILE;
    QLabel *label_input;
    QTextEdit *text_more;
    QWidget *page_output;
    QVBoxLayout *verticalLayout_2;
    QMdiArea *mdiArea;
    QTextBrowser *box_outputLog;
    QWidget *page_3D;
    QMdiArea *mdiArea_VTK;
    QProgressBar *progressBar;
    QLabel *label_hour;
    QLabel *hour_display;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuABM;
    QMenu *menuGraphs;
    QMenu *menuPlayer;
    QMenu *menuSnapshot;
    QStatusBar *statusbar;
    QToolBar *toolBar1;
    QButtonGroup *buttonGroup_CXCL13_BDRY;
    QButtonGroup *buttonGroup_SPECIES;
    QButtonGroup *buttonGroup_S1P_BDRY;
    QButtonGroup *buttonGroup_OXY_BDRY;
    QButtonGroup *buttonGroup_CCL21_BDRY;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1319, 977);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/ABM.png"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        MainWindow->setStyleSheet(QString::fromUtf8("#tabs{\n"
"color: white;\n"
"background-color: QLinearGradient( x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #88d, stop: 0.1 #99e, stop: 0.49 #77c, stop: 0.5 #66b, stop: 1 #77c);\n"
"border-width: 1px;\n"
"border-color: #339;\n"
"border-style: solid;\n"
"border-radius: 15;\n"
"}\n"
"\n"
"#text_more{\n"
"padding: 8px\n"
"}"));
        action_saveAs = new QAction(MainWindow);
        action_saveAs->setObjectName(QString::fromUtf8("action_saveAs"));
        action_save = new QAction(MainWindow);
        action_save->setObjectName(QString::fromUtf8("action_save"));
        action_open_input = new QAction(MainWindow);
        action_open_input->setObjectName(QString::fromUtf8("action_open_input"));
        action_stop = new QAction(MainWindow);
        action_stop->setObjectName(QString::fromUtf8("action_stop"));
        action_stop->setCheckable(false);
        action_stop->setEnabled(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/001_29.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_stop->setIcon(icon1);
        action_run = new QAction(MainWindow);
        action_run->setObjectName(QString::fromUtf8("action_run"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/icons/001_59.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_run->setIcon(icon2);
        action_inputs = new QAction(MainWindow);
        action_inputs->setObjectName(QString::fromUtf8("action_inputs"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/icons/001_45.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_inputs->setIcon(icon3);
        action_outputs = new QAction(MainWindow);
        action_outputs->setObjectName(QString::fromUtf8("action_outputs"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/icons/Display2.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_outputs->setIcon(icon4);
        action_VTK = new QAction(MainWindow);
        action_VTK->setObjectName(QString::fromUtf8("action_VTK"));
        action_VTK->setEnabled(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/icons/cell2.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_VTK->setIcon(icon5);
        action_pause = new QAction(MainWindow);
        action_pause->setObjectName(QString::fromUtf8("action_pause"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/icons/001_07.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_pause->setIcon(icon6);
        action_load_results = new QAction(MainWindow);
        action_load_results->setObjectName(QString::fromUtf8("action_load_results"));
        action_add_graph = new QAction(MainWindow);
        action_add_graph->setObjectName(QString::fromUtf8("action_add_graph"));
        action_remove_graph = new QAction(MainWindow);
        action_remove_graph->setObjectName(QString::fromUtf8("action_remove_graph"));
        action_remove_all = new QAction(MainWindow);
        action_remove_all->setObjectName(QString::fromUtf8("action_remove_all"));
        action_play_VTK = new QAction(MainWindow);
        action_play_VTK->setObjectName(QString::fromUtf8("action_play_VTK"));
        action_set_speed = new QAction(MainWindow);
        action_set_speed->setObjectName(QString::fromUtf8("action_set_speed"));
        action_save_snapshot = new QAction(MainWindow);
        action_save_snapshot->setObjectName(QString::fromUtf8("action_save_snapshot"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_5 = new QGridLayout(centralwidget);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        stackedWidget = new QStackedWidget(centralwidget);
        stackedWidget->setObjectName(QString::fromUtf8("stackedWidget"));
        page_input = new QWidget();
        page_input->setObjectName(QString::fromUtf8("page_input"));
        verticalLayout = new QVBoxLayout(page_input);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(-1, 0, -1, -1);
        tabs = new QTabWidget(page_input);
        tabs->setObjectName(QString::fromUtf8("tabs"));
        tabs->setEnabled(true);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(tabs->sizePolicy().hasHeightForWidth());
        tabs->setSizePolicy(sizePolicy);
        tabs->setMinimumSize(QSize(0, 0));
        tabs->setAutoFillBackground(false);
        tabs->setTabPosition(QTabWidget::West);
        tabs->setElideMode(Qt::ElideNone);
        tab_B = new QWidget();
        tab_B->setObjectName(QString::fromUtf8("tab_B"));
        tab_B->setEnabled(true);
        label_BC_AVIDITY_MEDIAN = new QMyLabel(tab_B);
        label_BC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("label_BC_AVIDITY_MEDIAN"));
        label_BC_AVIDITY_MEDIAN->setGeometry(QRect(760, 80, 61, 20));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_BC_AVIDITY_MEDIAN->sizePolicy().hasHeightForWidth());
        label_BC_AVIDITY_MEDIAN->setSizePolicy(sizePolicy1);
        label_BC_AVIDITY_MEDIAN->setCursor(QCursor(Qt::ArrowCursor));
        label_BC_AVIDITY_MEDIAN->setMouseTracking(true);
        label_BC_AVIDITY_MEDIAN->setWordWrap(false);
        label_BC_AVIDITY_MEDIAN->setTextInteractionFlags(Qt::LinksAccessibleByMouse);
        slider_TC_AVIDITY_MEDIAN = new QSlider(tab_B);
        slider_TC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("slider_TC_AVIDITY_MEDIAN"));
        slider_TC_AVIDITY_MEDIAN->setGeometry(QRect(823, 80, 61, 20));
        slider_TC_AVIDITY_MEDIAN->setOrientation(Qt::Horizontal);
        line_BC_AVIDITY_MEDIAN = new QLineEdit(tab_B);
        line_BC_AVIDITY_MEDIAN->setObjectName(QString::fromUtf8("line_BC_AVIDITY_MEDIAN"));
        line_BC_AVIDITY_MEDIAN->setEnabled(false);
        line_BC_AVIDITY_MEDIAN->setGeometry(QRect(910, 80, 51, 20));
        line_BC_AVIDITY_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_BC_AVIDITY_SHAPE = new QMyLabel(tab_B);
        label_BC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("label_BC_AVIDITY_SHAPE"));
        label_BC_AVIDITY_SHAPE->setGeometry(QRect(760, 120, 51, 20));
        sizePolicy1.setHeightForWidth(label_BC_AVIDITY_SHAPE->sizePolicy().hasHeightForWidth());
        label_BC_AVIDITY_SHAPE->setSizePolicy(sizePolicy1);
        label_BC_AVIDITY_SHAPE->setWordWrap(false);
        slider_TC_AVIDITY_SHAPE = new QSlider(tab_B);
        slider_TC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("slider_TC_AVIDITY_SHAPE"));
        slider_TC_AVIDITY_SHAPE->setGeometry(QRect(823, 120, 61, 20));
        slider_TC_AVIDITY_SHAPE->setOrientation(Qt::Horizontal);
        line_BC_AVIDITY_SHAPE = new QLineEdit(tab_B);
        line_BC_AVIDITY_SHAPE->setObjectName(QString::fromUtf8("line_BC_AVIDITY_SHAPE"));
        line_BC_AVIDITY_SHAPE->setEnabled(false);
        line_BC_AVIDITY_SHAPE->setGeometry(QRect(910, 120, 51, 20));
        line_BC_AVIDITY_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        layoutWidget = new QWidget(tab_B);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(0, 0, 391, 251));
        gridLayout = new QGridLayout(layoutWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_BC_COGNATE_FRACTION = new QMyLabel(layoutWidget);
        label_BC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("label_BC_COGNATE_FRACTION"));
        sizePolicy1.setHeightForWidth(label_BC_COGNATE_FRACTION->sizePolicy().hasHeightForWidth());
        label_BC_COGNATE_FRACTION->setSizePolicy(sizePolicy1);
        label_BC_COGNATE_FRACTION->setWordWrap(false);

        gridLayout->addWidget(label_BC_COGNATE_FRACTION, 0, 0, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 0, 1, 1, 1);

        slider_TC_COGNATE_FRACTION = new QSlider(layoutWidget);
        slider_TC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("slider_TC_COGNATE_FRACTION"));
        slider_TC_COGNATE_FRACTION->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_COGNATE_FRACTION, 0, 2, 1, 1);

        line_BC_COGNATE_FRACTION = new QLineEdit(layoutWidget);
        line_BC_COGNATE_FRACTION->setObjectName(QString::fromUtf8("line_BC_COGNATE_FRACTION"));
        line_BC_COGNATE_FRACTION->setMaximumSize(QSize(90, 16777215));
        line_BC_COGNATE_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_BC_COGNATE_FRACTION, 0, 3, 1, 1);

        label_BC_STIM_RATE_CONSTANT = new QMyLabel(layoutWidget);
        label_BC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("label_BC_STIM_RATE_CONSTANT"));
        sizePolicy1.setHeightForWidth(label_BC_STIM_RATE_CONSTANT->sizePolicy().hasHeightForWidth());
        label_BC_STIM_RATE_CONSTANT->setSizePolicy(sizePolicy1);
        label_BC_STIM_RATE_CONSTANT->setMouseTracking(false);
        label_BC_STIM_RATE_CONSTANT->setWordWrap(false);

        gridLayout->addWidget(label_BC_STIM_RATE_CONSTANT, 1, 0, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_4, 1, 1, 1, 1);

        slider_TC_STIM_RATE_CONSTANT = new QSlider(layoutWidget);
        slider_TC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("slider_TC_STIM_RATE_CONSTANT"));
        slider_TC_STIM_RATE_CONSTANT->setEnabled(true);
        slider_TC_STIM_RATE_CONSTANT->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_STIM_RATE_CONSTANT, 1, 2, 1, 1);

        line_BC_STIM_RATE_CONSTANT = new QLineEdit(layoutWidget);
        line_BC_STIM_RATE_CONSTANT->setObjectName(QString::fromUtf8("line_BC_STIM_RATE_CONSTANT"));
        line_BC_STIM_RATE_CONSTANT->setEnabled(false);
        line_BC_STIM_RATE_CONSTANT->setMaximumSize(QSize(90, 16777215));
        line_BC_STIM_RATE_CONSTANT->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_BC_STIM_RATE_CONSTANT, 1, 3, 1, 1);

        label_BC_STIM_HALFLIFE = new QMyLabel(layoutWidget);
        label_BC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("label_BC_STIM_HALFLIFE"));
        sizePolicy1.setHeightForWidth(label_BC_STIM_HALFLIFE->sizePolicy().hasHeightForWidth());
        label_BC_STIM_HALFLIFE->setSizePolicy(sizePolicy1);
        label_BC_STIM_HALFLIFE->setMouseTracking(false);
        label_BC_STIM_HALFLIFE->setWordWrap(false);

        gridLayout->addWidget(label_BC_STIM_HALFLIFE, 2, 0, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_5, 2, 1, 1, 1);

        line_BC_STIM_HALFLIFE = new QLineEdit(layoutWidget);
        line_BC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("line_BC_STIM_HALFLIFE"));
        line_BC_STIM_HALFLIFE->setEnabled(false);
        line_BC_STIM_HALFLIFE->setMaximumSize(QSize(90, 16777215));
        line_BC_STIM_HALFLIFE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_BC_STIM_HALFLIFE, 2, 3, 1, 1);

        units_TC_STIM_HALFLIFE = new QLabel(layoutWidget);
        units_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("units_TC_STIM_HALFLIFE"));

        gridLayout->addWidget(units_TC_STIM_HALFLIFE, 2, 4, 1, 1);

        label_MOTILITY_BETA = new QMyLabel(layoutWidget);
        label_MOTILITY_BETA->setObjectName(QString::fromUtf8("label_MOTILITY_BETA"));
        sizePolicy1.setHeightForWidth(label_MOTILITY_BETA->sizePolicy().hasHeightForWidth());
        label_MOTILITY_BETA->setSizePolicy(sizePolicy1);
        label_MOTILITY_BETA->setMouseTracking(false);
        label_MOTILITY_BETA->setWordWrap(false);

        gridLayout->addWidget(label_MOTILITY_BETA, 3, 0, 1, 1);

        horizontalSpacer_10 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_10, 3, 1, 1, 1);

        line_MOTILITY_BETA = new QLineEdit(layoutWidget);
        line_MOTILITY_BETA->setObjectName(QString::fromUtf8("line_MOTILITY_BETA"));
        line_MOTILITY_BETA->setMaximumSize(QSize(90, 16777215));
        line_MOTILITY_BETA->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_MOTILITY_BETA, 3, 3, 1, 1);

        label_MOTILITY_RHO = new QMyLabel(layoutWidget);
        label_MOTILITY_RHO->setObjectName(QString::fromUtf8("label_MOTILITY_RHO"));
        sizePolicy1.setHeightForWidth(label_MOTILITY_RHO->sizePolicy().hasHeightForWidth());
        label_MOTILITY_RHO->setSizePolicy(sizePolicy1);
        label_MOTILITY_RHO->setMouseTracking(false);
        label_MOTILITY_RHO->setWordWrap(false);

        gridLayout->addWidget(label_MOTILITY_RHO, 4, 0, 1, 1);

        horizontalSpacer_11 = new QSpacerItem(13, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_11, 4, 1, 1, 1);

        line_MOTILITY_RHO = new QLineEdit(layoutWidget);
        line_MOTILITY_RHO->setObjectName(QString::fromUtf8("line_MOTILITY_RHO"));
        line_MOTILITY_RHO->setMaximumSize(QSize(90, 16777215));
        line_MOTILITY_RHO->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(line_MOTILITY_RHO, 4, 3, 1, 1);

        slider_TC_STIM_HALFLIFE = new QSlider(layoutWidget);
        slider_TC_STIM_HALFLIFE->setObjectName(QString::fromUtf8("slider_TC_STIM_HALFLIFE"));
        slider_TC_STIM_HALFLIFE->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_TC_STIM_HALFLIFE, 2, 2, 1, 1);

        slider_MOTILITY_BETA = new QSlider(layoutWidget);
        slider_MOTILITY_BETA->setObjectName(QString::fromUtf8("slider_MOTILITY_BETA"));
        slider_MOTILITY_BETA->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_MOTILITY_BETA, 3, 2, 1, 1);

        slider_MOTILITY_RHO = new QSlider(layoutWidget);
        slider_MOTILITY_RHO->setObjectName(QString::fromUtf8("slider_MOTILITY_RHO"));
        slider_MOTILITY_RHO->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(slider_MOTILITY_RHO, 4, 2, 1, 1);

        qwtPlot_TC_AVIDITY = new QwtPlot(tab_B);
        qwtPlot_TC_AVIDITY->setObjectName(QString::fromUtf8("qwtPlot_TC_AVIDITY"));
        qwtPlot_TC_AVIDITY->setGeometry(QRect(420, 50, 300, 170));
        qwtPlot_DIVIDE1 = new QwtPlot(tab_B);
        qwtPlot_DIVIDE1->setObjectName(QString::fromUtf8("qwtPlot_DIVIDE1"));
        qwtPlot_DIVIDE1->setGeometry(QRect(420, 250, 300, 170));
        label_DIVIDE1_MEDIAN = new QMyLabel(tab_B);
        label_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("label_DIVIDE1_MEDIAN"));
        label_DIVIDE1_MEDIAN->setGeometry(QRect(760, 280, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE1_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DIVIDE1_MEDIAN->setSizePolicy(sizePolicy1);
        label_DIVIDE1_MEDIAN->setMouseTracking(false);
        label_DIVIDE1_MEDIAN->setWordWrap(false);
        slider_DIVIDE1_MEDIAN = new QSlider(tab_B);
        slider_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("slider_DIVIDE1_MEDIAN"));
        slider_DIVIDE1_MEDIAN->setGeometry(QRect(820, 280, 61, 16));
        slider_DIVIDE1_MEDIAN->setOrientation(Qt::Horizontal);
        line_DIVIDE1_MEDIAN = new QLineEdit(tab_B);
        line_DIVIDE1_MEDIAN->setObjectName(QString::fromUtf8("line_DIVIDE1_MEDIAN"));
        line_DIVIDE1_MEDIAN->setGeometry(QRect(910, 280, 51, 20));
        line_DIVIDE1_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DIVIDE1_SHAPE = new QMyLabel(tab_B);
        label_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("label_DIVIDE1_SHAPE"));
        label_DIVIDE1_SHAPE->setGeometry(QRect(760, 320, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE1_SHAPE->sizePolicy().hasHeightForWidth());
        label_DIVIDE1_SHAPE->setSizePolicy(sizePolicy1);
        label_DIVIDE1_SHAPE->setMouseTracking(false);
        label_DIVIDE1_SHAPE->setWordWrap(false);
        slider_DIVIDE1_SHAPE = new QSlider(tab_B);
        slider_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("slider_DIVIDE1_SHAPE"));
        slider_DIVIDE1_SHAPE->setGeometry(QRect(820, 320, 61, 16));
        slider_DIVIDE1_SHAPE->setOrientation(Qt::Horizontal);
        line_DIVIDE1_SHAPE = new QLineEdit(tab_B);
        line_DIVIDE1_SHAPE->setObjectName(QString::fromUtf8("line_DIVIDE1_SHAPE"));
        line_DIVIDE1_SHAPE->setGeometry(QRect(910, 320, 51, 20));
        line_DIVIDE1_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        alabel_dist = new QLabel(tab_B);
        alabel_dist->setObjectName(QString::fromUtf8("alabel_dist"));
        alabel_dist->setGeometry(QRect(550, 10, 291, 20));
        QFont font;
        font.setPointSize(14);
        font.setBold(true);
        font.setWeight(75);
        alabel_dist->setFont(font);
        line_DIVIDE2_SHAPE = new QLineEdit(tab_B);
        line_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("line_DIVIDE2_SHAPE"));
        line_DIVIDE2_SHAPE->setGeometry(QRect(910, 540, 51, 20));
        line_DIVIDE2_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DIVIDE2_MEDIAN = new QSlider(tab_B);
        slider_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("slider_DIVIDE2_MEDIAN"));
        slider_DIVIDE2_MEDIAN->setGeometry(QRect(820, 500, 61, 16));
        slider_DIVIDE2_MEDIAN->setOrientation(Qt::Horizontal);
        line_DIVIDE2_MEDIAN = new QLineEdit(tab_B);
        line_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("line_DIVIDE2_MEDIAN"));
        line_DIVIDE2_MEDIAN->setGeometry(QRect(910, 500, 51, 20));
        line_DIVIDE2_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DIVIDE2_MEDIAN = new QMyLabel(tab_B);
        label_DIVIDE2_MEDIAN->setObjectName(QString::fromUtf8("label_DIVIDE2_MEDIAN"));
        label_DIVIDE2_MEDIAN->setGeometry(QRect(760, 500, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE2_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DIVIDE2_MEDIAN->setSizePolicy(sizePolicy1);
        label_DIVIDE2_MEDIAN->setMouseTracking(false);
        label_DIVIDE2_MEDIAN->setWordWrap(false);
        label_DIVIDE2_SHAPE = new QMyLabel(tab_B);
        label_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("label_DIVIDE2_SHAPE"));
        label_DIVIDE2_SHAPE->setGeometry(QRect(760, 540, 61, 20));
        sizePolicy1.setHeightForWidth(label_DIVIDE2_SHAPE->sizePolicy().hasHeightForWidth());
        label_DIVIDE2_SHAPE->setSizePolicy(sizePolicy1);
        label_DIVIDE2_SHAPE->setMouseTracking(false);
        label_DIVIDE2_SHAPE->setWordWrap(false);
        slider_DIVIDE2_SHAPE = new QSlider(tab_B);
        slider_DIVIDE2_SHAPE->setObjectName(QString::fromUtf8("slider_DIVIDE2_SHAPE"));
        slider_DIVIDE2_SHAPE->setGeometry(QRect(820, 540, 61, 16));
        slider_DIVIDE2_SHAPE->setOrientation(Qt::Horizontal);
        qwtPlot_DIVIDE2 = new QwtPlot(tab_B);
        qwtPlot_DIVIDE2->setObjectName(QString::fromUtf8("qwtPlot_DIVIDE2"));
        qwtPlot_DIVIDE2->setGeometry(QRect(420, 450, 300, 170));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Arial"));
        font1.setPointSize(4);
        qwtPlot_DIVIDE2->setFont(font1);
        tabs->addTab(tab_B, QString());
        tab_DC = new QWidget();
        tab_DC->setObjectName(QString::fromUtf8("tab_DC"));
        tab_DC->setEnabled(true);
        layoutWidget1 = new QWidget(tab_DC);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 20, 361, 181));
        gridLayout_DC = new QGridLayout(layoutWidget1);
        gridLayout_DC->setObjectName(QString::fromUtf8("gridLayout_DC"));
        gridLayout_DC->setContentsMargins(0, 0, 0, 0);
        label_DC_BIND_DELAY = new QMyLabel(layoutWidget1);
        label_DC_BIND_DELAY->setObjectName(QString::fromUtf8("label_DC_BIND_DELAY"));
        label_DC_BIND_DELAY->setEnabled(true);
        label_DC_BIND_DELAY->setMouseTracking(false);
        label_DC_BIND_DELAY->setWordWrap(false);

        gridLayout_DC->addWidget(label_DC_BIND_DELAY, 0, 0, 1, 2);

        units_DC_BIND_DELAY = new QLabel(layoutWidget1);
        units_DC_BIND_DELAY->setObjectName(QString::fromUtf8("units_DC_BIND_DELAY"));

        gridLayout_DC->addWidget(units_DC_BIND_DELAY, 0, 3, 1, 1);

        label_DC_DENS_HALFLIFE = new QMyLabel(layoutWidget1);
        label_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("label_DC_DENS_HALFLIFE"));
        label_DC_DENS_HALFLIFE->setMouseTracking(false);
        label_DC_DENS_HALFLIFE->setWordWrap(false);

        gridLayout_DC->addWidget(label_DC_DENS_HALFLIFE, 1, 0, 1, 2);

        line_DC_DENS_HALFLIFE = new QLineEdit(layoutWidget1);
        line_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("line_DC_DENS_HALFLIFE"));
        line_DC_DENS_HALFLIFE->setMaximumSize(QSize(90, 16777215));
        line_DC_DENS_HALFLIFE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_DC->addWidget(line_DC_DENS_HALFLIFE, 1, 2, 1, 1);

        label_MAX_TC_BIND = new QMyLabel(layoutWidget1);
        label_MAX_TC_BIND->setObjectName(QString::fromUtf8("label_MAX_TC_BIND"));
        label_MAX_TC_BIND->setMouseTracking(false);
        label_MAX_TC_BIND->setWordWrap(false);

        gridLayout_DC->addWidget(label_MAX_TC_BIND, 2, 0, 1, 2);

        label_MAX_COG_BIND = new QMyLabel(layoutWidget1);
        label_MAX_COG_BIND->setObjectName(QString::fromUtf8("label_MAX_COG_BIND"));
        label_MAX_COG_BIND->setMouseTracking(false);
        label_MAX_COG_BIND->setWordWrap(false);

        gridLayout_DC->addWidget(label_MAX_COG_BIND, 3, 0, 1, 2);

        units_DC_DENS_HALFLIFE = new QLabel(layoutWidget1);
        units_DC_DENS_HALFLIFE->setObjectName(QString::fromUtf8("units_DC_DENS_HALFLIFE"));

        gridLayout_DC->addWidget(units_DC_DENS_HALFLIFE, 1, 3, 1, 1);

        spin_MAX_TC_BIND = new QSpinBox(layoutWidget1);
        spin_MAX_TC_BIND->setObjectName(QString::fromUtf8("spin_MAX_TC_BIND"));
        spin_MAX_TC_BIND->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_DC->addWidget(spin_MAX_TC_BIND, 2, 2, 1, 1);

        spin_MAX_COG_BIND = new QSpinBox(layoutWidget1);
        spin_MAX_COG_BIND->setObjectName(QString::fromUtf8("spin_MAX_COG_BIND"));
        spin_MAX_COG_BIND->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_DC->addWidget(spin_MAX_COG_BIND, 3, 2, 1, 1);

        line_DC_BIND_DELAY = new QLineEdit(layoutWidget1);
        line_DC_BIND_DELAY->setObjectName(QString::fromUtf8("line_DC_BIND_DELAY"));
        line_DC_BIND_DELAY->setMaximumSize(QSize(90, 16777215));
        line_DC_BIND_DELAY->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_DC->addWidget(line_DC_BIND_DELAY, 0, 2, 1, 1);

        label_DC_ANTIGEN_SHAPE = new QMyLabel(tab_DC);
        label_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("label_DC_ANTIGEN_SHAPE"));
        label_DC_ANTIGEN_SHAPE->setGeometry(QRect(790, 150, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_ANTIGEN_SHAPE->sizePolicy().hasHeightForWidth());
        label_DC_ANTIGEN_SHAPE->setSizePolicy(sizePolicy1);
        label_DC_ANTIGEN_SHAPE->setMouseTracking(false);
        label_DC_ANTIGEN_SHAPE->setWordWrap(false);
        qwtPlot_DC_ANTIGEN = new QwtPlot(tab_DC);
        qwtPlot_DC_ANTIGEN->setObjectName(QString::fromUtf8("qwtPlot_DC_ANTIGEN"));
        qwtPlot_DC_ANTIGEN->setGeometry(QRect(430, 70, 300, 200));
        line_DC_ANTIGEN_SHAPE = new QLineEdit(tab_DC);
        line_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("line_DC_ANTIGEN_SHAPE"));
        line_DC_ANTIGEN_SHAPE->setGeometry(QRect(940, 150, 51, 20));
        line_DC_ANTIGEN_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DC_ANTIGEN_MEDIAN = new QMyLabel(tab_DC);
        label_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("label_DC_ANTIGEN_MEDIAN"));
        label_DC_ANTIGEN_MEDIAN->setGeometry(QRect(790, 110, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_ANTIGEN_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DC_ANTIGEN_MEDIAN->setSizePolicy(sizePolicy1);
        label_DC_ANTIGEN_MEDIAN->setMouseTracking(false);
        label_DC_ANTIGEN_MEDIAN->setWordWrap(false);
        line_DC_ANTIGEN_MEDIAN = new QLineEdit(tab_DC);
        line_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("line_DC_ANTIGEN_MEDIAN"));
        line_DC_ANTIGEN_MEDIAN->setGeometry(QRect(940, 110, 51, 20));
        line_DC_ANTIGEN_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DC_ANTIGEN_SHAPE = new QSlider(tab_DC);
        slider_DC_ANTIGEN_SHAPE->setObjectName(QString::fromUtf8("slider_DC_ANTIGEN_SHAPE"));
        slider_DC_ANTIGEN_SHAPE->setGeometry(QRect(850, 150, 61, 16));
        slider_DC_ANTIGEN_SHAPE->setOrientation(Qt::Horizontal);
        slider_DC_ANTIGEN_MEDIAN = new QSlider(tab_DC);
        slider_DC_ANTIGEN_MEDIAN->setObjectName(QString::fromUtf8("slider_DC_ANTIGEN_MEDIAN"));
        slider_DC_ANTIGEN_MEDIAN->setGeometry(QRect(850, 110, 61, 16));
        slider_DC_ANTIGEN_MEDIAN->setOrientation(Qt::Horizontal);
        label_DC_LIFETIME_SHAPE = new QMyLabel(tab_DC);
        label_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("label_DC_LIFETIME_SHAPE"));
        label_DC_LIFETIME_SHAPE->setGeometry(QRect(790, 430, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_LIFETIME_SHAPE->sizePolicy().hasHeightForWidth());
        label_DC_LIFETIME_SHAPE->setSizePolicy(sizePolicy1);
        label_DC_LIFETIME_SHAPE->setMouseTracking(false);
        label_DC_LIFETIME_SHAPE->setWordWrap(false);
        qwtPlot_DC_LIFETIME = new QwtPlot(tab_DC);
        qwtPlot_DC_LIFETIME->setObjectName(QString::fromUtf8("qwtPlot_DC_LIFETIME"));
        qwtPlot_DC_LIFETIME->setGeometry(QRect(430, 350, 300, 200));
        line_DC_LIFETIME_SHAPE = new QLineEdit(tab_DC);
        line_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("line_DC_LIFETIME_SHAPE"));
        line_DC_LIFETIME_SHAPE->setGeometry(QRect(940, 430, 51, 20));
        line_DC_LIFETIME_SHAPE->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_DC_LIFETIME_MEDIAN = new QMyLabel(tab_DC);
        label_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("label_DC_LIFETIME_MEDIAN"));
        label_DC_LIFETIME_MEDIAN->setGeometry(QRect(790, 390, 61, 20));
        sizePolicy1.setHeightForWidth(label_DC_LIFETIME_MEDIAN->sizePolicy().hasHeightForWidth());
        label_DC_LIFETIME_MEDIAN->setSizePolicy(sizePolicy1);
        label_DC_LIFETIME_MEDIAN->setMouseTracking(false);
        label_DC_LIFETIME_MEDIAN->setWordWrap(false);
        line_DC_LIFETIME_MEDIAN = new QLineEdit(tab_DC);
        line_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("line_DC_LIFETIME_MEDIAN"));
        line_DC_LIFETIME_MEDIAN->setGeometry(QRect(940, 390, 51, 20));
        line_DC_LIFETIME_MEDIAN->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        slider_DC_LIFETIME_SHAPE = new QSlider(tab_DC);
        slider_DC_LIFETIME_SHAPE->setObjectName(QString::fromUtf8("slider_DC_LIFETIME_SHAPE"));
        slider_DC_LIFETIME_SHAPE->setGeometry(QRect(850, 430, 61, 16));
        slider_DC_LIFETIME_SHAPE->setOrientation(Qt::Horizontal);
        slider_DC_LIFETIME_MEDIAN = new QSlider(tab_DC);
        slider_DC_LIFETIME_MEDIAN->setObjectName(QString::fromUtf8("slider_DC_LIFETIME_MEDIAN"));
        slider_DC_LIFETIME_MEDIAN->setGeometry(QRect(850, 390, 61, 16));
        slider_DC_LIFETIME_MEDIAN->setOrientation(Qt::Horizontal);
        alabel_dist_2 = new QLabel(tab_DC);
        alabel_dist_2->setObjectName(QString::fromUtf8("alabel_dist_2"));
        alabel_dist_2->setGeometry(QRect(580, 20, 291, 20));
        alabel_dist_2->setFont(font);
        tabs->addTab(tab_DC, QString());
        tab_chemo = new QWidget();
        tab_chemo->setObjectName(QString::fromUtf8("tab_chemo"));
        cbox_USE_S1PR1 = new QCheckBox(tab_chemo);
        cbox_USE_S1PR1->setObjectName(QString::fromUtf8("cbox_USE_S1PR1"));
        cbox_USE_S1PR1->setGeometry(QRect(80, 100, 81, 18));
        gridLayoutWidget_2 = new QWidget(tab_chemo);
        gridLayoutWidget_2->setObjectName(QString::fromUtf8("gridLayoutWidget_2"));
        gridLayoutWidget_2->setGeometry(QRect(80, 370, 211, 126));
        gridLayout_S1P_2 = new QGridLayout(gridLayoutWidget_2);
        gridLayout_S1P_2->setObjectName(QString::fromUtf8("gridLayout_S1P_2"));
        gridLayout_S1P_2->setContentsMargins(0, 0, 0, 0);
        label_CCL21_BDRY_CONC = new QMyLabel(gridLayoutWidget_2);
        label_CCL21_BDRY_CONC->setObjectName(QString::fromUtf8("label_CCL21_BDRY_CONC"));

        gridLayout_S1P_2->addWidget(label_CCL21_BDRY_CONC, 1, 0, 1, 2);

        line_CCL21_BDRY_CONC = new QLineEdit(gridLayoutWidget_2);
        line_CCL21_BDRY_CONC->setObjectName(QString::fromUtf8("line_CCL21_BDRY_CONC"));

        gridLayout_S1P_2->addWidget(line_CCL21_BDRY_CONC, 1, 2, 1, 1);

        label_CCL21_DIFF_COEFF = new QMyLabel(gridLayoutWidget_2);
        label_CCL21_DIFF_COEFF->setObjectName(QString::fromUtf8("label_CCL21_DIFF_COEFF"));

        gridLayout_S1P_2->addWidget(label_CCL21_DIFF_COEFF, 2, 0, 1, 2);

        line_CCL21_DIFF_COEFF = new QLineEdit(gridLayoutWidget_2);
        line_CCL21_DIFF_COEFF->setObjectName(QString::fromUtf8("line_CCL21_DIFF_COEFF"));

        gridLayout_S1P_2->addWidget(line_CCL21_DIFF_COEFF, 2, 2, 1, 1);

        label_CCL21_HALFLIFE = new QMyLabel(gridLayoutWidget_2);
        label_CCL21_HALFLIFE->setObjectName(QString::fromUtf8("label_CCL21_HALFLIFE"));

        gridLayout_S1P_2->addWidget(label_CCL21_HALFLIFE, 3, 0, 1, 2);

        line_CCL21_HALFLIFE = new QLineEdit(gridLayoutWidget_2);
        line_CCL21_HALFLIFE->setObjectName(QString::fromUtf8("line_CCL21_HALFLIFE"));

        gridLayout_S1P_2->addWidget(line_CCL21_HALFLIFE, 3, 2, 1, 1);

        units_CCL21conc = new QLabel(gridLayoutWidget_2);
        units_CCL21conc->setObjectName(QString::fromUtf8("units_CCL21conc"));

        gridLayout_S1P_2->addWidget(units_CCL21conc, 1, 3, 1, 1);

        units_CCL21diff = new QLabel(gridLayoutWidget_2);
        units_CCL21diff->setObjectName(QString::fromUtf8("units_CCL21diff"));

        gridLayout_S1P_2->addWidget(units_CCL21diff, 2, 3, 1, 1);

        units_CCL21life = new QLabel(gridLayoutWidget_2);
        units_CCL21life->setObjectName(QString::fromUtf8("units_CCL21life"));

        gridLayout_S1P_2->addWidget(units_CCL21life, 3, 3, 1, 1);

        label_CCL21_STRENGTH = new QLabel(gridLayoutWidget_2);
        label_CCL21_STRENGTH->setObjectName(QString::fromUtf8("label_CCL21_STRENGTH"));

        gridLayout_S1P_2->addWidget(label_CCL21_STRENGTH, 4, 0, 1, 2);

        line_CCL21_STRENGTH = new QLineEdit(gridLayoutWidget_2);
        line_CCL21_STRENGTH->setObjectName(QString::fromUtf8("line_CCL21_STRENGTH"));

        gridLayout_S1P_2->addWidget(line_CCL21_STRENGTH, 4, 2, 1, 1);

        label_CCL21_BDRY_RATE = new QLabel(gridLayoutWidget_2);
        label_CCL21_BDRY_RATE->setObjectName(QString::fromUtf8("label_CCL21_BDRY_RATE"));

        gridLayout_S1P_2->addWidget(label_CCL21_BDRY_RATE, 0, 0, 1, 2);

        line_CCL21_BDRY_RATE = new QLineEdit(gridLayoutWidget_2);
        line_CCL21_BDRY_RATE->setObjectName(QString::fromUtf8("line_CCL21_BDRY_RATE"));

        gridLayout_S1P_2->addWidget(line_CCL21_BDRY_RATE, 0, 2, 1, 1);

        units_CCL21rate = new QLabel(gridLayoutWidget_2);
        units_CCL21rate->setObjectName(QString::fromUtf8("units_CCL21rate"));

        gridLayout_S1P_2->addWidget(units_CCL21rate, 0, 3, 1, 1);

        cbox_USE_CCR7 = new QCheckBox(tab_chemo);
        cbox_USE_CCR7->setObjectName(QString::fromUtf8("cbox_USE_CCR7"));
        cbox_USE_CCR7->setGeometry(QRect(80, 340, 81, 18));
        cbox_USE_EBI2 = new QCheckBox(tab_chemo);
        cbox_USE_EBI2->setObjectName(QString::fromUtf8("cbox_USE_EBI2"));
        cbox_USE_EBI2->setGeometry(QRect(440, 100, 71, 21));
        cbox_USE_CXCR5 = new QCheckBox(tab_chemo);
        cbox_USE_CXCR5->setObjectName(QString::fromUtf8("cbox_USE_CXCR5"));
        cbox_USE_CXCR5->setGeometry(QRect(440, 340, 81, 18));
        cbox_USE_S1PR2 = new QCheckBox(tab_chemo);
        cbox_USE_S1PR2->setObjectName(QString::fromUtf8("cbox_USE_S1PR2"));
        cbox_USE_S1PR2->setGeometry(QRect(190, 100, 81, 18));
        rbut_CCL21_BDRY_0 = new QRadioButton(tab_chemo);
        buttonGroup_CCL21_BDRY = new QButtonGroup(MainWindow);
        buttonGroup_CCL21_BDRY->setObjectName(QString::fromUtf8("buttonGroup_CCL21_BDRY"));
        buttonGroup_CCL21_BDRY->addButton(rbut_CCL21_BDRY_0);
        rbut_CCL21_BDRY_0->setObjectName(QString::fromUtf8("rbut_CCL21_BDRY_0"));
        rbut_CCL21_BDRY_0->setGeometry(QRect(300, 370, 111, 21));
        gridLayoutWidget_5 = new QWidget(tab_chemo);
        gridLayoutWidget_5->setObjectName(QString::fromUtf8("gridLayoutWidget_5"));
        gridLayoutWidget_5->setGeometry(QRect(440, 370, 220, 126));
        gridLayout_S1P_5 = new QGridLayout(gridLayoutWidget_5);
        gridLayout_S1P_5->setObjectName(QString::fromUtf8("gridLayout_S1P_5"));
        gridLayout_S1P_5->setContentsMargins(0, 0, 0, 0);
        label_CXCL13_BDRY_CONC = new QMyLabel(gridLayoutWidget_5);
        label_CXCL13_BDRY_CONC->setObjectName(QString::fromUtf8("label_CXCL13_BDRY_CONC"));

        gridLayout_S1P_5->addWidget(label_CXCL13_BDRY_CONC, 1, 0, 1, 2);

        line_CXCL13_BDRY_CONC = new QLineEdit(gridLayoutWidget_5);
        line_CXCL13_BDRY_CONC->setObjectName(QString::fromUtf8("line_CXCL13_BDRY_CONC"));

        gridLayout_S1P_5->addWidget(line_CXCL13_BDRY_CONC, 1, 2, 1, 1);

        label_CXCL13_DIFF_COEFF = new QMyLabel(gridLayoutWidget_5);
        label_CXCL13_DIFF_COEFF->setObjectName(QString::fromUtf8("label_CXCL13_DIFF_COEFF"));

        gridLayout_S1P_5->addWidget(label_CXCL13_DIFF_COEFF, 2, 0, 1, 2);

        line_CXCL13_DIFF_COEFF = new QLineEdit(gridLayoutWidget_5);
        line_CXCL13_DIFF_COEFF->setObjectName(QString::fromUtf8("line_CXCL13_DIFF_COEFF"));

        gridLayout_S1P_5->addWidget(line_CXCL13_DIFF_COEFF, 2, 2, 1, 1);

        label_CXCL13_HALFLIFE = new QMyLabel(gridLayoutWidget_5);
        label_CXCL13_HALFLIFE->setObjectName(QString::fromUtf8("label_CXCL13_HALFLIFE"));

        gridLayout_S1P_5->addWidget(label_CXCL13_HALFLIFE, 3, 0, 1, 2);

        line_CXCL13_HALFLIFE = new QLineEdit(gridLayoutWidget_5);
        line_CXCL13_HALFLIFE->setObjectName(QString::fromUtf8("line_CXCL13_HALFLIFE"));

        gridLayout_S1P_5->addWidget(line_CXCL13_HALFLIFE, 3, 2, 1, 1);

        units_CXCL13conc = new QLabel(gridLayoutWidget_5);
        units_CXCL13conc->setObjectName(QString::fromUtf8("units_CXCL13conc"));

        gridLayout_S1P_5->addWidget(units_CXCL13conc, 1, 3, 1, 1);

        units_CXCL13diff = new QLabel(gridLayoutWidget_5);
        units_CXCL13diff->setObjectName(QString::fromUtf8("units_CXCL13diff"));

        gridLayout_S1P_5->addWidget(units_CXCL13diff, 2, 3, 1, 1);

        units_CXCL13life = new QLabel(gridLayoutWidget_5);
        units_CXCL13life->setObjectName(QString::fromUtf8("units_CXCL13life"));

        gridLayout_S1P_5->addWidget(units_CXCL13life, 3, 3, 1, 1);

        label_CXCL13_STRENGTH = new QLabel(gridLayoutWidget_5);
        label_CXCL13_STRENGTH->setObjectName(QString::fromUtf8("label_CXCL13_STRENGTH"));

        gridLayout_S1P_5->addWidget(label_CXCL13_STRENGTH, 4, 0, 1, 2);

        line_CXCL13_STRENGTH = new QLineEdit(gridLayoutWidget_5);
        line_CXCL13_STRENGTH->setObjectName(QString::fromUtf8("line_CXCL13_STRENGTH"));

        gridLayout_S1P_5->addWidget(line_CXCL13_STRENGTH, 4, 2, 1, 1);

        label_CXCL13_BDRY_RATE = new QLabel(gridLayoutWidget_5);
        label_CXCL13_BDRY_RATE->setObjectName(QString::fromUtf8("label_CXCL13_BDRY_RATE"));

        gridLayout_S1P_5->addWidget(label_CXCL13_BDRY_RATE, 0, 0, 1, 2);

        line_CXCL13_BDRY_RATE = new QLineEdit(gridLayoutWidget_5);
        line_CXCL13_BDRY_RATE->setObjectName(QString::fromUtf8("line_CXCL13_BDRY_RATE"));

        gridLayout_S1P_5->addWidget(line_CXCL13_BDRY_RATE, 0, 2, 1, 1);

        units_CXCL13rate = new QLabel(gridLayoutWidget_5);
        units_CXCL13rate->setObjectName(QString::fromUtf8("units_CXCL13rate"));

        gridLayout_S1P_5->addWidget(units_CXCL13rate, 0, 3, 1, 1);

        rbut_CXCL13_BDRY_0 = new QRadioButton(tab_chemo);
        buttonGroup_CXCL13_BDRY = new QButtonGroup(MainWindow);
        buttonGroup_CXCL13_BDRY->setObjectName(QString::fromUtf8("buttonGroup_CXCL13_BDRY"));
        buttonGroup_CXCL13_BDRY->addButton(rbut_CXCL13_BDRY_0);
        rbut_CXCL13_BDRY_0->setObjectName(QString::fromUtf8("rbut_CXCL13_BDRY_0"));
        rbut_CXCL13_BDRY_0->setGeometry(QRect(670, 370, 111, 21));
        gridLayoutWidget_4 = new QWidget(tab_chemo);
        gridLayoutWidget_4->setObjectName(QString::fromUtf8("gridLayoutWidget_4"));
        gridLayoutWidget_4->setGeometry(QRect(440, 130, 221, 126));
        gridLayout_S1P_4 = new QGridLayout(gridLayoutWidget_4);
        gridLayout_S1P_4->setObjectName(QString::fromUtf8("gridLayout_S1P_4"));
        gridLayout_S1P_4->setContentsMargins(0, 0, 0, 0);
        label_OXY_BDRY_CONC = new QMyLabel(gridLayoutWidget_4);
        label_OXY_BDRY_CONC->setObjectName(QString::fromUtf8("label_OXY_BDRY_CONC"));

        gridLayout_S1P_4->addWidget(label_OXY_BDRY_CONC, 1, 0, 1, 2);

        line_OXY_BDRY_CONC = new QLineEdit(gridLayoutWidget_4);
        line_OXY_BDRY_CONC->setObjectName(QString::fromUtf8("line_OXY_BDRY_CONC"));

        gridLayout_S1P_4->addWidget(line_OXY_BDRY_CONC, 1, 2, 1, 1);

        label_OXY_DIFF_COEFF = new QMyLabel(gridLayoutWidget_4);
        label_OXY_DIFF_COEFF->setObjectName(QString::fromUtf8("label_OXY_DIFF_COEFF"));

        gridLayout_S1P_4->addWidget(label_OXY_DIFF_COEFF, 2, 0, 1, 2);

        line_OXY_DIFF_COEFF = new QLineEdit(gridLayoutWidget_4);
        line_OXY_DIFF_COEFF->setObjectName(QString::fromUtf8("line_OXY_DIFF_COEFF"));

        gridLayout_S1P_4->addWidget(line_OXY_DIFF_COEFF, 2, 2, 1, 1);

        label_OXY_HALFLIFE = new QMyLabel(gridLayoutWidget_4);
        label_OXY_HALFLIFE->setObjectName(QString::fromUtf8("label_OXY_HALFLIFE"));

        gridLayout_S1P_4->addWidget(label_OXY_HALFLIFE, 3, 0, 1, 2);

        line_OXY_HALFLIFE = new QLineEdit(gridLayoutWidget_4);
        line_OXY_HALFLIFE->setObjectName(QString::fromUtf8("line_OXY_HALFLIFE"));

        gridLayout_S1P_4->addWidget(line_OXY_HALFLIFE, 3, 2, 1, 1);

        units_OXYconc = new QLabel(gridLayoutWidget_4);
        units_OXYconc->setObjectName(QString::fromUtf8("units_OXYconc"));

        gridLayout_S1P_4->addWidget(units_OXYconc, 1, 3, 1, 1);

        units_OXYdiff = new QLabel(gridLayoutWidget_4);
        units_OXYdiff->setObjectName(QString::fromUtf8("units_OXYdiff"));

        gridLayout_S1P_4->addWidget(units_OXYdiff, 2, 3, 1, 1);

        units_OXYlife = new QLabel(gridLayoutWidget_4);
        units_OXYlife->setObjectName(QString::fromUtf8("units_OXYlife"));

        gridLayout_S1P_4->addWidget(units_OXYlife, 3, 3, 1, 1);

        label_OXY_STRENGTH = new QLabel(gridLayoutWidget_4);
        label_OXY_STRENGTH->setObjectName(QString::fromUtf8("label_OXY_STRENGTH"));

        gridLayout_S1P_4->addWidget(label_OXY_STRENGTH, 4, 0, 1, 2);

        line_OXY_STRENGTH = new QLineEdit(gridLayoutWidget_4);
        line_OXY_STRENGTH->setObjectName(QString::fromUtf8("line_OXY_STRENGTH"));

        gridLayout_S1P_4->addWidget(line_OXY_STRENGTH, 4, 2, 1, 1);

        label_OXY_BDRY_RATE = new QLabel(gridLayoutWidget_4);
        label_OXY_BDRY_RATE->setObjectName(QString::fromUtf8("label_OXY_BDRY_RATE"));

        gridLayout_S1P_4->addWidget(label_OXY_BDRY_RATE, 0, 0, 1, 2);

        line_OXY_BDRY_RATE = new QLineEdit(gridLayoutWidget_4);
        line_OXY_BDRY_RATE->setObjectName(QString::fromUtf8("line_OXY_BDRY_RATE"));

        gridLayout_S1P_4->addWidget(line_OXY_BDRY_RATE, 0, 2, 1, 1);

        units_OXYrate = new QLabel(gridLayoutWidget_4);
        units_OXYrate->setObjectName(QString::fromUtf8("units_OXYrate"));

        gridLayout_S1P_4->addWidget(units_OXYrate, 0, 3, 1, 1);

        rbut_OXY_BDRY_0 = new QRadioButton(tab_chemo);
        buttonGroup_OXY_BDRY = new QButtonGroup(MainWindow);
        buttonGroup_OXY_BDRY->setObjectName(QString::fromUtf8("buttonGroup_OXY_BDRY"));
        buttonGroup_OXY_BDRY->addButton(rbut_OXY_BDRY_0);
        rbut_OXY_BDRY_0->setObjectName(QString::fromUtf8("rbut_OXY_BDRY_0"));
        rbut_OXY_BDRY_0->setGeometry(QRect(670, 130, 111, 21));
        gridLayoutWidget_3 = new QWidget(tab_chemo);
        gridLayoutWidget_3->setObjectName(QString::fromUtf8("gridLayoutWidget_3"));
        gridLayoutWidget_3->setGeometry(QRect(80, 130, 211, 152));
        gridLayout_S1P_3 = new QGridLayout(gridLayoutWidget_3);
        gridLayout_S1P_3->setObjectName(QString::fromUtf8("gridLayout_S1P_3"));
        gridLayout_S1P_3->setContentsMargins(0, 0, 0, 0);
        label_S1P_BDRY_CONC = new QMyLabel(gridLayoutWidget_3);
        label_S1P_BDRY_CONC->setObjectName(QString::fromUtf8("label_S1P_BDRY_CONC"));

        gridLayout_S1P_3->addWidget(label_S1P_BDRY_CONC, 1, 0, 1, 2);

        line_S1P_BDRY_CONC = new QLineEdit(gridLayoutWidget_3);
        line_S1P_BDRY_CONC->setObjectName(QString::fromUtf8("line_S1P_BDRY_CONC"));

        gridLayout_S1P_3->addWidget(line_S1P_BDRY_CONC, 1, 2, 1, 1);

        label_S1P_DIFF_COEFF = new QMyLabel(gridLayoutWidget_3);
        label_S1P_DIFF_COEFF->setObjectName(QString::fromUtf8("label_S1P_DIFF_COEFF"));

        gridLayout_S1P_3->addWidget(label_S1P_DIFF_COEFF, 2, 0, 1, 2);

        line_S1P_DIFF_COEFF = new QLineEdit(gridLayoutWidget_3);
        line_S1P_DIFF_COEFF->setObjectName(QString::fromUtf8("line_S1P_DIFF_COEFF"));

        gridLayout_S1P_3->addWidget(line_S1P_DIFF_COEFF, 2, 2, 1, 1);

        label_S1P_HALFLIFE = new QMyLabel(gridLayoutWidget_3);
        label_S1P_HALFLIFE->setObjectName(QString::fromUtf8("label_S1P_HALFLIFE"));

        gridLayout_S1P_3->addWidget(label_S1P_HALFLIFE, 3, 0, 1, 2);

        line_S1P_HALFLIFE = new QLineEdit(gridLayoutWidget_3);
        line_S1P_HALFLIFE->setObjectName(QString::fromUtf8("line_S1P_HALFLIFE"));

        gridLayout_S1P_3->addWidget(line_S1P_HALFLIFE, 3, 2, 1, 1);

        units_S1Pconc = new QLabel(gridLayoutWidget_3);
        units_S1Pconc->setObjectName(QString::fromUtf8("units_S1Pconc"));

        gridLayout_S1P_3->addWidget(units_S1Pconc, 1, 3, 1, 1);

        units_S1Pdiff = new QLabel(gridLayoutWidget_3);
        units_S1Pdiff->setObjectName(QString::fromUtf8("units_S1Pdiff"));

        gridLayout_S1P_3->addWidget(units_S1Pdiff, 2, 3, 1, 1);

        units_S1Plife = new QLabel(gridLayoutWidget_3);
        units_S1Plife->setObjectName(QString::fromUtf8("units_S1Plife"));

        gridLayout_S1P_3->addWidget(units_S1Plife, 3, 3, 1, 1);

        label_S1P_STRENGTH_POS = new QLabel(gridLayoutWidget_3);
        label_S1P_STRENGTH_POS->setObjectName(QString::fromUtf8("label_S1P_STRENGTH_POS"));

        gridLayout_S1P_3->addWidget(label_S1P_STRENGTH_POS, 4, 0, 1, 2);

        line_S1P_STRENGTH_POS = new QLineEdit(gridLayoutWidget_3);
        line_S1P_STRENGTH_POS->setObjectName(QString::fromUtf8("line_S1P_STRENGTH_POS"));

        gridLayout_S1P_3->addWidget(line_S1P_STRENGTH_POS, 4, 2, 1, 1);

        label_S1P_STRENGTH_NEG = new QLabel(gridLayoutWidget_3);
        label_S1P_STRENGTH_NEG->setObjectName(QString::fromUtf8("label_S1P_STRENGTH_NEG"));

        gridLayout_S1P_3->addWidget(label_S1P_STRENGTH_NEG, 5, 0, 1, 1);

        line_S1P_STRENGTH_NEG = new QLineEdit(gridLayoutWidget_3);
        line_S1P_STRENGTH_NEG->setObjectName(QString::fromUtf8("line_S1P_STRENGTH_NEG"));

        gridLayout_S1P_3->addWidget(line_S1P_STRENGTH_NEG, 5, 2, 1, 1);

        label_S1P_BDRY_RATE = new QLabel(gridLayoutWidget_3);
        label_S1P_BDRY_RATE->setObjectName(QString::fromUtf8("label_S1P_BDRY_RATE"));

        gridLayout_S1P_3->addWidget(label_S1P_BDRY_RATE, 0, 0, 1, 2);

        line_S1P_BDRY_RATE = new QLineEdit(gridLayoutWidget_3);
        line_S1P_BDRY_RATE->setObjectName(QString::fromUtf8("line_S1P_BDRY_RATE"));

        gridLayout_S1P_3->addWidget(line_S1P_BDRY_RATE, 0, 2, 1, 1);

        units_S1Prate = new QLabel(gridLayoutWidget_3);
        units_S1Prate->setObjectName(QString::fromUtf8("units_S1Prate"));

        gridLayout_S1P_3->addWidget(units_S1Prate, 0, 3, 1, 1);

        rbut_S1P_BDRY_0 = new QRadioButton(tab_chemo);
        buttonGroup_S1P_BDRY = new QButtonGroup(MainWindow);
        buttonGroup_S1P_BDRY->setObjectName(QString::fromUtf8("buttonGroup_S1P_BDRY"));
        buttonGroup_S1P_BDRY->addButton(rbut_S1P_BDRY_0);
        rbut_S1P_BDRY_0->setObjectName(QString::fromUtf8("rbut_S1P_BDRY_0"));
        rbut_S1P_BDRY_0->setGeometry(QRect(300, 130, 91, 21));
        label_chemokine = new QLabel(tab_chemo);
        label_chemokine->setObjectName(QString::fromUtf8("label_chemokine"));
        label_chemokine->setGeometry(QRect(280, 30, 221, 16));
        QFont font2;
        font2.setPointSize(16);
        label_chemokine->setFont(font2);
        rbut_S1P_BDRY_1 = new QRadioButton(tab_chemo);
        buttonGroup_S1P_BDRY->addButton(rbut_S1P_BDRY_1);
        rbut_S1P_BDRY_1->setObjectName(QString::fromUtf8("rbut_S1P_BDRY_1"));
        rbut_S1P_BDRY_1->setGeometry(QRect(300, 150, 111, 21));
        rbut_S1P_BDRY_1->setChecked(true);
        rbut_CCL21_BDRY_1 = new QRadioButton(tab_chemo);
        buttonGroup_CCL21_BDRY->addButton(rbut_CCL21_BDRY_1);
        rbut_CCL21_BDRY_1->setObjectName(QString::fromUtf8("rbut_CCL21_BDRY_1"));
        rbut_CCL21_BDRY_1->setGeometry(QRect(300, 390, 111, 21));
        rbut_CCL21_BDRY_1->setChecked(true);
        rbut_OXY_BDRY_1 = new QRadioButton(tab_chemo);
        buttonGroup_OXY_BDRY->addButton(rbut_OXY_BDRY_1);
        rbut_OXY_BDRY_1->setObjectName(QString::fromUtf8("rbut_OXY_BDRY_1"));
        rbut_OXY_BDRY_1->setGeometry(QRect(670, 150, 111, 21));
        rbut_OXY_BDRY_1->setChecked(true);
        rbut_CXCL13_BDRY_1 = new QRadioButton(tab_chemo);
        buttonGroup_CXCL13_BDRY->addButton(rbut_CXCL13_BDRY_1);
        rbut_CXCL13_BDRY_1->setObjectName(QString::fromUtf8("rbut_CXCL13_BDRY_1"));
        rbut_CXCL13_BDRY_1->setGeometry(QRect(670, 390, 111, 21));
        rbut_CXCL13_BDRY_1->setChecked(true);
        tabs->addTab(tab_chemo, QString());
        tab_TCR = new QWidget();
        tab_TCR->setObjectName(QString::fromUtf8("tab_TCR"));
        tab_TCR->setEnabled(false);
        layoutWidget_2 = new QWidget(tab_TCR);
        layoutWidget_2->setObjectName(QString::fromUtf8("layoutWidget_2"));
        layoutWidget_2->setGeometry(QRect(20, 40, 421, 191));
        gridLayout_3 = new QGridLayout(layoutWidget_2);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setContentsMargins(0, 0, 0, 0);
        label_IL2_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_IL2_THRESHOLD->setObjectName(QString::fromUtf8("label_IL2_THRESHOLD"));

        gridLayout_3->addWidget(label_IL2_THRESHOLD, 0, 0, 1, 1);

        label_ACTIVATION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_ACTIVATION_THRESHOLD->setObjectName(QString::fromUtf8("label_ACTIVATION_THRESHOLD"));

        gridLayout_3->addWidget(label_ACTIVATION_THRESHOLD, 1, 0, 1, 1);

        label_FIRST_DIVISION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_FIRST_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("label_FIRST_DIVISION_THRESHOLD"));

        gridLayout_3->addWidget(label_FIRST_DIVISION_THRESHOLD, 2, 0, 1, 1);

        label_DIVISION_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("label_DIVISION_THRESHOLD"));
        QSizePolicy sizePolicy2(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        label_DIVISION_THRESHOLD->setSizePolicy(sizePolicy2);
        label_DIVISION_THRESHOLD->setMouseTracking(false);
        label_DIVISION_THRESHOLD->setWordWrap(false);

        gridLayout_3->addWidget(label_DIVISION_THRESHOLD, 3, 0, 1, 1);

        label_EXIT_THRESHOLD = new QMyLabel(layoutWidget_2);
        label_EXIT_THRESHOLD->setObjectName(QString::fromUtf8("label_EXIT_THRESHOLD"));
        label_EXIT_THRESHOLD->setMouseTracking(false);
        label_EXIT_THRESHOLD->setWordWrap(false);

        gridLayout_3->addWidget(label_EXIT_THRESHOLD, 4, 0, 1, 1);

        line_IL2_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_IL2_THRESHOLD->setObjectName(QString::fromUtf8("line_IL2_THRESHOLD"));
        QSizePolicy sizePolicy3(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(120);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(line_IL2_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_IL2_THRESHOLD->setSizePolicy(sizePolicy3);
        line_IL2_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_IL2_THRESHOLD, 0, 1, 1, 1);

        line_ACTIVATION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_ACTIVATION_THRESHOLD->setObjectName(QString::fromUtf8("line_ACTIVATION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_ACTIVATION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_ACTIVATION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_ACTIVATION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_ACTIVATION_THRESHOLD, 1, 1, 1, 1);

        line_FIRST_DIVISION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_FIRST_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("line_FIRST_DIVISION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_FIRST_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_FIRST_DIVISION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_FIRST_DIVISION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_FIRST_DIVISION_THRESHOLD, 2, 1, 1, 1);

        line_DIVISION_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_DIVISION_THRESHOLD->setObjectName(QString::fromUtf8("line_DIVISION_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_DIVISION_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_DIVISION_THRESHOLD->setSizePolicy(sizePolicy3);
        line_DIVISION_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_DIVISION_THRESHOLD, 3, 1, 1, 1);

        line_EXIT_THRESHOLD = new QLineEdit(layoutWidget_2);
        line_EXIT_THRESHOLD->setObjectName(QString::fromUtf8("line_EXIT_THRESHOLD"));
        sizePolicy3.setHeightForWidth(line_EXIT_THRESHOLD->sizePolicy().hasHeightForWidth());
        line_EXIT_THRESHOLD->setSizePolicy(sizePolicy3);
        line_EXIT_THRESHOLD->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_EXIT_THRESHOLD, 4, 1, 1, 1);

        label_STIMULATION_LIMIT = new QMyLabel(layoutWidget_2);
        label_STIMULATION_LIMIT->setObjectName(QString::fromUtf8("label_STIMULATION_LIMIT"));
        label_STIMULATION_LIMIT->setMouseTracking(false);
        label_STIMULATION_LIMIT->setWordWrap(false);

        gridLayout_3->addWidget(label_STIMULATION_LIMIT, 5, 0, 1, 1);

        line_STIMULATION_LIMIT = new QLineEdit(layoutWidget_2);
        line_STIMULATION_LIMIT->setObjectName(QString::fromUtf8("line_STIMULATION_LIMIT"));
        sizePolicy3.setHeightForWidth(line_STIMULATION_LIMIT->sizePolicy().hasHeightForWidth());
        line_STIMULATION_LIMIT->setSizePolicy(sizePolicy3);
        line_STIMULATION_LIMIT->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(line_STIMULATION_LIMIT, 5, 1, 1, 1);

        label_40 = new QLabel(tab_TCR);
        label_40->setObjectName(QString::fromUtf8("label_40"));
        label_40->setGeometry(QRect(10, 10, 211, 16));
        QFont font3;
        font3.setBold(true);
        font3.setWeight(75);
        label_40->setFont(font3);
        tabs->addTab(tab_TCR, QString());
        tab_run = new QWidget();
        tab_run->setObjectName(QString::fromUtf8("tab_run"));
        layoutWidget2 = new QWidget(tab_run);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(10, 0, 431, 668));
        gridLayout_4 = new QGridLayout(layoutWidget2);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_4->setContentsMargins(0, 0, 0, 0);
        label_NX = new QMyLabel(layoutWidget2);
        label_NX->setObjectName(QString::fromUtf8("label_NX"));
        label_NX->setWordWrap(false);

        gridLayout_4->addWidget(label_NX, 0, 0, 1, 1);

        horizontalSpacer_12 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_12, 0, 1, 1, 1);

        spin_NX = new QSpinBox(layoutWidget2);
        spin_NX->setObjectName(QString::fromUtf8("spin_NX"));
        spin_NX->setMaximumSize(QSize(120, 16777215));
        spin_NX->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NX->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NX, 0, 2, 1, 1);

        label_BLOB_RADIUS = new QMyLabel(layoutWidget2);
        label_BLOB_RADIUS->setObjectName(QString::fromUtf8("label_BLOB_RADIUS"));
        label_BLOB_RADIUS->setWordWrap(false);

        gridLayout_4->addWidget(label_BLOB_RADIUS, 1, 0, 1, 1);

        horizontalSpacer_13 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_13, 1, 1, 1, 1);

        line_BLOB_RADIUS = new QLineEdit(layoutWidget2);
        line_BLOB_RADIUS->setObjectName(QString::fromUtf8("line_BLOB_RADIUS"));
        line_BLOB_RADIUS->setMaximumSize(QSize(120, 16777215));
        line_BLOB_RADIUS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_BLOB_RADIUS, 1, 2, 1, 1);

        label_BC_FRACTION = new QMyLabel(layoutWidget2);
        label_BC_FRACTION->setObjectName(QString::fromUtf8("label_BC_FRACTION"));
        label_BC_FRACTION->setWordWrap(false);

        gridLayout_4->addWidget(label_BC_FRACTION, 2, 0, 1, 1);

        horizontalSpacer_14 = new QSpacerItem(38, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_14, 2, 1, 1, 1);

        line_BC_FRACTION = new QLineEdit(layoutWidget2);
        line_BC_FRACTION->setObjectName(QString::fromUtf8("line_BC_FRACTION"));
        line_BC_FRACTION->setMaximumSize(QSize(120, 16777215));
        line_BC_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_BC_FRACTION, 2, 2, 1, 1);

        horizontalSpacer_15 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_15, 3, 1, 1, 1);

        line_FLUID_FRACTION = new QLineEdit(layoutWidget2);
        line_FLUID_FRACTION->setObjectName(QString::fromUtf8("line_FLUID_FRACTION"));
        line_FLUID_FRACTION->setMaximumSize(QSize(120, 16777215));
        line_FLUID_FRACTION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_FLUID_FRACTION, 3, 2, 1, 1);

        label_RESIDENCE_TIME = new QMyLabel(layoutWidget2);
        label_RESIDENCE_TIME->setObjectName(QString::fromUtf8("label_RESIDENCE_TIME"));
        label_RESIDENCE_TIME->setMouseTracking(false);
        label_RESIDENCE_TIME->setWordWrap(false);

        gridLayout_4->addWidget(label_RESIDENCE_TIME, 4, 0, 1, 1);

        horizontalSpacer_21 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_21, 4, 1, 1, 1);

        line_RESIDENCE_TIME = new QLineEdit(layoutWidget2);
        line_RESIDENCE_TIME->setObjectName(QString::fromUtf8("line_RESIDENCE_TIME"));
        line_RESIDENCE_TIME->setMaximumSize(QSize(120, 16777215));
        line_RESIDENCE_TIME->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_RESIDENCE_TIME, 4, 2, 1, 1);

        units_RESIDENCE_TIME = new QLabel(layoutWidget2);
        units_RESIDENCE_TIME->setObjectName(QString::fromUtf8("units_RESIDENCE_TIME"));

        gridLayout_4->addWidget(units_RESIDENCE_TIME, 4, 3, 1, 1);

        label_INFLAMM_DAYS1 = new QMyLabel(layoutWidget2);
        label_INFLAMM_DAYS1->setObjectName(QString::fromUtf8("label_INFLAMM_DAYS1"));

        gridLayout_4->addWidget(label_INFLAMM_DAYS1, 5, 0, 1, 1);

        horizontalSpacer_22 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_22, 5, 1, 1, 1);

        line_INFLAMM_DAYS1 = new QLineEdit(layoutWidget2);
        line_INFLAMM_DAYS1->setObjectName(QString::fromUtf8("line_INFLAMM_DAYS1"));
        line_INFLAMM_DAYS1->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_DAYS1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_DAYS1, 5, 2, 1, 1);

        label_INFLAMM_DAYS2 = new QMyLabel(layoutWidget2);
        label_INFLAMM_DAYS2->setObjectName(QString::fromUtf8("label_INFLAMM_DAYS2"));

        gridLayout_4->addWidget(label_INFLAMM_DAYS2, 6, 0, 1, 1);

        horizontalSpacer_23 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_23, 6, 1, 1, 1);

        line_INFLAMM_DAYS2 = new QLineEdit(layoutWidget2);
        line_INFLAMM_DAYS2->setObjectName(QString::fromUtf8("line_INFLAMM_DAYS2"));
        line_INFLAMM_DAYS2->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_DAYS2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_DAYS2, 6, 2, 1, 1);

        label_INFLAMM_LEVEL = new QMyLabel(layoutWidget2);
        label_INFLAMM_LEVEL->setObjectName(QString::fromUtf8("label_INFLAMM_LEVEL"));

        gridLayout_4->addWidget(label_INFLAMM_LEVEL, 7, 0, 1, 1);

        horizontalSpacer_24 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_24, 7, 1, 1, 1);

        line_INFLAMM_LEVEL = new QLineEdit(layoutWidget2);
        line_INFLAMM_LEVEL->setObjectName(QString::fromUtf8("line_INFLAMM_LEVEL"));
        line_INFLAMM_LEVEL->setMaximumSize(QSize(120, 16777215));
        line_INFLAMM_LEVEL->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_INFLAMM_LEVEL, 7, 2, 1, 1);

        label_CHEMO_RADIUS = new QMyLabel(layoutWidget2);
        label_CHEMO_RADIUS->setObjectName(QString::fromUtf8("label_CHEMO_RADIUS"));
        label_CHEMO_RADIUS->setEnabled(false);

        gridLayout_4->addWidget(label_CHEMO_RADIUS, 8, 0, 1, 1);

        horizontalSpacer_27 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_27, 8, 1, 1, 1);

        line_CHEMO_RADIUS = new QLineEdit(layoutWidget2);
        line_CHEMO_RADIUS->setObjectName(QString::fromUtf8("line_CHEMO_RADIUS"));
        line_CHEMO_RADIUS->setEnabled(false);
        line_CHEMO_RADIUS->setMaximumSize(QSize(120, 16777215));
        line_CHEMO_RADIUS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_CHEMO_RADIUS, 8, 2, 1, 1);

        label_BASE_EXIT_PROB = new QMyLabel(layoutWidget2);
        label_BASE_EXIT_PROB->setObjectName(QString::fromUtf8("label_BASE_EXIT_PROB"));
        label_BASE_EXIT_PROB->setMouseTracking(false);
        label_BASE_EXIT_PROB->setWordWrap(false);

        gridLayout_4->addWidget(label_BASE_EXIT_PROB, 9, 0, 1, 1);

        horizontalSpacer_28 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_28, 9, 1, 1, 1);

        line_BASE_EXIT_PROB = new QLineEdit(layoutWidget2);
        line_BASE_EXIT_PROB->setObjectName(QString::fromUtf8("line_BASE_EXIT_PROB"));
        line_BASE_EXIT_PROB->setMaximumSize(QSize(120, 16777215));
        line_BASE_EXIT_PROB->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_BASE_EXIT_PROB, 9, 2, 1, 1);

        label_NDAYS = new QMyLabel(layoutWidget2);
        label_NDAYS->setObjectName(QString::fromUtf8("label_NDAYS"));
        label_NDAYS->setMouseTracking(false);
        label_NDAYS->setWordWrap(false);

        gridLayout_4->addWidget(label_NDAYS, 10, 0, 1, 1);

        horizontalSpacer_29 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_29, 10, 1, 1, 1);

        label_SEED1 = new QMyLabel(layoutWidget2);
        label_SEED1->setObjectName(QString::fromUtf8("label_SEED1"));
        label_SEED1->setMouseTracking(false);
        label_SEED1->setWordWrap(false);

        gridLayout_4->addWidget(label_SEED1, 11, 0, 1, 1);

        horizontalSpacer_30 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_30, 11, 1, 1, 1);

        spin_SEED1 = new QSpinBox(layoutWidget2);
        spin_SEED1->setObjectName(QString::fromUtf8("spin_SEED1"));
        spin_SEED1->setMaximumSize(QSize(120, 16777215));
        spin_SEED1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_SEED1->setMaximum(999999999);

        gridLayout_4->addWidget(spin_SEED1, 11, 2, 1, 1);

        label_SEED2 = new QMyLabel(layoutWidget2);
        label_SEED2->setObjectName(QString::fromUtf8("label_SEED2"));
        label_SEED2->setMouseTracking(false);
        label_SEED2->setWordWrap(false);

        gridLayout_4->addWidget(label_SEED2, 12, 0, 1, 1);

        horizontalSpacer_31 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_31, 12, 1, 1, 1);

        spin_SEED2 = new QSpinBox(layoutWidget2);
        spin_SEED2->setObjectName(QString::fromUtf8("spin_SEED2"));
        spin_SEED2->setMaximumSize(QSize(120, 16777215));
        spin_SEED2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_SEED2->setMaximum(999999999);

        gridLayout_4->addWidget(spin_SEED2, 12, 2, 1, 1);

        line_NDAYS = new QLineEdit(layoutWidget2);
        line_NDAYS->setObjectName(QString::fromUtf8("line_NDAYS"));
        line_NDAYS->setMaximumSize(QSize(120, 16777215));
        line_NDAYS->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(line_NDAYS, 10, 2, 1, 1);

        label_FLUID_FRACTION = new QMyLabel(layoutWidget2);
        label_FLUID_FRACTION->setObjectName(QString::fromUtf8("label_FLUID_FRACTION"));
        label_FLUID_FRACTION->setMouseTracking(false);
        label_FLUID_FRACTION->setWordWrap(false);

        gridLayout_4->addWidget(label_FLUID_FRACTION, 3, 0, 1, 1);

        label_NT_ANIMATION = new QMyLabel(layoutWidget2);
        label_NT_ANIMATION->setObjectName(QString::fromUtf8("label_NT_ANIMATION"));

        gridLayout_4->addWidget(label_NT_ANIMATION, 14, 0, 1, 1);

        horizontalSpacer_32 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_32, 14, 1, 1, 1);

        spin_NT_ANIMATION = new QSpinBox(layoutWidget2);
        spin_NT_ANIMATION->setObjectName(QString::fromUtf8("spin_NT_ANIMATION"));
        spin_NT_ANIMATION->setMaximumSize(QSize(120, 16777215));
        spin_NT_ANIMATION->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NT_ANIMATION->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NT_ANIMATION, 14, 2, 1, 1);

        units_CHEMO_RADIUS = new QLabel(layoutWidget2);
        units_CHEMO_RADIUS->setObjectName(QString::fromUtf8("units_CHEMO_RADIUS"));

        gridLayout_4->addWidget(units_CHEMO_RADIUS, 8, 3, 1, 1);

        units_NDAYS = new QLabel(layoutWidget2);
        units_NDAYS->setObjectName(QString::fromUtf8("units_NDAYS"));

        gridLayout_4->addWidget(units_NDAYS, 10, 3, 1, 1);

        units_INFLAMM1 = new QLabel(layoutWidget2);
        units_INFLAMM1->setObjectName(QString::fromUtf8("units_INFLAMM1"));

        gridLayout_4->addWidget(units_INFLAMM1, 5, 3, 1, 1);

        units_INFLAMM2 = new QLabel(layoutWidget2);
        units_INFLAMM2->setObjectName(QString::fromUtf8("units_INFLAMM2"));

        gridLayout_4->addWidget(units_INFLAMM2, 6, 3, 1, 1);

        units_BLOB_RADIUS = new QLabel(layoutWidget2);
        units_BLOB_RADIUS->setObjectName(QString::fromUtf8("units_BLOB_RADIUS"));

        gridLayout_4->addWidget(units_BLOB_RADIUS, 1, 3, 1, 1);

        label_NCPU = new QMyLabel(layoutWidget2);
        label_NCPU->setObjectName(QString::fromUtf8("label_NCPU"));

        gridLayout_4->addWidget(label_NCPU, 13, 0, 1, 1);

        horizontalSpacer_34 = new QSpacerItem(40, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_34, 13, 1, 1, 1);

        spin_NCPU = new QSpinBox(layoutWidget2);
        spin_NCPU->setObjectName(QString::fromUtf8("spin_NCPU"));
        spin_NCPU->setMaximumSize(QSize(120, 16777215));
        spin_NCPU->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spin_NCPU->setMaximum(999999999);

        gridLayout_4->addWidget(spin_NCPU, 13, 2, 1, 1);

        cbox_savepos = new QCheckBox(layoutWidget2);
        cbox_savepos->setObjectName(QString::fromUtf8("cbox_savepos"));

        gridLayout_4->addWidget(cbox_savepos, 15, 2, 1, 1);

        cbox_IV_SHOW_NONCOGNATE = new QCheckBox(tab_run);
        cbox_IV_SHOW_NONCOGNATE->setObjectName(QString::fromUtf8("cbox_IV_SHOW_NONCOGNATE"));
        cbox_IV_SHOW_NONCOGNATE->setGeometry(QRect(450, 600, 161, 18));
        cbox_USE_TRAFFIC = new QCheckBox(tab_run);
        cbox_USE_TRAFFIC->setObjectName(QString::fromUtf8("cbox_USE_TRAFFIC"));
        cbox_USE_TRAFFIC->setEnabled(false);
        cbox_USE_TRAFFIC->setGeometry(QRect(450, 290, 101, 18));
        rbut_SPECIES_1 = new QRadioButton(tab_run);
        buttonGroup_SPECIES = new QButtonGroup(MainWindow);
        buttonGroup_SPECIES->setObjectName(QString::fromUtf8("buttonGroup_SPECIES"));
        buttonGroup_SPECIES->addButton(rbut_SPECIES_1);
        rbut_SPECIES_1->setObjectName(QString::fromUtf8("rbut_SPECIES_1"));
        rbut_SPECIES_1->setGeometry(QRect(470, 30, 61, 18));
        rbut_SPECIES_0 = new QRadioButton(tab_run);
        buttonGroup_SPECIES->addButton(rbut_SPECIES_0);
        rbut_SPECIES_0->setObjectName(QString::fromUtf8("rbut_SPECIES_0"));
        rbut_SPECIES_0->setGeometry(QRect(470, 10, 71, 18));
        rbut_SPECIES_0->setChecked(true);
        cbox_USE_EXIT_CHEMOTAXIS = new QCheckBox(tab_run);
        cbox_USE_EXIT_CHEMOTAXIS->setObjectName(QString::fromUtf8("cbox_USE_EXIT_CHEMOTAXIS"));
        cbox_USE_EXIT_CHEMOTAXIS->setEnabled(false);
        cbox_USE_EXIT_CHEMOTAXIS->setGeometry(QRect(450, 450, 131, 18));
        cbox_COMPUTED_OUTFLOW = new QCheckBox(tab_run);
        cbox_COMPUTED_OUTFLOW->setObjectName(QString::fromUtf8("cbox_COMPUTED_OUTFLOW"));
        cbox_COMPUTED_OUTFLOW->setEnabled(false);
        cbox_COMPUTED_OUTFLOW->setGeometry(QRect(450, 400, 141, 18));
        label_INPUT_FILE = new QMyLabel(tab_run);
        label_INPUT_FILE->setObjectName(QString::fromUtf8("label_INPUT_FILE"));
        label_INPUT_FILE->setGeometry(QRect(10, 680, 121, 16));
        text_INPUT_FILE = new QLineEdit(tab_run);
        text_INPUT_FILE->setObjectName(QString::fromUtf8("text_INPUT_FILE"));
        text_INPUT_FILE->setGeometry(QRect(130, 680, 151, 20));
        tabs->addTab(tab_run, QString());

        verticalLayout->addWidget(tabs);

        label_input = new QLabel(page_input);
        label_input->setObjectName(QString::fromUtf8("label_input"));
        label_input->setMaximumSize(QSize(16777215, 30));
        QFont font4;
        font4.setFamily(QString::fromUtf8("Segoe UI"));
        font4.setPointSize(16);
        font4.setBold(true);
        font4.setWeight(75);
        label_input->setFont(font4);
        label_input->setAlignment(Qt::AlignCenter);

        verticalLayout->addWidget(label_input);

        text_more = new QTextEdit(page_input);
        text_more->setObjectName(QString::fromUtf8("text_more"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(text_more->sizePolicy().hasHeightForWidth());
        text_more->setSizePolicy(sizePolicy4);
        text_more->setMinimumSize(QSize(480, 0));
        text_more->setMaximumSize(QSize(16777215, 100));
        text_more->setFrameShape(QFrame::StyledPanel);
        text_more->setFrameShadow(QFrame::Sunken);
        text_more->setReadOnly(true);

        verticalLayout->addWidget(text_more);

        stackedWidget->addWidget(page_input);
        page_output = new QWidget();
        page_output->setObjectName(QString::fromUtf8("page_output"));
        verticalLayout_2 = new QVBoxLayout(page_output);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        mdiArea = new QMdiArea(page_output);
        mdiArea->setObjectName(QString::fromUtf8("mdiArea"));
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(5);
        sizePolicy5.setHeightForWidth(mdiArea->sizePolicy().hasHeightForWidth());
        mdiArea->setSizePolicy(sizePolicy5);
        mdiArea->setMinimumSize(QSize(451, 541));
        mdiArea->setSizeIncrement(QSize(1, 1));
        mdiArea->setBaseSize(QSize(1, 1));
        mdiArea->setFrameShape(QFrame::NoFrame);
        mdiArea->setFrameShadow(QFrame::Plain);
        mdiArea->setTabPosition(QTabWidget::North);

        verticalLayout_2->addWidget(mdiArea);

        box_outputLog = new QTextBrowser(page_output);
        box_outputLog->setObjectName(QString::fromUtf8("box_outputLog"));
        sizePolicy4.setHeightForWidth(box_outputLog->sizePolicy().hasHeightForWidth());
        box_outputLog->setSizePolicy(sizePolicy4);
        box_outputLog->setMinimumSize(QSize(451, 0));
        box_outputLog->setMaximumSize(QSize(16777215, 200));

        verticalLayout_2->addWidget(box_outputLog);

        stackedWidget->addWidget(page_output);
        page_3D = new QWidget();
        page_3D->setObjectName(QString::fromUtf8("page_3D"));
        mdiArea_VTK = new QMdiArea(page_3D);
        mdiArea_VTK->setObjectName(QString::fromUtf8("mdiArea_VTK"));
        mdiArea_VTK->setGeometry(QRect(50, 30, 900, 900));
        progressBar = new QProgressBar(page_3D);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 30, 21, 391));
        progressBar->setValue(24);
        progressBar->setOrientation(Qt::Vertical);
        label_hour = new QLabel(page_3D);
        label_hour->setObjectName(QString::fromUtf8("label_hour"));
        label_hour->setGeometry(QRect(10, 0, 46, 20));
        QFont font5;
        font5.setPointSize(10);
        label_hour->setFont(font5);
        hour_display = new QLabel(page_3D);
        hour_display->setObjectName(QString::fromUtf8("hour_display"));
        hour_display->setGeometry(QRect(50, 0, 46, 20));
        hour_display->setFont(font5);
        stackedWidget->addWidget(page_3D);

        gridLayout_5->addWidget(stackedWidget, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1319, 18));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuEdit = new QMenu(menubar);
        menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
        menuABM = new QMenu(menubar);
        menuABM->setObjectName(QString::fromUtf8("menuABM"));
        menuGraphs = new QMenu(menubar);
        menuGraphs->setObjectName(QString::fromUtf8("menuGraphs"));
        menuPlayer = new QMenu(menubar);
        menuPlayer->setObjectName(QString::fromUtf8("menuPlayer"));
        menuSnapshot = new QMenu(menubar);
        menuSnapshot->setObjectName(QString::fromUtf8("menuSnapshot"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
        toolBar1 = new QToolBar(MainWindow);
        toolBar1->setObjectName(QString::fromUtf8("toolBar1"));
        toolBar1->setEnabled(true);
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar1);
        QWidget::setTabOrder(line_BC_COGNATE_FRACTION, line_BC_STIM_RATE_CONSTANT);
        QWidget::setTabOrder(line_BC_STIM_RATE_CONSTANT, line_BC_STIM_HALFLIFE);
        QWidget::setTabOrder(line_BC_STIM_HALFLIFE, line_MOTILITY_BETA);
        QWidget::setTabOrder(line_MOTILITY_BETA, line_MOTILITY_RHO);
        QWidget::setTabOrder(line_MOTILITY_RHO, line_DC_BIND_DELAY);
        QWidget::setTabOrder(line_DC_BIND_DELAY, line_DC_DENS_HALFLIFE);
        QWidget::setTabOrder(line_DC_DENS_HALFLIFE, spin_MAX_TC_BIND);
        QWidget::setTabOrder(spin_MAX_TC_BIND, spin_MAX_COG_BIND);
        QWidget::setTabOrder(spin_MAX_COG_BIND, spin_NX);
        QWidget::setTabOrder(spin_NX, line_BLOB_RADIUS);
        QWidget::setTabOrder(line_BLOB_RADIUS, line_BC_FRACTION);
        QWidget::setTabOrder(line_BC_FRACTION, line_FLUID_FRACTION);
        QWidget::setTabOrder(line_FLUID_FRACTION, line_RESIDENCE_TIME);
        QWidget::setTabOrder(line_RESIDENCE_TIME, line_INFLAMM_DAYS1);
        QWidget::setTabOrder(line_INFLAMM_DAYS1, line_INFLAMM_DAYS2);
        QWidget::setTabOrder(line_INFLAMM_DAYS2, line_INFLAMM_LEVEL);
        QWidget::setTabOrder(line_INFLAMM_LEVEL, line_CHEMO_RADIUS);
        QWidget::setTabOrder(line_CHEMO_RADIUS, line_BASE_EXIT_PROB);
        QWidget::setTabOrder(line_BASE_EXIT_PROB, line_NDAYS);
        QWidget::setTabOrder(line_NDAYS, spin_SEED1);
        QWidget::setTabOrder(spin_SEED1, spin_SEED2);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuEdit->menuAction());
        menubar->addAction(menuABM->menuAction());
        menubar->addAction(menuGraphs->menuAction());
        menubar->addAction(menuPlayer->menuAction());
        menubar->addAction(menuSnapshot->menuAction());
        menuFile->addAction(action_open_input);
        menuFile->addAction(action_load_results);
        menuFile->addSeparator();
        menuFile->addAction(action_save);
        menuFile->addAction(action_saveAs);
        menuABM->addAction(action_run);
        menuABM->addAction(action_pause);
        menuABM->addAction(action_stop);
        menuGraphs->addAction(action_add_graph);
        menuGraphs->addAction(action_remove_graph);
        menuGraphs->addAction(action_remove_all);
        menuPlayer->addAction(action_play_VTK);
        menuPlayer->addAction(action_set_speed);
        menuSnapshot->addAction(action_save_snapshot);
        toolBar1->addAction(action_run);
        toolBar1->addAction(action_pause);
        toolBar1->addAction(action_stop);
        toolBar1->addSeparator();
        toolBar1->addAction(action_inputs);
        toolBar1->addAction(action_outputs);
        toolBar1->addAction(action_VTK);

        retranslateUi(MainWindow);

        stackedWidget->setCurrentIndex(0);
        tabs->setCurrentIndex(2);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "B Cell Follicle ABM", 0, QApplication::UnicodeUTF8));
        action_saveAs->setText(QApplication::translate("MainWindow", "Save &As", 0, QApplication::UnicodeUTF8));
        action_save->setText(QApplication::translate("MainWindow", "&Save", 0, QApplication::UnicodeUTF8));
        action_open_input->setText(QApplication::translate("MainWindow", "&Open input", 0, QApplication::UnicodeUTF8));
        action_open_input->setIconText(QApplication::translate("MainWindow", "Open input", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_open_input->setToolTip(QApplication::translate("MainWindow", "Open input file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_stop->setText(QApplication::translate("MainWindow", "Stop", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_stop->setToolTip(QApplication::translate("MainWindow", "Stop ABM", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_stop->setShortcut(QApplication::translate("MainWindow", "F6", 0, QApplication::UnicodeUTF8));
        action_run->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_run->setToolTip(QApplication::translate("MainWindow", "Run ABM", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_inputs->setText(QApplication::translate("MainWindow", "Inputs", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_inputs->setToolTip(QApplication::translate("MainWindow", "Switch to the Inputs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_outputs->setText(QApplication::translate("MainWindow", "Outputs", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_outputs->setToolTip(QApplication::translate("MainWindow", "Switch to the Outputs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_VTK->setText(QApplication::translate("MainWindow", "Animation", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_VTK->setToolTip(QApplication::translate("MainWindow", "Switch to animation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_pause->setText(QApplication::translate("MainWindow", "Pause", 0, QApplication::UnicodeUTF8));
        action_load_results->setText(QApplication::translate("MainWindow", "&Load results", 0, QApplication::UnicodeUTF8));
        action_load_results->setIconText(QApplication::translate("MainWindow", "Load results", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        action_load_results->setToolTip(QApplication::translate("MainWindow", "Load result file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        action_add_graph->setText(QApplication::translate("MainWindow", "Add graph", 0, QApplication::UnicodeUTF8));
        action_remove_graph->setText(QApplication::translate("MainWindow", "Remove graph", 0, QApplication::UnicodeUTF8));
        action_remove_all->setText(QApplication::translate("MainWindow", "Remove all", 0, QApplication::UnicodeUTF8));
        action_play_VTK->setText(QApplication::translate("MainWindow", "Play cell animation", 0, QApplication::UnicodeUTF8));
        action_set_speed->setText(QApplication::translate("MainWindow", "Set speed", 0, QApplication::UnicodeUTF8));
        action_save_snapshot->setText(QApplication::translate("MainWindow", "Save snapshot", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_BC_AVIDITY_MEDIAN->setToolTip(QApplication::translate("MainWindow", "TCR avidity has a lognormal distribution, described by the mean and shape parameters.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        label_BC_AVIDITY_MEDIAN->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        label_BC_AVIDITY_MEDIAN->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        label_BC_AVIDITY_MEDIAN->setText(QApplication::translate("MainWindow", "label_BC_AVIDITY_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_BC_AVIDITY_SHAPE->setText(QApplication::translate("MainWindow", "label_BC_AVIDITY_SHAPE", 0, QApplication::UnicodeUTF8));
        label_BC_COGNATE_FRACTION->setText(QApplication::translate("MainWindow", "label_BC_COGNATE_FRACTION", 0, QApplication::UnicodeUTF8));
        label_BC_STIM_RATE_CONSTANT->setText(QApplication::translate("MainWindow", "label_BC_STIM_RATE_CONSTANT", 0, QApplication::UnicodeUTF8));
        label_BC_STIM_HALFLIFE->setText(QApplication::translate("MainWindow", "label_BC_STIM_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_TC_STIM_HALFLIFE->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_MOTILITY_BETA->setText(QApplication::translate("MainWindow", "label_MOTILITY_BETA", 0, QApplication::UnicodeUTF8));
        label_MOTILITY_RHO->setText(QApplication::translate("MainWindow", "label_MOTILITY_RHO", 0, QApplication::UnicodeUTF8));
        label_DIVIDE1_MEDIAN->setText(QApplication::translate("MainWindow", "label_DIVIDE1_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DIVIDE1_SHAPE->setText(QApplication::translate("MainWindow", "label_DIVIDE1_SHAPE", 0, QApplication::UnicodeUTF8));
        alabel_dist->setText(QApplication::translate("MainWindow", "Probability distributions", 0, QApplication::UnicodeUTF8));
        label_DIVIDE2_MEDIAN->setText(QApplication::translate("MainWindow", "label_DIVIDE1_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DIVIDE2_SHAPE->setText(QApplication::translate("MainWindow", "label_DIVIDE1_SHAPE", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_B), QApplication::translate("MainWindow", "B Cell", 0, QApplication::UnicodeUTF8));
        label_DC_BIND_DELAY->setText(QApplication::translate("MainWindow", "label_DC_BIND_DELAY", 0, QApplication::UnicodeUTF8));
        units_DC_BIND_DELAY->setText(QApplication::translate("MainWindow", "min", 0, QApplication::UnicodeUTF8));
        label_DC_DENS_HALFLIFE->setText(QApplication::translate("MainWindow", "label_DC_DENS_HALFLIFE", 0, QApplication::UnicodeUTF8));
        label_MAX_TC_BIND->setText(QApplication::translate("MainWindow", "label_MAX_TC_BIND", 0, QApplication::UnicodeUTF8));
        label_MAX_COG_BIND->setText(QApplication::translate("MainWindow", "label_MAX_COG_BIND", 0, QApplication::UnicodeUTF8));
        units_DC_DENS_HALFLIFE->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_DC_ANTIGEN_SHAPE->setText(QApplication::translate("MainWindow", "label_DC_ANTIGEN_SHAPE", 0, QApplication::UnicodeUTF8));
        label_DC_ANTIGEN_MEDIAN->setText(QApplication::translate("MainWindow", "label_DC_ANTIGEN_MEDIAN", 0, QApplication::UnicodeUTF8));
        label_DC_LIFETIME_SHAPE->setText(QApplication::translate("MainWindow", "label_DC_LIFETIME_SHAPE", 0, QApplication::UnicodeUTF8));
        label_DC_LIFETIME_MEDIAN->setText(QApplication::translate("MainWindow", "label_DC_LIFETIME_MEDIAN", 0, QApplication::UnicodeUTF8));
        alabel_dist_2->setText(QApplication::translate("MainWindow", "Probability distributions", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_DC), QApplication::translate("MainWindow", "DC", 0, QApplication::UnicodeUTF8));
        cbox_USE_S1PR1->setText(QApplication::translate("MainWindow", "Use S1PR1?", 0, QApplication::UnicodeUTF8));
        label_CCL21_BDRY_CONC->setText(QApplication::translate("MainWindow", "label_CCL21_BDRY_CONC", 0, QApplication::UnicodeUTF8));
        label_CCL21_DIFF_COEFF->setText(QApplication::translate("MainWindow", "label_CCL21_DIFF_COEFF", 0, QApplication::UnicodeUTF8));
        label_CCL21_HALFLIFE->setText(QApplication::translate("MainWindow", "label_CCL21_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_CCL21conc->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg.L</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_CCL21diff->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">m</span><span style=\" font-size:8pt; vertical-align:super;\">2</span><span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_CCL21life->setText(QApplication::translate("MainWindow", "h", 0, QApplication::UnicodeUTF8));
        label_CCL21_STRENGTH->setText(QApplication::translate("MainWindow", "label_CCL21_STRENGTH", 0, QApplication::UnicodeUTF8));
        label_CCL21_BDRY_RATE->setText(QApplication::translate("MainWindow", "label_CCL21_BDRY_RATE", 0, QApplication::UnicodeUTF8));
        units_CCL21rate->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg</span>.<span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        cbox_USE_CCR7->setText(QApplication::translate("MainWindow", "Use CCR7?", 0, QApplication::UnicodeUTF8));
        cbox_USE_EBI2->setText(QApplication::translate("MainWindow", "Use EBI2?", 0, QApplication::UnicodeUTF8));
        cbox_USE_CXCR5->setText(QApplication::translate("MainWindow", "Use CXCR5?", 0, QApplication::UnicodeUTF8));
        cbox_USE_S1PR2->setText(QApplication::translate("MainWindow", "Use S1PR2?", 0, QApplication::UnicodeUTF8));
        rbut_CCL21_BDRY_0->setText(QApplication::translate("MainWindow", "Use secretion", 0, QApplication::UnicodeUTF8));
        label_CXCL13_BDRY_CONC->setText(QApplication::translate("MainWindow", "label_CXCL13_BDRY_CONC", 0, QApplication::UnicodeUTF8));
        label_CXCL13_DIFF_COEFF->setText(QApplication::translate("MainWindow", "label_CXCL13_DIFF_COEFF", 0, QApplication::UnicodeUTF8));
        label_CXCL13_HALFLIFE->setText(QApplication::translate("MainWindow", "label_CXCL13_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_CXCL13conc->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg.L</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_CXCL13diff->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">m</span><span style=\" font-size:8pt; vertical-align:super;\">2</span><span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_CXCL13life->setText(QApplication::translate("MainWindow", "h", 0, QApplication::UnicodeUTF8));
        label_CXCL13_STRENGTH->setText(QApplication::translate("MainWindow", "label_CXCL13_STRENGTH", 0, QApplication::UnicodeUTF8));
        label_CXCL13_BDRY_RATE->setText(QApplication::translate("MainWindow", "label_CXCL13_BDRY_RATE", 0, QApplication::UnicodeUTF8));
        units_CXCL13rate->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg</span>.<span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        rbut_CXCL13_BDRY_0->setText(QApplication::translate("MainWindow", "Use secretion", 0, QApplication::UnicodeUTF8));
        label_OXY_BDRY_CONC->setText(QApplication::translate("MainWindow", "label_OXY_BDRY_CONC", 0, QApplication::UnicodeUTF8));
        label_OXY_DIFF_COEFF->setText(QApplication::translate("MainWindow", "label_OXY_DIFF_COEFF", 0, QApplication::UnicodeUTF8));
        label_OXY_HALFLIFE->setText(QApplication::translate("MainWindow", "labelOXYP_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_OXYconc->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg.L</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_OXYdiff->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">m</span><span style=\" font-size:8pt; vertical-align:super;\">2</span><span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_OXYlife->setText(QApplication::translate("MainWindow", "h", 0, QApplication::UnicodeUTF8));
        label_OXY_STRENGTH->setText(QApplication::translate("MainWindow", "label_OXY_STRENGTH", 0, QApplication::UnicodeUTF8));
        label_OXY_BDRY_RATE->setText(QApplication::translate("MainWindow", "label_OXY_BDRY_RATE", 0, QApplication::UnicodeUTF8));
        units_OXYrate->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg</span>.<span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        rbut_OXY_BDRY_0->setText(QApplication::translate("MainWindow", "Use secretion", 0, QApplication::UnicodeUTF8));
        label_S1P_BDRY_CONC->setText(QApplication::translate("MainWindow", "label_S1P_BDRY_CONC", 0, QApplication::UnicodeUTF8));
        label_S1P_DIFF_COEFF->setText(QApplication::translate("MainWindow", "label_S1P_DIFF_COEFF", 0, QApplication::UnicodeUTF8));
        label_S1P_HALFLIFE->setText(QApplication::translate("MainWindow", "label_S1P_HALFLIFE", 0, QApplication::UnicodeUTF8));
        units_S1Pconc->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg.L</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_S1Pdiff->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">m</span><span style=\" font-size:8pt; vertical-align:super;\">2</span><span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        units_S1Plife->setText(QApplication::translate("MainWindow", "h", 0, QApplication::UnicodeUTF8));
        label_S1P_STRENGTH_POS->setText(QApplication::translate("MainWindow", "label_S1P_STRENGTH_POS", 0, QApplication::UnicodeUTF8));
        label_S1P_STRENGTH_NEG->setText(QApplication::translate("MainWindow", "label_S1P_STRENGTH_NEG", 0, QApplication::UnicodeUTF8));
        label_S1P_BDRY_RATE->setText(QApplication::translate("MainWindow", "label_S1P_BDRY_RATE", 0, QApplication::UnicodeUTF8));
        units_S1Prate->setText(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">mg</span>.<span style=\" font-size:8pt;\">s</span><span style=\" font-size:8pt; vertical-align:super;\">-1</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        rbut_S1P_BDRY_0->setText(QApplication::translate("MainWindow", "Use secretion", 0, QApplication::UnicodeUTF8));
        label_chemokine->setText(QApplication::translate("MainWindow", "Chemokine Parameters", 0, QApplication::UnicodeUTF8));
        rbut_S1P_BDRY_1->setText(QApplication::translate("MainWindow", "Use concentration", 0, QApplication::UnicodeUTF8));
        rbut_CCL21_BDRY_1->setText(QApplication::translate("MainWindow", "Use concentration", 0, QApplication::UnicodeUTF8));
        rbut_OXY_BDRY_1->setText(QApplication::translate("MainWindow", "Use concentration", 0, QApplication::UnicodeUTF8));
        rbut_CXCL13_BDRY_1->setText(QApplication::translate("MainWindow", "Use concentration", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_chemo), QApplication::translate("MainWindow", "Chemokine", 0, QApplication::UnicodeUTF8));
        label_IL2_THRESHOLD->setText(QApplication::translate("MainWindow", "label_IL2_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_ACTIVATION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_ACTIVATION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_FIRST_DIVISION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_FIRST_DIVISION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_DIVISION_THRESHOLD->setText(QApplication::translate("MainWindow", "label_DIVISION_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_EXIT_THRESHOLD->setText(QApplication::translate("MainWindow", "label_EXIT_THRESHOLD", 0, QApplication::UnicodeUTF8));
        label_STIMULATION_LIMIT->setText(QApplication::translate("MainWindow", "label_STIMULATION_LIMIT", 0, QApplication::UnicodeUTF8));
        label_40->setText(QApplication::translate("MainWindow", "TCR activation thresholds", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_TCR), QApplication::translate("MainWindow", "TCR Activation", 0, QApplication::UnicodeUTF8));
        label_NX->setText(QApplication::translate("MainWindow", "label_NX", 0, QApplication::UnicodeUTF8));
        label_BLOB_RADIUS->setText(QApplication::translate("MainWindow", "label_BLOB_RADIUS", 0, QApplication::UnicodeUTF8));
        label_BC_FRACTION->setText(QApplication::translate("MainWindow", "label_BC_FRACTION", 0, QApplication::UnicodeUTF8));
        label_RESIDENCE_TIME->setText(QApplication::translate("MainWindow", "label_RESIDENCE_TIME", 0, QApplication::UnicodeUTF8));
        units_RESIDENCE_TIME->setText(QApplication::translate("MainWindow", "hours", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_DAYS1->setText(QApplication::translate("MainWindow", "label_INFLAMM_DAYS1", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_DAYS2->setText(QApplication::translate("MainWindow", "label_INFLAMM_DAYS2", 0, QApplication::UnicodeUTF8));
        label_INFLAMM_LEVEL->setText(QApplication::translate("MainWindow", "label_INFLAMM_LEVEL", 0, QApplication::UnicodeUTF8));
        label_CHEMO_RADIUS->setText(QApplication::translate("MainWindow", "label_CHEMO_RADIUS", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_BASE_EXIT_PROB->setToolTip(QApplication::translate("MainWindow", "Base B cell exit probability", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_BASE_EXIT_PROB->setText(QApplication::translate("MainWindow", "label_BASE_EXIT_PROB", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NDAYS->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_NDAYS->setText(QApplication::translate("MainWindow", "label_NDAYS", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_SEED1->setToolTip(QApplication::translate("MainWindow", "seed vector for RNGs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_SEED1->setText(QApplication::translate("MainWindow", "label_SEED1", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_SEED2->setToolTip(QApplication::translate("MainWindow", "seed vector for RNGs", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_SEED2->setText(QApplication::translate("MainWindow", "label_SEED2", 0, QApplication::UnicodeUTF8));
        label_FLUID_FRACTION->setText(QApplication::translate("MainWindow", "label_FLUID_FRACTION", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NT_ANIMATION->setToolTip(QApplication::translate("MainWindow", "animation update interval", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_NT_ANIMATION->setText(QApplication::translate("MainWindow", "label_NT_ANIMATION", 0, QApplication::UnicodeUTF8));
        units_CHEMO_RADIUS->setText(QApplication::translate("MainWindow", "\302\265m", 0, QApplication::UnicodeUTF8));
        units_NDAYS->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_INFLAMM1->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_INFLAMM2->setText(QApplication::translate("MainWindow", "days", 0, QApplication::UnicodeUTF8));
        units_BLOB_RADIUS->setText(QApplication::translate("MainWindow", "sites", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_NCPU->setToolTip(QApplication::translate("MainWindow", "animation update interval", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_NCPU->setText(QApplication::translate("MainWindow", "label_NCPU", 0, QApplication::UnicodeUTF8));
        cbox_savepos->setText(QApplication::translate("MainWindow", "Save cell paths", 0, QApplication::UnicodeUTF8));
        cbox_IV_SHOW_NONCOGNATE->setText(QApplication::translate("MainWindow", "Display non-cognate B cells ", 0, QApplication::UnicodeUTF8));
        cbox_USE_TRAFFIC->setText(QApplication::translate("MainWindow", "B cell trafficking? ", 0, QApplication::UnicodeUTF8));
        rbut_SPECIES_1->setText(QApplication::translate("MainWindow", "Human", 0, QApplication::UnicodeUTF8));
        rbut_SPECIES_0->setText(QApplication::translate("MainWindow", "Mouse", 0, QApplication::UnicodeUTF8));
        cbox_USE_EXIT_CHEMOTAXIS->setText(QApplication::translate("MainWindow", "Use exit chemotaxis?", 0, QApplication::UnicodeUTF8));
        cbox_COMPUTED_OUTFLOW->setText(QApplication::translate("MainWindow", "Compute B cell outflow?", 0, QApplication::UnicodeUTF8));
        label_INPUT_FILE->setText(QApplication::translate("MainWindow", "Auxiliary input data file", 0, QApplication::UnicodeUTF8));
        tabs->setTabText(tabs->indexOf(tab_run), QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        label_input->setText(QApplication::translate("MainWindow", "Inputs", 0, QApplication::UnicodeUTF8));
        text_more->setDocumentTitle(QString());
        label_hour->setText(QApplication::translate("MainWindow", "Hour:", 0, QApplication::UnicodeUTF8));
        hour_display->setText(QString());
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuEdit->setTitle(QApplication::translate("MainWindow", "Edit", 0, QApplication::UnicodeUTF8));
        menuABM->setTitle(QApplication::translate("MainWindow", "ABM", 0, QApplication::UnicodeUTF8));
        menuGraphs->setTitle(QApplication::translate("MainWindow", "Graphs", 0, QApplication::UnicodeUTF8));
        menuPlayer->setTitle(QApplication::translate("MainWindow", "Player", 0, QApplication::UnicodeUTF8));
        menuSnapshot->setTitle(QApplication::translate("MainWindow", "Snapshot", 0, QApplication::UnicodeUTF8));
        toolBar1->setWindowTitle(QApplication::translate("MainWindow", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BCELL_GUI_H
