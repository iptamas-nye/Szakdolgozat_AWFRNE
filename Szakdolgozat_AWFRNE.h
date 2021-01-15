#pragma once

#include "qcustomplot.h"
#include <QtWidgets/QMainWindow>
#include "ui_Szakdolgozat_AWFRNE.h"

class Szakdolgozat_AWFRNE : public QMainWindow
{
    Q_OBJECT

public:
    Szakdolgozat_AWFRNE(QWidget *parent = Q_NULLPTR);
    ~Szakdolgozat_AWFRNE();
private:
    Ui::Szakdolgozat_AWFRNEClass ui;
    QCustomPlot* customPlot;
    QCPItemText* textItem;
    double minAddedX;                    // min x value of added points
    double maxAddedX;                    // max x value of added points
    double minAddedY;                    // min y value of added points
    double maxAddedY;                    // max y value of added points
    QVector<qreal> addedPointsX;
    QVector<qreal> addedPointsY;
    QVector<qreal> calculatedPointsX;
    QVector<qreal> calculatedPointsY;
    double f_interpolate(double);
    void run();
    void labelMessage(QString, int);
    
private slots:
    void slot_initInterpolation();
    void slot_addPoint();
    void slot_interpolate();
    void slot_onMouseMove(QMouseEvent* event);
    void slot_showCoordiantes();
    void slot_pointSelected();
};
