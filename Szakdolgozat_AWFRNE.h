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
    QVector<qreal> addedPointsX;
    QVector<qreal> addedPointsY;
    QVector<qreal> calculatedPointsX;
    QVector<qreal> calculatedPointsY;
    double f_interpolate(double);
    void labelMessage(QString, int);
    double minAddedPointsX();
    double minAddedPointsY();
    double maxAddedPointsX();
    double maxAddedPointsY();
    void run();
    
private slots:
    void slot_initInterpolation();
    void slot_addPoint();
    void slot_interpolate();
    void slot_showCoordiantes();
    void slot_pointSelected();
    void slot_onMouseMove(QMouseEvent* event);
    void slot_onWheel(QWheelEvent* event);

};
