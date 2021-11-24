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
    QVector<qreal> lagrangePointsX;
    QVector<qreal> lagrangePointsY;
    QVector<qreal> hermitePointsX;
    QVector<qreal> hermitePointsY;
    QVector<qreal> derivativePointsX;
    QVector<qreal> derivativePointsY;
    QVector<qreal> hermiteDerivatives;
    double lagrange_interpolate(double);
    double hermite_interpolate(double);
    double calculateDerivative(double);
    void drawDerivative(double, double, double);
    double omegaDerivative(double);
    double omegaSecondDerivative(double);
    void labelMessage(QString, int);
    void checkZoomLimit();
    double minAddedPointsX();
    double minAddedPointsY();
    double maxAddedPointsX();
    double maxAddedPointsY();
    void run();
    
private slots:
    void slot_initInterpolation();
    void slot_addPoint();
    void slot_interpolate();
    void slot_rangeChanged(QCPRange);
    void slot_pointSelected();
    void slot_onMouseMove(QMouseEvent* event);
    void slot_onWheel(QWheelEvent* event);
    void slot_buttonClicked(QMouseEvent* event);

};
