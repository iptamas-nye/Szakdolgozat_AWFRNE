#include "Szakdolgozat_AWFRNE.h"
#include "stdafx.h"

#define ERROR 0
#define SUCCESS 1
#define INTERPOLATION 0
#define ADDEDPOINTS 1
#define DERIVATIVES 2
#define INTERPOLATION_STEPS 200


Szakdolgozat_AWFRNE::Szakdolgozat_AWFRNE(QWidget *parent) : 
    calculatedPointsX(INTERPOLATION_STEPS),
    calculatedPointsY(INTERPOLATION_STEPS),
    QMainWindow(parent)
{
    ui.setupUi(this);
    textItem = new QCPItemText(ui.customPlot);
    run(); 
}

Szakdolgozat_AWFRNE::~Szakdolgozat_AWFRNE()
{
    //delete customPlot;
}

void Szakdolgozat_AWFRNE::run()
{ 
    textItem->setVisible(false);
    ui.lineEditX->setValidator(new QRegExpValidator(QRegExp("[+-]?\\d*\\.?\\d+"), this));
    ui.lineEditY->setValidator(new QRegExpValidator(QRegExp("[+-]?\\d*\\.?\\d+"), this));

    ui.customPlot->addGraph();
    ui.customPlot->graph(INTERPOLATION)->setSelectable(QCP::SelectionType::stNone);
    ui.customPlot->graph(INTERPOLATION)->setLineStyle(QCPGraph::lsLine);
    ui.customPlot->graph(INTERPOLATION)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 0));
    ui.customPlot->graph(INTERPOLATION)->setPen(QPen(Qt::red));
    ui.customPlot->addGraph();
    ui.customPlot->graph(ADDEDPOINTS)->setSelectable(QCP::SelectionType::stSingleData);
    ui.customPlot->graph(ADDEDPOINTS)->setLineStyle(QCPGraph::lsNone);
    ui.customPlot->graph(ADDEDPOINTS)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    ui.customPlot->graph(ADDEDPOINTS)->setPen(QPen(Qt::blue));
    ui.customPlot->addGraph();
    ui.customPlot->graph(DERIVATIVES)->setSelectable(QCP::SelectionType::stNone);
    ui.customPlot->graph(DERIVATIVES)->setLineStyle(QCPGraph::lsLine);
    ui.customPlot->graph(DERIVATIVES)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 0));
    ui.customPlot->graph(DERIVATIVES)->setPen(QPen(Qt::green));
    ui.customPlot->xAxis->setLabel("x");
    ui.customPlot->yAxis->setLabel("y");
    ui.customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    derivativePointsX.append(0); derivativePointsX.append(0);
    derivativePointsY.append(0); derivativePointsY.append(0);

    connect(ui.actionNew, &QAction::triggered, this, &Szakdolgozat_AWFRNE::slot_initInterpolation);  
    connect(ui.pushButtonAddPoint, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_addPoint);
    connect(ui.pushButtonInterpolation, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_interpolate);
    connect(ui.customPlot, &QCustomPlot::selectionChangedByUser, this, &Szakdolgozat_AWFRNE::slot_pointSelected);
    connect(ui.customPlot, &QCustomPlot::mouseMove, this, &Szakdolgozat_AWFRNE::slot_onMouseMove);
    connect(ui.customPlot, &QCustomPlot::mouseWheel, this, &Szakdolgozat_AWFRNE::slot_onWheel);
    connect(ui.customPlot, &QCustomPlot::mousePress, this, &Szakdolgozat_AWFRNE::slot_buttonClicked);
}

void Szakdolgozat_AWFRNE::slot_initInterpolation()
{
    addedPointsX.clear();
    addedPointsY.clear();
    ui.customPlot->graph(ADDEDPOINTS)->setVisible(false);
    ui.customPlot->graph(INTERPOLATION)->setVisible(false); 
    ui.customPlot->xAxis->setRange(-5, 5);
    ui.customPlot->yAxis->setRange(-5, 5);
    ui.labelMessage->setText("");
    ui.lineEditX->setText("");
    ui.lineEditY->setText("");
    ui.customPlot->replot();
}

void Szakdolgozat_AWFRNE::slot_addPoint()
{
    if (addedPointsX.size() >= 10)
    {
        labelMessage("No more points allowed!", ERROR);
        return;
    }
    auto lineEditXDouble = ui.lineEditX->text().toDouble();
    auto lineEditYDouble = ui.lineEditY->text().toDouble();

    if(addedPointsX.contains(lineEditXDouble))
    {
        labelMessage("Error! X coordinate already added!", ERROR);
        return;
    }
    addedPointsX.append(lineEditXDouble);
    addedPointsY.append(lineEditYDouble);
    ui.customPlot->graph(ADDEDPOINTS)->setVisible(true);
    ui.customPlot->graph(ADDEDPOINTS)->setData(addedPointsX, addedPointsY);
    auto minAddedX = minAddedPointsX();
    auto maxAddedX = maxAddedPointsX();
    auto minAddedY = minAddedPointsY();
    auto maxAddedY = maxAddedPointsY();
    if (addedPointsX.size() == 1)
    {
        ui.customPlot->xAxis->setRange(minAddedX-2, minAddedX+2);
        ui.customPlot->yAxis->setRange(minAddedY-2, maxAddedY+2);
    }
    else
    {
        ui.customPlot->xAxis->setRange(minAddedX-(maxAddedX-minAddedX), maxAddedX + (maxAddedX - minAddedX));
        ui.customPlot->yAxis->setRange(minAddedY - (maxAddedY - minAddedY), maxAddedY + (maxAddedY - minAddedY));
    }  
    ui.customPlot->replot();
    labelMessage("Point added successfully!", SUCCESS);
}

void Szakdolgozat_AWFRNE::drawDerivative(double x, double y, double der)
{
    double half = (maxAddedPointsX() - minAddedPointsX()) / 4;
    double xMin = x - half, xMax = x + half;
    double yMin = y - half*der, yMax = y + half*der;
    derivativePointsX[0] = xMin; derivativePointsX[1] = xMax;
    derivativePointsY[0] = yMin; derivativePointsY[1] = yMax;
    ui.customPlot->graph(DERIVATIVES)->setData(derivativePointsX, derivativePointsY);
}

double Szakdolgozat_AWFRNE::f_interpolate(double x)
{ 
    auto nofPoints = addedPointsX.size();
    if (nofPoints == 1) return addedPointsY.at(0);
    double result = 0; // Initialize result 
    double minToInterpolate = minAddedPointsX();
    double maxToInterpolate = maxAddedPointsX();
    for (auto i = 0; i < nofPoints; i++)
    {
        double term = addedPointsY.at(i);
        for (auto j = 0; j < nofPoints; j++)
        {
            if (j != i)
                term = term * (x - addedPointsX.at(j)) / double(addedPointsX.at(i) - addedPointsX.at(j));
        }
        result += term;
    }   
    return result;
}

double Szakdolgozat_AWFRNE::calculateDerivative(double x0)
{
    const double epszilon = 1.0e-6; 
    double x1 = x0 - epszilon;
    double x2 = x0 + epszilon;
    double y1 = f_interpolate(x1);
    double y2 = f_interpolate(x2);
    return (y2 - y1) / (x2 - x1);
}

void Szakdolgozat_AWFRNE::slot_interpolate()
{
    double minToInterpolate = (addedPointsX.size() == 1) ? ui.customPlot->xAxis->range().lower : minAddedPointsX();
    double maxToInterpolate = (addedPointsX.size() == 1) ? ui.customPlot->xAxis->range().upper : maxAddedPointsX();
    double step = (maxToInterpolate - minToInterpolate) / (INTERPOLATION_STEPS);
    auto maxCalculatedPointsY = std::numeric_limits<double>::min();
    auto minCalculatedPointsY = std::numeric_limits<double>::max();
    for (auto i = 0; i < INTERPOLATION_STEPS-1; i++)  //elements of calculatedPointsX and Y
    {
        calculatedPointsX[i] = minToInterpolate + i * step;
        calculatedPointsY[i] = f_interpolate(calculatedPointsX.at(i));
        if (calculatedPointsY.at(i) > maxCalculatedPointsY) maxCalculatedPointsY = calculatedPointsY.at(i);
        if (calculatedPointsY.at(i) < minCalculatedPointsY) minCalculatedPointsY = calculatedPointsY.at(i);
    }
    calculatedPointsX.back() = minToInterpolate + INTERPOLATION_STEPS * step;
    calculatedPointsY.back() = f_interpolate(calculatedPointsX.back());
    if (calculatedPointsY.back() > maxCalculatedPointsY) maxCalculatedPointsY = calculatedPointsY.back();
    if (calculatedPointsY.back() < minCalculatedPointsY) minCalculatedPointsY = calculatedPointsY.back();

    ui.customPlot->graph(INTERPOLATION)->setVisible(true);
    ui.customPlot->graph(INTERPOLATION)->setData(calculatedPointsX, calculatedPointsY);
    ui.customPlot->graph(INTERPOLATION)->rescaleAxes();
    if (addedPointsX.size() == 1)
    {
        ui.customPlot->yAxis->setRange(minCalculatedPointsY + 2, maxCalculatedPointsY - 2);
    }
    else
    {
        ui.customPlot->yAxis->setRange(minCalculatedPointsY, maxCalculatedPointsY);
    }
    ui.customPlot->replot();      
    labelMessage("Successful interpolation!", SUCCESS);
}

void Szakdolgozat_AWFRNE::slot_pointSelected()
{
    auto graph = ui.customPlot->graph(ADDEDPOINTS);
    auto selection = graph->selection();
    foreach(QCPDataRange dataRange, selection.dataRanges())
    {     
        auto begin = graph->data()->at(dataRange.begin()); // get range begin iterator from index
        auto end = graph->data()->at(dataRange.end()); // get range end iterator from index
        for (QCPGraphDataContainer::const_iterator it = begin; it != end; ++it)
        {
            qDebug() << it->key;
            qDebug() << it->value;
        }
    }
}

void Szakdolgozat_AWFRNE::slot_onMouseMove(QMouseEvent* event)
{
    if (!ui.checkBoxShowDerivative->isChecked())
    {
        return;
    }
    double x = ui.customPlot->xAxis->pixelToCoord(event->pos().x());
    double y = ui.customPlot->yAxis->pixelToCoord(event->pos().y());
    double derivative = calculateDerivative(x);
    statusBar()->showMessage(QString("%1, %2").arg(x).arg(y));
    double buffer = (maxAddedPointsY() - minAddedPointsY()) / 200;
    if ((y > f_interpolate(x) - buffer) && (y < f_interpolate(x) + buffer))
    {
        textItem->setText(QString("%1").arg(derivative));
        textItem->position->setCoords(QPointF(x, y));
        textItem->setFont(QFont(font().family(), 10));
        textItem->setVisible(true);
        ui.customPlot->graph(DERIVATIVES)->setVisible(true);
        drawDerivative(x, y, derivative);
        ui.customPlot->replot();
    }
    else 
    {
        textItem->setVisible(false);
        ui.customPlot->graph(DERIVATIVES)->setVisible(false);
        ui.customPlot->replot();
    }

}

void Szakdolgozat_AWFRNE::slot_buttonClicked(QMouseEvent* event)
{

    if (event->button() == Qt::RightButton)
    {
        auto graph = ui.customPlot->graph(ADDEDPOINTS);
        auto selection = graph->selection();
        double xToDelete, yToDelete;
        bool toDelete = false;
        foreach(QCPDataRange dataRange, selection.dataRanges())
        {
            auto begin = graph->data()->at(dataRange.begin()); // get range begin iterator from index
            auto end = graph->data()->at(dataRange.end()); // get range end iterator from index
            
            for (QCPGraphDataContainer::const_iterator it = begin; it != end; ++it)
            {
                xToDelete = it->key;
                yToDelete = it->value;
                toDelete = true;
            }
        }
        if (toDelete && !addedPointsX.empty())
        {
            for (auto i = 0; i < addedPointsX.size(); i++)
            {
                if (addedPointsX.at(i) == xToDelete && addedPointsY.at(i) == yToDelete)
                {
                    QMessageBox::StandardButton reply;
                    reply = QMessageBox::question(this, "Warning!", "Would you like to delete this point?",
                        QMessageBox::Yes | QMessageBox::No);
                    if (reply == QMessageBox::Yes) {
                        ui.customPlot->graph(INTERPOLATION)->setVisible(false);
                        ui.customPlot->graph(ADDEDPOINTS)->data()->remove(addedPointsX.at(i));
                        addedPointsX.remove(i);
                        addedPointsY.remove(i);
                        ui.customPlot->replot();
                    }              
                }
            }            
        }
    }
}

void Szakdolgozat_AWFRNE::slot_onWheel(QWheelEvent* event)
{
    qDebug() << "wheelEvent";
}


void Szakdolgozat_AWFRNE::labelMessage(QString message, int error)
{
    ui.labelMessage->setText(message);
    ui.labelMessage->setStyleSheet(error ? "color: green;" : "color: red;");
}

double Szakdolgozat_AWFRNE::minAddedPointsX()
{
    if (addedPointsX.isEmpty()) return std::numeric_limits<double>::max();
    auto minX = std::numeric_limits<double>::max();
    for (auto i : addedPointsX)
    {
        if (i < minX) minX = i;
    }
    return minX;
}

double Szakdolgozat_AWFRNE::minAddedPointsY()
{
    if (addedPointsY.isEmpty()) return std::numeric_limits<double>::max();
    auto minY = std::numeric_limits<double>::max();
    for (auto i : addedPointsY)
    {
        if (i < minY) minY = i;
    }
    return minY;
}

double Szakdolgozat_AWFRNE::maxAddedPointsX()
{
    if (addedPointsX.isEmpty()) return std::numeric_limits<double>::min();
    auto maxX = std::numeric_limits<double>::min();
    for (auto i : addedPointsX)
    {
        if (i > maxX) maxX = i;
    }
    return maxX;
}

double Szakdolgozat_AWFRNE::maxAddedPointsY()
{
    if (addedPointsY.isEmpty()) return std::numeric_limits<double>::min();
    auto maxY = std::numeric_limits<double>::min();
    for (auto i : addedPointsY)
    {
        if (i > maxY) maxY = i;
    }
    return maxY;
}


