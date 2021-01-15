#include "Szakdolgozat_AWFRNE.h"
#include "stdafx.h"

#define ERROR 0
#define SUCCESS 1
#define INTERPOLATION 0
#define ADDEDPOINTS 1
#define INTERPOLATION_STEPS 200


Szakdolgozat_AWFRNE::Szakdolgozat_AWFRNE(QWidget *parent) : 
    calculatedPointsX(INTERPOLATION_STEPS),
    calculatedPointsY(INTERPOLATION_STEPS),
    minAddedX(std::numeric_limits<double>::max()),
    maxAddedX(std::numeric_limits<double>::min()),
    minAddedY(std::numeric_limits<double>::max()),
    maxAddedY(std::numeric_limits<double>::min()),
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
    ui.customPlot->xAxis->setLabel("x");
    ui.customPlot->yAxis->setLabel("y");
    ui.customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    connect(ui.actionNew, &QAction::triggered, this, &Szakdolgozat_AWFRNE::slot_initInterpolation);
    connect(ui.customPlot, &QCustomPlot::mouseMove, this, &Szakdolgozat_AWFRNE::slot_onMouseMove);
    connect(ui.pushButtonAddPoint, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_addPoint);
    connect(ui.pushButtonInterpolation, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_interpolate);
    connect(ui.checkBoxShowCoords, &QCheckBox::stateChanged, this, &Szakdolgozat_AWFRNE::slot_showCoordiantes);
    connect(ui.customPlot, &QCustomPlot::selectionChangedByUser, this, &Szakdolgozat_AWFRNE::slot_pointSelected);
}

void Szakdolgozat_AWFRNE::labelMessage(QString message, int error)
{
    ui.labelMessage->setText(message);
    ui.labelMessage->setStyleSheet(error ? "color: green;" : "color: red;");
}

void Szakdolgozat_AWFRNE::slot_initInterpolation()
{
    minAddedX = std::numeric_limits<double>::max();
    maxAddedX = std::numeric_limits<double>::min();
    minAddedY = std::numeric_limits<double>::max();
    maxAddedY = std::numeric_limits<double>::min();
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
    if (addedPointsX.back() < minAddedX) minAddedX = addedPointsX.back();
    if (addedPointsX.back() > maxAddedX) maxAddedX = addedPointsX.back();
    if (addedPointsY.back() < minAddedY) minAddedY = addedPointsY.back();
    if (addedPointsY.back() > maxAddedY) maxAddedY = addedPointsY.back();
    ui.customPlot->graph(ADDEDPOINTS)->setVisible(true);
    ui.customPlot->graph(ADDEDPOINTS)->setData(addedPointsX, addedPointsY);
    if (addedPointsX.size() == 1)
    {
        ui.customPlot->xAxis->setRange(minAddedX - 1, minAddedX + 1);
        ui.customPlot->yAxis->setRange(minAddedY - 1, minAddedY + 1);
    }
    else
    {
        ui.customPlot->xAxis->setRange(minAddedX-(maxAddedX-minAddedX), maxAddedX + (maxAddedX - minAddedX));
        ui.customPlot->yAxis->setRange(minAddedY - (maxAddedY - minAddedY), maxAddedY + (maxAddedY - minAddedY));
    }  
    ui.customPlot->replot();
    labelMessage("Point added successfully!", SUCCESS);
}

double Szakdolgozat_AWFRNE::f_interpolate(double x)
{ 
    auto nofPoints = addedPointsX.size();
    if (nofPoints == 1) return addedPointsY.at(0);
    double result = 0; // Initialize result 
    double minToInterpolate = minAddedX;
    double maxToInterpolate = maxAddedX;
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

void Szakdolgozat_AWFRNE::slot_interpolate()
{
    double minToInterpolate = (addedPointsX.size() == 1) ? ui.customPlot->xAxis->range().lower : minAddedX;
    double maxToInterpolate = (addedPointsX.size() == 1) ? ui.customPlot->xAxis->range().upper : maxAddedX;
    double step = (maxToInterpolate - minToInterpolate) / INTERPOLATION_STEPS;

    for (auto ip = 0; ip < INTERPOLATION_STEPS; ip++)  //elements of calculatedPointsX and Y
    {
        calculatedPointsX[ip] = minToInterpolate + ip * step;
        calculatedPointsY[ip] = f_interpolate(calculatedPointsX[ip]);
    }
    ui.customPlot->graph(INTERPOLATION)->setVisible(true);
    ui.customPlot->graph(INTERPOLATION)->setData(calculatedPointsX, calculatedPointsY);
    ui.customPlot->graph(INTERPOLATION)->rescaleAxes();
    ui.customPlot->yAxis->setRange(minAddedY, maxAddedY);
    ui.customPlot->replot();      
    labelMessage("Successful interpolation!", SUCCESS);
}

void Szakdolgozat_AWFRNE::slot_onMouseMove(QMouseEvent* event)
{
    double x = ui.customPlot->xAxis->pixelToCoord(event->pos().x());
    double y = ui.customPlot->yAxis->pixelToCoord(event->pos().y());
    statusBar()->showMessage(QString("%1, %2").arg(x).arg(y));
    if (textItem->visible())
    {
        textItem->setText(QString("%1, %2").arg(x).arg(y));
        textItem->position->setCoords(QPointF(x, y));
        textItem->setFont(QFont(font().family(), 10));
        ui.customPlot->replot();
    }
    
}

void Szakdolgozat_AWFRNE::slot_showCoordiantes()
{
    textItem->setVisible(ui.checkBoxShowCoords->isChecked());
    ui.customPlot->replot();
}

void Szakdolgozat_AWFRNE::slot_pointSelected()
{
    auto graph = ui.customPlot->graph(ADDEDPOINTS);
    auto selection = graph->selection();
    foreach(QCPDataRange dataRange, selection.dataRanges())
    {
        //QCPGraphDataContainer::const_iterator begin = graph->data()->at(dataRange.begin()); // get range begin iterator from index
        //QCPGraphDataContainer::const_iterator end = graph->data()->at(dataRange.end()); // get range end iterator from index
        auto begin = graph->data()->at(dataRange.begin()); // get range begin iterator from index
        auto end = graph->data()->at(dataRange.end()); // get range end iterator from index
        for (QCPGraphDataContainer::const_iterator it = begin; it != end; ++it)
        {
            // iterator "it" will go through all selected data points, as an example, we calculate the value average
            qDebug() << it->key;
            qDebug() << it->value;
        }
    }
}
