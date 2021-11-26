#include "Szakdolgozat_AWFRNE.h"
#include "stdafx.h"

#define ERROR 0
#define SUCCESS 1
#define LAGRANGE 0
#define HERMITE 1
#define ADDEDPOINTS 2
#define DERIVATIVES 3
#define INTERPOLATION_STEPS 200
#define MAX_ADDED_POINTS 10
#define EMPTY_DERIVATIVE 654321


Szakdolgozat_AWFRNE::Szakdolgozat_AWFRNE(QWidget *parent) : 
    lagrangePointsX(INTERPOLATION_STEPS),
    lagrangePointsY(INTERPOLATION_STEPS),
    hermitePointsX(INTERPOLATION_STEPS),
    hermitePointsY(INTERPOLATION_STEPS),
    QMainWindow(parent)
{
    ui.setupUi(this);
    run(); 
}

Szakdolgozat_AWFRNE::~Szakdolgozat_AWFRNE()
{
    
}

void Szakdolgozat_AWFRNE::run()
{ 
    ui.lineEditX->setValidator(new QRegExpValidator(QRegExp("[+-]?\\d*\\.?\\d+"), this));
    ui.lineEditY->setValidator(new QRegExpValidator(QRegExp("[+-]?\\d*\\.?\\d+"), this));
    ui.lineEditDerivative->setValidator(new QRegExpValidator(QRegExp("[+-]?\\d*\\.?\\d+"), this));

    ui.customPlot->addGraph();
    ui.customPlot->graph(LAGRANGE)->setSelectable(QCP::SelectionType::stNone);
    ui.customPlot->graph(LAGRANGE)->setLineStyle(QCPGraph::lsLine);
    ui.customPlot->graph(LAGRANGE)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 0));
    ui.customPlot->graph(LAGRANGE)->setPen(QPen(Qt::red));
    ui.customPlot->addGraph();
    ui.customPlot->graph(HERMITE)->setSelectable(QCP::SelectionType::stNone);
    ui.customPlot->graph(HERMITE)->setLineStyle(QCPGraph::lsLine);
    ui.customPlot->graph(HERMITE)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 0));
    ui.customPlot->graph(HERMITE)->setPen(QPen(Qt::cyan));
    ui.customPlot->addGraph();
    ui.customPlot->graph(ADDEDPOINTS)->setSelectable(QCP::SelectionType::stSingleData);
    ui.customPlot->graph(ADDEDPOINTS)->setLineStyle(QCPGraph::lsNone);
    ui.customPlot->graph(ADDEDPOINTS)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    ui.customPlot->graph(ADDEDPOINTS)->setPen(QPen(Qt::blue));
    ui.customPlot->xAxis->setLabel("x");
    ui.customPlot->yAxis->setLabel("y");
    ui.customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    connect(ui.actionNew, &QAction::triggered, this, &Szakdolgozat_AWFRNE::slot_initInterpolation);  
    connect(ui.pushButtonAddPoint, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_addPoint);
    connect(ui.pushButtonInterpolate, &QPushButton::clicked, this, &Szakdolgozat_AWFRNE::slot_interpolate);
    connect(ui.customPlot, &QCustomPlot::mouseMove, this, &Szakdolgozat_AWFRNE::slot_onMouseMove);
    connect(ui.customPlot, &QCustomPlot::mousePress, this, &Szakdolgozat_AWFRNE::slot_buttonClicked);
    connect(ui.customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), this, SLOT(slot_rangeChanged(QCPRange)));
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

void Szakdolgozat_AWFRNE::checkZoomLimit()
{
    auto minX = minAddedPointsX();
    auto maxX = maxAddedPointsX();
    auto diff = maxX - minX;
    auto minLimit = minX - diff * 2;
    auto maxLimit = maxX + diff * 2;
}

void Szakdolgozat_AWFRNE::slot_initInterpolation()
{
    addedPointsX.clear();
    addedPointsY.clear();
    ui.customPlot->graph(ADDEDPOINTS)->setVisible(false);
    ui.customPlot->graph(LAGRANGE)->setVisible(false); 
    ui.customPlot->graph(HERMITE)->setVisible(false);
    ui.customPlot->xAxis->blockSignals(true);  // block signals
    ui.customPlot->xAxis->setRange(-5, 5);
    ui.customPlot->yAxis->setRange(-5, 5);
    ui.customPlot->xAxis->blockSignals(false);  // block signals
    ui.labelMessage->setText("");
    ui.lineEditX->setText("");
    ui.lineEditY->setText("");
    ui.lineEditDerivative->setText("");
    ui.customPlot->replot();
}

void Szakdolgozat_AWFRNE::slot_addPoint()
{
    if (ui.lineEditX->text().isEmpty() || ui.lineEditY->text().isEmpty())
    {
        labelMessage("Please add numbers!", ERROR);
        return;
    }
    if (addedPointsX.size() >= MAX_ADDED_POINTS)
    {
        labelMessage("No more points allowed!", ERROR);
        return;
    }
    auto lineEditXDouble = ui.lineEditX->text().toDouble();
    auto lineEditYDouble = ui.lineEditY->text().toDouble();
    double lineEditDerivative = 0.0;
    if (ui.lineEditDerivative->text().isEmpty())
    {
        lineEditDerivative = double(EMPTY_DERIVATIVE);
    }
    else 
    {
        lineEditDerivative = ui.lineEditDerivative->text().toDouble();
    }
    
    if(addedPointsX.contains(lineEditXDouble))
    {
        labelMessage("Error! X coordinate already added!", ERROR);
        return;
    }
    addedPointsX.append(lineEditXDouble);
    addedPointsY.append(lineEditYDouble);
    hermiteDerivatives.append(lineEditDerivative);
    ui.customPlot->graph(ADDEDPOINTS)->setVisible(true);
    ui.customPlot->graph(ADDEDPOINTS)->setData(addedPointsX, addedPointsY);
    auto minAddedX = minAddedPointsX();
    auto maxAddedX = maxAddedPointsX();
    auto minAddedY = minAddedPointsY();
    auto maxAddedY = maxAddedPointsY();
    ui.customPlot->xAxis->blockSignals(true);  // block signals
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
    ui.customPlot->xAxis->blockSignals(false); // block signals
    ui.customPlot->replot();
    labelMessage("Point added successfully!", SUCCESS);
}

double Szakdolgozat_AWFRNE::omegaDerivative(double x)
{
    auto size = addedPointsX.size();
    if (2 == size)
        return 2 * x - addedPointsX.at(0) - addedPointsX.at(1);
    auto sum = 0.0, prod = 1.0;
    for (auto i = 0; i < size; i++)
    {
        prod = 1;
        for (auto j = 0; j < size; j++)
        {
            if (i != j)
            {
                prod *= (x - addedPointsX.at(j));
            }
        }
        sum += prod;
    }
    return sum;
}

double Szakdolgozat_AWFRNE::omegaSecondDerivative(double x)
{
    auto size = addedPointsX.size();
    if (2 == size)
        return 2;
    auto sum = 0.0, prod = 1.0;
    for (auto i = 0; i < size - 1; i++)
    {    
        for (auto j = i + 1; j < size; j++)
        {
            prod = 1;
            for (auto k = 0; k < size; k++)
            {
                if ((i != k) && (j != k))
                {
                    prod *= (x - addedPointsX.at(k));
                }
            }
            sum += 2 * prod;
        }
        
    }
    return sum;
}

double Szakdolgozat_AWFRNE::lagrange_interpolate(double x)
{ 
    auto nofPoints = addedPointsX.size();
    if (nofPoints == 1) return addedPointsY.at(0);
    double result = 0; // Initialize result 
    for (auto i = 0; i < nofPoints; i++)
    {
        double Li = addedPointsY.at(i);
        for (auto j = 0; j < nofPoints; j++)
        {
            if (j != i)
                Li = Li * (x - addedPointsX.at(j)) / double(addedPointsX.at(i) - addedPointsX.at(j));
        }
        result += Li;
    }   
    return result;
}

double Szakdolgozat_AWFRNE::hermite_interpolate(double x)
{
    auto nofPoints = addedPointsX.size();
    if (nofPoints == 1) return hermiteDerivatives.at(0);
    double result, Li, h0=0, h1 = 0;
    for (auto i = 0; i < nofPoints; i++)
    {
        Li = 1;
        for (auto j = 0; j < nofPoints; j++)
        {
            if (j != i)
                Li = Li * (x - addedPointsX.at(j)) / double(addedPointsX.at(i) - addedPointsX.at(j));
        }
        h0 += addedPointsY.at(i) * Li * Li * (1 - (x - addedPointsX.at(i)) * (omegaSecondDerivative(addedPointsX.at(i))) / omegaDerivative(addedPointsX.at(i)));
        h1 += hermiteDerivatives.at(i) * Li * Li * (x - addedPointsX.at(i));
    }
    
    return h0 + h1;
}

void Szakdolgozat_AWFRNE::slot_interpolate()
{
    double minToInterpolate = ui.customPlot->xAxis->range().lower;
    double maxToInterpolate = ui.customPlot->xAxis->range().upper;
    double step = (maxToInterpolate - minToInterpolate) / (INTERPOLATION_STEPS);
    auto maxInterpolationPointsY = std::numeric_limits<double>::min();
    auto minInterpolationPointsY = std::numeric_limits<double>::max();
    if (ui.checkBoxLagrange->isChecked()) // ---------------- Lagrange -------------------------
    {
        for (auto i = 0; i < INTERPOLATION_STEPS-1; i++)  //elements of lagrangePointsX and Y
        {
            lagrangePointsX[i] = minToInterpolate + i * step;
            lagrangePointsY[i] = lagrange_interpolate(lagrangePointsX.at(i));
            if (lagrangePointsY.at(i) > maxInterpolationPointsY) maxInterpolationPointsY = lagrangePointsY.at(i);
            if (lagrangePointsY.at(i) < minInterpolationPointsY) minInterpolationPointsY = lagrangePointsY.at(i);
        }
        lagrangePointsX.back() = minToInterpolate + INTERPOLATION_STEPS * step;
        lagrangePointsY.back() = lagrange_interpolate(lagrangePointsX.back());
        if (lagrangePointsY.back() > maxInterpolationPointsY) maxInterpolationPointsY = lagrangePointsY.back();
        if (lagrangePointsY.back() < minInterpolationPointsY) minInterpolationPointsY = lagrangePointsY.back();

        ui.customPlot->graph(LAGRANGE)->setVisible(true);
        ui.customPlot->graph(LAGRANGE)->setData(lagrangePointsX, lagrangePointsY);
        ui.customPlot->graph(LAGRANGE)->rescaleAxes();
        labelMessage("Successful interpolation", SUCCESS);
    }
    else
    {
        ui.customPlot->graph(LAGRANGE)->setVisible(false);
    }
    if (ui.checkBoxHermite->isChecked())  // ---------------- Hermit -------------------------
    {
        if (hermiteDerivatives.contains(EMPTY_DERIVATIVE))
        {
            ui.customPlot->graph(HERMITE)->setVisible(false);
            labelMessage("Hermite interpolation needs derivatives at all points!", ERROR);
        }
        else
        {
            for (auto i = 0; i < INTERPOLATION_STEPS - 1; i++)  //elements of hermitePointsX and Y
            {
                hermitePointsX[i] = minToInterpolate + i * step;
                hermitePointsY[i] = hermite_interpolate(hermitePointsX.at(i));
                if (hermitePointsY.at(i) > maxInterpolationPointsY) maxInterpolationPointsY = hermitePointsY.at(i);
                if (hermitePointsY.at(i) < minInterpolationPointsY) minInterpolationPointsY = hermitePointsY.at(i);
            }
            hermitePointsX.back() = minToInterpolate + INTERPOLATION_STEPS * step;
            hermitePointsY.back() = hermite_interpolate(hermitePointsX.back());
            if (hermitePointsY.back() > maxInterpolationPointsY) maxInterpolationPointsY = hermitePointsY.back();
            if (hermitePointsY.back() < minInterpolationPointsY) minInterpolationPointsY = hermitePointsY.back();

            ui.customPlot->graph(HERMITE)->setVisible(true);
            ui.customPlot->graph(HERMITE)->setData(hermitePointsX, hermitePointsY);
            ui.customPlot->graph(HERMITE)->rescaleAxes();
            labelMessage("Successful interpolation", SUCCESS);
        } 
    }
    else 
    {
        ui.customPlot->graph(HERMITE)->setVisible(false);
    }
    ui.customPlot->xAxis->blockSignals(true);  // block signals
    if (addedPointsX.size() == 1)
    {
        ui.customPlot->yAxis->setRange(minInterpolationPointsY + 2, maxInterpolationPointsY - 2);
    }
    else
    {
        ui.customPlot->yAxis->setRange(minInterpolationPointsY, maxInterpolationPointsY);
    }
    ui.customPlot->xAxis->blockSignals(false);  // block signals
    ui.customPlot->replot();
}

void Szakdolgozat_AWFRNE::slot_rangeChanged(QCPRange range)
{
    if (ui.customPlot->graph(LAGRANGE)->visible())
    {
        slot_interpolate();
    }
    if (ui.customPlot->graph(HERMITE)->visible())
    {
        slot_interpolate();
    }
    checkZoomLimit();
}


void Szakdolgozat_AWFRNE::slot_onMouseMove(QMouseEvent* event)
{
    double x = ui.customPlot->xAxis->pixelToCoord(event->pos().x());
    double y = ui.customPlot->yAxis->pixelToCoord(event->pos().y());
    statusBar()->showMessage(QString("%1, %2").arg(x).arg(y));
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
                    QString message = QString::fromUtf8("Szeretné törölni a kijelölt pontot?");
                    QMessageBox::StandardButton reply;
                    reply = QMessageBox::question(this, "Figyelem!", message,
                        QMessageBox::Yes | QMessageBox::No);
                    if (reply == QMessageBox::Yes) {
                        ui.customPlot->graph(LAGRANGE)->setVisible(false);
                        ui.customPlot->graph(HERMITE)->setVisible(false);
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



