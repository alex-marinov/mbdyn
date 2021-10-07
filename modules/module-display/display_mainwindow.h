/* Copyright (C) 2008-2021
 *
 * Matteo Daniele <matteo.daniele@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */
// Author: Matteo Daniele <matteo.daniele@polimi.it>
// simulation data display for mbdyn
#ifndef DISPLAYMAINWINDOW_H
#define DISPLAYMAINWINDOW_H
#include "userelem.h"
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QApplication>
#include <QtCharts/QtCharts>

class MainWindow : public QMainWindow
{
    Q_OBJECT

    public:
        MainWindow(QWidget *parent = nullptr);
        ~MainWindow();

    private:
       
};

#endif // DISPLAYMAINWINDOW_H