# $Header$
###############################################################################
#
# Add module dependencies here; e.g., if in dir 'module-X' there's 
# a file 'myfunc.cc' write:
#
# MODULE_DEPENDENCIES=myfunc.o
#
# Then put this file in directory 'module-X' as 'Makefile.inc'
#
##############################################################################
QTINCLUDE=/usr/include/qt6

MODULE_DEPENDENCIES=display_mainwindow.lo

# MODULE_LINK=-lQtCore \
#             -lQtGui \
#             -lQMainWindow \
#             -lQApplication \
#             -lQtCharts \
#             -ldisplay_mainwindow
MODULE_INCLUDE=-I$(QTINCLUDE) \
               -I$(QTINCLUDE)/QtCore/ \
               -I$(QTINCLUDE)/QtGui/ \
               -I$(QTINCLUDE)/QtWidgets/ \
               -I$(QTINCLUDE)/QtCharts/ -std=c++11
