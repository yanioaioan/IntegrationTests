TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    OneDSolver.cpp

include(deployment.pri)
qtcAddDeployment()


#in

HEADERS += \
    OneDSolver.h
CONFIG +=c++11
