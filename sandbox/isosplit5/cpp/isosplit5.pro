QT += core
QT -= gui

CONFIG += c++11

DESTDIR = ../bin
OBJECTS_DIR = ../build
TARGET = isosplit5
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    isocut5.cpp \
    jisotonic.cpp \
    isosplit5.cpp

HEADERS += \
    isocut5.h \
    jisotonic.h \
    isosplit5.h
