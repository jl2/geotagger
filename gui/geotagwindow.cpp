/*
  geotagwindow.h
  
  Copyright (c) 2011, Jeremiah LaRocco jeremiah.larocco@gmail.com

  Permission to use, copy, modify, and/or distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <QtGui>
#include <stdexcept>
#include "geotagwindow.h"

void GeoTagWindow::readSettings() {
  settings = new QSettings(QSettings::IniFormat, QSettings::UserScope,
                           "GeoTag", "GeoTag");
    
}

GeoTagWindow::GeoTagWindow(QWidget *parent) : QMainWindow(parent) {
    readSettings();

    setupActions();
    setupToolBar();
    setupMenuBar();
    setupStatusBar();

    setWindowTitle(tr("GeoTag"));
}

GeoTagWindow::~GeoTagWindow() {
    settings->sync();
}

void GeoTagWindow::setupActions() {

    // Exit
    exitAction = new QAction(tr("Exit"), this);
    exitAction->setIcon(QIcon(":/images/exit.png"));
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));//exit()));

    // About
    aboutAction = new QAction(tr("About"), this);
    aboutAction->setIcon(QIcon(":/images/about.png"));
    aboutAction->setStatusTip(tr("About this geotagging program"));
    connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));

}

void GeoTagWindow::setupToolBar() {
    QToolBar *tb = new QToolBar("Main", this);
    
    addToolBar(tb);
}

void GeoTagWindow::setupMenuBar() {
    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(exitAction);

    menuBar()->addSeparator();
  
    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAction);
}

void GeoTagWindow::setupStatusBar() {
    statusBar()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

}

void GeoTagWindow::about() {
    QMessageBox::about(this,
                       tr("About"),
                       tr("<h2>GeoTag!</h2>"
                          "<p>By Jeremiah LaRocco.</p>"));
}
