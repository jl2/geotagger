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

#ifndef GEO_TAG_WINDOW_INCLUDE_H
#define GEO_TAG_WINDOW_INCLUDE_H

#include <QMainWindow>

class QToolBar;
class QMenu;
class QAction;
class QLabel;
class QFrame;
class QSettings;

class GeoTagWindow : public QMainWindow {
    Q_OBJECT;

    QToolBar *toolBar;
    QMenu *mainMenu;

public:
    GeoTagWindow(QWidget *parent = 0);
    ~GeoTagWindow();
  
public slots:
    void about();

/* private slots: */
/*     void resetView(); */
    // void pluginChanged();
/*     void startEvolving(); */
/*     void stopEvolving(); */
/*     void reset(); */
private:
    void readSettings();
    void setupActions();
    void setupFileActions();
    void setupEditActions();
    void setupToolsActions();
    void setupHelpActions();

    void setupToolBar();
    void setupMenuBar();
    void setupStatusBar();

    // void loadPlugins();

    // LifeWidget *life;
  
    // QAction *startAction;
    // QAction *stopAction;
    // QAction *resetAction;
    // QAction *configureAction;

    // QAction *resetViewAction;
    QAction *exitAction;
    QAction *aboutAction;

    // QLabel *curIterLabel;
    // QLabel *curPluginLabel;

    // QMap<QString, LifePlugin *> plugins;
    // QString curPlugin;

    QSettings *settings;
};

#endif
