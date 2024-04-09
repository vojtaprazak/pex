#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:29:18 2017

@author: pratsch
"""
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from Vector import *
from MapParser import MapParser
from PEETPicker import orient_pcle_to_point 
# ususal stuff
from copy import deepcopy
import itertools
import matplotlib.pyplot as plt
from PyQt4 import QtGui # Import the PyQt4 module we'll need
import sys # We need sys so that we can pass argv to QApplication
from PyQt4.QtCore import QThread, SIGNAL # for Threading
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import numpy as np
import SphericalViewer # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer



#Class for the Threading

# XYZ Slicers
class MyWorker_SliceViewerX(QThread):
    # Thread that wraps the x,y,z view
    def __init__(self,myslice, PrepViews):
        QThread.__init__(self)
        self.slice=myslice  
        self. PrepViews= PrepViews

    def __del__(self):
        self.wait()

    def run(self):
        # your logic here        
        self.emit(SIGNAL('draw_xslice(PyQt_PyObject)'),self.PrepViews.get_xslice(self.slice))
        self.emit(SIGNAL('checkforredrawX(PyQt_PyObject)'),self.slice)
        
class MyWorker_SliceViewerY(QThread):
    # Thread that wraps the x,y,z view
    def __init__(self,myslice, PrepViews):
        QThread.__init__(self)
        self.slice=myslice  
        self. PrepViews= PrepViews

    def __del__(self):
        self.wait()

    def run(self):
        # your logic here        
        self.emit(SIGNAL('draw_yslice(PyQt_PyObject)'),self.PrepViews.get_yslice(self.slice))
        self.emit(SIGNAL('checkforredrawY(PyQt_PyObject)'),self.slice)
        
class MyWorker_SliceViewerZ(QThread):
    # Thread that wraps the x,y,z view
    def __init__(self,myslice, PrepViews):
        QThread.__init__(self)
        self.slice=myslice  
        self. PrepViews= PrepViews

    def __del__(self):
        self.wait()

    def run(self):
        # your logic here        
        self.emit(SIGNAL('draw_zslice(PyQt_PyObject)'),self.PrepViews.get_zslice(self.slice))
        self.emit(SIGNAL('checkforredrawZ(PyQt_PyObject)'),self.slice)
        
# Chart viewers        
class MyWorker_S2Viewer(QThread):
    # chart viewers
    def __init__(self,radius,centerx,centery,centerz,PrepViews):
        QThread.__init__(self)
        self.PrepareChartView=PrepareChartView()
        self.radius=radius
        self.center=np.array([centerx,centery,centerz])
        self.PrepViews=PrepViews

    def __del__(self):
        self.wait()

    def run(self):
        # your logic here
        #print('radius:'+str(self.radius))
        X,Y,Z=self.PrepareChartView.getcoordinatesToPolarMesh(self.center,self.radius)
        Chartimage=self.PrepViews.getpoint_xyz(X,Y,Z)
        print('here')
        self.emit(SIGNAL('draw_chart(PyQt_PyObject)'),Chartimage)
        
         
        
       
        
        
class MainAppSphericalViewer(QtGui.QMainWindow, SphericalViewer.Ui_MainWindow):
    # Class for the GUI
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
                            # It sets up layout and widgets that are defined
        
        # Here the interactions are defined        
        self.pushButton_mod.clicked.connect(self.file_open_mod)
        self.pushButton_mrc.clicked.connect(self.file_open_mrc)
        self.pushButton_Start.clicked.connect(self.start_my_calc)
        
        
    # Now the functions used for the interactions are defined    
    def file_open_mod(self):
        name =QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        print(name)        
        self.lineEdit_mod.setText(name)
    
    def file_open_mrc(self):
        name =QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        print(name)        
        self.lineEdit_mrc.setText(name)
        
        
    def start_my_calc(self):
        self.textBrowser_OutputInfo.append('Started computation, please wait \n')
        self.textBrowser_OutputInfo.append('Loading mod-file... \n')
        try: 
            sphericalpts=PEETmodel(self.lineEdit_mod.text()).get_all_points()
        except:
            self.textBrowser_OutputInfo.append('Unexpected error while loading the mod-file... \n')
            raise  
        #self.textBrowser_OutputInfo.append(str(sphericalpts.shape))
        self.textBrowser_OutputInfo.append('Started estimation of sperical parameters... \n')
        try:
            center, radius = self.CHP_lsqSphereFit(sphericalpts)
          
        except:
            self.textBrowser_OutputInfo.append('Unexpected error in estimation of sperical parameters... \n')
            raise
        self.textBrowser_OutputInfo.append('Loading mrc-file... \n')
        try: 
            self.Prepview=PrepareViews(self.lineEdit_mrc.text(),center, radius)
        except:
            self.textBrowser_OutputInfo.append('Unexpected error while loading the mrc-file... \n')
            raise  
        # setting the limits
        max_xyz=self.Prepview.getsize()
        self.doubleSpinBox_centerx.setMaximum(max_xyz[0])
        self.doubleSpinBox_centerx.setValue(center[0])
        self.doubleSpinBox_centery.setMaximum(max_xyz[1])
        self.doubleSpinBox_centery.setValue(center[1])
        self.doubleSpinBox_centerz.setMaximum(max_xyz[2])
        self.doubleSpinBox_centerz.setValue(center[2])
        self.doubleSpinBox_radius.setMaximum(np.max([2.0*radius,200.]))
        self.doubleSpinBox_radius.setValue(radius)    
        #same for the slicer
        self.spinBox_slicex.setMaximum(max_xyz[0])
        self.horizontalScrollBar_slicex.setMaximum(max_xyz[0])
        self.spinBox_slicey.setMaximum(max_xyz[1])
        self.horizontalScrollBar_slicey.setMaximum(max_xyz[1])
        self.spinBox_slicez.setMaximum(max_xyz[2])
        self.horizontalScrollBar_slicez.setMaximum(max_xyz[2])
        self.spinBox_slicex.setValue(int(center[0]))
        self.spinBox_slicey.setValue(int(center[1]))
        self.spinBox_slicez.setValue(int(center[2]))
        
        
        self.textBrowser_OutputInfo.append('Init of grafics... \n')
        # TODO The surface view
        self.myChartThread =  MyWorker_S2Viewer(self.doubleSpinBox_radius.value(),self.doubleSpinBox_centerx.value(),
                                                self.doubleSpinBox_centery.value(),self.doubleSpinBox_centerz.value(),self.Prepview)                                       
        self.connect(self.myChartThread,SIGNAL('draw_chart(PyQt_PyObject)'),self.draw_chart)
        self.myChartThread.start()
        # These threads generate the slices
        # These variables ensure that only one thread of each kind is active
        self.myxSliceThread =  MyWorker_SliceViewerX(self.spinBox_slicex.value(),self.Prepview)                                       
        self.connect(self.myxSliceThread, SIGNAL("draw_xslice(PyQt_PyObject)"), self.draw_xslice) # Does the add_post function and is called by signal from thread
        self.connect(self.myxSliceThread, SIGNAL('checkforredrawX(PyQt_PyObject)'),self.restartThreadX)
        self.connect(self.myxSliceThread, SIGNAL('finished()'),self.restartThreadX_fin)
        self.myxSliceThread.start()
        self.myySliceThread =  MyWorker_SliceViewerY(self.spinBox_slicey.value(),self.Prepview)                                       
        self.connect(self.myySliceThread, SIGNAL("draw_yslice(PyQt_PyObject)"), self.draw_yslice) # Does the add_post function and is called by signal from thread
        self.connect(self.myySliceThread, SIGNAL('checkforredrawY(PyQt_PyObject)'),self.restartThreadY)
        self.connect(self.myySliceThread, SIGNAL('finished()'),self.restartThreadY_fin)
        self.myySliceThread.start()
        self.myzSliceThread =  MyWorker_SliceViewerZ(self.spinBox_slicez.value(),self.Prepview)                                       
        self.connect(self.myzSliceThread, SIGNAL("draw_zslice(PyQt_PyObject)"), self.draw_zslice) # Does the add_post function and is called by signal from thread
        self.connect(self.myzSliceThread, SIGNAL('checkforredrawZ(PyQt_PyObject)'),self.restartThreadZ)
        self.connect(self.myzSliceThread, SIGNAL('finished()'),self.restartThreadZ_fin)
        self.myzSliceThread.start()
        self.spinBox_slicex.valueChanged.connect(self.updateGx)
        self.spinBox_slicey.valueChanged.connect(self.updateGy)
        self.spinBox_slicez.valueChanged.connect(self.updateGz)
        self.doubleSpinBox_centerx.valueChanged.connect(self.updateChart)
        self.doubleSpinBox_centery.valueChanged.connect(self.updateChart)
        self.doubleSpinBox_centerz.valueChanged.connect(self.updateChart)
        self.doubleSpinBox_radius.valueChanged.connect(self.updateChart)
        self.textBrowser_OutputInfo.append('Have fun... \n')
        
        
   
    def rmmpl(self,): #not used
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()
        
          
    def draw_xslice(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        fig1 = Figure()
        self.ax1fx = fig1.add_subplot(111)
        self.ax1fx.imshow(imagarray,'gray')                     
        self.canvasx = FigureCanvas(fig1)
        self.verticalLayout_grafic_x.addWidget(self.canvasx)
        self.canvasx.draw()
        self.toolbarx = NavigationToolbar(self.canvasx, 
                self.widget_grafic_x, coordinates=True)
        self.verticalLayout_grafic_x.addWidget(self.toolbarx)
    
    def draw_xslice_rerun(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        self.ax1fx.imshow(imagarray,'gray')     
        self.canvasx.draw()        
    
    def draw_yslice(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        figy = Figure()
        self.ax1fy = figy.add_subplot(111)
        self.ax1fy.imshow(imagarray,'gray')                     
        self.canvasy = FigureCanvas( figy)
        self.verticalLayout_grafic_y.addWidget(self.canvasy)
        self.canvasy.draw()
        self.toolbary = NavigationToolbar(self.canvasy, 
                self.widget_grafic_y, coordinates=True)
        self.verticalLayout_grafic_y.addWidget(self.toolbary)

    def draw_yslice_rerun(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        self.ax1fy.imshow(imagarray,'gray')     
        self.canvasy.draw()        
        
    def draw_zslice(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        fig1 = Figure()
        self.ax1fz = fig1.add_subplot(111)
        self.ax1fz.imshow(imagarray,'gray')                     
        self.canvasz = FigureCanvas(fig1)
        self.verticalLayout_grafic_z.addWidget(self.canvasz)
        self.canvasz.draw()
        self.toolbarz = NavigationToolbar(self.canvasz, 
                self.widget_grafic_z, coordinates=True)
        self.verticalLayout_grafic_z.addWidget(self.toolbarz)
        
    def draw_zslice_rerun(self, imagarray):
        """
        This function is called by the thread emited signal
        it gets a string argument back from the thread! 
        """
        self.ax1fz.imshow(imagarray,'gray')     
        self.canvasz.draw()            

    # Things that interact with the thread or draw due to threads
    def updateGx(self):
        if self.myxSliceThread.isRunning()==False:
            self.myxSliceThread =  MyWorker_SliceViewerX(self.spinBox_slicex.value(),self.Prepview)                                       
            self.connect(self.myxSliceThread, SIGNAL("draw_xslice(PyQt_PyObject)"), self.draw_xslice_rerun) # Does the add_post function and is called by signal from thread
            self.myxSliceThread.start()   
        else:
            pass

    def updateGy(self):
        if self.myySliceThread.isRunning()==False:
            self.myySliceThread =  MyWorker_SliceViewerY(self.spinBox_slicey.value(),self.Prepview)                                       
            self.connect(self.myySliceThread, SIGNAL("draw_yslice(PyQt_PyObject)"), self.draw_yslice_rerun) # Does the add_post function and is called by signal from thread
            self.myySliceThread.start()        
        else:
            pass
    def updateGz(self):
        if self.myzSliceThread.isRunning()==False:
            self.myzSliceThread =  MyWorker_SliceViewerZ(self.spinBox_slicez.value(),self.Prepview)                                       
            self.connect(self.myzSliceThread, SIGNAL("draw_zslice(PyQt_PyObject)"), self.draw_zslice_rerun) # Does the add_post function and is called by signal from thread
            self.myzSliceThread.start()   
        else:
            pass
    
    def restartThreadX(self, xslice):
        self.TODOx=(xslice!=self.spinBox_slicex.value())
  
    def restartThreadX_fin(self):
        if self.TODOx:
            self.updateGx()
        else:
            pass
        
    def restartThreadY(self, yslice):
            self.TODOy=(yslice!=self.spinBox_slicey.value())
        
    def restartThreadY_fin(self):
        if self.TODOy:
            self.updateGy()
        else:
            pass
        
    def restartThreadZ(self, zslice):
            self.TODOz=(zslice!=self.spinBox_slicez.value())
        
    def restartThreadZ_fin(self):
        if self.TODOz:
            self.updateGz()
        else:
            pass    
        
    #Chart related stuff
    def draw_chart(self, imagarray):
        print('washere')
        fig1 = Figure()
        self.ax1fS2C = fig1.add_subplot(111)
        self.ax1fS2C.imshow(imagarray,'gray')                     
        self.canvasS2C = FigureCanvas(fig1)
        self.verticalLayout_Chart.addWidget(self.canvasS2C)
        self.canvasS2C.draw()
        self.toolbarS2C = NavigationToolbar(self.canvasS2C, 
                self.myGraficwidget_Chart, coordinates=True)
        self.verticalLayout_Chart.addWidget(self.toolbarS2C)

    def redraw_chart(self, imagarray):
        print('washere')
        self.ax1fS2C.imshow(imagarray,'gray')                     
        self.canvasS2C.draw()
        
        
        
    def updateChart(self):
        if self.myChartThread.isRunning()==False:
                self.myChartThread =  MyWorker_S2Viewer(self.doubleSpinBox_radius.value(),self.doubleSpinBox_centerx.value(),
                                                self.doubleSpinBox_centery.value(),self.doubleSpinBox_centerz.value(),self.Prepview)                                       
                self.connect(self.myChartThread,SIGNAL('draw_chart(PyQt_PyObject)'),self.redraw_chart)
                self.myChartThread.start() 
        else:
            pass
        
    # Coordinates related Stuff
    def CHP_lsqSphereFit(self,startlist):
        """
        Idea stolen from Matlab Central
        """
        x=startlist[:,0]
        y=startlist[:,1]
        z=startlist[:,2]
        A=np.ones([np.shape(x)[0],4])
        A[:,0]=x
        A[:,1]=y
        A[:,2]=z 
        b=-(x**2 + y**2 + z**2)     
        outp=np.linalg.lstsq(A,b) #Matlabs: A \ b
        center = -outp[0][0:3]/2
        radius = np.sqrt(np.sum(center**2)-outp[0][3])
        return center, radius 
    
class PrepareChartView ():
    def __init__(self):
        pass
    
    
    def getcoordinatesToPolarMesh(self,center,radius,phi=0,theta=0):
        #print('here radius: '+str(radius))
        phigrid, thetagrid = np.mgrid[ 0:2 * np.pi:401j,0:np.pi:201j]
        x=np.sin(thetagrid)*np.cos(phigrid)*radius+center[0]
        y=np.sin(thetagrid)*np.sin(phigrid)*radius+center[1]
        z=np.cos(thetagrid)*radius+center[2]
        return x,y,z

class PrepareViews ():
    # Class that stores the tomogram and gets the data
    def __init__(self, filename_mrc,center, radius):
        # filename_mrc is the name of the mrc file to be viewed, center and radius are parameters of the object
        try: 
            self.MRC=MapParser.readMRC(filename_mrc).fullMap# reads in z,y, x            
        except:
            print(("Unexpected error in mrc file:", sys.exc_info()[0]))
            raise   
    
    def get_zslice(self,z):
        try:
            zSlice=self.MRC[z,:,:]
        except:
            print('zSlice out of range')
            raise
            pass
        return zSlice
        
    def get_yslice(self,y):
        try:
            ySlice=self.MRC[:,y,:]
        except:
            print('yzSlice out of range')
            raise            
            pass
        return ySlice
        
    def get_xslice(self,x):
        try:
            xSlice=self.MRC[:,:,x]
        except:
            print('xSlice out of range')
            raise
            pass
        return xSlice
        
    def getpoint_xyz(self,x,y,z):  
        try:
           # removing things outside of ROI
           z=np.int32(z)
           tz=z>=self.MRC.shape[0]
           z[tz]=0
           z[z<0]=0
           y=np.int32(y)
           ty=y>=self.MRC.shape[1]
           y[ty]=0
           y[y<0]=0
           x=np.int32(x)
           tx=x>=self.MRC.shape[2]
           x[tx]=0
           x[x<0]=0 
            
           pointxyz=self.MRC[np.int32(z),np.int32(y),np.int32(x)] 
           pointxyz[tz]=0
           pointxyz[ty]=0
           pointxyz[tx]=0
                   
           print('that did work')
        except:
                print('that did not work')
                raise
                pass
        return pointxyz
    
    def getsize(self):
        sxyz=self.MRC.shape
        return sxyz[2]-1,sxyz[1]-1,sxyz[0]-1
        
def main():
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication
    form = MainAppSphericalViewer()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    sys.exit(app.exec_())                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()                              # run the main function
