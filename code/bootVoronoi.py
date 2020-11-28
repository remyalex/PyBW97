import os.path
import sys
import logging
import math
import random
from qgis.core import *
from PyQt4.QtGui import *
app = QApplication([])
QgsApplication.setPrefixPath("/usr", True)
QgsApplication.initQgis()
sys.path.append('/usr/share/qgis/python/plugins')
from processing.core.Processing import Processing
Processing.initialize()
from processing.tools import *
import pdb

class bootVoronoi:
  def __init__ (self, fileMmi, nBoot):
    self.nBoot = int(nBoot)
    self.nInt = 0
    wsPath, iFname = os.path.split(fileMmi)
    print iFname
    self.Mmis = self.readMMIs(wsPath, iFname)
    print self.Mmis
    self.mmiFiles = self.bootstrapSamp(wsPath, iFname)
    print self.mmiFiles
    QgsApplication.exitQgis()
    QApplication.exit()

  def readMMIs(self, wsPath, ifName):
    logging.info('...Leyendo archivo de intensidades: ' + os.path.join(wsPath, ifName))
    if not os.path.exists(os.path.join(wsPath, 'bootLayers')):
        os.mkdir(os.path.join(wsPath, 'bootLayers'))
    pntLy = QgsVectorLayer('Point?crs=epsg:4326&field=x:float&field=y:float&field=i:float', 'intPoints' , "memory")
    f = open(os.path.join(wsPath, ifName), 'r')
    line = f.readline()
    parameters = line.split(';')
    x, y = parameters[0].split('  ')
    self.maxMmiPoint = {'x':float(x), 'y':float(y)}
    self.sqRadious = int(parameters[1])
    self.sqResolution = int(parameters[2])
    line = f.readline()
    pMmiS = list()
    while line:
      data = line.split('\t')
      feat = QgsFeature()
      feat.setFields(pntLy.pendingFields())
      feat.setAttributes([float(data[0]), float(data[1]), float(data[2])])
      feat.setGeometry(QgsGeometry.fromPoint(QgsPoint(float(data[0]), float(data[1]))))
      res, outFeats = pntLy.dataProvider().addFeatures([feat])
      p = {'x':float(data[0]), 'y':float(data[1]), 'i': float(data[2])}
      pMmiS.append(p)
      line = f.readline()
      self.nInt += 1
    save = QgsVectorFileWriter.writeAsVectorFormat(pntLy, os.path.join(wsPath, 'bootLayers', ifName) + ".shp", "UTF-8", pntLy.crs() , "ESRI Shapefile")
    self.pointsLyr =  QgsVectorLayer(os.path.join(wsPath, 'bootLayers', ifName) + ".shp", 'intPoints' , "ogr")
    return pMmiS

  def bootstrapSamp(self, wsPath, ifName):
    logging.info('Iniciando bootstrap con ' + str(self.nBoot) +' remuestreos...')
    print('Iniciando bootstrap con %s remuestreos...', str(self.nBoot))
    mmiFiles = list()
    bootMmiS = list()
    #pdb.set_trace()
    crs = self.pointsLyr.crs().toWkt()
    #voronoiOut = general.runalg("qgis:voronoipolygons", self.pointsLyr, 1, str(os.path.join(wsPath, 'bootLayers', ifName + "_vrn")))
    voronoiOut = general.runalg("qgis:voronoipolygons", self.pointsLyr, 1, None)
    voronoiLyr = QgsVectorLayer(voronoiOut['OUTPUT'], 'intVrn' , "ogr")
    saveVrn = QgsVectorFileWriter.writeAsVectorFormat(voronoiLyr, os.path.join(wsPath, 'bootLayers', ifName + "_vrn") + ".shp", "UTF-8", voronoiLyr.crs() , "ESRI Shapefile")
    for i in range(1,self.nBoot):
        bootFl = open(os.path.join(wsPath, 'bootLayers', ifName + '_boot' + str(i)), 'w')
        bootFl.write(str(self.maxMmiPoint['x'])+'  '+str(self.maxMmiPoint['y'])+';'+str(self.sqRadious)+';'+str(self.sqResolution)+"\n")
        for poly in voronoiLyr.getFeatures():
         bounds = poly.geometry().boundingBox()
         x_coord = random.uniform(bounds.xMinimum(), bounds.xMaximum())
         y_coord = random.uniform(bounds.yMinimum(), bounds.yMaximum())
         while not QgsGeometry.fromPoint(QgsPoint(x_coord, y_coord)).within(poly.geometry()):
             x_coord = random.uniform(bounds.xMinimum(), bounds.xMaximum())
             y_coord = random.uniform(bounds.yMinimum(), bounds.yMaximum())
         p = {'x':x_coord, 'y':y_coord, 'i':poly.attributes()[2]}
         bootMmiS.append(p)
        for p in bootMmiS:
          bootFl.write(str(p['x'])+'\t'+str(p['y'])+'\t'+str(p['i'])+"\n")
        bootFl.close()
        bootMmiS = list()
        mmiFiles.append(os.path.join(wsPath, 'bootLayers', ifName + '_boot' + str(i)))
    mmiFiles.append(os.path.join(wsPath, ifName))
    return mmiFiles
