import os.path
import logging
import numpy
import math
import argparse
from datetime import datetime
from bootAffine import bootAffine
from bootVoronoi import bootVoronoi

class mmiMesh:
    def __init__ (self, fileMmi, md):
        #Funcion leer parametros y datos
        #numero de puntos de intensidad
        self.nInt = 0
        self.minRms = 100
        self.Mmis = self.readMMIs(fileMmi)
        #Comprobar si Radio de malla es mayor que la resolucion
        if(self.sqRadious>self.sqResolution):
            #Calcula numero de elementos en ejes de malla
            self.n = int(2*math.ceil(self.sqRadious/self.sqResolution))
            #Calculo de la resolucion grados decimales
            self.Dg = float(self.sqResolution)/6378*float(180)/math.pi
            #Llamado a procedimentos
            self.mesh = self.makeMesh(self.maxMmiPoint)
            self.A, self.B, self.C, self.D = self.choseModel(md)
            self.magn = self.computeMagnitudes()
            self.rms = self.computeRMS()
            #self.dists = self.computeDistances(maxMmiPoint)
        else:
            logging.error('La resolucion excede el area de la malla.')
            print('La resolucion excede el area de la malla.')

    def makeMesh(self, p):
        #Calculo de los puntos inicial y final de cada eje (XY)
        x0n = p['x']+(self.Dg*self.n/2)
        xn0 = p['x']-(self.Dg*self.n/2)
        y0n = p['y']+(self.Dg*self.n/2)
        yn0 = p['y']-(self.Dg*self.n/2)
        #Construir ejes intervalos
        xaxis = numpy.linspace(xn0,x0n,self.n)
        yaxis = numpy.linspace(yn0,y0n,self.n)
        #Delvolver valores de malla
        return numpy.meshgrid(xaxis, yaxis)

    def choseModel(self, mod):
        if mod == 1:
            return 2.2237, 1.6684, 0.041214, 0
        elif mod == 2:
            return 1.92, 2.33, 0.0021, 3.68

    def computeDistances(self, p):
        #Solo una colatitud para punto centro
        cl_e = (90 - p['y']) * math.pi/180
        #Evalua todas las colatitudes de la malla
        cl_p = (90 - self.mesh[1]) * math.pi/180
        #Evalua todas las DeltaLong de la malla
        dLong = numpy.fabs(p['x'] - self.mesh[0]) * math.pi/180
        #Evalua distancias esfericas
        distEsf = numpy.arccos(numpy.cos(cl_e)*numpy.cos(cl_p)+numpy.sin(cl_e)*numpy.sin(cl_p)*numpy.cos(dLong))
        #Retorna distancia en Metros (Radio de la Tierra X Distancia Esferica)
        return distEsf*6378.137

    def singleMi(self, inte, pp):
        Dist = self.computeDistance(pp, inte)
        LogDist = math.log10(Dist)
        return (1/self.B * (inte['i'] + self.A + self.C*(Dist) + self.D*(LogDist)))

    def computeMagnitudes(self):
        a = numpy.zeros(shape=(self.n, self.n))
        logging.info('...Calculando magnitudes ...')
        print('...Calculando magnitudes ...')
        for i in range(0,self.n):
            for j in range(0,self.n):
                sumMagn = 0
                for inte in self.Mmis:
                    pp = {'x':self.mesh[0][i,j], 'y':self.mesh[1][i,j]}
                    sumMagn += self.singleMi(inte, pp)
                a[i,j] = sumMagn/self.nInt
        return a

    def computeWeight(self, dist):
        if dist < 150:
            weight = 0.1 + math.cos((dist/150)*math.pi/2)
        else:
            weight = 0.1
        return weight

    def computeRMS(self):
        rms = numpy.zeros(shape=(self.n, self.n))
        logging.info('...Calculando incertidumbres ...')
        print('...Calculando incertidumbres ...')
        for i in range(0,self.n):
            for j in range(0,self.n):
                sumRMS = 0
                sumWeight = 0
                for inte in self.Mmis:
                    pp = {'x':self.mesh[0][i,j], 'y':self.mesh[1][i,j]}
                    weight = self.computeWeight(self.computeDistance(pp, inte))
                    sumWeight += math.pow(weight,2)
                    sumRMS += weight*math.pow((self.singleMi(inte, pp) - self.magn[i,j]), 2)
                rms[i,j] = math.sqrt((1/sumWeight)*sumRMS)
        for i in range(0,self.n):
            for j in range(0,self.n):
                if self.minRms > rms[i,j]:
                    self.minRms = rms[i,j]
        for i in range(0,self.n):
            for j in range(0,self.n):
                rms[i,j] = rms[i,j] - self.minRms
        return rms

    def computeDistance(self, p1, p2):
        #Solo una colatitud para punto centro
        cl_e = (90 - p1['y']) * math.pi/180
        #Evalua todas las colatitudes de la malla
        cl_p = (90 - p2['y']) * math.pi/180
        #Evalua todas las DeltaLong de la malla
        dLong = math.fabs(p1['x'] - p2['x']) * math.pi/180
        #Evalua distancias esfericas
        distEsf = math.acos(math.cos(cl_e)*math.cos(cl_p)+math.sin(cl_e)*math.sin(cl_p)*math.cos(dLong))
        #Retorna distancia en Metros (Radio de la Tierra X Distancia Esferica)
        return distEsf*6378.137

    def readMMIs(self, path):
        logging.info('...Leyendo archivo de intensidades: %s...', path)
        print('...Leyendo archivo de intensidades: %s...', path)
        f = open(path, 'r')
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
            p = {'x':float(data[0]), 'y':float(data[1]), 'i': float(data[2])}
            pMmiS.append(p)
            line = f.readline()
            self.nInt += 1
        return pMmiS

def main(intsFile, modelNum, bootsGroups=0, magnFlName=None, rmsFlName=None, rasterName=None):
    if not magnFlName:
        magnFlName = 'magnitudes'
    if not rmsFlName:
        rmsFlName = 'rms'
    wsPath, iFname = os.path.split(intsFile)
    logging.basicConfig(filename=os.path.join(wsPath,'computeLog_'+datetime.now().strftime('%Y%m%d-%H%M%S')+'.log',), level=logging.DEBUG)
    if bootsGroups:
        if rasterName:
            mmiFiles = bootVoronoi(intsFile, bootsGroups).mmiFiles
        else:
            mmiFiles = bootAffine(intsFile, bootsGroups).mmiFiles
    else:
        mmiFiles = [intsFile]
    logging.info('Iniciando calculos...')
    meshList = list()
    mesh = None
    for ints in mmiFiles:
        meshList.append(mmiMesh(ints, int(modelNum)))
    if len(mmiFiles) > 1:
        bootFl = open(os.path.join(wsPath, "bootOutput_"+datetime.now().strftime('%Y%m%d-%H%M%S')), 'w')
        bootFl.write("Long;Lat;Magn;minRMS\n")
        magn = numpy.zeros(shape=(max(meshList).n, max(meshList).n))
        rms = numpy.zeros(shape=(max(meshList).n, max(meshList).n))
        for i in range(0,max(meshList).n):
            for j in range(0,max(meshList).n):
                sumMagn = 0
                sumRms = 0
                for mesh in meshList:
                    if mesh.rms[i,j] == 0:
                        bootFl.write(str(mesh.mesh[0][i,j])+";"+str(mesh.mesh[1][i,j])+";"+str(mesh.magn[i,j])+";"+str(mesh.minRms)+"\n")
                    sumMagn += mesh.magn[i,j]
                    sumRms += mesh.rms[i,j]
                magn[i,j] = sumMagn/len(meshList)
                rms[i,j] = sumRms/len(meshList)
        oMesh = max(meshList).mesh
        oMagn = magn
        oRms = rms
        bootFl.close()
    else:
        mesh = min(meshList)
        oMesh = mesh.mesh
        oMagn = mesh.magn
        oRms = mesh.rms

    magnFl = open(os.path.join(wsPath, magnFlName), 'w')
    magnFl.write("Long;Lat;Magn\n")
    rmsFl = open(os.path.join(wsPath, rmsFlName), 'w')
    rmsFl.write("Long;Lat;RMS\n")
    for i in range(0,mesh.n):
        for j in range(0,mesh.n):
            magnFl.write(str(oMesh[0][i,j])+";"+str(oMesh[1][i,j])+";"+str(oMagn[i,j])+"\n")
            rmsFl.write(str(oMesh[0][i,j])+";"+str(oMesh[1][i,j])+";"+str(oRms[i,j])+"\n")
            if mesh.rms[i,j] == 0:
                logging.info("Epicentro: " + str(oMesh[0][i,j])+";"+str(oMesh[1][i,j])+". Magnitud:"+str(oMagn[i,j])+"+/-"+str(oRms[i,j]))
                print "Epicentro: " + str(oMesh[0][i,j])+";"+str(oMesh[1][i,j])+". Magnitud:"+str(oMagn[i,j])+"+/-"+str(oRms[i,j])
    magnFl.close()
    rmsFl.close()
    logging.info("Error de ajuste: " + str(mesh.minRms))
    print "Error de ajuste: " + str(mesh.minRms)

if __name__ == "__main__":
    epilog = ('***PyBW97 v1.0.0 [Pyhon 2.7]***\n'
              'Bakun & Wenworth Method for Earthquake Physical parameters estimation Program\n'
              ' Computes Physical parameters from Earthquake Intensities file in sample output format:\n'
              '     *** Intensities File Format ***\n'
              '         Head [Coordinates of Intensities Center and Grid Parameters]: <Long> <Lat>;<GridRadious(kms)>;<GridStep(kms)>\n'
              '         Next lines [Coordinates and intensity value (EMS98 scale)]: <Long> <Lat> <Intensity(ems98)>\n'
              '     *** End (<filename>.int)***\n'
              '--- Required Parameters ---\n'
              '     -if Path of intensities file: /<path>/<to>/<filename>.int\n'
              '     -m IPE Models (Intensity Prediction Equation):\n'
              '         <1> - Palme et al (2005)       :  I = 2.2 - 1.6*M + 4e-2*(Dist) + 0*(LogDist)\n'
              '         <2> - Gomez-Capera et al (2016):  I = 1.92 - 2.3*M + 2.1e-3*(Dist) + 3.68*(LogDist)\n'
              '--- End ---\n'
              '--- Optional Parameters ---\n'
              '     -mn Name of magnitudes file: <magnitudesName>\n'
              '     -rn Name of rms file: <rmsName>\n'
              '     -b  Number of Bootstrap resamples\n'
              '--- End ---\n\n'
              'Author: Remy Galan\n'
              'Licence: GNU GPLv3\n')
    parser = argparse.ArgumentParser(epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-if", "--intsFile",
                        required=True,
                        help="Entire path of EQ macro-intensities file in <int> format (Above sample - Path will be PyBW97 workspace)")
    parser.add_argument("-m", "--modelNum",
                        required=True,
                        help="Choose model number, from the above listed")
    parser.add_argument("-b", "--bootsGroups",
                        required=False,
                        help="Set number of bootstrap resamples")
    parser.add_argument("-mn", "--magnFlName",
                    required=False,
                    help="Name for magnitudes output-file")
    parser.add_argument("-rn", "--rmsFlName",
                        required=False,
                        help="Name for rms output-file")
    parser.add_argument("-br", "--rasterName",
                        required=False,
                        help="Name for input raster intensities file")
    parsed_args = parser.parse_args()
    main(**vars(parsed_args))
