import os.path
import logging
import math
import random
import copy
import pdb

class bootAffine:
  def __init__ (self, fileMmi, nBoot):
    self.nBoot = int(nBoot)
    self.nInt = 0
    self.Mmis = self.readMMIs(fileMmi)
    wsPath, iFname = os.path.split(fileMmi)
    self.mmiFiles = self.bootstrapSamp(wsPath, iFname)

  def readMMIs(self, path):
      logging.info('...Leyendo archivo de intensidades: ' + str(path) + '...')
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

  def bootstrapSamp(self, wsPath, ifName):
      logging.info('Iniciando bootstrap con ' + str(self.nBoot) +' remuestreos...')
      print('Iniciando bootstrap con %s remuestreos...', str(self.nBoot))
      if not os.path.exists(os.path.join(wsPath, 'bootSamples')):
          os.mkdir(os.path.join(wsPath, 'bootSamples'))
      mmiFiles = list()
      bootMmiS = list()
      for i in range(1,self.nBoot):
          bootFl = open(os.path.join(wsPath, 'bootSamples', ifName + '_boot' + str(i)), 'w')
          bootFl.write(str(self.maxMmiPoint['x'])+'  '+str(self.maxMmiPoint['y'])+';'+str(self.sqRadious)+';'+str(self.sqResolution)+"\n")
          while len(bootMmiS) < len(self.Mmis):
              p = random.choice(self.Mmis)
              copyP = copy.copy(p)
              if bootMmiS.count(p) == 0:
                  bootMmiS.append(p)
              else:
                  dL = self.maxMmiPoint['x'] - p['x']
                  dP = self.maxMmiPoint['y'] - p['y']
                  while not bootMmiS.count(copyP) == 0:
                      #pdb.set_trace()
                      t = random.uniform(0, 360)
                      dP_ = dL*math.sin(math.radians(t)) + dP*math.cos(math.radians(t))
                      dL_ = dL*math.sin(math.radians(t)) - dP*math.sin(math.radians(t))
                      copyP['x'] = self.maxMmiPoint['x'] - dL_
                      copyP['y'] = self.maxMmiPoint['y'] - dP_
                  bootMmiS.append(copyP)
          for p in bootMmiS:
            bootFl.write(str(p['x'])+'\t'+str(p['y'])+'\t'+str(p['i'])+"\n")
          bootFl.close()
          bootMmiS = list()
          mmiFiles.append(os.path.join(wsPath, 'bootSamples', ifName + '_boot' + str(i)))
      mmiFiles.append(os.path.join(wsPath, ifName))
      return mmiFiles
