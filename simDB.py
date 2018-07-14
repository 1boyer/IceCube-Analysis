#! /usr/bin/python
#      simDB
# This program simulates the occurrence of double
# bang flasher events within IceCube. It requires a geometry
# file be provided by simDBGeometry.py. This program stores
# its results in simDB.pickle for further use.
#
#   Version 1

import sys
import math
import random
import time
import pickle

  
try:
  random.seed(time.time())
except:
  print "Can't seed random number generator! Exiting!"
  exit(255)
class module:
  def __init__(self, Name, omRate=10):
    self.name = Name
    self.rate = omRate
    self._period = (10**9)/omRate
    self.next = random.randint(150, int(self._period)) 
    self.prev = 0
    if omRate == 10:
      self._meanDev = 4857757 # ns
      self._stdDev = 85 # ns
    else:
      self._meanDev = 0
      self._stdDev = 0
    self._period = self._period + random.gauss(self._meanDev, self._stdDev)
    
  def update(self, newTime): #Update list in time window
    if newTime > self.next:
      self.prev = self.next
      self.next = self.next + self._period


class daq:
  def __init__(self, geometry, hours=1, omNumber=10, updateRate=10, logfilename=None):
    
    self.maxTime = hours*3600*(10**9)
    self.omNum = omNumber
    self.mimics = 0
    self._mimicStamps = []
    self.daqTime = 0
    self.updatePeriod = 5e7
    self.omList = []
    self.positions = geometry
    self._window = 110.0
    
    if logfilename:
      try:
        self.logFile = open(logfilename, 'w')
      except:
        self.logFile = None
    else:
      self.logFile = None
      
    for key in random.sample( self.positions.keys() , self.omNum):
      self.omList.append( module( key, 10) )

  def tick(self):
    self.daqTime = self.daqTime + self.updatePeriod
    for om in self.omList:
      om.update(self.daqTime)
    if self.daqTime > self.maxTime:
      return False
    else: 
      return True

  def findMimics(self):
    for i in range( self.omNum ):
      for j in range( i ):
        self.mimicCheck(self.omList[i], self.omList[j])

  def recordMimic(self, stamp):
    if stamp not in self._mimicStamps:
      self._mimicStamps.append( stamp )
      self.mimics = self.mimics + 1
      if self.logFile:
        self.logFile.write( "Mimic between %s and %s @ %d\n" % (str(self.omList[i]), str(self.omList[j]), self.daqTime) )

  def mimicCheck(self, om1, om2 ):
    X = self.positions[om1.name][0] - self.positions[om2.name][0]
    Y = self.positions[om1.name][1] - self.positions[om2.name][1]
    Z = self.positions[om1.name][2] - self.positions[om2.name][2]
    spaceDistance = math.sqrt(X**2 + Y**2 + Z**2)
    spaceDistanceLow = spaceDistance - self._window
    spaceDistanceHigh = spaceDistance + self._window
    
    t00 = math.fabs(om1.next - om2.next)
    t10 = math.fabs(om1.prev - om2.next)
    t01 = math.fabs(om1.next - om2.prev)
    
    if spaceDistanceLow < t00 and t00 < spaceDistanceHigh:
      self.recordMimic( om1.name, om1.next, om2.name, om2.next )
    elif spaceDistanceLow < t01 and t01 < spaceDistanceHigh:
      self.recordMimic( om1.name, om1.prev, om2.name, om2.next )
    elif spaceDistanceLow < t10 and t10 < spaceDistanceHigh:
      self.recordMimic( om1.name, om1.next, om2.name, om2.prev )

  def run(self):
    run = True
    while run:
      self.findMimics()
      run = self.tick()
    if self.logFile:
      self.logFile.write( str(self) )

  def __str__(self):
    retStr = "Mimics found %d, Time Elapsed %d hrs\n" % \
      (self.mimics, self.daqTime/(3600*10**9))
    retStr = retStr + "Active OMs during Trial (%s)\n" % len(self.omList)
    for om in self.omList:
      retStr = retStr + "OM:" + om.name + "\tRate:" + str(om.rate) + "\n"
    return retStr


if __name__ == "__main__":

  if len(sys.argv) > 1 and sys.argv[1].isdigit():
    numOMs = int(sys.argv[1])
  else:
    numOMs = random.randint(1,80)

  try:  
    geoFile = open("simDB_geometry.pickle", 'r')
    geometry = pickle.load(geoFile)
    geoFile.close()
  except:
    print "Run simDBGeometry.py first! Exiting."
    sys.exit(255)
    
          #hours=1, omNumber=10, updateRate=10):
  simulation = daq(geometry, 4, numOMs, 10)
  print "Beginning simulation (%s OMs)...." % numOMs, 
  simulation.run()
  print "finished in %d processor seconds" % time.clock()
 
  try:
    infile = open("simDB.pickle", 'r')
    results = pickle.load(infile)
    infile.close()
  except:
    results = []
    print "No pickle file found... creating new pickle file"

  results.append( (simulation.mimics, simulation.omNum) )
  outfile = open("simDB.pickle", 'w')
  pickle.dump(results, outfile)
