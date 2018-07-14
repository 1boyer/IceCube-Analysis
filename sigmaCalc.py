#!/usr/bin/env python

# (c) 2010
# The IceCube Collaboration
# $Id: tauFinder $
#
# $Revision: 1$
# $Date: June 9, 2010$
# Lance Boyer
##############################

import sys
import os
import types
import math
from optparse import OptionParser, OptionValueError
from decimal import *

try:
  from icecube.icetray import *
  from icecube.dataclasses import *
  from icecube.dataio import *
except ImportError:
  print "IceTray is missing! This script needs IceTray!"
  sys.exit(1)


class digitalOpticalModule: #Tracks the Period Between Flashes
  def __init__(self, omKey, flashRate, StartTime=None):
    self.name = omKey
    self.rate = flashRate
    self._nominalPeriod = 10**9/flashRate
    self._lastFlashTime = StartTime
    self.periods = []
    self.deviations = []

  def addFlash(self, time):
    flashPeriod = time - self._lastFlashTime
    self._lastFlashTime = time
    self.periods.append( flashPeriod )
    self.deviations.append ( flashPeriod - self._nominalPeriod )

  def printStatistics(self):
    meanPeriod = self.getMeanPeriod()
    periodStandardDeviation = self.getPeriodStandardDeviation()
    meanDeviation = self.getMeanDeviation()
    deviationStandardDeviation = self.getDeviationStandardDeviation()
    
    print "Name: ", self.name, "  Nominal Period: ", self._nominalPeriod
    print "Mean Period: ", meanPeriod, "(%f)" % ( float(meanPeriod)/self._nominalPeriod)
    print "Period Standard Deviation: ", periodStandardDeviation
    print "Mean Deviation: ", meanDeviation, "(%f)" % (float(meanDeviation)/self._nominalPeriod)
    print "Deviation Standard Deviation: ", deviationStandardDeviation
    print 

  def getMeanPeriod(self):
    return sum(self.periods) / len(self.periods)
  def getPeriodStandardDeviation(self):
    return math.sqrt( math.fabs( sum([period**2 for period in self.periods])/len(self.periods) - self.getMeanPeriod()**2) )
  def getMeanDeviation(self):
    return sum(self.deviations) / len(self.deviations)
  def getDeviationStandardDeviation(self):
    return math.sqrt( math.fabs( sum([deviation**2 for deviation in self.deviations])/len(self.deviations) - self.getMeanDeviation()**2) )


class flashReader:

  def __init__(self, i3FilenameList=None):
    self.omList = {}
    self.startTime = None

    if type(i3FilenameList) == types.NoneType:
      i3FilenameList = []
    elif type(i3FilenameList) == types.StringType:
      i3FilenameList = [i3FilenameList]

    for i3Filename in i3FilenameList:
      self.processFile(i3Filename)

  def processFile(self,i3Filename):
    try:
      file = I3File(i3Filename)
    except:
      (type, val, tback) = sys.exc_info()
      print >> sys.stderr, "%s: Error Opening I3 file:" % sys.argv[0]
      print >> sys.stderr, val
      try:
        errno_val = val.errno
      except:
        return
    try:
      frame = file.pop_frame()
      while frame:
        self.frameRead(frame)
        frame = file.pop_frame()
    finally:
      file.close()


  def frameRead(self, frame):
    if "I3EventHeader" and "flasher" in frame:
      eventHeader = frame["I3EventHeader"]
      frameStartTime = eventHeader.StartTime.GetUTCDaqTime()/10

      flasherMap = frame["flasher"]
      for flash in flasherMap:
        flashTime = flash.GetFlashTime() + frameStartTime
        omKey =  flash.GetFlashingOM()
        rate = flash.GetRate()
        if omKey not in self.omList.keys():
          self.omList[omKey] = digitalOpticalModule(omKey, rate, flashTime)
        else:
          self.omList[omKey].addFlash(flashTime)

  def printModuleStats(self):
    tenHzList = []
    threeHzList = []
    
    for module in self.omList.values():
      module.printStatistics()
    for module in self.omList.values():
      if module.rate == 10:
        tenHzList.append( module )
      elif module.rate == 3:
        threeHzList.append( module )
      else:
        print "Module found with non-standard rate (not 3 or 10)"
        
    print "General DOM Statistics"
    if len( tenHzList):
      print "Average DOM Period (10 hz): ", sum( [module.getMeanPeriod() for module in tenHzList] ) / len(tenHzList)
      print "STD DOM Period (10 HZ):", math.sqrt( sum( [module.getMeanPeriod()**2 for module in tenHzList] ) / len(tenHzList) \
                                                                              - (sum([module.getMeanPeriod() for module in tenHzList] ) / len(tenHzList) )**2 )
      print "Average Period ST Deviation (10 hz):", sum( [module.getPeriodStandardDeviation() for module in tenHzList] ) / len(tenHzList)
    if len(threeHzList):
      print "Average DOM Period (3 hz): ", sum( [module.getMeanPeriod() for module in threeHzList] ) / len(threeHzList)
      print "STD DOM Period (10 HZ):", math.sqrt( sum( [module.getMeanPeriod()**2 for module in threeHzList] ) / len(threeHzList) \
                                                                              - (sum([module.getMeanPeriod() for module in threeHzList] ) / len(threeHzList) )**2 )
      print "Average Period ST Deviation (3 hz):", sum( [module.getPeriodStandardDeviation() for module in threeHzList] ) / len(threeHzList)

if __name__ == "__main__":

  if len( sys.argv ) < 2:
    print "Usage: %s I3file" % sys.argv[0]
    sys.exit( 0 )

  flashes = flashReader(sys.argv[1:])
  flashes.printModuleStats()
  
  sys.exit( 0 )
