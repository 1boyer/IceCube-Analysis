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

try:
  from icecube.icetray import *
  from icecube.dataclasses import *
  from icecube.dataio import *
except ImportError:
  print "IceTray is missing! This script needs IceTray!"
  sys.exit(1)


class flashEvent:
  def __init__(self, Time, Position, RecordKeeper, FrameID, Frame, Omkey=None):
    self.time = Time
    self.position = Position # I3 OM Position Class
    self.record = RecordKeeper # Record Keeping Information
    self.frame = Frame # Hold onto frame in memory
    self.omKey = Omkey
    self.frameID = FrameID

  def __repr__(self):
    return str(self.eventType) + " @ " + str( (self.position.X, self.position.Y, self.position.Z, int(self.time)) )

  def __str__(self):
    return self.__repr__()

class recordKeeper:
  def __init__(self, i3filename="No Name"):
    self.startTime = 0
    self.endTime = 0
    self.filename = i3filename
    self.frameCount = 0
    self.flashCount = 0
    self.pairCount = 0
    self.mimeCount = 0
    self.omCount = 0
    self.omRates = {}
    self.omFlashes = {}
    self.mimeList = []

  def updateTime(self, currentTime):
    if self.endTime:
      self.endTime = currentTime
    else:
      self.startTime = currentTime
      self.endTime = currentTime

    self.frameCount = self.frameCount + 1

  def addFlash(self, omKey):
    if omKey in self.omFlashes:
      self.omFlashes[omKey] = self.omFlashes[omKey] + 1
    else:
      self.omFlashes[omKey] = 1

    self.flashCount = self.flashCount + 1

  def addOM(self, omKey, rate):
    if omKey not in self.omRates:
      self.omRates[omKey] = rate
      self.omCount = self.omCount + 1
      self.omFlashes[omKey] = 1
    else:
      self.omFlashes[omKey] = self.omFlashes[omKey] + 1
      
    self.flashCount = self.flashCount + 1


  def addPair(self):
    self.pairCount = self.pairCount + 1

  def addMimic(self, mime):
    self.mimeList.append(mime)
    self.mimeCount = self.mimeCount + 1

  def __str__(self):
    timeElapsed = (self.endTime - self.startTime)/10**9
  
    fillString = "$&$&$&$&$&$&$&$&$&$ [ %s ] $&$&$&$&$&$&$&$&$&$" % self.filename 
    fillString = fillString + "\n%d Frames Mimicing Tau Event, %d Frames with Two Flashers, %d Total Frames [Time Elapsed %d s]\n"  \
        % (self.mimeCount, self.pairCount, self.frameCount, timeElapsed)
    fillString = fillString + "Number of active DOMS: %d (%d flashes)\n" % (self.omCount, self.flashCount)
    for om in self.omRates.keys():
      fillString = fillString + "OM:" + str(om) + " OM Rate:" + str(self.omRates[om]) + " Times Flashed:" + str(self.omFlashes[om]) + "\n"

    fillString = fillString + "\n Mimic Events:\n"
    for mime in self.mimeList:
      fillString = fillString + mime + "\n"

    return fillString

class eventList:
  def __init__(self, initialList = [], tauRadius=100.0, pairRadius = 5775.0):
    self._eventList = initialList
    self._pairWindowRadius = pairRadius
    self._tauWindowRadius = tauRadius
    self._keepers = []
  
  def _pruneOldEvents(self):
    referenceTime = self._eventList[0].time
    
    newList = []
    for event in self._eventList:
      if math.fabs(referenceTime - event.time) <= self._pairWindowRadius:
        newList.append(event)
    
    self._eventList = newList

  def append(self, Time, Position, RecordKeeper, FrameID, Frame, omKey=None):
    self._eventList = [flashEvent(Time, Position, RecordKeeper, FrameID, Frame, omKey)] + self._eventList
    self._pruneOldEvents()
    self._findPairs()
    
  def duration(self, flash1, flash2): # time separation (ns)
    return math.fabs(flash1 - flash2)

  def length(self, pos1, pos2): # Spatial separation (ns)
    return pos1.CalcDistance(pos2)/I3Constants.c

  def _findPairs(self):
    keep = False

    primaryEvent = self._eventList[0]
    ptime = primaryEvent.time
    ppos = primaryEvent.position
    pframeID = primaryEvent.frameID
    pomKey = primaryEvent.omKey
    precord = primaryEvent.record
    

    # Checking all events in current time window for 
    # a double pulse relationship

    for compareEvent in self._eventList[1:]:
      keep = True

      ctime = compareEvent.time
      cpos = compareEvent.position
      cframeID = compareEvent.frameID
      comKey = compareEvent.omKey

      if cframeID is not pframeID:
        print "Warning: Comparing different frames (%s and %s)! Timing errors likely!" % (cframeID, pframeID)
        continue
        
      precord.addPair()
      
      properTimeElapsed = math.fabs(self.duration(ctime,ptime) - self.length(cpos, ppos))
      if properTimeElapsed < self._tauWindowRadius:
        mime = pframeID + "[%2d%2d:%2d%2d:%d ns]" % ( pomKey.GetString(),\
              pomKey.GetOM(), comKey.GetString(), comKey.GetOM(), properTimeElapsed)
        precord.addMimic(mime)
        
    if keep and primaryEvent.frame not in self._keepers:
      self._keepers.append(primaryEvent.frame)

  def __repr__(self):
    return self._keepers

  def __str__(self):
    return str([i for i in self._keepers])

  def __getitem__(self, i):
    return self._keepers[i]
    
### EOC ###



class tauMimicFinder:

  def __init__(self, i3FilenameList=None, Uncertainty=10, OutputFile=None):
    self._positions = {}
    self._events = eventList([], Uncertainty)
    self._output = None
    self._records = []
    self._noGeometry = True

    if OutputFile:
      try:
        self._output = I3File()
        self._output.open_file(OutputFile, I3File.Writing)
      except:
        print "Error opening output file (%s)" % OutputFile
        self._output = None

    if type(i3FilenameList) == types.NoneType:
      i3FilenameList = []
    elif type(i3FilenameList) == types.StringType:
      i3FilenameList = [i3FilenameList]

    for i3Filename in i3FilenameList:
      self.processFile(i3Filename)

  def processFile(self,i3Filename):
    try:
      file = I3File(i3Filename)
      self._record = recordKeeper(i3Filename)
    except:
      (type, val, tback) = sys.exc_info()
      print >> sys.stderr, "tauFinder: Error Opening I3 file:"
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
      self._records.append(self._record)
      if self._output:
        for event in self._events: # Representation of self._events returns
          self._output.push(event)  # the I3 frames marked for keeping
      file.close()


  def frameRead(self, frame):

    if "I3Geometry" in frame and self._noGeometry:
      geometry = frame["I3Geometry"]
      omgeo = geometry.omgeo
      self._noGeometry = False
      for omgeoEntry in omgeo:
        omKey = omgeoEntry.key()
        pos = omgeoEntry.data().position
        self._positions[omKey] = pos
      if self._output:
        self._output.push(frame)

    if "I3EventHeader" and "flasher" in frame:
      eventHeader = frame["I3EventHeader"]
      frameStartTime = eventHeader.StartTime.GetUTCDaqTime()
      frameID = "%d.%d.%d" % (eventHeader.RunID, eventHeader.SubRunID, eventHeader.EventID)
      self._record.updateTime(frameStartTime)

      flasherMap = frame["flasher"]
      for flash in flasherMap:
        flashtime = flash.GetFlashTime()
        omKey =  flash.GetFlashingOM()
        self._record.addOM(omKey, flash.GetRate())
        self._events.append(flashtime + frameStartTime, self._positions[omKey], self._record, frameID, frame, omKey)
    
  def printStatistics(self):
    for record in self._records:
      print record


  ### EOC ###


if __name__ == "__main__":


  def OptparseCallback_windowRadius( option, opt_str, value, parser ):
        # If the value given for the option is negative...
    if value < 0.0:
        raise OptionValueError( "Value %g given for %s option. Value must be "
                                "non-negative"
                                % (value, opt_str) )

    value = round( value, 1 )
        # Store the value that passed inspection
    setattr( parser.values, option.dest, value )


  parser = OptionParser( usage = "%prog [options] i3datafile ..." )
  parser.add_option( "--windowRadius", action = "callback",
                       callback = OptparseCallback_windowRadius, type = "float",
                       dest = "windowRadius", default = 110.0,
                       metavar = "Nanoseconds",
                       help = "How far from a perfect flash separation to "
                              "allow an actual flash separation to count as a "
                              "Tau Double Bang mimic. Units are nanoseconds "
                              "rounded to the nearest tenth of a nanosecond. "
                              "Positive values only. Default is 110.0 ns." )
  parser.add_option("-f", "--file", action="store", type="string", dest="outfile", 
                       default=None, help="Output Tau Mimic Events to a file")

  (options, args) = parser.parse_args()

  if not len( args ):
    print "Usage: %s I3file" % sys.argv[0]
    sys.exit( 0 )

  pairs = tauMimicFinder(args, options.windowRadius, options.outfile)
  pairs.printStatistics()
  
  sys.exit( 0 )
