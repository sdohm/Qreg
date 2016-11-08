#!/usr/bin/python



#
# qreg2vmd.cpp
#
#  Created on: Nov 11,2015
#      Author: Sebastian Dohm <sebastian.dohm@uni-ulm.de>
#



#
# Copyright 2015 Sebastian Dohm
# 
# This file is part of Qreg2VMD.
#
# Qreg2VMD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Qreg2VMD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Qreg2VMD.  If not, see <http://www.gnu.org/licenses/>.
#



import argparse
import getopt
import os.path
import sys

def main(argv):

  inputFileName  = ''
  outputFileName = ''

  maxAtomCount = int(-1)
  currentAtomCount = int(-1)
  linesToGo = int(-1)

  try:
    opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])

  except getopt.GetoptError:
    print 'qreg2vmd.py -i <inputfile> -o <outputfile>'
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h' or opt == '--help':
      print 'qreg2vmd.py -i <inputfile> -o <outputfile>'
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputFileName = arg
    elif opt in ("-o", "--ofile"):
      outputFileName = arg

  if not os.path.isfile(inputFileName):
    print 'Input file ',inputFileName,' does not exist.'
    sys.exit(2)

  if os.path.isfile(outputFileName):
    print 'Cowardly declining to overwrite existing file ', outputFileName, '.'
    sys.exit(2)

  inputFile  = open(inputFileName , 'r')
  outputFile = open(outputFileName, 'w')

  # Get maxAtomCount.
  for line in iter(inputFile):
    nWords = line.split()
    wordCount = len(nWords)
    if wordCount == 1 and int(nWords[0]) > maxAtomCount:
      maxAtomCount = int(nWords[0])
      # print maxAtomCount
  inputFile.seek(0)

  # Add dummy Atoms and change currentAtomCount.
  for line in inputFile:
    linesToGo -= 1
    line = line.replace("O", "C")
    line = line.replace("H", "C")
    # print linesToGo
    nWords = line.split()
    wordCount = len(nWords)
    if wordCount == 1:
      outputFile.write(str(maxAtomCount) + '\n')
      currentAtomCount = int(nWords[0])
      # print currentAtomCount
      dummyCount = maxAtomCount - currentAtomCount
      # print dummyCount
      linesToGo = 2
    else:
      outputFile.write(line) 
    if linesToGo == 0:
      # print line
      currentDummy = line
      # print currentDummy
      while dummyCount != 0:
        outputFile.write(currentDummy)
        dummyCount -= 1
        
  inputFile.close
  outputFile.close

if __name__ == "__main__":
  main(sys.argv[1:])

