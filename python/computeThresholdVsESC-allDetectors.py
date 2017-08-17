#! /usr/bin/env python
'''
    Run computeThresholdVsESC on all detectors defined in a config file

    Created: 2017-08-17
    Author: Louis Moureaux
'''

import json
import os
from collections import namedtuple
from optparse import OptionParser

import ROOT as r
from gempython.utils.wrappers import runCommand

r.gROOT.SetBatch(True)

parser = OptionParser()
parser.add_option("-d", "--debug", action="store_true", dest="debug",
                    help="print extra debugging information", metavar="debug")
parser.add_option("-c", "--config", type="string", dest="config",
                    help="Specify configuration directory", metavar="path")
parser.add_option("--gemdata", type="string", dest="gemdata",
                    help="Specify location of the gemdata directory structure",
                    metavar="/gemdata", default="/gemdata")
parser.add_option("-t", "--time", type="string", dest="time",
                    help="Specify scan date and time", metavar="YYYY.MM.DD.HH.MM")
parser.add_option("--ztrim", type="float", dest="ztrim", default=0.0,
                    help="Specify the p value of the trim", metavar="ztrim")
(options, args) = parser.parse_args()

print 'Loading config files from %s' % options.config
configFile = open(options.config + '/detectors.json')

# See 'How to convert JSON data into a Python object'
# https://stackoverflow.com/questions/6578986/how-to-convert-json-data-into-a-python-object
config = json.load(configFile,
                   object_hook=lambda d: namedtuple('X', d.keys())(*d.values()))

# <gemdata>/<GEMINI>/scurve/<date>/SCurveData/SCurveFitData.root
inputBasePath = '%s/%%s/scurve/%s/SCurveData/SCurveFitData.root' % (
        options.gemdata, options.time)
print 'Input files: ' + inputBasePath % '<location>'

# <config>/<config.cfg>
configBasePath = '%s/%%s' % options.config
print 'Config files: ' + configBasePath % '<config.cfg>'

# Script base directory, will be used to locate computeThresholdVsESC.py
pythonDir = os.path.dirname(os.path.abspath(__file__))

haddCommand = ['hadd', '-f', 'allDetectorsCombined.root']

for det in config.detectors:
    print '### ' + det.serialNumber
    command = ['python2.7', '%s/computeThresholdVsESC.py' % pythonDir]
    if options.debug:
        command.append('--debug')
    command.append('--config=%s' % (configBasePath % det.config))
    command.append('--infilename=%s' % (inputBasePath % det.location))
    command.append('--ztrim=%f' % options.ztrim)
    command.append('--hvpt=%f' % det.hvpt)
    command.append('--detNameOverride=%s' % det.serialNumber)
    command.append('--outprefix=%s' % det.serialNumber.replace('/', ''))
    runCommand(command)
    haddCommand.append('%s.root' % det.serialNumber.replace('/', ''))

print '### Summing...'
runCommand(haddCommand)

print '### Extracting combined plots...'

def drawAndSave(histogram, filename, title):
    canvas = r.TCanvas('canv', 'canv', 500, 500)
    r.gStyle.SetOptStat(0)
    histogram.SetTitle(title)
    histogram.Draw('colz')
    canvas.Update()
    canvas.SaveAs(filename)

sumFile = r.TFile('allDetectorsCombined.root')
drawAndSave(sumFile.hot, 'allDetectorsCombined-hot.png', 'All Detectors Combined (hot)')
drawAndSave(sumFile.notHot, 'allDetectorsCombined-notHot.png', 'All Detectors Combined (not hot)')
allHits = sumFile.hot
allHits.Add(sumFile.notHot)
drawAndSave(allHits, 'allDetectorsCombined-all.png', 'All Detectors Combined (hot and not hot)')
