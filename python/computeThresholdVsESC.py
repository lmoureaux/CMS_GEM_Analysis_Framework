#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
    Created on Monday April 24
    
    @author: Brian L. Dorney, Louis Moureaux
    """

#Imports
import csv
import sys, os
import numpy as np
import root_numpy as rp
#from scipy import interpolate

#Python Toolkit Imports
from AnalysisSuiteEfficiencyPredictor import *
from AnalysisSuiteGainMap import *
from AnalysisSuiteClusterCharge import *
from AnalysisSuiteVFATNoise import *
from Utilities import *

#ROOT Imports
from ROOT import gROOT, TArrayD, TF1, TFile, TGraph, TGraph2D, TGraphErrors, TH2F, TRandom2
import ROOT as r

class AnalysisSuiteThresholdVsESC(AnalysisSuiteEfficiencyPredictor):
    #Here inputfilename should be a tab delimited file in the form of:
    #   field_name  value
    #
    #The field_name is not case sensitive but the value is case sensitive (as it may be a filepaths!)
    def __init__(self, inputfilename="config/configEffPredictor.cfg", debug=False):
        # Unfortunately we have to do a partial initialization of AnalysisSuiteEfficiencyPredictor.
        # It's not modular enough but we still want to use some of its methods.

        #Set the debug flag
        self.DEBUG  = debug

        #HV Points To Perform Calculations for
        self.LIST_HVPTS = []

        #Detector Names
        self.NAME_DET_DUT               = "Detector"
        self.NAME_DET_CLUSTQ            = "GE11-III-FIT-0001"
        self.NAME_DET_CLUSTSIZE         = "GE11-IV-CERN-0001"

        #Declare Detector Container
        self.PARAMS_DET_DUT             = PARAMS_DET()

        self.SECTOR_IETA_CLUSTSIZENORM  = 5
        self.SECTOR_IETA_QC5            = 4

        self.SECTOR_IPHI_CLUSTSIZENORM  = 2
        self.SECTOR_IPHI_QC5            = 2

        self.ANA_UNI_GRANULARITY	= 32

        #Declare Gain Containers
        self.PARAMS_GAIN_DET_DUT        = PARAMS_GAIN()
        self.PARAMS_GAIN_DET_CLUSTQ     = PARAMS_GAIN()
        self.PARAMS_GAIN_DET_CLUSTSIZE  = PARAMS_GAIN()

        #File Names
        self.FILE_CLUSTQ_MEAN           = "" #File having ordered triplet of (V_Drift, ClustSize, Mean) data
        self.FILE_CLUSTQ_MPV            = "" #As above but for MPV
        self.FILE_CLUSTQ_SIGMA          = "" #As above but for Sigma

        self.FILE_MAPPING_DUT           = "" #Mapping config file for detector under test
        self.FILE_MAPPING_DUT_VFAT2DET	= "" #Mapping config file for detectur under test which translates VFAT position to (ieta,iphi) coordinate

        self.FILE_MIP_AVG_CLUST_SIZE    = "" #File having a TF1 which is a fit of MIP <ClustSize> vs. Gain

        self.FILE_OUTPUT                = "" #Output filename to be created

        self.FILE_QC5_RESP_UNI          = "" #Framework output file
        self.FILE_SCURVE_DATA		= "" #TFile containing the scurveFitTree TTree with SCurve data

        #TObject Names
        self.TOBJ_NAME_FUNC_AVGCLUSTSIZE= ""

        #Load the input file and set all variables
        list_strLines = []
        with open(inputfilename) as inputFile:
            list_strLines = inputFile.readlines()

        inputFile.close()

        #strip new line character ('\n') from the file
        list_strLines = [x.strip() for x in list_strLines]

        if self.DEBUG:
            print "Field_Name\tValue"

        #Set Variables
        for iPos in range(0,len(list_strLines)):
            #Get this (field_name,value) pair
            strLine = list_strLines[iPos].split("\t")

            #Skip if commented (e.g. first character of first member is "#")
            if strLine[0][0] == "#":
                continue

            #Transform to upper case
            strLine[0] = strLine[0].upper()

            #Print (field_name, value) pair
            if self.DEBUG:
                print "{0}\t{1}".format(strLine[0],strLine[1])

            #Set values
            if strLine[0] == "FILE_FRAMEWORK_OUTPUT":
                #self.ANASUITEGAIN.openInputFile(strLine[1])
                self.FILE_QC5_RESP_UNI = strLine[1]
            elif strLine[0] == "FILE_CLUSTQ_MEAN":
                self.FILE_CLUSTQ_MEAN = strLine[1]
            elif strLine[0] == "FILE_CLUSTQ_MPV":
                self.FILE_CLUSTQ_MPV = strLine[1]
            elif strLine[0] == "FILE_CLUSTQ_SIGMA":
                self.FILE_CLUSTQ_SIGMA = strLine[1]
            elif strLine[0] == "FILE_CLUSTSIZE":
                self.FILE_MIP_AVG_CLUST_SIZE = strLine[1]
            elif strLine[0] == "FILE_OUTPUT":
                #self.FILE_OUTPUT = strLine[1]
                self.FILE_OUTPUT= TFile(strLine[1],"RECREATE","",1)
                self.DIR_CHARGE	= self.FILE_OUTPUT.mkdir("ChargeData")
                self.DIR_MAP_EFF= self.FILE_OUTPUT.mkdir("EffMaps")
                #self.DIR_MAP_GAIN=self.FILE_OUTPUT.mkdir("GainMaps")
            elif strLine[0] == "FILE_DUT_MAPPING_GEO":
                self.FILE_MAPPING_DUT = strLine[1]
            elif strLine[0] == "FILE_DUT_MAPPING_VFATPOS2IETAIPHI":
                self.FILE_MAPPING_DUT_VFAT2DET = strLine[1]
            elif strLine[0] == "FILE_DUT_SCURVEDATA":
                self.FILE_SCURVE_DATA = strLine[1]
            elif strLine[0] == "DUT_GAIN_P0":
                self.PARAMS_GAIN_DET_DUT.GAIN_CURVE_P0 = float(strLine[1])
            elif strLine[0] == "DUT_GAIN_P0_ERR":
                self.PARAMS_GAIN_DET_DUT.GAIN_CURVE_P0_ERR = float(strLine[1])
            elif strLine[0] == "DUT_GAIN_P1":
                self.PARAMS_GAIN_DET_DUT.GAIN_CURVE_P1 = float(strLine[1])
            elif strLine[0] == "DUT_GAIN_P1_ERR":
                self.PARAMS_GAIN_DET_DUT.GAIN_CURVE_P1_ERR = float(strLine[1])
            elif strLine[0] == "DUT_IETA_CLUST_SIZE_NORM":
                self.SECTOR_IETA_CLUSTSIZENORM = int(strLine[1])
            elif strLine[0] == "DUT_IETA_QC5_GAIN_CAL":
                self.SECTOR_IETA_QC5 = int(strLine[1])
            elif strLine[0] == "DUT_IPHI_CLUST_SIZE_NORM":
                self.SECTOR_IPHI_CLUSTSIZENORM = int(strLine[1])
            elif strLine[0] == "DUT_IPHI_QC5_GAIN_CAL":
                self.SECTOR_IPHI_QC5 = int(strLine[1])
            elif strLine[0] == "DUT_NUM_SIM_PTS_PER_RO":
                self.NUMSIMPTS_PER_RO = int(strLine[1])
            elif strLine[0] == "DUT_QC5_RESP_UNI_HVPT":
                self.PARAMS_DET_DUT.DET_IMON_QC5_RESP_UNI = float(strLine[1])
            elif strLine[0] == "DUT_SERIAL_NUMBER":
                self.NAME_DET_DUT = strLine[1]
            elif strLine[0] == "DUT_SLICE_GRANULARITY":
                self.ANA_UNI_GRANULARITY = int(strLine[1])
            elif strLine[0] == "DET_CLUSTQ_SERIAL_NUMBER":
                self.NAME_DET_CLUSTQ = strLine[1]
            elif strLine[0] == "DET_CLUSTQ_GAIN_P0":
                self.PARAMS_GAIN_DET_CLUSTQ.GAIN_CURVE_P0 = float(strLine[1])
            elif strLine[0] == "DET_CLUSTQ_GAIN_P0_ERR":
                self.PARAMS_GAIN_DET_CLUSTQ.GAIN_CURVE_P0_ERR = float(strLine[1])
            elif strLine[0] == "DET_CLUSTQ_GAIN_P1":
                self.PARAMS_GAIN_DET_CLUSTQ.GAIN_CURVE_P1 = float(strLine[1])
            elif strLine[0] == "DET_CLUSTQ_GAIN_P1_ERR":
                self.PARAMS_GAIN_DET_CLUSTQ.GAIN_CURVE_P1_ERR = float(strLine[1])
            elif strLine[0] == "DET_CLUSTSIZE_SERIAL_NUMBER":
                self.NAME_DET_CLUSTSIZE = strLine[1]
            elif strLine[0] == "DET_CLUSTSIZE_TF1_TNAME":
                self.TOBJ_NAME_FUNC_AVGCLUSTSIZE = strLine[1]
            elif strLine[0] == "EFF_HVPT_LIST":
                self.LIST_HVPTS = strLine[1].split(",")
            else:
                print "Input Field Name:"
                print strLine[0]
                print "Not recognized, please cross-check input file:"
                print inputfilename

        del list_strLines

        #Declare Cluster Charge Analysis Suite
        self.ANASUITECLUSTQ = AnalysisSuiteClusterCharge(file_out=self.FILE_OUTPUT,
                                                         params_gain=self.PARAMS_GAIN_DET_CLUSTQ,
                                                         calcGain=True,
                                                         debug=self.DEBUG)

        #Load the Landau Cluster Charge Data
        if self.DEBUG:
            print "Loading Landau Cluster Charge Data"

        if len(self.FILE_CLUSTQ_MEAN) == 0:
            print "Landau Cluster Charge Mean Filename not found, problem!"
        if len(self.FILE_CLUSTQ_MPV) == 0:
            print "Landau Cluster Charge MPV Filename not found, problem!"
        if len(self.FILE_CLUSTQ_SIGMA) == 0:
            print "Landau Cluster Charge Sigma Filename not found, problem!"

        self.ANASUITECLUSTQ.loadData(inputfilename_MEAN=self.FILE_CLUSTQ_MEAN,
                                     inputfilename_MPV=self.FILE_CLUSTQ_MPV,
                                     inputfilename_SIGMA=self.FILE_CLUSTQ_SIGMA)

        #Interpolate the Landau Cluster Charge Data
        if self.DEBUG:
            print "Interpolating Landau Cluster Charge Data"

        self.ANASUITECLUSTQ.interpolateData()

        #Declare Gain Map Analysis Suite
        self.PARAMS_DET_DUT.DETPOS_IETA = self.SECTOR_IETA_QC5
        self.PARAMS_DET_DUT.DETPOS_IPHI = self.SECTOR_IPHI_QC5

        if self.DEBUG:
            print "Loading Mapping"

        self.PARAMS_DET_DUT.loadMapping(self.FILE_MAPPING_DUT, self.DEBUG)

        self.ANASUITEGAIN   = AnalysisSuiteGainMap(file_out=self.FILE_OUTPUT,
                                                   inputfilename=self.FILE_QC5_RESP_UNI,
                                                   #outputfilename=self.FILE_OUTPUT,
                                                   #outputfileoption="RECREATE",
                                                   params_gain=self.PARAMS_GAIN_DET_DUT,
                                                   params_det=self.PARAMS_DET_DUT,
                                                   debug=self.DEBUG)

        #Load the MIP average cluster size vs. gain formulat
        fileMIPAvgClustSize = TFile(self.FILE_MIP_AVG_CLUST_SIZE, "READ","", 1)

        if self.DEBUG:
            print "Loading function for MIP <Clust Size> vs. Gain"

        self.FUNC_AVG_MIP_CLUSTSIZE = fileMIPAvgClustSize.Get(self.TOBJ_NAME_FUNC_AVGCLUSTSIZE)

        return

    def run(self, fHVPt, inputFile, outprefix, zTrim=0.):
        #Calculate the original gain map of the detector, needed to determine gain at fHVOrGain
        self.calcGainMap()

        #Get the gain map at fHVPt
        if self.DEBUG:
            print "Calculating Gain Map at Input HV Pt: {0}".format(fHVPt)
        
        g2D_Map_Gain_HVPt = self.ANASUITEGAIN.calcGainMapHV(self.NAME_DET_DUT, fHVPt)

        #Calculate the original normalized cluster map of the detector
        self.calcNormAvgClustSizeMap()

        #Now interesting problem N_pts in g2D_Map_Gain_HVPt <= N_pts in self.ANASUITEGAIN.G2D_MAP_AVG_CLUST_SIZE_NORM
        #We can make some numpy arrays
        #Then the (clustPos,y) members from G2D_MAP_AVG_CLUST_SIZE_NORM that don't appear in g2D_Map_Gain_HVPt
        
        #Get the gain and normalized average clustersize data
        data_Gain        = getDataTGraph2D(g2D_Map_Gain_HVPt)
        data_NormAvgCS   = getDataTGraph2D(self.ANASUITEGAIN.G2D_MAP_AVG_CLUST_SIZE_NORM)

        #Close TFiles that have been opened by self.ANASUITEGAIN
	
        #Create a dictionary where the coordinate point (x,y) is mapped to the (gain, Norm <CS>)
        dict_Coords_GainAndNormAvgCS = {(idx[0],idx[1]):[idx[2]] for idx in data_Gain}
        for idx in data_NormAvgCS:
            if (idx[0],idx[1]) in dict_Coords_GainAndNormAvgCS:
            	dict_Coords_GainAndNormAvgCS[(idx[0],idx[1])].append(idx[2])

        # Read (ieta,iphi) to vfatN map
        iEtaiPhiToVfatN = np.zeros((8, 3), dtype=int)
        with open(self.FILE_MAPPING_DUT_VFAT2DET) as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel-tab')
            reader.next() # Skip header
            for row in reader:
                iEtaiPhiToVfatN[int(row['ieta'])-1][int(row['iphi'])-1] = int(row['vfatPos'])

        landauMPV = np.zeros((24, 128))

        #Loop over all points of the detector
        for coordPt in dict_Coords_GainAndNormAvgCS:
            #Get Coordinate Info
            coordPt_iEtaiPhi	= self.PARAMS_DET_DUT.getiEtaiPhIndex(coordPt[0],coordPt[1])
            coordPt_iStripRange	= self.PARAMS_DET_DUT.LIST_DET_GEO_PARAMS[coordPt_iEtaiPhi[0]-1].getStripRange(coordPt[0], self.ANA_UNI_GRANULARITY)
            coordPt_iThisSlice	= self.PARAMS_DET_DUT.LIST_DET_GEO_PARAMS[coordPt_iEtaiPhi[0]-1].getSliceIdx(coordPt[0], self.ANA_UNI_GRANULARITY)

            #Get the gain and MIP cluster size
            fGain 	= dict_Coords_GainAndNormAvgCS[coordPt][0]
            fNormAvgCS	= dict_Coords_GainAndNormAvgCS[coordPt][1]
            fMIPAvgCS	= fNormAvgCS * self.FUNC_AVG_MIP_CLUSTSIZE.Eval(fGain)
            
            #Get the Landau Parameters
            fLandauMean	= self.ANASUITECLUSTQ.getInterpolatedMean(fMIPAvgCS, fGain)
            fLandauMPV	= self.ANASUITECLUSTQ.getInterpolatedMPV(fMIPAvgCS, fGain)
            fLandauSigma= self.ANASUITECLUSTQ.getInterpolatedSigma(fMIPAvgCS, fGain)

            #Save Landau MPV for relevant channels
            vfat = iEtaiPhiToVfatN[coordPt_iEtaiPhi[0]-1][coordPt_iEtaiPhi[1]-1]
            for chan in range(coordPt_iStripRange[0], coordPt_iStripRange[1]):
                landauMPV[vfat][chan] = fLandauMPV
                pass
            pass

        # DAC units -> fC conversion
        vToQb = -0.8
        vToQm = 0.05

        detName = self.NAME_DET_DUT

        notHotHistogram = r.TH2D('%s_notHot'%detName, '%s (not masked);Expected Signal Charge MPV [fC];mu-ztrim*sigma [fC]'%detName,64,4,8,64,vToQm*-0.5+vToQb,vToQm*255.5+vToQb)
        hotHistogram = r.TH2D('%s_hot'%detName, '%s (hot);Expected Signal Charge MPV [fC];mu-ztrim*sigma [fC]'%detName,64,4,8,64,vToQm*-0.5+vToQb,vToQm*255.5+vToQb)

        # Read S-curve fit tree and fill histograms
        inF = r.TFile(inputFile, 'READ')
        for event in inF.scurveFitTree:
            vfat = event.vfatN
            chan = event.vfatCH
            mpv = landauMPV[vfat][chan]
            opThreshold = event.threshold - zTrim * event.noise
            if not event.mask:
                notHotHistogram.Fill(mpv, vToQm*opThreshold+vToQb)
            elif event.maskReason == 0x01:
                hotHistogram.Fill(mpv, vToQm*opThreshold+vToQb)
            pass

        # Write
        outF = r.TFile(outprefix + '.root', 'RECREATE')
        notHotHistogram.Write()
        hotHistogram.Write()
        outF.Close()

        canvas = r.TCanvas('canv', 'canv', 500, 500)
        r.gStyle.SetOptStat(0)

        notHotHistogram.Draw('colz')
        canvas.Update()
        canvas.SaveAs(outprefix + '-notHot.png')

        hotHistogram.Draw('colz')
        canvas.Update()
        canvas.SaveAs(outprefix + '-hot.png')

        hotHistogram.SetTitle('%s (hot and not hot channels);Expected Signal Charge MPV [fC];mu-ztrim*sigma [fC]'%detName)
        hotHistogram.Add(notHotHistogram)
        hotHistogram.Draw('colz')
        canvas.Update()
        canvas.SaveAs(outprefix + '-all.png')

if __name__ == "__main__":
    r.gROOT.SetBatch(True)

    from optparse import OptionParser
    parser = OptionParser()

    parser.add_option("-d", "--debug", action="store_true", dest="debug",
                      help="print extra debugging information", metavar="debug")
    parser.add_option("-c", "--config", type="string", dest="config",
                      help="Specify Configuration Filename", metavar="config.cfg")
    parser.add_option("-i", "--infilename", type="string", dest="infilename",
                      help="Specify Input Filename", metavar="infilename")
    parser.add_option("-o", "--outprefix", type="string", dest="outprefix",
                      help="Specify Output Prefix", metavar="outprefix",
                      default="ESCAnalysis")
    parser.add_option("--ztrim", type="float", dest="ztrim", default=0.0,
                      help="Specify the p value of the trim", metavar="ztrim")
    parser.add_option("--hvpt", type="float", dest="hvpt", default=655.,
                      help="Specify HV point", metavar="hvpt")
    (options, args) = parser.parse_args()

    predictor = AnalysisSuiteThresholdVsESC(options.config, debug=options.debug)
    predictor.run(options.hvpt,
                  options.infilename,
                  options.outprefix,
                  zTrim = options.ztrim)
