//
//  SelectorHit.h
//  
//
//  Created by Brian L Dorney on 11/03/16.
//
//

#ifndef ____SelectorHit__
#define ____SelectorHit__

//C++ Includes
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

//Framework Includes
#include "DetectorMPGD.h"
#include "Selector.h"
#include "UniformityUtilityTypes.h"
#include "TimingUtilityFunctions.h"

//ROOT Includes
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<short>+;
//#endif

namespace QualityControl {
    namespace Uniformity {
        class SelectorHit : public Selector {
            
        public:
            //Constructors
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Default
            SelectorHit();
            
            //Actions - Methods that Do Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Filters an input vector<Uniformity::Cluster object based on the stored Uniformity::AnalysisSetupUniformity attribute
            virtual std::vector<Uniformity::Hit> filterHits(std::vector<Uniformity::Hit> vec_inputHits);
            
            //Filters an input vector<Uniformity::Cluster> object based on the provided Uniformity::AnalysisSetupUniformity object
            virtual std::vector<Uniformity::Hit> filterHits(std::vector<Uniformity::Hit> vec_inputHits, Uniformity::AnalysisSetupUniformity inputSetup){
                setAnalysisParameters(inputSetup);
                return filterHits(vec_inputHits);
            };
            
            //Check if Hit Passes selection stored in aSetupUniformity? True -> Passes; False -> Fails
            bool hitPassesSelection(Uniformity::Hit &inputHit);
            
            //Getters - Methods that Get (i.e. Return) Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Printers - Methods that Print Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Setters - Methods that Set Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            //Given an output ROOT file from amoreSRS with hits
            //Applies the hit selection and stores those selected hits in inputDet
            //Input is a std::string storing the physical filename
            virtual void setHits(std::string &strInputRootFileName, Uniformity::DetectorMPGD &inputDet);
            
            //Given an output ROOT file from amoreSRS with hits
            //Applies the hit selection and stores those selected hits in inputDet
            //Input is a TFile *
            virtual void setHits(TFile * file_InputRootFile, Uniformity::DetectorMPGD &inputDet);
            
            
            //As above but overwrites the stored AnalysisSetupUniformity object
            //Input is a std::string storing the physical filename
            virtual void setHits(std::string &strInputRootFileName, Uniformity::DetectorMPGD &inputDet, Uniformity::AnalysisSetupUniformity inputSetup){
                setAnalysisParameters(inputSetup);
                setHits(strInputRootFileName, inputDet);
                return;
            };
            
            //As above but overwrites the stored AnalysisSetupUniformity object
            //Input is a TFile *
            virtual void setHits(TFile * file_InputRootFile, Uniformity::DetectorMPGD &inputDet, Uniformity::AnalysisSetupUniformity inputSetup){
                setAnalysisParameters(inputSetup);
                setHits(file_InputRootFile, inputDet);
                return;
            };
            
        private:
            //Actions - Methods that Do Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Getters - Methods that Get (i.e. Return) Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Printers - Methods that Print Something
            //------------------------------------------------------------------------------------------------------------------------------------------
            
            //Setters - Methods that Set Something
            //------------------------------------------------------------------------------------------------------------------------------------------
        }; //End class SelectorHit
    } //End namespace Uniformity
} //End namespace QualityControl

#endif /* defined(____SelectorHit__) */
