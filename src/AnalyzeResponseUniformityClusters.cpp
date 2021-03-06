//
//  AnalyzeResponseUniformityClusters.cpp
//  
//
//  Created by Brian L Dorney on 28/01/16.
//
//

#include "DetectorMPGD.h"
#include "AnalyzeResponseUniformityClusters.h"

using std::cout;
using std::endl;
using std::make_shared;
using std::map;
using std::shared_ptr;
using std::string;
using std::vector;

using QualityControl::Timing::getString;
using QualityControl::Timing::printROOTFileStatus;
using QualityControl::Timing::HistoSetup;
using QualityControl::Timing::stofSafe;

using namespace QualityControl::Uniformity;

//Default Constructor
AnalyzeResponseUniformityClusters::AnalyzeResponseUniformityClusters(){
    strAnalysisName = "analysis";
} //End Default Constructor

//Set inputs at construction
AnalyzeResponseUniformityClusters::AnalyzeResponseUniformityClusters(AnalysisSetupUniformity inputSetup){
    strAnalysisName = "analysis";
    
    //Store Analysis Parameters
    aSetup = inputSetup;
    
    //Store Detector
    //detMPGD = inputDet;
} //End Constructor

//Loops over all stored clusters in an input DetectorMPGD object and fills histograms for the full detector
void AnalyzeResponseUniformityClusters::fillHistos(DetectorMPGD & inputDet){
    //Variable Declaration
    std::multimap<int, Cluster> map_clusters;
    vector<int> vec_iEvtList;
    
    //Determine Cluster Multiplicity - Detector Level
    map_clusters = inputDet.getClusters();
    vec_iEvtList = getVectorOfKeys( map_clusters );
    for (auto iterEvt = vec_iEvtList.begin(); iterEvt != vec_iEvtList.end(); ++iterEvt) {
        inputDet.hMulti_Clust->Fill( map_clusters.count( (*iterEvt) ) );
    }
    
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        
        //Determine Cluster Multiplicity - iEta Level
        map_clusters = inputDet.getClusters( (*iterEta).first );
        vec_iEvtList = getVectorOfKeys( map_clusters );
        for (auto iterEvt = vec_iEvtList.begin(); iterEvt != vec_iEvtList.end(); ++iterEvt) {
            (*iterEta).second.clustHistos.hMulti->Fill( map_clusters.count( (*iterEvt) ) );
        }
        
        //Loop Over Stored iPhi Sectors
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over iPhi Sectors
            
            //Determine Cluster Multiplicity - iPhi Level
            vec_iEvtList = getVectorOfKeys( (*iterPhi).second.map_clusters );
            for (auto iterEvt = vec_iEvtList.begin(); iterEvt != vec_iEvtList.end(); ++iterEvt) {
                (*iterPhi).second.clustHistos.hMulti->Fill( (*iterPhi).second.map_clusters.count( (*iterEvt) ) );
            }
            
            //Loop Over Stored Clusters
            for (auto iterClust = (*iterPhi).second.map_clusters.begin(); iterClust != (*iterPhi).second.map_clusters.end(); ++iterClust) { //Loop Over Stored Clusters
                //Fill iEta Histograms
                (*iterEta).second.clustHistos.hADC->Fill( (*iterClust).second.fADC );
                (*iterEta).second.clustHistos.hPos->Fill( (*iterClust).second.fPos_X );
                (*iterEta).second.clustHistos.hSize->Fill( (*iterClust).second.iSize );
                (*iterEta).second.clustHistos.hTime->Fill( (*iterClust).second.iTimeBin );
                
                if ( (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run.count(iNum_Run) > 0 ) {
                    (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run[iNum_Run]->Fill( (*iterClust).first, (*iterClust).second.fADC );
                    (*iterEta).second.clustHistos.map_hTime_v_EvtNum_by_Run[iNum_Run]->Fill( (*iterClust).first, (*iterClust).second.iTimeBin );
                }
                
                (*iterEta).second.clustHistos.hADC_v_Pos->Fill( (*iterClust).second.fPos_X, (*iterClust).second.fADC );
                (*iterEta).second.clustHistos.hADC_v_Size->Fill( (*iterClust).second.iSize, (*iterClust).second.fADC );
                (*iterEta).second.clustHistos.hADC_v_Time->Fill( (*iterClust).second.iTimeBin, (*iterClust).second.fADC );
                (*iterEta).second.clustHistos.hSize_v_Pos->Fill( (*iterClust).second.fPos_X, (*iterClust).second.iSize );
                
                //Fill iPhi Histograms
                (*iterPhi).second.clustHistos.hADC->Fill( (*iterClust).second.fADC );
                (*iterPhi).second.clustHistos.hSize->Fill( (*iterClust).second.iSize);
                (*iterPhi).second.clustHistos.hTime->Fill( (*iterClust).second.iTimeBin);
                
                if ( (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run.count(iNum_Run) > 0 ) {
                    (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run[iNum_Run]->Fill( (*iterClust).first, (*iterClust).second.fADC );
                    (*iterPhi).second.clustHistos.map_hTime_v_EvtNum_by_Run[iNum_Run]->Fill( (*iterClust).first, (*iterClust).second.iTimeBin );
                }
                
                (*iterPhi).second.clustHistos.hADC_v_Pos->Fill( (*iterClust).second.fPos_X, (*iterClust).second.fADC );
                (*iterPhi).second.clustHistos.hADC_v_Size->Fill( (*iterClust).second.iSize, (*iterClust).second.fADC );
                (*iterPhi).second.clustHistos.hADC_v_Time->Fill( (*iterClust).second.iTimeBin, (*iterClust).second.fADC );
                (*iterPhi).second.clustHistos.hSize_v_Pos->Fill( (*iterClust).second.fPos_X, (*iterClust).second.iSize );
            } //End Loop Over Stored Clusters
        } //End Loop Over iPhi Sectors
    } //End Loop Over iEta Sectors
    
    return;
} //End AnalyzeResponseUniformityClusters::fillHistos() - Full Detector

//Assumes Histos have been filled already (obviously)
void AnalyzeResponseUniformityClusters::fitHistos(DetectorMPGD & inputDet){
    //Variable Declaration
    Double_t *dPeakPos;
    
    float fMin = -1e12, fMax = 1e12;
    float fNormChi2;
    float fPkPos = 0, fPkPosErr = 0;        //Peak Position
    float fPkWidth = 0, fPkWidthErr = 0;    //Peak Width
    
    int iBinMin, iBinMax;	//Bins in histogram encapsulating fMin to fMax
    int iIdxPk, iIdxWidth;	//Position in fit parameter meaning container of the peak and the width parameters
    
    TSpectrum specADC(1,2);    //One peak; 2 sigma away from any other peak
    
    vector<float> vec_fFitRange;
    
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        
        //Loop Over Stored iPhi Sectors
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over iPhi Sectors
            
            //Loop Over Stored Slices
            for (auto iterSlice = (*iterPhi).second.map_slices.begin(); iterSlice != (*iterPhi).second.map_slices.end(); ++iterSlice ) { //Loop Over Slices
                cout<<"=======================================================================\n";
                cout<<"Attempting to Fit (iEta, iPhi, iSlice) = (" << (*iterEta).first << ", " << (*iterPhi).first << ", " << (*iterSlice).first << ")\n";
                
                //Check if the slice histogram does not exist, get it if it doesn't
                if ( (*iterSlice).second.hSlice_ClustADC == nullptr ){
                    if ( (*iterPhi).second.clustHistos.hADC_v_Pos == nullptr ) continue;
                    
                    (*iterSlice).second.hSlice_ClustADC = make_shared<TH1F>( *( (TH1F*) (*iterPhi).second.clustHistos.hADC_v_Pos->ProjectionY( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "Slice" + getString((*iterSlice).first) + "_clustADC").c_str(),(*iterSlice).first,(*iterSlice).first,"") ) );
                }
                
                //Skip this slice if the histogram has zero entries
                if ( !( (*iterSlice).second.hSlice_ClustADC->GetEntries() > 0) ) continue;
                
                //Find peak & store it's position
                specADC.Search( (*iterSlice).second.hSlice_ClustADC.get(), 2, "nobackground", 0.5 );
                dPeakPos = specADC.GetPositionX();
                
                //Initialize Fit
                (*iterSlice).second.fitSlice_ClustADC = make_shared<TF1>( getFit( (*iterEta).first, (*iterPhi).first, (*iterSlice).first, aSetup.fitSetup_clustADC, (*iterSlice).second.hSlice_ClustADC, specADC) );
                
                //Clear the calculated fit range from the previous slice
                vec_fFitRange.clear();
                for (auto iterRange = aSetup.fitSetup_clustADC.m_vec_strFit_Range.begin(); iterRange != aSetup.fitSetup_clustADC.m_vec_strFit_Range.end(); ++iterRange) { //Loop Over Fit Range
                    vec_fFitRange.push_back( getParsedInput( (*iterRange), (*iterSlice).second.hSlice_ClustADC, specADC ) );
                } //End Loop Over Fit Range
                
                //Perform Fit & Store the Result
                TFitResult fitRes_ADC;
                
                TVirtualFitter::SetMaxIterations(10000);
                
                if (vec_fFitRange.size() > 1) { //Case: Fit within the user specific range
                    fMin = (*std::min_element(vec_fFitRange.begin(), vec_fFitRange.end() ) );
                    fMax = (*std::max_element(vec_fFitRange.begin(), vec_fFitRange.end() ) );
                    
                    iBinMin = std::floor( (fMin - aSetup.histoSetup_clustADC.fHisto_xLower) / aSetup.histoSetup_clustADC.fHisto_BinWidth );
                    iBinMax = std::ceil( (fMax - aSetup.histoSetup_clustADC.fHisto_xLower) / aSetup.histoSetup_clustADC.fHisto_BinWidth ) + 1;
                    
                    //Skip this histo if integral over fit range is zero
                    if( !( (*iterSlice).second.hSlice_ClustADC->Integral(iBinMin, iBinMax) > 0 ) ){
                        (*iterSlice).second.fitSlice_ClustADC.reset();
                        continue;
                    }
                    
                    fitRes_ADC = *((*iterSlice).second.hSlice_ClustADC->Fit( (*iterSlice).second.fitSlice_ClustADC.get(),aSetup.fitSetup_clustADC.m_strFit_Option.c_str(),"", fMin, fMax) );
                } //End Case: Fit within the user specific range
                else{ //Case: No range to use
                    fitRes_ADC = *((*iterSlice).second.hSlice_ClustADC->Fit( (*iterSlice).second.fitSlice_ClustADC.get(),aSetup.fitSetup_clustADC.m_strFit_Option.c_str(),"") );
                } //End Case: No range to use
                
                //Determine which point in the TGraphs this is
                int iPoint = std::distance( (*iterPhi).second.map_slices.begin(), iterSlice) + aSetup.iUniformityGranularity * std::distance((*iterEta).second.map_sectorsPhi.begin(), iterPhi);
                
                //Store info from spectrum
                //Store - Number of Peaks (from spectrum)
                (*iterEta).second.gEta_ClustADC_Spec_NumPks->SetPoint(iPoint, (*iterSlice).second.fPos_Center, specADC.GetNPeaks() );
                (*iterEta).second.gEta_ClustADC_Spec_NumPks->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, 0 );
                
                //Store - Peak Position (from spectrum)
                (*iterEta).second.gEta_ClustADC_Spec_PkPos->SetPoint(iPoint, (*iterSlice).second.fPos_Center, dPeakPos[0] );
                (*iterEta).second.gEta_ClustADC_Spec_PkPos->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, 0. );
                
                //Find the index of the PEAK parameter and one of the width parameters (HWHM, FWHM, SIGMA)
                auto iterParamFWHM = std::find(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end(), "FWHM");
                auto iterParamHWHM = std::find(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end(), "HWHM");
                auto iterParamPEAK = std::find(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end(), "PEAK");
                auto iterParamSigma= std::find(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end(), "SIGMA");
                
                //Get the Peak Position
                iIdxPk = std::distance(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), iterParamPEAK);
                fPkPos      = (*iterSlice).second.fitSlice_ClustADC->GetParameter(iIdxPk);
                fPkPosErr   = (*iterSlice).second.fitSlice_ClustADC->GetParError(iIdxPk);
                
                //Get the Peak Width
                if ( iterParamHWHM != aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end() ){ //Case: Fit Parameter List Has Meaning HWHM
                    iIdxWidth= std::distance(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), iterParamHWHM);
                    fPkWidth      = 2. * (*iterSlice).second.fitSlice_ClustADC->GetParameter(iIdxWidth);
                    fPkWidthErr   = 2. * (*iterSlice).second.fitSlice_ClustADC->GetParError(iIdxWidth);
                } //End Case: Fit Parameter List Has Meaning HWHM
                else if( iterParamFWHM != aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end() ){ //Case: Fit Parameter List Has Meaning FWHM
                    iIdxWidth= std::distance(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), iterParamFWHM);
                    fPkWidth      = (*iterSlice).second.fitSlice_ClustADC->GetParameter(iIdxWidth);
                    fPkWidthErr   = (*iterSlice).second.fitSlice_ClustADC->GetParError(iIdxWidth);
                } //End Case: Fit Parameter List Has Meaning FWHM
                else if( iterParamSigma != aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.end() ){ //Case: Fit Parameter List Has Meaning SIGMA
                    iIdxWidth= std::distance(aSetup.fitSetup_clustADC.m_vec_strFit_ParamMeaning.begin(), iterParamSigma);
                    fPkWidth      = 2. * sqrt( 2. * log( 2. ) ) * (*iterSlice).second.fitSlice_ClustADC->GetParameter(iIdxWidth);
                    fPkWidthErr   = 2. * sqrt( 2. * log( 2. ) ) * (*iterSlice).second.fitSlice_ClustADC->GetParError(iIdxWidth);
                } //End Case: Fit Parameter List Has Meaning SIGMA
                
                //Get NormChi2 value
                fNormChi2 = (*iterSlice).second.fitSlice_ClustADC->GetChisquare() / (*iterSlice).second.fitSlice_ClustADC->GetNDF();
                
                //Was the Fit Valid?
                //i.e. did the minimizer succeed in finding the minimm
                (*iterSlice).second.iMinuitStatus  = fitRes_ADC.Status();
                if ( fitRes_ADC.IsValid()
                    && isQualityFit( (*iterSlice).second.fitSlice_ClustADC, iIdxPk )
                    && !std::isinf(fNormChi2)
                    && !std::isnan(fNormChi2) ) { //Case: Valid Fit!!!
                    (*iterPhi).second.fNFitSuccess++;
                    (*iterSlice).second.bFitAccepted = true;
                    
                    //Store Fit parameters - NormChi2
                    (*iterEta).second.gEta_ClustADC_Fit_NormChi2->SetPoint(iPoint, (*iterSlice).second.fPos_Center,  fNormChi2  );
                    (*iterEta).second.gEta_ClustADC_Fit_NormChi2->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, 0. );
                    
                    //Store Fit parameters - Peak Position (from fit)
                    (*iterEta).second.gEta_ClustADC_Fit_PkPos->SetPoint(iPoint, (*iterSlice).second.fPos_Center, fPkPos );
                    (*iterEta).second.gEta_ClustADC_Fit_PkPos->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, fPkPosErr );
                    
                    //Store Fit parameters - Peak Resolution (from fit)
                    (*iterEta).second.gEta_ClustADC_Fit_PkRes->SetPoint(iPoint, (*iterSlice).second.fPos_Center, fPkWidth / fPkPos );
                    (*iterEta).second.gEta_ClustADC_Fit_PkRes->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, sqrt( pow( fPkWidthErr / fPkPos, 2) + pow( ( fPkPosErr * fPkWidth ) / ( fPkPos * fPkPos), 2 ) ) );
                    
                    //Record observables for the summary stat (Used for checking uniformity)
                    inputDet.mset_fClustADC_Fit_NormChi2.insert( fNormChi2 );
                    inputDet.mset_fClustADC_Fit_PkPos.insert( fPkPos );
                    inputDet.mset_fClustADC_Fit_PkRes.insert( fPkWidth / fPkPos );
                } //End Case: Valid Fit!!!
                else{ //Case: Invalid Fit (minimizer did not find minumum)
                    //Store Fit parameters - Peak Position (from fit); when failing
                    (*iterEta).second.gEta_ClustADC_Fit_Failures->SetPoint(iPoint, (*iterSlice).second.fPos_Center, fPkPos );
                    (*iterEta).second.gEta_ClustADC_Fit_Failures->SetPointError(iPoint, 0.5 * (*iterSlice).second.fWidth, fPkPosErr );
                } //End Case: Invalid Fit (minimizer did not find minum)
            } //End Loop Over Slices
        } //End Loop Over iPhi Sectors
    } //End Loop Over iEta Sectors
    
    //Calculate statistics
    if ( inputDet.mset_fClustADC_Fit_NormChi2.size() > 0 ) { //Check if stored fit positions exist
        calcStatistics( inputDet.statClustADC_Fit_NormChi2, inputDet.mset_fClustADC_Fit_NormChi2, "ResponseFitNormChi2" );
    } //End Check if stored fit positions exist
    if ( inputDet.mset_fClustADC_Fit_PkPos.size() > 0 ) { //Check if stored fit positions exist
        calcStatistics( inputDet.statClustADC_Fit_PkPos, inputDet.mset_fClustADC_Fit_PkPos, "ResponseFitPkPos" );
    } //End Check if stored fit positions exist
    if ( inputDet.mset_fClustADC_Fit_PkRes.size() > 0 ) { //Check if stored fit positions exist
        calcStatistics( inputDet.statClustADC_Fit_PkRes, inputDet.mset_fClustADC_Fit_PkRes, "ResponseFitPkRes" );
    } //End Check if stored fit positions exist
    
    return;
} //End AnalyzeResponseUniformityClusters::fitHistos()

//Loops through the detector and initializes all cluster graphs
void AnalyzeResponseUniformityClusters::initGraphsClusters(DetectorMPGD & inputDet){
    //Variable Declaration
    
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        //Initialize Response uniformity graphs - Fit norm Chi2
        (*iterEta).second.gEta_ClustADC_Fit_NormChi2 = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Fit_NormChi2->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Fit_NormChi2" ) ).c_str() );
        
        //Initialize Response uniformity graphs - Fit peak pos
        (*iterEta).second.gEta_ClustADC_Fit_PkPos = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Fit_PkPos->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Fit_PkPos" ) ).c_str() );
        
        //Initialize Response uniformity graphs - Fit peak resolution
        (*iterEta).second.gEta_ClustADC_Fit_PkRes = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Fit_PkRes->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Fit_PkRes" ) ).c_str() );
        
        //Initialize Response uniformity graphs - Positions Were Fit Fails
        (*iterEta).second.gEta_ClustADC_Fit_Failures = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Fit_Failures->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Fit_Failures" ) ).c_str() );
        
        //Initialize Response uniformity graphs - Spec Number of Peaks
        (*iterEta).second.gEta_ClustADC_Spec_NumPks = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Spec_NumPks->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Spec_NumPks" ) ).c_str() );
        
        //Initialize Response uniformity graphs - Spec Peak Pos
        (*iterEta).second.gEta_ClustADC_Spec_PkPos = make_shared<TGraphErrors>( TGraphErrors( aSetup.iUniformityGranularity * (*iterEta).second.map_sectorsPhi.size() ) );
        (*iterEta).second.gEta_ClustADC_Spec_PkPos->SetName( ( getNameByIndex( (*iterEta).first, -1, -1, "g", "clustADC_Spec_PkPos" ) ).c_str() );
    } //End Loop Over iEta Sectors
    
    return;
} //End AnalyzeResponseUniformityClusters::initGraphsClusters()

//Loops through the detector and initalizes all cluster histograms
void AnalyzeResponseUniformityClusters::initHistosClusters(DetectorMPGD & inputDet){
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        
        //Grab Eta Sector width (for clustPos Histo)
        aSetup.histoSetup_clustPos.iHisto_nBins  = (*iterEta).second.map_sectorsPhi.size() * aSetup.iUniformityGranularity;
        aSetup.histoSetup_clustPos.fHisto_xLower = -0.5*(*iterEta).second.fWidth;
        aSetup.histoSetup_clustPos.fHisto_xUpper = 0.5*(*iterEta).second.fWidth;
        
        //Initialize iEta Histograms - 1D
        (*iterEta).second.clustHistos.hADC = make_shared<TH1F>(getHistogram((*iterEta).first, -1, aSetup.histoSetup_clustADC ) );
        (*iterEta).second.clustHistos.hMulti = make_shared<TH1F>(getHistogram((*iterEta).first, -1, aSetup.histoSetup_clustMulti ) );
        (*iterEta).second.clustHistos.hPos = make_shared<TH1F>(getHistogram((*iterEta).first, -1, aSetup.histoSetup_clustPos ) );
        (*iterEta).second.clustHistos.hSize = make_shared<TH1F>(getHistogram((*iterEta).first, -1, aSetup.histoSetup_clustSize ) );
        (*iterEta).second.clustHistos.hTime = make_shared<TH1F>(getHistogram((*iterEta).first, -1, aSetup.histoSetup_clustTime ) );
        
        //Initialize iEta Histograms - 2D
        (*iterEta).second.clustHistos.hADC_v_Pos = make_shared<TH2F>( getHistogram2D((*iterEta).first, -1, aSetup.histoSetup_clustPos, aSetup.histoSetup_clustADC ) );
        (*iterEta).second.clustHistos.hADC_v_Size = make_shared<TH2F>( getHistogram2D((*iterEta).first, -1, aSetup.histoSetup_clustSize, aSetup.histoSetup_clustADC ) );
        (*iterEta).second.clustHistos.hADC_v_Time = make_shared<TH2F>( getHistogram2D((*iterEta).first, -1, aSetup.histoSetup_clustTime, aSetup.histoSetup_clustADC ) );
        (*iterEta).second.clustHistos.hSize_v_Pos = make_shared<TH2F>( getHistogram2D((*iterEta).first, -1, aSetup.histoSetup_clustPos, aSetup.histoSetup_clustSize ) );
        
        //Loop Over Stored iPhi Sectors
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over iPhi Sectors
            //Special case for Cluster Position, the number of bins here should be factor 3 less than requested (eta case)
            Timing::HistoSetup setupClustPosPhi	= aSetup.histoSetup_clustPos;
            setupClustPosPhi.iHisto_nBins		= aSetup.iUniformityGranularity;
            setupClustPosPhi.fHisto_xLower		= (*iterPhi).second.fPos_Xlow;
            setupClustPosPhi.fHisto_xUpper 		= (*iterPhi).second.fPos_Xhigh;
            
            //Initialize iPhi Histograms - 1D
            (*iterPhi).second.clustHistos.hADC = make_shared<TH1F>(getHistogram( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustADC ) );
            (*iterPhi).second.clustHistos.hMulti = make_shared<TH1F>(getHistogram( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustMulti ) );
            (*iterPhi).second.clustHistos.hSize = make_shared<TH1F>(getHistogram( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustSize ) );
            (*iterPhi).second.clustHistos.hTime = make_shared<TH1F>(getHistogram( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustTime ) );
            
            //Initialize iPhi Histograms - 2D
            (*iterPhi).second.clustHistos.hADC_v_Pos = make_shared<TH2F>( getHistogram2D( (*iterEta).first, (*iterPhi).first, setupClustPosPhi, aSetup.histoSetup_clustADC ) );
            (*iterPhi).second.clustHistos.hADC_v_Size = make_shared<TH2F>( getHistogram2D( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustSize, aSetup.histoSetup_clustADC ) );
            (*iterPhi).second.clustHistos.hADC_v_Time = make_shared<TH2F>( getHistogram2D( (*iterEta).first, (*iterPhi).first, aSetup.histoSetup_clustTime, aSetup.histoSetup_clustADC ) );
            (*iterPhi).second.clustHistos.hSize_v_Pos = make_shared<TH2F>( getHistogram2D( (*iterEta).first, (*iterPhi).first, setupClustPosPhi, aSetup.histoSetup_clustSize ) );
            
            //Setup the Slices
            for (int i=1; i<= aSetup.iUniformityGranularity; ++i) { //Loop Over Slices
                //Create the slice
                SectorSlice slice;
                
                //Store position information for this slice
                slice.fPos_Center = (*iterPhi).second.clustHistos.hADC_v_Pos->GetXaxis()->GetBinCenter(i);
                slice.fWidth = (*iterPhi).second.clustHistos.hADC_v_Pos->GetXaxis()->GetBinWidth(i);
                
                //Store the slice
                (*iterPhi).second.map_slices[i] = slice;
            } //End Loop Over Slices
        } //End Loop Over iPhi Sectors
    } //End Loop Over iEta Sectors
    
    //Initialize histograms over the entire detector
    inputDet.hMulti_Clust = make_shared<TH1F>(getHistogram( -1, -1, aSetup.histoSetup_clustMulti ) );
    
    return;
} //End AnalyzeResponseUniformityClusters::initHistosClusters()

//Loops through the detector and initalizes all cluster histograms
void AnalyzeResponseUniformityClusters::initHistosClustersByRun(int iInputRunNo, DetectorMPGD & inputDet){
	//cout<<"AnalyzeResponseUniformityClusters::initHistosClustersByRun() - iInputRunNo = " << iInputRunNo << endl;

    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        //Initialize iEta Histograms - 1D
        
            //Placeholder
        
        //Initialize iEta Histograms - 2D
        (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run[iInputRunNo] = make_shared<TH2F>( TH2F(getNameByIndex( (*iterEta).first, -1, -1, "h", "clustADC_v_EvtNum_Run" + getString(iInputRunNo) ).c_str(), "", 50, 0, 500000, aSetup.histoSetup_clustADC.iHisto_nBins, aSetup.histoSetup_clustADC.fHisto_xLower, aSetup.histoSetup_clustADC.fHisto_xUpper ) );
        (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run[iInputRunNo]->Sumw2();
        
        (*iterEta).second.clustHistos.map_hTime_v_EvtNum_by_Run[iInputRunNo] = make_shared<TH2F>( TH2F(getNameByIndex( (*iterEta).first, -1, -1, "h", "clustTime_v_EvtNum_Run" + getString(iInputRunNo) ).c_str(), "", 50, 0, 500000, aSetup.histoSetup_clustTime.iHisto_nBins, aSetup.histoSetup_clustTime.fHisto_xLower, aSetup.histoSetup_clustTime.fHisto_xUpper ) );
        (*iterEta).second.clustHistos.map_hTime_v_EvtNum_by_Run[iInputRunNo]->Sumw2();
        
        //Loop Over Stored iPhi Sectors
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over iPhi Sectors
            //Initialize iPhi Histograms - 1D
            
                //Placeholder
            
            //Initialize iPhi Histograms - 2D
            (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run[iInputRunNo] = make_shared<TH2F>( TH2F(getNameByIndex( (*iterEta).first, (*iterPhi).first, -1, "h", "clustADC_v_EvtNum_Run" + getString(iInputRunNo) ).c_str(), "", 50, 0, 500000, aSetup.histoSetup_clustADC.iHisto_nBins, aSetup.histoSetup_clustADC.fHisto_xLower, aSetup.histoSetup_clustADC.fHisto_xUpper ) );
            (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run[iInputRunNo]->Sumw2();
            
            (*iterPhi).second.clustHistos.map_hTime_v_EvtNum_by_Run[iInputRunNo] = make_shared<TH2F>( TH2F(getNameByIndex( (*iterEta).first, (*iterPhi).first, -1, "h", "clustTime_v_EvtNum_Run" + getString(iInputRunNo) ).c_str(), "", 50, 0, 500000, aSetup.histoSetup_clustTime.iHisto_nBins, aSetup.histoSetup_clustTime.fHisto_xLower, aSetup.histoSetup_clustTime.fHisto_xUpper ) );
            (*iterPhi).second.clustHistos.map_hTime_v_EvtNum_by_Run[iInputRunNo]->Sumw2();
            
            //Setup the Slices
            /*for (int i=1; i<= aSetup.iUniformityGranularity; ++i) { //Loop Over Slices
                
                //Placeholder
                
            }*/ //End Loop Over Slices
        } //End Loop Over iPhi Sectors
    } //End Loop Over iEta Sectors
    
    //Initialize histograms over the entire detector
    
        //Placeholder
    
    return;
} //End AnalyzeResponseUniformityClusters::initHistosClustersByRun()

//Loads a ROOT file previously created by an instance of AnalyzeResponseUniformityClusters
//Loads all TObjects found in the input ROOT file into detMPGD;
//Any previously stored information in detMPGD is lost.
//Takes a std::string which stores the physical filename as input
void AnalyzeResponseUniformityClusters::loadHistosFromFile( std::string & strInputMappingFileName, std::string & strInputROOTFileName ){
    //Variable Declaration
    TFile *ptr_fileInput = nullptr;
    
    //This method will be called when the user wants to re-run the fitting on a previously created batch of histograms
    //The user will directly supply an AMORE mapping file, this will make the DetectorMPGD structure so that it matches the one created when the histograms where first booked
    //The user will indirectly supply an analysis config file because they want to re-run the fits and need to give new information
    //Use previously existing framework code to setup the detector MPGD, then this method behaves as the reverse of storeHistos()
    
    //Open the requested ROOT file
    //------------------------------------------------------
    ptr_fileInput = new TFile(strInputROOTFileName.c_str(),"READ","",1);
    
    //Check to see if data file opened successfully, if so load the tree
    //------------------------------------------------------
    if ( !ptr_fileInput->IsOpen() || ptr_fileInput->IsZombie() ) { //Case: failed to load ROOT file
        perror( ("Uniformity::AnalyzeResponseUniformityClusters::loadHistosFromFile() - error while opening file: " + strInputROOTFileName ).c_str() );
        Timing::printROOTFileStatus(ptr_fileInput);
        std::cout << "Exiting!!!\n";
        
        return;
    } //End Case: failed to load ROOT file
    
    //Call loadHistosFromFile below
    loadHistosFromFile(strInputMappingFileName, ptr_fileInput);
    
    //Close the file
    //------------------------------------------------------
    ptr_fileInput->Close();
    
    return;
} //End AnalyzeResponseUniformityClusters::loadHistosFromFile()

//Loads a ROOT file previously created by an instance of AnalyzeResponseUniformityClusters
//Loads all TObjects found in the input ROOT file into detMPGD;
//Any previously stored information in detMPGD is lost.
//Takes a TFile * which the histograms are written to as input
void AnalyzeResponseUniformityClusters::loadHistosFromFile(std::string & strInputMappingFileName, TFile * file_InputRootFile ){
    //Variable Declaration
    
    //This method will be called when the user wants to re-run the fitting on a previously created batch of histograms
    //The user will directly supply an AMORE mapping file, this will make the DetectorMPGD structure so that it matches the one created when the histograms where first booked
    //The user will indirectly supply an analysis config file because they want to re-run the fits and need to give new information
    //Use previously existing framework code to setup the detector MPGD, then this method behaves as the reverse of storeHistos()
    
    //Setup the MPGD object
    //------------------------------------------------------
    ParameterLoaderDetector loadDetector;
    loadDetector.loadAmoreMapping(strInputMappingFileName);
    detMPGD = loadDetector.getDetector();

    //Check to see if data file opened successfully, if so load the tree
    //------------------------------------------------------
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie() ) { //Case: failed to load ROOT file
        perror( ("Uniformity::AnalyzeResponseUniformityClusters::loadHistosFromFile() - error while opening file: " + (string) file_InputRootFile->GetName() ).c_str() );
        Timing::printROOTFileStatus(file_InputRootFile);
        std::cout << "Exiting!!!\n";
        
        return;
    } //End Case: failed to load ROOT file
    
    //Loop Through the file and load all stored TObjects
    //------------------------------------------------------
    //Loop Over Stored iEta Sectors
    for (auto iterEta = detMPGD.map_sectorsEta.begin(); iterEta != detMPGD.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        
        //Get Directory
        //-------------------------------------
        //Check to see if the directory exists already
        TDirectory *dir_SectorEta = file_InputRootFile->GetDirectory( ( "SectorEta" + getString( (*iterEta).first ) ).c_str(), false, "GetDirectory" );
        
        //If the above pointer is null the directory does NOT exist, skip this Eta Sector
        if (dir_SectorEta == nullptr) continue;
        
        //Load Histograms - ReadoutSectorEta Level
        //-------------------------------------
        dir_SectorEta->cd();        
        (*iterEta).second.clustHistos.hADC = make_shared<TH1F>( *((TH1F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustADC").c_str() ) ) );
        (*iterEta).second.clustHistos.hMulti = make_shared<TH1F>( *((TH1F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustMulti").c_str() ) ) );
        (*iterEta).second.clustHistos.hPos = make_shared<TH1F>( *((TH1F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustPos").c_str() ) ) );
        (*iterEta).second.clustHistos.hSize = make_shared<TH1F>( *((TH1F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustSize").c_str() ) ) );
        (*iterEta).second.clustHistos.hTime = make_shared<TH1F>( *((TH1F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustTime").c_str() ) ) );
        (*iterEta).second.clustHistos.hADC_v_Pos    = make_shared<TH2F>( *((TH2F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustADC_v_clustPos").c_str() ) ) );
        (*iterEta).second.clustHistos.hADC_v_Size   = make_shared<TH2F>( *((TH2F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustADC_v_clustSize").c_str() ) ) );
        (*iterEta).second.clustHistos.hADC_v_Time   = make_shared<TH2F>( *((TH2F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustADC_v_clustTime").c_str() ) ) );
        (*iterEta).second.clustHistos.hSize_v_Pos   = make_shared<TH2F>( *((TH2F*) dir_SectorEta->Get( ("h_iEta" + getString( (*iterEta).first ) +  "_clustSize_v_clustPos").c_str() ) ) );
        
        //Loop Over Stored iPhi Sectors within this iEta Sector
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over Stored iPhi Sectors
            //Get Directory
            //-------------------------------------
            //Check to see if the directory exists already
            TDirectory *dir_SectorPhi = dir_SectorEta->GetDirectory( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str(), false, "GetDirectory"  );
            
            //If the above pointer is null the directory does NOT exist, skip this Phi Sector
            if (dir_SectorPhi == nullptr) continue;
            
            //Load Histograms - ReadoutSectorPhi Level
            //-------------------------------------
            dir_SectorPhi->cd();
            (*iterPhi).second.clustHistos.hADC = make_shared<TH1F>( *((TH1F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustADC").c_str() ) ) );
            (*iterPhi).second.clustHistos.hMulti = make_shared<TH1F>( *((TH1F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustMulti").c_str() ) ) );	    
            (*iterPhi).second.clustHistos.hSize = make_shared<TH1F>( *((TH1F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustSize").c_str() ) ) );
            (*iterPhi).second.clustHistos.hTime = make_shared<TH1F>( *((TH1F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustTime").c_str() ) ) );
            (*iterPhi).second.clustHistos.hADC_v_Pos    = make_shared<TH2F>( *((TH2F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustADC_v_clustPos").c_str() ) ) );
            (*iterPhi).second.clustHistos.hADC_v_Size   = make_shared<TH2F>( *((TH2F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustADC_v_clustSize").c_str() ) ) );
            (*iterPhi).second.clustHistos.hADC_v_Time   = make_shared<TH2F>( *((TH2F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustADC_v_clustTime").c_str() ) ) );
            (*iterPhi).second.clustHistos.hSize_v_Pos    = make_shared<TH2F>( *((TH2F*) dir_SectorPhi->Get( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "_clustSize_v_clustPos").c_str() ) ) );
            
            //Check to see if 2D histo retrieved successfully
            if ( (*iterPhi).second.clustHistos.hADC_v_Pos == nullptr) continue;
            
            //Set to Global Directory - ReadoutSectorPhi Level
            //-------------------------------------
            //Prevents seg faulting when closing ptr_fileInput
            (*iterPhi).second.clustHistos.hADC_v_Pos->SetDirectory(gROOT);
            
            //Slices
            //Trickery, detMPGD now has eta and phi sectors setup.
            //Slices are NOT setup
            //Here loop from 1 to aSetup.iUniformityGranularity and either:
            //      Option 1: Load slices from the file (slower?)
            //      Option 2: project out from ReadoutSectorPhi::clustHistos::hADC_v_Pos
            //Implemented Option 2; faster many large number of I/O operations??
            for (int i=1; i <= aSetup.iUniformityGranularity; ++i ) { //Loop Over Slices
                //Set Histograms - Slice Level
                //-------------------------------------
                //Creat the slice
                SectorSlice slice;
                
                slice.hSlice_ClustADC = make_shared<TH1F>( *( (TH1F*) (*iterPhi).second.clustHistos.hADC_v_Pos->ProjectionY( ("h_iEta" + getString( (*iterEta).first ) + "iPhi" + getString( (*iterPhi).first ) + "Slice" + getString(i) + "_clustADC").c_str(),i,i,"") ) );
                
                //Make sure to set this histo to the global directory
                slice.hSlice_ClustADC->SetDirectory(gROOT);
                
                //Store position information for this slice
                slice.fPos_Center = (*iterPhi).second.clustHistos.hADC_v_Pos->GetXaxis()->GetBinCenter(i);
                slice.fWidth = (*iterPhi).second.clustHistos.hADC_v_Pos->GetXaxis()->GetBinWidth(i);
                
                //Store the slice
                (*iterPhi).second.map_slices[i] = slice;
            } //End Loop Over Slices
        } //End Loop Over Stored iPhi Sectors
    } //End Loop Over Stored iEta Sectors

    //Load Summary Case Histograms (Special Case for hMulti?)
    //-------------------------------------
    //Check to see if dir_Summary exists, if not create it
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory( "Summary", false, "GetDirectory" );
    if (dir_Summary != nullptr){
		detMPGD.hMulti_Clust = make_shared<TH1F>( *((TH1F*) dir_Summary->Get("h_Summary_clustMulti" ) ) );
	}

    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End AnalyzeResponseUniformityClusters::loadHistosFromFile()

//Stores booked histograms (for those histograms that are non-null)
//Takes a std::string which stores the physical filename as input
void AnalyzeResponseUniformityClusters::storeHistos( string & strOutputROOTFileName, std::string strOption, DetectorMPGD & inputDet){
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the store histos sequence
    storeHistos(ptr_fileOutput, inputDet);
    //storeHistos(ptr_fileOutput);
    
    //Close the ROOT file
    ptr_fileOutput->Close();
    
    return;
} //End storeHistos()

//Stores booked histograms (for those histograms that are non-null)
//Takes a TFile * which the histograms are written to as input
void AnalyzeResponseUniformityClusters::storeHistos( TFile * file_InputRootFile, DetectorMPGD & inputDet){
    //Variable Declaration
    
    //Check if File Failed to Open Correctly
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeHistos", "\tExiting; No Histograms have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Setup the summary histograms
    TH1F hclustADC_All( getHistogram(-1, -1, aSetup.histoSetup_clustADC) );
    TH1F hclustPos_All( getHistogram(-1, -1, aSetup.histoSetup_clustPos) );
    TH1F hclustSize_All( getHistogram(-1, -1, aSetup.histoSetup_clustSize) );
    TH1F hclustTime_All( getHistogram(-1, -1, aSetup.histoSetup_clustTime) );
    
    //Get/Make the Summary Directory
    //Check to see if the directory exists already
    TDirectory *dir_Summary = file_InputRootFile->GetDirectory("Summary", false, "GetDirectory" );
    
    //If the above pointer is null the directory does NOT exist, create it
    if (dir_Summary == nullptr) { //Case: Directory did not exist in file, CREATE
        dir_Summary = file_InputRootFile->mkdir("Summary");
    } //End Case: Directory did not exist in file, CREATE
    
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        //Get Directory
        //-------------------------------------
        //Check to see if dir_SectorEta directory exists already, if not create it
        TDirectory *dir_SectorEta = file_InputRootFile->GetDirectory( ( "SectorEta" + getString( (*iterEta).first ) ).c_str(), false, "GetDirectory" );
        if (dir_SectorEta == nullptr) { //Case: Directory did not exist in file, CREATE
            dir_SectorEta = file_InputRootFile->mkdir( ( "SectorEta" + getString( (*iterEta).first ) ).c_str() );
        } //End Case: Directory did not exist in file, CREATE
        
        //Check to see if dir_RunHistory_Eta directory exists already, if not create it
        TDirectory *dir_RunHistory_Eta = dir_SectorEta->GetDirectory( "RunHistory", false, "GetDirectory" );
        if (dir_RunHistory_Eta == nullptr) {
            dir_RunHistory_Eta = dir_SectorEta->mkdir( "RunHistory" );
        }
        
        //Check to see if dir_RunHistory_Eta_ADC directory exists already, if not create it
        TDirectory *dir_RunHistory_Eta_ADC = dir_RunHistory_Eta->GetDirectory( "ADC", false, "GetDirectory" );
        if (dir_RunHistory_Eta_ADC == nullptr) {
            dir_RunHistory_Eta_ADC = dir_RunHistory_Eta->mkdir( "ADC" );
        }
        
        //Check to see if dir_RunHistory_Eta_Time directory exists already, if not create it
        TDirectory *dir_RunHistory_Eta_Time = dir_RunHistory_Eta->GetDirectory( "Time", false, "GetDirectory" );
        if (dir_RunHistory_Eta_Time == nullptr) {
            dir_RunHistory_Eta_Time = dir_RunHistory_Eta->mkdir( "Time" );
        }
        
        //Add this sector to the summary histogram
        hclustADC_All.Add((*iterEta).second.clustHistos.hADC.get() );
        hclustPos_All.Add((*iterEta).second.clustHistos.hPos.get() );
        hclustSize_All.Add((*iterEta).second.clustHistos.hSize.get() );
        hclustTime_All.Add((*iterEta).second.clustHistos.hTime.get() );
        
        //Store Histograms - ReadoutSectorEta Level
        //-------------------------------------
        dir_SectorEta->cd();
        (*iterEta).second.clustHistos.hADC->Write();
        (*iterEta).second.clustHistos.hMulti->Write();
        (*iterEta).second.clustHistos.hPos->Write();
        (*iterEta).second.clustHistos.hSize->Write();
        (*iterEta).second.clustHistos.hTime->Write();
        (*iterEta).second.clustHistos.hADC_v_Pos->Write();
        (*iterEta).second.clustHistos.hADC_v_Size->Write();
        (*iterEta).second.clustHistos.hADC_v_Time->Write();
        (*iterEta).second.clustHistos.hSize_v_Pos->Write();
        
        auto iterHistoTime = (*iterEta).second.clustHistos.map_hTime_v_EvtNum_by_Run.begin();
        for(auto iterHistoADC = (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run.begin(); iterHistoADC != (*iterEta).second.clustHistos.map_hADC_v_EvtNum_by_Run.end(); ++iterHistoADC){
            dir_RunHistory_Eta_ADC->cd();
            (*iterHistoADC).second->Write();

            dir_RunHistory_Eta_Time->cd();
            (*iterHistoTime).second->Write();
            ++iterHistoTime;
        }
        
        //Loop Over Stored iPhi Sectors within this iEta Sector
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over Stored iPhi Sectors
            //Get Directory
            //-------------------------------------
            //Check to see if dir_SectorPhi directory exists already, if not create it
            TDirectory *dir_SectorPhi = dir_SectorEta->GetDirectory( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str(), false, "GetDirectory"  );
            if (dir_SectorPhi == nullptr) { //Case: Directory did not exist in file, CREATE
                dir_SectorPhi = dir_SectorEta->mkdir( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str() );
            } //End Case: Directory did not exist in file, CREATE
            
            //Check to see if dir_RunHistory_Phi directory exists already, if not create it
            TDirectory *dir_RunHistory_Phi = dir_SectorPhi->GetDirectory( "RunHistory", false, "GetDirectory" );
            if (dir_RunHistory_Phi == nullptr) {
                dir_RunHistory_Phi = dir_SectorPhi->mkdir( "RunHistory" );
            }
            
            //Check to see if dir_RunHistory_Phi_ADC directory exists already, if not create it
            TDirectory *dir_RunHistory_Phi_ADC = dir_RunHistory_Phi->GetDirectory( "ADC", false, "GetDirectory" );
            if (dir_RunHistory_Phi_ADC == nullptr) {
                dir_RunHistory_Phi_ADC = dir_RunHistory_Phi->mkdir( "ADC" );
            }
            
            //Check to see if dir_RunHistory_Phi_Time directory exists already, if not create it
            TDirectory *dir_RunHistory_Phi_Time = dir_RunHistory_Phi->GetDirectory( "Time", false, "GetDirectory" );
            if (dir_RunHistory_Phi_Time == nullptr) {
                dir_RunHistory_Phi_Time = dir_RunHistory_Phi->mkdir( "Time" );
            }
            
            //Store Histograms - ReadoutSectorPhi Level
            //-------------------------------------
            dir_SectorPhi->cd();
            (*iterPhi).second.clustHistos.hADC->Write();
            (*iterPhi).second.clustHistos.hMulti->Write();
            (*iterPhi).second.clustHistos.hSize->Write();
            (*iterPhi).second.clustHistos.hTime->Write();
            (*iterPhi).second.clustHistos.hADC_v_Pos->Write();
            (*iterPhi).second.clustHistos.hADC_v_Size->Write();
            (*iterPhi).second.clustHistos.hADC_v_Time->Write();
            (*iterPhi).second.clustHistos.hSize_v_Pos->Write();
            
            iterHistoTime = (*iterPhi).second.clustHistos.map_hTime_v_EvtNum_by_Run.begin();
            for(auto iterHistoADC = (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run.begin(); iterHistoADC != (*iterPhi).second.clustHistos.map_hADC_v_EvtNum_by_Run.end(); ++iterHistoADC){
                dir_RunHistory_Phi_ADC->cd();
                (*iterHistoADC).second->Write();

                dir_RunHistory_Phi_Time->cd();
                (*iterHistoTime).second->Write();
                ++iterHistoTime;
            }
        } //End Loop Over Stored iPhi Sectors
    } //End Loop Over Stored iEta Sectors
    
    //Store the Summary Histograms
    dir_Summary->cd();
    hclustADC_All.Write();
    inputDet.hMulti_Clust->Write();
    hclustPos_All.Write();
    hclustSize_All.Write();
    hclustTime_All.Write();
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End storeHistos()

//Stores booked cluster fits (for those fits that are non-null)
//Takes a std::string which stores the physical filename as input
void AnalyzeResponseUniformityClusters::storeFits( string & strOutputROOTFileName, std::string strOption, DetectorMPGD & inputDet ){
    //Variable Declaration
    TFile * ptr_fileOutput = new TFile(strOutputROOTFileName.c_str(), strOption.c_str(),"",1);
    
    //Check if File Failed to Open Correctly
    if ( !ptr_fileOutput->IsOpen() || ptr_fileOutput->IsZombie()  ) {
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits","Error: File I/O");
        printROOTFileStatus(ptr_fileOutput);
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits", "\tExiting; No Fits have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Call the store fits sequence
    storeFits(ptr_fileOutput, inputDet);
    //storeFits(ptr_fileOutput);
    
    //Close the ROOT file
    ptr_fileOutput->Close();
    
    return;
} //End storeFits()

//Stores booked cluster fits (for those fits that are non-null)
//Takes a TFile * which the histograms are written to as input
void AnalyzeResponseUniformityClusters::storeFits( TFile * file_InputRootFile, DetectorMPGD & inputDet){
    //TFile does not manage objects
    //TH1::AddDirectory(kFALSE);

    //Variable Declaration
    
    //Check if File Failed to Open Correctly
    if ( !file_InputRootFile->IsOpen() || file_InputRootFile->IsZombie()  ) {
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits","Error: File I/O");
        printROOTFileStatus(file_InputRootFile);
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits", "\tPlease cross check input file name, option, and the execution directory\n" );
        printClassMethodMsg("AnalyzeResponseUniformityClusters","storeFits", "\tExiting; No Fits have been stored!\n" );
        
        return;
    } //End Check if File Failed to Open Correctly
    
    //Loop Over Stored iEta Sectors
    for (auto iterEta = inputDet.map_sectorsEta.begin(); iterEta != inputDet.map_sectorsEta.end(); ++iterEta) { //Loop Over iEta Sectors
        
        //Get Directory
        //-------------------------------------
        //Check to see if dir_SectorEta exists already, if not create it
        TDirectory *dir_SectorEta = file_InputRootFile->GetDirectory( ( "SectorEta" + getString( (*iterEta).first ) ).c_str(), false, "GetDirectory" );
        if (dir_SectorEta == nullptr) { //Case: Directory did not exist in file, CREATE
            dir_SectorEta = file_InputRootFile->mkdir( ( "SectorEta" + getString( (*iterEta).first ) ).c_str() );
        } //End Case: Directory did not exist in file, CREATE
        
        //Store Fits - ReadoutSectorEta Level
        //-------------------------------------
        dir_SectorEta->cd();
        (*iterEta).second.gEta_ClustADC_Fit_NormChi2->Write();
        ( (*iterEta).second.gEta_ClustADC_Fit_PkPos.get() )->Write();
        (*iterEta).second.gEta_ClustADC_Fit_PkRes->Write();
        (*iterEta).second.gEta_ClustADC_Fit_Failures->Write();
        
        (*iterEta).second.gEta_ClustADC_Spec_NumPks->Write();
        (*iterEta).second.gEta_ClustADC_Spec_PkPos->Write();
        
        //Loop Over Stored iPhi Sectors within this iEta Sector
        for (auto iterPhi = (*iterEta).second.map_sectorsPhi.begin(); iterPhi != (*iterEta).second.map_sectorsPhi.end(); ++iterPhi) { //Loop Over Stored iPhi Sectors
            //Get Directory
            //-------------------------------------
            //Check to see if dir_SectorPhi exists already, if not create it
            TDirectory *dir_SectorPhi = dir_SectorEta->GetDirectory( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str(), false, "GetDirectory"  );
            if (dir_SectorPhi == nullptr) { //Case: Directory did not exist in file, CREATE
                dir_SectorPhi = dir_SectorEta->mkdir( ( "SectorPhi" + getString( (*iterPhi).first ) ).c_str() );
            } //End Case: Directory did not exist in file, CREATE
            
            //Store Fits - ReadoutSectorPhi Level
            //-------------------------------------
            dir_SectorPhi->cd();
            
            //No Fits defined at this level - yet
            
            //Slices
            //Now that all clusters have been analyzed we extract the slices
            for (auto iterSlice = (*iterPhi).second.map_slices.begin(); iterSlice != (*iterPhi).second.map_slices.end(); ++iterSlice ) { //Loop Over Slices
                //if ((*iterSlice).second.fitSlice_ClustADC == nullptr) { continue; }
                
                //Get Directory
                //-------------------------------------
                //Check to see if dir_Slice exists already, if not create it
                TDirectory *dir_Slice = dir_SectorPhi->GetDirectory( ( "Slice" + getString( (*iterSlice).first ) ).c_str(), false, "GetDirectory"  );
                if (dir_Slice == nullptr) { //Case: Directory did not exist in file, CREATE
                    dir_Slice = dir_SectorPhi->mkdir( ( "Slice" + getString( (*iterSlice).first ) ).c_str() );
                } //End Case: Directory did not exist in file, CREATE
                
                //Store Fits - Slice Level
                //-------------------------------------
                dir_Slice->cd();
                if ((*iterSlice).second.hSlice_ClustADC != nullptr) (*iterSlice).second.hSlice_ClustADC->Write();
                if ((*iterSlice).second.fitSlice_ClustADC != nullptr) (*iterSlice).second.fitSlice_ClustADC->Write();
            } //End Loop Over Slices
        } //End Loop Over Stored iPhi Sectors
    } //End Loop Over Stored iEta Sectors
    
    //Do not close file_InputRootFile it is used elsewhere
    
    return;
} //End storeFits()
