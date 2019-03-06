# Module to propagate to the MET calculation the re-calibration of the jets:
# MET is Type-1 PF MET and Type-1 corrections are re-calculated following criteria at https://twiki.cern.ch/twiki/bin/view/CMS/METType1Type2Formulae#3_The_Type_I_correction

import ROOT
import math, os,re
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator

class jetRecalib(Module):
    def __init__(self,  globalTag, jetType = "AK4PFchs"):

        if "AK4" in jetType : 
            self.jetBranchName = "Jet"
        elif "AK8" in jetType :
            self.jetBranchName = "FatJet"
            self.subJetBranchName = "SubJet"
        else:
            raise ValueError("ERROR: Invalid jet type = '%s'!" % jetType)
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        self.lenVar = "n" + self.jetBranchName
        # TODO: these lines are useless, please remove
        # To do : change to real values
        self.jmsVals = [1.00, 0.99, 1.01]
        
        self.jesInputFilePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/jme/"

        type1METParams={'jetPtThreshold':15., 'skipEMfractionThreshold':0.9, 'skipMuons':True} 
        #self.jetReCalibrator = JetReCalibrator(globalTag, jetType, True, self.jesInputFilePath, calculateSeparateCorrections = False, calculateType1METCorrection  = False)
        self.jetReCalibrator = JetReCalibrator(globalTag, jetType, True, self.jesInputFilePath, calculateSeparateCorrections=False, calculateType1METCorrection=True, type1METParams=type1METParams)
	
        # load libraries for accessing JES scale factors and uncertainties from txt files
        for library in [ "libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools" ]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '%s'" % library.replace("lib", ""))
                ROOT.gSystem.Load(library)

    def beginJob(self):
	pass

    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("%s_pt_nom" % self.jetBranchName, "F", lenVar=self.lenVar)
        self.out.branch("MET_pt_nom" , "F")
        self.out.branch("MET_phi_nom", "F")
        self.out.branch("MET_pt_naive" , "F")
        self.out.branch("MET_phi_naive", "F")
        self.out.branch("MET_pt_shift" , "F")
        self.out.branch("MET_phi_shift", "F")
            
                        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = Collection(event, self.jetBranchName )
        met = Object(event, "MET") # NOTE: Do not use to re-calibrate 2017 data
        rawmet = Object(event, "RawMET") 
        muons = Collection(event, "Muon")

        # "_nom" collections and branches will store re-calibrated quantities
        jets_pt_nom = []
        (met_px, met_py) = (met.pt*math.cos(met.phi), met.pt*math.sin(met.phi))
        (met_px_nom, met_py_nom) = (met_px, met_py)
        (met_px_naive,met_py_naive) = (met_px, met_py) # a' la stephane
        #met_px_nom = met_px useless line
        #met_py_nom = met_py useless line
                
        rho = getattr(event, self.rhoBranchName)
        if rho==0 or rho==None: print 'WARNING: rho not correctly set for this event'
        #print 'rho={}'.format(rho)

        metShift=[0,0] # met correction accumulator
        type1METCorr=[0,0] # 
        
        for jet in jets:
	    #jet_pt = jet.pt
	    jet_pt = self.jetReCalibrator.correct(jet=jet,muons=muons,rho=rho,delta=0,addCorr=False,addShifts=False, metShift=metShift, type1METCorr=type1METCorr)
            jet_pt_nom = jet_pt # don't smear resolution in data
            if jet_pt_nom < 0.0:
                jet_pt_nom *= -1.0
            jets_pt_nom.append(jet_pt_nom) # concerning jets I am done here ...

            # now comes the met with naive correction
            if jet_pt_nom > 15.:
              jet_cosPhi = math.cos(jet.phi)
              jet_sinPhi = math.sin(jet.phi)
              met_px_naive = met_px_naive - (jet_pt_nom - jet.pt)*jet_cosPhi
              met_py_naive = met_py_naive - (jet_pt_nom - jet.pt)*jet_sinPhi

        # at this point x,y components are already summed in metShift, type1METCorr
        met_px_nom = rawmet.pt*math.cos(rawmet.phi) + type1METCorr[0]
        met_py_nom = rawmet.pt*math.sin(rawmet.phi) + type1METCorr[1]
        met_pt_nom = math.sqrt(met_px_nom**2 + met_py_nom**2)
        met_phi_nom = math.atan2(met_py_nom, met_px_nom) 
        # 
        met_px_shift = met_px + metShift[0] 
        met_py_shift = met_py + metShift[1] 
        met_pt_shift = math.sqrt(met_px_shift**2 + met_py_shift**2)
        met_phi_shift = math.atan2(met_py_shift, met_px_shift) 

        met_pt_naive = math.sqrt(met_px_naive**2 + met_py_naive**2)
        met_phi_naive = math.atan2(met_py_naive,met_px_naive)
        # debug 

        self.out.fillBranch("%s_pt_nom" % self.jetBranchName, jets_pt_nom)
        self.out.fillBranch("MET_pt_nom", met_pt_nom)
        self.out.fillBranch("MET_phi_nom", met_phi_nom)        
        self.out.fillBranch("MET_pt_naive", met_pt_naive)
        self.out.fillBranch("MET_phi_naive", met_phi_naive)        
        self.out.fillBranch("MET_pt_shift", met_pt_shift)
        self.out.fillBranch("MET_phi_shift", met_pt_shift)        

################ my stuff to debug        

        #print 'orig_met_pt={:.1f}, orig_met_phi={:.2f}, nom_met_pt={:.1f} nom_met_phi={:.2f}, shift_met_pt={:.1f}, shift_met_phi={:.1f}'.format(met.pt, met.phi, met_pt_nom, met_phi_nom, met_pt_shift, met_phi_shift) 
        print 'orig_met_pt={:.1f}, orig_met_phi={:.2f}, raw_met_pt={:.1f}, raw_met_phi={:.2f} nom_met_pt={:.1f} nom_met_phi={:.2f}'.format(met.pt, met.phi, rawmet.pt, rawmet.phi, met_pt_nom, met_phi_nom) 

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

jetRecalib2017B = lambda : jetRecalib("Fall17_17Nov2017B_V6_DATA")
jetRecalib2017C = lambda : jetRecalib("Fall17_17Nov2017C_V6_DATA")
jetRecalib2017D = lambda : jetRecalib("Fall17_17Nov2017D_V6_DATA")
jetRecalib2017E = lambda : jetRecalib("Fall17_17Nov2017E_V6_DATA")
jetRecalib2017F = lambda : jetRecalib("Fall17_17Nov2017F_V6_DATA")
