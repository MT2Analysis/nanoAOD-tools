import ROOT
import os, types
from math import *
from PhysicsTools.HeppyCore.utils.deltar import *

class JetReCalibrator:
    def __init__(self,globalTag,jetFlavour,doResidualJECs,jecPath,upToLevel=3,
                 calculateSeparateCorrections=False,
                 calculateType1METCorrection=False, type1METParams={'jetPtThreshold':15., 'skipEMfractionThreshold':0.9, 'skipMuons':True} ):
        """Create a corrector object that reads the payloads from the text dumps of a global tag under
            CMGTools/RootTools/data/jec  (see the getJec.py there to make the dumps).
           It will apply the L1,L2,L3 and possibly the residual corrections to the jets.
           If configured to do so, it will also compute the type1 MET corrections."""
        self.globalTag = globalTag
        self.jetFlavour = jetFlavour
        self.doResidualJECs = doResidualJECs
        self.jecPath = jecPath
        self.upToLevel = upToLevel
        self.calculateType1METCorr = calculateType1METCorrection
        self.type1METParams  = type1METParams
        # Make base corrections
        path = os.path.expandvars(jecPath) #"%s/src/CMGTools/RootTools/data/jec" % os.environ['CMSSW_BASE'];
        self.L1JetPar  = ROOT.JetCorrectorParameters("%s/%s_L1FastJet_%s.txt" % (path,globalTag,jetFlavour),"");
        self.L2JetPar  = ROOT.JetCorrectorParameters("%s/%s_L2Relative_%s.txt" % (path,globalTag,jetFlavour),"");
        self.L3JetPar  = ROOT.JetCorrectorParameters("%s/%s_L3Absolute_%s.txt" % (path,globalTag,jetFlavour),"");
        self.vPar = ROOT.vector(ROOT.JetCorrectorParameters)()
        self.vPar.push_back(self.L1JetPar);
        if upToLevel >= 2: self.vPar.push_back(self.L2JetPar);
        if upToLevel >= 3: self.vPar.push_back(self.L3JetPar);
        # Add residuals if needed
        if doResidualJECs : 
            self.ResJetPar = ROOT.JetCorrectorParameters("%s/%s_L2L3Residual_%s.txt" % (path,globalTag,jetFlavour))
            self.vPar.push_back(self.ResJetPar);
        #Step3 (Construct a FactorizedJetCorrector object) 
        self.JetCorrector = ROOT.FactorizedJetCorrector(self.vPar)
        if os.path.exists("%s/%s_Uncertainty_%s.txt" % (path,globalTag,jetFlavour)):
            self.JetUncertainty = ROOT.JetCorrectionUncertainty("%s/%s_Uncertainty_%s.txt" % (path,globalTag,jetFlavour));
        elif os.path.exists("%s/Uncertainty_FAKE.txt" % path):
            self.JetUncertainty = ROOT.JetCorrectionUncertainty("%s/Uncertainty_FAKE.txt" % path);
        else:
            print 'Missing JEC uncertainty file "%s/%s_Uncertainty_%s.txt", so jet energy uncertainties will not be available' % (path,globalTag,jetFlavour)
            self.JetUncertainty = None
        self.separateJetCorrectors = {}
        if calculateSeparateCorrections or calculateType1METCorrection:
            self.vParL1 = ROOT.vector(ROOT.JetCorrectorParameters)()
            self.vParL1.push_back(self.L1JetPar)
            self.separateJetCorrectors["L1"] = ROOT.FactorizedJetCorrector(self.vParL1)
            if upToLevel >= 2 and calculateSeparateCorrections:
                self.vParL2 = ROOT.vector(ROOT.JetCorrectorParameters)()
                for i in [self.L1JetPar,self.L2JetPar]: self.vParL2.push_back(i)
                self.separateJetCorrectors["L1L2"] = ROOT.FactorizedJetCorrector(self.vParL2)
            if upToLevel >= 3 and calculateSeparateCorrections:
                self.vParL3 = ROOT.vector(ROOT.JetCorrectorParameters)()
                for i in [self.L1JetPar,self.L2JetPar,self.L3JetPar]: self.vParL3.push_back(i)
                self.separateJetCorrectors["L1L2L3"] = ROOT.FactorizedJetCorrector(self.vParL3)
            if doResidualJECs and calculateSeparateCorrections:
                self.vParL3Res = ROOT.vector(ROOT.JetCorrectorParameters)()
                for i in [self.L1JetPar,self.L2JetPar,self.L3JetPar,self.ResJetPar]: self.vParL3Res.push_back(i)
                self.separateJetCorrectors["L1L2L3Res"] = ROOT.FactorizedJetCorrector(self.vParL3Res)

    def getCorrection(self,jet,rho,delta=0,corrector=None):
        if not corrector: corrector = self.JetCorrector
        if corrector != self.JetCorrector and delta!=0: raise RuntimeError('Configuration not supported')
        corrector.setJetEta(jet.eta)
        corrector.setJetPt(jet.pt*(1.-jet.rawFactor))
        corrector.setJetA(jet.area)
        corrector.setRho(rho)
        corr = corrector.getCorrection()
        if delta != 0:
            if not self.JetUncertainty: raise RuntimeError("Jet energy scale uncertainty shifts requested, but not available")
            self.JetUncertainty.setJetEta(jet.eta())
            self.JetUncertainty.setJetPt(corr * jet.pt() * (1-jet.rawFactor()))
            try:
                jet.jetEnergyCorrUncertainty = self.JetUncertainty.getUncertainty(True) 
            except RuntimeError as r:
                print "Caught %s when getting uncertainty for jet of pt %.1f, eta %.2f\n" % (r,corr * jet.pt() * (1-jet.rawFactor()),jet.eta())
                jet.jetEnergyCorrUncertainty = 0.5
            #print "   jet with corr pt %6.2f has an uncertainty %.2f " % (jet.pt()*jet.rawFactor()*corr, jet.jetEnergyCorrUncertainty)
            corr *= max(0, 1+delta*jet.jetEnergyCorrUncertainty)
        #print 'orig_pt={}, raw_factor={}, raw_pt={}, rho={}, corr_factor={}, final_pt={}'.format(jet.pt, (1.-jet.rawFactor), jet.pt*(1.-jet.rawFactor), rho, corr, jet.pt*(1.-jet.rawFactor)*corr)
        return corr


    def rawP4forType1MET_(self, jet, muons):
        """Return the raw 4-vector, after subtracting the muons (if requested),
           or None if the jet fails the EMF cut."""
        p4 = jet.p4() * (1-jet.rawFactor) # note different definition of rawfactor wrt heppy
        #muons = Collection(event, "Muon")
        #if not jet.hasPFSpecific():
            # return raw 4-vector if there are no PF details (e.g. AK8 jet below threshold) Not possible on nanoAOD
        #    return p4
        #emf = ( jet.physObj.neutralEmEnergy() + jet.physObj.chargedEmEnergy() )/p4.E()
        #print 'debug in rawP4forType1MET_'
        emf = (jet.neEmEF + jet.chEmEF) / p4.E()
        if emf > self.type1METParams['skipEMfractionThreshold']:
            return None
        if self.type1METParams['skipMuons']:
            #print 'debug if skipMuons'
            #for idau in xrange(jet.numberOfDaughters()):
            #    pfcand = jet.daughter(idau)
            #    if pfcand.isGlobalMuon() or pfcand.isStandAloneMuon(): 
            #        p4 -= pfcand.p4()
            # NOTE: above cannot be implemented in nanoAOD 
            # -> check instead overlap with available muons
            if jet.muonIdx1!=-1:
              overlap_mu = muons[jet.muonIdx1]
              if overlap_mu.isGlobal or not overlap_mu.isTracker: # equivalent to global || standalone
            #    print 'going to subtract first muon' 
                p4 = - overlap_mu.p4()
            if jet.muonIdx2!=-1:
              overlap_mu = muons[jet.muonIdx2]
              if overlap_mu.isGlobal or not overlap_mu.isTracker: # equivalent to global || standalone
            #    print 'going to subtract second muon' 
                p4 = - overlap_mu.p4()
        #print p4.Pt()
        return p4



# TODO: add metShift, and T1 met correction business
    def correct(self,jet,muons,rho,delta=0,addCorr=False,addShifts=False, metShift=None, type1METCorr=None):
        """Corrects a jet energy (optionally shifting it also by delta times the JEC uncertainty)

           If addCorr, set jet.corr to the correction.
           If addShifts, set also the +1 and -1 jet shifts 

           The metShift vector will accumulate the x and y changes to the MET from the JEC, i.e. the 
           negative difference between the new and old jet momentum, for jets eligible for type1 MET 
           corrections, and after subtracting muons. The pt cut is applied to the new corrected pt.
           This shift can be applied on top of the *OLD TYPE1 MET*, but only if there was no change 
           in the L1 corrections nor in the definition of the type1 MET (e.g. jet pt cuts).

           The type1METCorr vector, will accumulate the x, y, sumEt type1 MET corrections, to be
           applied to the *RAW MET*
        """
        raw = 1.-jet.rawFactor  
        corr = self.getCorrection(jet,rho,delta)
        if corr <= 0:
            return jet.pt
        newpt = jet.pt*raw*corr

        # MET recalibration business
        if not metShift:
            metShift = [0,0]
        if not type1METCorr:
            type1METCorr = [0,0]
        if newpt > self.type1METParams['jetPtThreshold']:
            rawP4forT1 = self.rawP4forType1MET_(jet,muons)
            if rawP4forT1 and rawP4forT1.Pt()*corr > self.type1METParams['jetPtThreshold']:
                metShift[0] -= rawP4forT1.Px() * (corr - 1.0/raw)
                metShift[1] -= rawP4forT1.Py() * (corr - 1.0/raw)
                if self.calculateType1METCorr:
                    l1corr = self.getCorrection(jet,rho,delta=0,corrector=self.separateJetCorrectors["L1"])
                    print 'debug corr={}, L1corr={}, rawP4forT1.Px()={}, rawP4forT1.Py()={}'.format(corr, l1corr, rawP4forT1.Px(), rawP4forT1.Py())
                    #print "\tfor jet with raw pt %.5g, eta %.5g, dpx = %.5g, dpy = %.5g" % (
                    #            jet.pt()*raw, jet.eta(), 
                    #            rawP4forT1.Px()*(corr - l1corr), 
                    #            rawP4forT1.Py()*(corr - l1corr))
                    type1METCorr[0] -= rawP4forT1.Px() * (corr - l1corr) 
                    type1METCorr[1] -= rawP4forT1.Py() * (corr - l1corr) 
                    #type1METCorr[2] += rawP4forT1.Et() * (corr - l1corr) 
        print 'orig_pt={}, raw_factor={}, raw_pt={}, rho={}, corr_factor={}, final_pt={}'.format(jet.pt, raw, jet.pt*raw, rho, corr, jet.pt*raw*corr)
        print 'debug, metShiftX={}, metShiftY={}, type1METCorrX={}, type1METCorrY={}'.format(metShift[0],metShift[1],type1METCorr[0],type1METCorr[1])
        return newpt


