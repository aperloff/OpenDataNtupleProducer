#ifndef JetAnalyzer_h
#define JetAnalyzer_h

#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Common/interface/ValueMap.h"

// File service for saving the ROOT files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// MiniAOD PAT libraries
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// ROOT libraries
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;
using namespace reco;
using namespace pat;

class JetAnalyzer : public edm::EDAnalyzer {
    public:
        explicit JetAnalyzer(const edm::ParameterSet&);
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        ~JetAnalyzer();

    private:
        // Configuration flags
        bool saveJetTree_;

        // Tokens
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> extraTracksToken_;
        edm::EDGetTokenT<reco::GenParticleCollection> genPartToken_;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
        edm::EDGetTokenT<reco::GenJetCollection> EDMGenJetsToken_;
        edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo >> pileupInfoToken_;

        edm::EDGetTokenT<double> pfRhoAllToken_;
        edm::EDGetTokenT<double> pfRhoCentralToken_;
        edm::EDGetTokenT<double> pfRhoCentralNeutralToken_;
        edm::EDGetTokenT<double> pfRhoCentralChargedPileUpToken_;

        edm::EDGetTokenT<edm::ValueMap<float> > qglToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > ptDToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > axis2Token_;
        edm::EDGetTokenT<edm::ValueMap<int> > multToken_;

        // Configurable vertex parameters
        double goodVtxNdof;
        double goodVtxZ;
        double goodVtxRho;

        // Output
        edm::Service<TFileService> fs;
        std::string t_qgtagger;

        TFile* outputFile;
        TTree* jetTree;
        TTree* evtTree;

        // -------------------------
        // TTree variables
        // -------------------------
        static const UInt_t kMaxPF = 5000;
        static const UInt_t kMaxTrk = kMaxPF;
        static const UInt_t kMaxVtx = kMaxTrk;

        // Jet variables
        Float_t jetPt;
        Float_t jetEta;
        Float_t jetPhi;
        Float_t jetMass;
        Float_t jetGirth;
        Float_t jetArea;

        Float_t jetRawPt;
        Float_t jetRawMass;

        UInt_t jetLooseID;
        UInt_t jetTightID;
        UInt_t jetGenMatch;

        // Quark-gluon likelihood variables
        Float_t jetQGl;
        Float_t QG_ptD;
        Float_t QG_axis2;
        UInt_t QG_mult;

        // Jet flavor variables
        Int_t partonFlav;
        Int_t hadronFlav;
        Int_t physFlav;

        UInt_t isPhysUDS;
        UInt_t isPhysG;
        UInt_t isPhysOther;
        UInt_t isPartonUDS;
        UInt_t isPartonG;
        UInt_t isPartonOther;

        // Jet composition
        UInt_t jetChargedHadronMult;
        UInt_t jetNeutralHadronMult;
        UInt_t jetChargedMult;
        UInt_t jetNeutralMult;
        UInt_t jetMult;

        // PF variables
        UInt_t nPF;
        Float_t PF_pT[kMaxPF];
        Float_t PF_dR[kMaxPF];
        Float_t PF_dTheta[kMaxPF];
        Float_t PF_dPhi[kMaxPF];
        Float_t PF_dEta[kMaxPF];
        Float_t PF_mass[kMaxPF];
        Int_t PF_id[kMaxPF];
        UInt_t PF_fromPV[kMaxPF];
        UInt_t PF_fromAK4Jet[kMaxPF];

        // Track variables
        UInt_t nTrk;
        Int_t Trk_charge[kMaxTrk];
        Float_t Trk_pT[kMaxTrk];
        Float_t Trk_eta[kMaxTrk];
        Float_t Trk_phi[kMaxTrk];
        Float_t Trk_px[kMaxTrk];
        Float_t Trk_py[kMaxTrk];
        Float_t Trk_pz[kMaxTrk];
        Float_t Trk_x[kMaxTrk];
        Float_t Trk_y[kMaxTrk];
        Float_t Trk_z[kMaxTrk];
        Float_t Trk_dxypv[kMaxTrk];
        Float_t Trk_dxyerrorpv[kMaxTrk];
        Float_t Trk_dzpv[kMaxTrk];
        Float_t Trk_dzerrorpv[kMaxTrk];
        UInt_t Trk_isLost[kMaxTrk];

        // Generator level primary vertex variables
        UInt_t nGenVtx;
        Float_t GenVtx_x[kMaxVtx];
        Float_t GenVtx_y[kMaxVtx];
        Float_t GenVtx_z[kMaxVtx];
        Float_t GenVtx_sumPt[kMaxVtx];
        Float_t GenVtx_sumPt2[kMaxVtx];

        // Primary vertex variables
        UInt_t nVtx;
        Float_t Vtx_x[kMaxVtx];
        Float_t Vtx_y[kMaxVtx];
        Float_t Vtx_z[kMaxVtx];
        Float_t Vtx_xError[kMaxVtx];
        Float_t Vtx_yError[kMaxVtx];
        Float_t Vtx_zError[kMaxVtx];
        UInt_t Vtx_ndof[kMaxVtx];
        Float_t Vtx_chi2[kMaxVtx];
        UInt_t Vtx_ntrks[kMaxVtx];
        UInt_t Vtx_isValid[kMaxVtx];
        UInt_t Vtx_isFake[kMaxVtx];
        UInt_t Vtx_isGood[kMaxVtx];

        // Generator level jet variables
        Float_t genJetPt;
        Float_t genJetEta;
        Float_t genJetPhi;
        Float_t genJetMass;

        // Generator level jet's particle variables
        UInt_t nGenJetPF;
        Float_t genJetPF_pT[kMaxPF];
        Float_t genJetPF_dR[kMaxPF];
        Float_t genJetPF_dTheta[kMaxPF];
        Float_t genJetPF_mass[kMaxPF];
        Int_t genJetPF_id[kMaxPF];

        // Misc. jet variables
        UInt_t eventJetMult;
        UInt_t jetPtOrder;

        // Selection variables
        Float_t dPhiJetsLO;
        Float_t dEtaJetsLO;
        Float_t alpha;

        // Event variables
        ULong64_t event;
        UInt_t run;
        UInt_t lumi;

        // MC variables
        Float_t pthat;
        Float_t eventWeight;

        // Pileup variables
        Float_t rhoAll;
        Float_t rhoCentral;
        Float_t rhoCentralNeutral;
        Float_t rhoCentralChargedPileUp;
        UInt_t PV_npvsGood;
        UInt_t Pileup_nPU;
        Float_t Pileup_nTrueInt;
};

// Define a struct for storing a jet with its index within the event (needed for QG likelihood variables)
struct JetIndexed {
	pat::Jet jet;
	unsigned int eventIndex;
	JetIndexed(pat::Jet j, unsigned int eIdx) : jet(j), eventIndex(eIdx) {}
};

// Define a sort function for JetIndexed pT-ordering
struct higher_pT_sort
{
	inline bool operator() (const JetIndexed& jet1, const JetIndexed& jet2)
	{
		return ( jet1.jet.pt() > jet2.jet.pt() );
	}
};

#endif
