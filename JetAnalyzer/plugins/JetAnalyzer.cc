//  Jet tuple producer for 13 TeV Run2 MC samples
//  Specifically aimed at studies of gluon and light quark jets
//  Data is saved to file on a jet-by-jet basis, resulting in almost flat tuples
//
//  Author: Kimmo Kallonen
//  Based on previous work by: Petra-Maria Ekroos

#include "JetAnalyzer.h"

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig):
    saveJetTree_(iConfig.getParameter<bool>("saveJetTree")),
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    extraTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("extraTracks"))),
    genPartToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    EDMGenJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genEventInfoToken_(consumes <GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    pileupInfoToken_(consumes <std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("pileupInfo"))),
    pfRhoAllToken_(consumes <double> (iConfig.getParameter<edm::InputTag>("pfRhoAll"))),
    pfRhoCentralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("pfRhoCentral"))),
    pfRhoCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("pfRhoCentralNeutral"))),
    pfRhoCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("pfRhoCentralChargedPileUp"))),
    qglToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    ptDToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
    axis2Token_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multToken_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult")))
{
    goodVtxNdof = iConfig.getParameter<double>("confGoodVtxNdof");
    goodVtxZ = iConfig.getParameter<double>("confGoodVtxZ");
    goodVtxRho = iConfig.getParameter<double>("confGoodVtxRho");
}

JetAnalyzer::~JetAnalyzer()
{}

void JetAnalyzer::beginJob()
{
    if(saveJetTree_) {
      // Create the ROOT tree and add all the branches to it
      jetTree = fs->make<TTree>("jetTree", "jetTree");

      jetTree->Branch("jetPt", &jetPt, "jetPt/F");
      jetTree->Branch("jetEta", &jetEta, "jetEta/F");
      jetTree->Branch("jetPhi", &jetPhi, "jetPhi/F");
      jetTree->Branch("jetMass", &jetMass, "jetMass/F");
      jetTree->Branch("jetGirth", &jetGirth, "jetGirth/F");
      jetTree->Branch("jetArea", &jetArea, "jetArea/F");

      jetTree->Branch("jetRawPt", &jetRawPt, "jetRawPt/F");
      jetTree->Branch("jetRawMass", &jetRawMass, "jetRawMass/F");

      jetTree->Branch("jetLooseID", &jetLooseID, "jetLooseID/I");
      jetTree->Branch("jetTightID", &jetTightID, "jetTightID/I");
      jetTree->Branch("jetGenMatch", &jetGenMatch, "jetGenMatch/I");

      jetTree->Branch("jetQGl", &jetQGl, "jetQGl/F");
      jetTree->Branch("QG_ptD", &QG_ptD, "QG_ptD/F");
      jetTree->Branch("QG_axis2", &QG_axis2, "QG_axis2/F");
      jetTree->Branch("QG_mult", &QG_mult, "QG_mult/I");

      jetTree->Branch("partonFlav", &partonFlav, "partonFlav/I");
      jetTree->Branch("hadronFlav", &hadronFlav, "hadronFlav/I");
      jetTree->Branch("physFlav", &physFlav, "physFlav/I");

      jetTree->Branch("isPhysUDS", &isPhysUDS, "isPhysUDS/I");
      jetTree->Branch("isPhysG", &isPhysG, "isPhysG/I");
      jetTree->Branch("isPhysOther", &isPhysOther, "isPhysOther/I");
      jetTree->Branch("isPartonUDS", &isPartonUDS, "isPartonUDS/I");
      jetTree->Branch("isPartonG", &isPartonG, "isPartonG/I");
      jetTree->Branch("isPartonOther", &isPartonOther, "isPartonOther/I");

      jetTree->Branch("jetChargedHadronMult", &jetChargedHadronMult, "jetChargedHadronMult/I");
      jetTree->Branch("jetNeutralHadronMult", &jetNeutralHadronMult, "jetNeutralHadronMult/I");
      jetTree->Branch("jetChargedMult", &jetChargedMult, "jetChargedMult/I");
      jetTree->Branch("jetNeutralMult", &jetNeutralMult, "jetNeutralMult/I");
      jetTree->Branch("jetMult", &jetMult, "jetMult/I");

      jetTree->Branch("nPF", &nPF, "nPF/I");
      jetTree->Branch("PF_pT", &PF_pT, "PF_pT[nPF]/F");
      jetTree->Branch("PF_dR", &PF_dR, "PF_dR[nPF]/F");
      jetTree->Branch("PF_dTheta", &PF_dTheta, "PF_dTheta[nPF]/F");
      jetTree->Branch("PF_dPhi", &PF_dPhi, "PF_dPhi[nPF]/F");
      jetTree->Branch("PF_dEta", &PF_dEta, "PF_dEta[nPF]/F");
      jetTree->Branch("PF_mass", &PF_mass, "cPF_mass[nPF]/F");
      jetTree->Branch("PF_id", &PF_id, "PF_id[nPF]/I");
      jetTree->Branch("PF_fromPV", &PF_fromPV, "PF_fromPV[nPF]/I");
      jetTree->Branch("PF_fromAK4Jet", &PF_fromAK4Jet, "PF_fromAK4Jet[nPF]/I");

      jetTree->Branch("genJetPt", &genJetPt, "genJetPt/F");
      jetTree->Branch("genJetEta", &genJetEta, "genJetEta/F");
      jetTree->Branch("genJetPhi", &genJetPhi, "genJetPhi/F");
      jetTree->Branch("genJetMass", &genJetMass, "genJetMass/F");

      jetTree->Branch("nGenJetPF",&nGenJetPF,"nGenJetPF/I");
      jetTree->Branch("genJetPF_pT", &genJetPF_pT, "genJetPF_pT[nGenJetPF]/F");
      jetTree->Branch("genJetPF_dR", &genJetPF_dR, "genJetPF_dR[nGenJetPF]/F");
      jetTree->Branch("genJetPF_dTheta", &genJetPF_dTheta, "genJetPF_dTheta[nGenJetPF]/F");
      jetTree->Branch("genJetPF_mass", &genJetPF_mass, "genJetPF_mass[nGenJetPF]/F");
      jetTree->Branch("genJetPF_id", &genJetPF_id, "genJetPF_id[nGenJetPF]/I");

      jetTree->Branch("eventJetMult", &eventJetMult, "eventJetMult/I");
      jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");

      jetTree->Branch("dPhiJetsLO", &dPhiJetsLO, "dPhiJetsLO/F");
      jetTree->Branch("dEtaJetsLO", &dEtaJetsLO, "dEtaJetsLO/F");
      jetTree->Branch("alpha", &alpha, "alpha/F");

      jetTree->Branch("event", &event, "event/l");
      jetTree->Branch("run", &run, "run/I");
      jetTree->Branch("lumi", &lumi, "lumi/I");

      // Add descriptive comments to all of the branches
      jetTree->GetBranch("jetPt")->SetTitle("Transverse momentum of the jet");
      jetTree->GetBranch("jetEta")->SetTitle("Pseudorapidity of the jet");
      jetTree->GetBranch("jetPhi")->SetTitle("Azimuthal angle of the jet");
      jetTree->GetBranch("jetMass")->SetTitle("Mass of the jet");
      jetTree->GetBranch("jetGirth")->SetTitle("Girth of the jet (as defined in arXiv:1106.3076 [hep-ph])");
      jetTree->GetBranch("jetArea")->SetTitle("Catchment area of the jet; used for jet energy corrections");
      jetTree->GetBranch("jetRawPt")->SetTitle("Transverse momentum of the jet before energy corrections");
      jetTree->GetBranch("jetRawMass")->SetTitle("Mass of the jet before energy corrections");
      jetTree->GetBranch("jetLooseID")->SetTitle("Indicates if the jet passes loose selection criteria; used for dismissing fake jets");
      jetTree->GetBranch("jetTightID")->SetTitle("Indicates if the jet passes tight selection criteria; used for dismissing fake jets");
      jetTree->GetBranch("jetGenMatch")->SetTitle("1: if a matched generator level jet exists; 0: if no match was found");
      jetTree->GetBranch("jetQGl")->SetTitle("Quark vs Gluon likelihood discriminator");
      jetTree->GetBranch("QG_ptD")->SetTitle("Transverse momentum distribution among particle flow candidates within the jet; defined as sqrt(Sum(pt^2))/Sum(pt), where the sum is over the particle flow candidates of the jet");
      jetTree->GetBranch("QG_axis2")->SetTitle("Minor axis of the jet, calculated from the particle flow candidates");
      jetTree->GetBranch("QG_mult")->SetTitle("Multiplicity of jet constituents with additional cuts: 1 GeV transverse momentum threshold for neutral particles; charged particles required to be associated with the primary interaction vertex (PF_fromPV == 3)");
      jetTree->GetBranch("partonFlav")->SetTitle("Flavor of the jet; parton definition");
      jetTree->GetBranch("hadronFlav")->SetTitle("Flavor of the jet; hadron definition");
      jetTree->GetBranch("physFlav")->SetTitle("Flavor of the jet; physics definition");
      jetTree->GetBranch("isPhysUDS")->SetTitle("Indicates a light quark jet; |physFlav| == 1, 2, 3");
      jetTree->GetBranch("isPhysG")->SetTitle("Indicates a gluon jet; physFlav == 21");
      jetTree->GetBranch("isPhysOther")->SetTitle("Indicates a non-light quark/gluon jet; |physFlav| != 1, 2, 3, 21");
      jetTree->GetBranch("isPartonUDS")->SetTitle("Indicates a light quark jet; |physFlav| == 1, 2, 3");
      jetTree->GetBranch("isPartonG")->SetTitle("Indicates a gluon jet; physFlav == 21");
      jetTree->GetBranch("isPartonOther")->SetTitle("Indicates a non-light quark/gluon jet; |physFlav| != 1, 2, 3, 21");
      jetTree->GetBranch("jetChargedHadronMult")->SetTitle("Multiplicity of charged hadron jet constituents");
      jetTree->GetBranch("jetNeutralHadronMult")->SetTitle("Multiplicity of neutral hadron jet constituents");
      jetTree->GetBranch("jetChargedMult")->SetTitle("Multiplicity of charged jet constituents");
      jetTree->GetBranch("jetNeutralMult")->SetTitle("Multiplicity of neutral jet constituents");
      jetTree->GetBranch("jetMult")->SetTitle("Multiplicity of jet constituents");
      jetTree->GetBranch("nPF")->SetTitle("Number of particle flow candidates (particles reconstructed by the particle flow algorithm); contains all particles within |deltaPhi| < 1 && |deltaEta| < 1 from the center of the jet");
      jetTree->GetBranch("PF_pT")->SetTitle("Transverse momentum of a particle flow candidate");
      jetTree->GetBranch("PF_dR")->SetTitle("Distance of a particle flow candidate to the center of the jet");
      jetTree->GetBranch("PF_dTheta")->SetTitle("Polar angle of a particle flow candidate");
      jetTree->GetBranch("PF_dPhi")->SetTitle("Azimuthal angle of a particle flow candidate");
      jetTree->GetBranch("PF_dEta")->SetTitle("Pseudorapidity of a particle flow candidate");
      jetTree->GetBranch("PF_mass")->SetTitle("Mass of a particle flow candidate");
      jetTree->GetBranch("PF_id")->SetTitle("Generator level particle identifier for the particle flow candidates, as defined in the PDG particle numbering scheme");
      jetTree->GetBranch("PF_fromPV")->SetTitle("Indicates how tightly the particle is associated with the primary vertex; ranges from 3 to 0");
      jetTree->GetBranch("PF_fromAK4Jet")->SetTitle("1: if the particle flow candidate is a constituent of the reconstructed AK4 jet; 0: if it is not a constituent of the jet");
      jetTree->GetBranch("genJetPt")->SetTitle("Transverse momentum of the matched generator level jet");
      jetTree->GetBranch("genJetEta")->SetTitle("Pseudorapidity of the matched generator level jet");
      jetTree->GetBranch("genJetPhi")->SetTitle("Azimuthal angle of the matched generator level jet");
      jetTree->GetBranch("genJetMass")->SetTitle("Mass of the matched generator level jet");
      jetTree->GetBranch("nGenJetPF")->SetTitle("Number of particles in the matched generator level jet");
      jetTree->GetBranch("genJetPF_pT")->SetTitle("Transverse momentum of a particle in the matched generator level jet");
      jetTree->GetBranch("genJetPF_dR")->SetTitle("Distance of a particle to the center of the matched generator level jet ");
      jetTree->GetBranch("genJetPF_dTheta")->SetTitle("Polar angle of a particle in the matched generator level jet");
      jetTree->GetBranch("genJetPF_mass")->SetTitle("Mass of a particle in the matched generator level jet");
      jetTree->GetBranch("genJetPF_id")->SetTitle("Generator level particle identifier for the particles in the matched generator level jet, as defined in the PDG particle numbering scheme");
      jetTree->GetBranch("eventJetMult")->SetTitle("Multiplicity of jets in the event");
      jetTree->GetBranch("jetPtOrder")->SetTitle("Indicates the ranking number of the jet, as the jets are ordered by their transverse momenta within the event");
      jetTree->GetBranch("dPhiJetsLO")->SetTitle("Phi difference of the two leading jets");
      jetTree->GetBranch("dEtaJetsLO")->SetTitle("Eta difference of the two leading jets");
      jetTree->GetBranch("alpha")->SetTitle("If there are at least 3 jets in the event, alpha is the third jet's transverse momentum divided by the average transverse momentum of the two leading jets");
      jetTree->GetBranch("event")->SetTitle("Event number");
      jetTree->GetBranch("run")->SetTitle("Run number");
      jetTree->GetBranch("lumi")->SetTitle("Luminosity block");
    }

    evtTree = fs->make<TTree>("evtTree", "evtTree");

    evtTree->Branch("event", &event, "event/l");
    evtTree->Branch("run", &run, "run/I");
    evtTree->Branch("lumi", &lumi, "lumi/I");

    evtTree->Branch("pthat", &pthat, "pthat/F");
    evtTree->Branch("eventWeight", &eventWeight, "eventWeight/F");

    evtTree->Branch("rhoAll", &rhoAll, "rhoAll/F");
    evtTree->Branch("rhoCentral", &rhoCentral, "rhoCentral/F");
    evtTree->Branch("rhoCentralNeutral", &rhoCentralNeutral, "rhoCentralNeutral/F");
    evtTree->Branch("rhoCentralChargedPileUp", &rhoCentralChargedPileUp, "rhoCentralChargedPileUp/F");
    evtTree->Branch("PV_npvsGood", &PV_npvsGood, "PV_npvsGood/I");
    evtTree->Branch("Pileup_nPU", &Pileup_nPU, "Pileup_nPU/I");
    evtTree->Branch("Pileup_nTrueInt", &Pileup_nTrueInt, "Pileup_nTrueInt/F");

    evtTree->Branch("nTrk", &nTrk, "nTrk/I");
    evtTree->Branch("Trk_charge", &Trk_charge, "Trk_charge[nTrk]/I");
    evtTree->Branch("Trk_pT", &Trk_pT, "Trk_pT[nTrk]/F");
    evtTree->Branch("Trk_eta", &Trk_eta, "Trk_eta[nTrk]/F");
    evtTree->Branch("Trk_phi", &Trk_phi, "Trk_phi[nTrk]/F");
    evtTree->Branch("Trk_px", &Trk_px, "Trk_px[nTrk]/F");
    evtTree->Branch("Trk_py", &Trk_py, "Trk_py[nTrk]/F");
    evtTree->Branch("Trk_pz", &Trk_pz, "Trk_pz[nTrk]/F");
    evtTree->Branch("Trk_x", &Trk_x, "Trk_x[nTrk]/F");
    evtTree->Branch("Trk_y", &Trk_y, "Trk_y[nTrk]/F");
    evtTree->Branch("Trk_z", &Trk_z, "Trk_z[nTrk]/F");
    evtTree->Branch("Trk_dxypv", &Trk_dxypv, "Trk_dxypv[nTrk]/F");
    evtTree->Branch("Trk_dxyerrorpv", &Trk_dxyerrorpv, "Trk_dxyerrorpv[nTrk]/F");
    evtTree->Branch("Trk_dzpv", &Trk_dzpv, "Trk_dzpv[nTrk]/F");
    evtTree->Branch("Trk_dzerrorpv", &Trk_dzerrorpv, "Trk_dzerrorpv[nTrk]/F");
    evtTree->Branch("Trk_isLost", &Trk_isLost, "Trk_isLost[nTrk]/I");

    evtTree->Branch("nGenVtx", &nGenVtx, "nGenVtx/I");
    evtTree->Branch("GenVtx_x", &GenVtx_x, "GenVtx_x[nGenVtx]/F");
    evtTree->Branch("GenVtx_y", &GenVtx_y, "GenVtx_y[nGenVtx]/F");
    evtTree->Branch("GenVtx_z", &GenVtx_z, "GenVtx_z[nGenVtx]/F");
    evtTree->Branch("GenVtx_sumPt", &GenVtx_sumPt, "GenVtx_sumPt[nGenVtx]/F");
    evtTree->Branch("GenVtx_sumPt2", &GenVtx_sumPt2, "GenVtx_sumPt2[nGenVtx]/F");

    evtTree->Branch("nVtx", &nVtx, "nVtx/I");
    evtTree->Branch("Vtx_x", &Vtx_x, "Vtx_x[nVtx]/F");
    evtTree->Branch("Vtx_y", &Vtx_y, "Vtx_y[nVtx]/F");
    evtTree->Branch("Vtx_z", &Vtx_z, "Vtx_z[nVtx]/F");
    evtTree->Branch("Vtx_xError", &Vtx_xError, "Vtx_xError[nVtx]/F");
    evtTree->Branch("Vtx_yError", &Vtx_yError, "Vtx_yError[nVtx]/F");
    evtTree->Branch("Vtx_zError", &Vtx_zError, "Vtx_zError[nVtx]/F");
    evtTree->Branch("Vtx_ndof", &Vtx_ndof, "Vtx_ndof[nVtx]/I");
    evtTree->Branch("Vtx_chi2", &Vtx_chi2, "Vtx_chi2[nVtx]/F");
    evtTree->Branch("Vtx_ntrks", &Vtx_ntrks, "Vtx_ntrks[nVtx]/I");
    evtTree->Branch("Vtx_isValid", &Vtx_isValid, "Vtx_isValid[nVtx]/I");
    evtTree->Branch("Vtx_isFake", &Vtx_isFake, "Vtx_isFake[nVtx]/I");
    evtTree->Branch("Vtx_isGood", &Vtx_isGood, "Vtx_isGood[nVtx]/I");

    // Add descriptive comments to all of the branches
    evtTree->GetBranch("event")->SetTitle("Event number");
    evtTree->GetBranch("run")->SetTitle("Run number");
    evtTree->GetBranch("lumi")->SetTitle("Luminosity block");
    evtTree->GetBranch("pthat")->SetTitle("Transverse momentum of the generated hard process");
    evtTree->GetBranch("eventWeight")->SetTitle("Monte Carlo generator weight");
    evtTree->GetBranch("rhoAll")->SetTitle("The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates");
    evtTree->GetBranch("rhoCentral")->SetTitle("The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates with |eta|<2.5");
    evtTree->GetBranch("rhoCentralNeutral")->SetTitle("The median density (in GeV/A) of pile-up contamination per event; computed from all neutral PF candidates with |eta| < 2.5");
    evtTree->GetBranch("rhoCentralChargedPileUp")->SetTitle("The median density (in GeV/A) of pile-up contamination per event; computed from all PF charged hadrons associated to pileup vertices and with |eta| < 2.5");
    evtTree->GetBranch("PV_npvsGood")->SetTitle("The number of good reconstructed primary vertices; selection: !isFake && ndof > 4 && abs(z) <= 24 && position.Rho < 2");
    evtTree->GetBranch("Pileup_nPU")->SetTitle("The number of pileup interactions that have been added to the event in the current bunch crossing");
    evtTree->GetBranch("Pileup_nTrueInt")->SetTitle("The true mean number of the poisson distribution for this event from which the number of interactions in each bunch crossing has been sampled");

    evtTree->GetBranch("nTrk")->SetTitle("Number of tracks; contains all tracks matched to particle flow candidates and the lost tracks (tracks not matched to particle flow candidates)");
    evtTree->GetBranch("Trk_charge")->SetTitle("Charge of a track");
    evtTree->GetBranch("Trk_pT")->SetTitle("Transverse momentum of a track [GeV]");
    evtTree->GetBranch("Trk_eta")->SetTitle("Pseudorapidity of a track");
    evtTree->GetBranch("Trk_phi")->SetTitle("Azimuthal angle of a track");
    evtTree->GetBranch("Trk_px")->SetTitle("x coordinate of the track mometum vector [cm]");
    evtTree->GetBranch("Trk_py")->SetTitle("y coordinate of the track mometum vector [cm]");
    evtTree->GetBranch("Trk_pz")->SetTitle("z coordinate of the track mometum vector [cm]");
    evtTree->GetBranch("Trk_x")->SetTitle("Position of the track along the x direction at its reference point [cm]");
    evtTree->GetBranch("Trk_y")->SetTitle("Position of the track along the y direction at its reference point [cm]");
    evtTree->GetBranch("Trk_z")->SetTitle("Position of the track along the beamline at its reference point [cm]");
    evtTree->GetBranch("Trk_dxypv")->SetTitle("Transverse impact parameter of the track [cm]");
    evtTree->GetBranch("Trk_dxyerrorpv")->SetTitle("Error on the transverse impact parameter of the track [cm]");
    evtTree->GetBranch("Trk_dzpv")->SetTitle("Longitudinal impact parameter of the track [cm]");
    evtTree->GetBranch("Trk_dzerrorpv")->SetTitle("Error on the longitudinal impact parameter of the track [cm]");
    evtTree->GetBranch("Trk_isLost")->SetTitle("Indicates if the track was part of the lostTracks collections (i.e. not matched to a PF candidate)");

    evtTree->GetBranch("nGenVtx")->SetTitle("Number of generator level primary vertices");
    evtTree->GetBranch("GenVtx_x")->SetTitle("x coordinate of a generator level primary vertex [cm]");
    evtTree->GetBranch("GenVtx_y")->SetTitle("y coordinate of a generator level primary vertex [cm]");
    evtTree->GetBranch("GenVtx_z")->SetTitle("z coordinate of a generator level primary vertex [cm]");
    evtTree->GetBranch("GenVtx_sumPt")->SetTitle("Sum of the pT of the generator particles associated to this vertex [GeV]");
    evtTree->GetBranch("GenVtx_sumPt2")->SetTitle("Sum of the pT^2 of the generator particles associated to this vertex [GeV^2]");


    evtTree->GetBranch("nVtx")->SetTitle("Number of primary vertices");
    evtTree->GetBranch("Vtx_x")->SetTitle("x coordinate of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_y")->SetTitle("y coordinate of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_z")->SetTitle("z coordinate of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_xError")->SetTitle("Error on the x positon of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_yError")->SetTitle("Error on the y positon of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_zError")->SetTitle("Error on the z positon of a primary vertex [cm]");
    evtTree->GetBranch("Vtx_ndof")->SetTitle("Number of degrees of freedom of the primary vertex fit");
    evtTree->GetBranch("Vtx_chi2")->SetTitle("Chi square of the primary vertex fit");
    evtTree->GetBranch("Vtx_ntrks")->SetTitle("Number of tracks associated to a primary vertex");
    evtTree->GetBranch("Vtx_isValid")->SetTitle("Indicates if a primary vertex is valid");
    evtTree->GetBranch("Vtx_isFake")->SetTitle("Indicates if a primary vertex is fake");
    evtTree->GetBranch("Vtx_isGood")->SetTitle("Indicates if a primary vertex is good");
}

void JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
    edm::Handle<pat::PackedCandidateCollection> extraTrks;
    iEvent.getByToken(extraTracksToken_, extraTrks);
    edm::Handle<reco::GenParticleCollection> genParts;
    iEvent.getByToken(genPartToken_, genParts);
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(EDMGenJetsToken_, genJets);
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    edm::Handle<std::vector< PileupSummaryInfo >>  puInfo;
    iEvent.getByToken(pileupInfoToken_, puInfo);

    edm::Handle<double> pfRhoAllHandle;
    iEvent.getByToken(pfRhoAllToken_, pfRhoAllHandle);
    edm::Handle<double> pfRhoCentralHandle;
    iEvent.getByToken(pfRhoCentralToken_, pfRhoCentralHandle);
    edm::Handle<double> pfRhoCentralNeutralHandle;
    iEvent.getByToken(pfRhoCentralNeutralToken_, pfRhoCentralNeutralHandle);
    edm::Handle<double> pfRhoCentralChargedPileUpHandle;
    iEvent.getByToken(pfRhoCentralChargedPileUpToken_, pfRhoCentralChargedPileUpHandle);

    edm::Handle<edm::ValueMap<float>> qglHandle;
    iEvent.getByToken(qglToken_, qglHandle);
    edm::Handle<edm::ValueMap<float>> ptDHandle;
    iEvent.getByToken(ptDToken_, ptDHandle);
    edm::Handle<edm::ValueMap<float>> axis2Handle;
    iEvent.getByToken(axis2Token_, axis2Handle);
    edm::Handle<edm::ValueMap<int>> multHandle;
    iEvent.getByToken(multToken_, multHandle);

    if (saveJetTree_) {

      // Create vectors for the jets
      // sortedJets include all jets of the event, while selectedJets have pT and eta cuts
      vector<JetIndexed> sortedJets;
      vector<JetIndexed> selectedJets;

      // Loop over the jets to save them to the jet vectors for pT-ordering
      int iJetR = -1;
      for(pat::JetCollection::const_iterator jetIt = jets->begin(); jetIt!=jets->end(); ++jetIt) {
        const pat::Jet &jet = *jetIt;
        ++iJetR;
        sortedJets.push_back( JetIndexed( jet, iJetR) );
        // Select
        if ( (jet.pt() > 30) && (fabs(jet.eta()) < 2.5) ) {
	  selectedJets.push_back( JetIndexed( jet, iJetR) );
        }
      }

      // Sort the jets in pT order
      std::sort(sortedJets.begin(), sortedJets.end(), higher_pT_sort());
      std::sort(selectedJets.begin(), selectedJets.end(), higher_pT_sort());

      // Loop over the pT-ordered selected jets and save them to file
      for (unsigned int ptIdx = 0; ptIdx < selectedJets.size(); ++ptIdx) {
        // Make selective cuts on the event level
        if (sortedJets.size() < 2) continue;
        if (fabs(sortedJets[0].jet.eta()) > 2.5 || fabs(sortedJets[1].jet.eta()) > 2.5) continue;
        if (fabs(sortedJets[0].jet.pt()) < 30 || fabs(sortedJets[1].jet.pt()) < 30) continue;

        JetIndexed idxJet = selectedJets[ptIdx];
        const pat::Jet j = idxJet.jet;
        int iJetRef = idxJet.eventIndex;

        // Jet variables
        jetPt = j.pt();
        jetEta = j.eta();
        jetPhi = j.phi();
        jetMass = j.mass();
        jetArea = j.jetArea();

        jetRawPt = j.correctedJet("Uncorrected").pt();
        jetRawMass = j.correctedJet("Uncorrected").mass();

        jetChargedHadronMult = j.chargedHadronMultiplicity();
        jetNeutralHadronMult = j.neutralHadronMultiplicity();
        jetChargedMult = j.chargedMultiplicity();
        jetNeutralMult = j.neutralMultiplicity();

        jetPtOrder = ptIdx;

        // Determine jet IDs
        jetLooseID = 0;
        jetTightID = 0;

        Float_t nhf = j.neutralHadronEnergyFraction();
        Float_t nemf = j.neutralEmEnergyFraction();
        Float_t chf = j.chargedHadronEnergyFraction();
        Float_t cemf = j.chargedEmEnergyFraction();
        unsigned int numconst = j.chargedMultiplicity() + j.neutralMultiplicity();
        unsigned int chm = j.chargedMultiplicity();

        if (abs(j.eta())<=2.7 && (numconst>1 && nhf<0.99 && nemf<0.99) && ((abs(j.eta())<=2.4 && chf>0 && chm>0 && cemf<0.99) || abs(j.eta())>2.4)) {
		  jetLooseID = 1;
		  if (nhf<0.90 && nemf<0.90) {
			jetTightID = 1;
		  }
        }

        // Add variables for deltaPhi and deltaEta for the two leading jets of the event
        dPhiJetsLO = deltaPhi(sortedJets[0].jet.phi(), sortedJets[1].jet.phi());
        dEtaJetsLO = sortedJets[0].jet.eta() - sortedJets[1].jet.eta();

        // The alpha variable is the third jet's pT divided by the average of the two leading jets' pT
        alpha = 0;
        // Make sure that there are at least 3 jets in the event
        if(sortedJets.size() > 2) {
		  Float_t leadingPtAvg = (sortedJets[0].jet.pt() + sortedJets[1].jet.pt()) * 0.5;
		  alpha = sortedJets[2].jet.pt() / leadingPtAvg;
        }

        // Assign flavors for each jet using three different flavor definitions
        partonFlav = j.partonFlavour();
        hadronFlav = j.hadronFlavour();

        physFlav = 0;
        if (j.genParton()) physFlav = j.genParton()->pdgId();

        // For convenience, save variables distinguishing gluon, light quark and other jets
        isPartonUDS = 0;
        isPartonG = 0;
        isPartonOther = 0;
        isPhysUDS = 0;
        isPhysG = 0;
        isPhysOther = 0;

        // Physics definition for flavors
        if(abs(physFlav) == 1 || abs(physFlav) == 2 || abs(physFlav) == 3) {
		  isPhysUDS = 1;
        } else if(abs(physFlav) == 21) {
		  isPhysG = 1;
        } else {
		  isPhysOther = 1;
        }

        // Parton definition for flavors
        if(abs(partonFlav) == 1 || abs(partonFlav) == 2 || abs(partonFlav) == 3) {
		  isPartonUDS = 1;
        } else if(abs(partonFlav) == 21) {
		  isPartonG = 1;
        } else {
		  isPartonOther = 1;
        }

        // Quark-gluon likelihood variables
        edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets, iJetRef));
        jetQGl = (*qglHandle)[jetRef];
        QG_ptD = (*ptDHandle)[jetRef];
        QG_axis2 = (*axis2Handle)[jetRef];
        QG_mult = (*multHandle)[jetRef];

        eventJetMult = selectedJets.size();

        // Loop over the PF candidates contained inside the jet, first sorting them in pT order
        std::vector<reco::CandidatePtr> pfCands = j.daughterPtrVector();
        std::sort(pfCands.begin(), pfCands.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt(); });
        int njetpf = 0;

        // Create a PF map for easier matching later
        std::map<const pat::PackedCandidate*, const pat::PackedCandidate> pfMap;

        // Here the jet girth is also calculated
        jetGirth = 0;

        unsigned int pfCandsSize = pfCands.size();
        for (unsigned int i = 0; i < pfCandsSize; ++i) {
		  const pat::PackedCandidate &pf = dynamic_cast<const pat::PackedCandidate &>(*pfCands[i]);
		  const pat::PackedCandidate* pfPointer = &pf;
		  pfMap.insert(std::pair <const pat::PackedCandidate*, const pat::PackedCandidate> (pfPointer, pf));

		  float dPhi = deltaPhi(pf.phi(), j.phi());
		  float dY = pf.rapidity() - j.rapidity();

		  jetGirth += sqrt(dY*dY + dPhi*dPhi) * pf.pt()/j.pt();
		  ++njetpf;
        }
        jetMult = njetpf;

        // Generator level jet variables and its constituents
        jetGenMatch = 0;
        genJetPt = 0;
        genJetEta = 0;
        genJetPhi = 0;
        genJetMass = 0;
        int ng = 0;

        // Check if the jet has a matching generator level jet
        if(j.genJet()) {
		  jetGenMatch = 1;

		  const reco::GenJet* gj = j.genJet();
		  genJetPt = gj->pt();
		  genJetEta = gj->eta();
		  genJetPhi = gj->phi();
		  genJetMass = gj->mass();

		  // Loop over the genjet's constituents
		  std::vector<const pat::PackedGenParticle*> genParticles;
		  for (unsigned int i = 0; i < gj->numberOfDaughters(); ++i) {
			const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(gj->daughter(i));
			genParticles.push_back( genParticle );
		  }

		  // Sort the constituents in pT order
		  std::sort(genParticles.begin(), genParticles.end(), [](const pat::PackedGenParticle* p1, const pat::PackedGenParticle* p2) {return p1->pt() > p2->pt(); });

		  unsigned int genParticlesSize = genParticles.size();
		  for (unsigned int i = 0; i != genParticlesSize; ++i) {
			const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(genParticles[i]);

			float dEta = (genParticle->eta()-gj->eta());
			float dPhi = deltaPhi(genParticle->phi(), gj->phi());

			genJetPF_pT[ng] = genParticle->pt();
			genJetPF_dR[ng] = deltaR(gj->eta(), gj->phi(), genParticle->eta(), genParticle->phi());
			genJetPF_dTheta[ng] = std::atan2(dPhi, dEta);
			genJetPF_mass[ng] = genParticle->mass();
			genJetPF_id[ng] = genParticle->pdgId();
			++ng;
		  }
		  nGenJetPF = ng;
        }

        // Loop over all the PF candidates in an event and save those which are
        //  within the area of |deltaEta| < 1 & |deltaPhi| < 1 from the center of the jet
        if (!(kMaxPF < pfs->size()))
		  assert(kMaxPF > pfs->size());
        int npfs = 0;

        unsigned int pfsSize = pfs->size();
        for (unsigned int i = 0; i != pfsSize; ++i) {
		  const pat::PackedCandidate &pf = (*pfs)[i];
		  const pat::PackedCandidate* pfPointer = &pf;

		  // Check if the PF was contained in the AK4 jet
		  if (pfMap.count(pfPointer)) {
			PF_fromAK4Jet[npfs] = 1;
		  } else {
			PF_fromAK4Jet[npfs] = 0;
		  }

		  float dEta = (pf.eta()-j.eta());
		  float dPhi = deltaPhi(pf.phi(),j.phi());

		  // Only save the PF candidates within the desired area
		  if ( (fabs(dEta) > 1.0) || (fabs(dPhi) > 1.0) ) continue;
		  PF_pT[npfs] = pf.pt();
		  PF_dR[npfs] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
		  PF_dTheta[npfs] = std::atan2(dPhi, dEta);
		  PF_dPhi[npfs] = dPhi;
		  PF_dEta[npfs] = dEta;
		  PF_mass[npfs] = pf.mass();
		  PF_id[npfs] = pf.pdgId();
		  PF_fromPV[npfs] =  pf.fromPV();
		  ++npfs;
        }
        nPF = npfs;

        // Save the jet in the tree
        jetTree->Fill();
      }
    } // saveJetTree

    // Add event information to the jet-based tree
    event = iEvent.id().event();
    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();

    // MC variables
    pthat = -1;
    if (genEventInfo->hasBinningValues()) {
      pthat = genEventInfo->binningValues()[0];
    }
    eventWeight = genEventInfo->weight();

    // Pileup info
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    Pileup_nTrueInt = -1;
    Pileup_nPU = -1;
    for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) {
		Pileup_nTrueInt = PVI->getTrueNumInteractions();
		Pileup_nPU = PVI->getPU_NumInteractions();
		continue;
      }
    }

    // Number of good primary vertices
    int vtxGood = 0;
    for(VertexCollection::const_iterator i_vtx = vertices->begin(); i_vtx != vertices->end(); i_vtx++) {
      if (!(i_vtx->isFake()) && i_vtx->ndof() > goodVtxNdof && fabs(i_vtx->z()) <= goodVtxZ && fabs(i_vtx->position().Rho()) <= goodVtxRho) {
		vtxGood++;
      }
    }
    PV_npvsGood = vtxGood;

    // Generator particle vertex info
    if (!(kMaxVtx < genParts->size())) {
      assert(kMaxVtx > genParts->size());
	}
    std::vector<math::XYZPoint> vertices_vec;
    std::vector<double> sumPt_vec;
    std::vector<double> sumPt2_vec;
    for(const auto& iPart : *genParts) {
      const math::XYZPoint& vertex = iPart.vertex();
      auto it = std::find(vertices_vec.begin(), vertices_vec.end(), vertex);
      if (it == vertices_vec.end()) {
		vertices_vec.push_back(vertex);
		sumPt_vec.push_back(0.0);
		sumPt2_vec.push_back(0.0);
		it = std::prev(vertices_vec.end());
      }
      else {
		int idx = std::distance(vertices_vec.begin(),it);
		sumPt_vec[idx] += iPart.pt();
		sumPt2_vec[idx] += std::pow(iPart.pt(), 2.);
      }
    }
    int ngenvtxs = 0;
    for(const auto& vtx : vertices_vec) {
      GenVtx_x[ngenvtxs] = vtx.x();
      GenVtx_y[ngenvtxs] = vtx.y();
      GenVtx_z[ngenvtxs] = vtx.z();
      GenVtx_sumPt[ngenvtxs] = sumPt_vec[ngenvtxs];
      GenVtx_sumPt2[ngenvtxs] = sumPt2_vec[ngenvtxs];
      ++ngenvtxs;
    }
    nGenVtx = ngenvtxs;

    // Rhos
    rhoAll = *pfRhoAllHandle;
    rhoCentral = *pfRhoCentralHandle;
    rhoCentralNeutral = *pfRhoCentralNeutralHandle;
    rhoCentralChargedPileUp = *pfRhoCentralChargedPileUpHandle;

    // Loop over all of the primary vertices in an event and save the pertinant information
    const reco::Vertex* primaryVertex = nullptr;
    if (!(kMaxVtx < pfs->size())) {
      assert(kMaxVtx > pfs->size());
	}
    int nvtxs = 0;
    for (size_t i=0; i<vertices->size();i++) {
      const reco::Vertex& vertex = (*vertices)[i];

      Vtx_x[nvtxs] = vertex.x();
      Vtx_y[nvtxs] = vertex.y();
      Vtx_z[nvtxs] = vertex.z();
      Vtx_xError[nvtxs] = vertex.xError();
      Vtx_yError[nvtxs] = vertex.yError();
      Vtx_zError[nvtxs] = vertex.zError();
      Vtx_ndof[nvtxs] = vertex.ndof();
      Vtx_chi2[nvtxs] = vertex.chi2();
      Vtx_ntrks[nvtxs] = vertex.nTracks();
      Vtx_isValid[nvtxs] = vertex.isValid();
      Vtx_isFake[nvtxs] = vertex.isFake();
      if (!(vertex.isFake()) && vertex.ndof() > goodVtxNdof && fabs(vertex.z()) <= goodVtxZ && fabs(vertex.position().Rho()) <= goodVtxRho) {
		Vtx_isGood[nvtxs] = 1;
		if (primaryVertex == nullptr) {
		  primaryVertex = &vertices->at(i);
		}
      }
      ++nvtxs;
    }
    nVtx = nvtxs;

    // Loop over all the PF candidates in an event and save the tracks
    if (!(kMaxTrk < pfs->size() + extraTrks->size())) {
      assert(kMaxTrk > (pfs->size() + extraTrks->size()));
	}
    int ntrks = 0;

    unsigned int pfsSize = pfs->size();
    for (unsigned int i = 0; i != pfsSize; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
      const reco::Track &trk = pf.pseudoTrack();

      // Save the track information
      Trk_charge[ntrks] = trk.charge();
      Trk_pT[ntrks] = trk.pt();
      Trk_eta[ntrks] = trk.eta();
      Trk_phi[ntrks] = trk.phi();
      Trk_px[ntrks] = trk.px();
      Trk_py[ntrks] = trk.py();
      Trk_pz[ntrks] = trk.pz();
      Trk_x[ntrks] = trk.vx();;
      Trk_y[ntrks] = trk.vy();
      Trk_z[ntrks] = trk.vz();
      if (primaryVertex) {
		const auto& primaryVertexPos = primaryVertex->position();
		Trk_dxypv[ntrks] = trk.dxy(primaryVertexPos);
		Trk_dzpv[ntrks] = trk.dz(primaryVertexPos);
      }
      else {
		Trk_dxypv[ntrks] = -100;
		Trk_dzpv[ntrks] = -100;
      }
      Trk_dxyerrorpv[ntrks] = trk.dxyError();
      Trk_dzerrorpv[ntrks] = trk.dzError();
      Trk_isLost[ntrks] = 0;
      ++ntrks;
    }

    // Loop over all the PF candidates in an event and save the tracks
    unsigned int extraTrksSize = extraTrks->size();
    for (unsigned int i = 0; i != extraTrksSize; ++i) {
      const pat::PackedCandidate &pf = (*extraTrks)[i];
      const reco::Track &trk = pf.pseudoTrack();

      // Save the track information
      Trk_charge[ntrks] = trk.charge();
      Trk_pT[ntrks] = trk.pt();
      Trk_eta[ntrks] = trk.eta();
      Trk_phi[ntrks] = trk.phi();
      Trk_px[ntrks] = trk.px();
      Trk_py[ntrks] = trk.py();
      Trk_pz[ntrks] = trk.pz();
      Trk_x[ntrks] = trk.vx();;
      Trk_y[ntrks] = trk.vy();
      Trk_z[ntrks] = trk.vz();
      if (primaryVertex) {
		const auto& primaryVertexPos = primaryVertex->position();
		Trk_dxypv[ntrks] = trk.dxy(primaryVertexPos);
		Trk_dzpv[ntrks] = trk.dz(primaryVertexPos);
      }
      else {
		Trk_dxypv[ntrks] = -100;
		Trk_dzpv[ntrks] = -100;
      }
      Trk_dxyerrorpv[ntrks] = trk.dxyError();
      Trk_dzerrorpv[ntrks] = trk.dzError();
      Trk_isLost[ntrks] = 1;
      ++ntrks;
    }
    nTrk = ntrks;

    // Save the tracks tree
    evtTree->Fill();

}

// Define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
