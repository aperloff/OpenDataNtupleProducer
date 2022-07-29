# OpenDataNtupleProducer - a tuple producer from CMS OpenData MiniAOD which saves jets, tracks, and vertices

This is a CMSSW module for producing mostly flat tuples from 13 TeV Run2 MC samples.

The code is intended to run inside the CMS Virtual Machine or Docker environments.

It is based on the CERN Open Data sample code here: https://opendata.cern.ch/record/12100

## Setting up

First install the CMS Virtual Machine as per instructions here: http://opendata.cern.ch/docs/cms-virtual-machine-2011

Alternatively you can use one of the Docker images located at: https://gitlab.cern.ch/cms-cloud/cmssw-docker/container_registry

Then open up the terminal and create the work area:
```
mkdir WorkingArea
cd WorkingArea
cmsrel CMSSW_8_0_26_patch2
cd CMSSW_8_0_26_patch2/src
cmsenv                  # activates the CMSSW environment
git clone https://github.com/aperloff/OpenDataNtupleProducer.git
scram b                 # compiles the code
cd OpenDataNtupleProducer/JetAnalyzer
```

## Running the code

One can change the conditions for the code execution (most notably, the amount of processed events) in the file ```python/ConfFile_cfg.py```.

To generate the ntuples, execute
```
cmsRun python/ConfFile_cfg.py
```

This configuration file runs the ```plugins/JetAnalyzer.cc``` code, which performs the data processing and extracts the jets, tracks, and vertices from collision events. There is an additional configuration file ```python/QGLikelihood_cfi.py```, which is there for the quark/gluon jet likelihood variables. It references the ```database/QGL_cmssw8020_v2.db``` file. There is no need to touch these QGL files.

## Data exploration
To quickly explore the created file, first load it up in ROOT
```
root -l OpenDataNtuple_RunIISummer16_13TeV_MC.root
```
Here you can create a ```TBrowser``` and examine the contents of the file and the distributions for each variable.
```
new TBrowser
```

## Content of the tuples

This program is capable of saving two different ROOT trees. One tree (`jetTree`) contains values relevant to a given jet and entries are saved on a jet-by-jet basis. The other tree (`evtTree`) contains track and vertex information and entries are saved on an event-by-event basis.

In `jetTree`, the reconstructed jets are AK4 jets clustered from Particle Flow candidates. The standard L1+L2+L3+residual energy corrections are applied to the jets and pileup is reduced using the CHS algorithm. The following table lists the branches available in that tree:

| Data variable | Type | Description |
| :---------------------- | -----------------: | :---------------------- |
| jetPt | Float_t | Transverse momentum of the jet |
| jetEta | Float_t | Pseudorapidity (η) of the jet |
| jetPhi | Float_t | Azimuthal angle (ϕ) of the jet |
| jetMass | Float_t | Mass of the jet |
| jetGirth | Float_t | Girth of the jet (as defined in arXiv:1106.3076 [hep-ph]) |
| jetArea | Float_t | Catchment area of the jet; used for jet energy corrections |
| jetRawPt | Float_t | Transverse momentum of the jet before the energy corrections |
| jetRawMass | Float_t | Mass of the jet before the energy corrections |
| jetLooseID | UInt_t | Binary variable indicating whether the jet passes 'loose' criteria for being a real jet |
| jetTightID | UInt_t | Binary variable indicating whether the jet passes 'tight' criteria for being a real jet |
| jetGenMatch | UInt_t | 1: if a matched generator level jet exists; 0: if no match was found |
| jetQGl | Float_t | Quark-Gluon jet likelihood discriminant variable built out of the three following variables (see the report CMS-PAS-JME-13-002 for more information) |
| QG_ptD | Float_t | Jet energy variable (see CMS-PAS-JME-13-002) |
| QG_axis2 | Float_t | Minor axis of the jet (see CMS-PAS-JME-13-002) |
| QG_mult | UInt_t | Jet constituent multiplicity with additional cuts (see CMS-PAS-JME-13-002) |
| partonFlav | Int_t | Flavour of the jet, as defined by the CMS parton-based definition |
| hadronFlav | Int_t | Flavour of the jet, as defined by the CMS hadron-based definition |
| physFlav | Int_t | Flavour of the jet, as defined by the CMS 'physics' definition (if in doubt, use this) |
| isPartonUDS | UInt_t | Indicates light quark (Up, Down, Strange) jets: partonFlav = 1, 2, 3 |
| isPartonG | UInt_t | Indicates gluon jets: partonFlav = 21 |
| isPartonOther | UInt_t | Indicates any other kind of jet: partonflav != 1, 2, 3, 21 |
| isPhysUDS | UInt_t | Indicates light quark (Up, Down, Strange) jets: physFlav = 1, 2, 3 |
| isPhysG | UInt_t | Indicates gluon jets: physFlav = 21 |
| isPhysOther | UInt_t | Indicates any other kind of jet: physFlav != 1, 2, 3, 21 |
| jetChargedHadronMult | UInt_t | Multiplicity of charged hadron jet constituents |
| jetNeutralHadronMult | UInt_t | Multiplicity of neutral hadron jet constituents |
| jetChargedMult | UInt_t | Multiplicity of charged jet constituents |
| jetNeutralMult | UInt_t | Multiplicity of neutral jet constituents |
| jetMult | UInt_t | Multiplicity of jet constituents |
| nPF | UInt_t | Number of particle flow (PF) candidates (particles reconstructed by the particle flow algorithm); contains all particles within \|Δϕ\| < 1 and \|Δη\| < 1 from the center of the jet |
| PF_pT[nPF] | Float_t | Transverse momentum of a PF candidate |
| PF_dR[nPF] | Float_t | Distance of a PF candidate to the center of the jet |
| PF_dTheta[nPF] | Float_t | Polar angle (θ) of a PF candidate |
| PF_dPhi[nPF] | Float_t | Azimuthal angle (ϕ) of a PF candidate |
| PF_dEta[nPF] | Float_t | Pseudorapidity (η) of a PF candidate |
| PF_mass[nPF] | Float_t | Mass of a PF candidate |
| PF_id[nPF] | Int_t | Generator level particle identifier for the particle flow candidates, as defined in the PDG particle numbering scheme |
| PF_fromPV[nPF] | UInt_t | A number indicating how tightly a particle is associated with the primary vertex (ranges from 3 to 0) |
| PF_fromAK4Jet[nPF] | UInt_t | 1: if the particle flow candidate is a constituent of the reconstructed AK4 jet; 0: if it is not a constituent of the jet |
| genJetPt | Float_t | Transverse momentum of the matched generator level jet |
| genJetEta | Float_t | Pseudorapidity (η) of the matched generator level jet |
| genJetPhi | Float_t | Azimuthal angle (ϕ) of the matched generator level jet |
| genJetMass | Float_t | Mass of the matched generator level jet |
| nGenJetPF | UInt_t | Number of particles in the matched generator level jet |
| genPF_pT[nGenJetPF] | Float_t | Transverse momentum of a particle in the matched generator level jet |
| genPF_dR[nGenJetPF] | Float_t | Distance of a particle to the center of the matched generator level jet |
| genPF_dTheta[nGenJetPF] | Float_t | Polar angle (θ) of a particle in the matched generator level jet |
| genPF_mass[nGenJetPF] | Float_t | Mass of a particle in the matched generator level jet |
| genPF_id[nGenJetPF] | Int_t | Generator level particle identifier for the particles in the matched generator level jet, as defined in the PDG particle numbering scheme |
| eventJetMult | UInt_t | Multiplicity of jets in the event |
| jetPtOrder | UInt_t | Indicates the ranking number of the jet, as the jets are ordered by their transverse momenta within a single event |
| dPhiJetsLO | Float_t | The phi difference of the two leading jets |
| dEtaJetsLO | Float_t | The eta difference of the two leading jets |
| alpha | Float_t | If there are at least 3 jets in the event, alpha is the third jet's transverse momentum divided by the average transverse momentum of the two leading jets |
| event | ULong64_t | Event number |
| run | UInt_t | Run number |
| lumi | UInt_t | Luminosity block |

The following table lists the branches available in `evtTree` tree:

| Data variable | Type | Description |
| :---------------------- | -----------------: | :---------------------- |
| event | ULong64_t | Event number |
| run | UInt_t | Run number |
| lumi | UInt_t | Luminosity block |
| pthat | Float_t | Transverse momentum of the generated hard process |
| eventWeight | Float_t | Monte Carlo generator weight |
| rhoAll | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates |
| rhoCentral | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates with |eta|<2.5 |
| rhoCentralNeutral | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all neutral PF candidates with |eta| < 2.5 |
| rhoCentralChargedPileUp | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all PF charged hadrons associated to pileup vertices and with |eta| < 2.5 |
| PV_npvsGood | UInt_t | The number of good reconstructed primary vertices; selection: !isFake && ndof > 4 && abs(z) <= 24 && position.Rho < 2 |
| Pileup_nPU | UInt_t | The number of pileup interactions that have been added to the event in the current bunch crossing |
| Pileup_nTrueInt | UInt_t | The true mean number of the poisson distribution for this event from which the number of interactions in each bunch crossing has been sampled |
| nTrk | UInt_t | Number of tracks; contains all tracks matched to particle flow candidates and the lost tracks (tracks not matched to particle flow candidates) |
| Trk_charge[nTrk] | Int_t | Charge of a track |
| Trk_pT[nTrk] | Float_t | Transverse momentum of a track [GeV] |
| Trk_eta[nTrk] | Float_t | Pseudorapidity of a track |
| Trk_phi[nTrk] | Float_t | Azimuthal angle of a track |
| Trk_px[nTrk] | Float_t | x coordinate of the track mometum vector [cm] |
| Trk_py[nTrk] | Float_t | y coordinate of the track mometum vector [cm] |
| Trk_pz[nTrk] | Float_t | z coordinate of the track mometum vector [cm] |
| Trk_x[nTrk] | Float_t | Position of the track along the x direction at its reference point [cm] |
| Trk_y[nTrk] | Float_t | Position of the track along the y direction at its reference point [cm] |
| Trk_z[nTrk] | Float_t | Position of the track along the beamline at its reference point [cm] |
| Trk_dxypv[nTrk] | Float_t | Transverse impact parameter of the track [cm] |
| Trk_dxyerrorpv[nTrk] | Float_t | Error on the transverse impact parameter of the track [cm] |
| Trk_dzpv[nTrk] | Float_t | Longitudinal impact parameter of the track [cm] |
| Trk_dzerrorpv[nTrk] | Float_t | Error on the longitudinal impact parameter of the track [cm] |
| Trk_isLost[nTrk] | UInt_t | Indicates if the track was part of the lostTracks collections (i.e. not matched to a PF candidate) |
| nGenVtx | UInt_t | Number of generator level primary vertices |
| GenVtx_x[nGenVtx] | Float_t | x coordinate of a generator level primary vertex [cm] |
| GenVtx_y[nGenVtx] | Float_t | y coordinate of a generator level primary vertex [cm] |
| GenVtx_z[nGenVtx] | Float_t | z coordinate of a generator level primary vertex [cm] |
| GenVtx_sumPt[nGenVtx] | Float_t | Sum of the pT of the generator particles associated to this vertex [GeV] |
| GenVtx_sumPt2[nGenVtx] | Float_t | Sum of the pT^2 of the generator particles associated to this vertex [GeV^2] |
| nVtx | UInt_t | Number of primary vertices |
| Vtx_x[nVtx] | Float_t | x coordinate of a primary vertex [cm] |
| Vtx_y[nVtx] | Float_t | y coordinate of a primary vertex [cm] |
| Vtx_z[nVtx] | Float_t | z coordinate of a primary vertex [cm] |
| Vtx_xError[nVtx] | Float_t | Error on the x positon of a primary vertex [cm] |
| Vtx_yError[nVtx] | Float_t | Error on the y positon of a primary vertex [cm] |
| Vtx_zError[nVtx] | Float_t | Error on the z positon of a primary vertex [cm] |
| Vtx_ndof[nVtx] | UInt_t | Number of degrees of freedom of the primary vertex fit |
| Vtx_chi2[nVtx] | Float_t | Chi square of the primary vertex fit |
| Vtx_ntrks[nVtx] | UInt_t | Number of tracks associated to a primary vertex |
| Vtx_isValid[nVtx] | UInt_t | Indicates if a primary vertex is valid |
| Vtx_isFake[nVtx] | UInt_t | Indicates if a primary vertex is fake |
| Vtx_isGood[nVtx] | UInt_t | Indicates if a primary vertex is good |