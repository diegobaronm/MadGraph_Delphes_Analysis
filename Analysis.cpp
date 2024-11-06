#include "Kinematics.h"

void PartonAnalysis(const std::string& fileName, const std::string& outputFileName, bool doWeights = false){
    ROOT::RDataFrame df("mytree", fileName);
    ROOT::RDF::RNode node = df; 
    
    // 1) Build the electron container
    // Status = 1 means final state particles
    auto myFilter = [](std::vector<int> pdgid, std::vector<int> status) { 
        std::vector<int> indices{};
        for (int i{} ; i<pdgid.size() ; i++){ 
            if (abs(pdgid.at(i)) == 11 && status.at(i) == 1) indices.push_back(i); 
        } 
        return indices; 
    };

    // Define the dummy column for the weight
    std::string weightColumn = "eventWeight";
    if (!doWeights){
        node = node.Define("dummyWeight",[](const std::vector<int>& pdgid) { return 1.0; },{"pdgID"});
        weightColumn = "dummyWeight";
    }

    node = node.Define("el_indices",myFilter,{"pdgID","pdgStatus"});

    auto getSize = [](std::vector<int> vec) { return vec.size(); };
    node = node.Define("el_indices_size",getSize,{"el_indices"});

    auto buildParticlesFromIndices = [](std::vector<int>& indices, std::vector<float>& px, std::vector<float>& py, std::vector<float>& pz, std::vector<float>& E){
        std::vector<TLorentzVector> particles;
        for (auto i : indices){
            TLorentzVector p;
            p.SetPxPyPzE(px.at(i), py.at(i), pz.at(i), E.at(i));
            particles.push_back(p);
        }

        // Order by pt
        std::sort(particles.begin(), particles.end(), [](const TLorentzVector& a, const TLorentzVector& b) { return a.Pt() > b.Pt(); });

        return particles;
    };
    node = node.Define("electrons",buildParticlesFromIndices,{"el_indices","px","py","pz","E"});

    // 2) Build jets container
    auto myFilterJets = [](std::vector<int> pdgid, std::vector<int> status) { 
        std::vector<int> indices{};
        for (int i{} ; i<pdgid.size() ; i++){ 
            if ( abs(pdgid.at(i)) < 7 && status.at(i) == 1) indices.push_back(i); 
        } 
        return indices; 
    };

    node = node.Define("jet_indices",myFilterJets,{"pdgID","pdgStatus"});
    node = node.Define("jet_indices_size",getSize,{"jet_indices"});
    node = node.Define("jets",buildParticlesFromIndices,{"jet_indices","px","py","pz","E"});

    // 3) Define jets and leptons pt
    node = node.Define("Electron_1_pt","electrons.at(0).Pt()");
    node = node.Define("Electron_2_pt","electrons.at(1).Pt()");
    node = node.Define("Jet_1_pt","jets.at(0).Pt()");
    node = node.Define("Jet_2_pt","jets.at(1).Pt()");

    // Some more kinematics
    // 4)
    auto deltaRapidity = [](std::vector<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    // 5)
    auto invariantMassTwoLeading = [](std::vector<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    // 6)
    auto dilepCentrality = [](std::vector<TLorentzVector> leptons, std::vector<TLorentzVector> jets) { 
        float lepton_xi = (leptons.at(0) + leptons.at(1)).Rapidity();
        float dijet_xi = jets.at(0).Rapidity() + jets.at(1).Rapidity();
        float delta_y = abs(jets.at(0).Rapidity() - jets.at(1).Rapidity());
        float centrality = abs((lepton_xi - 0.5*dijet_xi)/delta_y);
        return centrality;
    };
    // 7)
    auto dileptonPt = [](std::vector<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };
    // 8)
    auto deltaPhi = [](std::vector<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };

    // Create histogram container
    std::vector histograms = {node.Histo1D("el_indices", weightColumn)};
    histograms.push_back(node.Histo1D("jet_indices", weightColumn));
    histograms.push_back(node.Histo1D("el_indices_size", weightColumn));
    histograms.push_back(node.Histo1D("jet_indices_size", weightColumn));
    histograms.push_back(node.Histo1D({"Electron_1_pt_basic", "Electron_1_pt_basic", 500, 0, 500},"Electron_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Electron_2_pt_basic", "Electron_2_pt_basic", 500, 0, 500},"Electron_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Jet_1_pt_basic", "Jet_1_pt_basic", 1000, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Jet_2_pt_basic", "Jet_2_pt_basic", 1000, 0, 1000},"Jet_2_pt", weightColumn)); 

    
    node = node.Define("deltaPhi_electrons",deltaPhi,{"electrons"});
    node = node.Define("deltaPhi_jets",deltaPhi,{"jets"});
    node = node.Define("deltaRapidity",deltaRapidity,{"jets"});
    node = node.Define ("mjj",invariantMassTwoLeading,{"jets"});
    node = node.Define ("mll",invariantMassTwoLeading,{"electrons"});
    node = node.Define("dilepCentrality",dilepCentrality,{"electrons","jets"});
    node = node.Define("dileptonPt",dileptonPt,{"electrons"});

    histograms.push_back(node.Histo1D({"deltaPhi_electrons_basic", "deltaPhi_electrons_basic", 64, 0, 3.2},"deltaPhi_electrons", weightColumn));
    histograms.push_back(node.Histo1D({"deltaPhi_jets_basic", "deltaPhi_jets_basic", 64, 0, 3.2},"deltaPhi_jets", weightColumn));
    histograms.push_back(node.Histo1D({"deltaRapidity_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mjj_basic", "mjj_basic", 5000, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"mll_basic", "mll_basic", 1000, 0, 1000},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"dilepCentrality_basic", "dilepCentrality_basic", 350,0,3.5},"dilepCentrality", weightColumn));
    histograms.push_back(node.Histo1D({"dileptonPt_basic", "dileptonPt_basic", 1000,0,1000},"dileptonPt", weightColumn));


    // Save histograms
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    output->mkdir("histograms");
    output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}

void DelphesAnalysis(const std::string& fileName, const std::string& outputFileName, bool doWeights = false){
    ROOT::RDataFrame df("Delphes", fileName);
    ROOT::RDF::RNode node = df; 

    // Print column names
    g_LOG(LogLevel::INFO, "Columns available:");
    std::vector<std::string> columnNames = node.GetColumnNames();
    for (auto name : columnNames){
        std::cout << name << std::endl;
    }

    // Define the dummy column for the weight
    std::string weightColumn = "EvtW";
    float nEvents = 50000;
    float luminosity = 139.0;
    if (!doWeights){
        node = node.Define(weightColumn,[](const std::vector<float>& dummyV) { return 1.0; },{"Event.Weight"});
    } else {
        node = node.Define("SampleWeight",[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight) { return eventWeight[10] * 1000 * luminosity / (nEvents * eventWeight[10]); },{"Event.Weight"});
        node = node.Define(weightColumn,[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight, float sampleWeight) { return eventWeight[0] * sampleWeight; },{"Event.Weight","SampleWeight"});
    }


    // Define functions
    // 1)
    auto getSize = [](ROOT::VecOps::RVec<Float_t> vec) { return vec.size(); };
    // 2)
    auto buildParticlesZeroMass = [](ROOT::VecOps::RVec<Float_t>& pt){
        ROOT::VecOps::RVec<Float_t> m;
        for (std::size_t i{}; i < pt.size(); i++){
            m.push_back(0.0f);
        }
        return m;
    };
    // 3)
    auto buildParticlesFromContainer = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<Float_t>& eta, ROOT::VecOps::RVec<Float_t>& phi, ROOT::VecOps::RVec<Float_t>& m){
        ROOT::VecOps::RVec<TLorentzVector> particles;
        for (std::size_t i{}; i < pt.size(); i++){
            TLorentzVector p;
            p.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i));
            particles.push_back(p);
        }
        return particles;
    };
    // 4)
    auto deltaPhiMuons = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };
    // 5) 
    auto isPtOrdered = [](ROOT::VecOps::RVec<Float_t> pt) { 
        if (pt.size() < 2) g_LOG(LogLevel::ERROR, "Function called with less than two particles to compare!");
        return pt.at(0) > pt.at(1); 
    };
    // 6)
    auto numberOfHeavyJets = [](ROOT::VecOps::RVec<unsigned int> btagFlagVec) { 
        int numberBJets = 0;
        for (auto btag : btagFlagVec){
            if (btag == 1) numberBJets++;
        }
        return numberBJets;
    };
    // 7)
    auto deltaRapidity = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    // 8)
    auto invariantMassTwoLeading = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    // 9)
    auto isJetInGap = [](const TLorentzVector& j1,const TLorentzVector& j2, const TLorentzVector& testJet) {
        float delta_y_j1j2 = abs(j1.Rapidity()-j2.Rapidity());
        float delta_y_j1test = abs(j1.Rapidity()-testJet.Rapidity());
        float delta_y_j2test = abs(j2.Rapidity()-testJet.Rapidity());
        return (delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2);
    };

    auto nJetsInGap = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> jets) { 
        int gapJets = 0;
        TLorentzVector j1 = jets.at(0);
        TLorentzVector j2 = jets.at(1);

        if (jets.size() == 2) return 0;

        for (std::size_t i{2}; i < jets.size(); i++){
            TLorentzVector testJet = jets.at(i);
            if (isJetInGap(j1,j2,testJet)) gapJets++;
        }
        return gapJets;
    };
    // 9) 
    auto ptBalance = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float scalarPtSum{0};
        TLorentzVector vectorPtSum;
        for (std::size_t i{}; i < leptons.size(); i++){
            if ( i > 2) break;
            scalarPtSum += leptons.at(i).Pt();
            vectorPtSum += leptons.at(i);
        }

        // Jets
        for (std::size_t i{}; i < jets.size(); i++){

            if ( i>1 ){
                if ( isJetInGap(jets.at(0),jets.at(1),jets.at(i)) ){
                    scalarPtSum += jets.at(i).Pt();
                    vectorPtSum += jets.at(i);
                }
            } else {
                scalarPtSum += jets.at(i).Pt();
                vectorPtSum += jets.at(i);
            }
        }
        return vectorPtSum.Pt()/scalarPtSum;
    };
    // 10)
    auto dilepCentrality = [](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float lepton_xi = (leptons.at(0) + leptons.at(1)).Rapidity();
        float dijet_xi = jets.at(0).Rapidity() + jets.at(1).Rapidity();
        float delta_y = abs(jets.at(0).Rapidity() - jets.at(1).Rapidity());
        float centrality = abs((lepton_xi - 0.5*dijet_xi)/delta_y);
        return centrality;
    };
    // 11)
    auto dileptonPt = [](ROOT::VecOps::RVec<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };


    // Define new columns
    node = node.Define("Electron_mass",buildParticlesZeroMass,{"Electron.PT"});
    node = node.Define("electrons_p4",buildParticlesFromContainer, {"Electron.PT", "Electron.Eta", "Electron.Phi", "Electron_mass"});
    node = node.Define("jets_p4",buildParticlesFromContainer,{"Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass"});    
    node = node.Define("deltaPhi_electrons",deltaPhiMuons,{"electrons_p4"});
    node = node.Define("ElectronIsPtOrdered",isPtOrdered,{"Electron.PT"});
    node = node.Define("JetisPtOrdered",isPtOrdered,{"Jet.PT"});
    node = node.Define("Electron_1_pt","Electron.PT[0]");
    node = node.Define("Electron_2_pt","Electron.PT[1]");
    node = node.Define("Jet_1_pt","Jet.PT[0]");
    node = node.Define("Jet_2_pt","Jet.PT[1]");
    node = node.Define("Bjet_size",numberOfHeavyJets,{"Jet.BTag"});
    node = node.Define("deltaRapidity",deltaRapidity,{"jets_p4"});
    node = node.Define("mjj",invariantMassTwoLeading,{"jets_p4"});
    node = node.Define("mll",invariantMassTwoLeading,{"electrons_p4"});
    node = node.Define("nJetsInGap",nJetsInGap,{"jets_p4"});
    node = node.Define("ptBalance",ptBalance,{"electrons_p4","jets_p4"});
    node = node.Define("dilepCentrality",dilepCentrality,{"electrons_p4","jets_p4"});
    node = node.Define("dileptonPt",dileptonPt,{"electrons_p4"});
    
    
    // Filters
    node = node.Filter("Electron_size == 2");
    node = node.Filter("Muon_size == 0");
    node = node.Filter("Jet_size == 2 || Jet_size == 3");
    node = node.Filter("Electron.Charge[0] != Electron.Charge[1]");

    node = node.Filter("Electron.PT[0] >= 27");
    node = node.Filter("Electron.PT[1] >= 27");
    node = node.Filter("Jet.PT[0] >= 25");
    node = node.Filter("Jet.PT[1] >= 25");
    node = node.Filter("Jet_size == 3 ? Jet.PT[2] >= 25 : 1");
    node = node.Filter("mjj >= 500");
    node = node.Filter("mll > 40");

    std::vector histograms = {node.Histo1D("ElectronIsPtOrdered")};
    histograms.push_back(node.Histo1D({"lep1_pt_basic", "Electron_1_pt_basic", 500, 0, 500},"Electron_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"lep2_pt_basic", "Electron_2_pt_basic", 500, 0, 500},"Electron_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt_basic", "Jet_1_pt_basic", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt_basic", "Jet_2_pt_basic", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi_basic", "deltaPhi_electrons_basic", 64, 0, 3.2},"deltaPhi_electrons", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj_basic", "mjj_basic", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass_basic", "mll_basic", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"n_jets_interval_basic", "nJetsInGap_basic", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"pt_bal_basic", "ptBalance_basic", 100, 0, 1},"ptBalance", weightColumn));
    histograms.push_back(node.Histo1D({"Z_centrality_basic", "dilepCentrality_basic", 100, 0, 1},"dilepCentrality", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco_basic", "dileptonPt_basic", 100, 0, 1000},"dileptonPt", weightColumn));

    node = node.Filter("Bjet_size == 0");
    node = node.Filter("Electron.PT[0] >= 50");
    node = node.Filter("Electron.PT[1] >= 40");
    node = node.Filter("Jet.PT[0] >= 75");
    node = node.Filter("Jet.PT[1] >= 70");
    node = node.Filter("deltaRapidity >= 2.0");
    node = node.Filter("mjj >= 1000");
    node = node.Filter("mll > 101 && mll > 81");
    node = node.Filter("nJetsInGap == 0");
    node = node.Filter("ptBalance <= 0.15");
    node = node.Filter("dilepCentrality < 0.5");


    // More Histograms
    histograms.push_back(node.Histo1D("JetisPtOrdered", weightColumn));
    histograms.push_back(node.Histo1D("Electron_size", weightColumn));
    histograms.push_back(node.Histo1D("Jet_size", weightColumn));
    histograms.push_back(node.Histo1D("Bjet_size", weightColumn));
    histograms.push_back(node.Histo1D({"lep1_pt", "lep1_pt", 500, 0, 500},"Electron_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"lep2_pt", "lep2_pt", 500, 0, 500},"Electron_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt", "ljet0_pt", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt", "ljet1_pt", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi", "delta_phi", 64, 0, 3.2},"deltaPhi_electrons", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y", "delta_y", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj", "mass_jj", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass", "inv_mass", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"nJetsInGap", "nJetsInGap", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"ptBalance", "ptBalance", 100, 0, 1},"ptBalance", weightColumn));
    histograms.push_back(node.Histo1D({"Z_centrality", "Z_centrality", 100, 0, 1},"dilepCentrality", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco_basic_all", "Z_pt_reco_basic_all", 100, 0, 1000},"dileptonPt", weightColumn));


    // Save histograms
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    //output->mkdir("histograms");
    //output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}

void DelphesAnalysisTau(const std::string& fileName, const std::string& outputFileName, bool doWeights = false){
    ROOT::RDataFrame df("Delphes", fileName);
    ROOT::RDF::RNode node = df; 

    // Print column names
    g_LOG(LogLevel::INFO, "Columns available:");
    std::vector<std::string> columnNames = node.GetColumnNames();
    for (auto name : columnNames){
        std::cout << name << std::endl;
    }

    // Define the dummy column for the weight
    std::string weightColumn = "EvtW";
    float nEvents = 50000;
    float luminosity = 139.0;
    if (!doWeights){
        node = node.Define(weightColumn,[](const std::vector<float>& dummyV) { return 1.0; },{"Event.Weight"});
    } else {
        node = node.Define("SampleWeight",[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight) { return eventWeight[10] * 1000 * luminosity / (nEvents * eventWeight[10]); },{"Event.Weight"});
        node = node.Define(weightColumn,[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight, float sampleWeight) { return eventWeight[0] * sampleWeight; },{"Event.Weight","SampleWeight"});
    }


    // Define functions
    // 1)
    auto getSize = [](ROOT::VecOps::RVec<Float_t> vec) { return vec.size(); };
    auto getTauSize = [](ROOT::VecOps::RVec<Float_t> jetVec, ROOT::VecOps::RVec<unsigned int> tauTag) { 
        int tauSize = 0;
        for (std::size_t i{}; i < jetVec.size(); i++){
            if (tauTag.at(i) == 1) tauSize++;
        }
        return tauSize;
    };
    auto getNOTauSize = [](ROOT::VecOps::RVec<Float_t> jetVec, ROOT::VecOps::RVec<unsigned int> tauTag) { 
        int tauSize = 0;
        for (std::size_t i{}; i < jetVec.size(); i++){
            if (tauTag.at(i) != 1) tauSize++;
        }
        return tauSize;
    };
    // 2)
    auto buildParticlesZeroMass = [](ROOT::VecOps::RVec<Float_t>& pt){
        ROOT::VecOps::RVec<Float_t> m;
        for (std::size_t i{}; i < pt.size(); i++){
            m.push_back(0.0f);
        }
        return m;
    };
    // 3) Build particles functions
    auto buildParticlesFromContainer = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<Float_t>& eta, ROOT::VecOps::RVec<Float_t>& phi, ROOT::VecOps::RVec<Float_t>& m){
        ROOT::VecOps::RVec<TLorentzVector> particles;
        float minPt = 27.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            TLorentzVector p;
            p.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i));
            particles.push_back(p);
        }
        return particles;
    };
    auto buildTausFromContainer = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<Float_t>& eta, ROOT::VecOps::RVec<Float_t>& phi, ROOT::VecOps::RVec<Float_t>& m, ROOT::VecOps::RVec<unsigned int>& tauTag){
        ROOT::VecOps::RVec<TLorentzVector> particles;
        float minPt = 25.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            if (tauTag.at(i) == 1){
                TLorentzVector p;
                p.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i));
                particles.push_back(p);
            }
        }
        return particles;
    };
    auto buildNOTausFromContainer = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<Float_t>& eta, ROOT::VecOps::RVec<Float_t>& phi, ROOT::VecOps::RVec<Float_t>& m, ROOT::VecOps::RVec<unsigned int>& tauTag){
        ROOT::VecOps::RVec<TLorentzVector> particles;
        float minPt = 25.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            if (tauTag.at(i) != 1){
                TLorentzVector p;
                p.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i));
                particles.push_back(p);
            }
        }
        return particles;
    };
    auto buildTauCharge = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<unsigned int>& tauTag, ROOT::VecOps::RVec<int>& charge){
        ROOT::VecOps::RVec<int> charges;
        float minPt = 25.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            if (tauTag.at(i) == 1){
                charges.push_back(charge.at(i));
            }
        }
        return charges;
    };
    auto buildLeptonCharge = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<int>& charge){
        ROOT::VecOps::RVec<int> charges;
        float minPt = 27.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            charges.push_back(charge.at(i));
        }
        return charges;
    };


    // 4)
    auto deltaPhi = [](ROOT::VecOps::RVec<TLorentzVector> taus, const TLorentzVector& lepton) { return Kinematics::del_phi(taus.at(0).Phi(),lepton.Phi()); };
    auto deltaPhiLepMET = [](const TLorentzVector& p1, ROOT::VecOps::RVec<float> metPhi) { return Kinematics::del_phi(p1.Phi(),metPhi.at(0)); };
    auto deltaPhiTauMET = [](ROOT::VecOps::RVec<TLorentzVector> taus, ROOT::VecOps::RVec<float> metPhi) { return Kinematics::del_phi(taus.at(0).Phi(),metPhi.at(0)); };
    // 5) 
    auto isPtOrdered = [](ROOT::VecOps::RVec<Float_t> pt) { 
        if (pt.size() < 2) g_LOG(LogLevel::ERROR, "Function called with less than two particles to compare!");
        return pt.at(0) > pt.at(1); 
    };
    // 6)
    auto numberOfHeavyJets = [](ROOT::VecOps::RVec<unsigned int> btagFlagVec, ROOT::VecOps::RVec<Float_t> pt, ROOT::VecOps::RVec<unsigned int>& tauTag) { 
        int numberBJets = 0;
        float minPt = 25.0;
        for (std::size_t i{}; i < pt.size(); i++){
            if (pt.at(i) < minPt) continue;
            if (tauTag.at(i) == 1) continue;
            if (btagFlagVec.at(i) == 1) numberBJets++;
        }
        return numberBJets;
    };
    // 7)
    auto deltaRapidity = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    // 8)
    auto invariantMassTwoLeading = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    auto invariantMass = [](const TLorentzVector& p1, ROOT::VecOps::RVec<TLorentzVector> taus) { return (sqrt( 2.0*p1.Dot(taus[0])) ); };
    // 9)
    auto isJetInGap = [](const TLorentzVector& j1,const TLorentzVector& j2, const TLorentzVector& testJet) {
        float delta_y_j1j2 = abs(j1.Rapidity()-j2.Rapidity());
        float delta_y_j1test = abs(j1.Rapidity()-testJet.Rapidity());
        float delta_y_j2test = abs(j2.Rapidity()-testJet.Rapidity());
        return (delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2);
    };

    auto nJetsInGap = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> jets) { 
        int gapJets = 0;
        TLorentzVector j1 = jets.at(0);
        TLorentzVector j2 = jets.at(1);

        if (jets.size() == 2) return 0;

        for (std::size_t i{2}; i < jets.size(); i++){
            TLorentzVector testJet = jets.at(i);
            if (isJetInGap(j1,j2,testJet)) gapJets++;
        }
        return gapJets;
    };
    // 9) 
    auto ptBalance = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float scalarPtSum{0};
        TLorentzVector vectorPtSum;
        for (std::size_t i{}; i < leptons.size(); i++){
            if ( i > 2) break;
            scalarPtSum += leptons.at(i).Pt();
            vectorPtSum += leptons.at(i);
        }

        // Jets
        for (std::size_t i{}; i < jets.size(); i++){

            if ( i>1 ){
                if ( isJetInGap(jets.at(0),jets.at(1),jets.at(i)) ){
                    scalarPtSum += jets.at(i).Pt();
                    vectorPtSum += jets.at(i);
                }
            } else {
                scalarPtSum += jets.at(i).Pt();
                vectorPtSum += jets.at(i);
            }
        }
        return vectorPtSum.Pt()/scalarPtSum;
    };
    // 10)
    auto dilepCentrality = [](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float lepton_xi = (leptons.at(0) + leptons.at(1)).Rapidity();
        float dijet_xi = jets.at(0).Rapidity() + jets.at(1).Rapidity();
        float delta_y = abs(jets.at(0).Rapidity() - jets.at(1).Rapidity());
        float centrality = abs((lepton_xi - 0.5*dijet_xi)/delta_y);
        return centrality;
    };
    // 11)
    auto dileptonPt = [](ROOT::VecOps::RVec<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };
    auto taulepPt = [](ROOT::VecOps::RVec<TLorentzVector> taus, const TLorentzVector& lep) { return (taus[0] + lep).Pt(); };
    

    std::cout << "Hello world!" << std::endl;
    // Define new columns
    // Electron
    node = node.Define("Electron_mass",buildParticlesZeroMass,{"Electron.PT"});
    node = node.Define("electrons_p4",buildParticlesFromContainer, {"Electron.PT", "Electron.Eta", "Electron.Phi", "Electron_mass"});
    node = node.Define("Electron_charge",buildLeptonCharge,{"Electron.PT","Electron.Charge"});
    //Muons
    node = node.Define("Muon_mass",buildParticlesZeroMass,{"Muon.PT"});
    node = node.Define("muons_p4",buildParticlesFromContainer, {"Muon.PT", "Muon.Eta", "Muon.Phi", "Muon_mass"});
    node = node.Define("Muon_charge",buildLeptonCharge,{"Muon.PT","Muon.Charge"});
    // Jets
    node = node.Define("jets_p4",buildNOTausFromContainer,{"Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass","Jet.TauTag"});
    node = node.Define("NOTauJet_size",getNOTauSize,{"Jet.PT","Jet.TauTag"});
    // Taus
    node = node.Define("Tau_size",getTauSize,{"Jet.PT","Jet.TauTag"});
    node = node.Define("taus_p4",buildTausFromContainer,{"Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass","Jet.TauTag"});
    node = node.Define("Tau_charge",buildTauCharge,{"Jet.PT","Jet.TauTag","Jet.Charge"});

    std::cout << "Hello world!" << std::endl;

    // Filters object multiplicity
    node = node.Filter("(electrons_p4.size() == 1 && muons_p4.size() == 0) || (electrons_p4.size() == 0 && muons_p4.size() == 1)");
    node = node.Filter("Tau_size >= 1");
    node = node.Filter("jets_p4.size() == 2 || jets_p4.size() == 3");

    // Define the lepton in the event
    node = node.Define("lepton_p4",[](ROOT::VecOps::RVec<TLorentzVector> electrons_p4, ROOT::VecOps::RVec<TLorentzVector> muons_p4) { 
        // Safety check
        if (electrons_p4.size() > 0 &&  muons_p4.size() > 0) throw std::runtime_error("More than one lepton in the event!");
        if (electrons_p4.size() == 1) return electrons_p4.at(0);
        if (muons_p4.size() == 1) return muons_p4.at(0);
        return TLorentzVector();
    },{"electrons_p4","muons_p4"});
    node = node.Define("lepton_charge",[](ROOT::VecOps::RVec<int> electron_charge, ROOT::VecOps::RVec<int> muon_charge) { 
        // Safety check
        if (electron_charge.size() > 0 &&  muon_charge.size() > 0) throw std::runtime_error("More than one lepton in the event!");
        if (electron_charge.size() == 1) return electron_charge.at(0);
        if (muon_charge.size() == 1) return muon_charge.at(0);
        return 0;
    },{"Electron_charge","Muon_charge"});


    node = node.Filter("lepton_charge != Tau_charge[0]");


    std::cout << "Hello world!" << std::endl;
    
    node = node.Define("deltaPhi",deltaPhi,{"taus_p4","lepton_p4"});
    node = node.Define("deltaPhi_lepMET",deltaPhiLepMET,{"lepton_p4","MissingET.Phi"});
    node = node.Define("deltaPhi_tauMET",deltaPhiTauMET,{"taus_p4","MissingET.Phi"});
    node = node.Define("Lepton_pt","lepton_p4.Pt()");
    node = node.Define("Tau_pt","taus_p4[0].Pt()");
    node = node.Define("Jet_1_pt","jets_p4[0].Pt()");
    node = node.Define("Jet_2_pt","jets_p4[1].Pt()");
    node = node.Define("Bjet_size",numberOfHeavyJets,{"Jet.BTag","Jet.PT","Jet.TauTag"});
    node = node.Define("deltaRapidity",deltaRapidity,{"jets_p4"});
    node = node.Define("mjj",invariantMassTwoLeading,{"jets_p4"});
    node = node.Define("mll",invariantMass,{"lepton_p4","taus_p4"});
    node = node.Define("nJetsInGap",nJetsInGap,{"jets_p4"});
    node = node.Define("topology", Kinematics::getTauTauTopology, {"deltaPhi_lepMET", "deltaPhi_tauMET", "deltaPhi"});
    node = node.Define("omega", Kinematics::getOmega, {"topology", "deltaPhi_lepMET", "deltaPhi_tauMET", "deltaPhi"});
    node = node.Define("tauLepPt",taulepPt,{"taus_p4","lepton_p4"});
    
    node = node.Filter("mjj >= 500");
    node = node.Filter("mll > 40");

    std::cout << "Hello world!" << std::endl;
    
    std::vector histograms = {node.Histo1D({"Lepton_pt_basic", "Lepton_pt_basic", 500, 0, 500},"Lepton_pt", weightColumn)};
    histograms.push_back(node.Histo1D({"Tau_pt_basic", "Tau_pt_basic", 500, 0, 500},"Tau_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt_basic", "Jet_1_pt_basic", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt_basic", "Jet_2_pt_basic", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi_basic", "delta_phi_basic", 64, 0, 3.2},"deltaPhi", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj_basic", "mjj_basic", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass_basic", "inv_mass_basic", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"n_jets_interval_basic", "nJetsInGap_basic", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco_basic", "Z_pt_reco_basic", 100, 0, 1000},"tauLepPt", weightColumn));
    histograms.push_back(node.Histo1D({"topology_basic", "topology_basic", 5, 0, 5},"topology", weightColumn));
    histograms.push_back(node.Histo1D({"Tau_size_basic", "Tau_size_basic", 4, 0, 4},"Tau_size", weightColumn));
    histograms.push_back(node.Histo1D({"Bjet_size_basic", "Bjet_size_basic", 4, 0, 4},"Bjet_size", weightColumn));

    node = node.Filter("topology == 1 || topology == 2 || topology == 3"); 

    // Define MET vector and tau p4
    node = node.Define("met_p4",[](ROOT::VecOps::RVec<float> metPt, ROOT::VecOps::RVec<float> metPhi) { 
        TLorentzVector met;
        met.SetPtEtaPhiM(metPt.at(0),0,metPhi.at(0),0);
        return met;
    },{"MissingET.MET","MissingET.Phi"});
    node = node.Define("tau_p4",[](ROOT::VecOps::RVec<TLorentzVector> taus) { return taus.at(0); },{"taus_p4"});

    // Neutrinos
    node = node.Define("tau_nu_p4",Kinematics::getTauNeutrino,{"topology","met_p4","tau_p4","lepton_p4"});
    node = node.Define("lep_nu_p4",Kinematics::getLepNeutrino,{"topology","met_p4","tau_p4","lepton_p4"});

    node = node.Define("visible_leptons_pack", [](const TLorentzVector& tau, const TLorentzVector& lepton){
        ROOT::VecOps::RVec<TLorentzVector> leptons;
        leptons.push_back(tau);
        leptons.push_back(lepton);
        return leptons;
    },{"tau_p4","lepton_p4"});
    node = node.Define("leptons_pack", [](const TLorentzVector& tau, const TLorentzVector& lepton, const TLorentzVector& tauNu, const TLorentzVector& lepNu){
        ROOT::VecOps::RVec<TLorentzVector> leptons;
        leptons.push_back(tau);
        leptons.push_back(lepton);
        leptons.push_back(tauNu);
        leptons.push_back(lepNu);
        return leptons;
    },{"tau_p4","lepton_p4","tau_nu_p4","lep_nu_p4"});


    
    node = node.Define("reco_mass",Kinematics::getRecoMass,{"topology","tau_p4","lepton_p4","tau_nu_p4","lep_nu_p4"});

    histograms.push_back(node.Histo1D({"omega_basic", "omega_basic", 60, -3, 3},"omega", weightColumn));
    histograms.push_back(node.Histo1D({"reco_mass_basic", "reco_mass_basic", 1000, 0, 1000},"reco_mass", weightColumn));
    std::cout << "Hello world!" << std::endl;


    node = node.Define("ptBalance",ptBalance,{"leptons_pack","jets_p4"});
    node = node.Define("dilepCentrality",dilepCentrality,{"visible_leptons_pack","jets_p4"});

    node = node.Filter("Bjet_size == 0");
    node = node.Filter("Jet_1_pt >= 75");
    node = node.Filter("Jet_2_pt >= 70");
    node = node.Filter("deltaRapidity >= 2.0");
    node = node.Filter("mjj >= 1000");
    node = node.Filter("reco_mass > 110");
    node = node.Filter("omega > -0.2 && omega < 1.6");
    node = node.Filter("nJetsInGap == 0");
    node = node.Filter("ptBalance <= 0.15");
    node = node.Filter("dilepCentrality < 0.5");

  
    // More Histograms
    histograms.push_back(node.Histo1D("NOTauJet_size", weightColumn));
    histograms.push_back(node.Histo1D({"Lepton_pt", "Lepton_pt", 500, 0, 500},"Lepton_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Tau_pt", "Tau_pt", 500, 0, 500},"Tau_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt", "Jet_1_pt", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt", "Jet_2_pt", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi", "delta_phi", 64, 0, 3.2},"deltaPhi", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y", "deltaRapidity", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj", "mjj", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass", "inv_mass", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"reco_mass", "reco_mass", 1000, 0, 1000},"reco_mass", weightColumn));
    histograms.push_back(node.Histo1D({"n_jets_interval", "nJetsInGap", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco", "Z_pt_reco", 100, 0, 1000},"tauLepPt", weightColumn));
    histograms.push_back(node.Histo1D({"topology", "topology", 5, 0, 5},"topology", weightColumn));
    histograms.push_back(node.Histo1D({"Tau_size", "Tau_size", 4, 0, 4},"Tau_size", weightColumn));
    histograms.push_back(node.Histo1D({"Bjet_size", "Bjet_size", 4, 0, 4},"Bjet_size", weightColumn));
    histograms.push_back(node.Histo1D({"ptBalance", "ptBalance", 100, 0, 1},"ptBalance", weightColumn));
    histograms.push_back(node.Histo1D({"dilepCentrality", "dilepCentrality", 100, 0, 1},"dilepCentrality", weightColumn));
    

    // Save histograms
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    //output->mkdir("histograms");
    //output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}


void Analysis(){
    // Example of parton - level analysis
    //PartonAnalysis("DATA/CompleteDiagrams.root", "PartonLevel_Complete_Diagrams.root", false);
    //PartonAnalysis("DATA/EWKDiagrams.root", "PartonLevel_EWK_Diagrams.root", false);
    //PartonAnalysis("DATA/NO_EWKDiagrams.root", "PartonLevel_NO_EWK_Diagrams.root", false);
    

    // Example of Delphes analysis with light leptons
    //DelphesAnalysis("DATA/EWK_Delphes.root", "EWK_Delphes.root", true);

    // Example of Delphes analysis with tau leptons
    DelphesAnalysisTau("DATA/Zp_Delphes.root", "Zp_Weighted.root", true);
    DelphesAnalysisTau("DATA/Zp_plus_SM_Delphes.root", "Zp_Interference_SM_Weighted.root", true);
    DelphesAnalysisTau("DATA/HM_SM_Delphes.root", "HM_SM_Weighted.root", true);
}