#include "Kinematics.h"

void PartonAnalysis(std::string fileName){
    ROOT::RDataFrame df("mytree", fileName);
    ROOT::RDF::RNode node = df; 
    
    auto myFilter = [](std::vector<int> pdgid) { std::vector<int> indices{}; for (int i{} ; i<pdgid.size() ; i++){ if (pdgid.at(i) == 11 || pdgid.at(i) ==-11) indices.push_back(i); } return indices; };
    node = node.Define("el_indices",myFilter,{"pdgID"});

    auto getSize = [](std::vector<int> vec) { return vec.size(); };
    node = node.Define("el_indices_size",getSize,{"el_indices"});

    auto buildParticlesFromIndices = [](std::vector<int>& indices, std::vector<float>& px, std::vector<float>& py, std::vector<float>& pz, std::vector<float>& E){
        std::vector<TLorentzVector> particles;
        for (auto i : indices){
            TLorentzVector p;
            p.SetPxPyPzE(px.at(i), py.at(i), pz.at(i), E.at(i));
            particles.push_back(p);
        }
        return particles;
    };
    node = node.Define("electrons",buildParticlesFromIndices,{"el_indices","px","py","pz","E"});

    auto deltaPhi = [](std::vector<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };
    node = node.Define("deltaPhi",deltaPhi,{"electrons"});

    auto h = node.Histo1D("el_indices");
    auto h2 = node.Histo1D("el_indices_size");
    auto h3 = node.Histo1D("deltaPhi");

    // Save histograms
    auto output = TFile::Open("full_diagrams_genlevelcuts_output.root", "RECREATE");
    output->cd();
    output->mkdir("histograms");
    output->cd("histograms");
    h->Write();
    h2->Write();
    h3->Write();
    output->Close();
}

void DelphesAnalysis(std::string fileName){
    ROOT::RDataFrame df("Delphes", fileName);
    ROOT::RDF::RNode node = df; 

    // Print column names
    g_LOG(LogLevel::INFO, "Columns available:");
    std::vector<std::string> columnNames = node.GetColumnNames();
    for (auto name : columnNames){
        std::cout << name << std::endl;
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
    auto isPtOrdered = [](ROOT::VecOps::RVec<Float_t> pt) { return pt.at(0) > pt.at(1); };
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

                


    // Define new columns
    node = node.Define("Muon_mass",buildParticlesZeroMass,{"Muon.PT"});
    node = node.Define("muons_p4",buildParticlesFromContainer, {"Muon.PT", "Muon.Eta", "Muon.Phi", "Muon_mass"});
    node = node.Define("jets_p4",buildParticlesFromContainer,{"Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass"});    
    node = node.Define("deltaPhi_muons",deltaPhiMuons,{"muons_p4"});
    node = node.Define("MuonisPtOrdered",isPtOrdered,{"Muon.PT"});
    node = node.Define("JetisPtOrdered",isPtOrdered,{"Jet.PT"});
    node = node.Define("Muon_1_pt","Muon.PT[0]");
    node = node.Define("Muon_2_pt","Muon.PT[1]");
    node = node.Define("Jet_1_pt","Jet.PT[0]");
    node = node.Define("Jet_2_pt","Jet.PT[1]");
    node = node.Define("Bjet_size",numberOfHeavyJets,{"Jet.BTag"});
    node = node.Define("deltaRapidity",deltaRapidity,{"jets_p4"});
    node = node.Define("mjj",invariantMassTwoLeading,{"jets_p4"});
    node = node.Define("mll",invariantMassTwoLeading,{"muons_p4"});
    node = node.Define("nJetsInGap",nJetsInGap,{"jets_p4"});
    node = node.Define("ptBalance",ptBalance,{"muons_p4","jets_p4"});
    node = node.Define("dilepCentrality",dilepCentrality,{"muons_p4","jets_p4"});
    
    // Filters
    node = node.Filter("Muon_size == 2");
    node = node.Filter("Electron_size == 0");
    node = node.Filter("Jet_size == 2 || Jet_size == 3");
    node = node.Filter("Muon.Charge[0] != Muon.Charge[1]");
    node = node.Filter("Bjet_size == 0");
    node = node.Filter("Muon.PT[0] >= 50");
    node = node.Filter("Muon.PT[1] >= 40");
    node = node.Filter("Jet.PT[0] >= 75");
    node = node.Filter("Jet.PT[1] >= 70");
    node = node.Filter("deltaRapidity >= 2.0");
    node = node.Filter("mjj >= 1000");
    node = node.Filter("mll < 101 && mll > 81");
    node = node.Filter("nJetsInGap == 0");
    node = node.Filter("ptBalance <= 0.15");


    // Histograms
    std::vector histograms = {node.Histo1D("MuonisPtOrdered")};
    histograms.push_back(node.Histo1D("JetisPtOrdered"));
    histograms.push_back(node.Histo1D("Muon_size"));
    histograms.push_back(node.Histo1D("Jet_size"));
    histograms.push_back(node.Histo1D("Bjet_size"));
    histograms.push_back(node.Histo1D({"Muon_1_pt", "Muon_1_pt", 300, 0, 300},"Muon_1_pt"));
    histograms.push_back(node.Histo1D({"Muon_2_pt", "Muon_2_pt", 300, 0, 300},"Muon_2_pt"));
    histograms.push_back(node.Histo1D({"Jet_1_pt", "Jet_1_pt", 300, 0, 300},"Jet_1_pt"));
    histograms.push_back(node.Histo1D({"Jet_2_pt", "Jet_2_pt", 300, 0, 300},"Jet_2_pt"));   
    histograms.push_back(node.Histo1D({"deltaPhi_muons", "deltaPhi_muons", 64, 0, 32},"deltaPhi_muons"));
    histograms.push_back(node.Histo1D({"deltaRapidity", "deltaRapidity", 100, 0, 10},"deltaRapidity"));
    histograms.push_back(node.Histo1D({"mjj", "mjj", 500, 0, 5000},"mjj"));
    histograms.push_back(node.Histo1D({"mll", "mll", 300, 0, 300},"mll"));
    histograms.push_back(node.Histo1D({"nJetsInGap", "nJetsInGap", 2, 0, 2},"nJetsInGap"));
    histograms.push_back(node.Histo1D({"ptBalance", "ptBalance", 100, 0, 1},"ptBalance"));
    histograms.push_back(node.Histo1D({"dilepCentrality", "dilepCentrality", 100, 0, 1},"dilepCentrality"));


    // Save histograms
    std::string outputFileName = "output_" + fileName;
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    output->mkdir("histograms");
    output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}


void Analysis(){
    DelphesAnalysis("full_diagrams_delphes.root");
}

/******************************************************************************
*Tree    :Delphes   : Analysis tree                                          *
*Entries :    50000 : Total =     10463008826 bytes  File  Size = 3699377770 *
*        :          : Tree compression factor =   2.83                       *
******************************************************************************
*Br    0 :Event     : Int_t Event_                                           *
*Entries :    50000 : Total  Size=     468763 bytes  File Size  =      80531 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.09     *
*............................................................................*
*Br    1 :Event.fUniqueID : UInt_t fUniqueID[Event_]                         *
*Entries :    50000 : Total  Size=     413420 bytes  File Size  =      81567 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.04     *
*............................................................................*
*Br    2 :Event.fBits : UInt_t fBits[Event_]                                 *
*Entries :    50000 : Total  Size=     412944 bytes  File Size  =      81107 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.06     *
*............................................................................*
*Br    3 :Event.Number : Long64_t Number[Event_]                             *
*Entries :    50000 : Total  Size=     613063 bytes  File Size  =     179313 *
*Baskets :      115 : Basket Size=      12800 bytes  Compression=   3.40     *
*............................................................................*
*Br    4 :Event.ReadTime : Float_t ReadTime[Event_]                          *
*Entries :    50000 : Total  Size=     413301 bytes  File Size  =     231747 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.77     *
*............................................................................*
*Br    5 :Event.ProcTime : Float_t ProcTime[Event_]                          *
*Entries :    50000 : Total  Size=     413301 bytes  File Size  =     228514 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.80     *
*............................................................................*
*Br    6 :Event.ProcessID : Int_t ProcessID[Event_]                          *
*Entries :    50000 : Total  Size=     413420 bytes  File Size  =      82486 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.98     *
*............................................................................*
*Br    7 :Event.MPI : Int_t MPI[Event_]                                      *
*Entries :    50000 : Total  Size=     412706 bytes  File Size  =      80880 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br    8 :Event.Weight : Float_t Weight[Event_]                              *
*Entries :    50000 : Total  Size=     413063 bytes  File Size  =      82369 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.98     *
*............................................................................*
*Br    9 :Event.CrossSection : Float_t CrossSection[Event_]                  *
*Entries :    50000 : Total  Size=     413777 bytes  File Size  =     260998 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.57     *
*............................................................................*
*Br   10 :Event.CrossSectionError : Float_t CrossSectionError[Event_]        *
*Entries :    50000 : Total  Size=     414372 bytes  File Size  =      84063 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.90     *
*............................................................................*
*Br   11 :Event.Scale : Float_t Scale[Event_]                                *
*Entries :    50000 : Total  Size=     412944 bytes  File Size  =      82140 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.99     *
*............................................................................*
*Br   12 :Event.AlphaQED : Float_t AlphaQED[Event_]                          *
*Entries :    50000 : Total  Size=     413301 bytes  File Size  =      82255 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.99     *
*............................................................................*
*Br   13 :Event.AlphaQCD : Float_t AlphaQCD[Event_]                          *
*Entries :    50000 : Total  Size=     413301 bytes  File Size  =      82255 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   4.99     *
*............................................................................*
*Br   14 :Event.ID1 : Int_t ID1[Event_]                                      *
*Entries :    50000 : Total  Size=     412706 bytes  File Size  =      80877 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   15 :Event.ID2 : Int_t ID2[Event_]                                      *
*Entries :    50000 : Total  Size=     412706 bytes  File Size  =      80877 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   16 :Event.X1  : Float_t X1[Event_]                                     *
*Entries :    50000 : Total  Size=     412587 bytes  File Size  =      80876 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   17 :Event.X2  : Float_t X2[Event_]                                     *
*Entries :    50000 : Total  Size=     412587 bytes  File Size  =      80876 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   18 :Event.ScalePDF : Float_t ScalePDF[Event_]                          *
*Entries :    50000 : Total  Size=     413301 bytes  File Size  =      81339 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.05     *
*............................................................................*
*Br   19 :Event.PDF1 : Float_t PDF1[Event_]                                  *
*Entries :    50000 : Total  Size=     412825 bytes  File Size  =      80879 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   20 :Event.PDF2 : Float_t PDF2[Event_]                                  *
*Entries :    50000 : Total  Size=     412825 bytes  File Size  =      80879 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.07     *
*............................................................................*
*Br   21 :Event_size : Event_size/I                                          *
*Entries :    50000 : Total  Size=     211825 bytes  File Size  =      13789 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  15.17     *
*............................................................................*
*Br   22 :Weight    : Int_t Weight_                                          *
*Entries :    50000 : Total  Size=     422976 bytes  File Size  =      81106 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.05     *
*............................................................................*
*Br   23 :Weight.fUniqueID : UInt_t fUniqueID[Weight_]                       *
*Entries :    50000 : Total  Size=     613541 bytes  File Size  =      76519 *
*Baskets :      115 : Basket Size=      12800 bytes  Compression=   7.98     *
*............................................................................*
*Br   24 :Weight.fBits : UInt_t fBits[Weight_]                               *
*Entries :    50000 : Total  Size=     613065 bytes  File Size  =      75946 *
*Baskets :      115 : Basket Size=      12800 bytes  Compression=   8.04     *
*............................................................................*
*Br   25 :Weight.Weight : Float_t Weight[Weight_]                            *
*Entries :    50000 : Total  Size=     613184 bytes  File Size  =      77779 *
*Baskets :      115 : Basket Size=      12800 bytes  Compression=   7.85     *
*............................................................................*
*Br   26 :Weight_size : Weight_size/I                                        *
*Entries :    50000 : Total  Size=     211944 bytes  File Size  =      13904 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  15.05     *
*............................................................................*
*Br   27 :Particle  : Int_t Particle_                                        *
*Entries :    50000 : Total  Size=     555108 bytes  File Size  =     200874 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.04     *
*............................................................................*
*Br   28 :Particle.fUniqueID : UInt_t fUniqueID[Particle_]                   *
*Entries :    50000 : Total  Size=  278381871 bytes  File Size  =   86955126 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   3.20     *
*............................................................................*
*Br   29 :Particle.fBits : UInt_t fBits[Particle_]                           *
*Entries :    50000 : Total  Size=  417471085 bytes  File Size  =    2069364 *
*Baskets :      400 : Basket Size=    1694720 bytes  Compression= 201.73     *
*............................................................................*
*Br   30 :Particle.PID : Int_t PID[Particle_]                                *
*Entries :    50000 : Total  Size=  278380251 bytes  File Size  =   47682930 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   5.84     *
*............................................................................*
*Br   31 :Particle.Status : Int_t Status[Particle_]                          *
*Entries :    50000 : Total  Size=  278381061 bytes  File Size  =   23800597 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=  11.70     *
*............................................................................*
*Br   32 :Particle.IsPU : Int_t IsPU[Particle_]                              *
*Entries :    50000 : Total  Size=  278380521 bytes  File Size  =    1436787 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression= 193.75     *
*............................................................................*
*Br   33 :Particle.M1 : Int_t M1[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =   90178521 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   3.09     *
*............................................................................*
*Br   34 :Particle.M2 : Int_t M2[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =   12646317 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=  22.01     *
*............................................................................*
*Br   35 :Particle.D1 : Int_t D1[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  100452633 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   2.77     *
*............................................................................*
*Br   36 :Particle.D2 : Int_t D2[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  100435054 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   2.77     *
*............................................................................*
*Br   37 :Particle.Charge : Int_t Charge[Particle_]                          *
*Entries :    50000 : Total  Size=  278381061 bytes  File Size  =   18251264 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=  15.25     *
*............................................................................*
*Br   38 :Particle.Mass : Float_t Mass[Particle_]                            *
*Entries :    50000 : Total  Size=  278380521 bytes  File Size  =   57990992 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   4.80     *
*............................................................................*
*Br   39 :Particle.E : Float_t E[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =  249810347 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.11     *
*............................................................................*
*Br   40 :Particle.Px : Float_t Px[Particle_]                                *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  251911963 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.11     *
*............................................................................*
*Br   41 :Particle.Py : Float_t Py[Particle_]                                *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  251676182 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.11     *
*............................................................................*
*Br   42 :Particle.Pz : Float_t Pz[Particle_]                                *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  255739692 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.09     *
*............................................................................*
*Br   43 :Particle.P : Float_t P[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =    1435729 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression= 193.89     *
*............................................................................*
*Br   44 :Particle.PT : Float_t PT[Particle_]                                *
*Entries :    50000 : Total  Size=  278379981 bytes  File Size  =  246624359 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.13     *
*............................................................................*
*Br   45 :Particle.Eta : Float_t Eta[Particle_]                              *
*Entries :    50000 : Total  Size=  278380251 bytes  File Size  =  236720817 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.18     *
*............................................................................*
*Br   46 :Particle.Phi : Float_t Phi[Particle_]                              *
*Entries :    50000 : Total  Size=  278380251 bytes  File Size  =  241655058 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.15     *
*............................................................................*
*Br   47 :Particle.Rapidity : Float_t Rapidity[Particle_]                    *
*Entries :    50000 : Total  Size=  278381601 bytes  File Size  =  236872942 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   1.18     *
*............................................................................*
*Br   48 :Particle.T : Float_t T[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =   37868063 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   7.35     *
*............................................................................*
*Br   49 :Particle.X : Float_t X[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =   38263282 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   7.28     *
*............................................................................*
*Br   50 :Particle.Y : Float_t Y[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =   38264297 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   7.28     *
*............................................................................*
*Br   51 :Particle.Z : Float_t Z[Particle_]                                  *
*Entries :    50000 : Total  Size=  278379711 bytes  File Size  =   38356611 *
*Baskets :      266 : Basket Size=    1694720 bytes  Compression=   7.26     *
*............................................................................*
*Br   52 :Particle_size : Particle_size/I                                    *
*Entries :    50000 : Total  Size=     212182 bytes  File Size  =     128426 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   1.63     *
*............................................................................*
*Br   53 :Track     : Int_t Track_                                           *
*Entries :    50000 : Total  Size=     547632 bytes  File Size  =     165051 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.48     *
*............................................................................*
*Br   54 :Track.fUniqueID : UInt_t fUniqueID[Track_]                         *
*Entries :    50000 : Total  Size=   12309217 bytes  File Size  =    2852071 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   4.31     *
*............................................................................*
*Br   55 :Track.fBits : UInt_t fBits[Track_]                                 *
*Entries :    50000 : Total  Size=   18356688 bytes  File Size  =     255126 *
*Baskets :      117 : Basket Size=     279040 bytes  Compression=  71.94     *
*............................................................................*
*Br   56 :Track.PID : Int_t PID[Track_]                                      *
*Entries :    50000 : Total  Size=   12308497 bytes  File Size  =    2332131 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.28     *
*............................................................................*
*Br   57 :Track.Charge : Int_t Charge[Track_]                                *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =    1494411 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   8.23     *
*............................................................................*
*Br   58 :Track.P   : Float_t P[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =   11229184 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.10     *
*............................................................................*
*Br   59 :Track.PT  : Float_t PT[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =   11195380 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.10     *
*............................................................................*
*Br   60 :Track.Eta : Float_t Eta[Track_]                                    *
*Entries :    50000 : Total  Size=   12308497 bytes  File Size  =   11359165 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br   61 :Track.Phi : Float_t Phi[Track_]                                    *
*Entries :    50000 : Total  Size=   12308497 bytes  File Size  =   11368684 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br   62 :Track.CtgTheta : Float_t CtgTheta[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =   11434303 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br   63 :Track.C   : Float_t C[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =     214958 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.25     *
*............................................................................*
*Br   64 :Track.Mass : Float_t Mass[Track_]                                  *
*Entries :    50000 : Total  Size=   12308617 bytes  File Size  =    1494512 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   8.23     *
*............................................................................*
*Br   65 :Track.EtaOuter : Float_t EtaOuter[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =   11350234 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br   66 :Track.PhiOuter : Float_t PhiOuter[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =   11369810 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br   67 :Track.T   : Float_t T[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =    1748321 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   7.04     *
*............................................................................*
*Br   68 :Track.X   : Float_t X[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =    1774619 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   6.93     *
*............................................................................*
*Br   69 :Track.Y   : Float_t Y[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =    1774677 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   6.93     *
*............................................................................*
*Br   70 :Track.Z   : Float_t Z[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =    1775293 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   6.93     *
*............................................................................*
*Br   71 :Track.TOuter : Float_t TOuter[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =   10848049 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.13     *
*............................................................................*
*Br   72 :Track.XOuter : Float_t XOuter[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =   11271845 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.09     *
*............................................................................*
*Br   73 :Track.YOuter : Float_t YOuter[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =   11271245 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.09     *
*............................................................................*
*Br   74 :Track.ZOuter : Float_t ZOuter[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =    8865700 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.39     *
*............................................................................*
*Br   75 :Track.Xd  : Float_t Xd[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =    2632729 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   4.67     *
*............................................................................*
*Br   76 :Track.Yd  : Float_t Yd[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =    2354375 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.23     *
*............................................................................*
*Br   77 :Track.Zd  : Float_t Zd[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =    3301385 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.73     *
*............................................................................*
*Br   78 :Track.L   : Float_t L[Track_]                                      *
*Entries :    50000 : Total  Size=   12308257 bytes  File Size  =   10851322 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.13     *
*............................................................................*
*Br   79 :Track.D0  : Float_t D0[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =    3862627 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.19     *
*............................................................................*
*Br   80 :Track.DZ  : Float_t DZ[Track_]                                     *
*Entries :    50000 : Total  Size=   12308377 bytes  File Size  =    3301385 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.73     *
*............................................................................*
*Br   81 :Track.Nclusters : Float_t Nclusters[Track_]                        *
*Entries :    50000 : Total  Size=   12309217 bytes  File Size  =     215859 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.01     *
*............................................................................*
*Br   82 :Track.dNdx : Float_t dNdx[Track_]                                  *
*Entries :    50000 : Total  Size=   12308617 bytes  File Size  =     215108 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.21     *
*............................................................................*
*Br   83 :Track.ErrorP : Float_t ErrorP[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =     215783 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.03     *
*............................................................................*
*Br   84 :Track.ErrorPT : Float_t ErrorPT[Track_]                            *
*Entries :    50000 : Total  Size=   12308977 bytes  File Size  =     215884 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.00     *
*............................................................................*
*Br   85 :Track.ErrorPhi : Float_t ErrorPhi[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =     215498 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.11     *
*............................................................................*
*Br   86 :Track.ErrorCtgTheta : Float_t ErrorCtgTheta[Track_]                *
*Entries :    50000 : Total  Size=   12309697 bytes  File Size  =     216342 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.89     *
*............................................................................*
*Br   87 :Track.ErrorT : Float_t ErrorT[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =     215783 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.03     *
*............................................................................*
*Br   88 :Track.ErrorD0 : Float_t ErrorD0[Track_]                            *
*Entries :    50000 : Total  Size=   12308977 bytes  File Size  =     215884 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.00     *
*............................................................................*
*Br   89 :Track.ErrorDZ : Float_t ErrorDZ[Track_]                            *
*Entries :    50000 : Total  Size=   12308977 bytes  File Size  =     215884 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.00     *
*............................................................................*
*Br   90 :Track.ErrorC : Float_t ErrorC[Track_]                              *
*Entries :    50000 : Total  Size=   12308857 bytes  File Size  =     215783 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.03     *
*............................................................................*
*Br   91 :Track.ErrorD0Phi : Float_t ErrorD0Phi[Track_]                      *
*Entries :    50000 : Total  Size=   12309337 bytes  File Size  =     216218 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.92     *
*............................................................................*
*Br   92 :Track.ErrorD0C : Float_t ErrorD0C[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =     215498 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.11     *
*............................................................................*
*Br   93 :Track.ErrorD0DZ : Float_t ErrorD0DZ[Track_]                        *
*Entries :    50000 : Total  Size=   12309217 bytes  File Size  =     215859 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.01     *
*............................................................................*
*Br   94 :Track.ErrorD0CtgTheta : Float_t ErrorD0CtgTheta[Track_]            *
*Entries :    50000 : Total  Size=   12309937 bytes  File Size  =     216845 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.76     *
*............................................................................*
*Br   95 :Track.ErrorPhiC : Float_t ErrorPhiC[Track_]                        *
*Entries :    50000 : Total  Size=   12309217 bytes  File Size  =     215859 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.01     *
*............................................................................*
*Br   96 :Track.ErrorPhiDZ : Float_t ErrorPhiDZ[Track_]                      *
*Entries :    50000 : Total  Size=   12309337 bytes  File Size  =     216218 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.92     *
*............................................................................*
*Br   97 :Track.ErrorPhiCtgTheta : Float_t ErrorPhiCtgTheta[Track_]          *
*Entries :    50000 : Total  Size=   12310057 bytes  File Size  =     216448 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.86     *
*............................................................................*
*Br   98 :Track.ErrorCDZ : Float_t ErrorCDZ[Track_]                          *
*Entries :    50000 : Total  Size=   12309097 bytes  File Size  =     215498 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.11     *
*............................................................................*
*Br   99 :Track.ErrorCCtgTheta : Float_t ErrorCCtgTheta[Track_]              *
*Entries :    50000 : Total  Size=   12309817 bytes  File Size  =     216729 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.79     *
*............................................................................*
*Br  100 :Track.ErrorDZCtgTheta : Float_t ErrorDZCtgTheta[Track_]            *
*Entries :    50000 : Total  Size=   12309937 bytes  File Size  =     216845 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.76     *
*............................................................................*
*Br  101 :Track.Particle : TRef Particle[Track_]                             *
*Entries :    50000 : Total  Size=   36500801 bytes  File Size  =    7251613 *
*Baskets :      119 : Basket Size=     551424 bytes  Compression=   5.03     *
*............................................................................*
*Br  102 :Track.VertexIndex : Int_t VertexIndex[Track_]                      *
*Entries :    50000 : Total  Size=   12309457 bytes  File Size  =     216463 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.85     *
*............................................................................*
*Br  103 :Track_size : Track_size/I                                          *
*Entries :    50000 : Total  Size=     211825 bytes  File Size  =      90284 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   2.32     *
*............................................................................*
*Br  104 :Tower     : Int_t Tower_                                           *
*Entries :    50000 : Total  Size=     451562 bytes  File Size  =     177562 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.31     *
*............................................................................*
*Br  105 :Tower.fUniqueID : UInt_t fUniqueID[Tower_]                         *
*Entries :    50000 : Total  Size=   39519569 bytes  File Size  =   16576836 *
*Baskets :      120 : Basket Size=     591360 bytes  Compression=   2.38     *
*............................................................................*
*Br  106 :Tower.fBits : UInt_t fBits[Tower_]                                 *
*Entries :    50000 : Total  Size=   59172192 bytes  File Size  =     451169 *
*Baskets :      123 : Basket Size=     883200 bytes  Compression= 131.15     *
*............................................................................*
*Br  107 :Tower.ET  : Float_t ET[Tower_]                                     *
*Entries :    50000 : Total  Size=   39518701 bytes  File Size  =   36233897 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.09     *
*............................................................................*
*Br  108 :Tower.Eta : Float_t Eta[Tower_]                                    *
*Entries :    50000 : Total  Size=   39518825 bytes  File Size  =   36072832 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.10     *
*............................................................................*
*Br  109 :Tower.Phi : Float_t Phi[Tower_]                                    *
*Entries :    50000 : Total  Size=   39518825 bytes  File Size  =   36604186 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.08     *
*............................................................................*
*Br  110 :Tower.E   : Float_t E[Tower_]                                      *
*Entries :    50000 : Total  Size=   39518577 bytes  File Size  =   36111829 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.09     *
*............................................................................*
*Br  111 :Tower.T   : Float_t T[Tower_]                                      *
*Entries :    50000 : Total  Size=   39518577 bytes  File Size  =   31536211 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.25     *
*............................................................................*
*Br  112 :Tower.NTimeHits : Int_t NTimeHits[Tower_]                          *
*Entries :    50000 : Total  Size=   39519569 bytes  File Size  =     355820 *
*Baskets :      120 : Basket Size=     591360 bytes  Compression= 111.06     *
*............................................................................*
*Br  113 :Tower.Eem : Float_t Eem[Tower_]                                    *
*Entries :    50000 : Total  Size=   39518825 bytes  File Size  =   20741670 *
*Baskets :      120 : Basket Size=     590848 bytes  Compression=   1.91     *
*............................................................................*
*Br  114 :Tower.Ehad : Float_t Ehad[Tower_]                                  *
*Entries :    50000 : Total  Size=   39518949 bytes  File Size  =   15798232 *
*Baskets :      120 : Basket Size=     591360 bytes  Compression=   2.50     *
*............................................................................*
*Br  115 :Tower.Edges[4] : Float_t Edges[Tower_]                             *
*Entries :    50000 : Total  Size=  157437980 bytes  File Size  =   52712042 *
*Baskets :      136 : Basket Size=    1694720 bytes  Compression=   2.99     *
*............................................................................*
*Br  116 :Tower.Particles : TRefArray Particles[Tower_]                      *
*Entries :    50000 : Total  Size=  310614100 bytes  File Size  =   35103484 *
*Baskets :      271 : Basket Size=    1694720 bytes  Compression=   8.85     *
*............................................................................*
*Br  117 :Tower_size : Tower_size/I                                          *
*Entries :    50000 : Total  Size=     211825 bytes  File Size  =     108052 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   1.94     *
*............................................................................*
*Br  118 :EFlowTrack : Int_t EFlowTrack_                                     *
*Entries :    50000 : Total  Size=     549217 bytes  File Size  =     165955 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.47     *
*............................................................................*
*Br  119 :EFlowTrack.fUniqueID : UInt_t fUniqueID[EFlowTrack_]               *
*Entries :    50000 : Total  Size=   12309827 bytes  File Size  =    5618620 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   2.19     *
*............................................................................*
*Br  120 :EFlowTrack.fBits : UInt_t fBits[EFlowTrack_]                       *
*Entries :    50000 : Total  Size=   18357303 bytes  File Size  =     255481 *
*Baskets :      117 : Basket Size=     279040 bytes  Compression=  71.84     *
*............................................................................*
*Br  121 :EFlowTrack.PID : Int_t PID[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12309107 bytes  File Size  =    2509789 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   4.90     *
*............................................................................*
*Br  122 :EFlowTrack.Charge : Int_t Charge[EFlowTrack_]                      *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =    1502732 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   8.19     *
*............................................................................*
*Br  123 :EFlowTrack.P : Float_t P[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =   11228617 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.10     *
*............................................................................*
*Br  124 :EFlowTrack.PT : Float_t PT[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =   11201226 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.10     *
*............................................................................*
*Br  125 :EFlowTrack.Eta : Float_t Eta[EFlowTrack_]                          *
*Entries :    50000 : Total  Size=   12309107 bytes  File Size  =   11359525 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br  126 :EFlowTrack.Phi : Float_t Phi[EFlowTrack_]                          *
*Entries :    50000 : Total  Size=   12309107 bytes  File Size  =   11369249 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br  127 :EFlowTrack.CtgTheta : Float_t CtgTheta[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =   11435413 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br  128 :EFlowTrack.C : Float_t C[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =     215783 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.03     *
*............................................................................*
*Br  129 :EFlowTrack.Mass : Float_t Mass[EFlowTrack_]                        *
*Entries :    50000 : Total  Size=   12309227 bytes  File Size  =    1604488 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   7.67     *
*............................................................................*
*Br  130 :EFlowTrack.EtaOuter : Float_t EtaOuter[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =   11350616 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br  131 :EFlowTrack.PhiOuter : Float_t PhiOuter[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =   11371109 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.08     *
*............................................................................*
*Br  132 :EFlowTrack.T : Float_t T[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =    2222096 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.54     *
*............................................................................*
*Br  133 :EFlowTrack.X : Float_t X[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =    2248479 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.47     *
*............................................................................*
*Br  134 :EFlowTrack.Y : Float_t Y[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =    2249160 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.47     *
*............................................................................*
*Br  135 :EFlowTrack.Z : Float_t Z[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =    2249251 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   5.47     *
*............................................................................*
*Br  136 :EFlowTrack.TOuter : Float_t TOuter[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =   10847697 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.13     *
*............................................................................*
*Br  137 :EFlowTrack.XOuter : Float_t XOuter[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =   11272476 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.09     *
*............................................................................*
*Br  138 :EFlowTrack.YOuter : Float_t YOuter[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =   11271876 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.09     *
*............................................................................*
*Br  139 :EFlowTrack.ZOuter : Float_t ZOuter[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =    8000546 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.54     *
*............................................................................*
*Br  140 :EFlowTrack.Xd : Float_t Xd[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =    3025999 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   4.07     *
*............................................................................*
*Br  141 :EFlowTrack.Yd : Float_t Yd[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =    2774369 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   4.44     *
*............................................................................*
*Br  142 :EFlowTrack.Zd : Float_t Zd[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =    3629961 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.39     *
*............................................................................*
*Br  143 :EFlowTrack.L : Float_t L[EFlowTrack_]                              *
*Entries :    50000 : Total  Size=   12308867 bytes  File Size  =   10850748 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   1.13     *
*............................................................................*
*Br  144 :EFlowTrack.D0 : Float_t D0[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =    4084630 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.01     *
*............................................................................*
*Br  145 :EFlowTrack.DZ : Float_t DZ[EFlowTrack_]                            *
*Entries :    50000 : Total  Size=   12308987 bytes  File Size  =    3629961 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=   3.39     *
*............................................................................*
*Br  146 :EFlowTrack.Nclusters : Float_t Nclusters[EFlowTrack_]              *
*Entries :    50000 : Total  Size=   12309827 bytes  File Size  =     216729 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.79     *
*............................................................................*
*Br  147 :EFlowTrack.dNdx : Float_t dNdx[EFlowTrack_]                        *
*Entries :    50000 : Total  Size=   12309227 bytes  File Size  =     215859 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  57.01     *
*............................................................................*
*Br  148 :EFlowTrack.ErrorP : Float_t ErrorP[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =     216348 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.88     *
*............................................................................*
*Br  149 :EFlowTrack.ErrorPT : Float_t ErrorPT[EFlowTrack_]                  *
*Entries :    50000 : Total  Size=   12309587 bytes  File Size  =     216021 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.97     *
*............................................................................*
*Br  150 :EFlowTrack.ErrorPhi : Float_t ErrorPhi[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =     216342 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.89     *
*............................................................................*
*Br  151 :EFlowTrack.ErrorCtgTheta : Float_t ErrorCtgTheta[EFlowTrack_]      *
*Entries :    50000 : Total  Size=   12310307 bytes  File Size  =     217169 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.67     *
*............................................................................*
*Br  152 :EFlowTrack.ErrorT : Float_t ErrorT[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =     216348 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.88     *
*............................................................................*
*Br  153 :EFlowTrack.ErrorD0 : Float_t ErrorD0[EFlowTrack_]                  *
*Entries :    50000 : Total  Size=   12309587 bytes  File Size  =     216021 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.97     *
*............................................................................*
*Br  154 :EFlowTrack.ErrorDZ : Float_t ErrorDZ[EFlowTrack_]                  *
*Entries :    50000 : Total  Size=   12309587 bytes  File Size  =     216021 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.97     *
*............................................................................*
*Br  155 :EFlowTrack.ErrorC : Float_t ErrorC[EFlowTrack_]                    *
*Entries :    50000 : Total  Size=   12309467 bytes  File Size  =     216348 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.88     *
*............................................................................*
*Br  156 :EFlowTrack.ErrorD0Phi : Float_t ErrorD0Phi[EFlowTrack_]            *
*Entries :    50000 : Total  Size=   12309947 bytes  File Size  =     216845 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.76     *
*............................................................................*
*Br  157 :EFlowTrack.ErrorD0C : Float_t ErrorD0C[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =     216342 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.89     *
*............................................................................*
*Br  158 :EFlowTrack.ErrorD0DZ : Float_t ErrorD0DZ[EFlowTrack_]              *
*Entries :    50000 : Total  Size=   12309827 bytes  File Size  =     216729 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.79     *
*............................................................................*
*Br  159 :EFlowTrack.ErrorD0CtgTheta : Float_t ErrorD0CtgTheta[EFlowTrack_]  *
*Entries :    50000 : Total  Size=   12310547 bytes  File Size  =     216889 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.75     *
*............................................................................*
*Br  160 :EFlowTrack.ErrorPhiC : Float_t ErrorPhiC[EFlowTrack_]              *
*Entries :    50000 : Total  Size=   12309827 bytes  File Size  =     216729 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.79     *
*............................................................................*
*Br  161 :EFlowTrack.ErrorPhiDZ : Float_t ErrorPhiDZ[EFlowTrack_]            *
*Entries :    50000 : Total  Size=   12309947 bytes  File Size  =     216845 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.76     *
*............................................................................*
*Br  162 :EFlowTrack.ErrorPhiCtgTheta : Float_t ErrorPhiCtgTheta[EFlowTrack_]*
*Entries :    50000 : Total  Size=   12310667 bytes  File Size  =     217256 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.65     *
*............................................................................*
*Br  163 :EFlowTrack.ErrorCDZ : Float_t ErrorCDZ[EFlowTrack_]                *
*Entries :    50000 : Total  Size=   12309707 bytes  File Size  =     216342 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.89     *
*............................................................................*
*Br  164 :EFlowTrack.ErrorCCtgTheta : Float_t ErrorCCtgTheta[EFlowTrack_]    *
*Entries :    50000 : Total  Size=   12310427 bytes  File Size  =     217280 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.64     *
*............................................................................*
*Br  165 :EFlowTrack.ErrorDZCtgTheta : Float_t ErrorDZCtgTheta[EFlowTrack_]  *
*Entries :    50000 : Total  Size=   12310547 bytes  File Size  =     216889 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.75     *
*............................................................................*
*Br  166 :EFlowTrack.Particle : TRef Particle[EFlowTrack_]                   *
*Entries :    50000 : Total  Size=   36501426 bytes  File Size  =    7954915 *
*Baskets :      119 : Basket Size=     551424 bytes  Compression=   4.59     *
*............................................................................*
*Br  167 :EFlowTrack.VertexIndex : Int_t VertexIndex[EFlowTrack_]            *
*Entries :    50000 : Total  Size=   12310067 bytes  File Size  =     216629 *
*Baskets :      116 : Basket Size=     188416 bytes  Compression=  56.81     *
*............................................................................*
*Br  168 :EFlowTrack_size : EFlowTrack_size/I                                *
*Entries :    50000 : Total  Size=     212420 bytes  File Size  =      90859 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   2.31     *
*............................................................................*
*Br  169 :EFlowPhoton : Int_t EFlowPhoton_                                   *
*Entries :    50000 : Total  Size=     449276 bytes  File Size  =     166825 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.46     *
*............................................................................*
*Br  170 :EFlowPhoton.fUniqueID : UInt_t fUniqueID[EFlowPhoton_]             *
*Entries :    50000 : Total  Size=   22326939 bytes  File Size  =    9128593 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   2.45     *
*............................................................................*
*Br  171 :EFlowPhoton.fBits : UInt_t fBits[EFlowPhoton_]                     *
*Entries :    50000 : Total  Size=   33382784 bytes  File Size  =     332857 *
*Baskets :      119 : Basket Size=     495104 bytes  Compression= 100.28     *
*............................................................................*
*Br  172 :EFlowPhoton.ET : Float_t ET[EFlowPhoton_]                          *
*Entries :    50000 : Total  Size=   22326085 bytes  File Size  =   20371411 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.10     *
*............................................................................*
*Br  173 :EFlowPhoton.Eta : Float_t Eta[EFlowPhoton_]                        *
*Entries :    50000 : Total  Size=   22326207 bytes  File Size  =   20424799 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.09     *
*............................................................................*
*Br  174 :EFlowPhoton.Phi : Float_t Phi[EFlowPhoton_]                        *
*Entries :    50000 : Total  Size=   22326207 bytes  File Size  =   20665044 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.08     *
*............................................................................*
*Br  175 :EFlowPhoton.E : Float_t E[EFlowPhoton_]                            *
*Entries :    50000 : Total  Size=   22325963 bytes  File Size  =   20252784 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.10     *
*............................................................................*
*Br  176 :EFlowPhoton.T : Float_t T[EFlowPhoton_]                            *
*Entries :    50000 : Total  Size=   22325963 bytes  File Size  =   17910960 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.25     *
*............................................................................*
*Br  177 :EFlowPhoton.NTimeHits : Int_t NTimeHits[EFlowPhoton_]              *
*Entries :    50000 : Total  Size=   22326939 bytes  File Size  =     274955 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=  81.19     *
*............................................................................*
*Br  178 :EFlowPhoton.Eem : Float_t Eem[EFlowPhoton_]                        *
*Entries :    50000 : Total  Size=   22326207 bytes  File Size  =   20252913 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=   1.10     *
*............................................................................*
*Br  179 :EFlowPhoton.Ehad : Float_t Ehad[EFlowPhoton_]                      *
*Entries :    50000 : Total  Size=   22326329 bytes  File Size  =     274076 *
*Baskets :      118 : Basket Size=     332288 bytes  Compression=  81.45     *
*............................................................................*
*Br  180 :EFlowPhoton.Edges[4] : Float_t Edges[EFlowPhoton_]                 *
*Entries :    50000 : Total  Size=   88665060 bytes  File Size  =   27999474 *
*Baskets :      126 : Basket Size=    1308160 bytes  Compression=   3.17     *
*............................................................................*
*Br  181 :EFlowPhoton.Particles : TRefArray Particles[EFlowPhoton_]          *
*Entries :    50000 : Total  Size=  172630990 bytes  File Size  =   16910778 *
*Baskets :      138 : Basket Size=    1694720 bytes  Compression=  10.21     *
*............................................................................*
*Br  182 :EFlowPhoton_size : EFlowPhoton_size/I                              *
*Entries :    50000 : Total  Size=     212539 bytes  File Size  =      95202 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   2.20     *
*............................................................................*
*Br  183 :EFlowNeutralHadron : Int_t EFlowNeutralHadron_                     *
*Entries :    50000 : Total  Size=     449799 bytes  File Size  =     166598 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   2.47     *
*............................................................................*
*Br  184 :EFlowNeutralHadron.fUniqueID :                                     *
*         | UInt_t fUniqueID[EFlowNeutralHadron_]                            *
*Entries :    50000 : Total  Size=   13539883 bytes  File Size  =    5919345 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   2.29     *
*............................................................................*
*Br  185 :EFlowNeutralHadron.fBits : UInt_t fBits[EFlowNeutralHadron_]       *
*Entries :    50000 : Total  Size=   20201907 bytes  File Size  =     266504 *
*Baskets :      117 : Basket Size=     307200 bytes  Compression=  75.79     *
*............................................................................*
*Br  186 :EFlowNeutralHadron.ET : Float_t ET[EFlowNeutralHadron_]            *
*Entries :    50000 : Total  Size=   13539043 bytes  File Size  =   12259879 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.10     *
*............................................................................*
*Br  187 :EFlowNeutralHadron.Eta : Float_t Eta[EFlowNeutralHadron_]          *
*Entries :    50000 : Total  Size=   13539163 bytes  File Size  =   12099072 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.12     *
*............................................................................*
*Br  188 :EFlowNeutralHadron.Phi : Float_t Phi[EFlowNeutralHadron_]          *
*Entries :    50000 : Total  Size=   13539163 bytes  File Size  =   12520982 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.08     *
*............................................................................*
*Br  189 :EFlowNeutralHadron.E : Float_t E[EFlowNeutralHadron_]              *
*Entries :    50000 : Total  Size=   13538923 bytes  File Size  =   12248299 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.11     *
*............................................................................*
*Br  190 :EFlowNeutralHadron.T : Float_t T[EFlowNeutralHadron_]              *
*Entries :    50000 : Total  Size=   13538923 bytes  File Size  =   10213118 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.33     *
*............................................................................*
*Br  191 :EFlowNeutralHadron.NTimeHits : Int_t NTimeHits[EFlowNeutralHadron_]*
*Entries :    50000 : Total  Size=   13539883 bytes  File Size  =     225507 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=  60.03     *
*............................................................................*
*Br  192 :EFlowNeutralHadron.Eem : Float_t Eem[EFlowNeutralHadron_]          *
*Entries :    50000 : Total  Size=   13539163 bytes  File Size  =     224303 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=  60.35     *
*............................................................................*
*Br  193 :EFlowNeutralHadron.Ehad : Float_t Ehad[EFlowNeutralHadron_]        *
*Entries :    50000 : Total  Size=   13539283 bytes  File Size  =   12248463 *
*Baskets :      116 : Basket Size=     206848 bytes  Compression=   1.11     *
*............................................................................*
*Br  194 :EFlowNeutralHadron.Edges[4] : Float_t Edges[EFlowNeutralHadron_]   *
*Entries :    50000 : Total  Size=   53514796 bytes  File Size  =   14734838 *
*Baskets :      122 : Basket Size=     807424 bytes  Compression=   3.63     *
*............................................................................*
*Br  195 :EFlowNeutralHadron.Particles :                                     *
*         | TRefArray Particles[EFlowNeutralHadron_]                         *
*Entries :    50000 : Total  Size=  106848782 bytes  File Size  =   13741474 *
*Baskets :      129 : Basket Size=    1609728 bytes  Compression=   7.78     *
*............................................................................*
*Br  196 :EFlowNeutralHadron_size : EFlowNeutralHadron_size/I                *
*Entries :    50000 : Total  Size=     213372 bytes  File Size  =      91672 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   2.30     *
*............................................................................*
*Br  197 :GenJet    : Int_t GenJet_                                          *
*Entries :    50000 : Total  Size=     537374 bytes  File Size  =     128197 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   3.20     *
*............................................................................*
*Br  198 :GenJet.fUniqueID : UInt_t fUniqueID[GenJet_]                       *
*Entries :    50000 : Total  Size=    1185305 bytes  File Size  =     109608 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  199 :GenJet.fBits : UInt_t fBits[GenJet_]                               *
*Entries :    50000 : Total  Size=    1184829 bytes  File Size  =     109226 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.82     *
*............................................................................*
*Br  200 :GenJet.PT : Float_t PT[GenJet_]                                    *
*Entries :    50000 : Total  Size=    1184472 bytes  File Size  =     999938 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.18     *
*............................................................................*
*Br  201 :GenJet.Eta : Float_t Eta[GenJet_]                                  *
*Entries :    50000 : Total  Size=    1184591 bytes  File Size  =    1030623 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.15     *
*............................................................................*
*Br  202 :GenJet.Phi : Float_t Phi[GenJet_]                                  *
*Entries :    50000 : Total  Size=    1184591 bytes  File Size  =    1027061 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.15     *
*............................................................................*
*Br  203 :GenJet.T  : Float_t T[GenJet_]                                     *
*Entries :    50000 : Total  Size=    1184353 bytes  File Size  =    1038598 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.14     *
*............................................................................*
*Br  204 :GenJet.Mass : Float_t Mass[GenJet_]                                *
*Entries :    50000 : Total  Size=    1184710 bytes  File Size  =     992601 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.19     *
*............................................................................*
*Br  205 :GenJet.DeltaEta : Float_t DeltaEta[GenJet_]                        *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     571457 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   2.07     *
*............................................................................*
*Br  206 :GenJet.DeltaPhi : Float_t DeltaPhi[GenJet_]                        *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     561582 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   2.11     *
*............................................................................*
*Br  207 :GenJet.Flavor : UInt_t Flavor[GenJet_]                             *
*Entries :    50000 : Total  Size=    1184948 bytes  File Size  =     109452 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.80     *
*............................................................................*
*Br  208 :GenJet.FlavorAlgo : UInt_t FlavorAlgo[GenJet_]                     *
*Entries :    50000 : Total  Size=    1185424 bytes  File Size  =     109890 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.76     *
*............................................................................*
*Br  209 :GenJet.FlavorPhys : UInt_t FlavorPhys[GenJet_]                     *
*Entries :    50000 : Total  Size=    1185424 bytes  File Size  =     109890 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.76     *
*............................................................................*
*Br  210 :GenJet.BTag : UInt_t BTag[GenJet_]                                 *
*Entries :    50000 : Total  Size=    1184710 bytes  File Size  =     109213 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.82     *
*............................................................................*
*Br  211 :GenJet.BTagAlgo : UInt_t BTagAlgo[GenJet_]                         *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     109629 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  212 :GenJet.BTagPhys : UInt_t BTagPhys[GenJet_]                         *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     109629 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  213 :GenJet.TauTag : UInt_t TauTag[GenJet_]                             *
*Entries :    50000 : Total  Size=    1184948 bytes  File Size  =     109452 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.80     *
*............................................................................*
*Br  214 :GenJet.TauWeight : Float_t TauWeight[GenJet_]                      *
*Entries :    50000 : Total  Size=    1185305 bytes  File Size  =     109608 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  215 :GenJet.Charge : Int_t Charge[GenJet_]                              *
*Entries :    50000 : Total  Size=    1184948 bytes  File Size  =     379802 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   3.11     *
*............................................................................*
*Br  216 :GenJet.EhadOverEem : Float_t EhadOverEem[GenJet_]                  *
*Entries :    50000 : Total  Size=    1185543 bytes  File Size  =     110732 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.68     *
*............................................................................*
*Br  217 :GenJet.NCharged : Int_t NCharged[GenJet_]                          *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     387855 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   3.05     *
*............................................................................*
*Br  218 :GenJet.NNeutrals : Int_t NNeutrals[GenJet_]                        *
*Entries :    50000 : Total  Size=    1185305 bytes  File Size  =     395810 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   2.99     *
*............................................................................*
*Br  219 :GenJet.NeutralEnergyFraction :                                     *
*         | Float_t NeutralEnergyFraction[GenJet_]                           *
*Entries :    50000 : Total  Size=    1186733 bytes  File Size  =    1008523 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.17     *
*............................................................................*
*Br  220 :GenJet.ChargedEnergyFraction :                                     *
*         | Float_t ChargedEnergyFraction[GenJet_]                           *
*Entries :    50000 : Total  Size=    1186733 bytes  File Size  =     962504 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=   1.23     *
*............................................................................*
*Br  221 :GenJet.Beta : Float_t Beta[GenJet_]                                *
*Entries :    50000 : Total  Size=    1184710 bytes  File Size  =     109213 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.82     *
*............................................................................*
*Br  222 :GenJet.BetaStar : Float_t BetaStar[GenJet_]                        *
*Entries :    50000 : Total  Size=    1185186 bytes  File Size  =     109629 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  223 :GenJet.MeanSqDeltaR : Float_t MeanSqDeltaR[GenJet_]                *
*Entries :    50000 : Total  Size=    1185662 bytes  File Size  =     110062 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.75     *
*............................................................................*
*Br  224 :GenJet.PTD : Float_t PTD[GenJet_]                                  *
*Entries :    50000 : Total  Size=    1184591 bytes  File Size  =     108375 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.91     *
*............................................................................*
*Br  225 :GenJet.FracPt[5] : Float_t FracPt[GenJet_]                         *
*Entries :    50000 : Total  Size=    5072352 bytes  File Size  =     156526 *
*Baskets :      115 : Basket Size=      80896 bytes  Compression=  32.39     *
*............................................................................*
*Br  226 :GenJet.Tau[5] : Float_t Tau[GenJet_]                               *
*Entries :    50000 : Total  Size=    5071995 bytes  File Size  =     156355 *
*Baskets :      115 : Basket Size=      80896 bytes  Compression=  32.42     *
*............................................................................*
*Br  227 :GenJet.SoftDroppedJet : TLorentzVector SoftDroppedJet[GenJet_]     *
*Entries :    50000 : Total  Size=   15762598 bytes  File Size  =     288103 *
*Baskets :      117 : Basket Size=     243712 bytes  Compression=  54.70     *
*............................................................................*
*Br  228 :GenJet.SoftDroppedSubJet1 :                                        *
*         | TLorentzVector SoftDroppedSubJet1[GenJet_]                       *
*Entries :    50000 : Total  Size=   15763082 bytes  File Size  =     288605 *
*Baskets :      117 : Basket Size=     243712 bytes  Compression=  54.61     *
*............................................................................*
*Br  229 :GenJet.SoftDroppedSubJet2 :                                        *
*         | TLorentzVector SoftDroppedSubJet2[GenJet_]                       *
*Entries :    50000 : Total  Size=   15763082 bytes  File Size  =     288605 *
*Baskets :      117 : Basket Size=     243712 bytes  Compression=  54.61     *
*............................................................................*
*Br  230 :GenJet.TrimmedP4[5] : TLorentzVector TrimmedP4[GenJet_]            *
*Entries :    50000 : Total  Size=   77956179 bytes  File Size  =     828143 *
*Baskets :      125 : Basket Size=    1189888 bytes  Compression=  94.13     *
*............................................................................*
*Br  231 :GenJet.PrunedP4[5] : TLorentzVector PrunedP4[GenJet_]              *
*Entries :    50000 : Total  Size=   77956050 bytes  File Size  =     827999 *
*Baskets :      125 : Basket Size=    1189888 bytes  Compression=  94.15     *
*............................................................................*
*Br  232 :GenJet.SoftDroppedP4[5] : TLorentzVector SoftDroppedP4[GenJet_]    *
*Entries :    50000 : Total  Size=   77956695 bytes  File Size  =     828627 *
*Baskets :      125 : Basket Size=    1189888 bytes  Compression=  94.08     *
*............................................................................*
*Br  233 :GenJet.NSubJetsTrimmed : Int_t NSubJetsTrimmed[GenJet_]            *
*Entries :    50000 : Total  Size=    1186019 bytes  File Size  =     109684 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.79     *
*............................................................................*
*Br  234 :GenJet.NSubJetsPruned : Int_t NSubJetsPruned[GenJet_]              *
*Entries :    50000 : Total  Size=    1185900 bytes  File Size  =     110305 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.73     *
*............................................................................*
*Br  235 :GenJet.NSubJetsSoftDropped : Int_t NSubJetsSoftDropped[GenJet_]    *
*Entries :    50000 : Total  Size=    1186495 bytes  File Size  =     110113 *
*Baskets :      115 : Basket Size=      21504 bytes  Compression=  10.75     *
*............................................................................*
*Br  236 :GenJet.ExclYmerge23 : Double_t ExclYmerge23[GenJet_]               *
*Entries :    50000 : Total  Size=    2157426 bytes  File Size  =     117105 *
*Baskets :      115 : Basket Size=      36352 bytes  Compression=  18.40     *
*............................................................................*
*Br  237 :GenJet.ExclYmerge34 : Double_t ExclYmerge34[GenJet_]               *
*Entries :    50000 : Total  Size=    2157426 bytes  File Size  =     117105 *
*Baskets :      115 : Basket Size=      36352 bytes  Compression=  18.40     *
*............................................................................*
*Br  238 :GenJet.ExclYmerge45 : Double_t ExclYmerge45[GenJet_]               *
*Entries :    50000 : Total  Size=    2157426 bytes  File Size  =     117105 *
*Baskets :      115 : Basket Size=      36352 bytes  Compression=  18.40     *
*............................................................................*
*Br  239 :GenJet.ExclYmerge56 : Double_t ExclYmerge56[GenJet_]               *
*Entries :    50000 : Total  Size=    2157426 bytes  File Size  =     117105 *
*Baskets :      115 : Basket Size=      36352 bytes  Compression=  18.40     *
*............................................................................*
*Br  240 :GenJet.Constituents : TRefArray Constituents[GenJet_]              *
*Entries :    50000 : Total  Size=   24817568 bytes  File Size  =   11238019 *
*Baskets :      118 : Basket Size=     383488 bytes  Compression=   2.21     *
*............................................................................*
*Br  241 :GenJet.Particles : TRefArray Particles[GenJet_]                    *
*Entries :    50000 : Total  Size=   24817202 bytes  File Size  =   11237569 *
*Baskets :      118 : Basket Size=     383488 bytes  Compression=   2.21     *
*............................................................................*
*Br  242 :GenJet.Area : TLorentzVector Area[GenJet_]                         *
*Entries :    50000 : Total  Size=   15761388 bytes  File Size  =     286920 *
*Baskets :      117 : Basket Size=     243712 bytes  Compression=  54.92     *
*............................................................................*
*Br  243 :GenJet_size : GenJet_size/I                                        *
*Entries :    50000 : Total  Size=     211944 bytes  File Size  =      47784 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   4.38     *
*............................................................................*
*Br  244 :GenMissingET : Int_t GenMissingET_                                 *
*Entries :    50000 : Total  Size=     429203 bytes  File Size  =      81678 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.02     *
*............................................................................*
*Br  245 :GenMissingET.fUniqueID : UInt_t fUniqueID[GenMissingET_]           *
*Entries :    50000 : Total  Size=     414267 bytes  File Size  =      82259 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.00     *
*............................................................................*
*Br  246 :GenMissingET.fBits : UInt_t fBits[GenMissingET_]                   *
*Entries :    50000 : Total  Size=     413791 bytes  File Size  =      81799 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.03     *
*............................................................................*
*Br  247 :GenMissingET.MET : Float_t MET[GenMissingET_]                      *
*Entries :    50000 : Total  Size=     413553 bytes  File Size  =     288845 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.42     *
*............................................................................*
*Br  248 :GenMissingET.Eta : Float_t Eta[GenMissingET_]                      *
*Entries :    50000 : Total  Size=     413553 bytes  File Size  =     285734 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.44     *
*............................................................................*
*Br  249 :GenMissingET.Phi : Float_t Phi[GenMissingET_]                      *
*Entries :    50000 : Total  Size=     413553 bytes  File Size  =     287284 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.43     *
*............................................................................*
*Br  250 :GenMissingET_size : GenMissingET_size/I                            *
*Entries :    50000 : Total  Size=     212658 bytes  File Size  =      14594 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  14.39     *
*............................................................................*
*Br  251 :Jet       : Int_t Jet_                                             *
*Entries :    50000 : Total  Size=     536091 bytes  File Size  =     126283 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   3.24     *
*............................................................................*
*Br  252 :Jet.fUniqueID : UInt_t fUniqueID[Jet_]                             *
*Entries :    50000 : Total  Size=     817742 bytes  File Size  =     101357 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.04     *
*............................................................................*
*Br  253 :Jet.fBits : UInt_t fBits[Jet_]                                     *
*Entries :    50000 : Total  Size=     817266 bytes  File Size  =     100882 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.07     *
*............................................................................*
*Br  254 :Jet.PT    : Float_t PT[Jet_]                                       *
*Entries :    50000 : Total  Size=     816909 bytes  File Size  =     662282 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.23     *
*............................................................................*
*Br  255 :Jet.Eta   : Float_t Eta[Jet_]                                      *
*Entries :    50000 : Total  Size=     817028 bytes  File Size  =     680132 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.20     *
*............................................................................*
*Br  256 :Jet.Phi   : Float_t Phi[Jet_]                                      *
*Entries :    50000 : Total  Size=     817028 bytes  File Size  =     680976 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.20     *
*............................................................................*
*Br  257 :Jet.T     : Float_t T[Jet_]                                        *
*Entries :    50000 : Total  Size=     816790 bytes  File Size  =     634710 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.28     *
*............................................................................*
*Br  258 :Jet.Mass  : Float_t Mass[Jet_]                                     *
*Entries :    50000 : Total  Size=     817147 bytes  File Size  =     662013 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.23     *
*............................................................................*
*Br  259 :Jet.DeltaEta : Float_t DeltaEta[Jet_]                              *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     471539 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.73     *
*............................................................................*
*Br  260 :Jet.DeltaPhi : Float_t DeltaPhi[Jet_]                              *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     468760 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.74     *
*............................................................................*
*Br  261 :Jet.Flavor : UInt_t Flavor[Jet_]                                   *
*Entries :    50000 : Total  Size=     817385 bytes  File Size  =     235806 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   3.45     *
*............................................................................*
*Br  262 :Jet.FlavorAlgo : UInt_t FlavorAlgo[Jet_]                           *
*Entries :    50000 : Total  Size=     817861 bytes  File Size  =     100591 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.10     *
*............................................................................*
*Br  263 :Jet.FlavorPhys : UInt_t FlavorPhys[Jet_]                           *
*Entries :    50000 : Total  Size=     817861 bytes  File Size  =     100591 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.10     *
*............................................................................*
*Br  264 :Jet.BTag  : UInt_t BTag[Jet_]                                      *
*Entries :    50000 : Total  Size=     817147 bytes  File Size  =     111677 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.29     *
*............................................................................*
*Br  265 :Jet.BTagAlgo : UInt_t BTagAlgo[Jet_]                               *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     103429 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.88     *
*............................................................................*
*Br  266 :Jet.BTagPhys : UInt_t BTagPhys[Jet_]                               *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     103453 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.88     *
*............................................................................*
*Br  267 :Jet.TauTag : UInt_t TauTag[Jet_]                                   *
*Entries :    50000 : Total  Size=     817385 bytes  File Size  =     105942 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.69     *
*............................................................................*
*Br  268 :Jet.TauWeight : Float_t TauWeight[Jet_]                            *
*Entries :    50000 : Total  Size=     817742 bytes  File Size  =     101357 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.04     *
*............................................................................*
*Br  269 :Jet.Charge : Int_t Charge[Jet_]                                    *
*Entries :    50000 : Total  Size=     817385 bytes  File Size  =     239120 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   3.41     *
*............................................................................*
*Br  270 :Jet.EhadOverEem : Float_t EhadOverEem[Jet_]                        *
*Entries :    50000 : Total  Size=     817980 bytes  File Size  =     662356 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   1.23     *
*............................................................................*
*Br  271 :Jet.NCharged : Int_t NCharged[Jet_]                                *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     145055 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   5.62     *
*............................................................................*
*Br  272 :Jet.NNeutrals : Int_t NNeutrals[Jet_]                              *
*Entries :    50000 : Total  Size=     817742 bytes  File Size  =     293670 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   2.78     *
*............................................................................*
*Br  273 :Jet.NeutralEnergyFraction : Float_t NeutralEnergyFraction[Jet_]    *
*Entries :    50000 : Total  Size=     819170 bytes  File Size  =     195483 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   4.18     *
*............................................................................*
*Br  274 :Jet.ChargedEnergyFraction : Float_t ChargedEnergyFraction[Jet_]    *
*Entries :    50000 : Total  Size=     819170 bytes  File Size  =     194522 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   4.20     *
*............................................................................*
*Br  275 :Jet.Beta  : Float_t Beta[Jet_]                                     *
*Entries :    50000 : Total  Size=     817147 bytes  File Size  =     100644 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.09     *
*............................................................................*
*Br  276 :Jet.BetaStar : Float_t BetaStar[Jet_]                              *
*Entries :    50000 : Total  Size=     817623 bytes  File Size  =     101155 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.06     *
*............................................................................*
*Br  277 :Jet.MeanSqDeltaR : Float_t MeanSqDeltaR[Jet_]                      *
*Entries :    50000 : Total  Size=     818099 bytes  File Size  =     101579 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.03     *
*............................................................................*
*Br  278 :Jet.PTD   : Float_t PTD[Jet_]                                      *
*Entries :    50000 : Total  Size=     817028 bytes  File Size  =     100628 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.09     *
*............................................................................*
*Br  279 :Jet.FracPt[5] : Float_t FracPt[Jet_]                               *
*Entries :    50000 : Total  Size=    3235989 bytes  File Size  =     134593 *
*Baskets :      115 : Basket Size=      53248 bytes  Compression=  24.02     *
*............................................................................*
*Br  280 :Jet.Tau[5] : Float_t Tau[Jet_]                                     *
*Entries :    50000 : Total  Size=    3235632 bytes  File Size  =     133495 *
*Baskets :      115 : Basket Size=      53248 bytes  Compression=  24.22     *
*............................................................................*
*Br  281 :Jet.SoftDroppedJet : TLorentzVector SoftDroppedJet[Jet_]           *
*Entries :    50000 : Total  Size=    9886913 bytes  File Size  =     226995 *
*Baskets :      116 : Basket Size=     155648 bytes  Compression=  43.54     *
*............................................................................*
*Br  282 :Jet.SoftDroppedSubJet1 : TLorentzVector SoftDroppedSubJet1[Jet_]   *
*Entries :    50000 : Total  Size=    9887393 bytes  File Size  =     227448 *
*Baskets :      116 : Basket Size=     155648 bytes  Compression=  43.46     *
*............................................................................*
*Br  283 :Jet.SoftDroppedSubJet2 : TLorentzVector SoftDroppedSubJet2[Jet_]   *
*Entries :    50000 : Total  Size=    9887393 bytes  File Size  =     227448 *
*Baskets :      116 : Basket Size=     155648 bytes  Compression=  43.46     *
*............................................................................*
*Br  284 :Jet.TrimmedP4[5] : TLorentzVector TrimmedP4[Jet_]                  *
*Entries :    50000 : Total  Size=   48579330 bytes  File Size  =     582492 *
*Baskets :      121 : Basket Size=     751616 bytes  Compression=  83.39     *
*............................................................................*
*Br  285 :Jet.PrunedP4[5] : TLorentzVector PrunedP4[Jet_]                    *
*Entries :    50000 : Total  Size=   48579205 bytes  File Size  =     582361 *
*Baskets :      121 : Basket Size=     751616 bytes  Compression=  83.41     *
*............................................................................*
*Br  286 :Jet.SoftDroppedP4[5] : TLorentzVector SoftDroppedP4[Jet_]          *
*Entries :    50000 : Total  Size=   48579830 bytes  File Size  =     582996 *
*Baskets :      121 : Basket Size=     751616 bytes  Compression=  83.32     *
*............................................................................*
*Br  287 :Jet.NSubJetsTrimmed : Int_t NSubJetsTrimmed[Jet_]                  *
*Entries :    50000 : Total  Size=     818456 bytes  File Size  =     102043 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.99     *
*............................................................................*
*Br  288 :Jet.NSubJetsPruned : Int_t NSubJetsPruned[Jet_]                    *
*Entries :    50000 : Total  Size=     818337 bytes  File Size  =     101077 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   8.07     *
*............................................................................*
*Br  289 :Jet.NSubJetsSoftDropped : Int_t NSubJetsSoftDropped[Jet_]          *
*Entries :    50000 : Total  Size=     818932 bytes  File Size  =     102551 *
*Baskets :      115 : Basket Size=      16384 bytes  Compression=   7.96     *
*............................................................................*
*Br  290 :Jet.ExclYmerge23 : Double_t ExclYmerge23[Jet_]                     *
*Entries :    50000 : Total  Size=    1422663 bytes  File Size  =     103786 *
*Baskets :      115 : Basket Size=      25600 bytes  Compression=  13.68     *
*............................................................................*
*Br  291 :Jet.ExclYmerge34 : Double_t ExclYmerge34[Jet_]                     *
*Entries :    50000 : Total  Size=    1422663 bytes  File Size  =     103786 *
*Baskets :      115 : Basket Size=      25600 bytes  Compression=  13.68     *
*............................................................................*
*Br  292 :Jet.ExclYmerge45 : Double_t ExclYmerge45[Jet_]                     *
*Entries :    50000 : Total  Size=    1422663 bytes  File Size  =     103786 *
*Baskets :      115 : Basket Size=      25600 bytes  Compression=  13.68     *
*............................................................................*
*Br  293 :Jet.ExclYmerge56 : Double_t ExclYmerge56[Jet_]                     *
*Entries :    50000 : Total  Size=    1422663 bytes  File Size  =     103786 *
*Baskets :      115 : Basket Size=      25600 bytes  Compression=  13.68     *
*............................................................................*
*Br  294 :Jet.Constituents : TRefArray Constituents[Jet_]                    *
*Entries :    50000 : Total  Size=   12656444 bytes  File Size  =    5377502 *
*Baskets :      116 : Basket Size=     199680 bytes  Compression=   2.35     *
*............................................................................*
*Br  295 :Jet.Particles : TRefArray Particles[Jet_]                          *
*Entries :    50000 : Total  Size=   15440883 bytes  File Size  =    6880651 *
*Baskets :      117 : Basket Size=     243712 bytes  Compression=   2.24     *
*............................................................................*
*Br  296 :Jet.Area  : TLorentzVector Area[Jet_]                              *
*Entries :    50000 : Total  Size=    9885713 bytes  File Size  =     225828 *
*Baskets :      116 : Basket Size=     155648 bytes  Compression=  43.76     *
*............................................................................*
*Br  297 :Jet_size  : Jet_size/I                                             *
*Entries :    50000 : Total  Size=     211587 bytes  File Size  =      46607 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   4.48     *
*............................................................................*
*Br  298 :Electron  : Int_t Electron_                                        *
*Entries :    50000 : Total  Size=     466660 bytes  File Size  =      81078 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.06     *
*............................................................................*
*Br  299 :Electron.fUniqueID : UInt_t fUniqueID[Electron_]                   *
*Entries :    50000 : Total  Size=     213899 bytes  File Size  =      15279 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.82     *
*............................................................................*
*Br  300 :Electron.fBits : UInt_t fBits[Electron_]                           *
*Entries :    50000 : Total  Size=     213423 bytes  File Size  =      14819 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.22     *
*............................................................................*
*Br  301 :Electron.PT : Float_t PT[Electron_]                                *
*Entries :    50000 : Total  Size=     213066 bytes  File Size  =      14590 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.42     *
*............................................................................*
*Br  302 :Electron.Eta : Float_t Eta[Electron_]                              *
*Entries :    50000 : Total  Size=     213185 bytes  File Size  =      14706 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.31     *
*............................................................................*
*Br  303 :Electron.Phi : Float_t Phi[Electron_]                              *
*Entries :    50000 : Total  Size=     213185 bytes  File Size  =      14705 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.31     *
*............................................................................*
*Br  304 :Electron.T : Float_t T[Electron_]                                  *
*Entries :    50000 : Total  Size=     212947 bytes  File Size  =      14475 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.52     *
*............................................................................*
*Br  305 :Electron.Charge : Int_t Charge[Electron_]                          *
*Entries :    50000 : Total  Size=     213542 bytes  File Size  =      15012 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.04     *
*............................................................................*
*Br  306 :Electron.EhadOverEem : Float_t EhadOverEem[Electron_]              *
*Entries :    50000 : Total  Size=     214137 bytes  File Size  =      15509 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.63     *
*............................................................................*
*Br  307 :Electron.Particle : TRef Particle[Electron_]                       *
*Entries :    50000 : Total  Size=     214012 bytes  File Size  =      15316 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.79     *
*............................................................................*
*Br  308 :Electron.IsolationVar : Float_t IsolationVar[Electron_]            *
*Entries :    50000 : Total  Size=     214256 bytes  File Size  =      15680 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.49     *
*............................................................................*
*Br  309 :Electron.IsolationVarRhoCorr :                                     *
*         | Float_t IsolationVarRhoCorr[Electron_]                           *
*Entries :    50000 : Total  Size=     215089 bytes  File Size  =      16485 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  12.88     *
*............................................................................*
*Br  310 :Electron.SumPtCharged : Float_t SumPtCharged[Electron_]            *
*Entries :    50000 : Total  Size=     214256 bytes  File Size  =      15658 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.51     *
*............................................................................*
*Br  311 :Electron.SumPtNeutral : Float_t SumPtNeutral[Electron_]            *
*Entries :    50000 : Total  Size=     214256 bytes  File Size  =      15655 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.51     *
*............................................................................*
*Br  312 :Electron.SumPtChargedPU : Float_t SumPtChargedPU[Electron_]        *
*Entries :    50000 : Total  Size=     214494 bytes  File Size  =      15854 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  13.36     *
*............................................................................*
*Br  313 :Electron.SumPt : Float_t SumPt[Electron_]                          *
*Entries :    50000 : Total  Size=     213423 bytes  File Size  =      14870 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.17     *
*............................................................................*
*Br  314 :Electron.D0 : Float_t D0[Electron_]                                *
*Entries :    50000 : Total  Size=     213066 bytes  File Size  =      14484 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.52     *
*............................................................................*
*Br  315 :Electron.DZ : Float_t DZ[Electron_]                                *
*Entries :    50000 : Total  Size=     213066 bytes  File Size  =      14483 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.52     *
*............................................................................*
*Br  316 :Electron.ErrorD0 : Float_t ErrorD0[Electron_]                      *
*Entries :    50000 : Total  Size=     213661 bytes  File Size  =      15049 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.02     *
*............................................................................*
*Br  317 :Electron.ErrorDZ : Float_t ErrorDZ[Electron_]                      *
*Entries :    50000 : Total  Size=     213661 bytes  File Size  =      15049 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=  14.02     *
*............................................................................*
*Br  318 :Electron_size : Electron_size/I                                    *
*Entries :    50000 : Total  Size=     212182 bytes  File Size  =      13776 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  15.21     *
*............................................................................*
*Br  319 :Photon    : Int_t Photon_                                          *
*Entries :    50000 : Total  Size=     458115 bytes  File Size  =      93233 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   4.39     *
*............................................................................*
*Br  320 :Photon.fUniqueID : UInt_t fUniqueID[Photon_]                       *
*Entries :    50000 : Total  Size=     225185 bytes  File Size  =      26493 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   8.40     *
*............................................................................*
*Br  321 :Photon.fBits : UInt_t fBits[Photon_]                               *
*Entries :    50000 : Total  Size=     224709 bytes  File Size  =      26011 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   8.53     *
*............................................................................*
*Br  322 :Photon.PT : Float_t PT[Photon_]                                    *
*Entries :    50000 : Total  Size=     224352 bytes  File Size  =      37622 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   5.89     *
*............................................................................*
*Br  323 :Photon.Eta : Float_t Eta[Photon_]                                  *
*Entries :    50000 : Total  Size=     224471 bytes  File Size  =      37972 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   5.84     *
*............................................................................*
*Br  324 :Photon.Phi : Float_t Phi[Photon_]                                  *
*Entries :    50000 : Total  Size=     224471 bytes  File Size  =      37931 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   5.85     *
*............................................................................*
*Br  325 :Photon.E  : Float_t E[Photon_]                                     *
*Entries :    50000 : Total  Size=     224233 bytes  File Size  =      37520 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   5.90     *
*............................................................................*
*Br  326 :Photon.T  : Float_t T[Photon_]                                     *
*Entries :    50000 : Total  Size=     224233 bytes  File Size  =      37521 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   5.90     *
*............................................................................*
*Br  327 :Photon.EhadOverEem : Float_t EhadOverEem[Photon_]                  *
*Entries :    50000 : Total  Size=     225423 bytes  File Size  =      26741 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   8.33     *
*............................................................................*
*Br  328 :Photon.Particles : TRefArray Particles[Photon_]                    *
*Entries :    50000 : Total  Size=     303970 bytes  File Size  =      39097 *
*Baskets :      115 : Basket Size=       8192 bytes  Compression=   7.70     *
*............................................................................*
*Br  329 :Photon.IsolationVar : Float_t IsolationVar[Photon_]                *
*Entries :    50000 : Total  Size=     225542 bytes  File Size  =      33548 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   6.64     *
*............................................................................*
*Br  330 :Photon.IsolationVarRhoCorr : Float_t IsolationVarRhoCorr[Photon_]  *
*Entries :    50000 : Total  Size=     226375 bytes  File Size  =      34374 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   6.51     *
*............................................................................*
*Br  331 :Photon.SumPtCharged : Float_t SumPtCharged[Photon_]                *
*Entries :    50000 : Total  Size=     225542 bytes  File Size  =      31808 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   7.00     *
*............................................................................*
*Br  332 :Photon.SumPtNeutral : Float_t SumPtNeutral[Photon_]                *
*Entries :    50000 : Total  Size=     225542 bytes  File Size  =      30666 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   7.27     *
*............................................................................*
*Br  333 :Photon.SumPtChargedPU : Float_t SumPtChargedPU[Photon_]            *
*Entries :    50000 : Total  Size=     225780 bytes  File Size  =      27086 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   8.23     *
*............................................................................*
*Br  334 :Photon.SumPt : Float_t SumPt[Photon_]                              *
*Entries :    50000 : Total  Size=     224709 bytes  File Size  =      32697 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   6.79     *
*............................................................................*
*Br  335 :Photon.Status : Int_t Status[Photon_]                              *
*Entries :    50000 : Total  Size=     224828 bytes  File Size  =      26129 *
*Baskets :      115 : Basket Size=       7168 bytes  Compression=   8.50     *
*............................................................................*
*Br  336 :Photon_size : Photon_size/I                                        *
*Entries :    50000 : Total  Size=     211944 bytes  File Size  =      22242 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   9.41     *
*............................................................................*
*Br  337 :Muon      : Int_t Muon_                                            *
*Entries :    50000 : Total  Size=     463021 bytes  File Size  =     112194 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   3.65     *
*............................................................................*
*Br  338 :Muon.fUniqueID : UInt_t fUniqueID[Muon_]                           *
*Entries :    50000 : Total  Size=     548999 bytes  File Size  =     276237 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.98     *
*............................................................................*
*Br  339 :Muon.fBits : UInt_t fBits[Muon_]                                   *
*Entries :    50000 : Total  Size=     716373 bytes  File Size  =     100148 *
*Baskets :      115 : Basket Size=      14336 bytes  Compression=   7.13     *
*............................................................................*
*Br  340 :Muon.PT   : Float_t PT[Muon_]                                      *
*Entries :    50000 : Total  Size=     548166 bytes  File Size  =     410349 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.33     *
*............................................................................*
*Br  341 :Muon.Eta  : Float_t Eta[Muon_]                                     *
*Entries :    50000 : Total  Size=     548285 bytes  File Size  =     421720 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.29     *
*............................................................................*
*Br  342 :Muon.Phi  : Float_t Phi[Muon_]                                     *
*Entries :    50000 : Total  Size=     548285 bytes  File Size  =     422407 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.29     *
*............................................................................*
*Br  343 :Muon.T    : Float_t T[Muon_]                                       *
*Entries :    50000 : Total  Size=     548047 bytes  File Size  =     403959 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.35     *
*............................................................................*
*Br  344 :Muon.Charge : Int_t Charge[Muon_]                                  *
*Entries :    50000 : Total  Size=     548642 bytes  File Size  =     140893 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   3.87     *
*............................................................................*
*Br  345 :Muon.Particle : TRef Particle[Muon_]                               *
*Entries :    50000 : Total  Size=    1220280 bytes  File Size  =     334003 *
*Baskets :      115 : Basket Size=      22016 bytes  Compression=   3.65     *
*............................................................................*
*Br  346 :Muon.IsolationVar : Float_t IsolationVar[Muon_]                    *
*Entries :    50000 : Total  Size=     549356 bytes  File Size  =     323504 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.69     *
*............................................................................*
*Br  347 :Muon.IsolationVarRhoCorr : Float_t IsolationVarRhoCorr[Muon_]      *
*Entries :    50000 : Total  Size=     550189 bytes  File Size  =     324179 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.69     *
*............................................................................*
*Br  348 :Muon.SumPtCharged : Float_t SumPtCharged[Muon_]                    *
*Entries :    50000 : Total  Size=     549356 bytes  File Size  =     272904 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   2.00     *
*............................................................................*
*Br  349 :Muon.SumPtNeutral : Float_t SumPtNeutral[Muon_]                    *
*Entries :    50000 : Total  Size=     549356 bytes  File Size  =     255498 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   2.14     *
*............................................................................*
*Br  350 :Muon.SumPtChargedPU : Float_t SumPtChargedPU[Muon_]                *
*Entries :    50000 : Total  Size=     549594 bytes  File Size  =      92615 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   5.90     *
*............................................................................*
*Br  351 :Muon.SumPt : Float_t SumPt[Muon_]                                  *
*Entries :    50000 : Total  Size=     548523 bytes  File Size  =     318528 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   1.71     *
*............................................................................*
*Br  352 :Muon.D0   : Float_t D0[Muon_]                                      *
*Entries :    50000 : Total  Size=     548166 bytes  File Size  =     169571 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   3.22     *
*............................................................................*
*Br  353 :Muon.DZ   : Float_t DZ[Muon_]                                      *
*Entries :    50000 : Total  Size=     548166 bytes  File Size  =     143799 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   3.79     *
*............................................................................*
*Br  354 :Muon.ErrorD0 : Float_t ErrorD0[Muon_]                              *
*Entries :    50000 : Total  Size=     548761 bytes  File Size  =      91755 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   5.95     *
*............................................................................*
*Br  355 :Muon.ErrorDZ : Float_t ErrorDZ[Muon_]                              *
*Entries :    50000 : Total  Size=     548761 bytes  File Size  =      91755 *
*Baskets :      115 : Basket Size=      11776 bytes  Compression=   5.95     *
*............................................................................*
*Br  356 :Muon_size : Muon_size/I                                            *
*Entries :    50000 : Total  Size=     211706 bytes  File Size  =      35248 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=   5.93     *
*............................................................................*
*Br  357 :MissingET : Int_t MissingET_                                       *
*Entries :    50000 : Total  Size=     428780 bytes  File Size  =      80991 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.06     *
*............................................................................*
*Br  358 :MissingET.fUniqueID : UInt_t fUniqueID[MissingET_]                 *
*Entries :    50000 : Total  Size=     413904 bytes  File Size  =      82027 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.01     *
*............................................................................*
*Br  359 :MissingET.fBits : UInt_t fBits[MissingET_]                         *
*Entries :    50000 : Total  Size=     413428 bytes  File Size  =      81567 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.04     *
*............................................................................*
*Br  360 :MissingET.MET : Float_t MET[MissingET_]                            *
*Entries :    50000 : Total  Size=     413190 bytes  File Size  =     280464 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.46     *
*............................................................................*
*Br  361 :MissingET.Eta : Float_t Eta[MissingET_]                            *
*Entries :    50000 : Total  Size=     413190 bytes  File Size  =     279869 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.47     *
*............................................................................*
*Br  362 :MissingET.Phi : Float_t Phi[MissingET_]                            *
*Entries :    50000 : Total  Size=     413190 bytes  File Size  =     285929 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.44     *
*............................................................................*
*Br  363 :MissingET_size : MissingET_size/I                                  *
*Entries :    50000 : Total  Size=     212301 bytes  File Size  =      14249 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  14.71     *
*............................................................................*
*Br  364 :ScalarHT  : Int_t ScalarHT_                                        *
*Entries :    50000 : Total  Size=     423236 bytes  File Size  =      81218 *
*Baskets :      115 : Basket Size=      64000 bytes  Compression=   5.05     *
*............................................................................*
*Br  365 :ScalarHT.fUniqueID : UInt_t fUniqueID[ScalarHT_]                   *
*Entries :    50000 : Total  Size=     413783 bytes  File Size  =      81799 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.03     *
*............................................................................*
*Br  366 :ScalarHT.fBits : UInt_t fBits[ScalarHT_]                           *
*Entries :    50000 : Total  Size=     413307 bytes  File Size  =      81339 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   5.05     *
*............................................................................*
*Br  367 :ScalarHT.HT : Float_t HT[ScalarHT_]                                *
*Entries :    50000 : Total  Size=     412950 bytes  File Size  =     277883 *
*Baskets :      115 : Basket Size=       9728 bytes  Compression=   1.48     *
*............................................................................*
*Br  368 :ScalarHT_size : ScalarHT_size/I                                    *
*Entries :    50000 : Total  Size=     212182 bytes  File Size  =      14134 *
*Baskets :      115 : Basket Size=       3381 bytes  Compression=  14.83     *
*............................................................................*/