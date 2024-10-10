#pragma once
#include<TLorentzVector.h>


// Error message with colors!
const char* g_ERROR_MESSAGE = "\033[1;31mERROR:\033[0m";
const char* g_DEBUG_MESSAGE = "\033[1;94mDEBUG:\033[0m";
const char* g_INFO_MESSAGE = "\033[1;92mINFO:\033[0m";


// Enum for log levels.
enum LogLevel {
    ERROR = 1,
    INFO,
    DEBUG,
};

// Class for a logger object.
class Logger {
    public: 
        // Constructor
        Logger(LogLevel level) {
            m_logLevel = level;
        }

        // Destructor
        ~Logger() {}

        static const char* getLogMessage(LogLevel level){
            const char* logMessage = nullptr;
            if (level == LogLevel::INFO) {
                logMessage = g_INFO_MESSAGE;
            } else if (level == LogLevel::DEBUG) {
                logMessage = g_DEBUG_MESSAGE;
            } else if (level == LogLevel::ERROR) {
                logMessage = g_ERROR_MESSAGE;
            }
            return logMessage;
        }

        LogLevel getLogLevel() {
            return m_logLevel;
        }

        // Functions to print the log message.
        template <typename T>
        void operator () (LogLevel level, const T& message) {
            if (level <= getLogLevel()) std::cout << getLogMessage(level) << " " << message << std::endl;
        }

        template <typename T, typename V>
        void operator () (LogLevel level, const T& message1, const V& message2) {
            if (level <= getLogLevel()) std::cout << getLogMessage(level) << " " << message1 << message2 << std::endl;
        }

    private:
        LogLevel m_logLevel = LogLevel::ERROR;
};


Logger g_LOG = Logger(LogLevel::INFO);

namespace Kinematics {


/**
 * @brief Function to calculate delta phi between two angles.
 * @param phi_1 angle 1.
 * @param phi_2 angle 2.
 * 
 */
double del_phi(double phi_1, double phi_2){
    double pi=TMath::Pi();
    double phi_1_norm, phi_2_norm;
    if (phi_1<0.0){
        phi_1_norm=phi_1+2*pi;
    }else {
        phi_1_norm=phi_1;
    }

    if (phi_2<0.0){
        phi_2_norm=phi_2+2*pi;
    }else {
        phi_2_norm=phi_2;
    }
    double delta=std::abs(phi_1_norm-phi_2_norm);
    if (delta>pi){
        delta=2*pi-delta;
        delta=std::abs(delta);
    }

    return delta;
}

/**
 * @brief This function calculates the invariant mass of a list of particles
 * @param particles A list of TLorentzVector objects pointers
 * 
 */
double Mass(const std::vector<const TLorentzVector*>& particles)
{
    double sum{0};
    for (const auto& particle1 : particles)
    {
        for (const auto& particle2 : particles)
        {
            if (particle1 == particle2) { continue;} // Skip the same particle
            sum += particle1->Dot(*particle2);
        }
    }
    return sqrt(sum);
}

/**
 * @brief This function returns true if the charge of the particles is consistent with the region name.
 * @param q1
 * @param q2
 * @param regionName
 */
bool isChargeCorrect(const std::string& regionName, float q1, float q2){
    if (regionName.size() < 2){
        g_LOG(LogLevel::ERROR, "Region name is not valid.");
        exit(1);
    }

    bool isOS = regionName.find("OS") != std::string::npos;
    bool isSS = regionName.find("SS") != std::string::npos;

    if (isOS) return q1!=q2;
    if (isSS) return q1==q2;
    else {
        g_LOG(LogLevel::ERROR, "Region name is not valid. The region should end with OS or SS.");
        exit(1);
    }
}

/**
 * @brief This function calculates the transverse mass of a particle.
 * @param particle The particle to calculate the transverse mass for.
 * @param MET The missing transverse energy.
 * @return The transverse mass of the particle.
 */
double TransverseMass(const TLorentzVector* particle, const TLorentzVector* MET)
{
    double mt = sqrt(2*particle->Pt()*MET->Pt() * (1-cos(del_phi(particle->Phi(),MET->Phi()))));
    return mt;
}

/**
 * @brief Function to see if a jet is inside the rapidity gap spanned by two jets.
 * @param test_jet Jet to test.
 * @param j1 leading jet.
 * @param j2: subleading jet.
 */
int is_inside_jets(TLorentzVector* test_jet,TLorentzVector* j1, TLorentzVector* j2){
  double delta_y_j1j2 = abs(j1->Rapidity()-j2->Rapidity());
  double delta_y_j1test = abs(j1->Rapidity()-test_jet->Rapidity());
  double delta_y_j2test = abs(j2->Rapidity()-test_jet->Rapidity());
  if(delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2) return 0;
  return 1;
}

/**
 * @brief Function to calculate the minimum delta R between a test particle and a container of particles.
 * @param test_particle: particle to test.
 * @param bool_vector_container: container of booleans to select particles.
 * @param @param jet_container: container of particles to test against.
 */
double min_deltaR(TLorentzVector* test_particle, const std::vector<UInt_t>& bool_vector_container,const std::vector<TLorentzVector*>& jet_container){

  std::vector<double> delta_Rs{};

  for (size_t index{0};index<jet_container.size();index++){
    if (bool_vector_container[index]!=0){
      delta_Rs.push_back(jet_container[index]->DeltaR(*test_particle));
    }
    else {break;}
  }

  double min_dR=*std::min_element(delta_Rs.begin(),delta_Rs.end());
  return min_dR;
}

/**
 * @brief Get the number of jets in the gap defined by two tagging jets.
 * @param jetsVector The two first jets should be the tagging jets.
 * @param jetsBoolVector The boolean vector indicating if the jet is present in the event or not.
 * @return The number of jets in the gap.
 */
int getNumberOfGapJets(const std::vector<TLorentzVector*>& jetsVector, const std::vector<UInt_t>& jetsBoolVector){
    // First check that the two tagging jets are present in the event.
    if (jetsBoolVector[0] == 0 || jetsBoolVector[1] == 0){
        g_LOG(LogLevel::ERROR, "At least one of the two tagging jets are not present in the event.");
        exit(1);
    }

    // Count the number of jets in the gap.
    int nGapJets{0};
    for (size_t index{2};index<jetsVector.size();index++){
        if (jetsBoolVector[index] == 0) continue;
        nGapJets += is_inside_jets(jetsVector[index], jetsVector[0], jetsVector[1]);
    }
    return nGapJets;
}

/**
 * @brief Construct the pT balance variable.
 * @param particles The particle container.
 * @return The pT balance of the event.
 */
double getPtBalance(const std::vector<TLorentzVector*>& particles){
    double scalarSum{0};
    TLorentzVector vectorSum{};
    for (const auto& particle : particles){
        scalarSum += particle->Pt();
        vectorSum += *particle;
    }

    // Calculate the pT balance.
    double pt_bal= vectorSum.Pt()/scalarSum;
    return pt_bal;
}

/**
 * @brief Construct the signed-centrality variable.
 * @param taggingJet1 The leading tagging jet.
 * @param taggingJet2 The sub-leading tagging jet.
 * @param Tau The tau jet.
 * @param Lepton The lepton.
 * @return The centrality of the event.  
 */
double getSignedCentrality(const TLorentzVector* taggingJet1, const TLorentzVector* taggingJet2, const TLorentzVector* Tau, const TLorentzVector* Lepton){
    double lepton_xi = ((*Tau) + (*Lepton)).Rapidity();
    double dijet_xi = taggingJet1->Rapidity() + taggingJet2->Rapidity();
    double delta_y = (taggingJet1->Rapidity() - taggingJet2->Rapidity());
    double centrality = (lepton_xi - 0.5*dijet_xi)/delta_y;
    return centrality;
}

enum class Region
{
    DefaultNoRW,
    SR,
    CRa,
    CRb,
    CRc
};

/** 
 * @brief Get the region (QCDjj or EWKjj) of the event.
 * @param centrality The centrality of the event.
 * @param nJetsInGap The number of jets in the gap.
 * @return The region of the event.
 */
Region getRegion(double centrality, int nJetsInGap){
    Region region = Region::DefaultNoRW;
    if ((centrality<0.5 && centrality<=1) && nJetsInGap==0){region = Region::SR;}
    else if ((centrality<0.5 && centrality<=1) && nJetsInGap==1){region = Region::CRa;}
    else if ((centrality>=0.5 && centrality<=1) && nJetsInGap==1){region = Region::CRb;}
    else if ((centrality>=0.5 && centrality<=1) && nJetsInGap==0){region = Region::CRc;}
    return region;
}

} // namespace Kinematics