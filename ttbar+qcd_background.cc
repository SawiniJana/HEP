#include "TFile.h"
#include "TTree.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "TClonesArray.h"

#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

struct Observables {
    double phi1, phi2, phi3, phi4;
    int label;
};
void processfile(const string& filename, int label, int max_jets, ofstream &fout){ //max_jets is what we set
    //label: to distinguish QCD from top jet

    TFile* f = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)f->Get("Delphes");
    ExRootTreeReader* reader = new ExRootTreeReader(tree);
    Long64_t nEntries = reader->GetEntries();

    TClonesArray* branchEFT  = reader->UseBranch("EFlowTrack");
    // charged hadrons + matched calo
    TClonesArray* branchEFP  = reader->UseBranch("EFlowPhoton");
    //photons meaning only calo, no tracks
    TClonesArray* branchEFN  = reader->UseBranch("EFlowNeutralHadron");
    //Neutral hadrons, meaning only calo, no track
    
    JetDefinition jet_def(antikt_algorithm, 0.8);
    SoftDrop sd(0.0, 0.1); // beta=0, zcut=0.1

    int iEvent = 0;
    int accepted = 0;

    for (Long64_t ev = 0; ev < nEntries; ev++) {
        //A cut
        reader->ReadEntry(ev);

        // Collect all EFlow (particle-flow) objects as FastJet inputs
        vector<PseudoJet> particles;

        for (int i = 0; i < branchEFT->GetEntries(); i++) {
            Track* t = (Track*)branchEFT->At(i);

            double pt   = t->PT;
            double eta  = t->Eta;
            double phi  = t->Phi;
            double mass = t->Mass;

            double px = pt*cos(phi);
            double py = pt*sin(phi);
            double pz = sinh(eta) * sqrt(pow(pt,2) + pow(mass,2));
            double E = cosh(eta) * sqrt(pow(pt,2) + pow(mass,2));
            
            /*
            std::cout
            << " E = " << E
            << " Px = " << px
            << " Py = " << py
            << " Pz = " << pz
            << std::endl;
            */
            particles.push_back(PseudoJet(px,py,pz,E));
        }

        for (int i = 0; i < branchEFP->GetEntries(); i++) {
            Tower* t = (Tower*)branchEFP->At(i);
            TLorentzVector p4 = t->P4(); //Constructing from 4 momentum
            
            double px = p4.Px();
            double py = p4.Py();
            double pz = p4.Pz();
            double E = p4.E();
            particles.push_back(PseudoJet(px,py,pz,E));
        }

        for (int i = 0; i < branchEFN->GetEntries(); i++) {
            Tower* t = (Tower*)branchEFN->At(i);
            TLorentzVector p4 = t->P4(); //Constructing from 4 momentum
            
            double px = p4.Px();
            double py = p4.Py();
            double pz = p4.Pz();
            double E = p4.E();
            particles.push_back(PseudoJet(px,py,pz,E));
        }

        if (particles.empty()) continue;
        
        /*
        std::cout << "Event " << ev
          << " particles = "
          << particles.size()
          << std::endl;
        */

        iEvent++;
        if (iEvent >= max_jets) break;

        ClusterSequence cs(particles, jet_def);
        auto jets = sorted_by_pt(cs.inclusive_jets());
        //cout<< "Jet Size: " << jets.size() << endl;

        int nboosted = 0;

        for (auto& jet : jets) {   
            //cout<< "Jet momentum: "<< jet.pt() << endl;     
            if (jet.pt() < 300.0) continue;
            // Soft-drop mass
            PseudoJet sd_jet = sd(jet);
            double msd = sd_jet.m();
            if (msd < 140.0 || msd > 220.0) continue;
            nboosted++;

            auto constituents = sd_jet.constituents();
            sort(constituents.begin(), constituents.end(), [](const PseudoJet &a, const PseudoJet &b){
                return a.pt() > b.pt();
            });
            if (constituents.size() > 25) constituents.resize(25);
            //cout<<"Event: "<< iEvent << " Size of constituents: "<< constituents.size() <<endl;
            //cout<< "Jet momentum: "<< jet.pt() << endl;
            //cout << "Soft-Drop Mass: " << msd << endl;

            //double totalMomentum = 0.0;
            double totalEnergy = 0.0;
            
            accepted++;

            for (const auto &c : constituents)
                totalEnergy += c.e();
                //totalMomentum += c.pt();
        
            double EECnarrow = 0.0;
            double EECwide = 0.0;
            double e3 = 0.0;
            double e2 = 0.0;

            for(int k=0; k<constituents.size(); k++){
                double E1 = constituents[k].e()/totalEnergy;
                
                for(int l=k+1; l<constituents.size(); l++){
                    double E2 = constituents[l].e()/totalEnergy;
                    
                    //pairwise angular separation theta of particles
                    double dR1 = constituents[k].delta_R(constituents[l]);
                    if (0.06 < dR1 && dR1 <= 0.20)
                        EECnarrow += E1 * E2;
                    if (0.20 < dR1 && dR1 <= 0.80)
                        EECwide += E1 * E2;
                    
                    e2 += E1*E2 * dR1;

                    for(int m=l+1; m<constituents.size(); m++){
                        
                        double E3 = constituents[m].e()/totalEnergy;
                        double dR2 = constituents[k].delta_R(constituents[m]);
                        double dR3 = constituents[m].delta_R(constituents[l]);

                        e3 += E1*E2*E3 * dR1*dR2*dR3;          
                    }
                }
            }    

            Observables Obs;
            Obs.phi1 = EECnarrow;
            Obs.phi2 = EECwide;
            Obs.phi3 = e2;
            Obs.phi4 = e3;
            Obs.label = label;
            fout<<Obs.phi1<<","<<Obs.phi2<<","<<Obs.phi3<<","<<Obs.phi4<<","
                <<Obs.label<<"\n";
        }
        //if (nboosted != 0)
        //cout<< "Event: " << iEvent <<" Total nboosted: "<<nboosted<< endl<< endl;
    
    }
    cout<< "Total accepted jets: "<< accepted<< endl;
}

int main() {

    ofstream fout("observables.csv");
    fout << "phi1,phi2,phi3,phi4,label\n";

    const string filename1 = "/home/sawini-jana/hep/MG5_aMC_v3_7_1/qcd_output.root";
    const string filename2 = "/home/sawini-jana/hep/MG5_aMC_v3_7_1/ttbar_output.root";
    
    cout<< "Ttbar output: " << endl;
    processfile(filename2, 1, 50000, fout);
    cout<< "QCD output: " << endl;
    processfile(filename1, 0, 50000, fout);

    fout.close();

    return 0;
}
