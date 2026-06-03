#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include "fastjet/contrib/Nsubjettiness.hh"

using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;

// One entry PER SUBJET, grouped by jet:
struct SubjetFeatures {
    double E, px, py, pz, eta, phi, pt;
};

struct JetGraph {
    vector<SubjetFeatures> subjets;  // nodes
    double mt;                        // label
    int jet_id;                       // which jet this belongs to
    double tau21;
    double tau32;
    double tau43;  
    int multiplicity;  
};

void run(int TargetJet, double mt, vector <JetGraph>& dataset){
    //setup pythia for generating events
    Pythia pythia;

    //switch on processes
    int beamA = 11; //proton-proton collision
    int beamB = -11;
    double CM = 2000.;

    // e+ e- collision
    pythia.readString("Beams:idA = " + to_string(beamA));
    pythia.readString("Beams:idB = " + to_string(beamB));
    pythia.readString("Beams:eCM = " + to_string(CM));
    pythia.readString("6:m0 = " + to_string(mt));

    // Enable e+e- → γ*/Z, no 
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); 

    // Force Z/γ* → t tbar
    pythia.readString("23:onIfMatch = 6 -6");

    //forcing hadronic decays of W
    pythia.readString("24:onMode = off");
    pythia.readString("24:onIfAny = 1 2 3 4 5");//u,d,s,c,b

    //layers
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:FSR = on");

    pythia.readString("HadronLevel:Hadronize = on");

    //Unique seed
    int randomSeed = 0; //-> output for multiple runs are not same for 0
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + to_string(randomSeed));

    //supress event-by-event for cleaner output
    pythia.readString("Print:quiet = on");

    pythia.readString("Init:showAllSettings = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showAllParticleData = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");

    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    pythia.init();

    //setup jet finder
    double jetR = 1.2; 
    double jetRR = 0.22;
    JetDefinition jetdef_antikt(ee_genkt_algorithm,  //how to cluster
                            jetR,  
                            -1,            //how wide jets are
                            E_scheme);          //how to combine 4-momenta in final jets
                                         //how to run the algorithm efficiently
        
    JetDefinition jetdef_AC(ee_genkt_algorithm,
                            jetRR,
                            0.0,
                            E_scheme);
    int jetcount = 0;
    int iEvent = 0;

    //event loop
    while (jetcount < TargetJet){
        iEvent++;

        if(!pythia.next()) continue;
//        cout<<"Generating event: "<< iEvent<< " Event size: "<<  pythia.event.size() <<"\n";

        //collecting particles for FastJet
        vector<PseudoJet> fjParticles;
        double totalE = 0;
        for (int j = 0; j < pythia.event.size(); j++) {
            if (!pythia.event[j].isFinal()) continue;
            if (!pythia.event[j].isVisible()) continue;
            totalE += pythia.event[j].e();

        }

        // Check BEFORE collecting particles:
        if (totalE < 1500.) continue;  // skip entire event
        //particle loop
        for (int j = 0; j < pythia.event.size(); j++) {
            //TO note: pythia.event[j] -> datatype = Particle& 
            if (!pythia.event[j].isFinal()) continue;
            if (!pythia.event[j].isVisible()) continue;
           
            //to handle the warning
            if (pythia.event[j].pT() < 1e-15) continue;

            double eta = pythia.event[j].eta();

            double px = pythia.event[j].px();
            double py = pythia.event[j].py();
            double pz = pythia.event[j].pz();
            double E = pythia.event[j].e();

            //here is where we are linking-> particle properties in fastjet objects
            fjParticles.push_back(PseudoJet(px, py, pz, E));
        }
//        cout<< endl;

        //Clustering using anti-kt
        ClusterSequence cs_antikt(fjParticles, jetdef_antikt);
        vector<PseudoJet>jets = sorted_by_pt(cs_antikt.inclusive_jets(2));

//        cout<<"Number of jets: " << jets.size() << endl;
        for (auto & jet : jets) {
            //cout << "Multiplicity: " << jet.constituents().size()<< endl;
            if (jet.e() < 300.0) continue;
            double jet_mass = jet.m();
            
            vector<PseudoJet> constituents = jet.constituents();

            //Reclustering using CA
            ClusterSequence cs_ac(constituents, jetdef_AC);
            vector<PseudoJet>subjets = sorted_by_pt(cs_ac.inclusive_jets(1e-10));
//            cout << "Subjet multiplicity before: "   << subjets.size() <<endl;
            
// Then per jet (not per subjet):
            double beta = 1.0;
            Nsubjettiness tau1_calc(1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness tau2_calc(2, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness tau3_calc(3, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness tau4_calc(4, OnePass_KT_Axes(), UnnormalizedMeasure(beta));

            double tau1 = tau1_calc(jet);
            double tau2 = tau2_calc(jet);
            double tau3 = tau3_calc(jet);
            double tau4 = tau4_calc(jet);

            double tau21 = tau2/tau1;  // 2-prong vs 1-prong
            double tau32 = tau3/tau2;  // 3-prong vs 2-prong
            double tau43 = tau4/tau3;
            
            int n_subjets = subjets.size();
            // Remove soft subjets that create large ζ values:
            double jet_e = jet.e();
            vector<PseudoJet> hard_subjets;
            for (auto& sj : subjets) {
                if (sj.e()/jet_e >= 0.001)
                    hard_subjets.push_back(sj);
            }
            subjets = hard_subjets;
//            cout<< "EEE tuples: " << multiplicity << endl;

            JetGraph jg;
            jg.mt = mt;
            jg.jet_id = jetcount;
            jg.tau43 = tau43;
            jg.tau32 = tau32;
            jg.tau21 = tau21;
            jg.multiplicity = n_subjets;

            for (auto& sj : hard_subjets) {
                
                SubjetFeatures sf;
                sf.E   = sj.e();
                sf.px  = sj.px();
                sf.py  = sj.py();
                sf.pz  = sj.pz();
                sf.eta = sj.eta();
                sf.phi = sj.phi();
                sf.pt  = sj.pt();

                jg.subjets.push_back(sf);
            }

            dataset.push_back(jg);

            jetcount++;
            if(jetcount >= TargetJet) break;
//            cout << "Subjet multiplicity after: "   << subjets.size() <<endl;
        }
    }
    cout << "Done! Generated " << iEvent << " events for " 
         << jetcount << " jets\n";
    
}

int main(){
    
    int N = 10000;
    
    for (double mt = 170.0; mt <= 180.0; mt += 0.5){    
        //rounding off
        mt = round(mt * 10.0)/10.0;

        //storage
        vector<JetGraph> dataset;
        run(N, mt, dataset);

        //Clean name
        ostringstream oss;
        oss << fixed << setprecision(1) << mt;
        ofstream logfile("/home/sawini-jana/hep/Paper/datasetforGNN_"+ oss.str() + ".txt");
        
        //loop through dataset                 
        for (const auto& jg : dataset) {
            // Write all subjets for this jet:
            for (const auto& sf : jg.subjets) {
                logfile << jg.jet_id<< " "   // which jet
                        << sf.E    << " "
                        << sf.px   << " "
                        << sf.py   << " "
                        << sf.pz   << " "
                        << sf.eta  << " "
                        << sf.phi  << " "
                        << sf.pt   << " "
                        << jg.tau43<< " "
                        << jg.tau32<< " "
                        << jg.tau21<< " "
                        << jg.multiplicity<< " "
                        << jg.mt               // label
                        << "\n";
            }
        }

        logfile.close();
        cout << "\nTotal EEEC tuples for mt = " << oss.str()
             << " : " << dataset.size() << "\n\n";    
    }

    return 0;
}
