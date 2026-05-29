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

using namespace Pythia8;
using namespace fastjet;

//setup jet finding settings
pair<JetDefinition, AreaDefinition> setupJet(double r,JetAlgorithm algorithm){
    double jetR = r;
    JetDefinition jetdef(algorithm,  //how to cluster
                         jetR,              //how wide jets are
                         E_scheme,          //how to combine 4-momenta in final jets
                         Best);             //how to run the algorithm efficiently
    
    
    AreaDefinition areaDefinition(active_area);

    return {jetdef, areaDefinition};
    }

//function for calculating theta using dot pdt
double tau(const PseudoJet &p1,const PseudoJet &p2){
    double px1 = p1.px(),py1 = p1.py(),pz1 = p1.pz();
    double px2 = p2.px(),py2 = p2.py(),pz2 = p2.pz();
    double mag1 = sqrt(pow(px1,2)+ pow(py1,2)+ pow(pz1,2));
    double mag2 = sqrt(pow(px2,2)+ pow(py2,2)+ pow(pz2,2));
    double dot = px1*px2 + py1*py2 + pz1*pz2;
    double cosTheta = dot/ (mag1 * mag2);
    
    //to hold the value inside [-1,1]
    cosTheta = max(-1., min(1., cosTheta));
    return (1. - cosTheta)/2.;
}
//creating a storage element for 5-D vector
struct JetFeatures {
    double tau1;
    double tau2;
    double tau3;
    double EEEC;
    double mt;
};

void run(int TargetJet, double mt, vector <JetFeatures>& dataset){
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
    double jetR = 1.2; //-> if we change in setupJet() we change here too
    auto [jetdef_antikt, areadef_antikt] = setupJet(jetR, antikt_algorithm);
    double jetRR = 0.1;
    auto [jetdef_AC, areadef_AC] = setupJet(jetRR, cambridge_aachen_algorithm);
    
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
        if (totalE < 500.) continue;  // skip entire event
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
        vector<PseudoJet>jets = sorted_by_pt(cs_antikt.inclusive_jets(10.));

//        cout<<"Number of jets: " << jets.size() << endl;
        for (auto & jet : jets) {
            //cout << "Multiplicity: " << jet.constituents().size()<< endl;
            if (jet.e() < 100.0) continue;
            double jet_mass = jet.m();
            
            
            if (jet_mass < 100. || jet_mass > 250.) continue;  // ← top mass window
            cout << "jet mass = " << jet_mass << endl;
            
            vector<PseudoJet> constituents = jet.constituents();

            //Reclustering using CA
            ClusterSequence cs_ac(constituents, jetdef_AC);
            vector<PseudoJet>subjets = sorted_by_pt(cs_ac.inclusive_jets(1e-15));
            
//            cout << "Subjet multiplicity: "   << subjets.size() <<endl;
            
            int multiplicity = 0;
            for(int k=0; k<subjets.size(); k++){
                for(int l=k+1; l<subjets.size(); l++){
                    for(int m=l+1; m<subjets.size(); m++){
                        
                        //EEC
                        double E1 = subjets[k].e();
                        double E2 = subjets[l].e();
                        double E3 = subjets[m].e();
                        
                        double EEEC = E1*E2*E3/(pow(CM, 3));
                        if(EEEC < 1e-6) continue;

                        //pairwise angular separation theta of particles
                        double tau1 = tau(subjets[k],subjets[l]);
                        double tau2 = tau(subjets[l],subjets[m]);
                        double tau3 = tau(subjets[m],subjets[k]);
                        //sorting in ascending order                        
                        vector<double> t = {tau1, tau2, tau3};
                        sort(t.begin(), t.end());
                        tau1 = t.at(0); tau2 = t.at(1); tau3 = t.at(2);
                        if (tau1 > 0.1 || tau2 > 0.1 || tau3 > 0.1) continue;
                        
                        JetFeatures jf;
                        jf.tau1 = tau1; jf.tau2 = tau2; jf.tau3 = tau3;
                        jf.EEEC = EEEC; jf.mt = mt;
                        dataset.push_back(jf);

                        multiplicity+= 1;
                    }
                }
            }
//            cout<< "EEE tuples: " << multiplicity << endl;
            jetcount++;
            if(jetcount >= TargetJet) break;
        }
    }
    cout << "Done! Generated " << iEvent << " events for " 
         << jetcount << " jets\n";
    
}

int main(){
    
    int N = 1000;
    
    for (double mt = 170.0; mt <= 171.0; mt += 0.1){    
        //rounding off
        mt = round(mt * 10.0)/10.0;

        //storage
        vector<JetFeatures> dataset;
        run(N, mt, dataset);

        //Clean name
        ostringstream oss;
        oss << fixed << setprecision(1) << mt;
        ofstream logfile("/home/sawini-jana/hep/Paper/dataset_"+ oss.str() + ".txt");
        
        //loop through dataset
         
        
        for (const auto& jet : dataset) {

            logfile << jet.tau1 << " "
                    << jet.tau2 << " "
                    << jet.tau3 << " "
                    << jet.EEEC << " "
                    << jet.mt
                    << "\n";
        }

        logfile.close();
        cout << "\nTotal EEEC tuples for mt = " << oss.str()
             << " : " << dataset.size() << "\n\n";    
    }

    return 0;
}