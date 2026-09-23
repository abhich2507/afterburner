#include "Histogrammer.hh"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>


Histogrammer::Histogrammer(std::shared_ptr<TFile> file, std::string dir_name):
        _file(std::move(file)),
        _dirname(std::move(dir_name))
{
}

void Histogrammer::initialize() {
    // creates ROOT file to write things to
    //_file = std::make_shared<TFile>(_output_file_path.c_str(), "RECREATE");
    _file->cd();
    auto dir = _file->mkdir(_dirname.c_str());
    dir->cd();
    // Particle Quantities
    partPtHist = new TH1D("partPt","Final State Particle Pt",500,0.,50.);
    partEtaHist = new TH1D("partEta","Final State Particle Eta",400,-10.,10.);
    partPhiHist = new TH1D("partPhi","Final State Particle Phi",100,-1.0*TMath::Pi(),TMath::Pi());
    partPtVsEtaHist = new TH2D("partPtVsEta","Final State Particle Pt Vs Eta",400,-10.,10.,500,0.,50.);
    partPhiVsEtaHist = new TH2D("partPhiVsEta","Final State Particle Phi Vs Eta",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());
    /*
    // Beam Shape
  TH1D *eCM = new TH1D("eCM","Modified - Nominal CM Energy",10000,-1.0,1.0);
  TH2D *pXY1 = new TH2D("pXY1","Hadron Beam Py Vs Px",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd1 = new TH2D("pXZProd1","Hadron Beam Px Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pYZProd1 = new TH2D("pYZProd1","Hadron Beam Py Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pXY2 = new TH2D("pXY2","Lepton Beam Py Vs Px",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd2 = new TH2D("pXZProd2","Lepton Beam Px Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pYZProd2 = new TH2D("pYZProd2","Lepton Beam Py Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH1D *pZ1 = new TH1D("pZ1","Hadron Beam Pz",10000,0.0,280.0);
  TH1D *pZ2 = new TH1D("pZ2","Lepton Beam Pz",10000,-20.,0.);

     TH2D *lepVsHadPartZ = new TH2D("lepVsHadPartZ","Intrabunch Z Positions of Colliding Leptons Vs Hadrons",5000,-500.0,500.0,5000,-500.0,500.0);
*/
    vtxX = new TH1D("vtxX","Vertex x;[mm]",100,-5.0,5.0);
    vtxY = new TH1D("vtxY","Vertex y;[mm]",100,-.05,.05);
    vtxZ = new TH1D("vtxZ","Vertex z;[mm]",100,-400.0,400.0);
    vtxT = new TH1D("vtxT","Time;[mm]",100,-400.0,400.0);
    vtx2X = new TH1D("vtx2X","Vertex x;[mm]",500,-5.0,5.0);
    vtx2Y = new TH1D("vtx2Y","Vertex y;[mm]",500,-5.0,5.0);
    vtx2Z = new TH1D("vtx2Z","Vertex z;[mm]",500,-500.0,500.0);
    vtx2T = new TH1D("vtx2T","Time;[mm]",500,-500.0,500.0);
    vtxYvsX = new TH2D("vtxYvsX","Vertex Y vs X;X [mm];Y [mm]",5000,-5.0,5.0,5000,-5.0,5.0);
    vtxXvsT = new TH2D("vtxXvsT","Vertex X vs T;T [mm];X [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
    vtxXvsZ = new TH2D("vtxXvsZ","Vertex X vs Z;Z [mm];X [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
    vtxYvsZ = new TH2D("vtxYvsZ","Vertex Y vs Z;Z [mm];Y [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
    vtxTvsZ = new TH2D("vtxTvsZ","Interaction Time Vs Z-vertex;Z [mm];T [mm]",5000,-500.0,500.0,5000,-500.0,500.0);
    vtxXvsTZSum = new TH2D("vtxXvsTZSum","Vertex X vs T+Z;T+Z [mm];X [mm]",5000,-500.,500.,5000,-5.,5.);
    vtxXvsTZDiff = new TH2D("vtxXvsTZDiff","Vertex X vs T-Z;T-Z [mm];X [mm]",5000,-500.,500.,5000,-5.,5.);

    atan2PxPz1Hist   = new TH1D("atan2PxPz1","",2500,-0.05,0.05);
    atan2PyPz1Hist   = new TH1D("atan2PyPz1","",2500,-0.01,0.01);
    atan2PyPtot1Hist = new TH1D("atan2PyPtot1","",500,-0.001,0.001);
    atan2PxPz2Hist   = new TH1D("atan2PxPz2","",2500,-0.01,0.01);
    atan2PyPz2Hist   = new TH1D("atan2PyPz2","",2500,-0.01,0.01);
}


void Histogrammer::process_event(HepMC3::GenEvent &event) {

    std::vector<std::shared_ptr<const HepMC3::GenParticle>> beam_particles;

    for(auto& prt: event.particles()) {
        // Select beam particles for further analysis
        if(prt->status() == 4) {
            beam_particles.push_back(prt);
        }

        bool partFin = prt->status() == 1;
        double partPt = prt->momentum().pt();
        double partEta = prt->momentum().eta();
        double partPhi = prt->momentum().phi();

        double px =  prt->momentum().px();
        double py =  prt->momentum().py();
        double pz =  prt->momentum().pz();
        double E =   prt->momentum().e();

        //cout << i << " " << p8.event[i].id() << " " << partPt << " " << partEta << " " << partPhi << endl;

        if(partFin && partEta>-10.0 && partEta<10.0 ) // && y<0.95 && y>0.01 && i > 7)
        {
            partPtHist->Fill(partPt);
            partPhiHist->Fill(partPhi);
            partEtaHist->Fill(partEta);

            partPtVsEtaHist->Fill(partEta,partPt);
            partPhiVsEtaHist->Fill(partEta,partPhi);

//            partStatusVsEtaHist->Fill(partEta,p8.event[i].status());
//
//            if(partPt > 1.0)
//            {
//                partPtHiHist->Fill(partPt);
//                partPhiHiHist->Fill(partPhi);
//                partEtaHiHist->Fill(partEta);
//
//                partPhiVsEtaHiHist->Fill(partEta,partPhi);
//
//                partStatusVsEtaHiHist->Fill(partEta,p8.event[i].status());
//            }
        }
    }

    // Beam effects plots



    // Check particles
    if(beam_particles.size() != 2) {
        std::stringstream msg("Input file should have exactly 2 beam particles (status code = 4). File has ");
        msg<<beam_particles.size();
        throw std::runtime_error(msg.str());
    }

    // --- identify beams by species, not by record order ---
    auto is_lepton = [](int pdg) { return std::abs(pdg) == 11 || std::abs(pdg) == 13; };

    std::shared_ptr<const HepMC3::GenParticle> lepton = nullptr;
    std::shared_ptr<const HepMC3::GenParticle> hadron = nullptr;
    for (const auto &bp : beam_particles) {
        if (is_lepton(bp->pdg_id())) lepton = bp;
        else                        hadron = bp;
    }
    if (!lepton || !hadron) {
        throw std::runtime_error("Could not identify one lepton and one hadron beam (status==4)");
    }

    auto mom_had = hadron->momentum();
    auto mom_lep = lepton->momentum();
    auto vtx_had = hadron->end_vertex();
    auto vtx_lep = lepton->end_vertex();
    // Note here, for many HepMC3 files these vertices might be 2 different objects.
    // We assume here that their x,y,z are identical

    if(_verbose && vtx_had && vtx_lep) {
        printf("had vtx %7.1f %7.1f %7.1f    lep vtx %7.1f %7.1f %7.1f\n",
               vtx_had->position().x(), vtx_had->position().y(), vtx_had->position().z(),
               vtx_lep->position().x(), vtx_lep->position().y(), vtx_lep->position().z());

    }

    if(vtx_had) {
        vtxX->Fill(vtx_had->position().x());
        vtxY->Fill(vtx_had->position().y());
        vtxZ->Fill(vtx_had->position().z());
        vtxT->Fill(vtx_had->position().t());
    } else {
        vtxX->Fill(0.0);
        vtxY->Fill(0.0);
        vtxZ->Fill(0.0);
        vtxT->Fill(0.0);
    }

    if(vtx_lep) {
        vtx2X->Fill(vtx_lep->position().x());
        vtx2Y->Fill(vtx_lep->position().y());
        vtx2Z->Fill(vtx_lep->position().z());
        vtx2T->Fill(vtx_lep->position().t());
    } else {
        vtx2X->Fill(0.0);
        vtx2Y->Fill(0.0);
        vtx2Z->Fill(0.0);
        vtx2T->Fill(0.0);
    }

    // --- direction-independent divergence angles ---
    // atan2(p_t, |pz|) gives the small angle w.r.t. the beam axis regardless of +z/-z travel
    auto ang = [](double p_t, double pz) { return TMath::ATan2(p_t, std::abs(pz)); };

    atan2PxPz1Hist  ->Fill(ang(mom_had.px(), mom_had.pz()));
    atan2PyPz1Hist  ->Fill(ang(mom_had.py(), mom_had.pz()));
    atan2PyPtot1Hist->Fill(TMath::ATan2(mom_had.py(), mom_had.length()));

    atan2PxPz2Hist  ->Fill(ang(mom_lep.px(), mom_lep.pz()));
    atan2PyPz2Hist  ->Fill(ang(mom_lep.py(), mom_lep.pz()));

//    vtxYvsX->Fill(p8.process[0].xProd(),p8.process[0].yProd());
//    vtxXvsT->Fill(p8.process[0].tProd(),p8.process[0].xProd());
//    vtxXvsZ->Fill(p8.process[0].zProd(),p8.process[0].xProd());
//    vtxYvsZ->Fill(p8.process[0].zProd(),p8.process[0].yProd());
//    vtxTvsZ->Fill(p8.process[0].zProd(),p8.process[0].tProd());
//
//    vtxXvsTZSum->Fill(p8.process[0].tProd()+p8.process[0].zProd(),p8.process[0].xProd());
//    vtxXvsTZDiff->Fill(p8.process[0].tProd()-p8.process[0].zProd(),p8.process[0].xProd());
//
//    double hadZ = p8.process[0].zProd() - TMath::Cos(0.0125)*p8.process[0].tProd();
//    double lepZ = p8.process[0].zProd() + TMath::Cos(0.0125)*p8.process[0].tProd();
//    lepVsHadPartZ->Fill(hadZ,lepZ);
}


Histogrammer::~Histogrammer() {

}

