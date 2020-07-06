#define missing_cxx
#include "missing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>

string pdgid_tran(int pdg)
{
    int day = pdg;
    switch (day)
    {
    case 22:
        return "gamma";
    case -11:
        return "e+";
    case 11:
        return "e-";
    case 12:
        return "neutrino";
    case -13:
        return "mu+";
    case 13:
        return "mu-";
    case 111:
        return "pi0";
    case 211:
        return "pi+";
    case -211:
        return "pi-";
    case 130:
        return "K0_Long";
    case 321:
        return "Kaon+";
    case -321:
        return "Kaon-";
    case 2112:
        return "neutron";
    case 2212:
        return "proton";
    case -2212:
        return "antiproton";
    case 221:
        return "eta";
    case 3122:
        return "lambda";
    case 3222:
        return "sigma+";
    case 3212:
        return "sigma0";
    case 3112:
        return "sigma-";
    case -2112:
        return "AntiNeutron";
    case -15:
        return "Tau+";
    case 15:
        return "Tau-";
    case 1000010020:
        return "Deuteron";
    case 1000010030:
        return "Triton";
    case 1000020040:
        return "Alpha";
    }
    return "others";
}



void missing::Loop()
{
//   In a ROOT session, you can do:
//      root> .L missing.C
//      root> missing t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
TH2F *hmissing=new TH2F("hmissing","hmissing",700,-1,6, 700, -1, 6);
TH2F *hmissing2=new TH2F("hmissing2","hmissing2",700,-1,6, 700, -1, 6);

    TH1F *h_Elost=new TH1F("h_Elost","h_Elost",700,-1,6);
    
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   double Eksum=0;
   double Etsum=0;
   double Eava=0;
   double E_in=0;
   double E_out=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Eksum=0;
      Etsum=0;
      Eava=0;
      E_in=0;
      E_out=0;
      for(int i=0;i<primary_second_id->size();i++){
        Eksum+=primary_second_Ek->at(i);
        Etsum+=primary_second_Et->at(i);
        if(primary_second_id->at(i)==111 || primary_second_id->at(i)==211 || primary_second_id->at(i)==-211 || primary_second_id->at(i)==130 || primary_second_id->at(i)==321 ||primary_second_id->at(i)==-321 ) Eava+=primary_second_Et->at(i);
        else Eava+=primary_second_Ek->at(i);
        E_in += primary_second_mass->at(i); 
        E_out += primary_second_Et->at(i); 
      }
      E_in += ParentEtot1;
      E_out += ParentEtot2;
      
      string name;
      //if(Etsum<)
       //double Emaybedep=0;
       if(ParentEtot1!=-99){
      if(ParentEtot1 < Eava){
        cout<<"***********************************"<<endl;
        cout<<"** Etot_pre: "<<ParentEtot1<<" post: "<<ParentEtot2<<endl; 
        cout<<"Eloss: "<< ParentEtot1 - ParentEtot2<<endl;
        cout<<"** totEdep: "<<ParenttotEdep<<" non-ion: "<<ParentnonionEdep<<endl;
        cout<<"** secondary sum: Etot: "<<Etsum<<" Ek: "<<Eksum<<" Eavaliable: "<<Eava<<endl;
          //Emaybedep=0;
        for(int i=0;i<primary_second_id->size();i++){
                name = pdgid_tran(primary_second_id->at(i));
                //if(name=="others")
                //if(ParentEtot1 < Eava)
                //cout<<"id: "<<primary_second_id->at(i)<<" Ek: "<< primary_second_Ek->at(i)<<" Et: "<<primary_second_Et->at(i)<<endl;
                //else
                cout<<"id: "<<name<<" Ek: "<< primary_second_Ek->at(i)<<" Et: "<<primary_second_Et->at(i)<<endl;
            //if(primary_second_Ek->at(i)<ParenttotEdep) Emaybedep+=primary_second_Ek->at(i);
        }
          //cout<<"maybe dep: "<<Emaybedep<<endl;
      }
           if(ParentEK2==0) hmissing->Fill(ParentEtot1, Eava);
           else      hmissing->Fill(ParentEtot1 - ParentEtot2, Eava);
       hmissing2->Fill(ParentEtot1 - ParentEtot2, (ParentEtot1 - ParentEtot2 - Eksum) / (ParentEtot1 - ParentEtot2));
       //if(ParentEK2==0)
       h_Elost->Fill(E_in - E_out);
       //if(ParentEK2==0) h_Elost->Fill( (ParentEtot1 - Eava ) / (ParentEtot1));
       //else h_Elost->Fill( (ParentEtot1-ParentEtot2 - Eava ) / (ParentEtot1-ParentEtot2));
       //else h_Elost->Fill( (ParentEtot1-ParentEtot2 - Eksum) / (ParentEtot1-ParentEtot2));
       //h_Elost->Fill( (ParenttotEdep) / (ParentEtot1));
       }
   }
   hmissing->Draw("colz");
   TCanvas *c2=new TCanvas();
   hmissing2->Draw("colz");
   TCanvas *c3=new TCanvas();
   h_Elost->GetXaxis()->SetTitle("missing energy (GeV)");
   h_Elost->GetXaxis()->SetTitle("missing energy / energy loss");
   h_Elost->Draw("");

}

