//
//  test.cxx
//  SST test beam 2017 analysis
//
//  Created by Daniil Sukhonos on 21.01.20.
//

#include <stdio.h>

#include <iostream>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>

using namespace std;

int main() {
	gStyle->SetOptStat(1111111);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetLabelSize(0.045,"y");
	gStyle->SetLabelSize(0.045);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.3);
	gStyle->SetStatStyle(0);
	gStyle->SetTitleFontSize(0.060);

	TFile inputFile("sept2017_run_222_mod.root","update");
	TTree * tree;
	inputFile.GetObject("tree",tree);
	Int_t ncluster;
	Bool_t ClustIsOnTrk[600];
	tree->SetBranchAddress("ncluster", &ncluster);
	tree->SetBranchAddress("ClustIsOnTrk", ClustIsOnTrk);
	Int_t nclusterIsOnTrk, nclusterNotOnTrk;
	branch_1 = tree->Branch("nclusterIsOnTrk", &nclusterIsOnTrk, "nclusterIsOnTrk/I");
	branch_2 = tree->Branch("nclusterNotOnTrk", &nclusterNotOnTrk, "nclusterNotOnTrk/I");
	nEntries = tree->GetEntries();
	nclusterIsOnTrk = 0;
	nclusterNotOnTrk = 0;

	for (int i = 0; i < nEntries; ++i) {
		if (i%10000 == 0) cout<<"Event:\t"<<i<<"\r"<<flush;
		tree->GetEntry(i);
		if (ncluster > 0) {
			for (int j = 0; j < ncluster; ++j) {
				if (ClustIsOnTrk[j]) ++nclusterIsOnTrk;
				else ++nclusterNotOnTrk;
			}
			branch_1->Fill();
			nclusterIsOnTrk = 0;
			branch_2->Fill();
			nclusterNotOnTrk = 0;
		}
	}

	tree->Print();

	tree->Write();
	inputFile.Write();
	inputFile.Close();
	return 0;
}

// ########################################################################################

#include <iostream>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>

using namespace std;

TH1I ** test() {

	TFile inputFile("sept2017_run_222_mod.root","read");
	TTree * tree;
	inputFile.GetObject("tree",tree);
	Int_t ncluster;
	Bool_t ClustIsOnTrk[600];
	tree->SetBranchAddress("ncluster", &ncluster);
	tree->SetBranchAddress("ClustIsOnTrk", ClustIsOnTrk);
	Int_t nclusterIsOnTrk, nclusterNotOnTrk;
	nEntries = tree->GetEntries();
	nclusterIsOnTrk = 0;
	nclusterNotOnTrk = 0;

	TH1I ** hist_arr = new TH1I*[7];
	TH1I * h_clust_not_on_trk = new TH1I("h_clust_not_on_trk","nclusterNotOnTrk",120,0,120);
	hist_arr[0] = h_clust_not_on_trk;
	TH1I * h_clust_is_on_trk = new TH1I("h_clust_is_on_trk","nclusterIsOnTrk",120,0,120);
	hist_arr[1] = h_clust_is_on_trk;
	TH1I * h_clust_not_on_trk_ntrk_0 = new TH1I("h_clust_not_on_trk_ntrk_0","nclusterNotOnTrk {ntrk == 0}",120,0,120);
	hist_arr[2] = h_clust_not_on_trk_ntrk_0;
	TH1I * h_clust_is_on_trk_ntrk_0 = new TH1I("h_clust_is_on_trk_ntrk_0","nclusterIsOnTrk {ntrk == 0}",120,0,120);
	hist_arr[3] = h_clust_is_on_trk_ntrk_0;
	TH1I * h_clust_is_on_trk_ntrk_m1 = new TH1I("h_clust_is_on_trk_ntrk_m1","nclusterIsOnTrk {ntrk >= 1}",120,0,120);
	hist_arr[4] = h_clust_is_on_trk_ntrk_m1;
	TH1I * h_clust_not_on_trk_ntrk_m1 = new TH1I("h_clust_not_on_trk_ntrk_m1","nclusterNotOnTrk {ntrk >= 1}",120,0,120);
	hist_arr[5] = h_clust_not_on_trk_ntrk_m1;
	TH1I * h_clust_not_on_trk_ntrk_1 = new TH1I("h_clust_not_on_trk_ntrk_1","nclusterNotOnTrk {ntrk == 1}",120,0,120);
	hist_arr[6] = h_clust_not_on_trk_ntrk_1;

	Int_t ntrk;
	tree->SetBranchAddress("ntrk", &ntrk);
	Int_t nclusterIsOnTrk_ntrk_0(0), nclusterIsOnTrk_ntrk_m1(0), nclusterNotOnTrk_ntrk_0(0);
	Int_t nclusterNotOnTrk_ntrk_1(0), nclusterNotOnTrk_ntrk_m1(0);
	Int_t bad_tracks(0), good_tracks(0), total_tracks(0);

	//Int_t ClustLayer[600];
	//Double_t ClustPos[600];
	//Double_t chi2[100];
	//Double_t ndof[100];
	//Double_t ItpX[100];
	//Double_t ItpY[100];
	//Double_t SlpX[100];
	//Double_t SlpY[100];
	//tree->SetBranchAddress("ClustLayer", ClustLayer);
	//tree->SetBranchAddress("ClustPos", ClustPos);
	//tree->SetBranchAddress("chi2", chi2);
	//tree->SetBranchAddress("ndof", ndof);
	//tree->SetBranchAddress("ItpX", ItpX);
	//tree->SetBranchAddress("ItpY", ItpY);
	//tree->SetBranchAddress("SlpX", SlpX);
	//tree->SetBranchAddress("SlpY", SlpY);

	for (int i = 0; i < nEntries; ++i) {
		if (i%10000 == 0) cout<<"Event:\t"<<i<<"\r"<<flush;
		tree->GetEntry(i);
		if (ncluster > 0) {
			for (int j = 0; j < ncluster; ++j) {
				if (ClustIsOnTrk[j]) {
					++nclusterIsOnTrk;
					if (ntrk == 0) ++nclusterIsOnTrk_ntrk_0;
					else ++nclusterIsOnTrk_ntrk_m1;
				}
				else {
					++nclusterNotOnTrk;
					if (ntrk == 0) ++nclusterNotOnTrk_ntrk_0;
					else {
						++nclusterNotOnTrk_ntrk_m1;
						if (ntrk == 1) ++nclusterNotOnTrk_ntrk_1;
					}
				}
			}
	//		if (ntrk == 1 && nclusterIsOnTrk == 4 && ncluster == 8) {
	//			cout << "Event: " << i << endl;
	//			for (int j = 0; j < ncluster; ++j) {
	//				cout << "Cluster No: " << j << endl;
	//				cout << "Is on trk: " << ClustIsOnTrk[j] << endl;
	//				cout << "ClustLayer: " << ClustLayer[j] << endl;
	//				cout << "ClustPos: " << ClustPos[j] << endl;
	//				cout << "!!!!!!!!!!!" << endl;
	//			}
	//			cout << "Tracking info ......... " << endl;
	//			cout << "chi2: " << chi2[0] << endl;
	//			cout << "ndof: " << ndof[0] << endl;
	//			cout << "ItpX: " << ItpX[0] << endl;
	//			cout << "ItpY: " << ItpY[0] << endl;
	//			cout << "SlpX: " << SlpX[0] << endl;
	//			cout << "SlpY: " << SlpY[0] << endl;
	//			cout << "........................" << endl;
	//			cin.get();
	//		}
			if (ntrk == 1) {
				++total_tracks;
				if (nclusterIsOnTrk < 8) {
					++bad_tracks;
				}
				else {
					++good_tracks;
				}
			}
			if (nclusterIsOnTrk_ntrk_m1 > 0) {
				h_clust_is_on_trk_ntrk_m1->Fill(nclusterIsOnTrk_ntrk_m1);
				nclusterIsOnTrk_ntrk_m1 = 0;
			}
			h_clust_not_on_trk->Fill(nclusterNotOnTrk);
			h_clust_is_on_trk->Fill(nclusterIsOnTrk);
			h_clust_not_on_trk_ntrk_0->Fill(nclusterNotOnTrk_ntrk_0);
			h_clust_is_on_trk_ntrk_0->Fill(nclusterIsOnTrk_ntrk_0);
			h_clust_not_on_trk_ntrk_m1->Fill(nclusterNotOnTrk_ntrk_m1);
			h_clust_not_on_trk_ntrk_1->Fill(nclusterNotOnTrk_ntrk_1);
			nclusterIsOnTrk = 0;
			nclusterIsOnTrk_ntrk_0 = 0;
			nclusterNotOnTrk = 0;
			nclusterNotOnTrk_ntrk_0 = 0;
			nclusterNotOnTrk_ntrk_m1 = 0;
			nclusterNotOnTrk_ntrk_1 = 0;
		}
	}
	cout << "Bad tracks ratio: " << bad_tracks/(1.0 * total_tracks) << endl;
	cout << "Good tracks ratio: " << good_tracks/(1.0 * total_tracks) << endl;

	inputFile.Close();

	return hist_arr;

}

int main() {
	
	gStyle->SetOptStat(1111111);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetLabelSize(0.045,"y");
	gStyle->SetLabelSize(0.045);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.3);
	gStyle->SetStatStyle(0);
	gStyle->SetTitleFontSize(0.060);

	TH1I ** h_arr;
	h_arr = test();
	TCanvas * c1 = new TCanvas("c1","",1280,800);
	gPad->SetLogy();
	h_arr[4]->Draw();
	c1->SaveAs("ncluster_on_trk_ntrk_m1_mod.pdf");
	return 0;
}
