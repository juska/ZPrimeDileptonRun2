//Chad Harrington - 10/19/2017

void ttbar( TString fname="root://131.225.204.161:1094//store/user/cihar29/2017Analysis_trees/TTjets.root",
            TString outdir="/uscmst1b_scratch/lpc1/3DayLifetime/cihar29/", int xs = 831760 ) {

  TFile* f_incl = TFile::Open(fname);

  TString htotal_name = "mass_totalEvts", hmet_name = "mass_met", fullmu_name = "fullMu";
  TH1D* htotal = (TH1D*) f_incl->Get( htotal_name );
  TH1D* hmet = (TH1D*) f_incl->Get( hmet_name );
  double totalBins = htotal->GetNbinsX(), totalLow = htotal->GetBinLowEdge(1), totalHigh = htotal->GetBinLowEdge(totalBins+1);
  double metBins = hmet->GetNbinsX(), metLow = hmet->GetBinLowEdge(1), metHigh = hmet->GetBinLowEdge(metBins+1);
  delete htotal;
  delete hmet;

  TTree* tree = (TTree*) f_incl->Get("T");
  tree->SetBranchStatus("*",1);
  TH2F* mu_mass = (TH2F*) f_incl->Get(fullmu_name+"_mass");
  double mu700Bin = mu_mass->GetXaxis()->FindBin(700), mu1000Bin = mu_mass->GetXaxis()->FindBin(1000), muBins = mu_mass->GetNbinsX();

  TFile *f1 = new TFile(outdir+"TTbar0-700.root","recreate");
  TTree *t1 = tree->CloneTree(0);
  TH1D  *htotal1 = new TH1D(htotal_name,htotal_name,totalBins,totalLow,totalHigh);
  TH1D  *hmet1 = new TH1D(hmet_name,hmet_name,metBins,metLow,metHigh);
  TH1D  *mu1 = mu_mass->ProjectionY(fullmu_name, 1, mu700Bin-1);  //Want mu to be TH1F* and match other mu distributions

  TFile *f2 = new TFile(outdir+"TTbar700-1000_fromTTincl.root","recreate");
  TTree *t2 = tree->CloneTree(0);
  TH1D  *htotal2 = new TH1D(htotal_name,htotal_name,totalBins,totalLow,totalHigh);
  TH1D  *hmet2 = new TH1D(hmet_name,hmet_name,metBins,metLow,metHigh);
  TH1D  *mu2 = mu_mass->ProjectionY(fullmu_name, mu700Bin, mu1000Bin-1);

  TFile *f3 = new TFile(outdir+"TTbar1000-inf_fromTTincl.root","recreate");
  TTree *t3 = tree->CloneTree(0);
  TH1D  *htotal3 = new TH1D(htotal_name,htotal_name,totalBins,totalLow,totalHigh);
  TH1D  *hmet3 = new TH1D(hmet_name,hmet_name,metBins,metLow,metHigh);
  TH1D  *mu3 = mu_mass->ProjectionY(fullmu_name, mu1000Bin, muBins);

  delete mu_mass;

  TString intnames[] = { "totalEvts", "dilep_cut", "leppt_cut", "dilepmass_cut", "jetpteta_cut", "met_cut" };
  TString doublenames[] = { "topPtWeightNOM", "topPtWeightDN", "pdfUP", "pdfDN", "q2UP", "q2DN" };

  vector< vector<int> > ints1 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<int> > ints2 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<int> > ints3 = { {0}, {0}, {0}, {0}, {0}, {0} };
  vector< vector<double> > doubles1 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };
  vector< vector<double> > doubles2 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };
  vector< vector<double> > doubles3 = { {0.}, {0.}, {0.}, {0.}, {0.}, {0.} };

  for (int i=0; i<6; i++) {

    TString hiname = "mass_" + intnames[i];
    hiname.ReplaceAll("_cut", "");
    TString hdname = "mass_" + doublenames[i];

    TH1D* hi = (TH1D*) f_incl->Get(hiname);
    TH1D* hd = (TH1D*) f_incl->Get(hdname);
    double hi700Bin = hi->FindBin(700), hi1000Bin = hi->FindBin(1000), hiBins = hi->GetNbinsX();
    double hd700Bin = hd->FindBin(700), hd1000Bin = hd->FindBin(1000), hdBins = hd->GetNbinsX();

    ints1[i][0] = int( hi->Integral(1,hi700Bin-1) );
    ints2[i][0] = int( hi->Integral(hi700Bin,hi1000Bin-1) );
    ints3[i][0] = int( hi->Integral(hi1000Bin,hiBins) );
    double itotal = hi->Integral(1,hiBins);

    doubles1[i][0] = hd->Integral(1,hd700Bin-1);
    doubles2[i][0] = hd->Integral(hd700Bin,hd1000Bin-1);
    doubles3[i][0] = hd->Integral(hd1000Bin,hdBins);
    double dtotal = hd->Integral(1,hdBins);

    f1->WriteObject(&ints1[i], intnames[i]);
    f2->WriteObject(&ints2[i], intnames[i]);
    f3->WriteObject(&ints3[i], intnames[i]);

    f1->WriteObject(&doubles1[i], doublenames[i]);
    f2->WriteObject(&doubles2[i], doublenames[i]);
    f3->WriteObject(&doubles3[i], doublenames[i]);

    if (intnames[i] == "totalEvts") { //no toppt weighting is UP
      cout << "TTbar0-700_topPtWeightUP\t" << ints1[i][0]/itotal*xs << endl;
      cout << "TTbar700-1000_topPtWeightUP\t" << ints2[i][0]/itotal*xs << endl;
      cout << "TTbar1000-inf_topPtWeightUP\t" << ints3[i][0]/itotal*xs << endl;

      for (int i=1; i<=hiBins; i++) {
        double content = hi->GetBinContent(i);

        if      (i<hi700Bin)  htotal1->SetBinContent(i, content);
        else if (i<hi1000Bin) htotal2->SetBinContent(i, content);
        else                  htotal3->SetBinContent(i, content);
      }
      htotal1->SetEntries(ints1[i][0]);
      htotal2->SetEntries(ints2[i][0]);
      htotal3->SetEntries(ints3[i][0]);
    }
    cout << "TTbar0-700_" + doublenames[i] << "\t" << doubles1[i][0]/dtotal*xs << endl;
    cout << "TTbar700-1000_" + doublenames[i] << "\t" << doubles2[i][0]/dtotal*xs << endl;
    cout << "TTbar1000-inf_" + doublenames[i] << "\t" << doubles3[i][0]/dtotal*xs << endl;
  }

  float mass_ttbar;
  tree->SetBranchAddress("mass_ttbar", &mass_ttbar);

  Long64_t nEntries = tree->GetEntries();
  for (Long64_t n=0; n<nEntries; n++) {
    tree->GetEntry(n);

    if      (mass_ttbar < 700)  { t1->Fill(); hmet1->Fill(mass_ttbar); }
    else if (mass_ttbar < 1000) { t2->Fill(); hmet2->Fill(mass_ttbar); }
    else                        { t3->Fill(); hmet3->Fill(mass_ttbar); }
  }

  f1->Write();
  f2->Write();
  f3->Write();

  f_incl->Close();
  f1->Close();
  f2->Close();
  f3->Close();
  delete f_incl;
  delete f1;
  delete f2;
  delete f3;
  cout << "done" << endl;
}
