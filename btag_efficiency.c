// Chad Harrington and Bahareh Roozbahani 5/11/2017 - btagging efficiency

void setStyle();

void btag_efficiency(TString dir = "OFF/", double ptDR = 300) {

  map<TString, TString> channels = { {"mm","#mu#mu"}, {"ee","ee"}, {"em","e#mu"}, {"z","ll"} };  // channel, label
  map<TString, float>   flavors  = { {"b",1.4}, {"c",1.4}, {"udsg",0.5} };  // flavor, ymax
  map<TString, pair<TString, vector<double> > > vars = {
    {"Pt", {"Jet p_{T} (GeV)", {0, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1400} } },
    {"Eta",{"Jet #eta", {-2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4} } },
    {"DR", {"#DeltaR_{min}(lepton, jet)", {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5,
                                           2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4} } }
  };  // variable, {axis title, bins}

  vector <TString> files = { "wjets", "ww", "wz", "zz", "DYJetsToLL_M10to50", "DYJetsToLL_M50",
                             "TTbar0-700", "TTbar700-1000_total", "TTbar1000-inf_total", "TTsemilep",
                             "ST_t-channel_top_4f", "ST_t-channel_antitop_4f", "ST_tW_top_5f_inclusiveDecays_v2", "ST_tW_antitop_5f_v2" };

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetGrid();

  TLegend* leg = new TLegend(.6,.9-2*0.04,.8,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  TLatex text;
  text.SetNDC();

  map<TString, TGraphAsymmErrors*> effs;
  for (auto const& it_var : vars) {
    TString var = it_var.first, title = it_var.second.first;
    vector<double> bins = it_var.second.second;

    TH1D* h = new TH1D("h", "h", bins.size()-1, &bins[0]);
    h->GetXaxis()->SetTitle(title);
    h->GetYaxis()->SetTitle("b-tagging #varepsilon");
    h->GetYaxis()->SetTitleOffset(1.2);

    for (auto const& it_flav : flavors) {
      TString flav = it_flav.first;  float ymax = it_flav.second;

      TH1D* comb_denom = new TH1D("comb_denom", "comb_denom", bins.size()-1, &bins[0]);  comb_denom->Sumw2();
      TH1D* comb_LWP =   new TH1D("comb_LWP", "comb_LWP", bins.size()-1, &bins[0]);      comb_LWP->Sumw2();
      TH1D* comb_MWP =   new TH1D("comb_MWP", "comb_MWP", bins.size()-1, &bins[0]);      comb_MWP->Sumw2();

      for (auto const& it_chan : channels) {
        TString chan = it_chan.first, clabel = it_chan.second;
        chan = chan=="z" ? "ll" : chan;

        TString LWP = var + "_" + flav + "_LWP_" + chan;
        TString MWP = var + "_" + flav + "_MWP_" + chan;
        effs[LWP] = new TGraphAsymmErrors();
        effs[MWP] = new TGraphAsymmErrors();

        if (chan != "ll") {
          TH1D* chan_denom = new TH1D("chan_denom", "chan_denom", bins.size()-1, &bins[0]);  chan_denom->Sumw2();
          TH1D* chan_LWP =   new TH1D("chan_LWP", "chan_LWP", bins.size()-1, &bins[0]);      chan_LWP->Sumw2();
          TH1D* chan_MWP =   new TH1D("chan_MWP", "chan_MWP", bins.size()-1, &bins[0]);      chan_MWP->Sumw2();

          for (auto const& fname : files) {
            TFile* file = TFile::Open( dir + chan + "/" + fname + "_" + chan + ".root" );
            if (!file) return;

            TH1D* denom, *numL, *numM;
            if (var == "DR") {
              TH2D* denom2 = (TH2D*) file->Get( "jetPtDR_" + flav );
              denom = denom2->ProjectionX( "denomDR", 1, denom2->GetYaxis()->FindBin(ptDR)-1 );
              TH2D* numL2 = (TH2D*) file->Get( "jetPtDR_bTagL_" + flav );
              numL = numL2->ProjectionX( "numLDR", 1, numL2->GetYaxis()->FindBin(ptDR)-1 );
              TH2D* numM2 = (TH2D*) file->Get( "jetPtDR_bTagM_" + flav );
              numM = numM2->ProjectionX( "numMDR", 1, numM2->GetYaxis()->FindBin(ptDR)-1 );
            }
            else {
              denom = (TH1D*) file->Get( "jet" + var + "_" + flav );
              numL = (TH1D*) file->Get( "jet" + var + "_bTagL_" + flav );
              numM = (TH1D*) file->Get( "jet" + var + "_bTagM_" + flav );
            }
            denom->Sumw2();  denom = (TH1D*) denom->Rebin(bins.size()-1, "denom", &bins[0]);
            chan_denom->Add(denom);  comb_denom->Add(denom);

            numL->Sumw2();  numL = (TH1D*) numL->Rebin(bins.size()-1, "numL", &bins[0]);
            chan_LWP->Add(numL);  comb_LWP->Add(numL);

            numM->Sumw2();  numM = (TH1D*) numM->Rebin(bins.size()-1, "numM", &bins[0]);
            chan_MWP->Add(numM);  comb_MWP->Add(numM);
          }
          effs[LWP]->BayesDivide(chan_LWP, chan_denom);
          effs[MWP]->BayesDivide(chan_MWP, chan_denom);
        }
        else {
          effs[LWP]->BayesDivide(comb_LWP, comb_denom);
          effs[MWP]->BayesDivide(comb_MWP, comb_denom);
        }
        h->GetYaxis()->SetRangeUser(0, ymax);
        h->Draw();

        effs[LWP]->SetMarkerStyle(20);
        effs[LWP]->SetMarkerColor(1);
        effs[LWP]->SetLineColor(1);
        effs[LWP]->Draw("pzsame");

        effs[MWP]->SetMarkerStyle(24);
        effs[MWP]->SetMarkerColor(1);
        effs[MWP]->SetLineColor(1);
        effs[MWP]->Draw("pzsame");

        leg->AddEntry(effs[LWP], "DeepCSVL", "PL");
        leg->AddEntry(effs[MWP], "DeepCSVM", "PL");
        leg->Draw();

        text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
        text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
        text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

        text.SetTextSize(0.04);  text.SetTextFont(62);
        text.DrawLatex(0.2, 0.87, flav + "-Jets");
        text.DrawLatex(0.2, 0.83, clabel + " channel");

        c->Print("plots/btag_eff_" + (var=="DR" ? "DRpt"+to_string(int(ptDR)) : var) + "_" + flav + "_" + chan + ".pdf");
        leg->Clear();
        c->Clear();
      }
    }
  }

  // Btag root file //

  TFile* outFile = new TFile("btag_eff.root","RECREATE");

  for (auto const& it_eff : effs) {
    TString name = it_eff.first, chan = name(name.Last('_')+1, 2);

    if (!name.Contains("Pt_") || chan == "ll") continue;
    name.Remove(0, name.Index('_')+1);

    outFile->cd();
    outFile->mkdir(chan + "/");
    outFile->cd("btag_eff.root:/" + chan);
    it_eff.second->Write(name);
  }

  outFile->Write();
  delete outFile;
  outFile = 0;

  return;
  // Create individual plots //

  vector <TString> filesToPlot = { "TTbar0-700", "TTbar700-1000_total", "TTbar1000-inf_total", "TTsemilep",
                                   "zprime_M-3000_W-300", "zprime_M-3500_W-350", "zprime_M-4000_W-400", "zprime_M-4500_W-450" };

  //vector <TString> filesToPlot = { "STschannel", "STtWchannel", "STtchannel", "SaTtWchannel", "SaTtchannel",
  //                                 "zprime_M-3000_W-300", "zprime_M-3500_W-350", "zprime_M-4000_W-400", "zprime_M-4500_W-450" };

  leg->SetY1NDC(.9-4*0.04);

  for (auto const& it_var : vars) {
    TString var = it_var.first, title = it_var.second.first;
    vector<double> bins = it_var.second.second;

    TH1D* h = new TH1D("h2", "h2", bins.size()-1, &bins[0]);
    h->GetXaxis()->SetTitle(title);
    h->GetYaxis()->SetTitle("b-tagging #varepsilon");
    h->GetYaxis()->SetTitleOffset(1.2);

    for (auto const& it_flav : flavors) {
      TString flav = it_flav.first;  float ymax = it_flav.second;

      TH1D* comb_denom =  new TH1D("comb_denom_ind", "comb_denom_ind", bins.size()-1, &bins[0]);  comb_denom->Sumw2();
      TH1D* comb_LWP =    new TH1D("comb_LWP_ind", "comb_LWP_ind", bins.size()-1, &bins[0]);      comb_LWP->Sumw2();
      TH1D* comb_MWP =    new TH1D("comb_MWP_ind", "comb_MWP_ind", bins.size()-1, &bins[0]);      comb_MWP->Sumw2();

      TH1D* comb_sigdenom = new TH1D("comb_sigdenom", "comb_sigdenom", bins.size()-1, &bins[0]);  comb_sigdenom->Sumw2();
      TH1D* comb_sigLWP =   new TH1D("comb_sigLWP", "comb_sigLWP", bins.size()-1, &bins[0]);      comb_sigLWP->Sumw2();
      TH1D* comb_sigMWP =   new TH1D("comb_sigMWP", "comb_sigMWP", bins.size()-1, &bins[0]);      comb_sigMWP->Sumw2();

      for (auto const& it_chan : channels) {
        TString chan = it_chan.first, clabel = it_chan.second;
        chan = chan=="z" ? "ll" : chan;

        TGraphAsymmErrors* g_LWP = new TGraphAsymmErrors();
        TGraphAsymmErrors* g_MWP = new TGraphAsymmErrors();
        TGraphAsymmErrors* g_sigLWP = new TGraphAsymmErrors();
        TGraphAsymmErrors* g_sigMWP = new TGraphAsymmErrors();

        if (chan != "ll") {
          TH1D* chan_denom =  new TH1D("chan_denom_ind", "chan_denom_ind", bins.size()-1, &bins[0]);  chan_denom->Sumw2();
          TH1D* chan_LWP =    new TH1D("chan_LWP_ind", "chan_LWP_ind", bins.size()-1, &bins[0]);      chan_LWP->Sumw2();
          TH1D* chan_MWP =    new TH1D("chan_MWP_ind", "chan_MWP_ind", bins.size()-1, &bins[0]);      chan_MWP->Sumw2();

          TH1D* chan_sigdenom = new TH1D("chan_sigdenom", "chan_sigdenom", bins.size()-1, &bins[0]);  chan_sigdenom->Sumw2();
          TH1D* chan_sigLWP =   new TH1D("chan_sigLWP", "chan_sigLWP", bins.size()-1, &bins[0]);      chan_sigLWP->Sumw2();
          TH1D* chan_sigMWP =   new TH1D("chan_sigMWP", "chan_sigMWP", bins.size()-1, &bins[0]);      chan_sigMWP->Sumw2();

          for (auto const& fname : filesToPlot) {
            TFile* file = TFile::Open( dir + chan + "/" + fname + "_" + chan + ".root" );

            TH1D* denom, *numL, *numM;
            if (var == "DR") {
              TH2D* denom2 = (TH2D*) file->Get( "jetPtDR_" + flav );
              denom = denom2->ProjectionX( "denomDR", 1, denom2->GetYaxis()->FindBin(ptDR)-1 );
              TH2D* numL2 = (TH2D*) file->Get( "jetPtDR_bTagL_" + flav );
              numL = numL2->ProjectionX( "numLDR", 1, numL2->GetYaxis()->FindBin(ptDR)-1 );
              TH2D* numM2 = (TH2D*) file->Get( "jetPtDR_bTagM_" + flav );
              numM = numM2->ProjectionX( "numMDR", 1, numM2->GetYaxis()->FindBin(ptDR)-1 );
            }
            else {
              denom = (TH1D*) file->Get( "jet" + var + "_" + flav );
              numL = (TH1D*) file->Get( "jet" + var + "_bTagL_" + flav );
              numM = (TH1D*) file->Get( "jet" + var + "_bTagM_" + flav );
            }
            denom->Sumw2();  denom = (TH1D*) denom->Rebin(bins.size()-1, "denom", &bins[0]);

            numL->Sumw2();  numL = (TH1D*) numL->Rebin(bins.size()-1, "numL", &bins[0]);

            numM->Sumw2();  numM = (TH1D*) numM->Rebin(bins.size()-1, "numM", &bins[0]);

            if ( fname.Contains("zprime", TString::kIgnoreCase) ) {
              chan_sigdenom->Add(denom);  comb_sigdenom->Add(denom);
              chan_sigLWP->Add(numL);     comb_sigLWP->Add(numL);
              chan_sigMWP->Add(numM);     comb_sigMWP->Add(numM);
            }
            else {
              chan_denom->Add(denom);  comb_denom->Add(denom);
              chan_LWP->Add(numL);     comb_LWP->Add(numL);
              chan_MWP->Add(numM);     comb_MWP->Add(numM);
            }
          }
          g_LWP->BayesDivide(chan_LWP, chan_denom);
          g_MWP->BayesDivide(chan_MWP, chan_denom);
          g_sigLWP->BayesDivide(chan_sigLWP, chan_sigdenom);
          g_sigMWP->BayesDivide(chan_sigMWP, chan_sigdenom);
        }
        else {
          g_LWP->BayesDivide(comb_LWP, comb_denom);
          g_MWP->BayesDivide(comb_MWP, comb_denom);
          g_sigLWP->BayesDivide(comb_sigLWP, comb_sigdenom);
          g_sigMWP->BayesDivide(comb_sigMWP, comb_sigdenom);
        }
        h->GetYaxis()->SetRangeUser(0, ymax);
        h->Draw();

        g_LWP->SetMarkerStyle(20);
        g_LWP->SetMarkerColor(kBlack);
        g_LWP->SetLineColor(kBlack);
        g_LWP->Draw("pzsame");

        g_MWP->SetMarkerStyle(24);
        g_MWP->SetMarkerColor(kBlack);
        g_MWP->SetLineColor(kBlack);
        g_MWP->Draw("pzsame");

        g_sigLWP->SetMarkerStyle(21);
        g_sigLWP->SetMarkerColor(kRed);
        g_sigLWP->SetLineColor(kRed);
        g_sigLWP->Draw("pzsame");

        g_sigMWP->SetMarkerStyle(25);
        g_sigMWP->SetMarkerColor(kRed);
        g_sigMWP->SetLineColor(kRed);
        g_sigMWP->Draw("pzsame");

        leg->AddEntry(g_LWP, "DeepCSVL,  ttbar", "PL");
        leg->AddEntry(g_MWP, "DeepCSVM, ttbar", "PL");
        //leg->AddEntry(g_LWP, "DeepCSVL,  Single-Top", "PL");
        //leg->AddEntry(g_MWP, "DeepCSVM, Single-Top", "PL");
        leg->AddEntry(g_sigLWP, "DeepCSVL,  Z'", "PL");
        leg->AddEntry(g_sigMWP, "DeepCSVM, Z'", "PL");
        leg->Draw();

        text.SetTextSize(0.05);  text.SetTextFont(61); text.DrawLatex(0.18, 0.96, "CMS");
        text.SetTextSize(0.04);  text.SetTextFont(52); text.DrawLatex(0.29, 0.96, "Simulation");
        text.SetTextSize(0.035); text.SetTextFont(42); text.DrawLatex(0.85, 0.96, "(13 TeV)");

        text.SetTextSize(0.04);  text.SetTextFont(62);
        text.DrawLatex(0.2, 0.87, flav + "-Jets");
        text.DrawLatex(0.2, 0.83, clabel + " channel");

        c->Print("plots/btag_eff_" + (var=="DR" ? "DRpt"+to_string(int(ptDR)) : var) + "_" + flav + "_" + chan + "_individual.pdf");
        leg->Clear();
        c->Clear();
      }
    }
  }

}

void setStyle(){

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();

//End Style//
}
