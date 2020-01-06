/// run using ./run_ROC.sh signal_sample background_sample region
TGraph* Rocgraph(TH1* sig, TH1* bkg, bool sig_pop_upper);
TGraph* CRgraph(TH1* hist);
void setPars(TString parFile);
void setStyle();
TString bkgName, sigName;
TString bkgLabel, sigLabel;
TString varName;


void ROC(TString file="ROC_pars.txt"){

  setStyle();
  if (file.IsNull()){
    cout << "provide input parameter file" << endl;
    cin >> file;
  }
  setPars(file);

  TFile* bkgFile = TFile::Open(bkgName);
  TFile* sigFile = TFile::Open(sigName);
  if (!bkgFile || !sigFile) return -1;

  TH1D* bkgHist = (TH1D*) bkgFile->FindObjectAny(varName);
  TH1D* sigHist = (TH1D*) sigFile->FindObjectAny(varName);

  if (bkgHist->GetNbinsX() != sigHist->GetNbinsX()){
    cout << Form("background and signal have different %s binning", varName.Data()) << endl;
    return -1;
  }

  TString bkgSample = bkgName( bkgName.Last('/')+1, bkgName.Last('_')-bkgName.Last('/')-1);
  TString sigSample = sigName( sigName.Last('/')+1, sigName.Last('_')-sigName.Last('/')-1);
  varName = varName( varName.Index('_')+1, varName.Length());


  vector< pair<TString, TString> > labels = { {"ttbar", "t#bar{t}"}, {"dy", "Z/#gamma*#rightarrowl^{+}l^{-}"},
					      {"st", "Single-Top"}, {"vv", "VV"}, 
					      {"wjet", "W+Jets"}, {"gluon", "g_{kk}"}, {"zprime", "Z'"}};

  map<TString, TString> axis = {{"cosTheta1","cos #theta_{top,lep}"}, {"cosTheta2","cos #theta_{antitop,lep}"},
				{"cosTheta1r","cos #theta_{top,lep}"}, {"cosTheta2r","cos #theta_{antitop,lep}"},
				{"TOP_xl","#chi_{L}^{t}"}, {"ANTITOP_xl","#chi_{L}^{#bar{t}}"},
				{"T_xl","#chi_{L}^{t}"}, {"ANTIT_xl","#chi_{L}^{#bar{t}}"},
			        {"MT2s","MT2(GeV)"}, {"MT2r","MT2(GeV)"}, {"sT_met","S_{T}(GeV)"},
				{"rmin0","DeltaR_{min0}"}, {"rmin1","DeltaR_{min1}"}};

  map<TString, TString> var_labels = {{"cosTheta1","cos #theta_{top,lep}"}, {"cosTheta2","cos #theta_{antitop,lep}"},
				      {"cosTheta1r","cos #theta_{top,lep}"}, {"cosTheta2r","cos #theta_{antitop,lep}"},
				      {"TOP_xl","#chi_{L}^{t}"}, {"ANTITOP_xl","#chi_{L}^{#bar{t}}"},
				      {"T_xl","#chi_{L}^{t}"}, {"ANTIT_xl","#chi_{L}^{#bar{t}}"},
			              {"MT2s","MT2"}, {"MT2r","MT2"}, {"sT_met","S_{T}"},
				      {"rmin0","DeltaR_{min0}"}, {"rmin1","DeltaR_{min1}"}};

  // for signal accumulated at higher values, sig_density is true
  map<TString, bool> sig_density = {{"cosTheta1",true}, {"cosTheta2",true},
				    {"cosTheta1r",true}, {"cosTheta2r",true},
				    {"TOP_xl",false}, {"ANTITOP_xl",false},
				    {"T_xl",false}, {"ANTIT_xl",false},
			            {"MT2s",false}, {"MT2r",true}, {"sT_met",true},
				    {"rmin0",false}, {"rmin1",false}};

  for (auto const& i_sample : labels) {
    if(bkgSample.Contains(i_sample.first, TString::kIgnoreCase)) bkgLabel = i_sample.second;
    if(sigSample.Contains(i_sample.first, TString::kIgnoreCase)) sigLabel = i_sample.second;    
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->cd();
  TH1D* h_temp = new TH1D("h1", "h1", sigHist->GetNbinsX(), 0., 1.);
  h_temp->GetYaxis()->SetRangeUser(0,1);
  h_temp->GetXaxis()->SetTitle(Form("#epsilon_{ %s}",sigLabel.Data()));
  h_temp->GetYaxis()->SetTitle(Form("1-#epsilon_{ %s}",bkgLabel.Data()));
  h_temp->GetXaxis()->SetTitleOffset(0.65);
  h_temp->GetYaxis()->SetTitleOffset(0.65);
  h_temp->Draw();

  TGraph* gROC = Rocgraph(sigHist, bkgHist, sig_density[varName]); 
  double * X = gROC->GetX();
  TF1* f = new TF1("f",[&](double *X, double *){ return gROC->Eval(X[0]); },0,1,1);
  double integral = f->Integral(0,1);
  
  gROC->SetMarkerColor(kRed);
  gROC->SetMarkerStyle(22);  
  gROC->Draw("samep");
  TLine *l = new TLine(0.,1.,1,0.);
  l->Draw("same");

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.05);
  text.DrawLatex(0.14, 0.92, "CMS"); 
  text.SetTextSize(0.04);
  text.SetTextFont(42);
  text.DrawLatex(0.2, 0.2, Form("AUC = %4.3f",  integral));  
  text.DrawLatex(0.2, 0.25, Form("Discriminator: %s", var_labels[varName].Data())); 

  /*if(varName == "cosTheta1r"|| varName == "cosTheta2r"|| varName == "T_xl" ||varName == "ANTIT_xl" || 
    varName == "MT2r" || varName == "rmin0" || varName == "rmin1" || varName == "sT_met"){
    text.DrawLatex(0.6,0.84,"#bf{#geq 1 btag}");
    text.DrawLatex(0.6,0.74,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}");
    text.DrawLatex(0.6,0.64,"#bf{p_{T}^{miss} > 30 GeV}");
  }*/
  if (bkgName.Contains("ON/") && sigName.Contains("ON/")){
    text.DrawLatex(0.6,0.9,"#Delta R_{sum} < 2");
  }

  c->Print(Form("roc_%s.pdf",varName.Data()));
  c->Clear();
  c->cd();
  vector<double> var_bins;
   for (int i=0; i < sigHist->GetNbinsX(); i++) {
      var_bins.push_back(sigHist->GetBinLowEdge(i+1));
   }

  TH1D* h_temp2 = new TH1D("h2", "h2", var_bins.size()-1,  &var_bins[0]);
  h_temp2->GetYaxis()->SetRangeUser(0,1);
  h_temp2->GetXaxis()->SetTitle(axis[varName].Data());
  h_temp2->GetYaxis()->SetTitle(Form("#Rgothic_{%s}",axis[varName].Data()));
  h_temp2->GetXaxis()->SetTitleOffset(0.9);
  h_temp2->GetYaxis()->SetTitleOffset(0.9);
  h_temp2->Draw();

  TGraph* gCR_sig = CRgraph(sigHist);
  TGraph* gCR_bkg = CRgraph(bkgHist);
  gCR_sig->SetMarkerColor(kRed);
  gCR_sig->SetMarkerStyle(20);
  gCR_sig->SetMarkerSize(1);
  gCR_sig->Draw("samep");

  gCR_bkg->SetMarkerColor(kBlue);
  gCR_bkg->SetMarkerStyle(20);
  gCR_bkg->SetMarkerSize(1);
  gCR_bkg->Draw("samep");

  text.SetTextColor(4);
  text.DrawLatex(0.6, 0.3, Form("Background: %s", bkgLabel.Data())); 
  text.SetTextColor(2);
  text.DrawLatex(0.6, 0.25, Form("Signal: %s",     sigLabel.Data())); 
  text.SetTextColor(1);
  text.SetTextSize(0.05);
  text.DrawLatex(0.14, 0.92, "CMS"); 
  text.SetTextSize(0.04);
  text.SetTextFont(42);
  c->Print(Form("cm_%s.pdf",varName.Data()));
  delete h_temp;
  delete h_temp2;
  delete c;

}

TGraph* Rocgraph(TH1* sig, TH1* bkg, bool sig_pop_upper) {

    const int nbins = sig->GetNbinsX();
    double SIG[nbins];
    double BKG[nbins];
    for ( int i = 0; i != nbins; ++i ) {
      if(sig_pop_upper){
        SIG[i] = sig->Integral(i+1,nbins)/sig->Integral();
        BKG[i] = 1.0 - ( bkg->Integral(i+1,nbins)/bkg->Integral());
      }
      else {
        SIG[i] = sig->Integral(1,1+i)/sig->Integral();
        BKG[i] = 1.0 - ( bkg->Integral(1,1+i)/bkg->Integral());
      }
    }
    TGraph *gROC = new TGraph(nbins,SIG,BKG);    
    return gROC;
}

TGraph* CRgraph(TH1* hist) {

    const int nbins = hist->GetNbinsX();
    double CR[nbins];
    double var[nbins];
    for ( int i = 0; i != nbins; ++i ) {
      var[i] = hist->GetBinCenter(i+1);
      CR[i] = hist->Integral(1,i+1)/hist->Integral();
    }
    TGraph *gCR = new TGraph(nbins,var,CR);
    return gCR;
}

void setPars(TString parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)){

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if (var == "bkgName")         bkgName = line;
    else if (var == "sigName")    sigName = line;
    else if (var == "varName")    varName = line;

  }
  file.close();
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
  tdrStyle->SetGridWidth(0);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.12);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.15);
  tdrStyle->SetPadRightMargin(0.10);

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
  tdrStyle->SetTitleSize(0.05, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.0);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.026, "XYZ");

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
}
