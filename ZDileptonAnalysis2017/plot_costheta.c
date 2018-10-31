#include "TH1.h"
#include "TCanvas.h"
void setStyle();
void plot_costheta(){

  map<TString,TFile*> inFile;
  inFile["ttbar"] = TFile::Open("TTBar.root");
  inFile["gkk"] = TFile::Open("gkk_3000.root");
  inFile["zprime"] = TFile::Open("zprime_3tev_10.root");

  map<TString,TTree*> tree;

  int nGen=1;
  float lep_top_costheta[nGen];
  float lep_antitop_costheta[nGen];

  map<TString, TH1F*> h_top_costheta, h_antitop_costheta;
  map<TString, TH2F*> topvsantitop;


  for (map<TString,TFile*> ::iterator it = inFile.begin(); it != inFile.end(); it++) {
    tree[it->first]= (TTree*) inFile[it->first]->Get("T");
  }
  for (map<TString,TTree*> ::iterator it = tree.begin(); it != tree.end(); it++) {
    tree[it->first]->SetBranchStatus("*",1);
    tree[it->first]->SetBranchAddress("lep_top_costheta", lep_top_costheta);
    tree[it->first]->SetBranchAddress("lep_antitop_costheta", lep_antitop_costheta);
    Long64_t nEntries = tree[it->first]->GetEntries();
    cout << nEntries << " Events" << endl;
    h_top_costheta[it->first] =  new TH1F(Form("h_top_costheta_%s",it->first.Data()), Form("h_top_costheta_%s",it->first.Data()), 40 , -1,1);
    h_antitop_costheta[it->first] =  new TH1F(Form("h_antitop_costheta_%s",it->first.Data()), Form("h_antitop_costheta_%s",it->first.Data()), 40 , -1,1);
    topvsantitop[it->first] = new TH2F (Form("%s",it->first.Data()),Form("%s",it->first.Data()),40,-1,1,40,-1.0,1.0);
    for (Long64_t n=0; n<nEntries; n++) {
      tree[it->first]->GetEntry(n);
      for (int i=0; i<1; i++) {
        //cout<<lep_top_costheta[i]<<endl;
        h_top_costheta[it->first]->Fill(lep_top_costheta[i]);
        h_antitop_costheta[it->first]->Fill(lep_antitop_costheta[i]);
        topvsantitop[it->first] ->Fill(lep_top_costheta[i],lep_antitop_costheta[i],1);
      }  
    }
  }
  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);

  TLegend* leg = new TLegend(.4,.65,.8,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  TH1F* h_cos = new TH1F ("","",40,-1.,1.);
  h_cos->GetYaxis()->SetRangeUser(0,0.06);
  c->cd();
  h_cos->GetXaxis()->SetTitle("cos #theta");
  h_cos->GetXaxis()->SetTitle("cos #theta_{top and lep}");
  h_cos->Draw();
  h_top_costheta["ttbar"] -> SetLineColor(kBlue);
  h_top_costheta["gkk"] -> SetLineColor(kRed);
  h_top_costheta["zprime"] -> SetLineColor(kGreen);
  h_top_costheta["ttbar"] -> Scale(1/h_top_costheta["ttbar"]->Integral());
  h_top_costheta["gkk"] -> Scale(1/h_top_costheta["gkk"]->Integral());
  h_top_costheta["zprime"] -> Scale(1/h_top_costheta["zprime"]->Integral());
  for (map<TString,TH1F*> ::iterator it = h_top_costheta.begin(); it != h_top_costheta.end(); it++) {  
    leg->SetTextSize(0.036);
    //leg->AddEntry(h_top_costheta[it->first], it->first.Data() , "l"); 
    h_top_costheta[it->first] -> Draw("sameHIST");
  }
  //leg->Draw();
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.02);

  text.SetTextColor(4);
  text.DrawLatex(0.2, 0.9, Form("TTBar mean = %4.3f , stdDev = %4.3f", h_top_costheta["ttbar"]->GetMean(), h_top_costheta["ttbar"]->GetStdDev() ) );
  text.SetTextColor(3);
  text.SetTextColor(2);
  text.DrawLatex(0.2, 0.85, Form("g_{kk} mean = %4.3f , stdDev = %4.3f", h_top_costheta["gkk"]->GetMean(), h_top_costheta["gkk"]->GetStdDev() ) );
  text.SetTextColor(3);
  text.SetTextColor(3);
  text.DrawLatex(0.2, 0.8, Form("Z' mean = %4.3f , stdDev = %4.3f", h_top_costheta["zprime"]->GetMean(), h_top_costheta["zprime"]->GetStdDev() ) );
  c->Print("costop.pdf");
  leg->Clear();
  c->Clear();

  c->cd();
  h_cos->Draw();
  h_antitop_costheta["ttbar"] -> SetLineColor(kBlue);
  h_antitop_costheta["gkk"] -> SetLineColor(kRed);
  h_antitop_costheta["zprime"] -> SetLineColor(kGreen);
  h_antitop_costheta["ttbar"] -> Scale(1/h_antitop_costheta["ttbar"]->Integral());
  h_antitop_costheta["gkk"] -> Scale(1/h_antitop_costheta["gkk"]->Integral());
  h_antitop_costheta["zprime"] -> Scale(1/h_antitop_costheta["zprime"]->Integral());
  for (map<TString,TH1F*> ::iterator it = h_antitop_costheta.begin(); it != h_antitop_costheta.end(); it++) {  
    leg->SetTextSize(0.036);
    leg->AddEntry(h_antitop_costheta[it->first], it->first.Data() , "l"); 
    h_antitop_costheta[it->first] -> Draw("sameHIST");
  }
  //leg->Draw();
  text.SetTextColor(4);
  text.DrawLatex(0.2, 0.9, Form("TTBar mean = %4.3f , stdDev = %4.3f", h_antitop_costheta["ttbar"]->GetMean(), h_antitop_costheta["ttbar"]->GetStdDev() ) );
  text.SetTextColor(2);
  text.DrawLatex(0.2, 0.85, Form("g_{kk} mean = %4.3f , stdDev = %4.3f", h_antitop_costheta["gkk"]->GetMean(), h_antitop_costheta["gkk"]->GetStdDev() ) );
  text.SetTextColor(3);
  text.DrawLatex(0.2, 0.8, Form("Z' mean = %4.3f , stdDev = %4.3f", h_antitop_costheta["zprime"]->GetMean(), h_antitop_costheta["zprime"]->GetStdDev() ) );
  text.SetTextColor(3);
  c->Print("cosantitop.pdf");
  leg->Clear();
  c->Clear();

  c->cd();
  h_cos->GetXaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetYaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetZaxis()->SetRangeUser(0,10);
  h_cos->Draw();
  topvsantitop["ttbar"] -> Draw("colz");
  topvsantitop["ttbar"]->GetXaxis()->SetTitle("cos #theta_{top,lep}");
  topvsantitop["ttbar"]->GetYaxis()->SetTitle("cos #theta_{antitop,lep}");
  text.SetTextSize(0.04);

  text.SetTextColor(1);
  text.DrawLatex(0.2, 0.97, "ttbar" );
  c->Print("ttbar_2D.pdf");
  leg->Clear();
  c->Clear();

  c->cd();
  h_cos->GetXaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetYaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetZaxis()->SetRangeUser(0,10);
  h_cos->Draw();
  topvsantitop["gkk"] -> Draw("colz");
  topvsantitop["gkk"]->GetXaxis()->SetTitle("cos #theta_{top,lep}");
  topvsantitop["gkk"]->GetYaxis()->SetTitle("cos #theta_{antitop,lep}");
  text.SetTextSize(0.04);

  text.SetTextColor(1);
  text.DrawLatex(0.2, 0.97, "g_{kk}" );
  c->Print("gkk_2D.pdf");
  leg->Clear();
  c->Clear();

  c->cd();
  h_cos->GetXaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetYaxis()->SetRangeUser(-1.1,1.1);
  h_cos->GetZaxis()->SetRangeUser(0,10);
  h_cos->Draw();
  topvsantitop["zprime"] -> Draw("colz");
  topvsantitop["zprime"]->GetXaxis()->SetTitle("cos #theta_{top,lep}");
  topvsantitop["zprime"]->GetYaxis()->SetTitle("cos #theta_{antitop,lep}");
  text.SetTextSize(0.04);

  text.SetTextColor(1);
  text.DrawLatex(0.2, 0.97, "zprime" );
  c->Print("zprime_2D.pdf");
  leg->Clear();
  c->Clear();
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

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.12);

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
}
