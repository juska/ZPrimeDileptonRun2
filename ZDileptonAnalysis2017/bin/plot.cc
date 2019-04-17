//charles harrington and bahareh roozbahani
//execute as plot plot_pars.txt

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TKey.h"
#include "TLine.h"

#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double delta(double delta1, double delta2, string pos_neg);
void setPars(const string& parFile);
void setStyle();
void drawText();
bool isBkg( const TString& set ) { return set=="ttbar" || set=="dy" || set=="wjet" || set=="st" || set=="vv"; }

//parameters- edit in plot_pars.txt
vector<TString> fileNames, systematics;
vector<double> rebin;
map<TString, float> sys_norm;
TString dataFileName, dir, outName, channel, theta, postfilename, region, hname, rightText, subplot, xtitle, zprime, gluon;
float xmin, xmax, ymin, ymax, subymin, subymax, sigScale;
bool logx=false, logy=false, plotData=false, plotImpact=false, fit=false;
static Int_t MLightGray = TColor::GetColor( "#dddddd" );
static Int_t LightGray  = TColor::GetColor( "#aaaaaa" );

int main(int argc, char* argv[]) {

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  setStyle();

  vector< pair<TString, TString> > labels = { {"ttbar", "t#bar{t}"}, {"dy", "Z/#gamma*#rightarrowl^{+}l^{-}"}, {"st", "Single-Top"}, {"vv", "VV"}, {"wjet", "W+Jets"},
                                              {"gluon", "g_{kk} 3 TeV (#sigma=10 pb)"}, {"zprime", "Z'   3 TeV (#sigma=10 pb)"}, {"bkg", "Background"} };

  map<TString, TH1D*> m_NOM;
  map<TString, map<TString, TH1D*> > m_UP, m_DN;

  fileNames.push_back( dataFileName );
  for (auto const& fname : fileNames) {

    TString set = "";

    if      ( fname == dataFileName )                                                set = "DATA";
    else if ( fname.Contains("gluon", TString::kIgnoreCase) ||
              fname.Contains("zprime", TString::kIgnoreCase) ) {

      if (theta == "gkk") {
        if ( !fname.Contains("gluon", TString::kIgnoreCase) ) continue;
        TString mass = fname(fname.Index("M-")+2, fname.Length()-fname.Index("M-")-2);
        set = "gkk" + mass;
      }
      else if (theta == "zp1" || theta == "zp10" || theta == "zp30") {
        if ( !fname.Contains("zprime", TString::kIgnoreCase) ) continue;

        TString mass = fname(fname.Index("M-")+2, fname.Index("_W")-fname.Index("M-")-2);
        TString width = fname(fname.Index("W-")+2, fname.Length()-fname.Index("W-")-2);
        set = "zp" + mass + "_" + width;

        if ( width.Contains('p', TString::kIgnoreCase) ) width.ReplaceAll('p', '.');
        int percent = stof( width.Data() ) / stof( mass.Data() ) * 100;

        if      (theta == "zp1"  && percent != 1)  continue;
        else if (theta == "zp10" && percent != 10) continue;
        else if (theta == "zp30" && percent != 30) continue;
      }
      else {
        if      ( fname == zprime )                                                  set = "zprime";
        else if ( fname == gluon )                                                   set = "gluon";
        else continue;
      }
    }

    else if ( fname.Contains("ttbar", TString::kIgnoreCase) )                        set = "ttbar";
    else if ( fname.Contains("dy", TString::kIgnoreCase) )                           set = "dy";
    else if ( fname.Contains("wjet", TString::kIgnoreCase) )                         set = "wjet";
    else if ( fname.Contains("st", TString::kIgnoreCase) ||
              fname.Contains("sat", TString::kIgnoreCase) )                          set = "st";
    else if ( fname.Contains("ww", TString::kIgnoreCase) ||
              fname.Contains("wz", TString::kIgnoreCase) ||
              fname.Contains("zz", TString::kIgnoreCase) )                           set = "vv";
    else if ( fname.Contains("qcd", TString::kIgnoreCase) &&
              ( (channel == "mm" &&  fname.Contains("Mu", TString::kIgnoreCase)) ||
                (channel == "ee" && !fname.Contains("Mu", TString::kIgnoreCase)) ||
                (channel == "em" &&  fname.Contains("Mu", TString::kIgnoreCase)) ) ) set = "wjet";
    else                                                                             continue;

    TFile* file = TFile::Open( dir + fname + "_" + channel + ".root" );
    if (!file) return -1;
    TH1D* h = (TH1D*) file->FindObjectAny(hname);
    if (!h) { cout << hname << " not found!" << endl;  return -1; }
    if ( theta!="gkk" && theta!="zp1" && theta!="zp10" && theta!="zp30" && (set=="zprime" || set=="gluon") ) h->Scale(sigScale);

    if (rebin.size() == 1) h->Rebin(rebin[0]);
    else                   h = (TH1D*) h->Rebin(rebin.size()-1, "h", &rebin[0]);

    if ( m_NOM.find(set) == m_NOM.end() ) m_NOM[set] = h;
    else                                  m_NOM[set]->Add(h);

    // Add up MC that are a part of bkg
    if ( isBkg(set) ) {
      if ( m_NOM.find("bkg") == m_NOM.end() ) m_NOM["bkg"] = (TH1D*) h->Clone("h_bkg"); // Need new hist since we use this pointer in two places
      else                                    m_NOM["bkg"]->Add(h);
    }

    if ( set=="DATA" || (postfilename != "" && postfilename != "false") ) continue;
    for (auto const& sys : systematics) {

      TH1D* h_UP, *h_DN;
      if ( ( (sys == "muTrigSys" || sys == "muIdSys") && channel == "ee" ) ||
           ( sys == "eleIdSys"   && channel == "mm" ) ||
           ( sys == "eleTrigSys" && channel != "ee" ) ||
           ( (sys == "topPtWeight" || sys == "pdf" || sys == "q2ttbar") && set != "ttbar") ) {

        h_UP = (TH1D*) h->Clone("h_UP");
        h_DN = (TH1D*) h->Clone("h_DN");
      }
      else {
        TFile* fileUP = TFile::Open( dir + fname + "_" + channel + "_" + sys + "UP.root" );
        h_UP = (TH1D*) fileUP->FindObjectAny(hname);
        if ( theta!="gkk" && theta!="zp1" && theta!="zp10" && theta!="zp30" && (set=="zprime" || set=="gluon") ) h_UP->Scale(sigScale);

        if (rebin.size() == 1) h_UP->Rebin(rebin[0]);
        else                   h_UP = (TH1D*) h_UP->Rebin(rebin.size()-1, "h_UP", &rebin[0]);

        if (sys == "topPtWeight" && set == "ttbar") {
          h_DN = (TH1D*) h->Clone("h_DN");
          h_DN->Scale(2);
          h_DN->Add(h_UP, -1);
        }
        else {
          TFile* fileDN = TFile::Open( dir + fname + "_" + channel + "_" + sys + "DOWN.root" );
          h_DN = (TH1D*) fileDN->FindObjectAny(hname);
          if ( theta!="gkk" && theta!="zp1" && theta!="zp10" && theta!="zp30" && (set=="zprime" || set=="gluon") ) h_DN->Scale(sigScale);

          if (rebin.size() == 1) h_DN->Rebin(rebin[0]);
          else                   h_DN = (TH1D*) h_DN->Rebin(rebin.size()-1, "h_DN", &rebin[0]);
        }
      }

      if ( m_UP[set].find(sys) == m_UP[set].end() ) m_UP[set][sys] = h_UP;
      else                                          m_UP[set][sys]->Add(h_UP);
      if ( m_DN[set].find(sys) == m_DN[set].end() ) m_DN[set][sys] = h_DN;
      else                                          m_DN[set][sys]->Add(h_DN);

      if ( isBkg(set) ) {
        if ( m_UP["bkg"].find(sys) == m_UP["bkg"].end() ) m_UP["bkg"][sys] = (TH1D*) h_UP->Clone("h_bkgUP");
        else                                              m_UP["bkg"][sys]->Add(h_UP);
        if ( m_DN["bkg"].find(sys) == m_DN["bkg"].end() ) m_DN["bkg"][sys] = (TH1D*) h_DN->Clone("h_bkgDN");
        else                                              m_DN["bkg"][sys]->Add(h_DN);
      }
    }
  }
  int nBins = m_NOM["DATA"]->GetNbinsX();

  /// THETA output ///

  if (theta == "zp1" || theta == "zp10" || theta == "zp30" || theta == "gkk") {
    TFile* outFile = new TFile(channel + region + "__" + theta + "__" + hname( hname.Index('_')+1, hname.Length()-hname.Index('_')-1 ) + ".root","RECREATE");
    outFile->cd();

    for (auto const& i_set : m_NOM) {
      TString set = i_set.first;
      if (set == "bkg") continue;
      i_set.second->Write( channel + region + "__" + set );
      for (int i=1; i<=nBins; i++) {
        if (i_set.second->GetBinError(i) / i_set.second->GetBinContent(i) > 0.2) cout << "Warning: " << set << " " << i_set.second->GetBinCenter(i) << endl;
      }

      if (set == "DATA") continue;
      for (auto const& sys : systematics) {
        if ( ( (sys == "muTrigSys" || sys == "muIdSys") && channel == "ee" ) ||
             ( sys == "eleIdSys"   && channel == "mm" ) ||
             ( sys == "eleTrigSys" && channel != "ee" ) ||
             ( (sys == "topPtWeight" || sys == "pdf" || sys == "q2ttbar") && set != "ttbar") ) continue;

        m_UP[set][sys]->Write( channel + region + "__" + set + "__" + sys + "__" + "up" );
        m_DN[set][sys]->Write( channel + region + "__" + set + "__" + sys + "__" + "down" );
      }
    }

    delete outFile;
    outFile = 0;
    return 1;
  }

  /// POST-FIT SYSTEMATICS ///

  if (postfilename != "" && postfilename != "false") {
    TFile* postfile = TFile::Open(postfilename + ".root");
    sys_norm.clear(); systematics.clear();

    TIter nextkey(postfile->GetListOfKeys());
    TKey* key;
    while ( (key = (TKey*)nextkey()) ) {
      TString keyname = key->GetName();

      // use systematics from post-fit file
      if (keyname.Contains("UP")) {
        int index = keyname.Index("__", 2, keyname.Index("__")+2, TString::kExact); // index of second "__"
        TString sys = keyname(index+2, keyname.Index("UP")-index-2);
        unsigned int idx;
        for (idx=0; idx<systematics.size(); idx++) { if (systematics[idx] == sys) break; }
        if (idx == systematics.size()) systematics.push_back(sys);
      }
    }
    // clone the binning
    for (auto const& sys : systematics) {
      m_UP["bkg"][sys] = (TH1D*) m_NOM["bkg"]->Clone(sys+"UP");
      m_DN["bkg"][sys] = (TH1D*) m_NOM["bkg"]->Clone(sys+"DN");
    }
    for (int i=1; i<=nBins; i++) {
      double bkg_content = 0;

      for (auto const& i_set : m_NOM) {
        TString set = i_set.first;
        if ( !isBkg(set) ) continue;

        TH1D* h = (TH1D*) postfile->FindObjectAny(channel + region + "__" + set);
        i_set.second->SetBinContent(i, h->GetBinContent(i));
        bkg_content += h->GetBinContent(i);
      }
      m_NOM["bkg"]->SetBinContent(i, bkg_content);

      for (auto const& sys : systematics) {
        double bkgUP_content=0, bkgDN_content=0;

        for (auto const& i_set : m_NOM) {
          TString set = i_set.first;
          if ( !isBkg(set) ) continue;

          TH1D* h_UP = (TH1D*) postfile->FindObjectAny(channel + region + "__" + set + "__" + sys + "UP");
          TH1D* h_DN = (TH1D*) postfile->FindObjectAny(channel + region + "__" + set + "__" + sys + "DN");
          bkgUP_content += h_UP->GetBinContent(i);  bkgDN_content += h_DN->GetBinContent(i);
        }
        m_UP["bkg"][sys]->SetBinContent(i, bkgUP_content);
        m_DN["bkg"][sys]->SetBinContent(i, bkgDN_content);
      }
    }
  }

  /// Calculate background systematics ///

  TGraphAsymmErrors* background = new TGraphAsymmErrors(); // TGraph background with full systematics for plotting         (UP and DOWN error)
  m_NOM["sys"] = (TH1D*) m_NOM["bkg"]->Clone("h_sys");     // TH1    background with full systematics for KS and Chi2 test (average error)
  for (auto const& i_norm : sys_norm) systematics.push_back( i_norm.first );

  for (int pt=0; pt<nBins; pt++) {
    int bin = pt+1;
    double nom = m_NOM["bkg"]->GetBinContent(bin), staterr = m_NOM["bkg"]->GetBinError(bin), errorUP=0, errorDN=0, errorAVG=0;

    background->SetPoint(pt, m_NOM["bkg"]->GetBinCenter(bin), nom);
    if (nom == 0) continue;

    for (auto const& sys : systematics) {
      double deltaUP, deltaDN;

      if (sys_norm.find(sys) != sys_norm.end()) {  // normalization-only systematics
        double norm = sys_norm[sys];

        if      (sys == "lumi") norm *= nom;
        else if (sys=="sig_dy") norm *= m_NOM["dy"]->GetBinContent(bin);
        else if (sys=="sig_st") norm *= m_NOM["st"]->GetBinContent(bin);
        else if (sys=="sig_db") norm *= m_NOM["vv"]->GetBinContent(bin);
        else                    norm = 0;

        deltaUP = norm;
        deltaDN = -1 * norm;
      }
      else {                                       // shape systematics
        deltaUP = m_UP["bkg"][sys]->GetBinContent(bin) - nom;
        deltaDN = m_DN["bkg"][sys]->GetBinContent(bin) - nom;
      }
      double d1 = delta(deltaUP, deltaDN, "+");
      double d2 = delta(deltaUP, deltaDN, "-");

      errorUP  += d1*d1;
      errorDN  += d2*d2;

      errorAVG += 0.5*(fabs(d1)+fabs(d2)) * 0.5*(fabs(d1)+fabs(d2));
    }
    background->SetPointEXhigh( pt, m_NOM["bkg"]->GetBinWidth(bin)/2 );
    background->SetPointEXlow(  pt, m_NOM["bkg"]->GetBinWidth(bin)/2 );
    background->SetPointEYhigh( pt, sqrt( errorUP + staterr*staterr ) );
    background->SetPointEYlow(  pt, sqrt( errorDN + staterr*staterr ) );

    m_NOM["sys"]->SetBinError( bin, sqrt( errorAVG + staterr*staterr ) );
  }

  /// PLOTTING ///

  m_NOM["DATA"]->SetMarkerStyle(20);
  m_NOM["DATA"]->SetLineColor(kBlack);
  m_NOM["DATA"]->SetMarkerSize(0.45);

  background->SetFillStyle(3245);
  background->SetFillColor(LightGray);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(1);

  m_NOM["ttbar"]->SetLineColor(2);
  m_NOM["ttbar"]->SetFillColor(2);
  m_NOM["dy"]->SetLineColor(8);
  m_NOM["dy"]->SetFillColor(8);
  m_NOM["wjet"]->SetLineColor(4);
  m_NOM["wjet"]->SetFillColor(4);
  m_NOM["vv"]->SetLineColor(6);
  m_NOM["vv"]->SetFillColor(6);
  m_NOM["st"]->SetLineColor(28);
  m_NOM["st"]->SetFillColor(28);

  m_NOM["zprime"]->SetLineColor(12);
  m_NOM["zprime"]->SetLineWidth(2);
  m_NOM["gluon"]->SetLineColor(9);
  m_NOM["gluon"]->SetLineWidth(2);
  m_NOM["gluon"]->SetLineStyle(2);

  THStack* mcStack = new THStack();
  for (auto const& i_set : labels) {
    TString set = i_set.first;
    if ( isBkg(set) ) mcStack->Add( m_NOM[set] );
  }

  TCanvas* c = new TCanvas("c", "c", 400, 400);
  float b_scale = 0.34, t_scale = 0.99 - b_scale;
  TPad* top = new TPad("top", "top", 0.01, b_scale, 0.99, 0.99);
  TPad* bottom = new TPad("bottom", "bottom", 0.01, 0.0, 0.99, b_scale);

  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    top->SetTopMargin(0.065);
    top->SetBottomMargin(0.03);
    top->SetLeftMargin(0.19);
    top->SetRightMargin(0.05);
    top->Draw();
    bottom->SetTopMargin(0.03);
    bottom->SetBottomMargin(0.35);
    bottom->SetLeftMargin(0.19);
    bottom->SetRightMargin(0.05);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }
  if (logx) gPad->SetLogx();
  if (logy) gPad->SetLogy();

  TH1D* hist=0;
  if (rebin.size() == 1) hist = new TH1D("hist", "hist", nBins, m_NOM["DATA"]->GetBinLowEdge(1), m_NOM["DATA"]->GetBinLowEdge(nBins+1));
  else                   hist = new TH1D("hist", "hist", rebin.size()-1, &rebin[0]);

  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    hist->GetXaxis()->SetTickLength(0.03/t_scale);
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.05/t_scale);
    hist->GetYaxis()->SetTitleOffset(1.5*t_scale);
    hist->GetYaxis()->SetLabelSize(0.03/t_scale);
  }
  else {
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetXaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetTitleOffset(1.2);

    if (logx) { hist->GetXaxis()->SetNoExponent(); hist->GetXaxis()->SetMoreLogLabels(); }
  }

  hist->GetXaxis()->SetNdivisions(5, 5, 0);
  hist->GetXaxis()->SetRangeUser(xmin, xmax);
  hist->GetYaxis()->SetTitle("Events");

  ymax = mcStack->GetMaximum();
  if (logy) { ymax *= 100; ymin = 0.1; }
  else        ymax *= 1.5;
  hist->GetYaxis()->SetRangeUser(ymin, int(ymax) );

  hist->Draw();
  int legEntries = labels.size(); // don't include bkg uncert in legend
  if (!plotData) legEntries--;

  TLegend* leg = 0;
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    leg = new TLegend(.65,.9-.055*legEntries,.85,.9);
    leg->SetTextSize(0.04);
  }
  else {
    leg = new TLegend(.65,.9-.04*legEntries,.85,.9);
    leg->SetTextSize(0.03);
  }
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  if (plotData) leg->AddEntry(m_NOM["DATA"], "Data", "PLE");

  for (auto const& i_set : labels) {
    TString set = i_set.first;
    if (set=="bkg") continue;

    leg->AddEntry( m_NOM[set], i_set.second.Data(), (set=="gluon" || set=="zprime" ? "L" : "F") );
  }

  mcStack->Draw("samehist");
  background->Draw("sameE2");
  m_NOM["gluon"]->Draw("samehist");
  m_NOM["zprime"]->Draw("samehist");
  if (plotData) m_NOM["DATA"]->Draw("samePE");

  gPad->RedrawAxis();
  leg->Draw();

  drawText();

  TH1D* bhist = 0;
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    bottom->cd();
    if (logx) gPad->SetLogx();

    bhist = (TH1D*) hist->Clone("bhist");

    TH1D* hsubplot = (TH1D*) m_NOM["DATA"]->Clone("hsubplot");

    TGraphAsymmErrors* bkg_env = (TGraphAsymmErrors*) background->Clone("bkg_env");
    TGraphAsymmErrors* bkg_env_stat = (TGraphAsymmErrors*) background->Clone("bkg_env_stat");
    double *yarray = background->GetY(), *eyhigh = background->GetEYhigh(), *eylow = background->GetEYlow();

    TString subytitle = "Data / Bkg";
    if (subplot=="ratio") {
      for (int pt=0; pt<nBins; pt++) {
        if (yarray[pt] == 0) continue;
        int bin=pt+1;

        hsubplot->SetBinContent(bin, m_NOM["DATA"]->GetBinContent(bin)/yarray[pt]);
        hsubplot->SetBinError(  bin, m_NOM["DATA"]->GetBinError(bin)  /yarray[pt]);

        bkg_env->SetPoint(       pt, m_NOM["DATA"]->GetBinCenter(bin), 1.);
        bkg_env->SetPointEYhigh( pt, eyhigh[pt]/yarray[pt] );
        bkg_env->SetPointEYlow(  pt, eylow[pt] /yarray[pt] );

        bkg_env_stat->SetPoint(       pt, m_NOM["DATA"]->GetBinCenter(bin), 1.);
        bkg_env_stat->SetPointEYhigh( pt, m_NOM["bkg"]->GetBinError(bin)/yarray[pt] );
        bkg_env_stat->SetPointEYlow(  pt, m_NOM["bkg"]->GetBinError(bin)/yarray[pt] );
      }
    }
    else {
      subytitle = "Pull";
      for (int pt=0; pt<nBins; pt++) {
        if (yarray[pt] == 0) continue;
        int bin=pt+1;

        double diff = m_NOM["DATA"]->GetBinContent(bin)-yarray[pt];
        double sigma = diff>0 ? sqrt( m_NOM["DATA"]->GetBinError(bin)*m_NOM["DATA"]->GetBinError(bin) + eyhigh[pt]*eyhigh[pt] )
                              : sqrt( m_NOM["DATA"]->GetBinError(bin)*m_NOM["DATA"]->GetBinError(bin) + eylow[pt]*eylow[pt] );

        if (sigma == 0) continue;
        hsubplot->SetBinContent(bin, diff/sigma);
      }
    }

    bhist->GetXaxis()->SetLabelSize(0.03/b_scale);
    bhist->GetXaxis()->SetTickLength(0.03/b_scale);
    bhist->GetXaxis()->SetTitleSize(0.05/b_scale);
    bhist->GetXaxis()->SetTitleOffset(1);
    bhist->GetXaxis()->SetTitle(xtitle);
    bhist->GetXaxis()->SetRangeUser(xmin, xmax);
    if (logx) { bhist->GetXaxis()->SetNoExponent(); bhist->GetXaxis()->SetMoreLogLabels(); }

    bhist->GetYaxis()->SetRangeUser(subymin, subymax);
    bhist->GetYaxis()->SetNdivisions(5, 3, 0);
    bhist->GetYaxis()->SetLabelSize(0.03/b_scale);
    bhist->GetYaxis()->SetTitleOffset(1.5*b_scale);
    bhist->GetYaxis()->SetTitle(subytitle);
    bhist->GetYaxis()->SetTitleSize(0.05/b_scale);
    bhist->GetYaxis()->CenterTitle(true);

    bhist->Draw();

    if (subplot=="ratio") {
      bkg_env_stat->SetLineColor(MLightGray);
      bkg_env_stat->SetFillColor(MLightGray);
      bkg_env_stat->SetFillStyle(1001);

      bkg_env->SetLineColor(LightGray);
      bkg_env->SetFillColor(LightGray);
      bkg_env->SetFillStyle(1001);

      bkg_env->Draw("sameE2");
      bkg_env_stat->Draw("sameE2");
      hsubplot->Draw("samePE");
      TLine* hline = new TLine(xmin,1,xmax,1);
      hline->SetLineStyle(7);
      hline->Draw();
    }
    else {
      hsubplot->SetFillColor(kGray+2);
      hsubplot->SetFillStyle(3144);
      hsubplot->Draw("samehist");
    }
    gPad->RedrawAxis();
  }

  c->Print( outName + ".pdf" );

  if (fit) {
    cout << channel << " " << hname << endl;
    cout << "KS Test\t\t"        << m_NOM["DATA"]->KolmogorovTest(m_NOM["bkg"])      << "\t" << m_NOM["DATA"]->KolmogorovTest(m_NOM["sys"])      << endl;
    cout << "KS Test (with X)\t" << m_NOM["DATA"]->KolmogorovTest(m_NOM["bkg"], "X") << "\t" << m_NOM["DATA"]->KolmogorovTest(m_NOM["sys"], "X") << endl;
    m_NOM["DATA"]->Chi2Test(m_NOM["bkg"], "UWP");
    m_NOM["DATA"]->Chi2Test(m_NOM["sys"], "UWP");
  }

  if (plotImpact) {

    delete leg;
    TLegend* leg = new TLegend(.65,.9-.06*4,.85,.9);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    bhist->GetYaxis()->SetTitleOffset(0.53);
    if (subplot=="ratio") bhist->GetYaxis()->SetTitle("#frac{Up/Down}{Nom}");
    else                  bhist->GetYaxis()->SetTitle("#frac{Up/Down - Nom}{Nom}");

    for (auto const& sys : systematics) {
      if (sys_norm.find(sys) != sys_norm.end()) continue;

      for (auto const& i_set : labels) {
        TString set = i_set.first;

        c->Clear("D");  leg->Clear();
        top->cd();  hist->Draw();  drawText();

        TH1D* h_NM = m_NOM[set];
        TH1D* h_UP = m_UP[set][sys];
        TH1D* h_DN = m_DN[set][sys];

        if (logy) hist->GetYaxis()->SetRangeUser(0.1, int(h_NM->GetMaximum()*100) );
        else      hist->GetYaxis()->SetRangeUser(ymin, int(h_NM->GetMaximum()*1.5) );

        h_NM->ResetAttFill();  h_NM->ResetAttLine();

        h_NM->SetLineColor(kBlack);  h_NM->SetLineWidth(2);  h_NM->Draw("histsame");
        h_UP->SetLineColor(8);       h_UP->SetLineWidth(2);  h_UP->Draw("histsame");
        h_DN->SetLineColor(kRed+1);  h_DN->SetLineWidth(2);  h_DN->Draw("histsame");

        leg->SetHeader( i_set.second.Data() );  leg->AddEntry(h_NM, "Nominal", "L");
        leg->AddEntry(h_UP, sys+" Up", "L");    leg->AddEntry(h_DN, sys+" Down", "L");
        leg->Draw();

        bottom->cd();  bhist->Draw();

        TH1D* hsubUP = (TH1D*) h_UP->Clone("hsubUP");
        TH1D* hsubDN = (TH1D*) h_DN->Clone("hsubDN");
        if (subplot=="ratio") { hsubUP->Divide(h_NM); hsubDN->Divide(h_NM); }
        else                  { hsubUP->Add(h_NM, -1); hsubUP->Divide(h_NM); hsubDN->Add(h_NM, -1); hsubDN->Divide(h_NM); }
        hsubUP->Draw("histsame");  hsubDN->Draw("histsame");

        c->Print("./impact/" + outName + "_" + set + "_" + sys + ".pdf");
      }
    }
  }

}

void drawText() {
  TLatex text;
  text.SetNDC();

  text.SetTextSize(0.045);
  text.SetTextFont(42);
  if ( plotData && (subplot=="ratio" || subplot=="pull") )
    text.DrawLatex(1-rightText.Length()/82., 0.95, rightText);
  else {
    text.SetTextSize(0.035);
    text.DrawLatex(1-rightText.Length()/75., 0.96, rightText);
  }

  text.SetTextSize(0.06);
  text.SetTextFont(61);
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) text.DrawLatex(0.2, 0.95, "CMS");
  else {
    text.SetTextSize(0.05);
    text.DrawLatex(0.18, 0.96, "CMS");
  }

  //text.SetTextSize(0.04);
  //text.SetTextFont(52);
  //text.DrawLatex(0.29, 0.96, "Simulation Preliminary"); //make bool

  text.SetTextSize(0.04);
  text.SetTextFont(42);
  float textposx, textposy, gap;
  if ( plotData && (subplot=="ratio" || subplot=="pull") ) {
    textposx = 0.25;
    textposy = 0.88;
    gap = 0.06;
  }
  else {
    text.SetTextSize(0.03);
    textposx = 0.22;
    textposy = 0.9;
    gap = 0.05;
  }

  if      (channel == "mm")  text.DrawLatex(textposx,textposy,"#bf{#mu#mu}");
  else if (channel == "ee")  text.DrawLatex(textposx,textposy,"#bf{ee}");
  else                       text.DrawLatex(textposx,textposy,"#bf{e#mu}");

  if (hname.Contains("0_") || hname.Contains("2_") || hname.Contains("4_") || hname.Contains("6_"))
    text.DrawLatex(textposx,textposy-gap,"#bf{= 0 btags}");
  else
    text.DrawLatex(textposx,textposy-gap,"#bf{#geq 1 btag}");

  if (hname.Contains("0_") || hname.Contains("1_") || hname.Contains("4_") || hname.Contains("5_"))
    text.DrawLatex(textposx,textposy-2*gap,"#bf{p_{T}^{j0}>100 GeV}");
  else
    text.DrawLatex(textposx,textposy-2*gap,"#bf{p_{T}^{j0}>100 GeV, p_{T}^{j1}>50 GeV}");

  if (hname.Contains("0_") || hname.Contains("1_") || hname.Contains("2_") || hname.Contains("3_"))
    text.DrawLatex(textposx,textposy-3*gap,"#bf{p_{T}^{miss} < 30 GeV}");
  else
    text.DrawLatex(textposx,textposy-3*gap,"#bf{p_{T}^{miss} > 30 GeV}");
}

void setPars(const string& parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)) {

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if (var == "fileNames") {
      while ( (delim_pos = line.find(' ')) != -1) {
        fileNames.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      fileNames.push_back( line.data() );
    }
    else if (var == "systematics") {
      while ( (delim_pos = line.find(' ')) != -1) {
        systematics.push_back( line.substr(0, delim_pos).data() );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      systematics.push_back( line.data() );
    }
    else if (var == "sys_norm") {
      while ( (delim_pos = line.find(' ')) != -1) {
        int col = line.find(':');
        sys_norm[line.substr(0, col).data()] = stof( line.substr(col+1, delim_pos-col-1) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      int col = line.find(':');
      sys_norm[line.substr(0, col).data()] = stof( line.substr(col+1, delim_pos-col-1) );
    }
    else if (var == "rebin") {
      while ( (delim_pos = line.find(' ')) != -1) {
        rebin.push_back( stod( line.substr(0, delim_pos) ) );
        line.erase(0, delim_pos + 1);
        while (line.at(0) == ' ') line.erase(0, 1);
      }
      rebin.push_back( stod(line) );
    }
    else if (var == "dataFileName") dataFileName = line.data();
    else if (var == "outName")      outName = line.data();
    else if (var == "dir")          dir = line.data();
    else if (var == "channel")      channel = line.data();
    else if (var == "theta")        theta = line.data();
    else if (var == "postfilename") postfilename = line.data();
    else if (var == "hname")        hname = line.data();
    else if (var == "rightText")    rightText = line.data();
    else if (var == "subplot")      subplot = line.data();
    else if (var == "xtitle")       xtitle = line.data();
    else if (var == "zprime")       zprime = line.data();
    else if (var == "gluon")        gluon = line.data();
    else if (var == "sigScale")     sigScale = stof(line);
    else if (var == "xmin")         xmin = stof(line);
    else if (var == "xmax")         xmax = stof(line);
    else if (var == "ymin")         ymin = stof(line);
    else if (var == "ymax")         ymax = stof(line);
    else if (var == "subymin")      subymin = stof(line);
    else if (var == "subymax")      subymax = stof(line);
    else if (var == "region")       { if (line != "false") region = line.data(); }
    else if (var == "logx")         { if (line == "true")  logx = true; }
    else if (var == "logy")         { if (line == "true")  logy = true; }
    else if (var == "plotData")     { if (line == "true")  plotData = true; }
    else if (var == "plotImpact")   { if (line == "true")  plotImpact = true; }
    else if (var == "fit")          { if (line == "true")  fit = true; }
  }
  file.close();
}

void setStyle(){
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(0.9);

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

double delta(double delta1, double delta2, string plus_minus) {
// pos_neg should be either "+" or "-"

  double del = 0. ;
  if (plus_minus=="+") {
    double dd = (delta1 > delta2) ? delta1 : delta2 ;
    if (dd > 0.) del = dd ;
  }
  else if (plus_minus=="-") {
    double dd = (delta1 < delta2) ? delta1 : delta2 ;
    if (dd < 0.) del = dd ;
  }
  return del;
}
