//Only statistical error:    createLog OFF/logs/
//Systematic and stat error: createLog dir cutname {systematics, ...}
//Example: createLog logs/ \>=_2_Jets,\_\>=_1_btag lumi sig_dy sig_st sig_db jec jer btagSF mistagSF pileup topPtWeight pdf q2ttbar muTrigSys muIdSys eleTrigSys eleIdSys

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <iomanip>
#include <string>

#include "TString.h"

using namespace std;

void readFile(const string& fileName, map<string, map<string, pair<double, double> > >& cuts);
double delta(double delta1, double delta2, string pos_neg);

string boldline   = "====================================================================================================================================================================";
string singleline = "--------------------------------------------------------------------------------------------------------------------------------------------------------------------";

string zprime = "ZprimeToTT_M3000_W300_total", gluon = "RSGluonToTT_M-3000";   //specify which signals to use in cutflow summary
double xs_zprime = 0.272788, xs_gluon = 0.16757;

//use vector to define order of samples in output tables
vector< pair<string, string> > set_labels = { {"ttbar","ttbar"}, {"dy","Drell-Yan"}, {"st","Single-Top"}, {"vv","Diboson"}, {"wjet","W+Jets"},
                                              {"bkg","Background"}, {zprime,zprime}, {gluon,gluon} };

vector< pair<string, string> > set_latex = { {"ttbar","t$\\bar{\\textrm{t}}$"}, {"dy","Z/$\\gamma^{*}\\rightarrow l^{+}l^{-}$"}, {"st","Single-Top"}, {"vv","VV"}, {"wjet","W+Jets"},
                                             {"bkg","Total Bkg"}, {zprime,"Z$'$ (10\\%, 3 TeV)"}, {gluon,"$\\textrm{g}_{\\textrm{kk}}$ (3 TeV)"} };

map<string, string> sys_labels = { {"lumi","luminosity"}, {"sig_dy","\u03C3(drell-yan)"}, {"sig_st","\u03C3(single-top)"}, {"sig_db","\u03C3(diboson)"},
                                   {"jec","JEC"}, {"jer","JER"}, {"btagSF","b-tagging"}, {"mistagSF","mis-tagging"}, {"pileup","pileup"},
                                   {"topPtWeight","top pT modeling"}, {"pdf","PDF"}, {"q2ttbar","Q2 ttbar"},
                                   {"muTrigSys","muon trigger"}, {"muIdSys","muon id"}, {"eleTrigSys","electron trigger"}, {"eleIdSys","electron id"} };

map<string, string> sys_latex = { {"lumi","luminosity"}, {"sig_dy","$\\sigma(drell-yan)$"}, {"sig_st","$\\sigma(single-top)$"}, {"sig_db","$\\sigma(diboson)$"},
                                  {"jec","JEC"}, {"jer","JER"}, {"btagSF","b-tagging"}, {"mistagSF","mis-tagging"}, {"pileup","pileup"},
                                  {"topPtWeight","top $p_{T}$ modeling"}, {"pdf","PDF"}, {"q2ttbar","Q2 ttbar"},
                                  {"muTrigSys","muon trigger"}, {"muIdSys","muon id"}, {"eleTrigSys","electron trigger"}, {"eleIdSys","electron id"} };

// Normalization-only systematics
map<string, double> sys_norm = { {"lumi",0.025}, {"sig_dy",0.2}, {"sig_st",0.16}, {"sig_db",0.15} };

int main(int argc, char* argv[]) {

  if      (argc == 1) { cout << "Please provide a directory." << endl; return -1; }
  else if (argc == 3) { cout << "Please provide a cut name (using '_' in place of spaces) followed by systematics." << endl; return -1; }

  string dir = argv[1];

  TString summary_cutname;
  vector<string> systematics;
  if (argc > 3) {
    summary_cutname = argv[2];  summary_cutname.ReplaceAll("_", " ");

    //shape and norm systematics: Read from command line
    for (int i=3; i<argc; i++) systematics.push_back( argv[i] );
  }

  //map(channel, map(cut, map(dataset, (N, stat error) ) ) )          map(channel, map(cut, map(dataset, (sysUP, sysDN) ) ) )
  map<string, map<string, map<string, pair<double, double> > > > nom, totalSys;

  //map(channel, map(sys, map(cut, map(dataset, (N, stat error) ) ) ) )
  map<string, map<string, map<string, map<string, pair<double, double> > > > > delUP, delDN;

  vector<string> channels = { "mm", "ee", "em", "ll" };

  //read in event yields
  for (auto const& chan : channels) {
    if (chan == "ll") continue;

    readFile(dir + "/log_" + chan + ".txt", nom[chan]);

    for (auto const& sys : systematics) {
      if      ( (sys == "muTrigSys" || sys == "muIdSys") && chan == "ee" ) continue;
      else if (  sys == "eleIdSys"   && chan == "mm" ) continue;
      else if (  sys == "eleTrigSys" && chan != "ee" ) continue;

      if (sys_norm.find(sys) == sys_norm.end()) {  // shape systematics
        readFile(dir + "/log_" + chan + "_" + sys + "UP.txt",   delUP[chan][sys]);
        readFile(dir + "/log_" + chan + "_" + sys + "DOWN.txt", delDN[chan][sys]);
      }
      for (auto const& it_cut : nom[chan]) {
        string cut = it_cut.first;

        for (auto const& it_set : nom[chan][cut]) {
          string set = it_set.first;

          if (sys_norm.find(sys) == sys_norm.end()) {  // shape systematics

            if (sys == "topPtWeight" || sys == "pdf" || sys == "q2ttbar") {
              if (set == "bkg" || set == "ttbar") {    // this works because "bkg" comes before "ttbar" in set map

                if ( sys == "topPtWeight" ) delDN[chan][sys][cut][set].first = nom[chan][cut]["ttbar"].first - delUP[chan][sys][cut]["ttbar"].first;
                else                        delDN[chan][sys][cut][set].first = delDN[chan][sys][cut]["ttbar"].first - nom[chan][cut]["ttbar"].first;

                delUP[chan][sys][cut][set].first = delUP[chan][sys][cut]["ttbar"].first - nom[chan][cut]["ttbar"].first;
              }
              else {
                delUP[chan][sys][cut][set].first = 0;
                delDN[chan][sys][cut][set].first = 0;
              }
            }
            else {
              delUP[chan][sys][cut][set].first -= nom[chan][cut][set].first;
              delDN[chan][sys][cut][set].first -= nom[chan][cut][set].first;
            }
          }
          else {                                       // normalization-only systematics
            double perEvent_sys = sys_norm[sys];

            if      (sys == "lumi")
              perEvent_sys *= nom[chan][cut][set].first;
            else if ( sys=="sig_dy" && (set == "dy" || set == "bkg") )
              perEvent_sys *= nom[chan][cut]["dy"].first;
            else if ( sys=="sig_st" && (set == "st" || set == "bkg") )
              perEvent_sys *= nom[chan][cut]["st"].first;
            else if ( sys=="sig_db" && (set == "vv" || set == "bkg") )
              perEvent_sys *= nom[chan][cut]["vv"].first;
            else perEvent_sys = 0.;

            delUP[chan][sys][cut][set].first = perEvent_sys;
            delDN[chan][sys][cut][set].first = perEvent_sys==0 ? 0 : -1 * perEvent_sys;
          }
        }
      }
    }
  }

  //fill combined channel
  for (auto const& it_cut : nom["mm"]) {
    string cut = it_cut.first;

    for (auto const& it_set : nom["mm"][cut]) {
      string set = it_set.first;

      if (cut <= "FMET Filters") {
        if (set == "data") {
          nom["ll"][cut][set].first =  nom["mm"][cut][set].first + nom["ee"][cut][set].first;
          nom["ll"][cut][set].second = sqrt( nom["mm"][cut][set].second*nom["mm"][cut][set].second + nom["ee"][cut][set].second*nom["ee"][cut][set].second );
        }
        else {
          nom["ll"][cut][set].first =  nom["mm"][cut][set].first;
          nom["ll"][cut][set].second = nom["mm"][cut][set].second;
        }
      }
      else {
        nom["ll"][cut][set].first =  nom["mm"][cut][set].first + nom["ee"][cut][set].first + nom["em"][cut][set].first;
        nom["ll"][cut][set].second = sqrt( nom["mm"][cut][set].second*nom["mm"][cut][set].second + nom["ee"][cut][set].second*nom["ee"][cut][set].second +
                                           nom["em"][cut][set].second*nom["em"][cut][set].second );
      }
      for (auto const& sys : systematics) {

        if (cut <= "FMET Filters") {
          if ( sys == "eleTrigSys" || sys == "eleIdSys" ) {
            delUP["ll"][sys][cut][set].first = delUP["ee"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["ee"][sys][cut][set].first;
          }
          else {
            delUP["ll"][sys][cut][set].first = delUP["mm"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["mm"][sys][cut][set].first;
          }
        }
        else {
          if      ( sys == "eleTrigSys" ) {
            delUP["ll"][sys][cut][set].first = delUP["ee"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["ee"][sys][cut][set].first;
          }
          else if ( sys == "eleIdSys" ) {
            delUP["ll"][sys][cut][set].first = delUP["ee"][sys][cut][set].first + delUP["em"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["ee"][sys][cut][set].first + delDN["em"][sys][cut][set].first;
          }
          else if ( sys == "muTrigSys" || sys == "muIdSys" ) {
            delUP["ll"][sys][cut][set].first = delUP["mm"][sys][cut][set].first + delUP["em"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["mm"][sys][cut][set].first + delDN["em"][sys][cut][set].first;
          }
          else {
            delUP["ll"][sys][cut][set].first = delUP["mm"][sys][cut][set].first + delUP["ee"][sys][cut][set].first + delUP["em"][sys][cut][set].first;
            delDN["ll"][sys][cut][set].first = delDN["mm"][sys][cut][set].first + delDN["ee"][sys][cut][set].first + delDN["em"][sys][cut][set].first;
          }
        }
      }
    }
  }
/*
  // diff between CR6 and CR5
  for (auto const& chan : channels) {

    for (auto const& it_set : nom[chan]["W>= 2 Jets, = 0 btags"]) {
      string set = it_set.first;
      nom[chan]["YCR6 - CR5"][set].first = nom[chan]["W>= 2 Jets, = 0 btags"][set].first - nom[chan]["V= 1 Jet, >= 1 btag"][set].first;

      nom[chan]["YCR6 - CR5"][set].second = sqrt( nom[chan]["W>= 2 Jets, = 0 btags"][set].second*nom[chan]["W>= 2 Jets, = 0 btags"][set].second
                                                + nom[chan]["V= 1 Jet, >= 1 btag"][set].second*nom[chan]["V= 1 Jet, >= 1 btag"][set].second );

      for (auto const& it_sys : delUP[chan]) {
        string sys = it_sys.first;

        delUP[chan][sys]["YCR6 - CR5"][set].first = delUP[chan][sys]["W>= 2 Jets, = 0 btags"][set].first - delUP[chan][sys]["V= 1 Jet, >= 1 btag"][set].first;
        delDN[chan][sys]["YCR6 - CR5"][set].first = delDN[chan][sys]["W>= 2 Jets, = 0 btags"][set].first - delDN[chan][sys]["V= 1 Jet, >= 1 btag"][set].first;
      }
    }
  }
*/
  string cut_initial = nom["mm"].begin()->first, cut_end = (--nom["mm"].end())->first;

  // create cutflow tables
  for (auto const& chan : channels) {
    ofstream file(dir + chan + "_cutflow.txt");

    /// SUMMARY TABLE ///
    file << boldline << endl;
    file << "                                              Cut Flow Table: Summary\n" ;
    file << boldline << endl << Form("%27s", "");

    for (auto const& p : set_labels) {
      string set = p.first;
      const char* label = p.second.data();  int len = strlen(label);

      if      (set == "bkg")                  file << "|||            Data           |||         Background        |||      Data/Background      ";
      else if (set == zprime || set == gluon) file << Form("|||%*s%*s",14+len/2,label,13-len/2,"");
    }
    for (auto const& it_cut : nom[chan]) {
      string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);

      file << endl << Form("%-27s", cutname.data());
      for (auto const& p : set_labels) {
        string set = p.first;

        if      (set == "bkg")
          file << Form("|||  %12.0f (%1.6f)  |||  %12.1f (%1.6f)  |||           %1.3f           ",
                        nom[chan][cut]["data"].first,  nom[chan][cut]["data"].first / nom[chan][cut_initial]["data"].first,
                        nom[chan][cut]["bkg"].first,   nom[chan][cut]["bkg"].first  / nom[chan][cut_initial]["bkg"].first,
                        nom[chan][cut]["data"].first / nom[chan][cut]["bkg"].first);
        else if (set == zprime || set == gluon)
          file << Form("|||  %12.1f (%1.6f)  ", nom[chan][cut][set].first, nom[chan][cut][set].first / nom[chan][cut_initial][set].first);
      }
      if (cutname == "MET Filters" || cutname == ">= 1 jet") file << endl << singleline;
    } // end cut loop 1 (Summary table)

    /// BACKGROUND TABLE ///
    file << endl << endl << boldline << endl;
    file << "                                              Cut Flow Table: Background\n" ;
    file << boldline << endl << Form("%27s", "");

    for (auto const& p : set_labels) {
      string set = p.first;
      if (set == "bkg" || set == zprime || set == gluon) continue;
      const char* label = p.second.data();  int len = strlen(label);

      file << Form("|||%*s%*s",14+len/2,label,13-len/2,"");
    }
    for (auto const& it_cut : nom[chan]) {
      string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);

      file << endl << Form("%-27s", cutname.data());
      for (auto const& p : set_labels) {
        string set = p.first;
        if (set == "bkg" || set == zprime || set == gluon) continue;

        file << Form("|||  %12.1f (%1.6f)  ", nom[chan][cut][set].first, nom[chan][cut][set].first / nom[chan][cut_initial][set].first);
      }
      if (cutname == "MET Filters" || cutname == ">= 1 jet") file << endl << singleline;
    } // end cut loop 2 (Background table)

    /// INDIVIDUAL SYSTEMATICS ///
    for (auto const& sys : systematics) {
      if ( delUP[chan].find(sys) == delUP[chan].end() ) continue;

      file << endl << endl << boldline << endl;
      file << "                                    Systematics: " << sys_labels[sys] <<  "       [absolute UP DOWN on first line and relative UP DOWN in \% on second line]\n";
      file << boldline << endl << Form("%27s", "");

      for (auto const& p : set_labels) {
        const char* label = p.second.data();  int len = strlen(label);
        file << Form("||%*s%*s",11+len/2,label,10-len/2,"");
      }

      /// CALCULATE TOTAL SYS ///
      for (auto const& it_cut : nom[chan]) {
        string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);

        for (auto const& it_set : nom[chan][cut]) {
          string set = it_set.first;

          double d1 = delta( delUP[chan][sys][cut][set].first, delDN[chan][sys][cut][set].first, "+" );
          double d2 = delta( delUP[chan][sys][cut][set].first, delDN[chan][sys][cut][set].first, "-" );

          totalSys[chan][cut][set].first  += d1*d1;
          totalSys[chan][cut][set].second += d2*d2;
        }

        file << endl << Form("%-27s", cutname.data());
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", delUP[chan][sys][cut][set].first, delDN[chan][sys][cut][set].first);
        }
        file << endl << Form("%-27s", "");
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", 100*delUP[chan][sys][cut][set].first/nom[chan][cut][set].first, 100*delDN[chan][sys][cut][set].first/nom[chan][cut][set].first);
        }
        if (cutname == "MET Filters" || cutname == ">= 1 jet") file << endl << singleline;
      }  // end cut loop 3
    } // end systematics loop 1 (Individual systematics)

    if (systematics.size() > 0) {
      string cutkey;

      /// TOTAL SYSTEMATICS CUTFLOW ///
      file << endl << endl << boldline << endl;
      file << "                                    Total Systematics:        [absolute +- on first line and relative +- in \% on second line]\n";
      file << boldline << endl << Form("%27s", "");

      for (auto const& p : set_labels) {
        const char* label = p.second.data();  int len = strlen(label);
        file << Form("||%*s%*s",11+len/2,label,10-len/2,"");
      }
      for (auto const& it_cut : nom[chan]) {
        string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);
        if (summary_cutname == cutname) cutkey = cut;

        file << endl << Form("%-27s", cutname.data());
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", sqrt( totalSys[chan][cut][set].first ), sqrt( totalSys[chan][cut][set].second ) );
        }
        file << endl << Form("%-27s", "");
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", 100*sqrt(totalSys[chan][cut][set].first)/nom[chan][cut][set].first, 100*sqrt(totalSys[chan][cut][set].second)/nom[chan][cut][set].first);
        }
        if (cutname == "MET Filters" || cutname == ">= 1 jet") file << endl << singleline;
      } // end cut loop 4 (Total Systematics)

      /// TOTAL SYSTEMATICS SUMMARY ///
      string cutname = cutkey;  cutname.erase(0, 1);
      file << endl << endl << boldline << endl;
      file << "           Total Systematics: " + cutname + "        [absolute +- on first line and relative +- in \% on second line]\n";
      file << boldline << endl << Form("%27s", "");

      for (auto const& p : set_labels) {
        const char* label = p.second.data();  int len = strlen(label);
        file << Form("||%*s%*s",11+len/2,label,10-len/2,"");
      }

      for (auto const& sys : systematics) {
        if ( delUP[chan].find(sys) == delUP[chan].end() ) continue;

        file << endl << Form("%-27s", sys_labels[sys].data());
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", delUP[chan][sys][cutkey][set].first, delDN[chan][sys][cutkey][set].first);
        }
        file << endl << Form("%-27s", "");
        for (auto const& p : set_labels) {
          string set = p.first;
          file << Form("|| %9.2f %9.2f ", 100*delUP[chan][sys][cutkey][set].first/nom[chan][cutkey][set].first, 100*delDN[chan][sys][cutkey][set].first/nom[chan][cutkey][set].first);
        }
      } // end systematics loop 2 (Summary of systematics)

      file << endl << Form("%-27s", "Total Sys");
      for (auto const& p : set_labels) {
        string set = p.first;
        file << Form("|| %9.2f %9.2f ", sqrt( totalSys[chan][cutkey][set].first ), sqrt( totalSys[chan][cutkey][set].second ) );
      }
      file << endl << Form("%-27s", "");
      for (auto const& p : set_labels) {
        string set = p.first;
        file << Form("|| %9.2f %9.2f ", 100*sqrt(totalSys[chan][cutkey][set].first)/nom[chan][cutkey][set].first, 100*sqrt(totalSys[chan][cutkey][set].second)/nom[chan][cutkey][set].first);
      }

      /// LATEX TOTAL SYSTEMATICS ///
      TString cut_latex = cutname;  cut_latex.ReplaceAll(">=", "$\\geq$");

      file << endl << endl << "\\renewcommand{\\arraystretch}{2}\n";
      file << "\\begin{sidewaystable}\n";
      file << "\\resizebox{\\textheight}{!}{\n";
      file << "\\setlength\\tabcolsep{2pt}\n";
      file << "\\fontsize{3mm}{3mm} \\selectfont \n";
      file << "  \\begin{tabular}{ |c||c|c|c|c|c|c|c|c| }\n";
      file << " \\multicolumn{9}{c}{channel : $" + chan + "$} \\\\\n";
      file << "  \\multicolumn{9}{c}{" + cut_latex + "} \\\\\n";
      file << "  \\hline\n";
      file << " Systematic";
      for (auto const& p : set_latex) file << " & " + p.second;
      file << " \\\\\n  \\hline\\hline\n";

      for (auto const& sys : systematics) {
        if ( delUP[chan].find(sys) == delUP[chan].end() ) continue;

        file << " " + sys_latex[sys];
        for (auto const& p : set_latex) {
          string set = p.first;
          file << Form(" & %.2f %.2f", 100*delUP[chan][sys][cutkey][set].first/nom[chan][cutkey][set].first, 100*delDN[chan][sys][cutkey][set].first/nom[chan][cutkey][set].first);
        }
        file << " \\\\\n \\hline\n";
      } // end systematics loop 3 (Summary of systematics latex)

      file << " Total Systematics";
      for (auto const& p : set_labels) {
        string set = p.first;
        file << Form(" & %.2f -%.2f", 100*sqrt(totalSys[chan][cutkey][set].first)/nom[chan][cutkey][set].first, 100*sqrt(totalSys[chan][cutkey][set].second)/nom[chan][cutkey][set].first);
      }
      file << " \\\\\n \\hline\n";
      file << "  \\end{tabular}}\n";
      file << "\\end{sidewaystable}" << endl;
    }

    /// LATEX FINAL EVENT CUTFLOW ///
    file << endl << " \\newpage " << endl;
    file << "\\begin{sidewaystable}\n";
    file << "\\resizebox{\\textheight}{!}{\n";
    file << "\\setlength\\tabcolsep{2pt}\n";
    file << "\\fontsize{3mm}{3mm} \\selectfont \n";
    file << "  \\begin{tabular}{ |c||c|c|c|c|c|c|c|c|c|c| }\n";
    file << " \\multicolumn{11}{c}{channel : $" + chan + "$} \\\\\n";
    file << "  \\hline\n";
    file << " Cut";
    for (auto const& p : set_latex) {
      file << " & " + p.second;
      if (p.first == "bkg") file << " & Data & Data/Bkg";
    }
    file << " \\\\\n  \\hline\\hline\n";

    for (auto const& it_cut : nom[chan]) {
      string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);  TString cut_latex = cutname;  cut_latex.ReplaceAll(">=", "$\\geq$");

      double ratio = nom[chan][cut]["data"].first/nom[chan][cut]["bkg"].first;
      double bkgUP = sqrt( nom[chan][cut]["bkg"].second*nom[chan][cut]["bkg"].second + totalSys[chan][cut]["bkg"].first );
      double bkgDN = sqrt( nom[chan][cut]["bkg"].second*nom[chan][cut]["bkg"].second + totalSys[chan][cut]["bkg"].second );

      file << " " + cut_latex;
      for (auto const& p : set_latex) {
        string set = p.first;

        if (systematics.size() > 0) {
          file << Form(" & $%.1f^{+%.1f}_{-%.1f}$", nom[chan][cut][set].first,
                                                    sqrt( nom[chan][cut][set].second*nom[chan][cut][set].second + totalSys[chan][cut][set].first ),
                                                    sqrt( nom[chan][cut][set].second*nom[chan][cut][set].second + totalSys[chan][cut][set].second ) );
          if (set == "bkg")
            file << Form(" & $%.0f \\pm %.1f$ & $%.2f^{+%.2f}_{-%.2f}$", nom[chan][cut]["data"].first, nom[chan][cut]["data"].second,
                          ratio, ratio * sqrt( 1/nom[chan][cut]["data"].first + bkgUP*bkgUP/nom[chan][cut]["bkg"].first/nom[chan][cut]["bkg"].first ),
                                 ratio * sqrt( 1/nom[chan][cut]["data"].first + bkgDN*bkgDN/nom[chan][cut]["bkg"].first/nom[chan][cut]["bkg"].first ) );
        }
        else {
          file << Form(" & %.1f $\\pm$ %.1f", nom[chan][cut][set].first, nom[chan][cut][set].second);
          if (set == "bkg")
            file << Form(" & %.0f $\\pm$ %.1f & %.2f $\\pm$ %.2f", nom[chan][cut]["data"].first, nom[chan][cut]["data"].second,
                          ratio, ratio * sqrt( 1/nom[chan][cut]["data"].first + nom[chan][cut]["bkg"].second*nom[chan][cut]["bkg"].second
                                                                              / nom[chan][cut]["bkg"].first/nom[chan][cut]["bkg"].first ) );
        }
      }
      file << " \\\\\n \\hline\n";
      if (cutname == "MET Filters" || cutname == ">= 1 jet") file << " \\hline\n";
    } // end cut loop 5 (Final Summary Latex)

    file << "  \\end{tabular}}\n";
    file << "\\end{sidewaystable}" << endl;

    file.close();

    if (chan == "ll") continue;
    // Efficiency File //
    ofstream efile(dir + chan + "_efficiencies.txt");

    ifstream file2(dir + "/log_" + chan + ".txt");

    efile << file2.rdbuf();
    file2.close();

    efile << endl << boldline << endl;
    efile << "                                              Cut Flow Table: Summary\n";
    efile << boldline << endl << Form("%27s", "");
    for (auto const& p : set_labels) {
      if (p.first == "bkg") continue;
      const char* label = p.second.data();  int len = strlen(label);

      efile << Form("|||%*s%*s",9+len/2,label,9-len/2,"");
    }
    for (auto const& it_cut : nom[chan]) {
      string cut = it_cut.first, cutname = it_cut.first;  cutname.erase(0, 1);

      efile << endl << Form("%-27s", cutname.data());
      for (auto const& p : set_labels) {
        string set = p.first;
        if (set == "bkg") continue;

        efile << Form("|||     %1.6f     ", nom[chan][cut][set].first / nom[chan][cut_initial][set].first);
      }
      if (cutname == "MET Filters" || cutname == ">= 1 jet") efile << endl << singleline;
    } // end cut loop 6 (Efficiency Summary)

    TString cut_latex = cut_end;  cut_latex.Remove(0, 1);  cut_latex.ReplaceAll(">=", "$\\geq$");

    efile << "\n\n\\begin{center}\n";
    efile << "  \\begin{tabular}{ |c|c| }\n";
    efile << "  \\multicolumn{2}{c}{channel "<< chan << " } \\\\\n";
    efile << "  \\multicolumn{2}{c}{" + cut_latex + "} \\\\\n";
    efile << "  \\hline\n";
    efile << "  Process & Efficiency \\\\\n";
    efile << "  \\hline\\hline\n";

    for (auto const& p : set_latex) {
      string set = p.first;
      if (set == "bkg") continue;

      efile << Form("  %s & %1.6f \\\\\n  \\hline", p.second.data(),  nom[chan][cut_end][set].first / nom[chan][cut_initial][set].first) << endl;
    }
    efile << "  \\end{tabular}\n";
    efile << "\\end{center}" << endl;

    efile.close();
  } // end channel loop

  // Combination Table Event Yields and Signal Efficiencies //
  ofstream combofile(dir + "combination_tables.txt");
  combofile << "//////////////////////////" << endl;
  combofile << "////// EVENT YIELDS //////" << endl;
  combofile << "//////////////////////////" << endl << endl;
  combofile << "\\renewcommand{\\arraystretch}{2}\n";

  // use correct cross section for signals
  for (auto const& chan : channels) {
    for (auto const& it_cut : nom[chan]) {
      string cut = it_cut.first;

      nom[chan][cut][zprime].first      *= xs_zprime;  nom[chan][cut][zprime].second      *= xs_zprime;
      totalSys[chan][cut][zprime].first *= xs_zprime;  totalSys[chan][cut][zprime].second *= xs_zprime;
      nom[chan][cut][gluon].first       *= xs_gluon;   nom[chan][cut][gluon].second       *= xs_gluon;
      totalSys[chan][cut][gluon].first  *= xs_gluon;   totalSys[chan][cut][gluon].second  *= xs_gluon;
    }
  }

  // re-arrange order of set vectors (move signals to beginning)
  set_labels.insert( set_labels.begin(), *(--set_labels.end()) );  set_labels.pop_back();
  set_labels.insert( set_labels.begin(), *(--set_labels.end()) );  set_labels.pop_back();
  set_latex. insert( set_latex.begin(),  *(--set_latex.end()) );   set_latex. pop_back();
  set_latex. insert( set_latex.begin(),  *(--set_latex.end()) );   set_latex. pop_back();

  for (auto const& it_cut : nom["mm"]) {
    string cut = it_cut.first;

    if (cut < "Q= 1 Jet, = 0 btags, metR") continue;
    TString cutname = cut;  cutname.Remove(0, 1);  cutname.ReplaceAll(">=", "$\\geq$");

    combofile << "\\begin{tabular}{ |c|c|c|c|c| }\n";
    combofile << "\\multicolumn{5}{c}{" + cutname + "} \\\\\n";
    combofile << "\\hline\n";
    combofile << "Sample & $\\mu\\mu$ Channel & ee Channel & e$\\mu$ Channel & Combined \\\\\n";
    combofile << "\\hline\n";

    for (auto const& p : set_latex) {
      string set = p.first;
      combofile << p.second;
      for (auto const& chan : channels)
        combofile << Form(" & $%.1f^{+%.1f}_{-%.1f}$", nom[chan][cut][set].first, sqrt( nom[chan][cut][set].second*nom[chan][cut][set].second + totalSys[chan][cut][set].first ),
                                                                                  sqrt( nom[chan][cut][set].second*nom[chan][cut][set].second + totalSys[chan][cut][set].second ) );
      combofile << " \\\\\n" << endl;
      if (set == gluon || set == "wjet") combofile << "\\hline\n";
    }
    combofile << "Data";
    for (auto const& chan : channels)
      combofile << Form(" & $%.0f \\pm %.1f$", nom[chan][cut]["data"].first, nom[chan][cut]["data"].second);
    combofile << " \\\\\n" << endl;

    combofile << "Data/Bkg";
    for (auto const& chan : channels) {
      double ratio = nom[chan][cut]["data"].first / nom[chan][cut]["bkg"].first;
      double bkgUP = sqrt( nom[chan][cut]["bkg"].second*nom[chan][cut]["bkg"].second + totalSys[chan][cut]["bkg"].first );
      double bkgDN = sqrt( nom[chan][cut]["bkg"].second*nom[chan][cut]["bkg"].second + totalSys[chan][cut]["bkg"].second );

      combofile << Form(" & $%.2f^{+%.2f}_{-%.2f}$", ratio,
                         ratio * sqrt( 1/nom[chan][cut]["data"].first + bkgUP*bkgUP/nom[chan][cut]["bkg"].first/nom[chan][cut]["bkg"].first ),
                         ratio * sqrt( 1/nom[chan][cut]["data"].first + bkgDN*bkgDN/nom[chan][cut]["bkg"].first/nom[chan][cut]["bkg"].first ) );
    }
    combofile << " \\\\\n" << endl;
    combofile << "\\hline\n";

    combofile << "S(Z$'$)/Bkg";
    for (auto const& chan : channels)
      combofile << Form(" & %.5f", nom[chan][cut][zprime].first / nom[chan][cut]["bkg"].first);
    combofile << " \\\\\n" << endl;

    combofile << "S($\\textrm{g}_{\\textrm{kk}}$)/Bkg";
    for (auto const& chan : channels)
      combofile << Form(" & %.5f", nom[chan][cut][gluon].first / nom[chan][cut]["bkg"].first);
    combofile << " \\\\\n" << endl;
    combofile << "\\hline\n";

    combofile << "S(Z$'$)/$\\sqrt{\\textrm{S(Z$'$) + Bkg}}$";
    for (auto const& chan : channels)
      combofile << Form(" & %.3f", nom[chan][cut][zprime].first / sqrt( nom[chan][cut][zprime].first + nom[chan][cut]["bkg"].first ) );
    combofile << " \\\\\n" << endl;

    combofile << "S($\\textrm{g}_{\\textrm{kk}})/\\sqrt{\\textrm{S(g}_{\\textrm{kk}}\\textrm{) + Bkg}}$";
    for (auto const& chan : channels)
      combofile << Form(" & %.3f", nom[chan][cut][gluon].first / sqrt( nom[chan][cut][gluon].first + nom[chan][cut]["bkg"].first ) );
    combofile << " \\\\\n" << endl;
    combofile << "\\hline\n";

    combofile << "\\end{tabular}\n" << endl;
  } //end cut loop

  combofile << "/////////////////////////////////" << endl;
  combofile << "////// SIGNAL EFFICIENCIES //////" << endl;
  combofile << "/////////////////////////////////" << endl << endl;

  vector<string> v_gkk = { "RSGluonToTT_M-500", "RSGluonToTT_M-750", "RSGluonToTT_M-1000", "RSGluonToTT_M-1250", "RSGluonToTT_M-1500", "RSGluonToTT_M-2000",
                           "RSGluonToTT_M-2500", "RSGluonToTT_M-3000", "RSGluonToTT_M-3500", "RSGluonToTT_M-4000", "RSGluonToTT_M-4500", "RSGluonToTT_M-5000",
                           "RSGluonToTT_M-5500", "RSGluonToTT_M-6000", "RSGluonToTT_M-6500", "RSGluonToTT_M-7000", "RSGluonToTT_M-7500", "RSGluonToTT_M-8000" };

  vector<string> v_zp1 = { "ZprimeToTT_M500_W5_total", "ZprimeToTT_M750_W7p5_total", "ZprimeToTT_M1000_W10_total", "ZprimeToTT_M1250_W12p5_ext",
                           "ZprimeToTT_M1500_W15_total", "ZprimeToTT_M2000_W20_total", "ZprimeToTT_M2500_W25_total", "ZprimeToTT_M3000_W30_total",
                           "ZprimeToTT_M3500_W35_total", "ZprimeToTT_M4000_W40_total", "ZprimeToTT_M5000_W50_total", "ZprimeToTT_M6000_W60_total",
                           "ZprimeToTT_M7000_W70_total", "ZprimeToTT_M8000_W80_total" };

  vector<string> v_zp10 = { "ZprimeToTT_M500_W50_total", "ZprimeToTT_M750_W75_total", "ZprimeToTT_M1000_W100_ext", "ZprimeToTT_M1250_W125_total",
                            "ZprimeToTT_M1500_W150_total", "ZprimeToTT_M2000_W200_total", "ZprimeToTT_M2500_W250_total", "ZprimeToTT_M3000_W300_total",
                            "ZprimeToTT_M3500_W350_total", "ZprimeToTT_M4000_W400_total", "ZprimeToTT_M5000_W500_total", "ZprimeToTT_M6000_W600_total",
                            "ZprimeToTT_M7000_W700_total", "ZprimeToTT_M8000_W800_total", "ZprimeToTT_M9000_W900_total" };

  vector<string> v_zp30 = { "ZprimeToTT_M500_W150_total", "ZprimeToTT_M750_W225_total", "ZprimeToTT_M1000_W300_ext", "ZprimeToTT_M1250_W375_total",
                            "ZprimeToTT_M1500_W450_total", "ZprimeToTT_M2000_W600_total", "ZprimeToTT_M2500_W750", "ZprimeToTT_M3000_W900_total",
                            "ZprimeToTT_M3500_W1050_total", "ZprimeToTT_M4000_W1200_total", "ZprimeToTT_M5000_W1500_total", "ZprimeToTT_M6000_W1800_total",
                            "ZprimeToTT_M7000_W2100_total", "ZprimeToTT_M8000_W2400_total", "ZprimeToTT_M9000_W2700_total" };

  map<string, vector<string> > m_sigs = { {"$\\textrm{g}_{\\textrm{kk}}$", v_gkk}, {"Z$'$ (1\\%)", v_zp1}, {"Z$'$ (10\\%)", v_zp10}, {"Z$'$ (30\\%)", v_zp30} };

  combofile << "\\begin{tabular}{ |c|c|c|c|c|c| }\n";
  combofile << "\\hline\n";
  combofile << "\\multicolumn{2}{|c|}{Signal} & $\\mu\\mu$ Channel & ee Channel & e$\\mu$ Channel & Combined \\\\\n";
  combofile << "\\hline\n";

  for (auto const& i_sig : m_sigs) {
    combofile << Form("\\multirow{%i}{*}{%s} ", int(i_sig.second.size()), i_sig.first.data());
    for (auto const& sig : i_sig.second) {

      TString tsig = sig;
      int index = tsig.Contains("gluon") ? tsig.Length() : tsig.Last('_');
      TString mass = tsig(tsig.Index("M-")+2, index-tsig.Index("M-")-2);

      double mm_eff = nom["mm"][cut_end][sig].first / nom["mm"][cut_initial][sig].first * 100;
      double ee_eff = nom["ee"][cut_end][sig].first / nom["ee"][cut_initial][sig].first * 100;
      double em_eff = nom["em"][cut_end][sig].first / nom["em"][cut_initial][sig].first * 100;

      combofile << Form("& %s & %.4f & %.4f & %.4f & %.4f \\\\", mass.Data(), mm_eff, ee_eff, em_eff, mm_eff+ee_eff+em_eff) << endl;
    }
    combofile << "\\hline\n";
  }
  combofile << "\\end{tabular}\n" << endl;

  combofile.close();
}

//map(cut, map(dataset, (N, stat error) ) )
void readFile(const string& fileName, map<string, map<string, pair<double, double> > >& cuts) {

  ifstream file(fileName);
  if ( !file.is_open() ) cout << fileName + " not found!" << endl;
  //cout << " =====> Reading file:  " << fileName << endl ;

  string line, channel;
  TString dataset;
  double weight=-1;
  char cut_order = 'A';  //cheesy way to keep cut order in map
  while (getline(file, line)) {

    if (line.length() > 0) { while (line.at(0) == ' ') line.erase(0, 1); }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;
    string str = line.substr(0, delim_pos);

    if (str == "Weight") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      weight = stod(line);
    }
    else if (str == "Processing") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      TString l = line;
      dataset = l( l.Last('/')+1, l.Last('.')-l.Last('/')-1 );

      //cout << "Filling " << dataset << " using weight " << weight << endl;
    }
    else if (str == "Channel:") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      channel = line;
    }
    else if (str == "Cut") {
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      //dataset = line.data();
      //cout << "Filling " << dataset << " using weight " << weight << endl;
      cut_order = 'A';  //reset cut order
    }
    else {
      delim_pos = line.find("|||");
      if (delim_pos == -1) continue;

      str = line.substr(0, delim_pos);
      if (str == "") continue;

      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());
      string cut = str;
      cut = cut_order + cut;  cut_order++;

      line.erase(0, delim_pos + 3);
      while (line.at(0) == ' ') line.erase(0, 1);

      str = line.substr(0, line.find("|||"));
      while (str.at(str.length()-1) == ' ') str.erase(str.length()-1, str.length());

      double N = stod(str);

      //data
      if ( dataset.Contains("Muon", TString::kIgnoreCase) || dataset.Contains("DoubleEG", TString::kIgnoreCase) )
        cuts[cut]["data"] = make_pair( N, sqrt(N) );

      //signal or background
      else {
        string key;

        if      ( dataset.Contains("ttbar", TString::kIgnoreCase) ||
                  dataset.Contains("TTsemilep", TString::kIgnoreCase) )                    key = "ttbar";
        else if ( dataset.Contains("dy", TString::kIgnoreCase) )                           key = "dy";
        else if ( dataset.Contains("wjet", TString::kIgnoreCase) )                         key = "wjet";
        else if ( dataset.Contains("st_", TString::kIgnoreCase) )                          key = "st";
        else if ( dataset.Contains("ww", TString::kIgnoreCase) ||
                  dataset.Contains("wz", TString::kIgnoreCase) ||
                  dataset.Contains("zz", TString::kIgnoreCase) )                           key = "vv";
        else if ( dataset.Contains("qcd", TString::kIgnoreCase) &&
                  ( (channel == "mm" &&  dataset.Contains("Mu", TString::kIgnoreCase)) ||
                    (channel == "ee" && !dataset.Contains("Mu", TString::kIgnoreCase)) ||
                    (channel == "em" &&  dataset.Contains("Mu", TString::kIgnoreCase)) ) ) key = "qcd";
        else                                                                               key = dataset.Data(); //signal

        if ( cuts[cut].find(key) == cuts[cut].end() ) cuts[cut][key] = make_pair( N, weight * sqrt(N/weight) );
        else {
          cuts[cut][key].first += N;
          cuts[cut][key].second = sqrt( cuts[cut][key].second*cuts[cut][key].second + weight*N );
        }

        //total background
        if ( !dataset.Contains("zprime", TString::kIgnoreCase) && !dataset.Contains("gluon", TString::kIgnoreCase) ) {
          if ( dataset.Contains("qcd", TString::kIgnoreCase) && !( (channel == "mm" &&  dataset.Contains("Mu", TString::kIgnoreCase)) ||
                                                                   (channel == "ee" && !dataset.Contains("Mu", TString::kIgnoreCase)) ||
                                                                   (channel == "em" &&  dataset.Contains("Mu", TString::kIgnoreCase)) ) ) continue;

          key = "bkg";
          if ( cuts[cut].find(key) == cuts[cut].end() ) cuts[cut][key] = make_pair( N, weight * sqrt(N/weight) );
          else {
            cuts[cut][key].first += N;
            cuts[cut][key].second = sqrt( cuts[cut][key].second*cuts[cut][key].second + weight*N );
          }
        }
      }
    }
  }
  file.close();
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
