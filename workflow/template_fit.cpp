

#include "tTree.h"
#include "binning_histos_small.h"
#include "result_paths.h"   // obsProdTag(): which MC production carries which observable
#include "Help_Functions.h"
#include "Draw_EEC.h"
#include "TMarker.h"   // best (A, f) point on the chi2 map, do_template_fit_chi2scan()
//#include "../CMSStyle.C" // CMS style impored from Matthew


// NOT CALLED since 2026-10-27: template_fit() runs do_template_fit_chi2scan() (AN-style chi2 scan,
// 2B vs 1B, no 0B) instead. Kept intact so the RooFit fit with 0B can be switched back on.
void do_template_fit_combined(const TString &HighEGdata_name, const TString &LowEGdata_name, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG, bool& also_bjet,  Variation ivar = NOMINAL){


    // Test CMS style 
    // setCMSStyle();

    // All plots go into one flat directory (sDirname_www); root files stay in sDirname.
    // Only the nominal variation is plotted; the others are still fitted because
    // draw_variation_uncertainity() needs their fractions for the 0B systematic.
    const bool save_png = (ivar == NOMINAL);
    TString sDir_canvas = sDirname_www;
        if (save_png) gSystem->mkdir(sDir_canvas, kTRUE);


    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    // The names carry the observable suffix from gFitObs(): "" for dR, so the dR fit reads
    // exactly the histograms it always read, and "_B" for the momentum balance.
    const ObsDef &obs = gFitObs();
    TString namehData = obs.n("h3D_data"); // h3D_data with eec on the weight 
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = obs.n("h3D_bb");
    TString nameh1B = obs.n("h3D_b");

    TString nameh0B = obs.n("h3D_0b");
    // -- Additional histograms
    TString namehmore2B = "";
 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(templates, "read");
            if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
        TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
        TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
            if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
            if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
        TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone("h3D_0b");  // clone to the
            // CANONICAL name, not nameh0B: the input histogram carries the observable
            // suffix ("h3D_0b_B"), but everything downstream -- Draw_EEC.h in particular
            // -- looks the templates up in the fit output by their unsuffixed names, the
            // same way h3D_b and h3D_bb are cloned just above. For dR the two spellings
            // coincide, which is why only the balance fit tripped over it.
            if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}

    //-- Bjet sample:
      TH3D *h3D_b_bjet = nullptr;
      TH3D *h3D_bb_bjet  = nullptr;
    if (also_bjet){
        TFile *file_bjet = new TFile(templates_bjet, "read");
        if (!file_bjet || file_bjet->IsZombie() ) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        
        TH3D* tmp_b  = dynamic_cast<TH3D*>(file_bjet->Get(nameh1B));
        TH3D* tmp_bb = dynamic_cast<TH3D*>(file_bjet->Get(nameh2B));
        if (!tmp_b) {Warning("Get", "Histogram %s not found", nameh1B.Data());}
        else{
            h3D_b_bjet = static_cast<TH3D*>(tmp_b->Clone("h3D_b_bjet"));
            h3D_b_bjet->SetDirectory(nullptr);
        }
        if (!tmp_bb) {Warning("Get", "Histogram %s not found", nameh2B.Data());}
        else {
            h3D_bb_bjet = static_cast<TH3D*>(tmp_bb->Clone("h3D_bb_bjet"));
            h3D_bb_bjet->SetDirectory(nullptr);
        }

        file_bjet->Close(); // NEW
    }

    //Open dataset:
    TFile *file_data = new TFile(HighEGdata_name, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
            cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
            h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
            if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
            if(alsoLowEG){
                TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
                    if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
                    cout << "file data name " << file_data_LowEG->GetName() << endl;
                TH3D *h3D_dataLowEG;
                    h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
                    if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
                h3D_data->Add(h3D_dataLowEG);
            }
        // -- data: set style
        h3D_data->SetTitle("Data");
        styleData(h3D_data);

        // -- bjet histograms
        if(also_bjet){
            h3D_b_bjet->SetTitle("1B (bjet)");h3D_bb_bjet->SetTitle("2B (bjet)");
            h3D_b_bjet->GetXaxis()->SetTitle("m_{2B} [GeV]"); // Name is wrong in MC templates (EEC)
            // Same category colours as the qcd templates -- 1B is 1B whichever sample it
            // came from -- with a hatched fill as the only thing saying "bjet sample".
            styleTemplate(h3D_b_bjet,  TFColor::c1B()); h3D_b_bjet ->SetFillStyle(3354);
            styleTemplate(h3D_bb_bjet, TFColor::c2B()); h3D_bb_bjet->SetFillStyle(3345);
        }

        // -- qcd
        h3D_b->SetTitle("1B (qcd)");h3D_bb->SetTitle("2B (qcd)"); h3D_nob ->SetTitle("0B (qcd)");
        styleTemplate(h3D_b,   TFColor::c1B());
        styleTemplate(h3D_bb,  TFColor::c2B());
        styleTemplate(h3D_nob, TFColor::c0B());

    // -- Define used bins. The binning is owned by binning_histos_small.h and the loops
    // below index jtpt_binsVector / dr_binsVector, so verify the inputs match it rather
    // than adopting whatever they contain.
    if (!CheckInputBinning(h3D_data)) return;
    if (!CheckInputBinning(h3D_b))    return;
    if (!CheckInputBinning(h3D_bb))   return;
    if (!CheckInputBinning(h3D_nob))  return;
    if (also_bjet) {
        if (!CheckInputBinning(h3D_b_bjet))  return;
        if (!CheckInputBinning(h3D_bb_bjet)) return;
    }

    bins_pt = jtpt_binsVectorSize - 1; // == jtpt_bins, and == h3D_data->GetNbinsZ()
    bins_dr = obs.nbins;   // the MEASURED OBSERVABLE's bin count, not necessarily dR's
    mb_bins = mb_binsVectorSize   - 1; // updated below if the mass axis is rebinned
        cout << "-- Data hist initial binning" << endl;
        cout << "pt bins = "<< bins_pt << endl;
        cout << "dr bins = "<< bins_dr << endl;
        cout << "mb bins = "<< mb_bins << endl;

    // -- Write used inputs 
    fout->cd();
    h3D_data->Write();
    if (also_bjet) { h3D_b_bjet->Write(); h3D_bb_bjet->Write();}
    h3D_b->Write(); h3D_bb->Write(); h3D_nob->Write();

    // ---- Weighted qcd + bjet combination (2026-10-24, Zoe) ----------------------------------
    // The b-enriched sample's event weights are normalised like an inclusive sample: its 2B
    // yield is ~14.5x qcd's, for the same physical b-jet cross section. So, before anything is
    // summed: (1) scale bjet to qcd's 2B yield, k = S_qcd / S_bjet; (2) mix the two with
    // statistics weights, w_bjet = Neff_bjet / (Neff_qcd + Neff_bjet), w_qcd = 1 - w_bjet. The
    // 1B and 2B templates become w_qcd*qcd + w_bjet*k*bjet -- at qcd normalisation, so the 0B
    // template (qcd only: the b-enriched light jets are a filter-biased subset) stays
    // comparable -- and the fractions below are taken from this combination too.
    // Same weights as setMixWeights() in apply_unfolding_2d.C (same 2B templates).
    // 2026-10-26: 1B gets its OWN k1 = S1_qcd / S1_bjet. One k from 2B was fine for Pythia
    // (bjet/qcd 14.41 for 2B, 14.43 for 1B, flat in pT, dR and A) but not for Herwig (11.9 vs
    // 12.7): there the combined 1B sat 6% above qcd normalisation, against a qcd-only 0B, which
    // shifted the 1B:0B mix of the background. The mixing weights stay the 2B ones.
    // Done AFTER the inputs are written (moved 2026-10-26), so h3D_b / h3D_bb in the fit file are
    // the raw qcd templates their names say -- Draw_EEC.h reads them as such.
    if (also_bjet && h3D_bb_bjet && h3D_b_bjet) {
        auto sums = [](const TH3D *h, double &S, double &Q) {
            S = 0.; Q = 0.;
            // Cells are 0 .. GetNcells()-1; "<=" read the last one twice (ROOT clamps the index).
            // for (int i = 0; i <= h->GetNcells(); ++i) { S += h->GetBinContent(i); Q += std::pow(h->GetBinError(i), 2); }   // until 2026-10-26
            for (int i = 0; i < h->GetNcells(); ++i) { S += h->GetBinContent(i); Q += std::pow(h->GetBinError(i), 2); }
        };
        double Sq, Qq, Sb, Qb, S1q, Q1q, S1b, Q1b;
        sums(h3D_bb, Sq, Qq); sums(h3D_bb_bjet, Sb, Qb);
        sums(h3D_b, S1q, Q1q); sums(h3D_b_bjet, S1b, Q1b);
        if (Sq > 0. && Sb > 0. && Qq > 0. && Qb > 0. && S1q > 0. && S1b > 0.) {
            const double k = Sq / Sb, k1 = S1q / S1b, nq = Sq * Sq / Qq, nb = Sb * Sb / Qb;
            const double wb = nb / (nq + nb), wq = 1. - wb;
            printf("qcd+bjet weighted combination: k = %.5f (bjet 2B / qcd 2B = %.2f), k1 = %.5f (bjet 1B / qcd 1B = %.2f),"
                   " Neff qcd %.0f, bjet %.0f -> w_qcd %.4f, w_bjet %.4f (bjet 2B scaled by %.5f, 1B by %.5f)\n",
                   k, 1. / k, k1, 1. / k1, nq, nb, wq, wb, wb * k, wb * k1);
            h3D_bb->Scale(wq);  h3D_b->Scale(wq);
            h3D_bb_bjet->Scale(wb * k); h3D_b_bjet->Scale(wb * k1);
            // h3D_bb_bjet->Scale(wb * k); h3D_b_bjet->Scale(wb * k);   // before 2026-10-26: 1B with the 2B k
        } else {
            Error("weighted combination", "empty 2B or 1B template -- qcd %g / %g, bjet %g / %g; not combining", Sq, S1q, Sb, S1b);
            return;
        }
    }

    // -- Choose the observable binning (dR, B) -- one source, gFitObs()
    const double* yBins = nullptr;
    int N_bins_dr = 0; N_bins_dr = obs.nbins; yBins = obs.bins;


    // -- Rebinning in mass axis when needed 
    if(ivar == FITRANGE_0_7)
    {
        // TH3D* h3_fitRange0to7 = MergeLastMassBinTo7GeV(h3_original);
        h3D_data = MergeLastMassBinTo7GeV(h3D_data);
                styleData(h3D_data);   // the merge hands back a fresh, unstyled histogram

        h3D_b    = MergeLastMassBinTo7GeV(h3D_b);
        h3D_bb   = MergeLastMassBinTo7GeV(h3D_bb);
        h3D_nob  = MergeLastMassBinTo7GeV(h3D_nob);

        if(also_bjet)
            {
                h3D_b_bjet  = MergeLastMassBinTo7GeV(h3D_b_bjet);
                h3D_bb_bjet = MergeLastMassBinTo7GeV(h3D_bb_bjet);
            }

        // update mass bins 
        mb_bins = h3D_data->GetNbinsX();
    }
    else if(ivar == FITRANGE_0_8)
    {
        h3D_data = MergeLastMassBinTo8GeV(h3D_data);
                styleData(h3D_data);   // the merge hands back a fresh, unstyled histogram

        h3D_b    = MergeLastMassBinTo8GeV(h3D_b);
        h3D_bb   = MergeLastMassBinTo8GeV(h3D_bb);
        h3D_nob  = MergeLastMassBinTo8GeV(h3D_nob);

        if(also_bjet)
            {
                h3D_b_bjet  = MergeLastMassBinTo8GeV(h3D_b_bjet);
                h3D_bb_bjet = MergeLastMassBinTo8GeV(h3D_bb_bjet);
            }

        // update mass bins 
        mb_bins = h3D_data->GetNbinsX();
    }
   

    // ---------------------------------------------------------
    /// -- For later drawing of S/B fractions, store true and fit result in the following histograms
        // Note that: since inetgaretd bins are pt 0 and dr 0, the hist of fractions should have #bins + 1 size 
        // x = dr, y = jetpt
        // Axis here is for the bin number instead of the values 
        TH2D *h_sig_fraction = new TH2D("h_sig_fraction", Form(";%s; jet pt", obs.axis.Data()), N_bins_dr+ 1,  1, N_bins_dr+ 2,
                                                                           h3D_data->GetNbinsZ() + 1 , 1, h3D_data->GetNbinsZ()+ 2 );                                                                                                                              
            h_sig_fraction->Reset();
            TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
            TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
            TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error"); 
            TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
            TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
            TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
            TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");

        // -- For unfolding: without integarted bins 
        TH2D *h_sig_fraction_fit = new TH2D("h_sig_fraction_fit", Form(";%s; jet pt", obs.axis.Data()), N_bins_dr,  yBins, jtpt_bins, jtpt_binsVector);                                                                                                                              
            h_sig_fraction_fit->Reset();
            TH2D *h_bkg_fraction_fit = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_fraction_fit");
            TH2D *h_sig_fraction_fit_error = (TH2D *) h_sig_fraction_fit->Clone("h_sig_fraction_fit_error");
            TH2D *h_bkg_fraction_fit_error = (TH2D *) h_bkg_fraction_fit->Clone("h_bkg_fraction_fit_error"); 
            // and corresponding true fractions for drawing purpose 
            TH2D *h_sig_frac_true_fitbins = (TH2D *) h_sig_fraction_fit ->Clone("h_sig_frac_true_fitbins");
            TH2D *h_sig_frac_true_error_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_sig_frac_true_error_fitbins");
            TH2D *h_bkg_frac_true_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_frac_true_fitbins");
            TH2D *h_bkg_frac_true_error_fitbins = (TH2D *) h_sig_fraction_fit->Clone("h_bkg_frac_true_error_fitbins");


    //-----------------------
    // -- CMS plot aesthetics, applied once for every canvas this function makes.
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);        // no ROOT title box: the CMS header + the dr/pT text carry it
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1);        // ticks mirrored on all four sides
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);

    // --- Vector to test the convergence
    std::vector <std::pair<int, int>> non_converge_bins;
 

    // Fitting - loop over dr and jtpt entries
    // Bin0: is integarted over the range 
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
    // for(Int_t ibin_pt = 1 ; ibin_pt <= 2; ibin_pt++){ // test 
        for(Int_t ibin_dr = 1; ibin_dr <= N_bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 0; ibin_dr <= 1; ibin_dr++){ // test 
            
            // define slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;
            Int_t SliceFirstbin_pt = ibin_pt;
            Int_t SliceLastbin_pt =  ibin_pt;
            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = N_bins_dr;}
            if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

            // Make projections 
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_data_mb->GetXaxis()->SetTitle("m_{2B} [GeV]");
                h_data_mb->SetTitle(h3D_data->GetTitle()); // upadte projection title
                styleData(h_data_mb);
            
            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
                h_bb->SetTitle( h3D_bb->GetTitle() );
                h_b->SetTitle( h3D_b->GetTitle() );
                h_nob->SetTitle( h3D_nob->GetTitle() );

            //  Make slices for bjet
            TH1D *h_bb_bjet = nullptr;
            TH1D *h_b_bjet = nullptr;

            if(also_bjet){
                h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr,  SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet->SetTitle(h3D_b_bjet->GetTitle());
                h_bb_bjet->SetTitle(h3D_bb_bjet->GetTitle()); 
            }
           
      
            // -- Deattach hists from root file 
            h_data_mb->SetDirectory(0);
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_nob->SetDirectory(0);
            if(also_bjet){
                h_b_bjet->SetDirectory(0);
                h_bb_bjet->SetDirectory(0);                
            }

            //Define common (pt, dr) canvas name  
            // Plot filenames name the observable. dR keeps the literal "deltaR" it has always
            // had, so its PNG names are unchanged; the balance gets "Bbin".
            TString sname_canvas = Form("ptbin_%d_%sbin_%d", ibin_pt,
                                        (obs.name == "dr") ? "deltaR" : obs.name.Data(), ibin_dr);

            // --- Compute Integrals
            // data:
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
            // --  Calculate true fractions to be used as initial values for the fit
            // double int2 = h_bb->Integral(1, mb_bins, "width");
            // double int1 = h_b ->Integral(1, mb_bins, "width");
            double int2 = h_bb->Integral(1, mb_bins, "width") + (also_bjet ? h_bb_bjet->Integral(1, mb_bins, "width") : 0.);
            double int1 = h_b ->Integral(1, mb_bins, "width") + (also_bjet ? h_b_bjet ->Integral(1, mb_bins, "width") : 0.);
            double int0 = h_nob->Integral(1, mb_bins, "width"); 
            double tot = int0 + int1 + int2;
                // compute the true fractions and their errors (for later comparison)           
                double sig_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : (int2 / (int0 + int1 + int2)); 
                double bkg_fraction_b_true = (int0  + int1 + int2 ) == 0 ? 0 : (int1 / (int0 + int1 + int2));
                double bkg_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : ( (int0 + int1)/(int0 + int1 + int2) );

                // 2026-10-24: the weighted qcd+bjet combination, like int2/int1 above (were qcd only)
                TH1D* h_b_comb  = (TH1D*) h_b ->Clone("h_b_comb");  if (also_bjet) h_b_comb ->Add(h_b_bjet);
                TH1D* h_bb_comb = (TH1D*) h_bb->Clone("h_bb_comb"); if (also_bjet) h_bb_comb->Add(h_bb_bjet);
                double True_bkg_b_err = 0.0;
                double True_bkg_b_integral = h_b_comb->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");

                double True_sig_err = 0.0;
                double True_sig_integral = h_bb_comb->IntegralAndError(1, mb_bins, True_sig_err, "width");

                // total bkg fraction: 0B + 1B (1B combined, 0B qcd)
                TH1D* h_b_nob_dijet = (TH1D*) h_b_comb->Clone("h_b_nob_dijet");  h_b_nob_dijet->Add(h_nob);
                double True_bkg_err = 0.0;
                double True_bkg_integral = h_b_nob_dijet->IntegralAndError(1, mb_bins, True_bkg_err, "width");

            // -- Compute other useful integrals 
                // From bjet sample 
                double int2_bjet = 0; 
                double int1_bjet =0;  
                if (also_bjet){
                    int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");
                    int2_bjet = h_bb_bjet->Integral(1, mb_bins, "width"); 
                }
      
                // -- check integrals in bjet and qcd sample
                // relative fraction of 1B: 2B in bjet and qcd samples is SAME. So you can combine the two samples with simple + (without reweighting).
                // cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int2 h_bb=" << int2 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "Dijet: int0 of 0B" << int0 << endl;

                // cout << "Dijet: 2B/(1B + 2B) = " << int2/(int2 + int1) << endl;
                // std::cout << "int2 h_bb_bjet=" << int2_bjet << std::endl;
                // std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                // cout << "Bjet: 2B/(2B + 1B) = " << int2_bjet/(int2_bjet + int1_bjet) << endl;
                // cout << " Dijet: 0B/(1B+2B) = " << int0/(int2 + int1) << endl;

            // --- Prepare PDFs for template fit
            // 1- Combine qcd + bjet samples 
            // Signal: 2B   
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){h_sumsig->Add(h_bb_bjet); h_sumsig->SetTitle("2B: qcd+bjet");}
                else{h_sumsig->SetTitle("2B: qcd");}
                Int_t h_sig_bins = h_sumsig->GetNbinsX();
                TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

            // Bkg: first 1B, and later added effectively 0B
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){h_sumbkg->Add(h_b_bjet); h_sumbkg->SetTitle("1B: qcd+bjet");}
                else{h_sumbkg->SetTitle("1B: qcd");}
                Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
                TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

                // for drawings
                    // Total sum bkg = 0B + 1B  
                    TH1D* h_sumbkg_0b_1b = (TH1D*) h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
                    h_sumbkg_0b_1b->Add(h_nob);
                    h_sumbkg_0b_1b->SetTitle("0B +1B (qcd+bjet)");


                // -- set the combined samples styles before fit
                styleTemplate(h_sumbkg, TFColor::c1B());   // 1B (qcd [+ bjet])
                styleTemplate(h_sumsig, TFColor::c2B());   // 2B (qcd [+ bjet])
                    // and for the bkg 1B + 0B: a sum of two categories, so the neutral
                    // colour rather than blue or green -- it is neither on its own.
                    styleTemplate(h_sumbkg_0b_1b, TFColor::bkg());
                    // 0B only
                    styleTemplate(h_nob, TFColor::c0B());

                // write to rootfile the used slices 
                fout->cd();
                h_data_mb->Write(); // data 
                h_sumsig->Write(); // 2B (qcd + bjet)
                h_sumbkg->Write(); // 1B (qcd + bjet)
                h_sumbkg_0b_1b->Write(); // 1B (qcd + bjet) + 0B 


/* ---- disabled (kept for reference): prefit control canvases (all contributions / raw templates) ----
    // -------- Draw prefits 
        // contirbutions seperated  
                // absolute yields: seperated contributions: qcd, bjet 
            THStack hstack_all_beforefit (Form("hstack_all_beforefit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms before fit");
                hstack_all_beforefit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                if(also_bjet) hstack_all_beforefit.Add(h_bb_bjet);
                hstack_all_beforefit.Add(h_bb);
                if(also_bjet) hstack_all_beforefit.Add(h_b_bjet); 
                hstack_all_beforefit.Add(h_b);
                hstack_all_beforefit.Add(h_nob);
                if (hstack_all_beforefit.GetMaximum() > h_data_mb->GetMaximum()/1e+04) { hstack_all_beforefit.SetMaximum(1.3* hstack_all_beforefit.GetMaximum());}
                else { hstack_all_beforefit.SetMaximum(1.3* h_data_mb->GetMaximum()/1e+04);}

                auto canva_beforefit = new TCanvas(Form("All_contributions_beforefit_%d_%d", ibin_dr, ibin_pt),"", 800, 800 );
                        canva_beforefit->cd();
                        // canva_beforefit->SetLogy();
                        hstack_all_beforefit.Draw("hist E");
                        TH1D* hd_norm_slice = (TH1D*) h_data_mb->Clone("hd_norm_slice");
                            hd_norm_slice->Scale(1./10000);
                            hd_norm_slice->Draw("HIST E same"); 
                            hd_norm_slice->SetTitle("data/1e+04");
                            gPad->Modified();   
                            gPad->Update();
                            canva_beforefit->Modified();
                            canva_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_beforefit->Write();
                            canva_beforefit->Print(Form("%s/%s_allcontributions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Draw Stack for combined qcd + bjet: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before normalization
            // absolute yields of Signal and bkg  
            THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_sumsig);
                    hstack_templatesforfit.Add(h_sumbkg_0b_1b);    
                    hstack_templatesforfit.SetMaximum(1.2 * hd_norm_slice->GetMaximum());
                
                    auto canva_sum_beforefit = new TCanvas(Form("templates_beforefit_%d_%d", ibin_dr, ibin_pt), "", 800, 800 );
                        canva_sum_beforefit->cd();
                        hstack_templatesforfit.Draw("hist E");
                        hd_norm_slice->Draw("HIST E same"); 
                            gPad->Modified();   
                            gPad->Update();
                            canva_sum_beforefit->Modified();
                            canva_sum_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_sum_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_sum_beforefit->Write();
                            canva_sum_beforefit->Print(Form("%s/%s_templates_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/ // ---- end disabled block ----

        // -- Safety for empty bins
            // Scale-free: "empty" must mean NO ENTRIES, not "integral below 1". The old
            // test was an absolute < 1 on the summed template integrals, which is only ever
            // a no-op for the EEC-weighted fit, where (pt1*pt2)^n makes those integrals
            // ~1e7. Run the same fit on UNWEIGHTED templates (EEC_WEIGHT_OFF, the yield
            // measurement) and the integrals are order 1 -- a whole QCD block sums to 0.86 --
            // so the sparser high-dr bins tripped this and were silently skipped, leaving
            // their signal fraction at exactly 0. Comparing against 0 asks the question the
            // comment always claimed to ask, and is a no-op for the EEC fit.
            if ((int0 + int1 + int2) <= 0. ) {
                cout << " ----------- empty bin -----------  ";
                // No MC but data: the fraction would stay 0 and empty this data bin (2026-10-26).
                if (integral_inputdata > 0.) non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
                // ----- To avoid empty bins if exist!, set them to eps value
                const double eps = 1e-6; 
                for (int i = 1; i <= h_sig_bins; i++){
                    if (h_sig->GetBinContent(i) <= 0){
                        h_sig->SetBinContent(i, eps);
                        cout << 
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
                        << endl; 
                    }
                }

                for (int i = 1; i <= h_bkg_bins; i++){
                    if (h_bkg->GetBinContent(i) <= 0){
                        h_bkg->SetBinContent(i, eps);
                        cout <<
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
                        << endl; 
                    }
                }
                /// Normalize safely for the Roofit
                if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
                {
                    cerr<< "Signal or bkg templates has zero integral"<< endl;
                    non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr)); // no fraction (2026-10-26)
                    continue;
                }
                // normalize signal and bkg
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                    // cout << "After Normalization (qcd +bjet) True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                    // cout << "After Normalization (qcd+bjet)  True 1B integral =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;

                // -- Effective bkg PDF. With a + b + c = 1 (a: 2B, b: 1B, c: 0B) and
                // b' + c' = 1, this gives b' = b/(b+c) and c' = c/(b+c), with b+c = 1-a.
                // NOTE: a, b, c come from int2/int1/int0. Since 2026-10-24 int2/int1 are the
                // WEIGHTED qcd+bjet combination (bjet normalised to qcd's 2B yield, see the top;
                // since 2026-10-26 its 1B to qcd's 1B yield),
                // int0 qcd 0B. Until then they were qcd alone, because the raw bjet weights
                // (~14.5x) would have biased c' downwards; normalised, they no longer do.
                    double eff_bkg0B = (1 - sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
                    double eff_bkg1B = 1. - eff_bkg0B;// b`
            // Build effective bkg hist: with new relaitve normalization, the integral should = 1
            // Normlaize h_nob to be added effectively to the normalized bkg distribution
                TH1D* norm_h_nob = (TH1D*) h_nob->Clone("norm_h_nob");
                // A slice without any 0B jet: 1/0 made this template NaN, and 0 x NaN made the
                // whole background NaN even at eff_bkg0B = 0 (2026-10-26). The background is then
                // 1B alone -- eff_bkg0B is 0 there anyway, int0 being 0.
                const double int_nob = norm_h_nob->Integral(1, h_bkg_bins, "width");
                if (int_nob > 0.) norm_h_nob->Scale(1./int_nob);
                else cout << Form("WARNING: (ptbin %d, deltaRbin %d) has no 0B jet -- background = 1B only", ibin_pt, ibin_dr) << endl;
                // 1B on its own, for the post-fit 1B slice (h_bkg is the 1B+0B mixture by then)
                TH1D* norm_h_1b = (TH1D*) h_sumbkg->Clone("norm_h_1b");
                const double int_1b = norm_h_1b->Integral(1, h_bkg_bins, "width");
                if (int_1b > 0.) norm_h_1b->Scale(1./int_1b);
            // Variation of 0B template contribution == the LIGHT-JET MISTAG systematic.
            // The 0B template is the jets with jtNbHad == 0 that still passed the b tag,
            // i.e. mistagged light and charm jets. How many of them sit under the data is
            // taken from MC alone -- nothing in this fit constrains it, because 0B is not a
            // free component: it is folded into the background PDF at the MC-predicted
            // ratio c' = 0B/(0B+1B). So the mistag rate enters the measurement only through
            // the SHAPE of the effective background, and the systematic is to rebuild that
            // shape with the 0B admixture scaled by w_var_0B.
            //
            // w = 2 / w = 0 is a deliberately conservative +/-100% on the mistag rate: it
            // brackets the measured light-flavour mistag scale factors (typically 20-50% at
            // a tight working point) by a wide margin, and it needs no external SF input.
            double w_var_0B = 1; 
                if (ivar == VARIED0B_UP) { w_var_0B = 2.0;}
                else if (ivar == VARIED0B_DOWN) {w_var_0B = 0.0;}
                else if (ivar == NOMINAL){w_var_0B = 1;}
                cout << "Before variation: eff 1B = "<< eff_bkg1B << endl;

                // change eff 1B accordingly                 
                eff_bkg0B *= w_var_0B;
                // Guard: eff_bkg1B = 1 - eff_bkg0B, so any eff_bkg0B > 1 makes the 1B term
                // in the Add() below NEGATIVE. The integral check further down would still
                // pass -- the two weights sum to 1 by construction -- but RooHistPdf would
                // silently clip the negative bins and fit a PDF that is not the one built
                // here. It cannot happen at the current WP (c' <= 0.19, so 2c' <= 0.37), but
                // a looser b tag or a coarser dR bin can reach it, so refuse loudly instead.
                if (eff_bkg0B > 1. || eff_bkg0B < 0.) {
                    cout << Form("WARNING: (ptbin %d, deltaRbin %d) 0B variation w = %.1f gives "
                                 "eff 0B = %.4f, outside [0,1] -- clipped. The variation is "
                                 "truncated in this bin, so its systematic is a LOWER bound.",
                                 ibin_pt, ibin_dr, w_var_0B, eff_bkg0B) << endl;
                    eff_bkg0B = std::min(1., std::max(0., eff_bkg0B));
                }
                eff_bkg1B = 1. - eff_bkg0B;
                    
                    cout << "eff 0B weight = "<< w_var_0B << endl;
                    cout << "eff 0B = "<< eff_bkg0B << endl;
                    cout << "eff 1B = "<< eff_bkg1B << endl;

                // --- Build effective Bkg template 
                h_bkg->Add(h_bkg, norm_h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist            
                // update combined normalized bkg integral and error  
                double err_int;
                double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
                    // cout << "Effective Bkg normalized hist integral +/- uncertainity = ? (should be 1):  "<< int_val << "+/-" << err_int << endl;
                    // cout << "effective 1B and 0B fractions = (sum should be 1): "<< eff_bkg1B << ", " << eff_bkg0B << endl;
                    // avoid float point failur 
                    // Written as !(<= tol) so a NaN fails too: NaN > tol is false (2026-10-26).
                    // Returning here writes no h_sig_fraction_fit, so the unfolding stops on it.
                    if (!(std::abs(int_val - 1) <= 1e-06)) {cerr << "ERROR: Effective BKG PDF is not normalized to 1 (integral " << int_val << ")!"<< endl; return;}
                    if (!(std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1) <= 1e-06)) {cerr << "ERROR: Signal PDF is not normalized to 1!"<< endl; return;}

            // -- Draw normalized PDFs of signal and effective background as inputs: prefit
            // set PDF titles 
            h_bkg->SetTitle("1B + 0B");
            h_sig->SetTitle("2B");
            styleTemplate(h_bkg, TFColor::bkg());   // 1B + 0B summed -> neutral
            styleTemplate(h_sig, TFColor::c2B());
                // Solid, not hatched: h_sig is cloned into the prefit and postfit stacks and
                // carries its fill style with it, so a hatch here desynced the two.


/* ---- disabled (kept for reference): normalized-PDF prefit canvas ----
            // normlaized data for comparison 
            TH1D* hnorm_data_self = (TH1D*)h_data_mb ->Clone("hnorm_data_self");
                hnorm_data_self->Scale(1./hnorm_data_self->Integral(1, mb_bins, "width"));
                hnorm_data_self->SetTitle("data self normalized");

                    auto canva_pdf_norm_beforefit = new TCanvas("Pdfs_norm_beforefit","", 800, 800 );
                        canva_pdf_norm_beforefit->cd();
                        h_bkg->Draw("hist E");
                        h_sig->Draw("hist E same");
                        hnorm_data_self->Draw("PE same");
                        canva_pdf_norm_beforefit->SetTitle("PDFs before fit");
                            gPad->Modified();   
                            gPad->Update();
                            canva_pdf_norm_beforefit->Modified();
                            canva_pdf_norm_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_pdf_norm_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            canva_pdf_norm_beforefit->Write();
                            canva_pdf_norm_beforefit->Print(Form("%s/%s_PDF_norm_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/ // ---- end disabled block ----
    /*        

            // what about a stack of PDFs before fit 
                // these are not normalized to 1, but normlaized such that the total PDFs are 1, using their qcd fractions  
                THStack hstack_pfds_scaledtoqcd (Form("hstack_pfds_scaledtoqcd_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true fractions, before fit");
                    hstack_pfds_scaledtoqcd.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true fractions in qcd 
                        TH1D* h_sig_scaled = (TH1D*) h_sig ->Clone("h_sig_scaled"); h_sig_scaled->Scale(sig_fraction_true);
                        TH1D* h_bkg_scaled = (TH1D*) h_bkg ->Clone("h_bkg_scaled"); h_bkg_scaled->Scale(1. - sig_fraction_true);
                            hstack_pfds_scaledtoqcd.Add(h_sig_scaled); 
                            hstack_pfds_scaledtoqcd.Add(h_bkg_scaled);
                            hstack_pfds_scaledtoqcd.SetMaximum(1.2 * hstack_pfds_scaledtoqcd.GetMaximum());

                            auto canva_pdfs_scaledtoqcd = new TCanvas("canva_pdfs_scaledtoqcd","", 800, 800 );
                                canva_pdfs_scaledtoqcd->cd();
                                hstack_pfds_scaledtoqcd.Draw("Hist E");
                                hnorm_data_self->Draw("HIST E same"); 
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd->Modified();
                                canva_pdfs_scaledtoqcd->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_scaledtoqcd->Write();
                                canva_pdfs_scaledtoqcd->Print(Form("%s/%s_pdfs_scaledtoqcdfractions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            

                // -- Prefit PDF: Signal, BKG
                THStack hstack_pfds_scaledtoqcd_int (Form("hstack_pfds_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        TH1D* h_sig_scaledint = (TH1D*) h_sig ->Clone("h_sig_scaledint"); h_sig_scaledint->Scale(int2);
                        TH1D* h_bkg_scaledint = (TH1D*) h_bkg ->Clone("h_bkg_scaledint"); h_bkg_scaledint->Scale(int1 + int0);
                            hstack_pfds_scaledtoqcd_int.Add(h_sig_scaledint); 
                            hstack_pfds_scaledtoqcd_int.Add(h_bkg_scaledint);
                            hstack_pfds_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_scaledtoqcd_int.GetMaximum());

                            auto canva_pdfs_scaledtoqcd_int = new TCanvas("canva_pdfs_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_scaledtoqcd_int->cd();
                                hstack_pfds_scaledtoqcd_int.Draw("Hist E");
                                TH1D* hd_scaledtoqcdint = (TH1D*) hnorm_data_self->Clone("hd_scaledtoqcdint"); hd_scaledtoqcdint->Scale(int0+int1+int2);
                                        hd_scaledtoqcdint->SetTitle("data scaled to qcd integral");
                                        hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd_int->Modified();
                                canva_pdfs_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_scaledtoqcd_int->Write();
                                canva_pdfs_scaledtoqcd_int->Print(Form("%s/%s_pdfs_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
                                    // test 0B: (1B + 2B) after scaling to dijet fractions: is it like the true fractions?
                                    cout << "----------------" << endl;
                                    cout << "After combine Dijet+ bjet, and scaling to qcd integrals " << endl;
                                    cout << "Integarl of scaled PDF(2B) = " <<   h_sig_scaledint->Integral() << endl;
                                    cout << "Integral of scaled PDF(1B+0B) = " << h_bkg_scaledint->Integral() << endl;

                // -- Prefit PDF: signal, BKG (seperated 0B and 1B)
                 THStack hstack_pfds_seperated_scaledtoqcd_int (Form("hstack_pfds_seperated_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_seperated_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        // 1B (bjet and dijet)
                        TH1D* hnorm_sumbkg_1B_scaledint = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledint");
                                hnorm_sumbkg_1B_scaledint->Scale(1./hnorm_sumbkg_1B_scaledint->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledint->Scale(int1);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_sig_scaledint); 
                        hstack_pfds_seperated_scaledtoqcd_int.Add(hnorm_sumbkg_1B_scaledint);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_nob);
                        hstack_pfds_seperated_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_seperated_scaledtoqcd_int.GetMaximum());
                            auto canva_pdfs_seperated_scaledtoqcd_int = new TCanvas("canva_pdfs_seperated_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_seperated_scaledtoqcd_int->cd();
                                hstack_pfds_seperated_scaledtoqcd_int.Draw("Hist E");
                                hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_seperated_scaledtoqcd_int->Modified();
                                canva_pdfs_seperated_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_seperated_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                canva_pdfs_seperated_scaledtoqcd_int->Write();
                                canva_pdfs_seperated_scaledtoqcd_int->Print(Form("%s/%s_pdfs_seperated_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
    */                             
            // Before fit: pdfs scaled to data 
             THStack hstack_pfds_scaledtoData (Form("hstack_pfds_scaledtoData_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to data, before fit");
                    hstack_pfds_scaledtoData.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the data inetgral 
                        // normalized then multiplied by their fraction from qcd, then scaled to data
                        TH1D* hnorm_sumbkg_1B_scaledtoData = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledtoData");
                                hnorm_sumbkg_1B_scaledtoData->Scale(1./hnorm_sumbkg_1B_scaledtoData->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledtoData->Scale(integral_inputdata * int1/tot);
                        TH1D* h_nob_scaledtoData = (TH1D*) norm_h_nob->Clone("h_nob_scaledtoData");
                                 h_nob_scaledtoData->Scale(integral_inputdata * int0/tot);
                        TH1D* h_2B_scaledtoData = (TH1D*) h_sig->Clone("h_2B_scaledtoData");
                                h_2B_scaledtoData->Scale(integral_inputdata * int2/tot);
                        hstack_pfds_scaledtoData.Add(h_2B_scaledtoData); // 2B 
                        hstack_pfds_scaledtoData.Add(hnorm_sumbkg_1B_scaledtoData);
                        hstack_pfds_scaledtoData.Add(h_nob_scaledtoData); // 0B
                        hstack_pfds_scaledtoData.SetMaximum(1.6 * h_data_mb->GetMaximum());

            // Modify prefit to add Ratio : Data/total MC 
               TCanvas* cRatio_prefit = new TCanvas(Form("Prefit_RatioPlot_%s", sname_canvas.Data()), "", 1100, 1100);
                    TPad* pad11 = new TPad("pad11","",0,0.2,1,1);//  0,0.2,1,1
                    TPad* pad22 = new TPad("pad21","",0,0,1,0.24);
                    pad11->SetBottomMargin(0.07); // 0.06 
                    pad11->SetLeftMargin(0.18); // for y axis title space 
                    pad11->SetTopMargin(0.14); // new to allow cms label
                    pad22->SetTopMargin(0.03);// 0.02     // bottom pad (very small)
                    pad22->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad22->SetLeftMargin(0.18);
                    pad11->SetTicks(1, 1); pad22->SetTicks(1, 1);   // ticks on all four sides
                    pad11->Draw();
                    pad22->Draw(); 
                    pad11->cd();
                    // Draw frame to control y axis name: frame needed to control y axis name 
                    TH1F *frame_pre = pad11->DrawFrame(h_data_mb->GetXaxis()->GetXmin(), 0, h_data_mb->GetXaxis()->GetXmax(), h_data_mb->GetMaximum() * 1.6 );
                        frame_pre->GetYaxis()->SetTitle("Counts/[GeV]");
                         frame_pre->GetXaxis()->SetLabelSize(0); // Remove labels of this Pad 
                        hstack_pfds_scaledtoData.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad11, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg_pre = CreateLegend(0.70, 0.6, 0.85, 0.75, // 0.63, 0.6, 0.85, 0.75, 
                            {h_data_mb, h_sig, h_sumbkg, norm_h_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "2B", "1B", "0B"} 
                        );
                        leg_pre->Draw("same");
                        // Real data on this canvas -> "Internal", not "Simulation Internal".
                        drawCMSHeader(pad11);
                        pad11->Modified(); // force refresh
                        pad11->Update();
                    pad22->cd(); 
                    // Draw ratio: data/total MC stack
                    TH1D* hsatckMC_total = (TH1D*) hstack_pfds_scaledtoData.GetStack()->Last()->Clone("hsatckMC_total");
                    AddRatioPlot(h_data_mb, hsatckMC_total, "Data/MC");
                    pad22->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        // cRatio_prefit->Write(); // canvas not stored in the root file
                        if (save_png) cRatio_prefit->Print(Form("%s/Prefit_%s.png", sDir_canvas.Data(), sname_canvas.Data()));
                        // was: a second copy into the per-variation _www subfolder, for non-integrated bins only


            ///// Fitting
            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     
            // Inputs 
            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
            RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
            RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
            RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);

            // Create list of templates
            RooArgList template_list(sig_template, bkg_template, "template_list");

            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
                                                                                                                                                                          //instead of sign histogram here we would put the data histogram
            Int_t status = result->status();
            if (status != 0) {
                // One retry before the bin is given up (2026-10-26): Strategy 2 (full Hessian at
                // every step), started from the middle of [0,1] instead of the MC fraction.
                cout << Form("WARNING: (ptbin %d, bin %d) fit status %d -- retrying with Strategy 2 from a = 0.5",
                             ibin_pt, ibin_dr, status) << endl;
                sig_fraction_val.setVal(0.5);
                delete result;
                result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(2), RooFit::Minos(false));
                status = result->status();
            }
            result->Print();

            std::cout << "covariance matrix:" << std::endl;
            (result->covarianceMatrix().Print());

            //Check if it converged for a dr and jtpt bin
            if (status != 0) {
                std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                         << ", ix = " << ibin_dr 
                         << " did not converge\n\n\n\n" << std::endl;
                non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();


            //Print some check
            std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a; // signal 2B
            // Total bkg: (1-a)
            p1 = (1-a)*eff_bkg1B; // 1B bkg 
            p2 = (1-a)*eff_bkg0B; // 0B bkg 
            errP0 = da;
            errP1 = da *eff_bkg1B;
            errP2 = da *eff_bkg0B;

            // std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
            // std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


            // std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
            // std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


            // std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
            // std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;


            // -- updated for the new hsit binning for the true and fitresult S/B fractions (for later drawings)
            // Fit result 
            h_sig_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, p0);
            h_sig_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr +1, ibin_pt+1, errP0);
            h_bkg_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, 1- p0);
            h_bkg_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0); // errorp1 = error (1-p0)
            h_bkg_fraction_error->SetBinContent(ibin_dr +1, ibin_pt +1, errP0);

                // For non integaretd bins 
            if(!(ibin_dr == 0 || ibin_pt == 0)){
                h_sig_fraction_fit->SetBinContent(ibin_dr , ibin_pt , p0);
                h_sig_fraction_fit->SetBinError(ibin_dr, ibin_pt , errP0);
                h_sig_fraction_fit_error->SetBinContent(ibin_dr, ibin_pt, errP0);
                h_bkg_fraction_fit->SetBinContent(ibin_dr, ibin_pt, 1- p0);
                h_bkg_fraction_fit->SetBinError(ibin_dr , ibin_pt , errP0);
                h_bkg_fraction_fit_error->SetBinContent(ibin_dr, ibin_pt, errP0);
            }



            // -- Compute siganl and bkg fractions uncertainity 
            //save the true fraction
            // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
            // 2026-10-26: in quadrature (S and B are independent samples); was the linear sum
            // ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/(S+B)^2.
            double err_true_frac = std::hypot((int0 + int1) * True_sig_err, int2 * True_bkg_err)/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B
            h_sig_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, sig_fraction_true);
            h_sig_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, bkg_fraction_true);
            h_bkg_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true_error->SetBinContent(ibin_dr +1, ibin_pt +1, err_true_frac);
            
            if(!(ibin_dr == 0 || ibin_pt == 0))
            {
                h_sig_frac_true_fitbins ->SetBinContent(ibin_dr , ibin_pt, sig_fraction_true);
                h_sig_frac_true_fitbins ->SetBinError(ibin_dr, ibin_pt , err_true_frac);
                h_bkg_frac_true_fitbins ->SetBinContent(ibin_dr , ibin_pt, bkg_fraction_true);
                h_bkg_frac_true_fitbins ->SetBinError(ibin_dr, ibin_pt , err_true_frac);
                h_bkg_frac_true_error_fitbins ->SetBinContent(ibin_dr , ibin_pt , err_true_frac);

            }   

            // -- After fits: save mass distribution re and post-fit 
            TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("2B");
                h_sig_fit->GetYaxis()->SetTitle("Counts/[GeV^{2}]");

                styleTemplate(h_sig_fit, TFColor::c2B());   // 2B -> red
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("1B+ 0B");
                styleTemplate(h_bkg_fit, TFColor::bkg());   // 1B + 0B summed -> neutral

            // 1B and 0B slices with their OWN shapes (2026-10-26): both were clones of h_bkg, the
            // 1B+0B mixture, so they were drawn and written with one and the same shape. They
            // still add up to h_bkg_fit: p1 + p2 = 1-a, mixed exactly as h_bkg was.
            TH1D *h_bkg_fit_1b = (TH1D*) norm_h_1b->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) norm_h_nob->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("0B");
                styleTemplate(h_bkg_fit_nob, TFColor::c0B());  // 0B -> green

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("Total fit");
                    // The sum of the three slices, not a category of its own: neutral dark
                    // outline and no fill, so it cannot be misread as one of them.
                    h_total_fit->SetFillColor(0);
                    h_total_fit->SetFillStyle(0);
                    h_total_fit->SetMarkerColor(TFColor::total());
                    h_total_fit->SetLineColor(TFColor::total());
                    h_total_fit->SetMarkerStyle(1);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();     // 2B
            h_bkg_fit_1b->Write();  // 1B
            h_bkg_fit_nob->Write(); // 0B
            h_bkg_fit->Write();     // 1B + 0B
            h_total_fit->Write();   // 2B + 1B + 0B

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
/* ---- disabled (kept for reference): after-fit overlay canvas (data vs. total fit vs. components) ----
            auto canva_afterfit = new TCanvas(Form("ALLHist_%s", sname_canvas_afterfit.Data()) ,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetTitle("Data");
                h_data_mb->SetLineWidth(2);
                if(also_bjet) h_bb_bjet->SetLineWidth(2);
                h_sig_fit->SetLineWidth(2);
                if(also_bjet) h_b_bjet->SetLineStyle(9);
                h_bkg_fit->SetLineWidth(2);
                if (h_data_mb->GetMaximum() < h_total_fit->GetMaximum() ) { h_data_mb->SetMaximum(1.3 * h_total_fit->GetMaximum());}
                else { h_data_mb->SetMaximum(1.3 * h_data_mb->GetMaximum()); }

                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");
                h_bkg_fit_1b->Draw("HIST E same");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                canva_afterfit->SetTitle("");
                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               
*/ // ---- end disabled block ----
                h_sig_fit->SetLineWidth(1); // stack-slice outline: 1, so it does not read as a curve

        // -- After fit: stacked: seperated contibutions 
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.Add(h_bkg_fit_nob);
                hstack_afterfit.SetMaximum(1.2 * hstack_afterfit.GetMaximum());

/* ---- disabled (kept for reference): after-fit control canvases (stacked / normalized stack) ----
            auto canva_stack_afterfit = new TCanvas(Form("All_templates_Data_stacked_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit: Data Vs. pdf signal Vs. bkg 
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_normstack = (TH1D*)  h_bkg_fit->Clone("h_bkg_fit_normstack");  h_bkg_fit_normstack ->Scale(1./integral_inputdata);
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_normstack);
                hstack_norm_afterfit.SetMaximum(1.3 * hnorm_data_self->GetMaximum());

                auto canva_stack_norm_afterfit = new TCanvas(Form("PDFs_Data_stacked_norm_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hnorm_data_self->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            //  -- Build legend and write plots     
                fout->cd();
                DrawCommonTextTopRight(canva_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_afterfit,ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_norm_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_norm_afterfit->GetName()));
*/ // ---- end disabled block ----

            //-- ratio plot: data /total fit
                TCanvas* c = new TCanvas(Form("RatioPlot_%s", sname_canvas_afterfit.Data()), "", 1100, 1100); // 900, 1100
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);//  0,0.2,1,1
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    pad1->SetBottomMargin(0.07);
                    pad1->SetLeftMargin(0.18); // for y axis title space 
                    pad1->SetTopMargin(0.14); // new to allow cms label
                    pad2->SetTopMargin(0.03);// 0.02     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2->SetLeftMargin(0.18);
                    pad1->SetTicks(1, 1); pad2->SetTicks(1, 1);   // ticks on all four sides
                    pad1->Draw();
                    pad2->Draw();
                    pad1->cd();
                    // Frame controls the y axis name. Range comes from h_data_mb, not from
                    // hstack_afterfit: THStack::GetXaxis() is null until the stack is painted,
                    // and the stack sums to h_total_fit, so this is equivalent.
                    Double_t ymax_postfit = 1.3 * std::max(h_data_mb->GetMaximum(), h_total_fit->GetMaximum());
                    TH1F *frame = pad1->DrawFrame(h_data_mb->GetXaxis()->GetXmin(), 0, h_data_mb->GetXaxis()->GetXmax(), ymax_postfit);
                        frame->GetYaxis()->SetTitle("Counts/[GeV]");
                        frame->GetXaxis()->SetLabelSize(0); // Remove labels of this Pad
                        hstack_afterfit.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // h_total_fit->Draw("P E SAME");
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad1, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg = CreateLegend(0.70, 0.6, 0.85, 0.75, // 0.63, 0.6, 0.85, 0.75, 
                            {h_data_mb, h_sig_fit, h_bkg_fit_1b, h_bkg_fit_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "", "", ""} // use default titles 
                        );
                        leg->Draw("same");
                        // Real data on this canvas -> "Internal", not "Simulation Internal".
                        drawCMSHeader(pad1);
                        pad1->Modified(); // force refresh
                        pad1->Update();
                    pad2->cd(); 
                // setCMSStyle(); // also for Pad2 
                    AddRatioPlot(h_data_mb, h_total_fit);
                    pad2->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        // c ->Write(); // canvas not stored in the root file
                        if (save_png) c ->Print(Form("%s/%s.png", sDir_canvas.Data(), c->GetName())); // RatioPlot_ptbin_N_deltaRbin_M_afterfit.png
                        // was: a second copy into the per-variation _www subfolder, for non-integrated bins only

                cout << "---------------------\n\n\n" << endl; 
            } // loop over deltaR bins 
    }

/* ---- disabled (kept for reference): duplicate TH3D write, inputs already saved before the fit loop ----
    // Save histograms
    // TH3D 
    for (auto h : {h3D_data, h3D_bb, h3D_b}){h->Write();}
    if (also_bjet){
        for (auto h : {h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    }   
*/ // ---- end disabled block ----



    // TH2D: signal / background fractions. Bin errors live inside these histograms, so the
    // separate *_error clones are no longer written.
    // 2026-10-26: a bin with no fit result (failed twice, or data but no usable template) would
    // keep signal fraction 0 in h_sig_fraction_fit, and apply_unfolding_2d.C multiplies the data
    // by it -- that data bin would be emptied without a word. So the two *_fraction_fit
    // histograms the unfolding reads are NOT written then: it stops on the missing histogram.
    const bool missing_fractions = !non_converge_bins.empty();
    fout->cd();
    for (auto h : {h_sig_fraction,
                   h_bkg_fraction,
                   h_sig_frac_true,
                   h_bkg_frac_true,

                   h_sig_fraction_fit,
                   h_bkg_fraction_fit,
                   h_sig_frac_true_fitbins,
                   h_bkg_frac_true_fitbins
                   }) {
                    if (missing_fractions && (h == h_sig_fraction_fit || h == h_bkg_fraction_fit)) continue;
                    h->Write();
    }

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
    if (missing_fractions)
        std::cerr << "ERROR: " << non_converge_bins.size() << " (pT, " << obs.name << ") bin(s) without a signal fraction ("
                  << varNames[ivar] << ") -- h_sig_fraction_fit NOT written to " << fout_name << std::endl;
        
    fout->Close();

}


void Draw_template_Run3(TString &templates, TString pT_selection, TString folder, TString &fout_name, Variation ivar = NOMINAL){
    // -- modified function. 
    // -- Draw only prefits. One input sample.

    bool also_bjet = false;

    // -- Make subdirectory for printed canvases only. The root files are in main directory 
    TString sDir_canvas = Form("%s/%s", sDirname.Data(), varNames[ivar].Data());
    gSystem->mkdir(sDir_canvas, kTRUE); // to be clean 

    // -- For output histograms 
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");

    // -- Using new naming convention
    // 2026-10-26: names and binning from gFitObs(), as in do_template_fit_combined() -- were
    // the dR ones, so a B call read dR histograms (and failed CheckInputBinning below).
    const ObsDef &obs = gFitObs();
    TString namehData = obs.n("h3D_data"); // h3D_data with eec on the weight
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = obs.n("h3D_bb");
    TString nameh1B = obs.n("h3D_b");

    TString nameh0B = obs.n("h3D_0b");
    // -- Additional histograms
    TString namehmore2B = "";
 
    //-- Dijet sample:
    TFile *file_dijet = new TFile(templates, "read");
            if (!file_dijet) {Error("Input File:", "File does not exist'%s'", file_dijet->GetName());}
        TH3D *h3D_b = (TH3D*)file_dijet->Get(nameh1B)->Clone("h3D_b");
        TH3D *h3D_bb = (TH3D*)file_dijet->Get(nameh2B)->Clone("h3D_bb");
            if(!h3D_b){Error("Get:", "histogram does not exist '%s' ",h3D_b->GetName()); return;}
            if(!h3D_bb){Error("Get:", "histogram does not exist '%s' ",h3D_bb->GetName()); return;}
        TH3D *h3D_nob = (TH3D*)file_dijet->Get(nameh0B)->Clone("h3D_0b");  // clone to the
            // CANONICAL name, not nameh0B: the input histogram carries the observable
            // suffix ("h3D_0b_B"), but everything downstream -- Draw_EEC.h in particular
            // -- looks the templates up in the fit output by their unsuffixed names, the
            // same way h3D_b and h3D_bb are cloned just above. For dR the two spellings
            // coincide, which is why only the balance fit tripped over it.
            if(!h3D_nob){Error("Get:", "histogram does not exist '%s' ",h3D_nob->GetName()); return;}
/*
    //-- Bjet sample:
    TFile *file_bjet = new TFile(templates_bjet, "read");
         if (!file_bjet) {Error("Input File:", "File does not exist'%s'", file_bjet->GetName());}
        TH3D *h3D_b_bjet  = (TH3D*)file_bjet->Get(nameh1B)->Clone("h3D_b_bjet");
        TH3D *h3D_bb_bjet = (TH3D*)file_bjet->Get(nameh2B)->Clone("h3D_bb_bjet");
        if(!h3D_b_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_b_bjet->GetName()); return;}
        if(!h3D_bb_bjet){Error("Get:", "histogram does not exist '%s' ",h3D_bb_bjet->GetName()); return;}

    //Open dataset:
    TFile *file_data = new TFile(HighEGdata_name, "read");
        if (!file_data) {Error("Input File:", "File does not exist'%s'", file_data->GetName());return;}
            cout << "file data name " << file_data->GetName() << endl;
        TH3D *h3D_data;
            h3D_data = (TH3D*)file_data->Get(namehData)->Clone("h3D_data");
            if(!h3D_data){Error("Get:", "histogram does not exist '%s' ",h3D_data->GetName()); return;}
            if(alsoLowEG){
                TFile *file_data_LowEG = new TFile(LowEGdata_name, "read");
                    if (!file_data_LowEG) {Error("Input File:", "File does not exist'%s'", file_data_LowEG->GetName());return;}
                    cout << "file data name " << file_data_LowEG->GetName() << endl;
                TH3D *h3D_dataLowEG;
                    h3D_dataLowEG = (TH3D*)file_data_LowEG->Get(namehData)->Clone("h3D_dataLowEG");
                    if(!h3D_dataLowEG){Error("Get:", "histogram does not exist '%s' ",h3D_dataLowEG->GetName()); return;}
                h3D_data->Add(h3D_dataLowEG);
            }
        // -- data: set style
        h3D_data->SetTitle("Data");
        styleData(h3D_data);

        // bjet
        h3D_b_bjet->SetTitle("1B (bjet)");h3D_bb_bjet->SetTitle("2B (bjet)");
        // its x axis name is not correct
        h3D_b_bjet->GetXaxis()->SetTitle("m_{2B} [GeV]"); // Name is wrong in MC templates (EEC)
        styleTemplate(h3D_b_bjet,  TFColor::c1B()); h3D_b_bjet ->SetFillStyle(3354);
        styleTemplate(h3D_bb_bjet, TFColor::c2B()); h3D_bb_bjet->SetFillStyle(3345);
*/
        // qcd
        h3D_b->SetTitle("1B (qcd)");h3D_bb->SetTitle("2B (qcd)"); h3D_nob ->SetTitle("0B (qcd)");
        styleTemplate(h3D_b,   TFColor::c1B());
        styleTemplate(h3D_bb,  TFColor::c2B());
        styleTemplate(h3D_nob, TFColor::c0B());

    // -- Define used bins. As in do_template_fit_combined(): validate the inputs against
    // binning_histos_small.h instead of adopting their axes (these are globals).
    if (!CheckInputBinning(h3D_bb))  return;
    if (!CheckInputBinning(h3D_b))   return;
    if (!CheckInputBinning(h3D_nob)) return;

        bins_pt = jtpt_binsVectorSize - 1;
        bins_dr = obs.nbins;   // was dr_binsVectorSize - 1 (2026-10-26)
        mb_bins = mb_binsVectorSize   - 1;

        cout << "-- Data hist initial binning" << endl;
        cout << "pt bins = "<< bins_pt << endl;
        cout << "dr bins = "<< bins_dr << endl;
        cout << "mb bins = "<< mb_bins << endl;

    // -- Write used inputs 
    fout->cd();
    // h3D_data->Write();
    // h3D_b_bjet->Write(); 
    // h3D_bb_bjet->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h3D_nob->Write(); 

    const double* yBins = nullptr; //dr array 
    int N_bins_dr = 0; 
    N_bins_dr = bins_dr;
    yBins = obs.bins;   // was dr_binsVector (2026-10-26)
/*
    // ---------------------------------------------------------
    /// -- For later drawing of S/B fractions, store true and fit result in the following histograms
        // Note that: since inetgaretd bins are pt 0 and dr 0, the hist of fractions should have #bins + 1 size 
        // x = dr, y = jetpt
        // Axis here is for the bin number instead of the values 
          TH2D *h_sig_fraction = new TH2D("h_sig_fraction", ";dr; jet pt", N_bins_dr+ 1,  1, N_bins_dr+ 2,
               h3D_data->GetNbinsZ() + 1 , 1, h3D_data->GetNbinsZ()+ 2 );                                                                                                                              
            h_sig_fraction->Reset();
            TH2D *h_bkg_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_fraction");
            TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
            TH2D *h_bkg_fraction_error = (TH2D *) h_bkg_fraction->Clone("h_bkg_fraction_error"); 
            TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
            TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");
            TH2D *h_bkg_frac_true = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true");
            TH2D *h_bkg_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_frac_true_error");
*/
    //-----------------------
    // -- general style: same CMS aesthetics as do_template_fit_combined()
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);

    // ---------------------------------------------------------
    // --- Vector to test the convergence
    std::vector <std::pair<int, int>> non_converge_bins;
 
    // Fitting - loop over dr and jtpt entries
    // Bin0: is integarted over the range 
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
    // for(Int_t ibin_pt = 1 ; ibin_pt <= 1; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= N_bins_dr; ibin_dr++){
        // for(Int_t ibin_dr = 1; ibin_dr <= 0; ibin_dr++){
            
            // define slice
            Int_t SliceFirstbin_dr = ibin_dr;
            Int_t SliceLastbin_dr =  ibin_dr;
            Int_t SliceFirstbin_pt = ibin_pt;
            Int_t SliceLastbin_pt =  ibin_pt;

            if(!ibin_dr){SliceFirstbin_dr = 1; SliceLastbin_dr = N_bins_dr;}
            if (!ibin_pt){ SliceFirstbin_pt = 1;  SliceLastbin_pt = bins_pt;}

            /*
            // Make projections 
            TH1D *h_data_mb = (TH1D *) h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
                h_data_mb->GetXaxis()->SetTitle("m_{2B} [GeV]");
                h_data_mb->SetTitle(h3D_data->GetTitle()); // upadte projection title
            */

            // Make slices for dijet
            TH1D *h_bb = (TH1D *) h3D_bb->ProjectionX(Form("h_bb_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_b = (TH1D *) h3D_b->ProjectionX(Form("h_b_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
            TH1D *h_nob = (TH1D *) h3D_nob ->ProjectionX(Form("h_nob_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, ibin_pt, ibin_pt);
                h_bb->SetTitle( h3D_bb->GetTitle() );
                h_b->SetTitle( h3D_b->GetTitle() );
                h_nob->SetTitle( h3D_nob->GetTitle() );

            // -- instead of the commented def.
                TH1D *h_bb_bjet; 
                TH1D *h_b_bjet;

/*
            //  Make slices for bjet 
            TH1D *h_bb_bjet = (TH1D *) h3D_bb_bjet ->ProjectionX(Form("h_bb_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr, SliceFirstbin_pt, SliceLastbin_pt);
            TH1D *h_b_bjet = (TH1D *) h3D_b_bjet ->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), SliceFirstbin_dr, SliceLastbin_dr,  SliceFirstbin_pt, SliceLastbin_pt);
                h_b_bjet->SetTitle(h3D_b_bjet->GetTitle());
                h_bb_bjet->SetTitle(h3D_bb_bjet->GetTitle());
*/



            // -- Deattach hists from root file 
            h_b->SetDirectory(0);
            h_bb->SetDirectory(0);
            h_nob->SetDirectory(0);
            /*
            h_data_mb->SetDirectory(0);
            h_b_bjet->SetDirectory(0);
            h_bb_bjet->SetDirectory(0);
            */

            // string (dr, pt) 
            // Plot filenames name the observable. dR keeps the literal "deltaR" it has always
            // had, so its PNG names are unchanged; the balance gets "Bbin".
            // Draw_template_Run3 has no local ObsDef, so read the driver-set global.
            TString sname_canvas = Form("ptbin_%d_%sbin_%d", ibin_pt,
                                        (gFitObs().name == "dr") ? "deltaR" : gFitObs().name.Data(), ibin_dr);

            // --  Calculate true fractions to be used as initial values for the fit (the true fractions are the qcd ones)            
            double int2 = h_bb->Integral(1, mb_bins, "width");
            double int1 = h_b ->Integral(1, mb_bins, "width");
            double int0 = h_nob->Integral(1, mb_bins, "width"); 
            double tot = int0 + int1 + int2;

                // compute the true fractions and their errors (for later comparison)           
                double sig_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : (int2 / (int0 + int1 + int2)); 
                double bkg_fraction_b_true = (int0  + int1 + int2 ) == 0 ? 0 : (int1 / (int0 + int1 + int2));
                double bkg_fraction_true = (int0  + int1 + int2 ) == 0 ? 0 : ( (int0 + int1)/(int0 + int1 + int2) );

                double True_bkg_b_err = 0.0;
                double True_bkg_b_integral = h_b->IntegralAndError(1, mb_bins, True_bkg_b_err, "width");

                double True_sig_err = 0.0;
                double True_sig_integral = h_bb->IntegralAndError(1, mb_bins, True_sig_err, "width");

                // total bkg fraction: 0B + 1B (qcd only)
                TH1D* h_b_nob_dijet = (TH1D*) h_b ->Clone("h_b_nob_dijet");  h_b_nob_dijet->Add(h_nob);
                double True_bkg_err = 0.0;
                double True_bkg_integral = h_b_nob_dijet->IntegralAndError(1, mb_bins, True_bkg_err, "width");


            // -- Compute other useful integrals 
            /*    
            // From bjet sample 
            double int2_bjet = h_bb_bjet->Integral(1, mb_bins, "width"); 
            double int1_bjet = h_b_bjet->Integral(1, mb_bins, "width");
            // data:
            double integral_inputdata = h_data_mb->Integral(1, mb_bins,"width");
            */

                // -- check integrals in bjet and qcd sample
                // relative fraction of 1B: 2B in bjet and qcd samples is SAME. So you can combine the two samples with simple + (without reweighting).
                // cout << "Total input data integral = "<< integral_inputdata << endl;
                // std::cout << "int2 h_bb=" << int2 << std::endl;
                // std::cout << "int1 h_b =" << int1 << std::endl;
                // cout << "Dijet: int0 of 0B" << int0 << endl;

                // cout << "Dijet: 2B/(1B + 2B) = " << int2/(int2 + int1) << endl;
                // std::cout << "int2 h_bb_bjet=" << int2_bjet << std::endl;
                // std::cout << "int1 h_b_bjet =" << int1_bjet << std::endl;
                // cout << "2B/(2B + 1B) in qcd sample = " << int2_bjet/(int2_bjet + int1_bjet) << endl;
                // cout << " Dijet: 0B/(1B+2B) = " << int0/(int2 + int1) << endl;

            // --- Prepare PDFs for template fit

            // 1- Combine qcd + bjet samples 
            // Signal: 2B   
            TH1D *h_sumsig = (TH1D*) h_bb->Clone(Form("h_sumsig_%d_%d", ibin_dr, ibin_pt));
                if(also_bjet){ h_sumsig->Add(h_bb_bjet);h_sumsig->SetTitle("2B: qcd+bjet");}
                Int_t h_sig_bins = h_sumsig->GetNbinsX();
                TH1D* h_sig = (TH1D*) h_sumsig->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

            // Bkg: first 1B, and later added effectively 0B
            TH1D *h_sumbkg = (TH1D*) h_b->Clone(Form("h_sumbkg_%d_%d", ibin_dr, ibin_pt));
                 if(also_bjet) { h_sumbkg->Add(h_b_bjet); h_sumbkg->SetTitle("1B: qcd+bjet");}
                 else{ h_sumbkg->SetTitle("1B: qcd"); }
                Int_t h_bkg_bins = h_sumbkg->GetNbinsX();
                TH1D* h_bkg = (TH1D*)  h_sumbkg->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt)); // to be normalized (next step)

                // for drawings
                    // Total sum bkg = 0B + 1B  
                    TH1D* h_sumbkg_0b_1b = (TH1D*) h_sumbkg->Clone(Form("h_sumbkg_0b_1b_%d_%d", ibin_dr, ibin_pt));
                    h_sumbkg_0b_1b->Add(h_nob);
                    if (also_bjet) h_sumbkg_0b_1b->SetTitle("0B +1B (qcd+bjet)");
                    else {h_sumbkg_0b_1b->SetTitle("0B +1B (qcd)"); }


                // -- set the combined samples styles before fit
                styleTemplate(h_sumbkg, TFColor::c1B());   // 1B
                styleTemplate(h_sumsig, TFColor::c2B());   // 2B
                    // and for the bkg 1B + 0B: a sum of two categories -> neutral
                    styleTemplate(h_sumbkg_0b_1b, TFColor::bkg());

                // write to rootfile the used slices 
                fout->cd();
                // h_data_mb->Write(); // data 
                h_sumsig->Write(); // 2B (qcd + bjet)
                h_sumbkg->Write(); // 1B (qcd + bjet)
                h_sumbkg_0b_1b->Write(); // 1B (qcd + bjet) + 0B 


    // -------- Draw prefits 
        // contirbutions seperated  
        // absolute yields: seperated contributions: qcd, bjet 
            THStack hstack_all_beforefit (Form("hstack_all_beforefit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms before fit");
                hstack_all_beforefit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                // hstack_all_beforefit.Add(h_bb_bjet);
                hstack_all_beforefit.Add(h_bb);
                // hstack_all_beforefit.Add(h_b_bjet); 
                hstack_all_beforefit.Add(h_b);
                hstack_all_beforefit.Add(h_nob);
                // if (hstack_all_beforefit.GetMaximum() > h_data_mb->GetMaximum()/1e+04) { hstack_all_beforefit.SetMaximum(1.3* hstack_all_beforefit.GetMaximum());}
                // else { hstack_all_beforefit.SetMaximum(1.3* h_data_mb->GetMaximum()/1e+04);}
                auto canva_beforefit = new TCanvas(Form("All_contributions_beforefit_%d_%d", ibin_dr, ibin_pt),"", 800, 800 );
                        canva_beforefit->cd();
                        // canva_beforefit->SetLogy();
                        hstack_all_beforefit.Draw("hist E");
                        /*TH1D* hd_norm_slice = (TH1D*) h_data_mb->Clone("hd_norm_slice");
                            hd_norm_slice->Scale(1./10000);
                            hd_norm_slice->Draw("HIST E same"); 
                            hd_norm_slice->SetTitle("data/1e+04");
                        */
                            gPad->Modified();   
                            gPad->Update();
                            canva_beforefit->Modified();
                            canva_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_beforefit, "Work In Progress");
                            canva_beforefit->Write();
                            canva_beforefit->Print(Form("%s/%s_allcontributions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Draw Stack for combined qcd + bjet: 1B, Combined 2B, and 0b : the real contrituions that is to be fitted before normalization
            // absolute yields of Signal and bkg  
            THStack hstack_templatesforfit(Form("hstack_templatesforfit_%d_%d", ibin_dr, ibin_pt),"Mass stacked histograms without normalization");
                    hstack_templatesforfit.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                    hstack_templatesforfit.Add(h_sumsig);
                    hstack_templatesforfit.Add(h_sumbkg_0b_1b);    
                    // hstack_templatesforfit.SetMaximum(1.2 * hd_norm_slice->GetMaximum());
                
                    auto canva_sum_beforefit = new TCanvas(Form("templates_beforefit_%d_%d", ibin_dr, ibin_pt), "", 800, 800 );
                        canva_sum_beforefit->cd();
                        hstack_templatesforfit.Draw("hist E");
                        // hd_norm_slice->Draw("HIST E same"); 
                            gPad->Modified();   
                            gPad->Update();
                            canva_sum_beforefit->Modified();
                            canva_sum_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_sum_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_sum_beforefit, "Work In Progress");
                            canva_sum_beforefit->Write();
                            canva_sum_beforefit->Print(Form("%s/%s_templates_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));

        // -- Safety for empty bins
            // Scale-free: "empty" must mean NO ENTRIES, not "integral below 1". The old
            // test was an absolute < 1 on the summed template integrals, which is only ever
            // a no-op for the EEC-weighted fit, where (pt1*pt2)^n makes those integrals
            // ~1e7. Run the same fit on UNWEIGHTED templates (EEC_WEIGHT_OFF, the yield
            // measurement) and the integrals are order 1 -- a whole QCD block sums to 0.86 --
            // so the sparser high-dr bins tripped this and were silently skipped, leaving
            // their signal fraction at exactly 0. Comparing against 0 asks the question the
            // comment always claimed to ask, and is a no-op for the EEC fit.
            if ((int0 + int1 + int2) <= 0. ) { cout << " ----------- empty bin -----------  "; continue;}   
                // ----- To avoid empty bins if exist!, set them to eps value
                const double eps = 1e-6; 
                for (int i = 1; i <= h_sig_bins; i++){
                    if (h_sig->GetBinContent(i) <= 0){
                        h_sig->SetBinContent(i, eps);
                        cout << 
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty  signal Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i) 
                        << endl; 
                    }
                }

                for (int i = 1; i <= h_bkg_bins; i++){
                    if (h_bkg->GetBinContent(i) <= 0){
                        h_bkg->SetBinContent(i, eps);
                        cout <<
                        Form("INFO: (ptbin %d, deltaRbin %d) has empty bkg  Bin %d, value set to 1e-06 to avoid fit failur", ibin_pt, ibin_dr, i)  
                        << endl; 
                    }
                }
                /// Normalize safely for the Roofit
                if(h_sig->Integral(1, h_sig_bins, "width") == 0.0 || h_bkg->Integral(1, h_bkg_bins, "width") == 0.0)
                {
                    cerr<< "Signal or bkg templates has zero integral"<< endl;
                    continue;
                }
                // normalize signal and bkg 
                h_sig->Scale(1/h_sig->Integral(1, h_sig_bins, "width"));
                h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins, "width"));
                    // cout << "After Normalization (qcd +bjet) True 2B integral = " << h_sig->Integral(1, h_sig_bins, "width") << endl;
                    // cout << "After Normalization (qcd+bjet)  True 1B integral =  "<<  h_bkg->Integral(1, h_bkg_bins, "width") << endl;

                // -- Effective bkg PDF. With a + b + c = 1 (a: 2B, b: 1B, c: 0B) and
                // b' + c' = 1, this gives b' = b/(b+c) and c' = c/(b+c), with b+c = 1-a.
                // NOTE: a, b, c come from int2/int1/int0, which are the DIJET (qcd) integrals
                // only -- h_bb/h_b/h_nob, never the _bjet ones. The shapes being mixed are
                // qcd+bjet when also_bjet. That is deliberate: qcd is the sample whose
                // flavour composition matches the data, the bjet sample is there for template
                // statistics and would bias c' downwards if it entered the ratio.
                    double eff_bkg0B = (1 - sig_fraction_true - bkg_fraction_b_true)/(1- sig_fraction_true);// c`
                    double eff_bkg1B = 1. - eff_bkg0B;// b`
            // Build effective bkg hist: with new relaitve normalization, the integral should = 1
            // Normlaize h_nob to be added effectively to the normalized bkg distribution
                TH1D* norm_h_nob = (TH1D*) h_nob->Clone("norm_h_nob"); norm_h_nob->Scale(1./norm_h_nob->Integral(1,  h_bkg_bins, "width"));
            // Variation of 0B template contribution == the LIGHT-JET MISTAG systematic.
            // The 0B template is the jets with jtNbHad == 0 that still passed the b tag,
            // i.e. mistagged light and charm jets. How many of them sit under the data is
            // taken from MC alone -- nothing in this fit constrains it, because 0B is not a
            // free component: it is folded into the background PDF at the MC-predicted
            // ratio c' = 0B/(0B+1B). So the mistag rate enters the measurement only through
            // the SHAPE of the effective background, and the systematic is to rebuild that
            // shape with the 0B admixture scaled by w_var_0B.
            //
            // w = 2 / w = 0 is a deliberately conservative +/-100% on the mistag rate: it
            // brackets the measured light-flavour mistag scale factors (typically 20-50% at
            // a tight working point) by a wide margin, and it needs no external SF input.
            double w_var_0B = 1; 
                if (ivar == VARIED0B_UP) { w_var_0B = 2.0;}
                else if (ivar == VARIED0B_DOWN) {w_var_0B = 0.0;}
                else if (ivar == NOMINAL){w_var_0B = 1;}
                cout << "Before variation: eff 1B = "<< eff_bkg1B << endl;

                // change eff 1B accordingly                 
                eff_bkg0B *= w_var_0B;
                // Guard: eff_bkg1B = 1 - eff_bkg0B, so any eff_bkg0B > 1 makes the 1B term
                // in the Add() below NEGATIVE. The integral check further down would still
                // pass -- the two weights sum to 1 by construction -- but RooHistPdf would
                // silently clip the negative bins and fit a PDF that is not the one built
                // here. It cannot happen at the current WP (c' <= 0.19, so 2c' <= 0.37), but
                // a looser b tag or a coarser dR bin can reach it, so refuse loudly instead.
                if (eff_bkg0B > 1. || eff_bkg0B < 0.) {
                    cout << Form("WARNING: (ptbin %d, deltaRbin %d) 0B variation w = %.1f gives "
                                 "eff 0B = %.4f, outside [0,1] -- clipped. The variation is "
                                 "truncated in this bin, so its systematic is a LOWER bound.",
                                 ibin_pt, ibin_dr, w_var_0B, eff_bkg0B) << endl;
                    eff_bkg0B = std::min(1., std::max(0., eff_bkg0B));
                }
                eff_bkg1B = 1. - eff_bkg0B;
                    
                    cout << "eff 0B weight = "<< w_var_0B << endl;
                    cout << "eff 0B = "<< eff_bkg0B << endl;
                    cout << "eff 1B = "<< eff_bkg1B << endl;

                // --- Build effective Bkg template 
                h_bkg->Add(h_bkg, norm_h_nob, eff_bkg1B, eff_bkg0B);//  eff_bkg1B x Normalized_1Bhist + eff_bkg0B x Normalized_0Bhist            
                // update combined normalized bkg integral and error  
                double err_int;
                double int_val = h_bkg->IntegralAndError(1,  h_bkg_bins, err_int ,"width");
                    // cout << "Effective Bkg normalized hist integral +/- uncertainity = ? (should be 1):  "<< int_val << "+/-" << err_int << endl;
                    // cout << "effective 1B and 0B fractions = (sum should be 1): "<< eff_bkg1B << ", " << eff_bkg0B << endl;
                    // avoid float point failur 
                    // Written as !(<= tol) so a NaN fails too: NaN > tol is false (2026-10-26).
                    // Returning here writes no h_sig_fraction_fit, so the unfolding stops on it.
                    if (!(std::abs(int_val - 1) <= 1e-06)) {cerr << "ERROR: Effective BKG PDF is not normalized to 1 (integral " << int_val << ")!"<< endl; return;}
                    if (!(std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1) <= 1e-06)) {cerr << "ERROR: Signal PDF is not normalized to 1!"<< endl; return;}

            // -- Draw normalized PDFs of signal and effective background as inputs: prefit
            // set PDF titles 
            h_bkg->SetTitle("1B + 0B");
            h_sig->SetTitle("2B");
            styleTemplate(h_bkg, TFColor::bkg());   // 1B + 0B summed -> neutral
            styleTemplate(h_sig, TFColor::c2B());


            // normlaized data for comparison
            /*TH1D* hnorm_data_self = (TH1D*)h_data_mb ->Clone("hnorm_data_self");
                hnorm_data_self->Scale(1./hnorm_data_self->Integral(1, mb_bins, "width"));
                hnorm_data_self->SetTitle("data self normalized");
              */
                    auto canva_pdf_norm_beforefit = new TCanvas("Pdfs_norm_beforefit","", 800, 800 );
                        canva_pdf_norm_beforefit->cd();
                        h_bkg->Draw("hist E");
                        h_sig->Draw("hist E same");
                        // hnorm_data_self->Draw("PE same");
                        canva_pdf_norm_beforefit->SetTitle("PDFs before fit");
                            gPad->Modified();   
                            gPad->Update();
                            canva_pdf_norm_beforefit->Modified();
                            canva_pdf_norm_beforefit->Update();
                            fout->cd();
                            DrawCommonTextTopRight(canva_pdf_norm_beforefit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                            // MC templates only on this canvas -> "Simulation Internal".
                            drawCMSHeader(canva_pdf_norm_beforefit, "Work In Progress");
                            canva_pdf_norm_beforefit->Write();
                            canva_pdf_norm_beforefit->Print(Form("%s/%s_PDF_norm_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
            // what about a stack of PDFs before fit 
                // these are not normalized to 1, but normlaized such that the total PDFs are 1, using their qcd fractions  
                THStack hstack_pfds_scaledtoqcd (Form("hstack_pfds_scaledtoqcd_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true fractions, before fit");
                    hstack_pfds_scaledtoqcd.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true fractions in qcd 
                        TH1D* h_sig_scaled = (TH1D*) h_sig ->Clone("h_sig_scaled"); h_sig_scaled->Scale(sig_fraction_true);
                        TH1D* h_bkg_scaled = (TH1D*) h_bkg ->Clone("h_bkg_scaled"); h_bkg_scaled->Scale(1. - sig_fraction_true);
                            hstack_pfds_scaledtoqcd.Add(h_sig_scaled); 
                            hstack_pfds_scaledtoqcd.Add(h_bkg_scaled);
                            hstack_pfds_scaledtoqcd.SetMaximum(1.2 * hstack_pfds_scaledtoqcd.GetMaximum());

                            auto canva_pdfs_scaledtoqcd = new TCanvas("canva_pdfs_scaledtoqcd","", 800, 800 );
                                canva_pdfs_scaledtoqcd->cd();
                                hstack_pfds_scaledtoqcd.Draw("Hist E");
                                // hnorm_data_self->Draw("HIST E same"); 
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd->Modified();
                                canva_pdfs_scaledtoqcd->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_scaledtoqcd, "Work In Progress");
                                canva_pdfs_scaledtoqcd->Write();
                                canva_pdfs_scaledtoqcd->Print(Form("%s/%s_pdfs_scaledtoqcdfractions_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            

                // -- Prefit PDF: Signal, BKG
                THStack hstack_pfds_scaledtoqcd_int (Form("hstack_pfds_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        TH1D* h_sig_scaledint = (TH1D*) h_sig ->Clone("h_sig_scaledint"); h_sig_scaledint->Scale(int2);
                        TH1D* h_bkg_scaledint = (TH1D*) h_bkg ->Clone("h_bkg_scaledint"); h_bkg_scaledint->Scale(int1 + int0);
                            hstack_pfds_scaledtoqcd_int.Add(h_sig_scaledint); 
                            hstack_pfds_scaledtoqcd_int.Add(h_bkg_scaledint);
                            hstack_pfds_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_scaledtoqcd_int.GetMaximum());

                            auto canva_pdfs_scaledtoqcd_int = new TCanvas("canva_pdfs_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_scaledtoqcd_int->cd();
                                hstack_pfds_scaledtoqcd_int.Draw("Hist E");
                               /* TH1D* hd_scaledtoqcdint = (TH1D*) hnorm_data_self->Clone("hd_scaledtoqcdint"); hd_scaledtoqcdint->Scale(int0+int1+int2);
                                        hd_scaledtoqcdint->SetTitle("data scaled to qcd integral");
                                        hd_scaledtoqcdint->Draw("Hist PE same");
                                */
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoqcd_int->Modified();
                                canva_pdfs_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_scaledtoqcd_int, "Work In Progress");
                                canva_pdfs_scaledtoqcd_int->Write();
                                canva_pdfs_scaledtoqcd_int->Print(Form("%s/%s_pdfs_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            
                                    // test 0B: (1B + 2B) after scaling to dijet fractions: is it like the true fractions?
                                    cout << "----------------" << endl;
                                    cout << "After combine Dijet+ bjet, and scaling to qcd integrals " << endl;
                                    cout << "Integarl of scaled PDF(2B) = " <<   h_sig_scaledint->Integral() << endl;
                                    cout << "Integral of scaled PDF(1B+0B) = " << h_bkg_scaledint->Integral() << endl;

                // -- Prefit PDF: signal, BKG (seperated 0B and 1B)
                 THStack hstack_pfds_seperated_scaledtoqcd_int (Form("hstack_pfds_seperated_scaledtoqcd_int_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to true integrals, before fit");
                    hstack_pfds_seperated_scaledtoqcd_int.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the true integrals in qcd 
                        // 1B (bjet and dijet)
                        TH1D* hnorm_sumbkg_1B_scaledint = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledint");
                                 hnorm_sumbkg_1B_scaledint->Scale(1./hnorm_sumbkg_1B_scaledint->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledint->Scale(int1);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_sig_scaledint); 
                        hstack_pfds_seperated_scaledtoqcd_int.Add(hnorm_sumbkg_1B_scaledint);
                        hstack_pfds_seperated_scaledtoqcd_int.Add(h_nob);
                        hstack_pfds_seperated_scaledtoqcd_int.SetMaximum(1.2 * hstack_pfds_seperated_scaledtoqcd_int.GetMaximum());


                            auto canva_pdfs_seperated_scaledtoqcd_int = new TCanvas("canva_pdfs_seperated_scaledtoqcd_int","", 800, 800 );
                                canva_pdfs_seperated_scaledtoqcd_int->cd();
                                hstack_pfds_seperated_scaledtoqcd_int.Draw("Hist E");
                                // hd_scaledtoqcdint->Draw("Hist PE same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_seperated_scaledtoqcd_int->Modified();
                                canva_pdfs_seperated_scaledtoqcd_int->Update();
                                fout->cd();
                                DrawCommonTextTopRight(canva_pdfs_seperated_scaledtoqcd_int, ibin_dr, ibin_pt, yBins, N_bins_dr);
                                // MC templates only on this canvas -> "Simulation Internal".
                                drawCMSHeader(canva_pdfs_seperated_scaledtoqcd_int, "Work In Progress");
                                canva_pdfs_seperated_scaledtoqcd_int->Write();
                                canva_pdfs_seperated_scaledtoqcd_int->Print(Form("%s/%s_pdfs_seperated_scaledtoqcd_int_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
            /*                     
            // Before fit: pdfs scaled to data 
             THStack hstack_pfds_scaledtoData (Form("hstack_pfds_scaledtoData_%d_%d", ibin_dr, ibin_pt),"PDFs scaled to data, before fit");
                    hstack_pfds_scaledtoData.SetTitle(Form("DeltaRBin_%d_PtBin_%d;m_{2B} [GeV];",  ibin_dr, ibin_pt));
                        // Create scaled PDFs to the data inetgral 
                        // normalized then multiplied by their fraction from qcd, then scaled to data
                        TH1D* hnorm_sumbkg_1B_scaledtoData = (TH1D*) h_sumbkg->Clone("hnorm_sumbkg_1B_scaledtoData");
                                hnorm_sumbkg_1B_scaledtoData->Scale(1./hnorm_sumbkg_1B_scaledtoData->Integral(1, mb_bins, "width"));
                                hnorm_sumbkg_1B_scaledtoData->Scale(integral_inputdata * int1/tot);
                        TH1D* h_nob_scaledtoData = (TH1D*) norm_h_nob->Clone("h_nob_scaledtoData");
                                 h_nob_scaledtoData->Scale(integral_inputdata * int0/tot);
                        TH1D* h_2B_scaledtoData = (TH1D*) h_sig->Clone("h_2B_scaledtoData");
                                h_2B_scaledtoData->Scale(integral_inputdata * int2/tot);
                        hstack_pfds_scaledtoData.Add(h_2B_scaledtoData); // 2B 
                        hstack_pfds_scaledtoData.Add(hnorm_sumbkg_1B_scaledtoData);
                        hstack_pfds_scaledtoData.Add(h_nob_scaledtoData); // 0B
                        hstack_pfds_scaledtoData.SetMaximum(1.2 * hstack_pfds_scaledtoData.GetMaximum());
                            auto canva_pdfs_scaledtoData = new TCanvas("canva_pdfs_scaledtoData","", 800, 800 );
                                canva_pdfs_scaledtoData->cd();
                                hstack_pfds_scaledtoData.Draw("Hist E");
                                h_data_mb->Draw("Hist PE same");
                                DrawCommonTextTopRight(canva_pdfs_scaledtoData, ibin_dr, ibin_pt,yBins, N_bins_dr, false);
                                TLegend* leg_before = CreateLegend(0.54, 0.6, 0.85, 0.8,
                                    {h_data_mb, h_sig, h_sumbkg, norm_h_nob},
                                    {"LPE", "LF", "LF", "LF"},
                                    {"Data", "2B", "1B", "0B"} // use default titles 
                                    );
                                    leg_before->Draw("same");
                                gPad->Modified();   
                                gPad->Update();
                                canva_pdfs_scaledtoData->Modified();
                                canva_pdfs_scaledtoData->Update();
                                fout->cd();
                                canva_pdfs_scaledtoData->Write();
                                canva_pdfs_scaledtoData->Print(Form("%s/%s_pdfs_scaledtoData_beforefit.png", sDir_canvas.Data(), sname_canvas.Data()));
*/

// return;
                                     
/*
            ///// Fitting
            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     
            // Inputs 
            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
            RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
            RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
            RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);

            // Create list of templates
            RooArgList template_list(sig_template, bkg_template, "template_list");

            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val", sig_fraction_true , 0., 1);// , // 1-bkg_fraction_b_true

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, false); // fasle: dont use extended fit: assume fraction not absolute yields  
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name            
                                                                                                                                                                          //instead of sign histogram here we would put the data histogram
            Int_t status = result->status();
            if (status != 0) {
                // One retry before the bin is given up (2026-10-26): Strategy 2 (full Hessian at
                // every step), started from the middle of [0,1] instead of the MC fraction.
                cout << Form("WARNING: (ptbin %d, bin %d) fit status %d -- retrying with Strategy 2 from a = 0.5",
                             ibin_pt, ibin_dr, status) << endl;
                sig_fraction_val.setVal(0.5);
                delete result;
                result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(2), RooFit::Minos(false));
                status = result->status();
            }
            result->Print();

            std::cout << "covariance matrix:" << std::endl;
            (result->covarianceMatrix().Print());

            //Check if it converged for a dr and jtpt bin
            if (status != 0) {
                std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                         << ", ix = " << ibin_dr 
                         << " did not converge\n\n\n\n" << std::endl;
                non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();


            //Print some check
            std::cout << "RooFit result: \n Signal PDF fraction: a = " << a << " \n its uncertainity: da = " << da << std::endl;
            /// Rescale the fit paraemters: a is for Signal + charm-light --> for now charmLight contribution is SET to ZERO I dont need it now.
            
            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a; // signal 2B
            // Total bkg: (1-a)
            p1 = (1-a)*eff_bkg1B; // 1B bkg 
            p2 = (1-a)*eff_bkg0B; // 0B bkg 
            errP0 = da;
            errP1 = da *eff_bkg1B;
            errP2 = da *eff_bkg0B;

            // std::cout << "a: 2B =" << sig_fraction_true << ", a': after fit=" << p0 << std::endl;            
            // std::cout << "a'/a  for 2B = " << p0/sig_fraction_true << std::endl;


            // std::cout << "b: 1B =" << bkg_fraction_b_true << ", b': 1B =" << p1 << std::endl;
            // std::cout << "b'/b for 1B = " << p1/bkg_fraction_b_true << std::endl;


            // std::cout << "c: 0B =" << (1 - sig_fraction_true - bkg_fraction_b_true) << ", c': 0B after fit =" << p2 << std::endl;
            // std::cout << "c'/c = " << p2/(1-sig_fraction_true - bkg_fraction_b_true) << std::endl;


            // -- updated for the new hsit binning for the true and fitresult S/B fractions (for later drawings)
            // Fit result 
            h_sig_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, p0);
            h_sig_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);

            h_bkg_fraction->SetBinContent(ibin_dr +1, ibin_pt +1, 1- p0);
            h_bkg_fraction->SetBinError(ibin_dr +1, ibin_pt +1, errP0); // errorp1 = error (1-p0)
            h_bkg_fraction_error->SetBinContent(ibin_dr +1, ibin_pt +1, errP0);


            // -- Compute siganl and bkg fractions uncertainity 
            //save the true fraction
            // -- sigma True for S or B  = (B x sigmaS + SxsigmaB)/(S+B)²
            // 2026-10-26: in quadrature (S and B are independent samples); was the linear sum
            // ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/(S+B)^2.
            double err_true_frac = std::hypot((int0 + int1) * True_sig_err, int2 * True_bkg_err)/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B
            h_sig_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, sig_fraction_true);
            h_sig_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            
            h_bkg_frac_true->SetBinContent(ibin_dr +1, ibin_pt +1, bkg_fraction_true);
            h_bkg_frac_true->SetBinError(ibin_dr +1, ibin_pt +1, err_true_frac);
            h_bkg_frac_true_error->SetBinContent(ibin_dr +1, ibin_pt +1, err_true_frac);
            

            // -- After fits: save mass distribution re and post-fit 
            TH1D *h_sig_fit = (TH1D*) h_sig->Clone(Form("h_sig_fit_%d_%d", ibin_dr, ibin_pt));
                h_sig_fit->Scale(p0 * integral_inputdata);
                h_sig_fit->SetTitle("2B");
                h_sig_fit->GetYaxis()->SetTitle("Counts/[GeV^{2}]");

                styleTemplate(h_sig_fit, TFColor::c2B());   // 2B -> red
                cout << "Signal integral after fit = " << h_sig_fit->Integral(1, mb_bins, "width") << endl;

            TH1D *h_bkg_fit = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_%d_%d", ibin_dr, ibin_pt)); // total bkg
                h_bkg_fit->Scale(integral_inputdata * (1.0 - p0));
                h_bkg_fit->SetTitle("1B+ 0B");
                styleTemplate(h_bkg_fit, TFColor::bkg());   // 1B + 0B summed -> neutral

            // 1B and 0B slices with their OWN shapes (2026-10-26): both were clones of h_bkg, the
            // 1B+0B mixture, so they were drawn and written with one and the same shape. They
            // still add up to h_bkg_fit: p1 + p2 = 1-a, mixed exactly as h_bkg was.
            TH1D *h_bkg_fit_1b = (TH1D*) norm_h_1b->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) norm_h_nob->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
                h_bkg_fit_nob->Scale(integral_inputdata * p2);
                h_bkg_fit_nob->SetTitle("0B");
                styleTemplate(h_bkg_fit_nob, TFColor::c0B());  // 0B -> green

                // And the post-fit template (fitted sig + bkg)
            TH1D* h_total_fit = (TH1D*) h_sig_fit->Clone(Form("h_total_fit_%d_%d", ibin_dr, ibin_pt));
                  h_total_fit->Add(h_bkg_fit);
                  h_total_fit->SetTitle("Total fit");
                    // The sum of the three slices, not a category of its own: neutral dark
                    // outline and no fill, so it cannot be misread as one of them.
                    h_total_fit->SetFillColor(0);
                    h_total_fit->SetFillStyle(0);
                    h_total_fit->SetMarkerColor(TFColor::total());
                    h_total_fit->SetLineColor(TFColor::total());
                    h_total_fit->SetMarkerStyle(1);
                    h_total_fit->SetLineWidth(2);


            // Save Signal and bkg mass distributions after the fit 
            fout->cd();
            h_sig_fit->Write();
            h_bkg_fit->Write();
            h_total_fit->Write();

            /// -- Draw useful canvas: Distiburions of Sig, Bkg, MC before and after fit 
            TString sname_canvas_afterfit = sname_canvas + "_afterfit";
            auto canva_afterfit = new TCanvas(Form("ALLHist_%s", sname_canvas_afterfit.Data()) ,Form("Templaets pre and post-fit, %s", sname_canvas.Data()), 800, 800 );
                canva_afterfit->cd();
                h_data_mb->SetTitle("Data");
                h_data_mb->SetLineWidth(2);
                h_bb_bjet->SetLineWidth(2);
                h_sig_fit->SetLineWidth(2);
                h_b_bjet->SetLineStyle(9);
                h_bkg_fit->SetLineWidth(2);
                if (h_data_mb->GetMaximum() < h_total_fit->GetMaximum() ) { h_data_mb->SetMaximum(1.3 * h_total_fit->GetMaximum());}
                else { h_data_mb->SetMaximum(1.3 * h_data_mb->GetMaximum()); }

                h_data_mb->Draw("P E");
                h_total_fit->Draw("P E SAME");
                h_bkg_fit_1b->Draw("HIST E same");
                h_bkg_fit_nob->Draw("HIST E SAME");
                h_sig_fit->Draw("HIST E SAME");
                canva_afterfit->SetTitle("");
                gPad->Modified();   
                gPad->Update();
                canva_afterfit->Modified();
                canva_afterfit->Update();
               

        // -- After fit: stacked: seperated contibutions 
            THStack hstack_afterfit("hstack_afterfit","Mass stacked histogram");
                hstack_afterfit.SetTitle(";m_{2B} [GeV];");
                hstack_afterfit.Add(h_sig_fit);
                hstack_afterfit.Add(h_bkg_fit_1b);
                hstack_afterfit.Add(h_bkg_fit_nob);
                hstack_afterfit.SetMaximum(1.2 * hstack_afterfit.GetMaximum());

            auto canva_stack_afterfit = new TCanvas(Form("All_templates_Data_stacked_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_afterfit->cd();
                hstack_afterfit.Draw("hist E"); 
                h_data_mb->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_afterfit->Modified();
                canva_stack_afterfit->Update();

            // normalized version: data and templates afterfit: Data Vs. pdf signal Vs. bkg 
            THStack hstack_norm_afterfit("hstack_norm_afterfit","Mass stacked histogram");
                hstack_norm_afterfit.SetTitle(";m_{2B} [GeV];");    
                TH1D* h_sig_fit_normstack = (TH1D*)  h_sig_fit->Clone("h_sig_fit_normstack");  h_sig_fit_normstack ->Scale(1./integral_inputdata);
                TH1D* h_bkg_fit_normstack = (TH1D*)  h_bkg_fit->Clone("h_bkg_fit_normstack");  h_bkg_fit_normstack ->Scale(1./integral_inputdata);
                hstack_norm_afterfit.Add(h_sig_fit_normstack);
                hstack_norm_afterfit.Add(h_bkg_fit_normstack);
                hstack_norm_afterfit.SetMaximum(1.3 * hnorm_data_self->GetMaximum());

                auto canva_stack_norm_afterfit = new TCanvas(Form("PDFs_Data_stacked_norm_%s", sname_canvas_afterfit.Data()),Form(""), 800, 800 );
                canva_stack_norm_afterfit->cd();
                hstack_norm_afterfit.Draw("hist E"); 
                hnorm_data_self->Draw("PE same"); 
                gPad->Modified();   
                gPad->Update();
                canva_stack_norm_afterfit->Modified();
                canva_stack_norm_afterfit->Update();

            //  -- Build legend and write plots     
                fout->cd();
                DrawCommonTextTopRight(canva_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_afterfit->Write();
                canva_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_afterfit,ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_afterfit->Write();
                canva_stack_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_afterfit->GetName()));


                DrawCommonTextTopRight(canva_stack_norm_afterfit, ibin_dr, ibin_pt, yBins, N_bins_dr);
                canva_stack_norm_afterfit->Write();
                canva_stack_norm_afterfit->Print(Form("%s/%s.png", sDir_canvas.Data(), canva_stack_norm_afterfit->GetName()));

            //-- ratio plot: data /total fit
                TCanvas* c = new TCanvas(Form("RatioPlot_%s", sname_canvas_afterfit.Data()), "", 900, 1100); // 800, 1100
                    TPad* pad1 = new TPad("pad1","",0,0.2,1,1);
                    TPad* pad2 = new TPad("pad2","",0,0,1,0.24);
                    // pad1->SetBottomMargin(0.13);
                    pad1->SetBottomMargin(0.06);
                    pad1->SetLeftMargin(0.18); // for y axis title space 
                    pad2->SetTopMargin(0.02);     // bottom pad (very small)
                    pad2->SetBottomMargin(0.40);  // keep space for x-axis labels
                    pad2->SetLeftMargin(0.18);  

                    pad1->Draw();
                    pad2->Draw(); 
                    pad1->cd();
                    // Draw frame to control y axis name: frame needed to control y axis name 
                    TH1F *frame = pad1->DrawFrame(hstack_afterfit.GetXaxis()->GetXmin(), 0, hstack_afterfit.GetXaxis()->GetXmax(), hstack_afterfit.GetMaximum()*1.3);
                        frame->GetYaxis()->SetTitle("Counts/[GeV]");
                        hstack_afterfit.Draw("hist E same"); 
                        h_data_mb->Draw("PE same"); 
                        // h_total_fit->Draw("P E SAME");
                        // Add (dr, pt) bins legend 
                        DrawCommonTextTopRight(pad1, ibin_dr, ibin_pt, yBins,N_bins_dr ,false); // without default bildlegend of other objects
                        // use new Legend for enties (withut hframe)
                        TLegend* leg = CreateLegend(0.54, 0.6, 0.85, 0.8,
                            {h_data_mb, h_sig_fit, h_bkg_fit_1b, h_bkg_fit_nob},
                            {"LPE", "LF", "LF", "LF"},
                            {"Data", "", "", ""} // use default titles 
                        );
                        leg->Draw("same");
                        pad1->Modified(); // force refresh 
                        pad1->Update();
                    pad2->cd(); 
                    AddRatioPlot(h_data_mb, h_total_fit);
                    pad2->SetTickx(1);// → draws ticks on both bottom and top
                        fout->cd();
                        c ->Write();
                        c ->Print(Form("%s/%s.png", sDir_canvas.Data(), c->GetName()));


                cout << "---------------------\n\n\n" << endl; 
*/
            } // loop over deltaR bins 
    }

/*
    // Save histograms
    // TH3D 
    // for (auto h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
    for (auto h : {h3D_data, h3D_bb_bjet, h3D_b_bjet}) {h->Write();}
        

        
    // TH2D     
    for (auto h : {h_sig_fraction, h_sig_fraction_error,
                   h_bkg_fraction, h_bkg_fraction_error,
                   h_sig_frac_true, h_sig_frac_true_error,
                   h_bkg_frac_true, h_bkg_frac_true_error 
                   }) {
                    h->Write();
    }
*/

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
    fout->Close();

}


// The MC templates for one sample and generator: the MCGEN file that
// create_files_for_template_fit.cpp wrote for it. Paths are listed rather than built,
// because the Pythia8 ones are Afnan's merged production (btagWP712, _80_9999_2) and the
// Herwig ones come from run_agg_ntuple_chunks.sh here (btagWP0712, _upartv2).
// track_eff_unc = true takes the templates from the 3%-track-drop production instead
// (TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh, then the per-block MCGEN files hadd'ed
// to the agg_ntuple_chunks top level). Only Pythia8 has that production today; asking for
// it with generator "herwig" returns "" and the caller stops with a message.
TString mcgenTemplates(const TString &sample, const TString &generator,
                       bool track_eff_unc = false, bool eec_weight_off = false,
                       const TString &observable = "dr")
{
    // The B (momentum-balance) histograms exist ONLY in the "_upartv2_B" productions --
    // the ones run from this working copy after the 2026-09-22 z -> B rename. Afnan's merged
    // Pythia files below contain dR only, and the pre-rename productions name the balance
    // axis "_z", so a B fit cannot use either and must not silently fall back to them.
    //
    // Every B production is Zoe's own, all under ONE naming scheme, so the path is built
    // rather than listed: the generator picks the directory, and the tag is the same one
    // apply_unfolding_2d.C reads its per-block files with (mcVarTag + obsProdTag), in the
    // order the macro and the run script write it:
    //     <sample>_fMCGEN [_trkdrop030] [_noeecw] _upartv2_B .root
    // That covers the three fits the B band needs -- Pythia nominal, Pythia 3% track drop
    // ("Tracking efficiency") and Herwig nominal ("MC template modeling") -- each from its
    // per-block files hadd'ed to the agg_ntuple_chunks top level. A combination nobody
    // produced resolves to a file that does not exist, and the driver stops on it.
    // The rho era (result_paths.h, 2026-09-28) sends dr through this generic branch too: every
    // template comes from Zoe's own _rho production, and the Afnan / pre-_rho fallbacks below
    // point at inputs that no longer exist upstream.
    if (observable != "dr" || rhoInputs()) {
        const TString base = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/";
        const bool herwig  = (generator == "herwig");
        TString subdir;
        if (sample == "qcd")  subdir = herwig ? "QCDHerwig"  : "QCD";
        if (sample == "bjet") subdir = herwig ? "bJetHerwig" : "bJet";
        if (subdir.Length() == 0) return "";
        // trkTag / eecWeightTag / obsProdTag are result_paths.h's, so the literals are shared
        // with apply_unfolding_2d.C and cannot drift from what it reads.
        // svsfTag(): the SV b scale-factor production when SV_SF is set (2026-10-12).
        const TString tag = trkTag(track_eff_unc) + svsfTag() + trackBdtSfTag()
                          + eecWeightTag(eec_weight_off)
                          + obsProdTag(observable);
        return base + subdir + "/agg_ntuple_chunks/"
               "Run3_btagWP0712_template_for_fit_histos_3D_" + sample + "_fMCGEN" + tag + ".root";
    }

    // The YIELD production (EEC weight off). Not a variation of the EEC fit -- the same fit
    // performed on unweighted templates, so the signal fraction it returns is the one a
    // yield measurement needs. Zoe's own blocks, hadd'ed to the agg_ntuple_chunks top level:
    // Afnan's merged Pythia production is EEC-weighted and has no unweighted counterpart,
    // so there is deliberately nothing to fall back on here.
    if (eec_weight_off) {
        if (generator != "pythia") return "";  // no Herwig yield production
        // The tracking variation of the yield run: 3% of reco tracks dropped AND the EEC
        // weight off. Tagged "_trkdrop030_noeecw" -- the macro appends TrkEffSyst::tag()
        // then EecWeight::tag(), in that order, and variationTag() in result_paths.h builds
        // the reader side the same way round.
        const TString tag = track_eff_unc ? "_trkdrop030_noeecw" : "_noeecw";
        if (sample == "qcd")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_qcd_fMCGEN"
                 + tag + "_upartv2.root";
        if (sample == "bjet")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_bjet_fMCGEN"
                 + tag + "_upartv2.root";
        return "";
    }
    if (track_eff_unc) {
        if (generator != "pythia") return "";   // no varied Herwig production
        if (sample == "qcd")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_qcd_fMCGEN_trkdrop030_upartv2.root";
        if (sample == "bjet")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_bjet_fMCGEN_trkdrop030_upartv2.root";
        return "";
    }
    if (generator == "pythia") {
        if (sample == "qcd")
            return "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2MCGEN_merged.root";
        if (sample == "bjet")
            return "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2MCGEN_merged.root";
    }
    if (generator == "herwig") {
        if (sample == "qcd")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCDHerwig/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_qcd_fMCGEN_upartv2.root";
        if (sample == "bjet")
            return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/bJetHerwig/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_bjet_fMCGEN_upartv2.root";
    }
    return "";
}

// ============================================================================
// do_template_fit_chi2scan() -- the template fit done the way CMS-AN-2025/114 does its
// prompt / non-prompt fit (Sec. 9.6 and App. B; Lorenzo Frosina's nonPrompt_Subtraction_chi2.C).
// 2026-10-27, Zoe: "do the same, but without the 0B for now". Replaces the RooFit fit of
// do_template_fit_combined() in template_fit(); that function is kept, not called.
//
// Our 2B / 1B play the AN's prompt / non-prompt. Per (pT, observable) bin:
//   model(m2B) = k * [ f * T_2B(m2B) + (1-f) * T_1B(m2B) ],   T = unit-area MC templates
//   (the 2B fraction is called k on plots and printouts -- Zoe, 2026-10-27; f / f_* in the code)
//   chi2(k, f) = sum_bins (data - model)^2 / (sigma_data^2 + sigma_model^2)
// sigma_model is the template MC error propagated through the model, so -- unlike the RooFit
// fit, whose likelihood ignored it -- the template statistics are in the fit. (k, f) are
// SCANNED on a grid, f in [0,1], and the pair with the lowest chi2 is taken. The error on f is
// the AN's chi2 profile at FIXED k*: the f where chi2 = chi2_min + 1 on either side, clipped to
// [0,1] -- asymmetric, and written as h_sig_fraction_fit_up / _down next to the nominal.
//
// Nominal (Zoe, 2026-10-27, final state of the day):
//   - background = 1B + 0B at the MC ratio (kWith0B), 0B from qcd only;
//   - templates REWEIGHTED IN JET pT: w(pT) = data/MC (2B+1B+0B) jet-pT spectrum in 5 GeV bins,
//     normalised within each coarse fit pT bin (block after the qcd+bjet combination);
//   - (kM2BReweight, now OFF): reweighting in m2B as in the AN, w(m2B) = data / MC sum on every
//     template. chi2 = 0 by construction; dropped by Zoe at the end of 2026-10-27.
// Kept for comparison: the fit on the plain templates (h_sig_fraction_fit_noRW) and on the jet-pT
// reweighted ones only (h_sig_fraction_fit_jetptRW). m2B: all 10 bins, 0-10 GeV.
// Needs the _ptfine histograms, i.e. the "...ptw" step-1 production (template_fit() reads it).
//
// Differences from the AN, all deliberate:
//   - no 0B: the background is the 1B template alone (Zoe, "for now"). A 0B admixture would be
//     a third template; nothing below assumes there are only two except T_1B.
//   - templates WITH the MC event weights (Zoe, "use the event weights for now"). The AN fits
//     unweighted templates (its App. H); the weighted-vs-unweighted template plots of
//     Lorenzo's macro need an unweighted production we do not have, so they are not made.
//   - the grid: k in [0, 2 x data integral] in 500 steps (the AN: [0, 5 x data maximum]) and
//     f in 1000 steps of 0.001 (the AN: 100 of 0.01). k* ~ data integral since the templates
//     are unit-area, and our fit errors on f are a few 0.01, so the AN's grid would quantise
//     the result at the size of its own error.
//   - chi2/ndf is quoted with ndf = bins used - 2 (the AN macro counts the bins only).
//
// Output (nominal_<fout_name> in sDirname): h_sig_fraction_fit (+_up/_down, error = larger
// side), h_sig_frac_true_fitbins (MC 2B/(2B+1B)), h_k_fit, h_chi2ndf_fit, and per bin the data
// slice, the two templates, the chi2 map and its profile. PNGs into sDirname_www.
// ============================================================================
// zeroBScale (2026-10-28): the light-jet mistag variation -- the 0B template (and its fine-pT copy,
// so the jet-pT weights see the same MC) scaled by 2 (var0B_2) or 0 (var0B_0) before anything uses
// it, as the RooFit fit's 0B x2 / x0. 1 = nominal.
void do_template_fit_chi2scan(const TString &data_name, const TString &templates,
                              const TString &templates_bjet, bool also_bjet, const TString &fout_name,
                              double zeroBScale = 1.)
{
    const ObsDef &obs = gFitObs();
    const TString sDir_canvas = sDirname_www;
    gSystem->mkdir(sDir_canvas, kTRUE);

    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    setBluePurpleRedPalette();   // Zoe: the usual red/blue scale on the chi2 and fraction maps (was ROOT's default)

    // ---- inputs: same files and names as do_template_fit_combined(), minus the 0B ----
    auto get3D = [](const TString &file, const TString &name, const char *clone) -> TH3D* {
        TFile *f = TFile::Open(file, "read");
        if (!f || f->IsZombie()) { Error("chi2scan", "cannot open %s", file.Data()); return nullptr; }
        TH3D *h = dynamic_cast<TH3D*>(f->Get(name));
        if (!h) { Error("chi2scan", "no %s in %s", name.Data(), file.Data()); f->Close(); return nullptr; }
        h = (TH3D*) h->Clone(clone); h->SetDirectory(nullptr); f->Close();
        return h;
    };
    TH3D *h3D_data = get3D(data_name, obs.n("h3D_data"), "h3D_data");
    TH3D *h3D_bb   = get3D(templates, obs.n("h3D_bb"), "h3D_bb");
    TH3D *h3D_b    = get3D(templates, obs.n("h3D_b"),  "h3D_b");
    // 0B (2026-10-27, Zoe "and if you add the 0B?"): with kWith0B the background is 1B + 0B at the MC
    // ratio, as in the RooFit fit -- 0B from qcd only (the b-enriched sample's light jets are
    // filter-biased), at qcd normalisation like the combined 1B. false = 2B vs 1B only.
    const bool kWith0B = true;
    // m2B reweighting of the MC sum (AN Sec. 9.6): OFF since the end of 2026-10-27 (Zoe: "you don't need
    // to do the mass reweighting"). Nominal = 0B + jet-pT reweighting only. An injection test showed the
    // m2B reweighting recovers only ~10-65% of a true 2B-fraction difference from MC (template overlap).
    const bool kM2BReweight = false;
    TH3D *h3D_0b   = get3D(templates, obs.n("h3D_0b"), "h3D_0b");
    if (!h3D_0b || !CheckInputBinning(h3D_0b)) return;
    if (zeroBScale != 1.) { h3D_0b->Scale(zeroBScale); printf("0B template scaled by %g (mistag variation)\n", zeroBScale); }
    TH3D *h3D_bb_bjet = also_bjet ? get3D(templates_bjet, obs.n("h3D_bb"), "h3D_bb_bjet") : nullptr;
    TH3D *h3D_b_bjet  = also_bjet ? get3D(templates_bjet, obs.n("h3D_b"),  "h3D_b_bjet")  : nullptr;
    if (!h3D_data || !h3D_bb || !h3D_b || (also_bjet && (!h3D_bb_bjet || !h3D_b_bjet))) return;
    for (TH3D *h : {h3D_data, h3D_bb, h3D_b, h3D_bb_bjet, h3D_b_bjet})
        if (h && !CheckInputBinning(h)) return;

    // qcd + bjet: the weighted combination of do_template_fit_combined() -- bjet scaled to qcd's
    // 2B (s2b) and 1B (s1b) yields, then the two mixed with Neff weights.
    // The factors are kept (c_q, c_bb, c_b) to combine the fine-pT templates the same way below.
    double c_q = 1., c_bb = 0., c_b = 0.;
    if (also_bjet) {
        auto sums = [](const TH3D *h, double &S, double &Q) {
            S = 0.; Q = 0.;
            for (int i = 0; i < h->GetNcells(); ++i) { S += h->GetBinContent(i); Q += std::pow(h->GetBinError(i), 2); }
        };
        double Sq, Qq, Sb, Qb, S1q, Q1q, S1b, Q1b;
        sums(h3D_bb, Sq, Qq); sums(h3D_bb_bjet, Sb, Qb); sums(h3D_b, S1q, Q1q); sums(h3D_b_bjet, S1b, Q1b);
        if (!(Sq > 0. && Sb > 0. && Qq > 0. && Qb > 0. && S1q > 0. && S1b > 0.)) {
            Error("chi2scan", "empty 2B or 1B template, not combining"); return;
        }
        const double s2b = Sq / Sb, s1b = S1q / S1b, nq = Sq * Sq / Qq, nb = Sb * Sb / Qb;
        const double wb = nb / (nq + nb), wq = 1. - wb;
        printf("qcd+bjet combination: bjet scale 2B %.5f, 1B %.5f, w_qcd %.4f, w_bjet %.4f\n", s2b, s1b, wq, wb);
        h3D_bb->Scale(wq); h3D_b->Scale(wq);
        h3D_bb->Add(h3D_bb_bjet, wb * s2b); h3D_b->Add(h3D_b_bjet, wb * s1b);
        c_q = wq; c_bb = wb * s2b; c_b = wb * s1b;
    }

    // ---- jet-pT reweighting of the MC (2026-10-27, Zoe) ----
    // The _ptfine histograms (step 1, 5 GeV jet-pT bins 80-400 GeV since 2026-10-28, 80-300 before) give the jet-pT spectrum of data and
    // of MC = 2B + 1B (qcd+bjet, combined as above) + 0B (qcd). Within each coarse pT bin of the fit
    // both are normalised, and w(pT) = data/MC per 5 GeV bin; so the weights only reshape the spectrum
    // inside a coarse bin, never move yield between bins. Summed over m2B and the observable (one curve
    // per coarse bin: stat-limited otherwise), the 0B included because data contain it. The 2B and 1B
    // templates are then rebuilt per coarse bin as sum_i w_i T(m2B, x, pT_i). pT is not the fit
    // variable, so unlike the m2B reweighting this cannot force f onto the MC.
    TH3D *h3D_data_pt = get3D(data_name, obs.n("h3D_data_ptfine"), "h3D_data_ptfine");
    TH3D *hF_bb = get3D(templates, obs.n("h3D_bb_ptfine"), "h3D_bb_ptfine");
    TH3D *hF_b  = get3D(templates, obs.n("h3D_b_ptfine"),  "h3D_b_ptfine");
    TH3D *hF_0b = get3D(templates, obs.n("h3D_0b_ptfine"), "h3D_0b_ptfine");
    TH3D *hF_bb_bjet = also_bjet ? get3D(templates_bjet, obs.n("h3D_bb_ptfine"), "h3D_bb_ptfine_bjet") : nullptr;
    TH3D *hF_b_bjet  = also_bjet ? get3D(templates_bjet, obs.n("h3D_b_ptfine"),  "h3D_b_ptfine_bjet")  : nullptr;
    if (!h3D_data_pt || !hF_bb || !hF_b || !hF_0b || (also_bjet && (!hF_bb_bjet || !hF_b_bjet))) {
        Error("chi2scan", "no _ptfine histograms -- they exist only in the ...ptw productions (step 1 of 2026-10-27)");
        return;
    }
    // 2026-10-28: data and MC _ptfine must come from the SAME step-1 binning (80-300 before 10-28,
    // 80-400 after). The first 80-400 rerun fitted new MC against old data: the loops below index
    // both with the MC's bins, so nothing failed -- the weights were just wrong. Refuse instead.
    for (TH3D *h : {hF_b, hF_0b, hF_bb_bjet, hF_b_bjet, h3D_data_pt}) {
        if (!h) continue;
        const TAxis *a = h->GetZaxis(), *r = hF_bb->GetZaxis();
        if (a->GetNbins() != r->GetNbins() || std::abs(a->GetXmin() - r->GetXmin()) > 1e-6 || std::abs(a->GetXmax() - r->GetXmax()) > 1e-6) {
            Error("chi2scan", "_ptfine binning differs: %s has %d bins %g-%g, 2B MC %d bins %g-%g -- "
                  "data and MC productions are from different step-1 versions; re-make the older one",
                  h->GetName(), a->GetNbins(), a->GetXmin(), a->GetXmax(), r->GetNbins(), r->GetXmin(), r->GetXmax());
            return;
        }
    }
    if (zeroBScale != 1.) hF_0b->Scale(zeroBScale);
    hF_bb->Scale(c_q); hF_b->Scale(c_q);
    if (also_bjet) { hF_bb->Add(hF_bb_bjet, c_bb); hF_b->Add(hF_b_bjet, c_b); }
    // const int nFine = hF_bb->GetNbinsZ();   // until 2026-10-28 evening: all fine bins (80-400 GeV)
    // 2026-10-28 (Zoe: "do the reweighting up til 300 GeV and throw everything above"): only the fine
    // bins below kJetPtMax are used -- for the weights AND for every fit input rebuilt from them below,
    // so jets above 300 GeV are out of the fit, data and MC alike. Step 1 still fills 80-400.
    const double kJetPtMax = 300.;
    const int nFine = hF_bb->GetZaxis()->FindBin(kJetPtMax - 1e-6);
    auto coarseOf = [&](double pt) { for (int c = 1; c < jtpt_bins; ++c) if (pt < jtpt_binsVector[c]) return c; return (int) jtpt_bins; };

    TH1D *h_ptspec_data = h3D_data_pt->ProjectionZ("h_jetpt_data", 1, h3D_data_pt->GetNbinsX(), 1, h3D_data_pt->GetNbinsY());
    TH3D *hF_mc = (TH3D*) hF_bb->Clone("h3D_mc_ptfine"); hF_mc->Add(hF_b); hF_mc->Add(hF_0b);
    TH1D *h_ptspec_mc = hF_mc->ProjectionZ("h_jetpt_mc", 1, hF_mc->GetNbinsX(), 1, hF_mc->GetNbinsY());
    TH1D *h_w_jetpt = (TH1D*) h_ptspec_data->Clone("h_w_jetpt"); h_w_jetpt->Reset();
    for (TH1D *h : {h_ptspec_data, h_ptspec_mc, h_w_jetpt}) h->SetDirectory(nullptr);
    // std::vector<double> Dc(jtpt_bins + 1, 0.), Mc(jtpt_bins + 1, 0.);   // until 2026-10-28: per coarse pT bin
    // for (int i = 1; i <= nFine; ++i) {
    //     const int c = coarseOf(h_w_jetpt->GetXaxis()->GetBinCenter(i));
    //     Dc[c] += h_ptspec_data->GetBinContent(i); Mc[c] += h_ptspec_mc->GetBinContent(i);
    // }
    // 2026-10-28 (Zoe: "you should integrate it only in the full range"): data and MC normalised ONCE,
    // over the whole fine range used (80-300 GeV since the evening, 80-400 before). The templates are
    // unit-area per fit bin, so f does not move.
    const double Dall = h_ptspec_data->Integral(1, nFine), Mall = h_ptspec_mc->Integral(1, nFine);
    for (int i = 1; i <= nFine; ++i) {
        // const int c = coarseOf(h_w_jetpt->GetXaxis()->GetBinCenter(i));
        // const double d = h_ptspec_data->GetBinContent(i) / Dc[c], ed = h_ptspec_data->GetBinError(i) / Dc[c];
        // const double m = h_ptspec_mc->GetBinContent(i) / Mc[c],   em = h_ptspec_mc->GetBinError(i) / Mc[c];
        const double d = h_ptspec_data->GetBinContent(i) / Dall, ed = h_ptspec_data->GetBinError(i) / Dall;
        const double m = h_ptspec_mc->GetBinContent(i) / Mall,   em = h_ptspec_mc->GetBinError(i) / Mall;
        // An MC-empty bin cannot be reweighted: w = 1 there (data in it, if any, stays unmodelled).
        h_w_jetpt->SetBinContent(i, m > 0. ? d / m : 1.);
        h_w_jetpt->SetBinError  (i, (m > 0. && d > 0.) ? (d / m) * std::hypot(ed / d, em / m) : 0.);
    }
    // The reweighted coarse templates, with the coarse axes of h3D_bb / h3D_b.
    // unit = true: w = 1, i.e. the coarse histogram re-made from the fine one (jets < 400 GeV only).
    auto rebuild = [&](const TH3D *fine, const TH3D *coarse, const char *name, bool unit = false) {
        TH3D *h = (TH3D*) coarse->Clone(name); h->Reset(); h->SetDirectory(nullptr);
        for (int im = 1; im <= fine->GetNbinsX(); ++im)
            for (int ix = 1; ix <= fine->GetNbinsY(); ++ix)
                for (int i = 1; i <= nFine; ++i) {
                    const int c = coarseOf(fine->GetZaxis()->GetBinCenter(i));
                    // const double w = h_w_jetpt->GetBinContent(i);   // until 2026-10-28
                    const double w = unit ? 1. : h_w_jetpt->GetBinContent(i);
                    const int b = h->GetBin(im, ix, c);
                    h->AddBinContent(b, w * fine->GetBinContent(im, ix, i));
                    h->SetBinError(b, std::hypot(h->GetBinError(b), w * fine->GetBinError(im, ix, i)));
                }
        return h;
    };
    TH3D *h3D_bb_rw = rebuild(hF_bb, h3D_bb, "h3D_bb_jetptrw");
    TH3D *h3D_b_rw  = rebuild(hF_b,  h3D_b,  "h3D_b_jetptrw");
    TH3D *h3D_0b_rw = rebuild(hF_0b, h3D_0b, "h3D_0b_jetptrw");
    // the fine histograms must add up to the coarse ones the fit used before
    // printf("jet-pT reweighting: fine/coarse 2B %.6f, 1B %.6f (should be 1); weights 80-100: %.3f-%.3f, >=100: %.3f-%.3f\n",
    //        hF_bb->Integral() / h3D_bb->Integral(), hF_b->Integral() / h3D_b->Integral(),
    //        h_w_jetpt->GetBinContent(1), h_w_jetpt->GetBinContent(4), h_w_jetpt->GetBinContent(5), h_w_jetpt->GetBinContent(nFine));
    // 2026-10-28 (Zoe: "dont keep jets with more than 400 GeV"): the _ptfine histograms stop at 400 GeV
    // (step 1 does not fill them above), so the data and the un-reweighted templates the fit uses are
    // re-made from them too -- every input of the fit without the jets >= 400 GeV. fine/coarse = the
    // share kept (< 1 by the jets above 400 GeV only).
    // printf("jet-pT reweighting: kept below 400 GeV (fine/coarse) ...   // until 2026-10-28 evening
    auto kept = [&](const TH3D *f) { return f->Integral(0, -1, 0, -1, 1, nFine); };
    printf("jet-pT reweighting: kept below %.0f GeV (fine/coarse) data %.6f, 2B %.6f, 1B %.6f, 0B %.6f; weights %.3f (first 5 GeV bin) ... %.3f (%.0f-%.0f GeV)\n",
           kJetPtMax, kept(h3D_data_pt) / h3D_data->Integral(), kept(hF_bb) / h3D_bb->Integral(),
           kept(hF_b) / h3D_b->Integral(), kept(hF_0b) / h3D_0b->Integral(),
           h_w_jetpt->GetBinContent(1), h_w_jetpt->GetBinContent(nFine),
           h_w_jetpt->GetXaxis()->GetBinLowEdge(nFine), h_w_jetpt->GetXaxis()->GetBinUpEdge(nFine));
    // h3D_data = rebuild(h3D_data_pt, h3D_data, "h3D_data_below400", true);   // first version: renamed the
    // h3D_bb   = rebuild(hF_bb, h3D_bb, "h3D_bb_below400", true);               // written h3D_data, which
    // h3D_b    = rebuild(hF_b,  h3D_b,  "h3D_b_below400",  true);               // plot_svsf_comparison.C reads
    // h3D_0b   = rebuild(hF_0b, h3D_0b, "h3D_0b_below400", true);
    h3D_data = rebuild(h3D_data_pt, h3D_data, h3D_data->GetName(), true);   // same names as before the cut
    h3D_bb   = rebuild(hF_bb, h3D_bb, h3D_bb->GetName(), true);
    h3D_b    = rebuild(hF_b,  h3D_b,  h3D_b->GetName(),  true);
    h3D_0b   = rebuild(hF_0b, h3D_0b, h3D_0b->GetName(), true);

    // the jet-pT spectra (unit area over the range used, 80-300 GeV) and the weights
    {
        TCanvas cv("c_jetpt_rw", "", 1100, 1100);
        TPad *p1 = new TPad("pjpt1", "", 0, 0.3, 1, 1), *p2 = new TPad("pjpt2", "", 0, 0, 1, 0.34);
        p1->SetBottomMargin(0.03); p1->SetLeftMargin(0.15); p1->SetTopMargin(0.10); p1->SetLogy();
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.30); p2->SetLeftMargin(0.15);
        p1->SetTicks(1, 1); p2->SetTicks(1, 1);
        p1->Draw(); p2->Draw(); p1->cd();
        // TH1D *dn = (TH1D*) h_ptspec_data->Clone("jpt_dn"); dn->Scale(1. / dn->Integral());   // until 2026-10-28 evening
        // TH1D *mn = (TH1D*) h_ptspec_mc->Clone("jpt_mn");   mn->Scale(1. / mn->Integral());
        TH1D *dn = (TH1D*) h_ptspec_data->Clone("jpt_dn"); dn->Scale(1. / dn->Integral(1, nFine));
        TH1D *mn = (TH1D*) h_ptspec_mc->Clone("jpt_mn");   mn->Scale(1. / mn->Integral(1, nFine));
        dn->GetXaxis()->SetRange(1, nFine); mn->GetXaxis()->SetRange(1, nFine);
        styleData(dn); dn->SetTitle("Data");
        mn->SetLineColor(kBlue + 1); mn->SetLineWidth(2); mn->SetFillStyle(0); mn->SetTitle("MC (2B+1B+0B)");
        mn->GetYaxis()->SetTitle("1/N dN/dp_{T}^{jet}"); mn->GetXaxis()->SetLabelSize(0);
        mn->Draw("hist E"); dn->Draw("PE same");
        TLegend *lg = CreateLegend(0.60, 0.70, 0.86, 0.84, {dn, mn}, {"LPE", "L"});
        lg->Draw();
        drawCMSHeader(p1);
        p2->cd();
        TH1D *wd = (TH1D*) h_w_jetpt->Clone("jpt_w");
        wd->SetMarkerStyle(20); wd->SetLineColor(kBlack); wd->SetMinimum(0.); wd->SetMaximum(2.5);
        wd->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]"); wd->GetYaxis()->SetTitle("w = Data/MC");
        wd->GetYaxis()->CenterTitle(true); wd->GetYaxis()->SetNdivisions(505);
        wd->GetYaxis()->SetTitleSize(0.08); wd->GetYaxis()->SetLabelSize(0.07); wd->GetYaxis()->SetTitleOffset(0.7);
        wd->GetXaxis()->SetTitleSize(0.09); wd->GetXaxis()->SetLabelSize(0.08);
        wd->GetXaxis()->SetRange(1, nFine);   // 2026-10-28 evening: only the bins used (80-300 GeV)
        wd->Draw("EP");
        // TLine one(wd->GetXaxis()->GetXmin(), 1., wd->GetXaxis()->GetXmax(), 1.); one.SetLineStyle(2); one.Draw();
        TLine one(wd->GetXaxis()->GetXmin(), 1., wd->GetXaxis()->GetBinUpEdge(nFine), 1.); one.SetLineStyle(2); one.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.06);
        // tx.DrawLatex(0.17, 0.85, "w normalised within each fit p_{T} bin (80-100, >100 GeV)");   // until 2026-10-28
        // tx.DrawLatex(0.17, 0.85, "w = data/MC, both normalised over 80-400 GeV (jets above 400 GeV dropped)");
        tx.DrawLatex(0.17, 0.85, Form("w = data/MC, both normalised over 80-%.0f GeV (jets above %.0f GeV dropped)", kJetPtMax, kJetPtMax));
        cv.Print(Form("%s/Chi2Scan_JetPtReweighting.png", sDir_canvas.Data()));
    }

    // 2026-10-27 (Zoe): fit m2B on 0-8 GeV, the last bin 7-8 GeV holding everything above 7 GeV
    // (the 10 GeV fold included) -- MergeMassTail(), same as the RooFit fit's fitRange0to8.
    // Zoe, later the same day: "keep the 10 mass bins for now" -> merge disabled, fit on 0-10 GeV.
    // h3D_data = MergeLastMassBinTo8GeV(h3D_data);
    // h3D_bb   = MergeLastMassBinTo8GeV(h3D_bb);
    // h3D_b    = MergeLastMassBinTo8GeV(h3D_b);

    // 2026-10-27, end of day (Zoe): merge the last two m2B bins, 8-9 and 9-10 GeV (the latter already
    // holding everything above 10), into one 8-10 GeV bin -> 9 mass bins. Done after the jet-pT
    // reweighting (rebuild() fills by mass-bin index of the 10-bin _ptfine histograms). The last bin is
    // 2 GeV wide: the fit uses bin contents, so only the plot labels ("per bin") care.
    auto mergeLast2Mass = [](TH3D *h) {
        std::vector<double> ex, ey, ez;
        for (int i = 1; i <= h->GetNbinsX() + 1; ++i) if (i != h->GetNbinsX()) ex.push_back(h->GetXaxis()->GetBinLowEdge(i));
        for (int i = 1; i <= h->GetNbinsY() + 1; ++i) ey.push_back(h->GetYaxis()->GetBinLowEdge(i));
        for (int i = 1; i <= h->GetNbinsZ() + 1; ++i) ez.push_back(h->GetZaxis()->GetBinLowEdge(i));
        TH3D *n = new TH3D(Form("%s_m8to10", h->GetName()), h->GetTitle(), (int) ex.size() - 1, ex.data(),
                           (int) ey.size() - 1, ey.data(), (int) ez.size() - 1, ez.data());
        n->Sumw2(); n->SetDirectory(nullptr);
        const int last = h->GetNbinsX();   // 9-10 GeV goes into the new last bin, with 8-9
        for (int ix = 0; ix <= h->GetNbinsX() + 1; ++ix)
            for (int iy = 0; iy <= h->GetNbinsY() + 1; ++iy)
                for (int iz = 0; iz <= h->GetNbinsZ() + 1; ++iz) {
                    const int jx = (ix >= last) ? ix - 1 : ix;   // last bin and overflow shift down by one
                    const int b = n->GetBin(jx, iy, iz);
                    n->AddBinContent(b, h->GetBinContent(ix, iy, iz));
                    n->SetBinError(b, std::hypot(n->GetBinError(b), h->GetBinError(ix, iy, iz)));
                }
        return n;
    };
    // for (TH3D **h : {&h3D_data, &h3D_bb, &h3D_b, &h3D_0b, &h3D_bb_rw, &h3D_b_rw, &h3D_0b_rw}) *h = mergeLast2Mass(*h);   // 8-10 GeV last bin
    (void) mergeLast2Mass;
    // Then (Zoe, same evening): cut at 8 GeV with 7-8 GeV as the overflow bin -- everything above 7 GeV
    // in it -> 8 mass bins, 0-8 GeV. MergeLastMassBinTo8GeV() = MergeMassTail(h, 8) from Help_Functions.h.
    for (TH3D **h : {&h3D_data, &h3D_bb, &h3D_b, &h3D_0b, &h3D_bb_rw, &h3D_b_rw, &h3D_0b_rw}) {
        TH3D *m = MergeLastMassBinTo8GeV(*h); m->SetDirectory(nullptr); *h = m;
    }

    // const int nPt = jtpt_bins, nX = obs.nbins, nM = mb_bins;   // 0-10 GeV, until the merge above
    const int nPt = jtpt_bins, nX = obs.nbins, nM = h3D_data->GetNbinsX();
    TFile *fout = new TFile(Form("%s/%s", sDirname.Data(), fout_name.Data()), "recreate");
    fout->cd();
    h3D_data->Write(); h3D_bb->Write(); h3D_b->Write();   // the 2B/1B here are the COMBINED templates
    h3D_bb_rw->Write(); h3D_b_rw->Write();                // ... and jet-pT reweighted (the nominal)
    h3D_0b->Write(); h3D_0b_rw->Write();                  // 0B: qcd only
    h_ptspec_data->Write(); h_ptspec_mc->Write(); h_w_jetpt->Write();

    auto book = [&](const char *name) {
        TH2D *h = new TH2D(name, Form(";%s;p_{T}^{jet} [GeV]", obs.axis.Data()), nX, obs.bins, nPt, jtpt_binsVector);
        h->SetDirectory(nullptr); return h;
    };
    TH2D *h_sig_fraction_fit      = book("h_sig_fraction_fit");
    TH2D *h_sig_fraction_fit_up   = book("h_sig_fraction_fit_up");
    TH2D *h_sig_fraction_fit_down = book("h_sig_fraction_fit_down");
    TH2D *h_bkg_fraction_fit      = book("h_bkg_fraction_fit");
    TH2D *h_sig_frac_true_fitbins = book("h_sig_frac_true_fitbins");
    TH2D *h_k_fit                 = book("h_k_fit");
    TH2D *h_chi2ndf_fit           = book("h_chi2ndf_fit");
    TH2D *h_sig_fraction_fit_noRW = book("h_sig_fraction_fit_noRW");   // fit on the templates before reweighting
    TH2D *h_sig_frac_true_rw      = book("h_sig_frac_true_rw");        // 2B fraction of the reweighted MC (jet pT + m2B)
    TH2D *h_sig_fraction_fit_jetptRW = book("h_sig_fraction_fit_jetptRW"); // fit on the jet-pT-reweighted templates only
    std::vector<std::pair<int,int>> failed;

    // grid (see the header for why it differs from the AN's)
    const int nk = 500, nF = 1000;
    // the chi2 maps are kept and drawn after the bin loop, on one common scale (2026-10-28)
    struct MapToDraw { TH2D *h; double k, f; int ix, ip; TString sname; double zmax, chi2min; };
    std::vector<MapToDraw> maps;

    for (int ip = 1; ip <= nPt; ++ip) {
        for (int ix = 1; ix <= nX; ++ix) {
            const TString sname = Form("ptbin_%d_%sbin_%d", ip, (obs.name == "dr") ? "deltaR" : obs.name.Data(), ix);

            TH1D *h_data = h3D_data->ProjectionX(Form("h_data_mb_%d_%d", ix, ip), ix, ix, ip, ip);
            // TH1D *h_2b   = h3D_bb  ->ProjectionX(Form("h_template_2b_%d_%d", ix, ip), ix, ix, ip, ip);   // before the jet-pT reweighting
            // TH1D *h_1b   = h3D_b   ->ProjectionX(Form("h_template_1b_%d_%d", ix, ip), ix, ix, ip, ip);
            // nominal: the jet-pT-reweighted templates; _noRW: the plain ones, for comparison
            TH1D *h_2b      = h3D_bb_rw->ProjectionX(Form("h_template_2b_%d_%d", ix, ip), ix, ix, ip, ip);
            TH1D *h_1b      = h3D_b_rw ->ProjectionX(Form("h_template_1b_%d_%d", ix, ip), ix, ix, ip, ip);
            TH1D *h_2b_noRW = h3D_bb   ->ProjectionX(Form("h_template_2b_noRW_%d_%d", ix, ip), ix, ix, ip, ip);
            TH1D *h_1b_noRW = h3D_b    ->ProjectionX(Form("h_template_1b_noRW_%d_%d", ix, ip), ix, ix, ip, ip);
            for (TH1D *h : {h_data, h_2b, h_1b, h_2b_noRW, h_1b_noRW}) { h->SetDirectory(nullptr); h->GetXaxis()->SetTitle("m_{2B} [GeV]"); }
            // 0B: the background template becomes 1B + 0B, summed BEFORE normalising, i.e. at the MC
            // ratio. h_1only / h_0only keep the separate (jet-pT reweighted) shapes for the plots.
            TH1D *h_0b      = h3D_0b_rw->ProjectionX(Form("h_template_0b_%d_%d", ix, ip), ix, ix, ip, ip);
            TH1D *h_0b_noRW = h3D_0b   ->ProjectionX(Form("h_template_0b_noRW_%d_%d", ix, ip), ix, ix, ip, ip);
            TH1D *h_1only   = (TH1D*) h_1b->Clone(Form("h_template_1bonly_%d_%d", ix, ip));
            TH1D *h_0only   = (TH1D*) h_0b->Clone(Form("h_template_0bonly_%d_%d", ix, ip));
            for (TH1D *h : {h_0b, h_0b_noRW, h_1only, h_0only}) h->SetDirectory(nullptr);
            const double I1only = h_1only->Integral(1, nM), I0only = kWith0B ? h_0only->Integral(1, nM) : 0.;
            if (kWith0B) { h_1b->Add(h_0b); h_1b_noRW->Add(h_0b_noRW); }
            const TString bkgLabel = kWith0B ? "1B+0B" : "1B";
            // c0 and the unit-area 1B / 0B shapes: computed after the m2B reweighting below
            // const double c0 = (I1only + I0only > 0.) ? I0only / (I1only + I0only) : 0.;   // 0B share of the background
            // if (I1only > 0.) h_1only->Scale(1. / I1only);
            // if (I0only > 0.) h_0only->Scale(1. / I0only);
            (void) I1only; (void) I0only;
            styleData(h_data); h_data->SetTitle("Data");
            styleTemplate(h_2b, TFColor::c2B()); h_2b->SetTitle("2B");
            // styleTemplate(h_1b, TFColor::c1B()); h_1b->SetTitle("1B");   // before the 0B
            styleTemplate(h_1b, kWith0B ? TFColor::bkg() : TFColor::c1B()); h_1b->SetTitle(bkgLabel);

            double e2 = 0., e1 = 0.;
            const double I_data = h_data->Integral(1, nM);
            // MC before the reweighting (f_mc, as before) and after it (f_mc_rw)
            const double I2 = h_2b_noRW->IntegralAndError(1, nM, e2), I1 = h_1b_noRW->IntegralAndError(1, nM, e1);
            const double Iw2 = h_2b->Integral(1, nM), Iw1 = h_1b->Integral(1, nM);
            if (!(I_data > 0.) || !(I2 > 0.) || !(I1 > 0.) || !(Iw2 > 0.) || !(Iw1 > 0.)) {
                cout << Form("WARNING: (ptbin %d, bin %d) empty data or template -- no fraction", ip, ix) << endl;
                failed.push_back({ip, ix});
                continue;
            }
            // MC expectation of the 2B fraction (no 0B: 2B / (2B + 1B))
            const double f_mc = I2 / (I1 + I2);   // I1 is 1B+0B when kWith0B: f_MC = 2B / all
            h_sig_frac_true_fitbins->SetBinContent(ix, ip, f_mc);
            h_sig_frac_true_fitbins->SetBinError(ix, ip, std::hypot(I1 * e2, I2 * e1) / std::pow(I1 + I2, 2));

            // unit-area templates (errors scale with them; the integral's own error is ignored, as in the AN)
            // h_2b->Scale(1. / I2); h_1b->Scale(1. / I1);   // m2B-reweighting era: h_2b/h_1b were the plain ones
            h_2b->Scale(1. / Iw2);      h_1b->Scale(1. / Iw1);        // jet-pT reweighted (nominal)
            h_2b_noRW->Scale(1. / I2);  h_1b_noRW->Scale(1. / I1);    // plain
            const double f_mc_rw = f_mc * (Iw2 / I2) / (f_mc * (Iw2 / I2) + (1. - f_mc) * (Iw1 / I1));   // = Iw2/(Iw1+Iw2)

            // the jet-pT-reweighted templates alone, kept for the comparison
            TH1D *h_2b_pt = (TH1D*) h_2b->Clone(Form("h_template_2b_jetptRW_%d_%d", ix, ip));
            TH1D *h_1b_pt = (TH1D*) h_1b->Clone(Form("h_template_1b_jetptRW_%d_%d", ix, ip));
            h_2b_pt->SetDirectory(nullptr); h_1b_pt->SetDirectory(nullptr);

            // ---- m2B reweighting of the SUM, on top of the jet-pT one (AN Sec. 9.6, Figs. 95-97) ----
            // Back as the nominal (Zoe, later 2026-10-27): "you reweight ONLY the shapes you use, and the
            // reweighting is done ONLY for the sum ... you are not fixing it, you are telling the MC that it
            // does not know the fraction, but the data are distributed that way". w(m) = data / MC sum,
            // both unit area, MC sum = f_MC T_2B + (1-f_MC) T_bkg (bkg = 1B+0B, so the sum is 2B+1B+0B like
            // the data), the same w(m) on every template, then renormalised. The fit then gives
            // f = f_MC <w>_2B / (f_MC <w>_2B + (1-f_MC) <w>_bkg), with chi2 = 0 by construction: f follows
            // where in m2B the data/MC excess sits, assuming each mass bin's excess is shared in the MC's
            // local 2B : background proportion.
            // const bool kM2BReweight = true;   // now set at the top, next to kWith0B
            TH1D *h_mc_tot = (TH1D*) h_2b->Clone(Form("h_mc_sum_%d_%d", ix, ip));
            h_mc_tot->Add(h_2b, h_1b, f_mc_rw, 1. - f_mc_rw);
            TH1D *h_rw = (TH1D*) h_data->Clone(Form("h_reweight_%d_%d", ix, ip));   // w(m) = data/MC sum
            h_rw->Scale(1. / I_data);
            for (int i = 1; i <= nM; ++i) {
                const double m = h_mc_tot->GetBinContent(i);
                h_rw->SetBinContent(i, (kM2BReweight && m > 0.) ? h_rw->GetBinContent(i) / m : 1.);
                h_rw->SetBinError  (i, (kM2BReweight && m > 0.) ? h_rw->GetBinError(i)   / m : 0.);
            }
            h_mc_tot->SetDirectory(nullptr); h_rw->SetDirectory(nullptr);
            for (TH1D *h : {h_2b, h_1b, h_1only, h_0only})
                for (int i = 1; i <= nM; ++i) {
                    const double w = h_rw->GetBinContent(i);
                    h->SetBinContent(i, w * h->GetBinContent(i)); h->SetBinError(i, w * h->GetBinError(i));
                }
            const double Im2 = h_2b->Integral(1, nM), Im1 = h_1b->Integral(1, nM);
            if (!(Im2 > 0.) || !(Im1 > 0.)) {
                cout << Form("WARNING: (ptbin %d, bin %d) a reweighted template is empty -- no fraction", ip, ix) << endl;
                failed.push_back({ip, ix});
                continue;
            }
            h_2b->Scale(1. / Im2); h_1b->Scale(1. / Im1);
            const double f_mc_m2b = f_mc_rw * Im2 / (f_mc_rw * Im2 + (1. - f_mc_rw) * Im1);   // 2B fraction of the fully reweighted MC
            const double I1w = h_1only->Integral(1, nM), I0w = kWith0B ? h_0only->Integral(1, nM) : 0.;
            const double c0 = (I1w + I0w > 0.) ? I0w / (I1w + I0w) : 0.;   // 0B share of the background
            if (I1w > 0.) h_1only->Scale(1. / I1w);
            if (I0w > 0.) h_0only->Scale(1. / I0w);

/* ---- disabled (kept for reference): the AN m2B reweighting, replaced by the jet-pT reweighting (Zoe, 2026-10-27) ----
   Circular: w(m2B) = data/MC in the FIT variable forces the fit onto the reweighted MC fraction, chi2 = 0.
            // ---- template reweighting to data (AN Sec. 9.6, Figs. 95-97), 2026-10-27, Zoe ----
            // Per (pT, observable) bin: w(m) = data(m) / MC(m), both unit-area, with MC = the 2B + 1B
            // sum at the MC fractions, f_MC T_2B + (1-f_MC) T_1B. The SAME w(m) multiplies both
            // templates, which are then renormalised. The scan below fits the reweighted ones; the
            // plain-template fit is kept as h_sig_fraction_fit_noRW. Template errors scale with w; the
            // data error on w is not propagated (as in the AN).
            // NOTE: at f = f_MC the reweighted model equals the data in every m2B bin, so the fit lands
            // on f_RW = f_MC * <w>_2B / (f_MC <w>_2B + (1-f_MC) <w>_1B) -- "the fractions of the
            // reweighted simulation", as the AN says -- with chi2 ~ 0. The data/MC shape difference is
            // shared between 2B and 1B in the MC's own local composition, bin by bin in m2B.
            TH1D *h_2b_noRW = (TH1D*) h_2b->Clone(Form("h_template_2b_noRW_%d_%d", ix, ip));
            TH1D *h_1b_noRW = (TH1D*) h_1b->Clone(Form("h_template_1b_noRW_%d_%d", ix, ip));
            TH1D *h_mc_tot  = (TH1D*) h_2b->Clone(Form("h_mc_2b1b_%d_%d", ix, ip));
            h_mc_tot->Add(h_2b, h_1b, f_mc, 1. - f_mc);
            TH1D *h_rw = (TH1D*) h_data->Clone(Form("h_reweight_%d_%d", ix, ip));   // w(m) = data/MC
            h_rw->Scale(1. / I_data);
            for (int i = 1; i <= nM; ++i) {
                const double m = h_mc_tot->GetBinContent(i);
                h_rw->SetBinContent(i, m > 0. ? h_rw->GetBinContent(i) / m : 0.);
                h_rw->SetBinError  (i, m > 0. ? h_rw->GetBinError(i)   / m : 0.);
            }
            for (TH1D *h : {h_2b_noRW, h_1b_noRW, h_mc_tot, h_rw}) h->SetDirectory(nullptr);
            double Iw2 = 0., Iw1 = 0.;
            for (int i = 1; i <= nM; ++i) {
                const double w = h_rw->GetBinContent(i);
                h_2b->SetBinContent(i, w * h_2b->GetBinContent(i)); h_2b->SetBinError(i, w * h_2b->GetBinError(i));
                h_1b->SetBinContent(i, w * h_1b->GetBinContent(i)); h_1b->SetBinError(i, w * h_1b->GetBinError(i));
            }
            Iw2 = h_2b->Integral(1, nM); Iw1 = h_1b->Integral(1, nM);
            if (!(Iw2 > 0.) || !(Iw1 > 0.)) {
                cout << Form("WARNING: (ptbin %d, bin %d) a reweighted template is empty -- no fraction", ip, ix) << endl;
                failed.push_back({ip, ix});
                continue;
            }
            h_2b->Scale(1. / Iw2); h_1b->Scale(1. / Iw1);
            const double f_mc_rw = f_mc * Iw2 / (f_mc * Iw2 + (1. - f_mc) * Iw1);   // 2B fraction of the reweighted MC
*/ // ---- end disabled block ----

            std::vector<double> d(nM), ed(nM);
            for (int i = 0; i < nM; ++i) { d[i] = h_data->GetBinContent(i + 1); ed[i] = h_data->GetBinError(i + 1); }
            int nUsed = 0;
            for (int i = 0; i < nM; ++i) if (ed[i] > 0.) ++nUsed;   // the AN uses only bins with a data error
            const int ndf = nUsed - 2;

            // One (k, f) scan + f profile for a given pair of unit-area templates.
            struct ScanResult { double k, f, f_left, f_right, chi2; TH2D *map; TH1D *prof; };
            auto scan = [&](const TH1D *hs, const TH1D *hb, const TString &tag) -> ScanResult {
                std::vector<double> s(nM), es(nM), b(nM), eb(nM);
                for (int i = 0; i < nM; ++i) {
                    s[i] = hs->GetBinContent(i + 1); es[i] = hs->GetBinError(i + 1);
                    b[i] = hb->GetBinContent(i + 1); eb[i] = hb->GetBinError(i + 1);
                }
                auto chi2At = [&](double k, double f) {
                    double c = 0.;
                    for (int i = 0; i < nM; ++i) {
                        if (!(ed[i] > 0.)) continue;
                        const double m  = k * (f * s[i] + (1. - f) * b[i]);
                        const double vm = k * k * (f * f * es[i] * es[i] + (1. - f) * (1. - f) * eb[i] * eb[i]);
                        c += (d[i] - m) * (d[i] - m) / (ed[i] * ed[i] + vm);
                    }
                    return c;
                };

                // ---- the (k, f) scan ----
                const double kmax = 2. * I_data, dk = kmax / nk, dF = 1. / nF;
                TH2D *hChi2 = new TH2D(Form("hChi2%s_%d_%d", tag.Data(), ix, ip), ";k;2B fraction f;#chi^{2}",
                                       nk + 1, -0.5 * dk, kmax + 0.5 * dk, nF + 1, -0.5 * dF, 1. + 0.5 * dF);
                hChi2->SetDirectory(nullptr);
                double best = 1e30; int iBest = -1;
                for (int i = 0; i <= nk; ++i) {
                    for (int j = 0; j <= nF; ++j) {
                        const double c = chi2At(i * dk, j * dF);
                        hChi2->SetBinContent(i + 1, j + 1, c);
                        if (c < best) { best = c; iBest = i; }
                    }
                }
                const double k_best = iBest * dk;

                // ---- chi2 profile in f at fixed k*, crossings of chi2_min + 1 (clipped to [0, 1]) ----
                // TH1D *h_prof = hChi2->ProjectionY(Form("hChi2Profile_%d_%d", ix, ip), iBest + 1, iBest + 1);   // AN: the map's own column
                // Computed afresh on a 10x finer f grid at k*: the error on f is ~0.004, so the map's
                // 0.001 steps would quantise the chi2_min + 1 crossings. Its minimum refines f*.
                const int nFp = 10000; const double dFp = 1. / nFp;
                TH1D *h_prof = new TH1D(Form("hChi2Profile%s_%d_%d", tag.Data(), ix, ip), ";2B fraction f;#chi^{2}", nFp + 1, -0.5 * dFp, 1. + 0.5 * dFp);
                h_prof->SetDirectory(nullptr);
                int jpBest = 0;
                for (int j = 0; j <= nFp; ++j) {
                    const double c = chi2At(k_best, j * dFp);
                    h_prof->SetBinContent(j + 1, c);
                    if (c < h_prof->GetBinContent(jpBest + 1)) jpBest = j;
                }
                best = h_prof->GetBinContent(jpBest + 1);
                const double c1 = best + 1.;
                double f_left = 0., f_right = 1.;
                for (int j = jpBest; j >= 0; --j)  if (h_prof->GetBinContent(j + 1) >= c1) { f_left  = j * dFp; break; }
                for (int j = jpBest; j <= nFp; ++j) if (h_prof->GetBinContent(j + 1) >= c1) { f_right = j * dFp; break; }
                return {k_best, jpBest * dFp, f_left, f_right, best, hChi2, h_prof};
            };

            const ScanResult r0 = scan(h_2b_noRW, h_1b_noRW, "_noRW");   // plain templates, for comparison
            delete r0.map;
            const ScanResult rpt = scan(h_2b_pt, h_1b_pt, "_jetptRW");   // jet-pT reweighting only, for comparison
            delete rpt.map;
            h_sig_fraction_fit_jetptRW->SetBinContent(ix, ip, rpt.f);
            h_sig_fraction_fit_jetptRW->SetBinError  (ix, ip, std::max(rpt.f - rpt.f_left, rpt.f_right - rpt.f));
            const ScanResult r  = scan(h_2b, h_1b, "");                   // jet pT + m2B reweighted: the nominal
            const double k_best = r.k, f_best = r.f, f_left = r.f_left, f_right = r.f_right, best = r.chi2, c1 = best + 1.;
            TH2D *hChi2 = r.map; TH1D *h_prof = r.prof;
            const double err_lo = f_best - f_left, err_hi = f_right - f_best;
            h_sig_fraction_fit_noRW->SetBinContent(ix, ip, r0.f);
            h_sig_fraction_fit_noRW->SetBinError  (ix, ip, std::max(r0.f - r0.f_left, r0.f_right - r0.f));
            // h_sig_frac_true_rw     ->SetBinContent(ix, ip, f_mc_rw);   // jet-pT reweighting only
            h_sig_frac_true_rw     ->SetBinContent(ix, ip, f_mc_m2b);

            printf("(ptbin %d, %s bin %d): f = %.3f +%.3f -%.3f, chi2/ndf = %.2f/%d (reweighted MC %.3f) | "
                   "jet-pT RW only: f = %.3f, chi2/ndf = %.1f | no RW: f = %.3f, chi2/ndf = %.1f | MC %.3f, 0B share of bkg %.3f%s\n",
                   ip, obs.name.Data(), ix, f_best, err_hi, err_lo, best, ndf, f_mc_m2b, rpt.f, rpt.chi2, r0.f, r0.chi2, f_mc, c0,
                   (f_best >= 1. - 1e-9 || f_best <= 1e-9) ? "  <-- f at the boundary" : "");

            h_sig_fraction_fit     ->SetBinContent(ix, ip, f_best);
            h_sig_fraction_fit     ->SetBinError  (ix, ip, std::max(err_lo, err_hi));
            h_sig_fraction_fit_up  ->SetBinContent(ix, ip, f_right);
            h_sig_fraction_fit_down->SetBinContent(ix, ip, f_left);
            h_bkg_fraction_fit     ->SetBinContent(ix, ip, 1. - f_best);
            h_bkg_fraction_fit     ->SetBinError  (ix, ip, std::max(err_lo, err_hi));
            h_k_fit                ->SetBinContent(ix, ip, k_best);
            h_chi2ndf_fit          ->SetBinContent(ix, ip, ndf > 0 ? best / ndf : 0.);

            // best model and its two components
            TH1D *h_fit_2b = (TH1D*) h_2b->Clone(Form("h_sig_fit_%d_%d", ix, ip)); h_fit_2b->Scale(k_best * f_best);
            // TH1D *h_fit_1b = (TH1D*) h_1b->Clone(Form("h_bkg_fit_1b_%d_%d", ix, ip)); h_fit_1b->Scale(k_best * (1. - f_best));   // before the 0B
            // background split back into 1B and 0B, each with its own shape, at the MC ratio c0
            TH1D *h_fit_1b = (TH1D*) h_1only->Clone(Form("h_bkg_fit_1b_%d_%d", ix, ip));  h_fit_1b->Scale(k_best * (1. - f_best) * (1. - c0));
            TH1D *h_fit_0b = (TH1D*) h_0only->Clone(Form("h_bkg_fit_nob_%d_%d", ix, ip)); h_fit_0b->Scale(k_best * (1. - f_best) * c0);
            styleTemplate(h_fit_1b, TFColor::c1B()); h_fit_1b->SetTitle("1B");
            styleTemplate(h_fit_0b, TFColor::c0B()); h_fit_0b->SetTitle("0B");
            TH1D *h_model  = (TH1D*) h_fit_2b->Clone(Form("h_total_fit_%d_%d", ix, ip)); h_model->Add(h_fit_1b); h_model->Add(h_fit_0b);
            h_model->SetFillStyle(0); h_model->SetLineColor(TFColor::total()); h_model->SetMarkerColor(TFColor::total());
            h_model->SetLineWidth(2); h_model->SetTitle(Form("Model (f = %.0f%%)", 100. * f_best));

            fout->cd();
            // hChi2 itself is not written: 501 x 1001 bins per slice would make the file ~70 MB.
            // h_template_2b/1b are the REWEIGHTED unit-area templates; _noRW the ones before.
            for (TObject *o : {(TObject*) h_data, (TObject*) h_2b, (TObject*) h_1b,
                               (TObject*) h_prof, (TObject*) h_fit_2b, (TObject*) h_fit_1b, (TObject*) h_model,
                               (TObject*) h_2b_noRW, (TObject*) h_1b_noRW, // (TObject*) h_mc_tot, (TObject*) h_rw,   // m2B reweighting, disabled
                               (TObject*) h_fit_0b, (TObject*) h_0only,
                               (TObject*) h_2b_pt, (TObject*) h_1b_pt, (TObject*) h_mc_tot, (TObject*) h_rw, (TObject*) rpt.prof,
                               (TObject*) r0.prof}) o->Write();

            // ---- (0) the jet-pT reweighting seen in m2B: 2B and 1B templates before (dashed) and after
            //          (solid), unit area, ratio after/before ----
            {
                TCanvas cv(Form("c_ptrw_%s", sname.Data()), "", 1100, 1100);
                TPad *p1 = new TPad(Form("pptrw1_%s", sname.Data()), "", 0, 0.2, 1, 1);
                TPad *p2 = new TPad(Form("pptrw2_%s", sname.Data()), "", 0, 0, 1, 0.24);
                p1->SetBottomMargin(0.07); p1->SetLeftMargin(0.18); p1->SetTopMargin(0.14);
                p2->SetTopMargin(0.03); p2->SetBottomMargin(0.40); p2->SetLeftMargin(0.18);
                p1->SetTicks(1, 1); p2->SetTicks(1, 1);
                p1->Draw(); p2->Draw(); p1->cd();
                // TH1D *a2 = (TH1D*) h_2b->Clone(Form("a2_%d_%d", ix, ip)),      *a1 = (TH1D*) h_1b->Clone(Form("a1_%d_%d", ix, ip));   // before the m2B reweighting came back
                TH1D *a2 = (TH1D*) h_2b_pt->Clone(Form("a2_%d_%d", ix, ip)),   *a1 = (TH1D*) h_1b_pt->Clone(Form("a1_%d_%d", ix, ip));
                TH1D *b2 = (TH1D*) h_2b_noRW->Clone(Form("b2_%d_%d", ix, ip)), *b1 = (TH1D*) h_1b_noRW->Clone(Form("b1_%d_%d", ix, ip));
                styleCurve(a2, TFColor::c2B(), true);  styleCurve(b2, TFColor::c2B(), false);
                styleCurve(a1, TFColor::c1B(), true);  styleCurve(b1, TFColor::c1B(), false);
                TH1F *fr = p1->DrawFrame(0., 0., h_data->GetXaxis()->GetXmax(), 1.4 * std::max({a2->GetMaximum(), a1->GetMaximum(), b2->GetMaximum(), b1->GetMaximum()}));
                fr->GetYaxis()->SetTitle("Fraction per bin"); fr->GetXaxis()->SetLabelSize(0);
                for (TH1D *h : {b2, b1, a2, a1}) h->Draw("hist same");
                TLegend *lg = CreateLegend(0.56, 0.44, 0.86, 0.64, {a2, b2, a1, b1}, {"L", "L", "L", "L"},
                                           {"2B, jet-p_{T} reweighted", "2B, not reweighted",
                                            Form("%s, jet-p_{T} reweighted", bkgLabel.Data()), Form("%s, not reweighted", bkgLabel.Data())});
                lg->Draw();
                DrawCommonTextTopRight(p1, ix, ip, obs.bins, nX, false,
                                       Form("f_{MC} = %.1f%% #rightarrow %.1f%%", 100. * f_mc, 100. * f_mc_rw));
                drawCMSHeader(p1);
                p2->cd();
                TH1D *r2 = (TH1D*) a2->Clone(Form("r2_%d_%d", ix, ip)); r2->Divide(b2);
                TH1D *r1 = (TH1D*) a1->Clone(Form("r1_%d_%d", ix, ip)); r1->Divide(b1);
                r2->SetMinimum(0.8); r2->SetMaximum(1.2);
                r2->GetYaxis()->SetTitle("rew./not"); r2->GetYaxis()->CenterTitle(true); r2->GetYaxis()->SetNdivisions(505);
                r2->GetYaxis()->SetTitleSize(0.10); r2->GetYaxis()->SetLabelSize(0.08); r2->GetYaxis()->SetTitleOffset(0.5);
                r2->GetXaxis()->SetTitleSize(0.12); r2->GetXaxis()->SetLabelSize(0.1);
                r2->Draw("hist"); r1->Draw("hist same");
                TLine one(0., 1., h_data->GetXaxis()->GetXmax(), 1.); one.SetLineStyle(2); one.Draw();
                cv.Print(Form("%s/Chi2Scan_TemplatesJetPtRW_%s.png", sDir_canvas.Data(), sname.Data()));
            }
            // ---- (0b) the m2B reweighting: data vs the MC sum (2B + bkg at f_MC, after the jet-pT step),
            //           ratio = w(m2B) (AN Figs. 95-97) ----
            if (kM2BReweight) {
                TCanvas cv(Form("c_rwm_%s", sname.Data()), "", 1100, 1100);
                TPad *p1 = new TPad(Form("prwm1_%s", sname.Data()), "", 0, 0.2, 1, 1);
                TPad *p2 = new TPad(Form("prwm2_%s", sname.Data()), "", 0, 0, 1, 0.24);
                p1->SetBottomMargin(0.07); p1->SetLeftMargin(0.18); p1->SetTopMargin(0.14);
                p2->SetTopMargin(0.03); p2->SetBottomMargin(0.40); p2->SetLeftMargin(0.18);
                p1->SetTicks(1, 1); p2->SetTicks(1, 1);
                p1->Draw(); p2->Draw(); p1->cd();
                TH1D *dn = (TH1D*) h_data->Clone(Form("dnrwm_%d_%d", ix, ip)); dn->Scale(1. / I_data);
                TH1D *mc = (TH1D*) h_mc_tot->Clone(Form("mcrwm_%d_%d", ix, ip));
                mc->SetFillColorAlpha(TFColor::bkg(), 0.5); mc->SetLineColor(TFColor::bkg());
                TH1F *fr = p1->DrawFrame(0., 0., h_data->GetXaxis()->GetXmax(), 1.4 * std::max(dn->GetMaximum(), mc->GetMaximum()));
                fr->GetYaxis()->SetTitle("Fraction per bin"); fr->GetXaxis()->SetLabelSize(0);
                mc->Draw("hist E same"); dn->Draw("PE same");
                TLegend *lg = CreateLegend(0.62, 0.50, 0.86, 0.62, {dn, mc}, {"LPE", "F"},
                                           {"Data", Form("MC (2B+%s), f_{MC} = %.0f%%", bkgLabel.Data(), 100. * f_mc_rw)});
                lg->Draw();
                DrawCommonTextTopRight(p1, ix, ip, obs.bins, nX, false, Form("reweighted MC: f = %.1f%%", 100. * f_mc_m2b));
                drawCMSHeader(p1);
                p2->cd();
                TH1D *wdraw = (TH1D*) h_rw->Clone(Form("wdrawm_%d_%d", ix, ip));
                wdraw->SetMarkerStyle(20); wdraw->SetLineColor(kBlack);
                wdraw->GetYaxis()->SetTitle("Data/MC = w"); wdraw->GetYaxis()->CenterTitle(true);
                wdraw->GetYaxis()->SetNdivisions(505); wdraw->GetYaxis()->SetTitleSize(0.10); wdraw->GetYaxis()->SetLabelSize(0.08);
                wdraw->GetYaxis()->SetTitleOffset(0.5); wdraw->GetXaxis()->SetTitleSize(0.12); wdraw->GetXaxis()->SetLabelSize(0.1);
                wdraw->SetMinimum(0.); wdraw->SetMaximum(std::min(3., 1.2 * wdraw->GetMaximum()));
                wdraw->Draw("EP");
                TLine one(0., 1., h_data->GetXaxis()->GetXmax(), 1.); one.SetLineStyle(2); one.Draw();
                cv.Print(Form("%s/Chi2Scan_Reweighting_%s.png", sDir_canvas.Data(), sname.Data()));
            }
/* ---- disabled (kept for reference): the m2B-reweighting plot (data vs MC, ratio = w(m2B)), replaced above ----
            // ---- (0) the reweighting: data vs MC (2B+1B at the MC fractions), ratio = w(m) (AN Figs. 95-97) ----
            {
                TCanvas cv(Form("c_rw_%s", sname.Data()), "", 1100, 1100);
                TPad *p1 = new TPad(Form("prw1_%s", sname.Data()), "", 0, 0.2, 1, 1);
                TPad *p2 = new TPad(Form("prw2_%s", sname.Data()), "", 0, 0, 1, 0.24);
                p1->SetBottomMargin(0.07); p1->SetLeftMargin(0.18); p1->SetTopMargin(0.14);
                p2->SetTopMargin(0.03); p2->SetBottomMargin(0.40); p2->SetLeftMargin(0.18);
                p1->SetTicks(1, 1); p2->SetTicks(1, 1);
                p1->Draw(); p2->Draw(); p1->cd();
                TH1D *dn = (TH1D*) h_data->Clone(Form("dnrw_%d_%d", ix, ip)); dn->Scale(1. / I_data);
                TH1D *mc = (TH1D*) h_mc_tot->Clone(Form("mcrw_%d_%d", ix, ip));
                mc->SetFillColorAlpha(TFColor::bkg(), 0.5); mc->SetLineColor(TFColor::bkg()); mc->SetTitle("MC (2B + 1B)");
                TH1F *fr = p1->DrawFrame(0., 0., h_data->GetXaxis()->GetXmax(), 1.4 * std::max(dn->GetMaximum(), mc->GetMaximum()));
                fr->GetYaxis()->SetTitle("1/N dN/dm_{2B}"); fr->GetXaxis()->SetLabelSize(0);
                mc->Draw("hist E same"); dn->Draw("PE same");
                TLegend *lg = CreateLegend(0.66, 0.50, 0.86, 0.62, {dn, mc}, {"LPE", "F"},
                                           {"Data", Form("MC, f_{MC} = %.0f%%", 100. * f_mc)});
                lg->Draw();
                DrawCommonTextTopRight(p1, ix, ip, obs.bins, nX, false, Form("reweighted MC: f = %.1f%%", 100. * f_mc_rw));
                drawCMSHeader(p1);
                p2->cd();
                TH1D *wdraw = (TH1D*) h_rw->Clone(Form("wdraw_%d_%d", ix, ip));
                wdraw->SetMarkerStyle(20); wdraw->SetLineColor(kBlack);
                wdraw->GetYaxis()->SetTitle("Data/MC = w"); wdraw->GetYaxis()->CenterTitle(true);
                wdraw->GetYaxis()->SetNdivisions(505); wdraw->GetYaxis()->SetTitleSize(0.10); wdraw->GetYaxis()->SetLabelSize(0.08);
                wdraw->GetYaxis()->SetTitleOffset(0.5); wdraw->GetXaxis()->SetTitleSize(0.12); wdraw->GetXaxis()->SetLabelSize(0.1);
                wdraw->SetMinimum(0.); wdraw->SetMaximum(std::min(3., 1.2 * wdraw->GetMaximum()));
                wdraw->Draw("EP");
                TLine one(0., 1., h_data->GetXaxis()->GetXmax(), 1.); one.SetLineStyle(2); one.Draw();
                cv.Print(Form("%s/Chi2Scan_Reweighting_%s.png", sDir_canvas.Data(), sname.Data()));
            }
*/ // ---- end disabled block ----

            // ---- plots (AN Figs. 40-42, 98-103) ----
            // (1) shapes: unit-area data vs the two templates
            {
                TCanvas cv(Form("c_shapes_%s", sname.Data()), "", 900, 800);
                cv.SetLeftMargin(0.16); cv.SetTopMargin(0.10);
                TH1D *dn = (TH1D*) h_data->Clone(Form("dnorm_%d_%d", ix, ip)); dn->Scale(1. / I_data);
                TH1D *s2 = (TH1D*) h_2b->Clone(Form("s2_%d_%d", ix, ip)); s2->SetFillStyle(0); s2->SetLineColor(TFColor::c2B()); s2->SetLineWidth(2);
                TH1D *s1 = (TH1D*) h_1b->Clone(Form("s1_%d_%d", ix, ip)); s1->SetFillStyle(0); s1->SetLineColor(TFColor::c1B()); s1->SetLineWidth(2);
                const double ymax = 1.5 * std::max({dn->GetMaximum(), s2->GetMaximum(), s1->GetMaximum()});
                TH1F *fr = cv.DrawFrame(0., 0., h_data->GetXaxis()->GetXmax(), ymax);
                fr->GetXaxis()->SetTitle("m_{2B} [GeV]"); fr->GetYaxis()->SetTitle("Fraction per bin");
                s2->Draw("hist E same"); s1->Draw("hist E same"); dn->Draw("PE same");
                TLegend *lg = CreateLegend(0.62, 0.62, 0.86, 0.76, {dn, s2, s1}, {"PE", "L", "L"}, {"Data (2B+1B+0B)", "2B", bkgLabel.Data()});
                lg->Draw();
                DrawCommonTextTopRight((TPad*) &cv, ix, ip, obs.bins, nX, false);
                drawCMSHeader((TPad*) &cv);
                cv.Print(Form("%s/Chi2Scan_Templates_%s.png", sDir_canvas.Data(), sname.Data()));
            }
            // (2) chi2 map in (k, f), best point = red star (AN Figs. 98-100)
            {
                TCanvas cv(Form("c_chi2map_%s", sname.Data()), "", 900, 800);
                cv.SetLeftMargin(0.13); cv.SetRightMargin(0.16); cv.SetTopMargin(0.10);
                // zoom on the valley: +-10% in k, +-10 errors in f, colours up to chi2_min + 25
                const double fw = 10. * std::max(0.005, std::max(err_lo, err_hi));
                // hChi2->SetMaximum(best + 25.); hChi2->SetMinimum(best);   // until 2026-10-28
                hChi2->GetXaxis()->SetRangeUser(0.9 * k_best, 1.1 * k_best);
                // hChi2->GetYaxis()->SetRangeUser(std::max(0., f_best - fw), std::min(1., f_best + fw));   // Zoe: f axis stays 0-1
                (void) fw;
                // 2026-10-28 (Zoe): every map on the scale 1 -> the largest chi2 in the drawn region, log z
                // (that maximum is ~1e6, the minimum ~1-100). A clone with chi2 < 1 raised to 1, so no
                // cell falls below the scale and is left unpainted.
                // TH2D *hMap = (TH2D*) hChi2->Clone(Form("%s_draw", hChi2->GetName())); hMap->SetDirectory(nullptr);
                // Then (Zoe, "with the same limits, to compare bin by bin"): x = k / N_data of the bin, so
                // every map has the same axes, k/N_data 0.9-1.1 and f 0-1 (k ~ N_data for unit-area templates).
                const TAxis *ka = hChi2->GetXaxis();
                TH2D *hMap = new TH2D(Form("%s_draw", hChi2->GetName()), ";k / N_{data};2B fraction f;#chi^{2}",
                                      ka->GetNbins(), ka->GetXmin() / I_data, ka->GetXmax() / I_data,
                                      hChi2->GetNbinsY(), hChi2->GetYaxis()->GetXmin(), hChi2->GetYaxis()->GetXmax());
                hMap->SetDirectory(nullptr);
                for (int i = 1; i <= hMap->GetNbinsX(); ++i)
                    for (int j = 1; j <= hMap->GetNbinsY(); ++j) hMap->SetBinContent(i, j, std::max(1., hChi2->GetBinContent(i, j)));
                hMap->GetXaxis()->SetRangeUser(0.9, 1.1);
                double zmax = 1.;
                for (int i = hMap->GetXaxis()->GetFirst(); i <= hMap->GetXaxis()->GetLast(); ++i)
                    for (int j = 1; j <= hMap->GetNbinsY(); ++j) zmax = std::max(zmax, hChi2->GetBinContent(i, j));
                // hMap->SetMinimum(1.); hMap->SetMaximum(zmax);
                // cv.SetLogz();
                // // hChi2->Draw("COLZ");
                // hMap->Draw("COLZ");
                // TMarker *mk = new TMarker(k_best, f_best, 29); mk->SetMarkerColor(kRed); mk->SetMarkerSize(2.); mk->Draw();
                // DrawCommonTextTopRight((TPad*) &cv, ix, ip, obs.bins, nX, false);
                // drawCMSHeader((TPad*) &cv);
                // cv.Print(Form("%s/Chi2Scan_Map_%s.png", sDir_canvas.Data(), sname.Data()));
                // Then (Zoe): the SAME scale on every map -> drawn after the bin loop, 1 -> the largest
                // chi2 of all maps of this fit (both pT bins, every observable bin).
                // maps.push_back({hMap, k_best, f_best, ix, ip, sname, zmax});
                maps.push_back({hMap, k_best / I_data, f_best, ix, ip, sname, zmax, best});
                hChi2->GetXaxis()->UnZoom(); hChi2->GetYaxis()->UnZoom();
            }
            // (3) chi2 profile in f at k* with the chi2_min + 1 crossings (AN Figs. 101-103)
            {
                TCanvas cv(Form("c_chi2prof_%s", sname.Data()), "", 900, 800);
                cv.SetLeftMargin(0.14); cv.SetTopMargin(0.10);
                h_prof->SetLineColor(kBlack); h_prof->SetLineWidth(2);
                h_prof->GetXaxis()->SetTitle("2B fraction f"); h_prof->GetYaxis()->SetTitle("#chi^{2}");
                h_prof->SetMinimum(best - 0.2); h_prof->SetMaximum(best + 2.);
                // const double fw = 4. * std::max(0.005, std::max(err_lo, err_hi));   // +-4 errors around f*
                // h_prof->GetXaxis()->SetRangeUser(std::max(0., f_best - fw), std::min(1., f_best + fw));   // Zoe: f axis stays 0-1
                h_prof->Draw("hist");
                TLine l0(0., c1, 1., c1); l0.SetLineColor(kRed); l0.SetLineStyle(2); l0.Draw();
                TLine lL(f_left, best - 0.2, f_left, c1 + 1.);  lL.SetLineColor(kRed); lL.SetLineStyle(2); lL.Draw();
                TLine lR(f_right, best - 0.2, f_right, c1 + 1.); lR.SetLineColor(kRed); lR.SetLineStyle(2); lR.Draw();
                TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.03);
                tx.DrawLatex(0.18, 0.80, Form("f = %.1f^{+%.1f}_{-%.1f} %%", 100. * f_best, 100. * err_hi, 100. * err_lo));
                tx.DrawLatex(0.18, 0.75, Form("k = %.4g", k_best));
                tx.DrawLatex(0.18, 0.70, Form("MC: f = %.1f %%", 100. * f_mc));
                DrawCommonTextTopRight((TPad*) &cv, ix, ip, obs.bins, nX, false);
                drawCMSHeader((TPad*) &cv);
                cv.Print(Form("%s/Chi2Scan_Profile_%s.png", sDir_canvas.Data(), sname.Data()));
            }
            // (4) data vs best model, with data/model (AN Figs. 40-42); the two components stacked
            {
                TCanvas cv(Form("c_model_%s", sname.Data()), "", 1100, 1100);
                TPad *p1 = new TPad(Form("p1_%s", sname.Data()), "", 0, 0.2, 1, 1);
                TPad *p2 = new TPad(Form("p2_%s", sname.Data()), "", 0, 0, 1, 0.24);
                p1->SetBottomMargin(0.07); p1->SetLeftMargin(0.18); p1->SetTopMargin(0.14);
                p2->SetTopMargin(0.03); p2->SetBottomMargin(0.40); p2->SetLeftMargin(0.18);
                p1->SetTicks(1, 1); p2->SetTicks(1, 1);
                p1->Draw(); p2->Draw(); p1->cd();
                THStack st(Form("st_%s", sname.Data()), "");
                st.Add(h_fit_2b); st.Add(h_fit_1b);
                if (kWith0B) st.Add(h_fit_0b);
                TH1F *fr = p1->DrawFrame(0., 0., h_data->GetXaxis()->GetXmax(), 1.4 * std::max(h_data->GetMaximum(), h_model->GetMaximum()));
                fr->GetYaxis()->SetTitle("Counts / bin"); fr->GetXaxis()->SetLabelSize(0);
                st.Draw("hist same"); h_model->Draw("hist same"); h_data->Draw("PE same");
                // TLegend *lg = CreateLegend(0.68, 0.42, 0.86, 0.61, {h_data, h_model, h_fit_2b, h_fit_1b},
                //                            {"LPE", "L", "F", "F"}, {"Data", "", "2B", "1B"});   // before the 0B
                TLegend *lg = kWith0B
                    ? CreateLegend(0.68, 0.38, 0.86, 0.61, {h_data, h_model, h_fit_2b, h_fit_1b, h_fit_0b},
                                   {"LPE", "L", "F", "F", "F"}, {"Data", "", "2B", "1B", "0B"})
                    : CreateLegend(0.68, 0.42, 0.86, 0.61, {h_data, h_model, h_fit_2b, h_fit_1b},
                                   {"LPE", "L", "F", "F"}, {"Data", "", "2B", "1B"});
                lg->Draw();
                DrawCommonTextTopRight(p1, ix, ip, obs.bins, nX, false, Form("#chi^{2}/ndf = %.1f/%d", best, ndf));
                drawCMSHeader(p1);
                p2->cd();
                AddRatioPlot(h_data, h_model, "Data/Model");
                cv.Print(Form("%s/Chi2Scan_DataVsModel_%s.png", sDir_canvas.Data(), sname.Data()));
            }
            delete hChi2;
        } // observable bins
    } // pT bins

    // ---- the chi2 maps, all on one colour scale (2026-10-28, Zoe) ----
    {
        // double zall = 1.;
        // for (const MapToDraw &m : maps) zall = std::max(zall, m.zmax);   // log z up to ~1e4: confusing (Zoe)
        // Then (Zoe): LINEAR, 1 -> just above the largest chi2_min of all bins (+25, rounded up to 10,
        // ~160); the rest of each map saturates.
        double cmax = 0.;
        for (const MapToDraw &m : maps) cmax = std::max(cmax, m.chi2min);
        const double zall = 10. * std::ceil((cmax + 25.) / 10.);
        for (MapToDraw &m : maps) {
            TCanvas cv(Form("c_chi2map_%s", m.sname.Data()), "", 900, 800);
            cv.SetLeftMargin(0.13); cv.SetRightMargin(0.16); cv.SetTopMargin(0.10);
            // cv.SetLogz();
            cv.SetRightMargin(0.18); m.h->GetZaxis()->SetTitleOffset(1.3);   // z title off the top label
            m.h->SetMinimum(1.); m.h->SetMaximum(zall);
            m.h->Draw("COLZ");
            TMarker *mk = new TMarker(m.k, m.f, 29); mk->SetMarkerColor(kRed); mk->SetMarkerSize(2.); mk->Draw();
            DrawCommonTextTopRight((TPad*) &cv, m.ix, m.ip, obs.bins, nX, false);
            drawCMSHeader((TPad*) &cv);
            cv.Print(Form("%s/Chi2Scan_Map_%s.png", sDir_canvas.Data(), m.sname.Data()));
            delete m.h;
        }
        printf("chi2 maps: common colour scale 1 - %.3g\n", zall);
    }

    // ---- 2B fraction maps, from data and from MC, in % (AN Fig. 43) ----
    for (int w = 0; w < 2; ++w) {
        TH2D *src = w ? h_sig_frac_true_fitbins : h_sig_fraction_fit;
        TH2D *pct = (TH2D*) src->Clone(Form("%s_percent", src->GetName()));
        pct->Scale(100.);
        TCanvas cv(Form("c_frac_%d", w), "", 1000, 800);
        cv.SetLeftMargin(0.12); cv.SetRightMargin(0.15); cv.SetTopMargin(0.10);
        gStyle->SetPaintTextFormat(".1f");
        pct->SetMinimum(0.); pct->SetMaximum(100.);
        pct->GetZaxis()->SetTitle(w ? "2B fraction, MC [%]" : "2B fraction, data [%]");
        pct->Draw("COLZ TEXT");
        drawCMSHeader((TPad*) &cv);
        cv.Print(Form("%s/Chi2Scan_SigFraction_%s.png", sDir_canvas.Data(), w ? "MC" : "Data"));
        delete pct;
    }

    // ---- the 2B fraction vs the observable from each method (AN Fig. 104) ----
    for (int ip = 1; ip <= nPt; ++ip) {
        TH1D *hm  = h_sig_frac_true_fitbins->ProjectionX(Form("fcmp_mc_%d", ip),   ip, ip);
        TH1D *hn  = h_sig_fraction_fit_noRW->ProjectionX(Form("fcmp_norw_%d", ip), ip, ip);
        TH1D *hr  = h_sig_fraction_fit     ->ProjectionX(Form("fcmp_rw_%d", ip),   ip, ip);
        TH1D *hp  = h_sig_fraction_fit_jetptRW->ProjectionX(Form("fcmp_ptrw_%d", ip), ip, ip);
        for (TH1D *h : {hm, hn, hr, hp}) { h->SetDirectory(nullptr); h->GetXaxis()->SetTitle(obs.axis); }
        styleCurve(hm, kBlue + 1, false);           hm->SetTitle("MC (before reweighting)");
        styleCurve(hn, TFColor::cmpB(), true);      hn->SetTitle("Fit, templates not reweighted");
        styleCurve(hp, TFColor::cmpA(), true);      hp->SetTitle("Fit, jet-p_{T} reweighting only");
        // styleCurve(hr, kRed + 1, true);             hr->SetTitle("Fit, reweighted templates (nominal)");   // m2B reweighting
        // styleCurve(hr, kRed + 1, true);             hr->SetTitle("Fit, jet-p_{T}-reweighted templates (nominal)");
        styleCurve(hr, kRed + 1, true);
        hr->SetTitle(kM2BReweight ? "Fit, jet-p_{T} + m_{2B} reweighting (nominal)" : "Fit, jet-p_{T} reweighting (nominal)");
        TCanvas cv(Form("c_fcmp_%d", ip), "", 1100, 1100);
        TPad *p1 = new TPad(Form("pfc1_%d", ip), "", 0, 0.2, 1, 1);
        TPad *p2 = new TPad(Form("pfc2_%d", ip), "", 0, 0, 1, 0.24);
        p1->SetBottomMargin(0.07); p1->SetLeftMargin(0.15); p1->SetTopMargin(0.12);
        p2->SetTopMargin(0.03); p2->SetBottomMargin(0.40); p2->SetLeftMargin(0.15);
        p1->SetTicks(1, 1); p2->SetTicks(1, 1);
        p1->Draw(); p2->Draw(); p1->cd();
        hm->SetMinimum(0.); hm->SetMaximum(1.35); hm->GetYaxis()->SetTitle("2B fraction f");
        hm->GetXaxis()->SetLabelSize(0);
        // without the m2B reweighting the jet-pT-only fit IS the nominal: drawn once
        hm->Draw("E1"); hn->Draw("E1 same"); if (kM2BReweight) hp->Draw("E1 same"); hr->Draw("E1 same");
        TLegend *lg = kM2BReweight ? CreateLegend(0.18, 0.66, 0.60, 0.86, {hm, hn, hp, hr}, {"LP", "LP", "LP", "LP"})
                                   : CreateLegend(0.18, 0.70, 0.60, 0.86, {hm, hn, hr}, {"LP", "LP", "LP"});
        lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.03);
        tx.DrawLatex(0.62, 0.82, pt_label(jtpt_binsVector[ip - 1], jtpt_binsVector[ip]));
        drawCMSHeader(p1);
        p2->cd();
        TH1D *rn = (TH1D*) hn->Clone(Form("fcmp_rn_%d", ip)); rn->Divide(hm);
        TH1D *rr = (TH1D*) hr->Clone(Form("fcmp_rr_%d", ip)); rr->Divide(hm);
        TH1D *rp = (TH1D*) hp->Clone(Form("fcmp_rp_%d", ip)); rp->Divide(hm);
        // rr->SetMinimum(0.6); rr->SetMaximum(1.4 > rr->GetMaximum() ? 1.4 : 2.5);   // cut off the noRW points
        rr->SetMinimum(0.6); rr->SetMaximum(std::min(3., std::max(1.4, 1.1 * std::max(rr->GetMaximum(), rn->GetMaximum()))));
        rr->GetYaxis()->SetTitle("Fit / MC"); rr->GetYaxis()->CenterTitle(true); rr->GetYaxis()->SetNdivisions(505);
        rr->GetYaxis()->SetTitleSize(0.10); rr->GetYaxis()->SetLabelSize(0.08); rr->GetYaxis()->SetTitleOffset(0.5);
        rr->GetXaxis()->SetTitleSize(0.12); rr->GetXaxis()->SetLabelSize(0.1);
        rr->Draw("E1"); rn->Draw("E1 same"); if (kM2BReweight) rp->Draw("E1 same");
        TLine one(obs.bins[0], 1., obs.bins[nX], 1.); one.SetLineStyle(2); one.Draw();
        cv.Print(Form("%s/Chi2Scan_FractionMethods_ptbin_%d.png", sDir_canvas.Data(), ip));
    }

    // ---- data yield vs the observable, split into 2B and 1B with the fitted f (AN Fig. 44) ----
    for (int ip = 1; ip <= nPt; ++ip) {
        TH1D *h_tot = h3D_data->ProjectionY(Form("h_yield_data_pt%d", ip), 1, nM, ip, ip);
        h_tot->SetDirectory(nullptr);
        TH1D *h_sig = (TH1D*) h_tot->Clone(Form("h_yield_2b_pt%d", ip));
        TH1D *h_bkg = (TH1D*) h_tot->Clone(Form("h_yield_1b_pt%d", ip));
        for (int ix = 1; ix <= nX; ++ix) {
            const double f = h_sig_fraction_fit->GetBinContent(ix, ip);
            h_sig->SetBinContent(ix, f * h_tot->GetBinContent(ix));        h_sig->SetBinError(ix, f * h_tot->GetBinError(ix));
            h_bkg->SetBinContent(ix, (1. - f) * h_tot->GetBinContent(ix)); h_bkg->SetBinError(ix, (1. - f) * h_tot->GetBinError(ix));
        }
        for (TH1D *h : {h_tot, h_sig, h_bkg}) { h->Scale(1., "width"); h->GetXaxis()->SetTitle(obs.axis); }
        styleData(h_tot);   h_tot->SetMarkerStyle(kOpenCircle);
        // styleCurve(h_sig, TFColor::c2B(), true);   // until 2026-10-28: red points
        // 2026-10-28 (Zoe): the 2B filled like the 1B+0B, same transparency, both drawn from 0 (not stacked)
        styleTemplate(h_sig, TFColor::c2B()); h_sig->SetFillColorAlpha(TFColor::c2B(), 0.35);
        styleTemplate(h_bkg, TFColor::c1B()); h_bkg->SetFillColorAlpha(TFColor::c1B(), 0.35);
        TCanvas cv(Form("c_sep_%d", ip), "", 1000, 1000);
        cv.SetLeftMargin(0.15); cv.SetTopMargin(0.10);
        // if (obs.name == "dr") cv.SetLogx();   // disabled: breaks if the first dR edge is 0
        const double ymax = 1.5 * h_tot->GetMaximum();
        h_tot->SetMaximum(ymax); h_tot->SetMinimum(0.);
        h_tot->GetYaxis()->SetTitle(Form("dN/d%s", obs.axis.Data()));
        // h_tot->Draw("PE1"); h_bkg->Draw("hist same"); h_sig->Draw("PE1 same"); h_tot->Draw("PE1 same");
        h_tot->Draw("PE1"); h_bkg->Draw("hist same"); h_sig->Draw("hist same"); h_tot->Draw("PE1 same");
        // TLegend *lg = CreateLegend(0.55, 0.68, 0.86, 0.84, {h_tot, h_sig, h_bkg}, {"LPE", "LPE", "F"},
        TLegend *lg = CreateLegend(0.55, 0.68, 0.86, 0.84, {h_tot, h_sig, h_bkg}, {"LPE", "F", "F"},
                                   {"Data, all b-tagged", "2B (f #times data)",
                                    kWith0B ? "1B+0B ((1-f) #times data)" : "1B ((1-f) #times data)"});
        lg->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.03);
        tx.DrawLatex(0.55, 0.64, pt_label(jtpt_binsVector[ip - 1], jtpt_binsVector[ip]));
        drawCMSHeader((TPad*) &cv);
        cv.Print(Form("%s/Chi2Scan_Separation_ptbin_%d.png", sDir_canvas.Data(), ip));
        fout->cd(); h_tot->Write(); h_sig->Write(); h_bkg->Write();
    }

    // Like do_template_fit_combined(): a bin without a fraction must stop the unfolding, not
    // empty its data bin, so h_sig_fraction_fit is then not written.
    fout->cd();
    for (TH2D *h : {h_sig_fraction_fit, h_sig_fraction_fit_up, h_sig_fraction_fit_down, h_bkg_fraction_fit,
                    h_sig_frac_true_fitbins, h_k_fit, h_chi2ndf_fit, h_sig_fraction_fit_noRW, h_sig_frac_true_rw,
                    h_sig_fraction_fit_jetptRW}) {
        if (!failed.empty() && (h == h_sig_fraction_fit || h == h_bkg_fraction_fit)) continue;
        h->Write();
    }
    for (auto p : failed) cout << "No fraction for (pT " << p.first << ", " << obs.name << " " << p.second << ")" << endl;
    if (!failed.empty())
        std::cerr << "ERROR: " << failed.size() << " bin(s) without a 2B fraction -- h_sig_fraction_fit NOT written" << std::endl;
    fout->Close();
}

// SAMPLE:    qcd | both
//            qcd  = fit with the dijet templates only (h3D_b, h3D_bb, h3D_0b)
//            both = also take h3D_b and h3D_bb from the bjet sample (also_bjet)
//            There is deliberately no "bjet": h3D_0b exists only in the qcd sample, the
//            bjet sample being filtered to b jets, so a bjet-only fit has no 0B template.
// GENERATOR: pythia | herwig -- which MC the templates come from. The data file is the
//            same either way. Flags apply to RunN 3; the Run2 branch is untouched.
//
// TRACK_EFF_UNC: false (nominal) | true -- refit with the templates from the 3%
//            track-drop production, so the signal fraction moves with the tracking
//            efficiency too. The DATA being fitted is the same file either way; only the
//            MC templates change. Writes to its own _trkdrop030 output directory, so the
//            nominal fit is never overwritten. Pythia only.
//
// e.g.  root -l -b -q 'template_fit.cpp("both","pythia")'
//       root -l -b -q 'template_fit.cpp("both","herwig")'
//       root -l -b -q 'template_fit.cpp("both","pythia",true)'   // tracking variation
void template_fit(TString SAMPLE = "both", TString GENERATOR = "pythia",
                  bool TRACK_EFF_UNC = false, bool EEC_WEIGHT_OFF = false,
                  TString OBSERVABLE = "dr"){

    // Which axis the fit slices in. "dr" is the original behaviour in every respect --
    // same input files, same histogram names, same output directory -- so an existing call
    // that does not pass this argument is unaffected.
    const ObsDef obs = obsByName(OBSERVABLE);
    if (obs.nbins == 0) {
        std::cerr << "ERROR: unknown OBSERVABLE '" << OBSERVABLE << "' (use dr | B)" << std::endl;
        return;
    }
    // Set BEFORE anything is read or booked: CheckInputBinning(), the histogram names and
    // the plot labels all read gFitObs().
    setFitObservable(obs);

    /* ---- disabled (kept for reference): B was Pythia-only, no track drop ----
    // True until 2026-09-24, when the B band grew to the full dr set. mcgenTemplates() now
    // builds the Herwig and track-drop B paths, and a production that was never made fails
    // on the missing-file check below, naming the file.
    if (OBSERVABLE != "dr" && (GENERATOR != "pythia" || TRACK_EFF_UNC)) {
        std::cerr << "ERROR: the " << OBSERVABLE << " templates exist only as Pythia "
                  << "productions with no "
                  << "track-drop variation. Asked for generator '" << GENERATOR
                  << "', track_eff_unc=" << TRACK_EFF_UNC << "." << std::endl;
        return;
    }
    ---- end disabled ---- */
    // For B, EEC_WEIGHT_OFF is not a variation -- it IS the measurement. The EEC weight
    // pT_b1*pT_b2 equals B(1-B)S², an analytic function of B itself, so an EEC-weighted B
    // distribution is the yield times a known kinematic factor. See README.md (repo root).
    if (nominalEecWeightOff(OBSERVABLE) && !EEC_WEIGHT_OFF)
        std::cout << "NOTE: fitting the EEC-WEIGHTED " << OBSERVABLE << " templates. The "
                  << "measurement is the YIELD run (EEC_WEIGHT_OFF=true); this one is the "
                  << "comparison." << std::endl;

    if (SAMPLE == "bjet") {
        std::cerr << "ERROR: SAMPLE 'bjet' is not a valid template fit: the 0B template "
                  << "(h3D_0b) exists only in the qcd sample, the bjet sample being "
                  << "filtered to b jets. Use qcd or both." << std::endl;
        return;
    }
    if (SAMPLE != "qcd" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | both)" << std::endl;
        return;
    }
    if (GENERATOR != "pythia" && GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown GENERATOR '" << GENERATOR << "' (use pythia | herwig)" << std::endl;
        return;
    }

    // One output directory per flag combination. sDirname / sDirname_www are the globals
    // from Help_Functions.h that every drawing function writes through, so set them here
    // before anything is created.
    // Same tag create_files_for_template_fit.cpp puts on the varied files, so the fit
    // that used them is identifiable from its directory name alone. "" when nominal, so
    // nominal output paths are unchanged.
    const TString trk_tag  = TRACK_EFF_UNC  ? "_trkdrop030" : "";
    // The yield fit gets its own directory for the same reason: it is a different
    // observable, not a variation, and must never overwrite the EEC fit.
    const TString eecw_tag = EEC_WEIGHT_OFF ? "_noeecw" : "";
    // Observable tag: "" for dR so the dR fit keeps its existing directory, "_B" otherwise.
    // A B fit must never land in the dR fit's folder -- apply_unfolding_2d.C reads
    // h_sig_fraction_fit by name, and a B-binned one there would be silently wrong.
    const TString obs_tag  = (OBSERVABLE == "dr") ? "" : ("_" + OBSERVABLE);
    // sDirname     = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/"
    //                "TemplateFits_" + SAMPLE + "_" + GENERATOR + trk_tag + svsfTag()
    //                + trackBdtSfTag() + eecw_tag + obs_tag + eraTag();   // the RooFit fit, until 2026-10-27
    // 2026-10-27: the AN-style chi2-scan fit (do_template_fit_chi2scan) writes to its own
    // "_chi2scan" directory, so the RooFit fits that apply_unfolding_2d.C reads stay untouched.
    sDirname     = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/"
                   "TemplateFits_" + SAMPLE + "_" + GENERATOR + trk_tag + svsfTag()
                   + trackBdtSfTag() + eecw_tag + obs_tag + eraTag() + "_chi2scan";
    sDirname_www = sDirname;

    // -- Output folder to save the result of the tests
    gSystem->mkdir(sDirname, kTRUE);// Predefined in Help.h -- holds the root files
    gSystem->mkdir(sDirname_www, kTRUE);// single flat folder holding every png
    // TString folder = Form("/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/TemplateFit_Run3/%s/", sDirname.Data()); // this is sDirname
    TString folder = Form("%s", sDirname.Data()); // this is sDirname
        cout << "Sample: " << SAMPLE << ", generator: " << GENERATOR
             << ", tracking-eff variation: " << (TRACK_EFF_UNC ? "ON (3% tracks dropped)" : "off")
             << ", SV b SF: " << (svsfTag().Length() ? svsfTag() : TString("off"))
             << ", track-BDT SF: " << (trackBdtSfTag().Length() ? trackBdtSfTag() : TString("off")) << endl;
        cout << "Output folder path: "<< folder << endl;

    Int_t RunN = 3; // 3;


    //Get data and mc labels
    TString pT_selection = "80_inf";
  
    // Add LowEG data 
    bool alsoLowEG = false; 
    bool also_bjet = false;

    // -- Input data/MC files names
    TString dataset_HG = ""; 
    TString dataset_LG = ""; 
    TString dataname = "All";
    TString templates_dijet = ""; 
    TString templates_bjet = "";
    TString fout_name = "TemplateFits_histos_3d_" + pT_selection +  ".root";

    // -- Set data/MC samples to use
    if (RunN == 3){
        alsoLowEG = false;
        // "both" adds the bjet sample's b/bb templates on top of the dijet ones.
        also_bjet = (SAMPLE == "both");
        // btagWP<NNN> follows BTAG_WP in the run scripts.
        // dataset_HG = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_upartv2.root"; // Does not exist!

        // Data is the same whichever generator the templates come from -- but NOT the same
        // whichever observable is being measured. h3D_data is filled with eec * weight_tree
        // like everything else, so the EEC-weighted data file would be fitted against
        // unweighted templates, which mixes two observables. The yield run therefore takes
        // its own data production (make_hardprobes_condor_scripts.sh with EEC_WEIGHT_OFF=true,
        // hadd'ed to the agg_template_chunks top level).
        // The B data file is the "_upartv2_B" production, for the same reason as the
        // templates: Afnan's data file has no h3D_data_B. Data and templates must come from
        // productions that agree on the observable, or the fit is meaningless.
        // ⚠️ h3D_data_<obs>, not h_count_data_<obs>: with the EEC weight off the fill
        // weight is 1 * prescale, so h3D_data carries the TRIGGER PRESCALE and h_count does
        // not (the prescale multiplies eec, never the count). The yield measurement wants
        // the prescale-corrected one -- 1.949M against 1.773M raw jets.
        dataset_HG = dataTemplateFile(EEC_WEIGHT_OFF, OBSERVABLE);
        templates_dijet = mcgenTemplates("qcd",  GENERATOR, TRACK_EFF_UNC, EEC_WEIGHT_OFF, OBSERVABLE);
        templates_bjet  = mcgenTemplates("bjet", GENERATOR, TRACK_EFF_UNC, EEC_WEIGHT_OFF, OBSERVABLE);
        fout_name = Form("Run%d_TemplateFits_histos_3d_%s.root", RunN, pT_selection.Data());

        if (templates_dijet.Length() == 0 ||
            (also_bjet && templates_bjet.Length() == 0)) {
            std::cerr << "ERROR: no MC templates for generator '" << GENERATOR
                      << "' -- add the paths to mcgenTemplates()" << std::endl;
            return;
        }
        for (const TString &f : { templates_dijet, also_bjet ? templates_bjet : templates_dijet }) {
            if (gSystem->AccessPathName(f)) {
                std::cerr << "ERROR: missing template file " << f << std::endl;
                return;
            }
        }
        cout << "Templates (dijet): " << templates_dijet << endl;
        if (also_bjet) cout << "Templates (bjet) : " << templates_bjet << endl;
    }
    else if (RunN == 2){
        alsoLowEG = true;
        also_bjet = true;
        dataset_HG = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_HighEG_btag_0990_small_bins.root"; 
        dataset_LG = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_LowEG_btag_0990_small_bins.root"; 
        templates_dijet = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_qcd_btag_0990_small_bins.root";
        templates_bjet = " /data_CMS/cms/zaidan/analysis_lise/pulido/small_bins/template_for_fit_histos_3D_bjet_btag_0990_small_bins.root";
        fout_name = Form("Run%d_TemplateFits_histos_3d_%s.root", RunN, pT_selection.Data());
    }


    // --- Start work from here -----
    // 2026-10-27: the AN-style chi2 scan, 2B vs 1B, no 0B, event-weighted templates. Everything
    // below it (the RooFit fit per variation, its summary plots, the 0B systematic) is the old
    // path, commented out: not used by this fit.
    if (RunN == 3) {
        // The jet-pT reweighting needs the _ptfine histograms, which only the "<era>ptw" step-1
        // production has (OUT_TAG=upartv2rhojecvetobenrptw[_B], 2026-10-27): same content otherwise.
        // for (TString *p : {&dataset_HG, &templates_dijet, &templates_bjet}) p->ReplaceAll(eraTag(), eraTag() + "ptw");
        // JET_PT_MAX chains (2026-10-29) use ONE production per variation: the _ptfine histograms are in the
        // plain productions since step 1 fills them always, so no separate "ptw" set.
        // Same for the SV_SEED chain (2026-11-05).
        // if (ptMaxTag().Length() == 0)   // until 2026-11-05
        if (ptMaxTag().Length() == 0 && svSeedTag().Length() == 0)
            for (TString *p : {&dataset_HG, &templates_dijet, &templates_bjet}) p->ReplaceAll(eraTag(), eraTag() + "ptw");
        cout << "chi2 scan reads the ptw production:\n  " << dataset_HG << "\n  " << templates_dijet
             << "\n  " << templates_bjet << endl;
        for (const TString &f : {dataset_HG, templates_dijet, templates_bjet})
            if (gSystem->AccessPathName(f)) { std::cerr << "ERROR: missing " << f << std::endl; return; }
        // 2026-10-28: the light-jet mistag variations, the same files the RooFit fit wrote
        // (var0B_2_ / var0B_0_, read by apply_unfolding_2d.C with TF_VARIATION), their plots in
        // subfolders. Run before the nominal so nothing of theirs can land on the nominal's plots.
        const TString www_nom = sDirname_www;
        for (const auto &v : std::vector<std::pair<TString, double>>{ {"var0B_2", 2.}, {"var0B_0", 0.} }) {
            sDirname_www = www_nom + "/" + v.first;
            gSystem->mkdir(sDirname_www, kTRUE);
            do_template_fit_chi2scan(dataset_HG, templates_dijet, templates_bjet, also_bjet, v.first + "_" + fout_name, v.second);
        }
        sDirname_www = www_nom;
        TString newfout_name = "nominal_" + fout_name;
        do_template_fit_chi2scan(dataset_HG, templates_dijet, templates_bjet, also_bjet, newfout_name);
        return;
    }

/* ---- disabled (kept for reference): RooFit fit with 0B, its variations, summary plots and 0B systematic (not used by the chi2 scan) ----
    TString sfoutputPlots_dijet = Form("Run%d_Summary_histo_templatefits.root", RunN);
    TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "RECREATE");
        if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {
            std::cout << "Error opening file!" << std::endl;
            return;
        }
    
       
        // -- Draw Prefit templates for Run3 MC (qcd sample)
            // Draw_template_Run3(templates_dijet, pT_selection, folder, fout_name);

        // --  Loop over variations on templates: one root file per variation 
        // --  other png drawings are on seperate directories, for simplisity.
        for (int ivar = 0; ivar < 4; ivar++) //
        // for (int ivar = 0; ivar < 1; ivar++) //
        {
            TString newfout_name = varNames[ivar]+ "_" + fout_name;
            do_template_fit_combined(dataset_HG,dataset_LG,templates_dijet, templates_bjet,  pT_selection, folder, newfout_name, alsoLowEG, also_bjet, (Variation) ivar); // default: NOMINAL variation 

            // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++) {  // test 
            for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
                /// Draw S/B fractions 
                // -- to test Draw fraction only: READ foutputPlots_dijet instead of RECREATE
                    // TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Update"); 
                    // if (! foutputPlots_dijet->IsOpen()){ foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Read"); }
                    // if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {std::cout << "Error opening file!" << std::endl; return;}
                
                draw_template_fit_result(newfout_name, foutputPlots_dijet, dataname, folder, pT_selection, ibin_pt, (Variation) ivar); 
                
                /// Draw EEC 
                draw_eec_simple(newfout_name, foutputPlots_dijet ,folder, also_bjet, ibin_pt, (Variation) ivar);

            } // end loop over ibin_pt

        } // end loop over ivar




        //-- Get 0B systematics 
        cout << "Calculate systematic uncertaintiy " << endl;
        TFile *fsys = new TFile( Form("%s/Result_syst_uncert_templatefit.root", sDirname.Data()),"recreate");
                    /// to test seperatly
                    // TFile *foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Update"); 
                    // if (! foutputPlots_dijet->IsOpen()){ foutputPlots_dijet = new TFile(Form("%s/%s", sDirname.Data(), sfoutputPlots_dijet.Data()), "Read"); }
                    // if (!foutputPlots_dijet || foutputPlots_dijet->IsZombie()) {std::cout << "Error opening file!" << std::endl; return;}
            for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++) {
            // for(Int_t ibin_pt = 1; ibin_pt <= 1; ibin_pt++) {  // test 

                    draw_variation_uncertainity(foutputPlots_dijet, fsys, ibin_pt);
            }

            fsys->Print();
            fsys->Close();
            delete fsys;
    
       


    foutputPlots_dijet->Print();
    foutputPlots_dijet->Close();
    delete foutputPlots_dijet;
*/ // ---- end disabled block ----

}


// ============================================================================
// plot_z_first_look() -- the pre-fit first look at an observable. Moved in from
// plot_z_first_look.C on 2026-09-25 to keep the number of files down; unchanged. It reads the
// same templates this file fits, and uses the same palette (TFColor, styleTemplate).
//
// First look at a measured observable, straight out of step 1 of the chain -- BEFORE any
// template fit, signal fraction, unfolding or correction. This is the pre-fit step of the
// analysis, the one that answers "is this observable measurable at all".
//
// ⚠️ THE NAME IS HISTORICAL. The momentum balance it was written for is called B, not z,
// since 2026-09-22; the function keeps its name so the command lines and the plot folders
// stay recognisable. Read "z" in plot_z_first_look as "the balance".
//
// It takes the observable as its LAST argument:
//
//   root -l -b -q -e '.L template_fit.cpp+' \
//                 -e 'plot_z_first_look("both","pythia","noeecw_upartv2_B","B")'
//
// The tag must be a production that actually filled that observable's axis -- checked
// against obsProdTag(), not assumed. ⚠️ Pre-rename productions ("_zfirst", "_3obs",
// "_lnfblin", "_fb") name the balance axis "_z" and will not be found under "B".
//
// WHAT THIS IS, AND WHAT IT IS NOT
// --------------------------------
// Nothing here is a measurement. The data curve is the raw EEC-weighted 2b-tagged yield per
// B bin: it still contains the mistagged and single-b background that the template fit
// exists to subtract, it has not been unfolded, and no efficiency, purity or EEC-weight
// correction has been applied. The MC curves are reco-level too, so data and MC are
// compared like for like -- but neither is a particle-level number.
//
// ⚠️ And unlike the dR chain, NO UParT scale factor is applied anywhere here: that
// calibration is binned in reco dR and has no B equivalent (2026-09-15 decision). The dR
// result has a correction that this does not.
//
// The point is to answer, before the template fit is written:
//
//   1. do the flavour templates SEPARATE in B? (plot 1)
//      If the 2b, 1b and 0b B shapes were identical, a fit per B bin would be pointless --
//      the B dependence of the signal fraction would be flat and could be taken from dR.
//   2. does m_2B still separate signal from background INSIDE each B bin? (plot 2)
//      This is the one that decides whether the fit will work: the fit discriminates on
//      m_2B, sliced in (B, pT), so each slice needs both shape separation and statistics.
//      This is the most important plot here.
//   3. is the response well behaved in B? (plot 3)
//      Purity, efficiency and the migration matrix on the chosen 4-bin B binning.
//   4. do data and MC look like each other in B at all? (plot 4)
//
// RUN (on LLR, after source setup_roounfold_env.sh):
//
//   root -l -b -q -e '.L template_fit.cpp+' -e 'plot_z_first_look("qcd","pythia","upartv2_B")'
//
// The third argument is the production tag (OUT_TAG in the run scripts). It defaults to the
// "_upartv2_B" production rather than the nominal, so this macro cannot silently read a
// nominal file that has no balance histograms and report empty plots.
// ============================================================================

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iomanip>
#include <iostream>
#include <vector>
#include "THStack.h"

namespace {

const char *kBase    = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
const char *kBtagTag = "btagWP0712";

// Analysis palette: red = measurement, blue = particle-level MC, green/purple for the rest.
// Data here is not a measurement yet (no fit, no unfolding), but it IS the data, so it keeps
// red; the MC flavour templates take the "rest" colours.
const int kColData = kRed + 1;
const int kCol2b   = kBlue + 1;    // true 2b -- the signal
const int kCol1b   = kGreen + 2;   // single b
const int kCol0b   = kMagenta + 2; // mistagged light/charm

TString sampleSubdir(const TString &sample, const TString &generator)
{
    const bool herwig = (generator == "herwig");
    if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
    if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
    return "";
}

// Sum one named histogram over every per-block file that exists. Same "take what is on
// disk" rule as aggChunkFiles() in apply_unfolding_2d.C -- block counts differ per sample.
template <class TH>
TH *sumOverBlocks(const std::vector<TString> &files, const char *hname, const char *newname)
{
    TH *out = nullptr;
    int used = 0, missing = 0;
    for (const TString &f : files) {
        TFile *fin = TFile::Open(f);
        if (!fin || fin->IsZombie()) { delete fin; continue; }
        TH *h = dynamic_cast<TH *>(fin->Get(hname));
        if (!h) { missing++; fin->Close(); delete fin; continue; }
        if (!out) { out = (TH *) h->Clone(newname); out->SetDirectory(nullptr); }
        else        out->Add(h);
        used++;
        fin->Close();
        delete fin;
    }
    if (!out)
        std::cerr << "ERROR: '" << hname << "' not found in any of the " << files.size()
                  << " input files. If this is a B histogram, the production predates the "
                  << "B observable (or predates the z -> B rename) -- re-run "
                  << "run_agg_ntuple_chunks.sh." << std::endl;
    else if (missing)
        std::cerr << "WARNING: '" << hname << "' missing from " << missing << " of "
                  << (used + missing) << " files" << std::endl;
    return out;
}

// "both" = qcd + bjet, the nominal configuration. The histogram NAMES are identical in the
// two samples' files (only the filename carries the sample tag), so summing across both
// file lists is all "both" needs -- the same operation apply_unfolding_2d.C performs.
//
// ⚠️ Physics caveat, inherited: "both" adds a b-enriched sample onto inclusive QCD, which
// double-counts b jets unless the sample weights already handle the overlap. The bjet
// sample is there for TEMPLATE STATISTICS. Comparing qcd against both shows how much it
// moves.
std::vector<TString> samplesIn(const TString &sample)
{
    if (sample == "both") return {"qcd", "bjet"};
    return {sample};
}

// prodPre (2026-11-09): a production tag that sits BEFORE the MC-SF tags in the filename, as
// run_agg_ntuple_chunks.sh writes it -- "_trkdrop030" for the 3% track-drop MC. Empty = nominal.
std::vector<TString> mcFiles(const TString &sample, const TString &generator,
                             const TString &tag, bool rmatrix, const TString &prodPre = "")
{
    std::vector<TString> files;
    for (const TString &s : samplesIn(sample)) {
        const TString subdir = sampleSubdir(s, generator);
        if (subdir.Length() == 0) return {};
        int found = 0;
        for (int b = 0; b < 50; ++b) {
            TString f = TString::Format(
                "%s/%s/agg_ntuple_chunks/block_%04d/%sRun3_%s_template_for_fit_histos_3D_%s_f%s%s_%s.root",
                kBase, subdir.Data(), b, rmatrix ? "RMatrix_" : "", kBtagTag, s.Data(),
                rmatrix ? "" : "MCGEN",
                // The SV b SF production (svsfTag(), 2026-10-12): MC only -- the same tag
                // also names the data files (dataFiles below), which carry no SF.
                (prodPre + svsfTag() + trackBdtSfTag()).Data(), tag.Data());
            if (gSystem->AccessPathName(f) == kFALSE) { files.push_back(f); ++found; }
        }
        if (found == 0)
            std::cerr << "WARNING: no " << (rmatrix ? "RMatrix" : "template")
                      << " files for sample '" << s << "' under " << subdir
                      << " with tag '" << tag << "'" << std::endl;
    }
    return files;
}

// Data is split over 5 primary datasets x 10 blocks.
std::vector<TString> dataFiles(const TString &tag)
{
    std::vector<TString> files;
    for (int pd = 0; pd < 5; ++pd)
        for (int b = 0; b < 10; ++b) {
            TString f = TString::Format(
                "%s/HardProbes/agg_template_chunks/HardProbes%d/block_%04d/"
                "Run3_%s_template_for_fit_histos_3D_data_fMCGEN_%s.root",
                kBase, pd, b, kBtagTag, tag.Data());
            if (gSystem->AccessPathName(f) == kFALSE) files.push_back(f);
        }
    return files;
}

// ---- jet-pT slicing -------------------------------------------------------------------
// jtpt_binsVector = {80, 100, 120} gives two bins, and they are NOT symmetric:
//   bin 1  80 < pT < 100     the lower bin
//   bin 2  pT > 100          OPEN-ENDED -- jtpt_fill() folds every jet above 120 into it,
//                            so it is "100 to infinity", not "100 to 120"
// The EEC results are quoted per pT bin (its plots carry _pt1/_pt2), so these do too.
// ptbin 0 means integrated over the whole range, 80 -> inf.
//
// There is no ROOT underflow/overflow to worry about on this axis: the pT cut in
// passRecoJetKinematics keeps everything below 80 out, and jtpt_fill() folds the high side
// into bin 2, so bins 1..N hold every filled jet (verified: both are exactly 0).
int ptLo(int ptbin) { return ptbin == 0 ? 1 : ptbin; }
int ptHi(int ptbin, TH1 *h) { return ptbin == 0 ? h->GetNbinsZ() : ptbin; }
// TH2D(observable, pT): the pT axis is Y, not Z.
int ptHiY(int ptbin, TH1 *h) { return ptbin == 0 ? h->GetNbinsY() : ptbin; }

TString ptLabel(int ptbin)
{
    if (ptbin == 0) return Form("p_{T} > %g GeV (all)", jtpt_binsVector[0]);
    // Open-ended last bin -> quote it as a threshold, the same rule pt_label() uses.
    if (ptbin == jtpt_bins) return Form("p_{T} > %g GeV", jtpt_binsVector[ptbin - 1]);
    return Form("%g < p_{T} < %g GeV", jtpt_binsVector[ptbin - 1], jtpt_binsVector[ptbin]);
}

TString ptSuffix(int ptbin) { return ptbin == 0 ? "_ptall" : Form("_pt%d", ptbin); }

// Project a TH3D(m_2B, X, pT) onto its X (observable) axis, over all m_2B, for one pT bin.
TH1D *projX(TH3D *h, const char *name, int ptbin)
{
    if (!h) return nullptr;
    TH1D *p = (TH1D *) h->ProjectionY(name, 1, h->GetNbinsX(), ptLo(ptbin), ptHi(ptbin, h));
    p->SetDirectory(nullptr);
    return p;
}

// Project onto m_2B for one bin of the observable axis, for one pT bin.
TH1D *projMB(TH3D *h, int ix, const char *name, int ptbin)
{
    if (!h) return nullptr;
    TH1D *p = (TH1D *) h->ProjectionX(name, ix, ix, ptLo(ptbin), ptHi(ptbin, h));
    p->SetDirectory(nullptr);
    return p;
}

void styleShape(TH1 *h, int col, int style = 1)
{
    if (!h) return;
    h->SetLineColor(col);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
    h->SetMarkerColor(col);
    h->SetStats(0);
}

void unitArea(TH1 *h)
{
    if (h && h->Integral() > 0) h->Scale(1. / h->Integral());
}

} // anonymous namespace


static void z_first_look_one(const TString &sample, const TString &generator,
                             const TString &tag, int ptbin,
                             const TString &observable)
{
    // A "_noeecw" production was filled with the EEC weight OFF, so h3D_* holds yields
    // (dN/dB), not an energy-weighted correlator. The labels have to say which, or the two
    // sets of plots are indistinguishable once they are out of their folders.
    const bool eecOff  = tag.Contains("noeecw");
    const TString wLab = eecOff ? "yields" : "EEC weighted";
    // Reject a bad sample/generator string rather than silently falling back to a default,
    // the same rule apply_unfolding_2d.C follows.
    if (sample != "qcd" && sample != "bjet" && sample != "both") {
        std::cerr << "ERROR: sample must be qcd | bjet | both (got '" << sample << "')" << std::endl;
        return;
    }
    if (generator != "pythia" && generator != "herwig") {
        std::cerr << "ERROR: generator must be pythia | herwig (got '" << generator << "')" << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Any observable, not just the balance: obsByName() supplies the binning, the axis title and the
    // "_<obs>" histogram-name suffix, so this macro reads whichever axis it is asked for.
    // Its own name is historical -- the balance was the first observable that needed a look before the
    // fit existed, and every other one needs the same four questions answered.
    const ObsDef obsd = obsByName(observable);
    if (obsd.nbins == 0) {
        std::cerr << "ERROR: unknown observable '" << observable << "' (use dr | B)"
                  << std::endl;
        return;
    }
    // The production tag and the observable have to agree, or this draws one observable's
    // histograms out of another observable's production -- which is exactly the mistake
    // obsProdTag() exists to prevent. It returns "_upartv2_3obs"; the tag here is an OUT_TAG
    // with no leading underscore and possibly a "noeecw_" in front, so compare on the stem.
    {
        TString stem = obsProdTag(observable);
        if (stem.BeginsWith("_")) stem.Remove(0, 1);
        // 2026-10-28: the "<era>ptw" productions (TF_METHOD=chi2scan) are the same observable's
        TString stemPtw = stem; stemPtw.ReplaceAll(eraTag()(1, eraTag().Length() - 1), eraTag()(1, eraTag().Length() - 1) + "ptw");
        // if (!tag.Contains(stem)) {   // until 2026-10-28
        if (!tag.Contains(stem) && !tag.Contains(stemPtw)) {
            std::cerr << "ERROR: tag '" << tag << "' is not a production of observable '"
                      << observable << "' (expected it to contain '" << stem << "'). "
                      << "Refusing to read one observable's histograms from another's "
                      << "production." << std::endl;
            return;
        }
    }

    std::cout << "\n================ " << observable << " first look: " << sample << " " << generator
              << ", tag '" << tag << "', " << ptLabel(ptbin)
              << " ================" << std::endl;

    std::vector<TString> f_tmpl = mcFiles(sample, generator, tag, false);
    std::vector<TString> f_rm   = mcFiles(sample, generator, tag, true);
    std::vector<TString> f_data = dataFiles(tag);
    std::cout << "  MC template files: " << f_tmpl.size()
              << ", MC RMatrix files: " << f_rm.size()
              << ", data files: " << f_data.size() << std::endl;
    if (f_tmpl.empty()) {
        std::cerr << "ERROR: no MC template files for tag '" << tag << "'. "
                  << "Run:  OUT_TAG=" << tag << " SAMPLE=" << sample
                  << " GENERATOR=" << generator << " ./run_agg_ntuple_chunks.sh" << std::endl;
        return;
    }

    // ---- MC flavour templates, in B and in dR for reference ----
    TH3D *h3_2b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_bb"), "s_2b_B");
    TH3D *h3_1b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_b"),  "s_1b_B");
    TH3D *h3_0b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h3D_0b"), "s_0b_B");
    if (!h3_2b) return;
    TH3D *h3_2b_dr = sumOverBlocks<TH3D>(f_tmpl, "h3D_bb", "s_2b_dr");

    TH3D *h3_data = f_data.empty() ? nullptr
                                   : sumOverBlocks<TH3D>(f_data, obsd.n("h3D_data"), "s_data_B");
    if (!h3_data)
        std::cout << "  NOTE: no data for this tag yet -- MC-only plots will be made."
                  << std::endl;

    // "<obs>_first_look_...": the folder is named from the observable, so the balance now
    // writes into "B_first_look_..." where it used to write "z_first_look_...". The old
    // folders are left where they are rather than moved.
    TString outdir = TString::Format("%s/results/%s_first_look_%s_%s%s_%s",
                                     kBase, observable.Data(), sample.Data(),
                                     generator.Data(), (svsfTag() + trackBdtSfTag()).Data(), tag.Data());
    gSystem->mkdir(outdir, true);

    // ---- 1. Do the flavour templates separate in B? ----
    {
        TH1D *p2b = projX(h3_2b, "p2b", ptbin), *p1b = projX(h3_1b, "p1b", ptbin), *p0b = projX(h3_0b, "p0b", ptbin);
        styleShape(p2b, kCol2b); styleShape(p1b, kCol1b); styleShape(p0b, kCol0b);

        std::cout << "\n--- 1. " << (eecOff ? "YIELD" : "EEC-weighted")
                  << " " << observable << " shapes by truth flavour (unit area) ---" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "    " << observable << " bin      |    2b    |    1b    |    0b" << std::endl;
        std::cout << "    -------------+----------+----------+---------" << std::endl;
        TH1D *n2b = (TH1D *) p2b->Clone("n2b"); unitArea(n2b);
        TH1D *n1b = p1b ? (TH1D *) p1b->Clone("n1b") : nullptr; unitArea(n1b);
        TH1D *n0b = p0b ? (TH1D *) p0b->Clone("n0b") : nullptr; unitArea(n0b);
        for (int i = 1; i <= n2b->GetNbinsX(); ++i)
            std::cout << "    " << std::setw(5) << n2b->GetXaxis()->GetBinLowEdge(i)
                      << " - " << std::setw(5) << n2b->GetXaxis()->GetBinUpEdge(i)
                      << " | " << std::setw(8) << n2b->GetBinContent(i)
                      << " | " << std::setw(8) << (n1b ? n1b->GetBinContent(i) : 0.)
                      << " | " << std::setw(8) << (n0b ? n0b->GetBinContent(i) : 0.) << std::endl;

        TCanvas c("c1", "", 900, 700);
        n2b->SetTitle(Form(";%s;normalised", obsd.axis.Data()));
        double mx = n2b->GetMaximum();
        if (n1b) mx = std::max(mx, n1b->GetMaximum());
        if (n0b) mx = std::max(mx, n0b->GetMaximum());
        n2b->SetMaximum(1.8 * mx);   // 1.8, not 1.45: the info lines now sit inside the frame
        n2b->SetMinimum(0.);
        n2b->Draw("hist");
        if (n1b) n1b->Draw("hist same");
        if (n0b) n0b->Draw("hist same");
        TLegend leg(0.16, 0.72, 0.60, 0.89);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(n2b, "true 2b (signal)", "l");
        if (n1b) leg.AddEntry(n1b, "1b", "l");
        if (n0b) leg.AddEntry(n0b, "0b (mistag)", "l");
        leg.Draw();
        // Inside the frame, under the legend: the space above the frame is the CMS header's.
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030);
        lat.DrawLatex(0.16, 0.68, Form("%s %s, reco level", sample.Data(), generator.Data()));
        lat.DrawLatex(0.16, 0.63, Form("%s, %s", ptLabel(ptbin).Data(), wLab.Data()));
        drawCMSHeader(&c);
        c.SaveAs(outdir + "/1_" + observable + "_templates_by_flavour" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/1_" + observable + "_templates_by_flavour" + ptSuffix(ptbin) + ".png");
    }

    // ---- 2. m_2B inside each B bin: what the template fit will actually fit ----
    {
        const int nz = h3_2b->GetNbinsY();
        TCanvas c("c2", "", 1300, 380);
        c.Divide(nz, 1, 0.001, 0.001);
        std::cout << Form("\n--- 2. m_2B separation inside each %s bin ---", observable.Data()) << std::endl;
        std::cout << "    " << observable << " bin      |   N(2b)   |   N(1b)   |   N(0b)   | S/(S+B)" << std::endl;
        std::cout << "    -------------+-----------+-----------+-----------+--------" << std::endl;
        // Keep the projections alive until the canvas is written.
        std::vector<TH1D *> keep;
        for (int iz = 1; iz <= nz; ++iz) {
            TH1D *m2b = projMB(h3_2b, iz, Form("m2b_%d", iz), ptbin);
            TH1D *m1b = projMB(h3_1b, iz, Form("m1b_%d", iz), ptbin);
            TH1D *m0b = projMB(h3_0b, iz, Form("m0b_%d", iz), ptbin);
            const double n2 = m2b ? m2b->Integral() : 0.;
            const double n1 = m1b ? m1b->Integral() : 0.;
            const double n0 = m0b ? m0b->Integral() : 0.;
            std::cout << "    " << std::setw(5) << h3_2b->GetYaxis()->GetBinLowEdge(iz)
                      << " - " << std::setw(5) << h3_2b->GetYaxis()->GetBinUpEdge(iz)
                      << " | " << std::setw(9) << n2 << " | " << std::setw(9) << n1
                      << " | " << std::setw(9) << n0
                      << " | " << std::setw(7)
                      << ((n2 + n1 + n0) > 0 ? n2 / (n2 + n1 + n0) : 0.) << std::endl;

            styleShape(m2b, kCol2b); styleShape(m1b, kCol1b); styleShape(m0b, kCol0b);
            unitArea(m2b); unitArea(m1b); unitArea(m0b);
            c.cd(iz);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14);
            double mx = m2b ? m2b->GetMaximum() : 1.;
            if (m1b) mx = std::max(mx, m1b->GetMaximum());
            if (m0b) mx = std::max(mx, m0b->GetMaximum());
            if (m2b) {
                m2b->SetTitle(";m_{2B} [GeV];normalised");
                m2b->SetMaximum(1.9 * mx); m2b->SetMinimum(0.);
                m2b->GetXaxis()->SetTitleSize(0.055); m2b->GetYaxis()->SetTitleSize(0.055);
                m2b->Draw("hist");
            }
            if (m1b) m1b->Draw("hist same");
            if (m0b) m0b->Draw("hist same");
            TLatex *lt = new TLatex(); lt->SetNDC(); lt->SetTextSize(0.055);
            // Name the observable, not "B": these panels are the fit's input and end up in
            // talks on their own, where "0.30 < B < 0.50" on a dr plot would be simply wrong.
            // 0.84, inside the frame: the margin above it carries the CMS header in panel 1.
            lt->DrawLatex(0.20, 0.84, Form("%.2f < %s < %.2f",
                                           h3_2b->GetYaxis()->GetBinLowEdge(iz),
                                           obsSymbol(observable).Data(),
                                           h3_2b->GetYaxis()->GetBinUpEdge(iz)));
            if (iz == 1) {
                TLegend *lg = new TLegend(0.45, 0.52, 0.92, 0.78);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
                lg->AddEntry(m2b, "2b", "l");
                if (m1b) lg->AddEntry(m1b, "1b", "l");
                if (m0b) lg->AddEntry(m0b, "0b", "l");
                lg->Draw();
            }
            // Panels are too narrow for the beam label, so the CMS mark alone, in panel 1.
            if (iz == 1) drawCMSHeader((TPad *) gPad, "Work In Progress", "");
            keep.push_back(m2b); keep.push_back(m1b); keep.push_back(m0b);
        }
        c.SaveAs(outdir + "/2_mB_in_" + observable + "_bins" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/2_mB_in_" + observable + "_bins" + ptSuffix(ptbin) + ".png");
    }

    // ---- 3. Response quality in B: purity, efficiency, migration ----
    if (!f_rm.empty()) {
        TH2D *pur_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_purity_numerator_tf"), "pn_B");
        TH2D *pud_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_purity_denominator_tf"), "pd_B");
        TH2D *efn_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_efficiency_numerator_tf"), "en_B");
        TH2D *efd_B = sumOverBlocks<TH2D>(f_rm, obsd.n("h_full_efficiency_denominator_tf"), "ed_B");
        // dR reference, same quantities on its own binning.
        TH2D *pur_d = sumOverBlocks<TH2D>(f_rm, "h_full_purity_numerator_tf", "pn_d");
        TH2D *pud_d = sumOverBlocks<TH2D>(f_rm, "h_full_purity_denominator_tf", "pd_d");

        // ⚠️ Ratios are RECOMPUTED from the summed numerators and denominators, never
        // summed themselves -- adding N stored ratios gives N x the true value. Same rule
        // as apply_unfolding_2d.C's "How 'both' is summed".
        auto ratio1D = [&](TH2D *num, TH2D *den, const char *name) -> TH1D * {
            if (!num || !den) return nullptr;
            TH1D *n = num->ProjectionX(Form("%s_n", name), ptLo(ptbin), ptHiY(ptbin, num));
            TH1D *d = den->ProjectionX(Form("%s_d", name), ptLo(ptbin), ptHiY(ptbin, den));
            TH1D *r = (TH1D *) n->Clone(name);
            r->Divide(n, d, 1., 1., "b");
            r->SetDirectory(nullptr);
            return r;
        };
        TH1D *purity_B = ratio1D(pur_B, pud_B, "purity_B");
        TH1D *eff_B    = ratio1D(efn_B, efd_B, "eff_B");
        TH1D *purity_d = ratio1D(pur_d, pud_d, "purity_d");

        std::cout << Form("\n--- 3. Response quality in %s (dR's own values for scale) ---", observable.Data()) << std::endl;
        if (purity_B) {
            std::cout << "    " << observable << " bin      |  purity  | efficiency" << std::endl;
            std::cout << "    -------------+----------+-----------" << std::endl;
            for (int i = 1; i <= purity_B->GetNbinsX(); ++i)
                std::cout << "    " << std::setw(5) << purity_B->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << purity_B->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(8) << purity_B->GetBinContent(i)
                          << " | " << std::setw(10) << (eff_B ? eff_B->GetBinContent(i) : 0.)
                          << std::endl;
        }
        if (purity_d) {
            std::cout << "    dR purity, for scale:";
            for (int i = 1; i <= purity_d->GetNbinsX(); ++i)
                std::cout << " " << std::setprecision(3) << purity_d->GetBinContent(i);
            std::cout << std::endl;
        }

        if (purity_B) {
            TCanvas c("c3", "", 900, 700);
            styleShape(purity_B, kCol2b);
            styleShape(eff_B, kCol1b);
            purity_B->SetTitle(Form(";%s;fraction", obsd.axis.Data()));
            purity_B->SetMinimum(0.); purity_B->SetMaximum(1.35);
            purity_B->Draw("hist");
            if (eff_B) eff_B->Draw("hist same");
            TLegend leg(0.16, 0.74, 0.62, 0.89);
            leg.SetBorderSize(0); leg.SetFillStyle(0);
            leg.AddEntry(purity_B, "purity (reco-binned)", "l");
            if (eff_B) leg.AddEntry(eff_B, "efficiency (gen-binned)", "l");
            leg.Draw();
            TLine l(obsd.bins[0], 1., obsd.bins[obsd.nbins], 1.);
            l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.Draw();
            drawCMSHeader(&c);
            c.SaveAs(outdir + "/3_" + observable + "_purity_efficiency" + ptSuffix(ptbin) + ".pdf");
            c.SaveAs(outdir + "/3_" + observable + "_purity_efficiency" + ptSuffix(ptbin) + ".png");
        }
    }

    // ---- 4. Data vs MC in B, reco level, shapes only ----
    if (h3_data) {
        TH1D *pdata = projX(h3_data, "pdata", ptbin);
        TH1D *pmc   = projX(h3_2b, "pmc_sum", ptbin);
        if (h3_1b) pmc->Add(projX(h3_1b, "pmc_1b", ptbin));
        if (h3_0b) pmc->Add(projX(h3_0b, "pmc_0b", ptbin));
        styleShape(pdata, kColData); styleShape(pmc, kCol2b);
        pdata->SetMarkerStyle(20);

        std::cout << Form("\n--- 4. Data vs MC in %s (reco level, unit area, NO fit/corrections) ---", observable.Data())
                  << std::endl;
        std::cout << "    " << observable << " bin      |    data  |    MC    |  data/MC" << std::endl;
        std::cout << "    -------------+----------+----------+---------" << std::endl;
        unitArea(pdata); unitArea(pmc);
        for (int i = 1; i <= pdata->GetNbinsX(); ++i) {
            const double d = pdata->GetBinContent(i), m = pmc->GetBinContent(i);
            std::cout << "    " << std::setw(5) << pdata->GetXaxis()->GetBinLowEdge(i)
                      << " - " << std::setw(5) << pdata->GetXaxis()->GetBinUpEdge(i)
                      << " | " << std::setw(8) << d << " | " << std::setw(8) << m
                      << " | " << std::setw(8) << (m > 0 ? d / m : 0.) << std::endl;
        }

        TCanvas c("c4", "", 900, 700);
        pmc->SetTitle(Form(";%s;normalised", obsd.axis.Data()));
        pmc->SetMinimum(0.);
        pmc->SetMaximum(1.8 * std::max(pmc->GetMaximum(), pdata->GetMaximum()));
        pmc->Draw("hist");
        pdata->Draw("e same");
        TLegend leg(0.16, 0.74, 0.65, 0.89);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(pdata, "data (2b-tagged, pre-fit)", "lep");
        leg.AddEntry(pmc, Form("%s MC, all flavours", sample.Data()), "l");
        leg.Draw();
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030);
        lat.DrawLatex(0.16, 0.69, Form("%s, reco level, %s", ptLabel(ptbin).Data(), wLab.Data()));
        lat.DrawLatex(0.16, 0.64, "no fit / unfolding / UParT SF");
        drawCMSHeader(&c);
        c.SaveAs(outdir + "/4_" + observable + "_data_vs_mc" + ptSuffix(ptbin) + ".pdf");
        c.SaveAs(outdir + "/4_" + observable + "_data_vs_mc" + ptSuffix(ptbin) + ".png");
    }

    // ---- 5. RAW YIELDS per B bin ----------------------------------------------------
    // The h_count_* histograms are the yields: same jets, same selection, same binning as
    // the h3D_* templates, but filled with the tree weight ONLY -- no EEC weight. That is
    // the difference between "how many jets" and "how much EEC", and it is why none of the
    // plots above is a yield (they are all EEC-weighted or unit-area normalised).
    //
    // For DATA the tree weight is 1, so h_count_data IS the raw jet count and its error is
    // sqrt(N). ⚠️ The trigger prescale is NOT in it -- the prescale multiplies `eec`, not
    // the count -- so this is the number of jets actually recorded, not a luminosity-scaled
    // yield. For MC the tree weight is the sample weight, so the count is a weighted
    // prediction; GetEntries() is printed next to it as the raw statistics behind it.
    {
        std::cout << Form("\n--- 5. RAW YIELDS per %s bin (h_count_*, NO EEC weight) ---", observable.Data()) << std::endl;

        TH3D *c_data = f_data.empty() ? nullptr
                                      : sumOverBlocks<TH3D>(f_data, obsd.n("h_count_data"), "c_data_B");
        TH3D *c_2b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_bb"), "c_2b_B");
        TH3D *c_1b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_b"),  "c_1b_B");
        TH3D *c_0b = sumOverBlocks<TH3D>(f_tmpl, obsd.n("h_count_0b"), "c_0b_B");

        TH1D *y_data = projX(c_data, "y_data", ptbin);
        TH1D *y_2b   = projX(c_2b,   "y_2b", ptbin);
        TH1D *y_1b   = projX(c_1b,   "y_1b", ptbin);
        TH1D *y_0b   = projX(c_0b,   "y_0b", ptbin);

        std::cout << std::fixed;
        if (y_data) {
            std::cout << "\n  DATA -- raw b-tagged 2-SV jet counts (prescale NOT applied):\n";
            std::cout << "    " << observable << " bin      |      jets  |  sqrt(N)  | rel. stat" << std::endl;
            std::cout << "    -------------+------------+-----------+----------" << std::endl;
            double tot = 0.;
            for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
                const double N = y_data->GetBinContent(i);
                tot += N;
                std::cout << "    " << std::setprecision(3)
                          << std::setw(5) << y_data->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << y_data->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(10) << std::setprecision(0) << N
                          << " | " << std::setw(9) << std::setprecision(1) << std::sqrt(N)
                          << " | " << std::setw(8) << std::setprecision(4)
                          << (N > 0 ? std::sqrt(N) / N : 0.) << std::endl;
            }
            std::cout << "    total: " << std::setprecision(0) << tot << " jets" << std::endl;
        }

        if (y_2b) {
            std::cout << "\n  MC (" << sample << " " << generator
                      << ") -- weighted counts, with raw entries behind them:\n";
            std::cout << "    " << observable << " bin      |   2b (wt) |   1b (wt) |   0b (wt) | S/(S+B) | 2b eff. entries"
                      << std::endl;
            std::cout << "    -------------+-----------+-----------+-----------+---------+----------------"
                      << std::endl;
            // Effective entries per bin, (sum w)^2 / sum w^2, from the content and its
            // error. This is the statistical power actually behind a weighted bin --
            // GetEntries() on the projection would give the TOTAL over all bins, the same
            // number four times over, which is why it is not used here.
            auto effEntries = [](TH1D *h, int i) -> double {
                const double c = h->GetBinContent(i), e = h->GetBinError(i);
                return (e > 0.) ? (c * c) / (e * e) : 0.;
            };
            for (int i = 1; i <= y_2b->GetNbinsX(); ++i) {
                const double n2 = y_2b->GetBinContent(i);
                const double n1 = y_1b ? y_1b->GetBinContent(i) : 0.;
                const double n0 = y_0b ? y_0b->GetBinContent(i) : 0.;
                std::cout << "    " << std::setprecision(3)
                          << std::setw(5) << y_2b->GetXaxis()->GetBinLowEdge(i)
                          << " - " << std::setw(5) << y_2b->GetXaxis()->GetBinUpEdge(i)
                          << " | " << std::setw(9) << std::setprecision(4) << n2
                          << " | " << std::setw(9) << n1
                          << " | " << std::setw(9) << n0
                          << " | " << std::setw(7)
                          << ((n2 + n1 + n0) > 0 ? n2 / (n2 + n1 + n0) : 0.)
                          << " | " << std::setw(15) << std::setprecision(0)
                          << effEntries(y_2b, i) << std::endl;
            }
            std::cout << "\n  ⚠️ MC counts are WEIGHTED (sample weight), so they are a prediction,"
                      << " not a number of\n     simulated jets, and they are NOT normalised to the"
                      << " data luminosity -- do not\n     compare the MC and data columns directly."
                      << std::endl;
        }

        // Yields plot: the three MC flavour components stacked, with data on top.
        //
        // Colours and styles come from Help_Functions.h -- TFColor::c2B/c1B/c0B and
        // styleTemplate/styleData -- i.e. exactly what the template fit draws its own plots
        // with: 2B red, 1B blue, 0B green, data black closed circles. Someone flipping
        // between this plot and the fit's m_2B plots is looking at the same categories, so
        // they have to be the same colours.
        //
        // Stack order matches the fit's: 2B at the bottom, then 1B, 0B on top.
        //
        // ⚠️ MC IS SCALED TO THE DATA. The MC counts are weighted predictions and are NOT
        // normalised to the data luminosity, so the raw MC and data numbers are not
        // comparable (the printout above says so too). Scaling the SUM to the data integral,
        // with the three components keeping their relative proportions, is the pre-fit
        // normalisation the template fit itself starts from -- what the plot then shows is
        // whether the SHAPES agree, which is the question worth asking before the fit.
        if (y_data) {
            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 14. * font_scale;

            const double n_jets = y_data->Integral();

            // Scale the MC sum to the data, preserving the flavour proportions.
            double mc_sum = 0.;
            for (TH1D *h : {y_2b, y_1b, y_0b}) if (h) mc_sum += h->Integral();
            const double mc_scale = (mc_sum > 0.) ? n_jets / mc_sum : 1.;
            for (TH1D *h : {y_2b, y_1b, y_0b}) if (h) h->Scale(mc_scale);

            std::cout << "    MC scaled to data by x" << mc_scale
                      << " (sum of 2b+1b+0b normalised to the data integral)" << std::endl;

            TCanvas c("c5", "", 800, 700);
            c.SetTicks(1, 1);
            c.SetLeftMargin(0.15); c.SetRightMargin(0.05);
            c.SetTopMargin(0.09);  c.SetBottomMargin(0.15);

            if (y_2b) styleTemplate(y_2b, TFColor::c2B());
            if (y_1b) styleTemplate(y_1b, TFColor::c1B());
            if (y_0b) styleTemplate(y_0b, TFColor::c0B());
            styleData(y_data);

            THStack st("st_raw_yields", "");
            if (y_2b) st.Add(y_2b);      // bottom
            if (y_1b) st.Add(y_1b);
            if (y_0b) st.Add(y_0b);      // top

            // Frame from the data histogram, so the axis is the observable's own binning and
            // the range covers whichever of data/stack is taller.
            TH1D *frame = (TH1D *) y_data->Clone("h_raw_frame");
            frame->SetDirectory(nullptr);
            frame->Reset();
            frame->SetTitle("");
            frame->SetMinimum(0.);
            frame->SetMaximum(1.45 * std::max(y_data->GetMaximum(), st.GetMaximum()));
            frame->GetYaxis()->SetTitle("2b jets per bin");
            frame->GetXaxis()->SetTitle(obsd.axis);
            for (TAxis *ax : { frame->GetXaxis(), frame->GetYaxis() }) {
                ax->CenterTitle(true);
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            frame->GetXaxis()->SetTitleOffset(1.15);
            frame->GetYaxis()->SetTitleOffset(1.9);
            frame->GetYaxis()->SetNdivisions(505);
            // ROOT pads an axis it was not given an explicit range for -- B would be drawn
            // out to 1.1, past where the observable can go. Same fix as the other macros.
            if (observable != "dr") frame->GetXaxis()->SetRange(1, frame->GetNbinsX());
            frame->Draw("AXIS");
            st.Draw("hist same");
            y_data->Draw("PE X0 same");

            TLatex cms; cms.SetNDC();
            cms.SetTextFont(62); cms.SetTextSize(0.045); cms.DrawLatex(0.15, 0.945, "CMS");
            cms.SetTextFont(52); cms.SetTextSize(0.036); cms.DrawLatex(0.255, 0.945, "Work In Progress");
            cms.SetTextFont(42); cms.SetTextSize(0.036); cms.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // Two lines only, by request: which pT bin, and what the plot is. The jet count,
            // the MC scale factor and the "no fit, no corrections" caveat were on the plot
            // and are now only in the printout above -- they have not stopped being true.
            // ⚠️ In particular the MC here IS scaled to the data (mc_scale, printed below),
            // so nothing on the canvas says so any more. Say it in the caption.
            TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.032);
            note.SetTextAlign(13);
            double y_note = 0.87;
            note.DrawLatex(0.19, y_note, ptLabel(ptbin));
            note.DrawLatex(0.19, y_note -= 0.048, "Raw yields");

            TLegend leg(0.62, 0.66, 0.94, 0.87);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetMargin(0.25);
            leg.SetTextFont(font_code); leg.SetTextSize(legend_size);
            leg.AddEntry(y_data, "Data", "pe1");
            if (y_2b) leg.AddEntry(y_2b, "MC 2B (signal)", "f");
            if (y_1b) leg.AddEntry(y_1b, "MC 1B", "f");
            if (y_0b) leg.AddEntry(y_0b, "MC 0B (mistag)", "f");
            leg.Draw();
            c.RedrawAxis();

            c.SaveAs(outdir + "/5_" + observable + "_raw_yields" + ptSuffix(ptbin) + ".pdf");
            c.SaveAs(outdir + "/5_" + observable + "_raw_yields" + ptSuffix(ptbin) + ".png");
        }
    }

    std::cout << "\nWritten to " << outdir << std::endl;
    std::cout << "(/data_CMS is not on the sshfs mount -- use sync_plots_to_eos.sh or scp)"
              << std::endl;
}

// Entry point: one pass per jet-pT slice.
//
// The EEC results are quoted per pT bin, so these are too. Two passes:
//   _pt1    80 < pT < 100   the lower bin
//   _pt2    pT > 100        the main bin, OPEN-ENDED (jtpt_fill folds >120 into it)
//
// The pT-INTEGRATED pass is deliberately NOT produced. The two pT bins have visibly
// different signal fractions (0.33-0.37 against 0.47-0.48 from the B template fit), so a
// pT-integrated shape is a blend of two populations rather than a measurement of either,
// and having it sitting in the same folder as the real ones only invites it being read as
// a result. z_first_look_one() still accepts ptbin = 0 for it if it is ever wanted for a
// statistics check -- it is just never called that way:
//
//   /* ---- disabled (kept for reference): the pT-integrated pass ----
//      z_first_look_one(sample, generator, tag, 0);   // writes *_ptall.*
//      ---- */
void plot_z_first_look(const TString &sample    = "qcd",
                       const TString &generator = "pythia",
                       const TString &tag       = "upartv2_B",
                       const TString &observable = "B")
{
    for (int ptbin = 1; ptbin <= jtpt_bins; ++ptbin)
        z_first_look_one(sample, generator, tag, ptbin, observable);
}


// ============================================================================
// RAW YIELDS: Pythia vs Herwig (added 2026-09-29)
// ============================================================================
// Plot 5 of the first look (the flavour-stacked raw yields) for BOTH generators on one
// canvas. Same yields as plot 5 -- the h_count_* histograms, tree weight only, no EEC
// weight, summed over all m_2B -- and the same normalisation: each generator's 2B+1B+0B
// total is scaled to the data integral, flavour proportions kept. So what differs between
// the two curves is the SHAPE of the reco-level yield, and (dashed) how much of it each
// generator calls signal.
//
// Top pad : data, the Pythia and Herwig totals (solid), their 2B components (dashed).
// Lower   : each total over the data, with the data's relative stat error as a grey band.
//
// Colours: orange (Pythia) vs purple (Herwig) -- TFColor::cmpA/cmpB, the pair
// Help_Functions.h reserves for exactly two variants of one quantity. Data black.
//
//   root -l -b -q -e '.L template_fit.cpp+' \
//        -e 'plot_raw_yields_generators("both","upartv2rho","dr")'
//   root -l -b -q -e '.L template_fit.cpp+' \
//        -e 'plot_raw_yields_generators("both","noeecw_upartv2rho_B","B")'
static void raw_yields_generators_one(const TString &sample, const TString &tag, int ptbin,
                                      const TString &observable)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const ObsDef obsd = obsByName(observable);
    if (obsd.nbins == 0) {
        std::cerr << "ERROR: unknown observable '" << observable << "' (use dr | B)" << std::endl;
        return;
    }
    // Same guard as z_first_look_one(): never read one observable out of another's production.
    {
        TString stem = obsProdTag(observable);
        if (stem.BeginsWith("_")) stem.Remove(0, 1);
        // 2026-10-28: the "<era>ptw" productions (TF_METHOD=chi2scan) are the same observable's
        TString stemPtw = stem; stemPtw.ReplaceAll(eraTag()(1, eraTag().Length() - 1), eraTag()(1, eraTag().Length() - 1) + "ptw");
        // if (!tag.Contains(stem)) {   // until 2026-10-28
        if (!tag.Contains(stem) && !tag.Contains(stemPtw)) {
            std::cerr << "ERROR: tag '" << tag << "' is not a production of observable '"
                      << observable << "' (expected it to contain '" << stem << "')." << std::endl;
            return;
        }
    }

    std::vector<TString> f_data = dataFiles(tag);
    TH3D *c_data = f_data.empty() ? nullptr
                                  : sumOverBlocks<TH3D>(f_data, obsd.n("h_count_data"), "cg_data");
    TH1D *y_data = projX(c_data, "yg_data", ptbin);
    if (!y_data) {
        std::cerr << "ERROR: no data yields for tag '" << tag << "'" << std::endl;
        return;
    }
    const double n_jets = y_data->Integral();

    // One generator's yields: total (2B+1B+0B) and the 2B part, both scaled by the SAME
    // factor, so the total matches the data and 2B keeps its share of it.
    struct GenYields { TH1D *tot = nullptr; TH1D *sig = nullptr; double scale = 1.; };
    auto load = [&](const TString &gen) -> GenYields {
        GenYields g;
        std::vector<TString> f = mcFiles(sample, gen, tag, false);
        if (f.empty()) return g;
        TH1D *y2 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_bb"), "cg_2b_" + gen), "yg_2b_" + gen, ptbin);
        TH1D *y1 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_b"),  "cg_1b_" + gen), "yg_1b_" + gen, ptbin);
        TH1D *y0 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_0b"), "cg_0b_" + gen), "yg_0b_" + gen, ptbin);
        if (!y2) return g;
        g.sig = y2;
        g.tot = (TH1D *) y2->Clone("yg_tot_" + gen);
        g.tot->SetDirectory(nullptr);
        if (y1) g.tot->Add(y1);
        if (y0) g.tot->Add(y0);
        g.scale = (g.tot->Integral() > 0.) ? n_jets / g.tot->Integral() : 1.;
        g.tot->Scale(g.scale);
        g.sig->Scale(g.scale);
        return g;
    };
    GenYields py = load("pythia");
    GenYields hw = load("herwig");
    if (!py.tot || !hw.tot) {
        std::cerr << "ERROR: missing " << (!py.tot ? "Pythia" : "Herwig") << " yields for tag '"
                  << tag << "'" << std::endl;
        return;
    }

    // ---- printout ----
    std::cout << "\n================ raw yields, Pythia vs Herwig: " << observable << ", "
              << sample << ", tag '" << tag << "', " << ptLabel(ptbin) << " ================"
              << std::endl;
    std::cout << "  MC totals scaled to the data (" << std::setprecision(0) << std::fixed
              << n_jets << " jets): Pythia x" << std::setprecision(4) << py.scale
              << ", Herwig x" << hw.scale << std::endl;
    std::cout << "    " << observable << " bin      |    data  |  Pythia  |  Herwig  | Py/data | Hw/data | Hw/Py  | 2B frac Py | 2B frac Hw"
              << std::endl;
    std::cout << "    -------------+----------+----------+----------+---------+---------+--------+------------+-----------"
              << std::endl;
    for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
        const double d = y_data->GetBinContent(i);
        const double p = py.tot->GetBinContent(i), h = hw.tot->GetBinContent(i);
        std::cout << "    " << std::setprecision(3)
                  << std::setw(5) << y_data->GetXaxis()->GetBinLowEdge(i)
                  << " - " << std::setw(5) << y_data->GetXaxis()->GetBinUpEdge(i)
                  << " | " << std::setw(8) << std::setprecision(0) << d
                  << " | " << std::setw(8) << p
                  << " | " << std::setw(8) << h
                  << " | " << std::setw(7) << std::setprecision(3) << (d > 0 ? p / d : 0.)
                  << " | " << std::setw(7) << (d > 0 ? h / d : 0.)
                  << " | " << std::setw(6) << (p > 0 ? h / p : 0.)
                  << " | " << std::setw(10) << (p > 0 ? py.sig->GetBinContent(i) / p : 0.)
                  << " | " << std::setw(9) << (h > 0 ? hw.sig->GetBinContent(i) / h : 0.)
                  << std::endl;
    }

    // ---- plot ----
    const Float_t font_scale  = 1200. / 800.;
    const Style_t font_code   = 43;
    const Float_t label_size  = 15. * font_scale;
    const Float_t title_size  = 15. * font_scale;
    const Float_t legend_size = 13. * font_scale;
    const Color_t col_py = TFColor::cmpA();
    const Color_t col_hw = TFColor::cmpB();

    TCanvas c("c_raw_gen", "", 800, 850);
    TPad *top = new TPad("rg_top", "", 0., 0.32, 1., 1.);
    TPad *bot = new TPad("rg_bot", "", 0., 0.,   1., 0.32);
    for (TPad *q : {top, bot}) { q->SetTicks(1, 1); q->SetFillColor(0); }
    top->SetMargin(0.15, 0.05, 0.0,  0.09);
    bot->SetMargin(0.15, 0.05, 0.26, 0.0);   // as the data-vs-gen ratio pad, which fits its title
    c.cd(); top->Draw(); bot->Draw();

    auto styleLine = [](TH1D *h, Color_t col, Style_t ls) {
        h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(ls);
        h->SetMarkerSize(0); h->SetFillStyle(0);
    };
    styleLine(py.tot, col_py, kSolid);  styleLine(py.sig, col_py, kDashed);
    styleLine(hw.tot, col_hw, kSolid);  styleLine(hw.sig, col_hw, kDashed);
    styleData(y_data);

    auto styleAxes = [&](TH1D *f) {
        for (TAxis *ax : { f->GetXaxis(), f->GetYaxis() }) {
            ax->CenterTitle(true);
            ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
        }
        if (observable != "dr") f->GetXaxis()->SetRange(1, f->GetNbinsX());
    };

    top->cd();
    TH1D *frame = (TH1D *) y_data->Clone("rg_frame");
    frame->SetDirectory(nullptr); frame->Reset(); frame->SetTitle("");
    frame->SetMinimum(0.001);
    frame->SetMaximum(1.5 * std::max({y_data->GetMaximum(), py.tot->GetMaximum(), hw.tot->GetMaximum()}));
    frame->GetYaxis()->SetTitle("2b jets per bin");
    styleAxes(frame);
    frame->GetYaxis()->SetTitleOffset(1.9);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("AXIS");
    py.sig->Draw("hist same");
    hw.sig->Draw("hist same");
    py.tot->Draw("hist same");
    hw.tot->Draw("hist same");
    y_data->Draw("PE X0 same");

    TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.037); note.SetTextAlign(13);
    note.DrawLatex(0.19, 0.86, ptLabel(ptbin));
    note.DrawLatex(0.19, 0.80, "Raw yields, reco level");
    note.DrawLatex(0.19, 0.74, "MC scaled to data");

    TLegend leg(0.52, 0.56, 0.94, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetMargin(0.22);
    leg.SetTextFont(font_code); leg.SetTextSize(legend_size);
    leg.AddEntry(y_data, "Data", "pe1");
    leg.AddEntry(py.tot, "Pythia 8, all flavours", "l");
    leg.AddEntry(py.sig, "Pythia 8, 2B only", "l");
    leg.AddEntry(hw.tot, "Herwig, all flavours", "l");
    leg.AddEntry(hw.sig, "Herwig, 2B only", "l");
    leg.Draw();
    drawCMSHeader(top);
    top->RedrawAxis();

    bot->cd();
    TH1D *rframe = (TH1D *) y_data->Clone("rg_rframe");
    rframe->SetDirectory(nullptr); rframe->Reset(); rframe->SetTitle("");
    TH1D *r_py = (TH1D *) py.tot->Clone("rg_r_py"); r_py->SetDirectory(nullptr);
    TH1D *r_hw = (TH1D *) hw.tot->Clone("rg_r_hw"); r_hw->SetDirectory(nullptr);
    TH1D *r_st = (TH1D *) y_data->Clone("rg_r_stat"); r_st->SetDirectory(nullptr);
    double rdev = 0.;
    for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
        const double d = y_data->GetBinContent(i);
        r_py->SetBinContent(i, d > 0 ? py.tot->GetBinContent(i) / d : 0.); r_py->SetBinError(i, 0.);
        r_hw->SetBinContent(i, d > 0 ? hw.tot->GetBinContent(i) / d : 0.); r_hw->SetBinError(i, 0.);
        r_st->SetBinContent(i, 1.);
        r_st->SetBinError(i, d > 0 ? y_data->GetBinError(i) / d : 0.);
        rdev = std::max({rdev, std::fabs(r_py->GetBinContent(i) - 1.), std::fabs(r_hw->GetBinContent(i) - 1.)});
    }
    const double rhalf = std::max(0.25, 1.2 * rdev);
    rframe->SetMinimum(1. - rhalf); rframe->SetMaximum(1. + rhalf);
    rframe->GetYaxis()->SetTitle("MC / Data");
    rframe->GetXaxis()->SetTitle(obsd.axis);
    styleAxes(rframe);
    rframe->GetYaxis()->SetTitleOffset(1.9);
    rframe->GetXaxis()->SetTitleOffset(1.1);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->Draw("AXIS");
    r_st->SetFillColor(kGray); r_st->SetFillStyle(1001); r_st->SetLineWidth(0); r_st->SetMarkerSize(0);
    r_st->Draw("E2 same");
    TLine one(obsd.bins[0], 1., obsd.bins[obsd.nbins], 1.);
    one.SetLineStyle(2); one.SetLineColor(kBlack); one.Draw();
    r_py->Draw("hist same");
    r_hw->Draw("hist same");
    bot->RedrawAxis();

    const TString outdir = TString::Format("%s/results/%s_raw_yields_generators_%s%s_%s",
                                           kBase, observable.Data(), sample.Data(),
                                           (svsfTag() + trackBdtSfTag()).Data(), tag.Data());
    gSystem->mkdir(outdir, true);
    const TString stem = outdir + "/" + observable + "_raw_yields_pythia_vs_herwig" + ptSuffix(ptbin);
    c.SaveAs(stem + ".pdf");
    c.SaveAs(stem + ".png");
}

void plot_raw_yields_generators(const TString &sample     = "both",
                                const TString &tag        = "upartv2rho",
                                const TString &observable = "dr")
{
    for (int ptbin = 1; ptbin <= jtpt_bins; ++ptbin)
        raw_yields_generators_one(sample, tag, ptbin, observable);
}


// ============================================================================
// RAW YIELDS: with vs without the 3% track drop (added 2026-11-09)
// ============================================================================
// The tracking-efficiency systematic at the reco level, before any fit: the same Pythia
// h_count_* yields as plot_raw_yields_generators (tree weight, summed over m_2B), nominal
// against the TRACK_EFF_UNC=true production (_trkdrop030, tracking_efficiency_syst.h).
//
// NORMALISATION: BOTH MC curves are scaled by the NOMINAL's factor to the data, not each to
// the data. Dropping tracks loses SVs and b tags, so the track-drop MC has FEWER jets; one
// common factor keeps that loss visible instead of normalising it away. Data is the same
// file for both (the drop is MC only) and is drawn for reference.
//
// Top pad : data, nominal and 3% drop totals (solid) and their 2B parts (dashed).
// Lower   : 3% drop / nominal, all flavours (solid) and 2B (dashed).
// Colours : orange (nominal) vs purple (3% drop), TFColor::cmpA/cmpB as for Pythia/Herwig.
//
//   root -l -b -q -e '.L template_fit.cpp+' \
//        -e 'plot_raw_yields_trackeff("both","upartv2rhojecvetobenrptw","dr")' \
//        -e 'plot_raw_yields_trackeff("both","noeecw_upartv2rhojecvetobenrptw_B","B")'
static void raw_yields_trackeff_one(const TString &sample, const TString &tag, int ptbin,
                                    const TString &observable)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const ObsDef obsd = obsByName(observable);
    if (obsd.nbins == 0) {
        std::cerr << "ERROR: unknown observable '" << observable << "' (use dr | B)" << std::endl;
        return;
    }
    {   // same guard as raw_yields_generators_one()
        TString stem = obsProdTag(observable);
        if (stem.BeginsWith("_")) stem.Remove(0, 1);
        TString stemPtw = stem; stemPtw.ReplaceAll(eraTag()(1, eraTag().Length() - 1), eraTag()(1, eraTag().Length() - 1) + "ptw");
        if (!tag.Contains(stem) && !tag.Contains(stemPtw)) {
            std::cerr << "ERROR: tag '" << tag << "' is not a production of observable '"
                      << observable << "' (expected it to contain '" << stem << "')." << std::endl;
            return;
        }
    }

    std::vector<TString> f_data = dataFiles(tag);
    TH3D *c_data = f_data.empty() ? nullptr
                                  : sumOverBlocks<TH3D>(f_data, obsd.n("h_count_data"), "ct_data");
    TH1D *y_data = projX(c_data, "yt_data", ptbin);
    if (!y_data) { std::cerr << "ERROR: no data yields for tag '" << tag << "'" << std::endl; return; }
    const double n_jets = y_data->Integral();

    struct Yields { TH1D *tot = nullptr; TH1D *sig = nullptr; };
    auto load = [&](const TString &pre, const TString &nm) -> Yields {
        Yields g;
        std::vector<TString> f = mcFiles(sample, "pythia", tag, false, pre);
        if (f.empty()) return g;
        TH1D *y2 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_bb"), "ct_2b_" + nm), "yt_2b_" + nm, ptbin);
        TH1D *y1 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_b"),  "ct_1b_" + nm), "yt_1b_" + nm, ptbin);
        TH1D *y0 = projX(sumOverBlocks<TH3D>(f, obsd.n("h_count_0b"), "ct_0b_" + nm), "yt_0b_" + nm, ptbin);
        if (!y2) return g;
        g.sig = y2;
        g.tot = (TH1D *) y2->Clone("yt_tot_" + nm);
        g.tot->SetDirectory(nullptr);
        if (y1) g.tot->Add(y1);
        if (y0) g.tot->Add(y0);
        return g;
    };
    Yields nom = load("", "nom");
    Yields drp = load("_trkdrop030", "drop");
    if (!nom.tot || !drp.tot) {
        std::cerr << "ERROR: missing " << (!nom.tot ? "nominal" : "3% track-drop") << " Pythia yields for tag '"
                  << tag << "'" << std::endl;
        return;
    }
    const double scale = (nom.tot->Integral() > 0.) ? n_jets / nom.tot->Integral() : 1.;
    for (TH1D *h : {nom.tot, nom.sig, drp.tot, drp.sig}) h->Scale(scale);   // ONE factor, see above

    // ---- printout ----
    std::cout << "\n================ raw yields, nominal vs 3% track drop: " << observable << ", "
              << sample << ", tag '" << tag << "', " << ptLabel(ptbin) << " ================" << std::endl;
    std::cout << "  both MC scaled by the nominal's factor x" << std::setprecision(4) << std::fixed << scale
              << "; integral drop / nominal: all " << drp.tot->Integral() / nom.tot->Integral()
              << ", 2B " << (nom.sig->Integral() > 0. ? drp.sig->Integral() / nom.sig->Integral() : 0.) << std::endl;
    std::cout << "    " << observable << " bin      |    data  |  nominal |  3% drop | drop/nom | 2B drop/nom | nom/data"
              << std::endl;
    for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
        const double d = y_data->GetBinContent(i), n = nom.tot->GetBinContent(i), p = drp.tot->GetBinContent(i);
        const double ns = nom.sig->GetBinContent(i), ps = drp.sig->GetBinContent(i);
        std::cout << "    " << std::setprecision(3)
                  << std::setw(5) << y_data->GetXaxis()->GetBinLowEdge(i)
                  << " - " << std::setw(5) << y_data->GetXaxis()->GetBinUpEdge(i)
                  << " | " << std::setw(8) << std::setprecision(0) << d
                  << " | " << std::setw(8) << n
                  << " | " << std::setw(8) << p
                  << " | " << std::setw(8) << std::setprecision(4) << (n > 0 ? p / n : 0.)
                  << " | " << std::setw(11) << (ns > 0 ? ps / ns : 0.)
                  << " | " << std::setw(8) << (d > 0 ? n / d : 0.) << std::endl;
    }

    // ---- plot ----
    const Float_t font_scale  = 1200. / 800.;
    const Style_t font_code   = 43;
    const Float_t label_size  = 15. * font_scale;
    const Float_t title_size  = 15. * font_scale;
    const Float_t legend_size = 13. * font_scale;
    const Color_t col_nom = TFColor::cmpA();
    const Color_t col_drp = TFColor::cmpB();

    TCanvas c("c_raw_trk", "", 800, 850);
    TPad *top = new TPad("rt_top", "", 0., 0.32, 1., 1.);
    TPad *bot = new TPad("rt_bot", "", 0., 0.,   1., 0.32);
    for (TPad *q : {top, bot}) { q->SetTicks(1, 1); q->SetFillColor(0); }
    top->SetMargin(0.15, 0.05, 0.0,  0.09);
    bot->SetMargin(0.15, 0.05, 0.26, 0.0);
    c.cd(); top->Draw(); bot->Draw();

    auto styleLine = [](TH1D *h, Color_t col, Style_t ls) {
        h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(ls);
        h->SetMarkerSize(0); h->SetFillStyle(0);
    };
    styleLine(nom.tot, col_nom, kSolid);  styleLine(nom.sig, col_nom, kDashed);
    styleLine(drp.tot, col_drp, kSolid);  styleLine(drp.sig, col_drp, kDashed);
    styleData(y_data);

    auto styleAxes = [&](TH1D *f) {
        for (TAxis *ax : { f->GetXaxis(), f->GetYaxis() }) {
            ax->CenterTitle(true);
            ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
        }
        if (observable != "dr") f->GetXaxis()->SetRange(1, f->GetNbinsX());
    };

    top->cd();
    TH1D *frame = (TH1D *) y_data->Clone("rt_frame");
    frame->SetDirectory(nullptr); frame->Reset(); frame->SetTitle("");
    frame->SetMinimum(0.001);
    frame->SetMaximum(1.5 * std::max({y_data->GetMaximum(), nom.tot->GetMaximum(), drp.tot->GetMaximum()}));
    frame->GetYaxis()->SetTitle("2b jets per bin");
    styleAxes(frame);
    frame->GetYaxis()->SetTitleOffset(1.9);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("AXIS");
    nom.sig->Draw("hist same");
    drp.sig->Draw("hist same");
    nom.tot->Draw("hist same");
    drp.tot->Draw("hist same");
    y_data->Draw("PE X0 same");

    TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.037); note.SetTextAlign(13);
    note.DrawLatex(0.19, 0.86, ptLabel(ptbin));
    note.DrawLatex(0.19, 0.80, "Raw yields, reco level, Pythia 8");
    note.DrawLatex(0.19, 0.74, "MC: one scale (nominal to data)");

    TLegend leg(0.52, 0.56, 0.94, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetMargin(0.22);
    leg.SetTextFont(font_code); leg.SetTextSize(legend_size);
    leg.AddEntry(y_data,  "Data", "pe1");
    leg.AddEntry(nom.tot, "Nominal, all flavours", "l");
    leg.AddEntry(nom.sig, "Nominal, 2B only", "l");
    leg.AddEntry(drp.tot, "3% tracks dropped, all flavours", "l");
    leg.AddEntry(drp.sig, "3% tracks dropped, 2B only", "l");
    leg.Draw();
    drawCMSHeader(top);
    top->RedrawAxis();

    bot->cd();
    TH1D *rframe = (TH1D *) y_data->Clone("rt_rframe");
    rframe->SetDirectory(nullptr); rframe->Reset(); rframe->SetTitle("");
    TH1D *r_tot = (TH1D *) drp.tot->Clone("rt_r_tot"); r_tot->SetDirectory(nullptr);
    TH1D *r_sig = (TH1D *) drp.sig->Clone("rt_r_sig"); r_sig->SetDirectory(nullptr);
    double rdev = 0.;
    for (int i = 1; i <= y_data->GetNbinsX(); ++i) {
        const double n = nom.tot->GetBinContent(i), ns = nom.sig->GetBinContent(i);
        // no ratio errors: the two MC samples are the same events, almost fully correlated
        r_tot->SetBinContent(i, n  > 0 ? drp.tot->GetBinContent(i) / n  : 0.); r_tot->SetBinError(i, 0.);
        r_sig->SetBinContent(i, ns > 0 ? drp.sig->GetBinContent(i) / ns : 0.); r_sig->SetBinError(i, 0.);
        if (n  > 0) rdev = std::max(rdev, std::fabs(r_tot->GetBinContent(i) - 1.));
        if (ns > 0) rdev = std::max(rdev, std::fabs(r_sig->GetBinContent(i) - 1.));
    }
    const double rhalf = std::max(0.05, 1.3 * rdev);
    rframe->SetMinimum(1. - rhalf); rframe->SetMaximum(1. + rhalf);
    rframe->GetYaxis()->SetTitle("3% drop / nominal");
    rframe->GetXaxis()->SetTitle(obsd.axis);
    styleAxes(rframe);
    rframe->GetYaxis()->SetTitleOffset(1.9);
    rframe->GetXaxis()->SetTitleOffset(1.1);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->Draw("AXIS");
    TLine one(obsd.bins[0], 1., obsd.bins[obsd.nbins], 1.);
    one.SetLineStyle(2); one.SetLineColor(kBlack); one.Draw();
    r_tot->Draw("hist same");
    r_sig->Draw("hist same");
    bot->RedrawAxis();

    const TString outdir = TString::Format("%s/results/%s_raw_yields_trackeff_%s%s_%s",
                                           kBase, observable.Data(), sample.Data(),
                                           (svsfTag() + trackBdtSfTag()).Data(), tag.Data());
    gSystem->mkdir(outdir, true);
    const TString stem = outdir + "/" + observable + "_raw_yields_nominal_vs_trkdrop" + ptSuffix(ptbin);
    c.SaveAs(stem + ".pdf");
    c.SaveAs(stem + ".png");
}

void plot_raw_yields_trackeff(const TString &sample     = "both",
                              const TString &tag        = "upartv2rhojecvetobenrptw",
                              const TString &observable = "dr")
{
    for (int ptbin = 1; ptbin <= jtpt_bins; ++ptbin)
        raw_yields_trackeff_one(sample, tag, ptbin, observable);
}


// ============================================================================
// DATA AFTER THE TEMPLATE FIT: Pythia vs Herwig templates (added 2026-09-29)
// ============================================================================
// The data times the signal fraction each fit extracts, bin by bin -- exactly the step
// apply_unfolding_2d.C does ("Multiplying by signal fraction"), before the UParT SF and
// before any unfolding. Read straight from the two template-fit output files, so the ONLY
// thing that differs between the two curves is which generator's templates were fitted:
//   data   = h3D_data of the fit file, projected over m_2B, one pT bin
//   f_sig  = h_sig_fraction_fit (observable x pT)
// Both files hold the same data (checked below). For dr the fit, and so the data here, is
// EEC-weighted; for B it is the yield fit (_noeecw) -- as in the measurement chain.
//
// Top pad : data before the fit (black), after the Pythia fit (orange), after the Herwig
//           fit (purple). Lower: Herwig / Pythia, which is the ratio of the signal fractions.
// Errors on the after-fit curves: data stat (+) the fit's signal-fraction error.
//
//   root -l -b -q -e '.L template_fit.cpp+' -e 'plot_data_after_fit_generators("dr")'
//   root -l -b -q -e '.L template_fit.cpp+' -e 'plot_data_after_fit_generators("B")'
static TString tfFitFile(const TString &generator, const TString &observable, bool eecOff)
{
    // Same naming as templateFitFile() in apply_unfolding_2d.C: sample "both", nominal fit.
    const TString eecw_tag = eecOff ? "_noeecw" : "";
    const TString obs_tag  = (observable == "dr") ? "" : ("_" + observable);
    return TString(kBase) + "/results/TemplateFit_Run3/TemplateFits_both_" + generator
         + svsfTag() + trackBdtSfTag() + eecw_tag + obs_tag + eraTag()
         + "/nominal_Run3_TemplateFits_histos_3d_80_inf.root";
}

static void data_after_fit_generators_one(const TString &observable, int ptbin)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const ObsDef obsd = obsByName(observable);
    if (obsd.nbins == 0) {
        std::cerr << "ERROR: unknown observable '" << observable << "' (use dr | B)" << std::endl;
        return;
    }
    const bool eecOff = nominalEecWeightOff(observable);   // B: yield fit; dr: EEC-weighted

    struct FitOut { TH1D *data = nullptr; TH1D *after = nullptr; TH1D *fsig = nullptr; };
    auto load = [&](const TString &gen) -> FitOut {
        FitOut o;
        const TString fn = tfFitFile(gen, observable, eecOff);
        TFile *f = TFile::Open(fn);
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open " << fn << " -- run template_fit.cpp for " << gen
                      << " first" << std::endl;
            return o;
        }
        TH3D *d3 = dynamic_cast<TH3D *>(f->Get("h3D_data"));
        TH2D *sf = dynamic_cast<TH2D *>(f->Get("h_sig_fraction_fit"));
        if (!d3 || !sf) {
            std::cerr << "ERROR: h3D_data / h_sig_fraction_fit missing in " << fn << std::endl;
            f->Close(); return o;
        }
        o.data = projX(d3, "af_data_" + gen, ptbin);
        o.fsig = sf->ProjectionX("af_fsig_" + gen, ptbin, ptbin);
        o.fsig->SetDirectory(nullptr);
        if (o.fsig->GetNbinsX() != o.data->GetNbinsX()) {
            std::cerr << "ERROR: signal fraction has " << o.fsig->GetNbinsX() << " bins, data "
                      << o.data->GetNbinsX() << " in " << fn << std::endl;
            f->Close(); o = FitOut(); return o;
        }
        o.after = (TH1D *) o.data->Clone("af_after_" + gen);
        o.after->SetDirectory(nullptr);
        o.after->Multiply(o.fsig);   // errors: data stat (+) fit error, in quadrature (relative)
        f->Close();
        return o;
    };
    FitOut py = load("pythia");
    FitOut hw = load("herwig");
    if (!py.after || !hw.after) return;

    // The two fits must have been run on the same data, or the comparison is meaningless.
    for (int i = 1; i <= py.data->GetNbinsX(); ++i) {
        const double a = py.data->GetBinContent(i), b = hw.data->GetBinContent(i);
        if (std::fabs(a - b) > 1e-6 * std::max(std::fabs(a), 1.)) {
            std::cerr << "ERROR: the Pythia and Herwig fits were run on DIFFERENT data (bin " << i
                      << ": " << a << " vs " << b << ") -- not comparing them." << std::endl;
            return;
        }
    }

    // ---- printout ----
    std::cout << "\n================ data after the template fit, Pythia vs Herwig: "
              << observable << ", " << ptLabel(ptbin) << (eecOff ? ", yields" : ", EEC-weighted")
              << " ================" << std::endl;
    std::cout << "    " << observable << " bin      |  f_sig Py      |  f_sig Hw      |  Hw/Py" << std::endl;
    std::cout << "    -------------+----------------+----------------+--------" << std::endl;
    for (int i = 1; i <= py.fsig->GetNbinsX(); ++i) {
        const double p = py.fsig->GetBinContent(i), h = hw.fsig->GetBinContent(i);
        std::cout << "    " << std::fixed << std::setprecision(3)
                  << std::setw(5) << py.data->GetXaxis()->GetBinLowEdge(i)
                  << " - " << std::setw(5) << py.data->GetXaxis()->GetBinUpEdge(i)
                  << " | " << std::setw(6) << p << " +- " << std::setw(5) << py.fsig->GetBinError(i)
                  << " | " << std::setw(6) << h << " +- " << std::setw(5) << hw.fsig->GetBinError(i)
                  << " | " << std::setw(6) << (p > 0 ? h / p : 0.) << std::endl;
    }
    std::cout << "    integral after fit: Pythia " << std::setprecision(4) << py.after->Integral()
              << ", Herwig " << hw.after->Integral() << "  (Hw/Py = "
              << (py.after->Integral() > 0 ? hw.after->Integral() / py.after->Integral() : 0.)
              << ")" << std::endl;

    // ---- plot ----
    const Float_t font_scale  = 1200. / 800.;
    const Style_t font_code   = 43;
    const Float_t label_size  = 15. * font_scale;
    const Float_t title_size  = 15. * font_scale;
    const Float_t legend_size = 13. * font_scale;
    const Color_t col_py = TFColor::cmpA();
    const Color_t col_hw = TFColor::cmpB();

    TCanvas c("c_af_gen", "", 800, 850);
    TPad *top = new TPad("af_top", "", 0., 0.32, 1., 1.);
    TPad *bot = new TPad("af_bot", "", 0., 0.,   1., 0.32);
    for (TPad *q : {top, bot}) { q->SetTicks(1, 1); q->SetFillColor(0); }
    top->SetMargin(0.15, 0.05, 0.0,  0.09);
    bot->SetMargin(0.15, 0.05, 0.26, 0.0);
    c.cd(); top->Draw(); bot->Draw();

    auto styleAfter = [](TH1D *h, Color_t col) {
        h->SetLineColor(col); h->SetLineWidth(2); h->SetMarkerColor(col);
        h->SetMarkerSize(0); h->SetFillStyle(0);
    };
    styleAfter(py.after, col_py);
    styleAfter(hw.after, col_hw);
    styleData(py.data);
    py.data->SetMarkerStyle(24);   // open: the data BEFORE the fit, not a result

    auto styleAxes = [&](TH1D *f) {
        for (TAxis *ax : { f->GetXaxis(), f->GetYaxis() }) {
            ax->CenterTitle(true);
            ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
        }
        if (observable != "dr") f->GetXaxis()->SetRange(1, f->GetNbinsX());
    };

    top->cd();
    TH1D *frame = (TH1D *) py.data->Clone("af_frame");
    frame->SetDirectory(nullptr); frame->Reset(); frame->SetTitle("");
    frame->SetMinimum(0.001);
    frame->SetMaximum(1.5 * py.data->GetMaximum());
    frame->GetYaxis()->SetTitle(eecOff ? "2b jets per bin" : "EEC-weighted data per bin");
    styleAxes(frame);
    frame->GetYaxis()->SetTitleOffset(1.9);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("AXIS");
    py.after->Draw("hist same");  py.after->Draw("E1 X0 same");
    hw.after->Draw("hist same");  hw.after->Draw("E1 X0 same");
    py.data->Draw("PE X0 same");

    TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.037); note.SetTextAlign(13);
    note.DrawLatex(0.19, 0.86, ptLabel(ptbin));
    note.DrawLatex(0.19, 0.80, "Reco level, no unfolding");
    note.DrawLatex(0.19, 0.74, eecOff ? "Yields" : "EEC weighted");

    TLegend leg(0.50, 0.66, 0.94, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetMargin(0.22);
    leg.SetTextFont(font_code); leg.SetTextSize(legend_size);
    leg.AddEntry(py.data, "Data, before fit", "pe1");
    leg.AddEntry(py.after, "Data #times f_{sig}, Pythia 8 fit", "l");
    leg.AddEntry(hw.after, "Data #times f_{sig}, Herwig fit", "l");
    leg.Draw();
    drawCMSHeader(top);
    top->RedrawAxis();

    bot->cd();
    TH1D *rframe = (TH1D *) py.data->Clone("af_rframe");
    rframe->SetDirectory(nullptr); rframe->Reset(); rframe->SetTitle("");
    // Herwig / Pythia of the after-fit data = f_sig(Hw) / f_sig(Py); the data cancels, so the
    // error is the two fits' errors only.
    TH1D *r = (TH1D *) hw.fsig->Clone("af_ratio"); r->SetDirectory(nullptr);
    r->Divide(py.fsig);
    double rdev = 0.;
    for (int i = 1; i <= r->GetNbinsX(); ++i)
        rdev = std::max(rdev, std::fabs(r->GetBinContent(i) - 1.) + r->GetBinError(i));
    const double rhalf = std::max(0.15, 1.2 * rdev);
    rframe->SetMinimum(1. - rhalf); rframe->SetMaximum(1. + rhalf);
    rframe->GetYaxis()->SetTitle("Herwig / Pythia");
    rframe->GetXaxis()->SetTitle(obsd.axis);
    styleAxes(rframe);
    rframe->GetYaxis()->SetTitleOffset(1.9);
    rframe->GetXaxis()->SetTitleOffset(1.1);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->Draw("AXIS");
    TLine one(obsd.bins[0], 1., obsd.bins[obsd.nbins], 1.);
    one.SetLineStyle(2); one.SetLineColor(kBlack); one.Draw();
    r->SetLineColor(col_hw); r->SetMarkerColor(col_hw); r->SetMarkerStyle(20); r->SetMarkerSize(1.0);
    r->SetLineWidth(2);
    r->Draw("PE1 X0 same");
    bot->RedrawAxis();

    const TString outdir = TString::Format("%s/results/%s_data_after_fit_generators_both%s%s",
                                           kBase, observable.Data(),
                                           (svsfTag() + trackBdtSfTag()).Data(),
                                           eraTag().Data());
    gSystem->mkdir(outdir, true);
    const TString stem = outdir + "/" + observable + "_data_after_fit_pythia_vs_herwig" + ptSuffix(ptbin);
    c.SaveAs(stem + ".pdf");
    c.SaveAs(stem + ".png");
}

void plot_data_after_fit_generators(const TString &observable = "dr")
{
    for (int ptbin = 1; ptbin <= jtpt_bins; ++ptbin)
        data_after_fit_generators_one(observable, ptbin);
}


// ============================================================================
// plot_jet_pt_spectrum() (2026-10-29, Zoe: "plot the jet pt spectrum of reco pt")
// Reco jet pT of the selected jets (2 SVs, b-tagged, |eta| < 2, veto map) from the _ptfine histograms
// of the A production (no EEC weight: one entry per jet), summed over m2B and A. Data against QCD Pythia
// split by flavour (2B / 1B / 0B), MC scaled to the data integral; ratio data/MC below. Lines at the
// measured-bin edges 100 and 200 GeV. Read from the open-ended chain (JET_PT_MAX unset), so the part
// the 200 GeV cut removes is visible. Needs the "<era>ptw" production (5 GeV bins, 80-400 GeV).
//   root -l -b -q -e '.L template_fit.cpp+' -e 'plot_jet_pt_spectrum()'
// ============================================================================
void plot_jet_pt_spectrum()
{
    const TString P = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
    TString era = eraTag(); era.Remove(0, 1);   // "upartv2rhojecvetobenr" (+ ptmax.. if JET_PT_MAX is set)
    const TString fdata = P + "/HardProbes/agg_template_chunks/Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_noeecw_" + era + "ptw_B.root";
    const TString fqcd  = P + "/QCD/agg_ntuple_chunks/Run3_btagWP0712_template_for_fit_histos_3D_qcd_fMCGEN" + svsfTag() + trackBdtSfTag() + "_noeecw_" + era + "ptw_B.root";
    auto proj = [](const TString &f, const char *h, const char *name) -> TH1D * {
        TFile *in = TFile::Open(f, "read");
        if (!in || in->IsZombie()) { std::cerr << "ERROR: cannot open " << f << std::endl; return nullptr; }
        TH3D *h3 = dynamic_cast<TH3D *>(in->Get(h));
        if (!h3) { std::cerr << "ERROR: no " << h << " in " << f << std::endl; in->Close(); return nullptr; }
        TH1D *p = h3->ProjectionZ(name, 1, h3->GetNbinsX(), 1, h3->GetNbinsY(), "e");
        p->SetDirectory(nullptr); in->Close();
        return p;
    };
    TH1D *hd  = proj(fdata, "h3D_data_ptfine_B", "jpt_data");
    TH1D *hbb = proj(fqcd,  "h3D_bb_ptfine_B",   "jpt_2b");
    TH1D *hb  = proj(fqcd,  "h3D_b_ptfine_B",    "jpt_1b");
    TH1D *h0  = proj(fqcd,  "h3D_0b_ptfine_B",   "jpt_0b");
    if (!hd || !hbb || !hb || !h0) return;
    TH1D *hmc = (TH1D *) hbb->Clone("jpt_mc"); hmc->Add(hb); hmc->Add(h0);
    const double s = hd->Integral() / hmc->Integral();
    for (TH1D *h : {hbb, hb, h0, hmc}) h->Scale(s);
    for (TH1D *h : {hd, hbb, hb, h0, hmc}) h->Scale(1., "width");

    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    TCanvas c("c_jetpt_spectrum", "", 1000, 1000);
    TPad *p1 = new TPad("pjs1", "", 0, 0.3, 1, 1), *p2 = new TPad("pjs2", "", 0, 0, 1, 0.3);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.14); p1->SetTopMargin(0.08); p1->SetLogy(); p1->SetTicks(1, 1);
    p2->SetTopMargin(0.03); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.14); p2->SetTicks(1, 1);
    p1->Draw(); p2->Draw(); p1->cd();
    THStack *st = new THStack("st_jpt", "");
    for (auto pr : std::vector<std::pair<TH1D *, Color_t>>{{h0, TFColor::c0B()}, {hb, TFColor::c1B()}, {hbb, TFColor::c2B()}}) {
        pr.first->SetFillColor(pr.second); pr.first->SetLineColor(pr.second); st->Add(pr.first);
    }
    hd->SetMarkerStyle(20); hd->SetMarkerColor(TFColor::data()); hd->SetLineColor(TFColor::data());
    hd->SetMaximum(5. * hd->GetMaximum()); hd->SetMinimum(0.5 * hd->GetMinimum(0.));
    hd->GetYaxis()->SetTitle("dN_{jet}/dp_{T}^{jet} [GeV^{-1}]"); hd->GetXaxis()->SetLabelSize(0);
    hd->GetYaxis()->SetTitleSize(0.05); hd->GetYaxis()->SetTitleOffset(1.25);
    hd->Draw("PE"); st->Draw("hist same"); hd->Draw("PE same");
    const double ylo = hd->GetMinimum(), yhi = hd->GetMaximum();
    for (double x : {100., 200.}) { TLine *l = new TLine(x, ylo, x, yhi); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
    TLegend *lg = new TLegend(0.58, 0.62, 0.92, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
    lg->AddEntry(hd, "Data", "pe"); lg->AddEntry(hbb, "2B (QCD Pythia)", "f"); lg->AddEntry(hb, "1B", "f"); lg->AddEntry(h0, "0B", "f");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.035);
    tx.SetTextSize(0.032);   // under the legend, clear of the stack
    tx.DrawLatex(0.50, 0.56, "2 SVs, b-tagged, |#eta_{jet}| < 2, one entry per jet");
    tx.DrawLatex(0.50, 0.51, "MC scaled to data; lines at 100 and 200 GeV");
    drawCMSHeader(p1);
    p2->cd();
    TH1D *r = (TH1D *) hd->Clone("jpt_ratio"); r->Divide(hmc);
    r->SetMinimum(0.4); r->SetMaximum(1.6);
    r->GetYaxis()->SetTitle("Data / MC"); r->GetYaxis()->CenterTitle(true); r->GetYaxis()->SetNdivisions(505);
    r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleOffset(0.55);
    r->GetXaxis()->SetTitle("reco p_{T}^{jet} [GeV]"); r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetLabelSize(0.10);
    r->GetXaxis()->SetLabelSize(0.10);
    r->Draw("PE");
    TLine *one = new TLine(r->GetXaxis()->GetXmin(), 1., r->GetXaxis()->GetXmax(), 1.); one->SetLineStyle(2); one->Draw();
    for (double x : {100., 200.}) { TLine *l = new TLine(x, 0.4, x, 1.6); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
    const TString out = P + "/results/jet_pt_spectrum_" + era;
    gSystem->mkdir(out, kTRUE);
    c.Print(out + "/jet_pt_reco_spectrum.png"); c.Print(out + "/jet_pt_reco_spectrum.pdf");
    const int b100 = hd->FindBin(100.01), b200 = hd->FindBin(200.01);
    printf("data jets: 80-100 %.0f, 100-200 %.0f, 200-400 %.0f (fraction above 200 of >100: %.3f)\n",
           hd->Integral(1, b100 - 1, "width"), hd->Integral(b100, b200 - 1, "width"), hd->Integral(b200, hd->GetNbinsX(), "width"),
           hd->Integral(b200, hd->GetNbinsX(), "width") / hd->Integral(b100, hd->GetNbinsX(), "width"));
    std::cout << "Wrote " << out << "/jet_pt_reco_spectrum.{png,pdf}" << std::endl;
}
