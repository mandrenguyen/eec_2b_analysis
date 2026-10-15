

#include "tTree.h"
#include "binning_histos_small.h"
#include "result_paths.h"   // obsProdTag(): which MC production carries which observable
#include "Help_Functions.h"
#include "Draw_EEC.h"
//#include "../CMSStyle.C" // CMS style impored from Matthew


void do_template_fit_combined(const TString &HighEGdata_name, const TString &LowEGdata_name, TString &templates, TString &templates_bjet, TString pT_selection, TString folder, TString &fout_name, bool& alsoLowEG, bool& also_bjet,  Variation ivar = NOMINAL){
    /*
        // ----- WORK in PROGRESS ---- 14/04/2026
        // 2B now covers 2B + more-than-2B in one histogram. 0B comes from the qcd sample
        // only: the bjet sample is filtered to b jets, so its 0B templates are not physical.
        // Argument add for possible varations on the fit
    */

    // Test CMS style 
    // setCMSStyle();

    // All plots go into one flat directory (sDirname_www); root files stay in sDirname.
    // Only the nominal variation is plotted; the others are still fitted because
    // draw_variation_uncertainity() needs their fractions for the 0B systematic.
    const bool save_png = (ivar == NOMINAL);
    TString sDir_canvas = sDirname_www;
        if (save_png) gSystem->mkdir(sDir_canvas, kTRUE);

/* ---- disabled (kept for reference): per-variation subfolders and the separate _www tree ----
    TString sDir_canvas = Form("%s/%s", sDirname.Data(), varNames[ivar].Data());
        gSystem->mkdir(sDir_canvas, kTRUE); // to be clean

    // Make subdirectory for webversion only
    TString sDir_canvas_www = Form("%s/%s", sDirname_www.Data(), varNames[ivar].Data());
        gSystem->mkdir(sDir_canvas_www, kTRUE);
*/ // ---- end disabled block ----

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
                    if ( std::abs(int_val - 1)  > 1e-06) {cout << "Effective BKG PDF is not normalized to 1!"<< endl; return;}
                    if ( std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1)  > 1e-06) {cout << "Signal PDF is not normalized to 1!"<< endl; return;}

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
            double err_true_frac = ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B 
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

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
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
                    h->Write();
    }

    // //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
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
    TString namehData = "h3D_data"; // h3D_data with eec on the weight 
    TString namehMC = "";

    // -- Input Histogram name: signal and bkg templates
    TString nameh2B = "h3D_bb";
    TString nameh1B = "h3D_b";

    TString nameh0B = "h3D_0b";
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
        bins_dr = dr_binsVectorSize   - 1;
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
    yBins = dr_binsVector;     
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
                    if ( std::abs(int_val - 1)  > 1e-06) {cout << "Effective BKG PDF is not normalized to 1!"<< endl; return;}
                    if ( std::abs( h_sig->Integral(1,  h_sig_bins, "width") - 1)  > 1e-06) {cout << "Signal PDF is not normalized to 1!"<< endl; return;}

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
            double err_true_frac = ( (int0 + int1) * True_sig_err + int2 * True_bkg_err )/TMath::Power(int0 + int1+ int2, 2);// 0B + 1B and 2B 
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

            TH1D *h_bkg_fit_1b = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_1b_%d_%d", ibin_dr, ibin_pt)); // 1B only
                h_bkg_fit_1b->Scale(integral_inputdata * p1);
                h_bkg_fit_1b->SetTitle("1B");
                styleTemplate(h_bkg_fit_1b, TFColor::c1B());   // 1B -> blue

            TH1D *h_bkg_fit_nob = (TH1D*) h_bkg->Clone(Form("h_bkg_fit_nob_%d_%d", ibin_dr, ibin_pt)); // 0B only
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
    sDirname     = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/"
                   "TemplateFits_" + SAMPLE + "_" + GENERATOR + trk_tag + svsfTag()
                   + trackBdtSfTag() + eecw_tag + obs_tag + eraTag();
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

std::vector<TString> mcFiles(const TString &sample, const TString &generator,
                             const TString &tag, bool rmatrix)
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
                (svsfTag() + trackBdtSfTag()).Data(), tag.Data());
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
        if (!tag.Contains(stem)) {
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
        if (!tag.Contains(stem)) {
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
