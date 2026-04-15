/*
 * ApplyRegularization.C
 * 
 * Description:
 *   - This macro imports the 'data_release.root' file and applies the regularization matrix to a truth predictions.
 * 
 * Usage:
 *   root -l ApplyRegularization.C
 *
 * Author: J.P. Detje
 * Date: 14 Aug 2025
 */

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include <string> // added

void ApplyRegularization()
{
    // Open the release file next to this macro
    TFile* file = TFile::Open("data_release.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open data_release.root" << std::endl;
        return;
    }

    // Fetch TMatrixD objects
    auto* mAc           = dynamic_cast<TMatrixD*>(file->Get("A_C_diff"));
    auto* vPrediction   = dynamic_cast<TMatrixD*>(file->Get("uBTune_diff"));
    auto* vData         = dynamic_cast<TMatrixD*>(file->Get("unfolded_data_diff"));
    auto* mCov          = dynamic_cast<TMatrixD*>(file->Get("cov_unfolded_diff"));
    auto* vWidth        = dynamic_cast<TMatrixD*>(file->Get("bin_widths_diff"));

    const bool applyBinWidth = true; // Normalise bins by width to match plots in paper

    if (!mAc || !vPrediction || !vData || !mCov || !vWidth) {
        std::cerr << "Error: Missing inputs (A_C_diff, uBTune_diff, unfolded_data_diff, cov_unfolded_diff, bin_widths_diff)." << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Prediction vector is concatenated cos(theta_mu), p_mu, cos(theta_pi), p_pi and theta_{mu pi} bins; Substitute with your own prediction.
    TMatrixD prediction    = *vPrediction;  // n x 1 prediction vector
    TMatrixD ac            = *mAc;          // n x n regularization matrix
    TMatrixD data          = *vData;        // n x 1 unfolded data vector
    TMatrixD cov           = *mCov;         // n x n unfolded covariance matrix
    TMatrixD width         = *vWidth;       // n x 1 bin widths vector
    file->Close();
    delete file;

    const int n = data.GetNrows();
    if (ac.GetNrows() != ac.GetNcols() || ac.GetNcols() != n ||
        prediction.GetNrows() != n || prediction.GetNcols() != 1 ||
        width.GetNrows() != n || width.GetNcols() != 1) {
        std::cerr << "Error: Dimension mismatch." << std::endl;
        return;
    }

    // Apply A_C smearing to the prediction: v_smear = A_C * v
    TMatrixD predictionSmeared(ac, TMatrixD::kMult, prediction);

    const std::string unitSuffix = applyBinWidth ? "(/1 or /rad or /GeV)" : ""; // When normalising by bin width, each bin has units of that variable
    const std::string titlePred  = std::string("A_{C} #times prediction;Bin index;10^{-38} cm^{2}/Ar") + unitSuffix;
    const std::string titleData  = std::string("A_{C} #times prediction vs unfolded data;Bin index;10^{-38} cm^{2}/Ar") + unitSuffix;
    TH1D* hPredictionSmeared = new TH1D("uBTune_diff_smeared", titlePred.c_str(), n, 0.0, n);
    TH1D* hData              = new TH1D("unfolded_data_diff_hist", titleData.c_str(), n, 0.0, n);
    hPredictionSmeared->SetDirectory(nullptr);
    hData->SetDirectory(nullptr);

    for (int i = 0; i < n; ++i) {
        const double w = applyBinWidth ? width(i, 0) : 1.0;
        const double pred = predictionSmeared(i, 0);
        const double val  = data(i, 0);
        const double err  = std::sqrt(cov(i, i));

        hPredictionSmeared->SetBinContent(i + 1, pred / w);
        hPredictionSmeared->SetBinError(i + 1, 0.0);

        hData->SetBinContent(i + 1, val / w);
        hData->SetBinError(i + 1, err / w);
    }

    TCanvas* c = new TCanvas("c_apply_ac", "A_C x uBTune vs unfolded", 900, 650);
    gStyle->SetOptStat(0);
    hData->Draw("E1");
    hPredictionSmeared->SetLineColor(kRed);
    hPredictionSmeared->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.10, 0.80, 0.45, 0.90);
    leg->AddEntry(hData, "Unfolded data", "lep");
    leg->AddEntry(hPredictionSmeared, "A_{C} #times prediction", "l");
    leg->Draw();
}