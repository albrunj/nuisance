
Data release for "Measurement of single charged pion production in charged-current νμ-Ar interactions with the MicroBooNE detector"

flux.root
    Contains:

    - Simulated muon neutrino energy distribution (GeV) at MicroBooNE:
      hEnumu_cv

    - Simulated muon antineutrino energy distribution (GeV) at MicroBooNE:
      hEnumubar_cv

    Create numu and numubar predictions according to these flux distributions, then add the histograms together:
        (fluxIntNuMuBar/(fluxIntNuMu + fluxIntNuMuBar))*numubarPrediction + fluxIntNuMu/(fluxIntNuMu + fluxIntNuMuBar)*numuPrediction,
    where fluxIntNuMu and fluxIntNuMuBar are the integrated distributions (e.g. via hEnumu_cv->Integral()).

data_release.root
    Contains:

    - Histogram of the unfolded data for the differential cross section with full uncertainty:
      unfolded_data_total (Units of 10^{-38} cm^{2}/Ar)

    - The unfolded data for the differential cross section bins (without uncertainties):
      unfolded_data_diff (without bin-width normalisation; units of 10^{-38} cm^{2}/Ar)

    - The covariance matrix with uncertainties on unfolded_data_diff:
      cov_unfolded_diff (Units of (10^{-38} cm^{2}/Ar)^2)

    - The bin widths. Element-wise dividing unfolded_data_diff by this is needed to reproduce the plots in the paper:
      bin_widths_diff (Units of 1, rad and GeV)

    - The 2D regularization matrix that needs to be applied to predictions to compare them with data (not needed for total cross-section):
      A_C_diff

    - As an example for a truth prediction, the GENIE v3 MicroBooNE tune prediction for the differential cross section bins is included:
      uBTune_diff (Units of 10^{-38} cm^{2}/Ar)

ApplyRegularization.C
    A ROOT macro that takes the square-root of the diagonal elements from the covariance matrix, sets them as uncertainties for the unfolded data. It then applies the regularization matrix to a truth prediction (uBTune_diff by default). Both histograms are scaled by their bin widths and plotted.
    The resulting histogram shows the concatenated differential cross section bins:
    Bins 0-10:  cos(theta_mu)
    Bins 11-15: p_mu [Units of GeV]
    Bins 16-22: cos(theta_pi)
    Bins 23-26: p_pi [Units of GeV]
    Bins 27-33: theta_{mu pi} [Units of rad]