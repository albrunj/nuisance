// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
*    This file is part of NUISANCE.
*
*    NUISANCE is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    NUISANCE is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include <fstream>

#include "InteractionModes.h"
#include "MicroBooNEBlockHandler.h"
#include "MicroBooNEHelper.h"
#include "MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu.h"
#include "TMatrixD.h"

// Anonymous namespace
namespace {
  constexpr int MASS_NUMBER_40AR = 40;
  constexpr int PROTON = 2212;
  constexpr int MU_MINUS = 13;
  constexpr int PION_PLUS = 211;
  constexpr int PION_MINUS = -211;

  constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
  constexpr double NEUTRON_MASS = 0.93956541; // GeV
  constexpr double PROTON_MASS = 0.93827208; // GeV
  constexpr double MUON_MASS = 0.10565837; // GeV
  constexpr double PION_MASS = 0.13957039; // GeV
  constexpr double BINDING_ENERGY = 0.02478; // 40Ar, GeV

}

MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu
  ::MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu( nuiskey samplekey )
{
  fSettings = LoadSampleSettings( samplekey );
  std::string name = fSettings.GetS( "name" );

  // The main histograms use the bin number on the x-axis
  fSettings.SetXTitle( "bin number" );

  fSettings.SetYTitle("#sigma (10^{-38} cm^{2}/^{40}Ar)");

  // Sample overview ---------------------------------------------------
  std::string descrip = name + " sample.\n" \
                        "Target: Ar\n" \
                        "Flux: BNB FHC numu\n" \
                        "Signal: CC1Pipm\n"
                        "Contact: microboone_info@fnal.gov\n"
                        "Reference: Phys. Rev. D 113, 032007 (2026)\n"
                        "DOI: https://doi.org/10.1103/t2cw-cdx2\n";

  fSettings.SetDescription( descrip );
  fSettings.SetTitle( name );
  fSettings.SetAllowedTypes( "FULL", "FIX/FULL" );
  fSettings.SetEnuRange( 0.0, 6.8 );
  fSettings.DefineAllowedTargets( "Ar" );
  fSettings.DefineAllowedSpecies( "numu" );
  FinaliseSampleSettings();

  // Scale factor for the flux-averaged total cross section
  // (10^{-38} cm^2 / Ar) in each bin
  fScaleFactor = GetEventHistogram()->Integral( "width" )
    * MASS_NUMBER_40AR / fNEvents / TotalIntegratedFlux();

  // Get bin definitions
  this->LoadBinDefinitions();

  std::string file_name_AC( FitPar::GetDataBase()
    + "/MicroBooNE/BNB_NumuCC1Pipm_2025/mat_table_add_smear.txt" );

  // Load the additional smearing matrix
  TMatrixD A_C = MicroBooNEHelper::load_matrix( file_name_AC );
  fAddSmear = std::make_shared< TMatrixD >( A_C );

  // Load the measured data points
  std::string file_name_data( FitPar::GetDataBase()
    + "/MicroBooNE/BNB_NumuCC1Pipm_2025/vec_table_unfolded_signal.txt" );
  TMatrixD temp_data_mat = MicroBooNEHelper::load_matrix( file_name_data );
  fDataHist = MicroBooNEHelper::to_histogram( temp_data_mat );

  fDataHist->SetNameTitle( (fSettings.GetName() + "_data").c_str(),
    fSettings.GetFullTitles().c_str() );

  // Also retrieve the total covariance matrix for the measurement
  std::string cov_file_name( FitPar::GetDataBase()
    + "/MicroBooNE/BNB_NumuCC1Pipm_2025/mat_table_cov_total.txt" );
  auto temp_cov_matrix = MicroBooNEHelper::load_matrix( cov_file_name );
  fFullCovar = MicroBooNEHelper::to_symmetric_matrix( temp_cov_matrix );

  // Now invert the covariance matrix and store the result
  temp_cov_matrix.Invert();
  covar = MicroBooNEHelper::to_symmetric_matrix( temp_cov_matrix );

  fDecomp = StatUtils::GetDecomp( fFullCovar );
  if ( !fDecomp ) {
    NUIS_ABORT( "Failed to compute Cholesky decomposition of covariance matrix" );
  }
  
  /*
  TDecompChol chol( *fFullCovar );
  chol.Decompose();
  fDecomp = new TMatrixDSym( fFullCovar->GetNrows(),
    chol.GetU().GetMatrixArray(), "" );
  */

  // Push the diagonals of fFullCovar onto the data histogram
  StatUtils::SetDataErrorFromCov( fDataHist, fFullCovar, 1.0, false );

  // Setup fMCHist from data
  fMCHist = dynamic_cast< TH1D* >( fDataHist->Clone() );
  fMCHist->SetNameTitle( (fSettings.GetName() + "_MC").c_str(),
   fSettings.GetFullTitles().c_str() );
  fMCHist->Reset();

  fMCStat = dynamic_cast< TH1D* >( fMCHist->Clone() );
  fMCStat->Reset();

  // Since we're using a 1D histogram with bin number along the x-axis, it
  // doesn't make sense to subdivide the bins. Rather than doing that
  // automatically, I just copy the original data histogram binning here. Thus,
  // MCFine ends up using the same bins as regular MC.
  fMCFine = dynamic_cast< TH1D* >( fMCHist->Clone() );
  fMCFine->SetNameTitle( (fSettings.GetName() + "_MC_FINE").c_str(),
    fSettings.GetFullTitles().c_str() );
  fMCFine->Reset();

  // Set up the MC modes histogram
  fMCHist_Modes = new TrueModeStack( (fSettings.GetName() + "_MODES").c_str(),
    "True Channels", fMCHist );
  fMCHist_Modes->SetTitleX( fDataHist->GetXaxis()->GetTitle() );
  fMCHist_Modes->SetTitleY( fDataHist->GetYaxis()->GetTitle() );
  this->SetAutoProcessTH1( fMCHist_Modes, kCMD_Reset, kCMD_Norm, kCMD_Write );

  //this->FinaliseMeasurement();
}


bool MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::isSignal( FitEvent* event ) {
  // PDG codes of interest
  const int MUON = 13;
  const int ANTI_MUON = -13;
  const int MUON_NEUTRINO = 14;
  const int ANTI_MUON_NEUTRINO = -14;
  const int PROTON = 2212;
  const int PION_PLUS = 211;
  const int PION_MINUS = -211;

  // Require the event to be a numu CC inclusive interaction
  if ( !SignalDef::isCCINC(event, MUON_NEUTRINO, EnuMin, EnuMax) && !SignalDef::isCCINC(event, ANTI_MUON_NEUTRINO, EnuMin, EnuMax) ) return false;

  // Require exactly one charged pion in the final state
  int npi = event->NumFSParticle(PION_PLUS)
        + event->NumFSParticle(PION_MINUS);

  if ( npi != 1 ) return false;

  // Reject events with neutral pions of any momenta
  if (event->NumFSParticle(111) != 0) return false;
  if (event->NumFSParticle(-111) != 0) return false;

  //reject events with kaons of any momenta
  if (event->NumFSParticle(321) != 0) return false;
  if (event->NumFSParticle(-321) != 0) return false;

  // Impose kinematic limits in the signal definition
  double p_mu = event->GetHMFSParticle( MUON )->fP.Vect().Mag(); // MeV
  double p_pi = 0.;
  if ( event->NumFSParticle( PION_PLUS ) == 1 ) {
    p_pi = event->GetHMFSParticle( PION_PLUS )->fP.Vect().Mag(); // MeV
  }
  else if ( event->NumFSParticle( PION_MINUS ) == 1 ) {
    p_pi = event->GetHMFSParticle( PION_MINUS )->fP.Vect().Mag(); // MeV
  }

  // The muon momentum must be at least 150 MeV/c
  if ( p_mu <= 150.) return false;

  // The pion momentum must be higher than 100 MeV/c
  if ( p_pi <= 100. ) return false;

  //angle between the muon and pion must be less than 2.65 radians
  TVector3 p3_mu;
  if ( event->NumFSParticle( MUON ) == 1 ) {
    p3_mu = event->GetHMFSParticle( MUON )->fP.Vect();
  }
  else if ( event->NumFSParticle( ANTI_MUON ) == 1 ) {
    p3_mu = event->GetHMFSParticle( ANTI_MUON )->fP.Vect();
  }
  TVector3 p3_pi;
    if ( event->NumFSParticle( PION_PLUS ) == 1 ) {
        p3_pi = event->GetHMFSParticle( PION_PLUS )->fP.Vect();
    }
    else if ( event->NumFSParticle( PION_MINUS ) == 1 ) {
        p3_pi = event->GetHMFSParticle( PION_MINUS )->fP.Vect();
    }
    double angle_mu_pi = p3_mu.Angle( p3_pi );
    if ( angle_mu_pi >= 2.65 ) return false;


  // If we've made it here, then the current event has passed all of the
  // requirements in the signal definition
  return true;
}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::FillEventVariables( FitEvent* event ) {

  // Clear out the vector of passing bins, which may have already been filled
  // for the previous event
  fPassingBins.clear();

  if ( event->NumFSParticle(MU_MINUS) == 0 ) return;
  if ( event->NumFSParticle(PION_PLUS) == 0 && event->NumFSParticle(PION_MINUS) == 0 ) return;

  // Loop over each of the bin definitions. Keep track of the bins that
  // pass all cuts (and should thus be filled in this event)
  for ( size_t b = 0u; b < fBinDefinitions.size(); ++b ) {

    // Start out by assuming that the event passes all cuts
    bool passed_cuts = true;

    // Apply all cuts
    const auto& my_cuts = fBinDefinitions.at( b );
    for ( const auto& cut : my_cuts ) {
      passed_cuts &= cut.evaluate( event );
    }

    // If all cuts were passed, add this bin to the vector of indices for bins
    // that should be filled
    if ( passed_cuts ) fPassingBins.push_back( b );

  } // loop over bins

}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::LoadBinDefinitions() {
  std::string binning_file_name( FitPar::GetDataBase()
    + "/MicroBooNE/BNB_NumuCC1Pipm_2025/bin_defs.txt" );

  std::ifstream bin_file( binning_file_name );
  if ( !bin_file.is_open() ) {
    NUIS_ABORT( "Could not open bin definitions file " + binning_file_name );
  }
  std::string dummy_str;
  int bin_type, dummy_int;
  size_t num_true_bins = 0u;

  if (!(bin_file >> dummy_str >> dummy_str >> num_true_bins)) {
    NUIS_ABORT( "Could not read number of true bins from bin definitions file "
      + binning_file_name );
  }

  const std::string delimiter( "&&" );
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    int bin_type = 0;
    int dummy_int = 0;
    if (!(bin_file >> bin_type >> dummy_int)) {
      NUIS_ABORT( "Could not read bin type and block index from bin definitions file "
        + binning_file_name );
    }

    // Use two calls to std::getline using a double quote delimiter
    // in order to get the contents of the next double-quoted string
    std::string bin_def;
    if ( !std::getline(bin_file >> std::ws, bin_def, '\"') ) {
      NUIS_ABORT( "Could not read bin definition from bin definitions file "
        + binning_file_name );
    }
    if (!std::getline(bin_file, bin_def, '\"')) {
      NUIS_ABORT("Failed to read quoted bin definition in: " + binning_file_name);
    }

    // Skip bins of type == 1 (background true bins)
    if ( bin_type == 1 ) continue;

    auto cut_pos = bin_def.find(delimiter);
    if (cut_pos == std::string::npos) {
      NUIS_ABORT( "Could not find delimiter '&&' in bin definition: " + bin_def );
      continue;
    }

    // Skip the text before the first "&&" (here we assume that it is a simple
    // bool for the signal definition)
    std::string all_cuts = bin_def.substr(
      bin_def.find(delimiter) + delimiter.length() );

    // Create an empty vector of MyCutPipm objects to start defining the new bin
    fBinDefinitions.emplace_back();

    // Loop over the rest of the cuts and check each one
    size_t pos = 0u;
    do {
      // Advance to the next individual cut, removing it from the temporary
      // copy of the list of cuts
      pos = all_cuts.find( delimiter );
      std::string cut = all_cuts.substr( 0, pos );
      all_cuts.erase( 0, pos + delimiter.length() );

      // Parse the cut
      std::string var_name;
      std::string comp_op;
      double cut_val;
      std::stringstream cut_ss( cut );

      cut_ss >> var_name >> comp_op >> cut_val;

      // Interpret the variable for this cut and set up a function
      // object that will calculate it from the input FitEvent
      std::function< double(FitEvent*) > getter;

      if ( var_name == "mc_p3_pi.Mag()" ) {
        getter = [=]( FitEvent* ev ) -> double {
            double x = 0.;
            if ( ev->NumFSParticle( PION_PLUS ) == 1 ) {
                x = ev->GetHMFSParticle( PION_PLUS )->fP.Vect().Mag() / 1e3; // GeV
            }
            else if ( ev->NumFSParticle( PION_MINUS ) == 1 ) {
                x = ev->GetHMFSParticle( PION_MINUS )->fP.Vect().Mag() / 1e3; // GeV
            }
          return x;
        };
      }
      else if ( var_name == "mc_p3_pi.CosTheta()" ) {
        getter = [=]( FitEvent* ev ) -> double {
            double x = 0.;
            if ( ev->NumFSParticle( PION_PLUS ) == 1 ) {
                x = ev->GetHMFSParticle( PION_PLUS )->fP.Vect().CosTheta();
            }
            else if ( ev->NumFSParticle( PION_MINUS ) == 1 ) {
                x = ev->GetHMFSParticle( PION_MINUS )->fP.Vect().CosTheta();
            }
          return x;
        };
      }
      else if ( var_name == "mc_p3_mu.Mag()" ) {
        getter = [=]( FitEvent* ev ) -> double {
          return ev->GetHMFSParticle( MU_MINUS )->fP.Vect().Mag() / 1e3; // GeV
        };
      }
      else if ( var_name == "mc_p3_mu.CosTheta()" ) {
        getter = [=]( FitEvent* ev ) -> double {
          return ev->GetHMFSParticle( MU_MINUS )->fP.Vect().CosTheta();
        };
      }
      else if ( var_name == "mc_theta_mu_pi" ) {
        getter = [=]( FitEvent* ev ) -> double {
          const TVector3& p3mu = ev->GetHMFSParticle( MU_MINUS )->fP.Vect();
          TVector3 p3pi;
            if ( ev->NumFSParticle( PION_PLUS ) == 1 ) {
                p3pi = ev->GetHMFSParticle( PION_PLUS )->fP.Vect();
            }
            else if ( ev->NumFSParticle( PION_MINUS ) == 1 ) {
                p3pi = ev->GetHMFSParticle( PION_MINUS )->fP.Vect();
            }
          double denom = p3mu.Mag() * p3pi.Mag();
          if ( denom == 0. ) return 0.;
          double cosang = p3mu.Dot(p3pi) / denom;
          cosang = std::max( -1., std::min( 1., cosang ) );
          double theta_mupi = std::acos( cosang );
          return theta_mupi;
        };
      }
      else NUIS_ABORT( "Unrecognized cut variable " + var_name );

      // Test the cut based on the appropriate comparison operator
      std::function< bool(double) > tester;
      if ( comp_op == ">=" ) {
        tester = [cut_val]( double x ) -> bool {
          return x >= cut_val;
        };
      }
      else if ( comp_op == "<" ) {
        tester = [cut_val]( double x ) -> bool {
          return x < cut_val;
        };
      }
      else NUIS_ABORT( "Unrecognized comparison operator " + comp_op );

      // Add the finished cut to the definition for the current bin
      fBinDefinitions.back().emplace_back( getter, tester );

    } while ( pos != std::string::npos );

  } // loop over true bins

}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::FillHistograms() {

  if ( !Signal ) return;

  // The histograms have bin number as the x-axis variable, so just fill them
  // using the bin index of each bin that passed all cuts
  for ( const auto& bin : fPassingBins ) {
    NUIS_LOG(DEB, "Fill MCHist: " << bin << ", " << Weight);

    fMCHist->Fill( bin, Weight );
    fMCStat->Fill( bin, 1.0 );
    if ( fMCHist_Modes ) fMCHist_Modes->Fill( Mode, bin, Weight );

    fMCFine->Fill( bin, Weight );
    if ( fMCFine_Modes ) fMCFine_Modes->Fill( Mode, bin, Weight );
  }

}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::ConvertEventRates() {

  // Do the standard conversion
  Measurement1D::ConvertEventRates();

  // Clone the binning from the MC histogram
  fMCHistWithAC.reset(
    dynamic_cast< TH1D* >( fMCHist->Clone() )
  );
  fMCHistWithAC->SetNameTitle( (fSettings.GetName() + "_MC_with_AC").c_str(),
   fSettings.GetFullTitles().c_str() );
  fMCHistWithAC->Reset();
  fMCHistWithAC->SetDirectory( nullptr );

  // Build a column vector using the predicted cross sections
  int num_bins = fMCHist->GetNbinsX();
  TMatrixD pred( num_bins, 1 );
  for ( int b = 0; b < num_bins; ++b ) {
    pred( b, 0 ) = fMCHist->GetBinContent( b + 1 );
  }

  // Apply the additional smearing matrix A_C to create a new prediction
  TMatrixD new_pred( *fAddSmear, TMatrixD::kMult, pred );

  // Store it in our clone of the original MC prediction histogram
  for ( int b = 0; b < num_bins; ++b ) {
    double xsec = new_pred( b, 0 );
    fMCHistWithAC->SetBinContent( b + 1, xsec );
  }

  // Prepare the slice histograms now that everything else is ready
  this->PrepareSlices();
}

double MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::GetLikelihood() {

  if ( fNoData || !fDataHist ) return 0.;

  // Apply Masking to MC if Required.
  if ( fIsMask and fMaskHist ) {
    NUIS_ABORT("Bin masks not yet supported by"
      " the MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu sample" );
    //PlotUtils::MaskBins(fMCHist, fMaskHist);
  }

  // Likelihood Calculation
  double stat = 0.;
  if ( fIsChi2 ) {
    if ( fIsNS ) {
      NUIS_ABORT("Norm-shape covariance not yet supported by"
        " the MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu sample" );
    }
    stat = StatUtils::GetChi2FromCov( fDataHist, fMCHistWithAC.get(),
      covar, NULL, 1.0, 1.0, fIsWriting ? fResidualHist : NULL );
  }

  fLikelihood = stat;

  return stat;
}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::Write( std::string drawOpt ) {
  // Write the standard information
  Measurement1D::Write( drawOpt );
  MicroBooNEHelper::WriteResults(fSettings.GetName(), fDataHist, fMCHistWithAC.get(),
    covar, fBlockHandler, fMCHist_Slices, fMCModeHists_Slices, fMCHist_Modes);

}

void MicroBooNE_BNB_NumuCC1Pipm_2025_XSec_nu::PrepareSlices() {
  // Load the definitions of the blocks of bins
  std::string block_file_name( FitPar::GetDataBase()
    + "/MicroBooNE/BNB_NumuCC1Pipm_2025/bin_blocks.txt" );

  fBlockHandler = std::make_shared< MicroBooNEBlockHandler >(
    fSettings.GetName(), block_file_name );

  // Clone the data histograms owned by the block handler to create
  // corresponding MC histograms
  for ( auto& hist : fBlockHandler->fHists ) {
    // Attach the MC suffix to the clone's name, then attach the data suffix
    // to the original
    auto mc_hist_name = hist->GetName() + std::string( "_MC" );
    hist->SetName( (hist->GetName() + std::string( "_data" )).c_str() );
    auto* temp_mc_hist = dynamic_cast< TH1D* >(
      hist->Clone( mc_hist_name.c_str() )
    );
    temp_mc_hist->SetDirectory( nullptr );
    temp_mc_hist->Reset();

    temp_mc_hist->SetLineColor( kRed );
    temp_mc_hist->SetLineStyle( 1 );
    temp_mc_hist->SetLineWidth( 1 );

    temp_mc_hist->SetFillColor( 0 );
    temp_mc_hist->SetFillStyle( 1001 );

    fMCHist_Slices.emplace_back( temp_mc_hist );

    // Make MC slice histograms for individual scattering modes
    for ( int m : fMCHist_Modes->fmodes ) {
      // Convert the interaction mode integer into a string label
      auto mode = static_cast< InputHandler::InteractionModes >( m );
      std::ostringstream oss;
      oss << '_' << mode;
      auto mode_hist_name = temp_mc_hist->GetName() + oss.str();
      // Set up the slice histogram for this mode
      auto* temp_mode_hist = dynamic_cast< TH1D* >(
        temp_mc_hist->Clone( mode_hist_name.c_str() )
      );
      temp_mode_hist->SetDirectory( nullptr );
      temp_mode_hist->Reset();
      fMCModeHists_Slices[ m ].emplace_back( temp_mode_hist );
    }

  }

  // Now that all slice histograms have been created, populate them with
  // corresponding entries from the full measurement histogram (that is
  // organized in terms of global bin number).
  // NOTE: Bins with one or more infinite edges are automatically skipped
  // by the block handler, so we don't have to worry about bins that have
  // no match in the slices
  for ( const auto& [mk, mv] : fBlockHandler->fBinMap ) {
    // The ROOT histogram used to store the global bin contents has one-based
    // indices, but the input table used to define the blocks is zero-based
    int global_bin = mk + 1;

    // Zero-based index of the slice histogram to which this global bin belongs
    int hist_idx = mv.histIdx; // position in the vector of slices

    // One-based bin index to be populated in the target histogram
    int local_bin = mv.rootBin;

    // Product of finite bin widths for this bin (used to convert to a
    // differential cross section)
    double widths = mv.width;

    // Retrieve the values of interest from the global histograms, dividing
    // by the bin width(s) to obtain a differential xsec each time
    double data_val = fDataHist->GetBinContent( global_bin ) / widths;
    double data_err = fDataHist->GetBinError( global_bin ) / widths;
    double mc_val = fMCHistWithAC->GetBinContent( global_bin ) / widths;
    std::map< int, double > mode_vals;
    for ( int m : fMCHist_Modes->fmodes ) {
      int mode_idx = fMCHist_Modes->ConvertModeToIndex( m );
      auto* mode_hist = fMCHist_Modes->GetHist( mode_idx );
      double value = mode_hist->GetBinContent( global_bin );
      mode_vals[ m ] = value / widths;
    }

    // Set the bin contents in the target slice histograms
    auto& data_hist = fBlockHandler->fHists.at( hist_idx );
    data_hist->SetBinContent( local_bin, data_val );
    data_hist->SetBinError( local_bin, data_err );

    auto& mc_hist = fMCHist_Slices.at( hist_idx );
    mc_hist->SetBinContent( local_bin, mc_val );

    for ( const auto& [ mode, val ] : mode_vals ) {
      auto& mode_hist = fMCModeHists_Slices.at( mode ).at( hist_idx );
      mode_hist->SetBinContent( local_bin, val );
    }

  }
}