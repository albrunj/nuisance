// MicroBooNEHelper.cxx

#include "MicroBooNEHelper.h"
#include "MicroBooNEBlockHandler.h"
#include "InteractionModes.h"
#include "StandardStacks.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "StatUtils.h"
#include "TNamed.h"

#include "TH1D.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

namespace MicroBooNEHelper {
  namespace {


  } // namespace

  // Loads a matrix from a text file following the format used
  // in the data release for this MicroBooNE measurement
  TMatrixD load_matrix( const std::string& input_file_name ) {

    // Get the table of matrix element values
    std::ifstream matrix_table_file( input_file_name );
    if ( !matrix_table_file.is_open() ) {
    throw std::runtime_error( "MicroBooNEHelper::load_matrix: cannot open "
      + input_file_name );
    }

    // Peek at the file contents to decide whether we're working with
    // a matrix or a column vector
    std::string dummy;
    matrix_table_file >> dummy >> dummy >> dummy;
    bool is_matrix = ( dummy == "numYbins" );

    // Return to the beginning of the file for parsing
    matrix_table_file.seekg( 0 );

    // Get the matrix or vector dimensions from the header line(s)
    int num_x_bins, num_y_bins;

    matrix_table_file >> dummy >> num_x_bins;
    if ( is_matrix ) {
      matrix_table_file >> dummy >> num_y_bins;

      // Skip the next header line which contains the data column names
      std::getline( matrix_table_file, dummy );
    }
    else {
      num_y_bins = 1;
    }

    // Create a TMatrixD with the correct dimensions
    TMatrixD matrix( num_x_bins, num_y_bins );

    // Parse its contents from the remaining lines
    std::string line;
    while ( std::getline(matrix_table_file, line) ) {
      int bin1, bin2;
      double element;

      std::stringstream temp_ss( line );
      temp_ss >> bin1;
      if ( is_matrix ) {
        temp_ss >> bin2;
      }
      else {
        bin2 = 0;
      }
      temp_ss >> element;

      if ( bin1 < num_x_bins && bin2 < num_y_bins ) {
        matrix( bin1, bin2 ) = element;
      }
    }

    return matrix;
  }

  // Helper function that converts a TMatrixD into the TMatrixDSym*
  // needed to initialize some class members inherited from Measurement1D.
  // Assumes that the input is a square matrix (TODO: add error handling).
  TMatrixDSym* to_symmetric_matrix( const TMatrixD& mat ) {
    int num_rows = mat.GetNrows();
    TMatrixDSym* sym = new TMatrixDSym( num_rows );
    for ( int a = 0; a < num_rows; ++a ) {
      for ( int b = 0; b < num_rows; ++b ) {
        sym->operator()( a, b ) = mat( a, b );
      }
    }
    return sym;
  }

  // Helper function that creates a TH1D from a TMatrixD. Currently
  // assumes that the input matrix has a single column (TODO: Add
  // error handling)
  TH1D* to_histogram( const TMatrixD& vec ) {
    int num_rows = vec.GetNrows();
    auto* hist = new TH1D( "vec_hist", "", num_rows, 0., num_rows );
    for ( int a = 0; a < num_rows; ++a ) {
      double value = vec( a, 0 );
      hist->SetBinContent( a + 1, value ); // ROOT bin indices are one-based
    }
    return hist;
  }

  //take in results and write out standard information
  void WriteResults(const std::string& sampleName, TH1D* data_hist, TH1D* mc_hist_with_ac, TMatrixDSym* covar,
    std::shared_ptr<MicroBooNEBlockHandler> block_handler, std::vector<std::shared_ptr<TH1D>>& mc_slice_hists, const std::map<int,
    std::vector<std::shared_ptr<TH1D>>>& mc_mode_slice_hists, TrueModeStack* mc_modes){

    // To start extra chi-squared calculations, create a mask that excludes
    // everything. Note that bins are kept if the mask entry is zero, and
    // discarded otherwise.
    int num_bins = data_hist->GetNbinsX();
    TH1I mask_all( "mask_all", "mask all", num_bins, 0., num_bins );
    mask_all.SetDirectory( nullptr );
    for ( int b = 1; b <= num_bins; ++b ) mask_all.SetBinContent( b, 1 );

    // Compute chi-squared scores for each block of bins separately using
    // masks. Store the results in TNamed objects in the output
    // file.
    std::shared_ptr< TH1I > block_mask;
    for ( const auto& [ block_idx, bin_vec ] : block_handler->fBlockBins ) {
      // Clone the exclude-everything mask and zero out bins that belong
      // to the current block (thus including them)
      block_mask.reset(
        dynamic_cast< TH1I* >( mask_all.Clone("block_mask") )
      );

      // Global bin indices are zero-based, while the ROOT histogram bins
      // are one-based, so we correct for this here
      for ( int b : bin_vec ) block_mask->SetBinContent( b + 1, 0 );

      // Compute the chi-squared statistic for the current block using the mask
      double chi2 = StatUtils::GetChi2FromCov( data_hist, mc_hist_with_ac,
        covar, block_mask.get(), 1.0, 1.0, nullptr );

      // Create a TNamed object to store the result. We use a TNamed so that we
      // can store the chi^2 value and number of bins together for convenient
      // viewing
      std::string param_name = sampleName + "_Block"
        + std::to_string( block_idx ) + "_Chi2";
      std::string chi2_per_bin_str = std::to_string( chi2 )
        + " / " + std::to_string( bin_vec.size() );
      TNamed chi2_param( param_name.c_str(), chi2_per_bin_str.c_str() );

      chi2_param.Write();
    }

    // Now we follow a similar procedure to compute and store chi-squared values
    // for each individual slice histogram in multi-slice blocks
    std::shared_ptr< TH1I > slice_mask;
    int num_slices = block_handler->fHists.size();
    for ( int s = 0u; s < num_slices; ++s ) {

      // Before doing anything else, write the histograms themselves for the
      // current slice to the output file
      block_handler->fHists.at( s )->Write(); // data slice histogram
      mc_slice_hists.at( s )->Write(); // total MC slice histogram
      for ( int m : mc_modes->fmodes ) {
        // MC slice histogram for mode m
        mc_mode_slice_hists.at( m ).at( s )->Write();
      }

      // Find the block and slice number within the block for the current slice
      // histogram
      int block_idx = -1; // dummy value
      int slice_idx = -1; // dummy value
      for ( const auto& [ bl, hist_idx_vec ] : block_handler->fBlockHists ) {
        for ( size_t h = 0u; h < hist_idx_vec.size(); ++h ) {
          int h_idx = hist_idx_vec.at( h );
          if ( h_idx == s ) {
            block_idx = bl;
            slice_idx = h;
          }
        }
      }

      // NOTE: This continue statement has been removed since some 1D blocks
      // include underflow or overflow bins that don't show up in their
      // slice histogram. In some cases, the user may want to see both
      // values, so it is easiest just to output all the information.
      //
      // OLD: If this is a 1D block, then there is only a single slice, and we
      // have already computed a suitable chi^2 score in the block-by-block
      // results above. We can therefore skip to the next slice.
      //if ( blockHandler->fBlockHists.size() <= 1u ) continue;

      // Clone the exclude-everything mask and zero out bins that belong
      // to the current slice (thus including them)
      slice_mask.reset(
        dynamic_cast< TH1I* >( mask_all.Clone("slice_mask") )
      );

      int num_unmasked_bins = 0;
      for ( const auto& [ mk, mv ] : block_handler->fBinMap ) {
        // Unmask only bins that belong to the current slice histogram
        int hist_index = mv.histIdx;
        if ( hist_index != s ) continue;
        // Do the unmasking, taking into account the one-based indexing
        // of the ROOT histograms and the zero-based bin indices from the table
        // of block definitions
        int global_bin = mk;
        slice_mask->SetBinContent( global_bin + 1, 0 );
        ++num_unmasked_bins;
      }

      // Compute the chi-squared statistic for the current block using the mask
      double chi2 = StatUtils::GetChi2FromCov( data_hist, mc_hist_with_ac,
        covar, slice_mask.get(), 1.0, 1.0, nullptr );

      // Create a TNamed object to store the result. We use a TNamed so that we
      // can store the chi^2 value and number of bins together for convenient
      // viewing
      std::string param_name = sampleName + "_Block"
        + std::to_string( block_idx ) + "_Slice"
        + std::to_string( slice_idx ) + "_Chi2";
      std::string chi2_per_bin_str = std::to_string( chi2 )
        + " / " + std::to_string( num_unmasked_bins );
      TNamed chi2_param( param_name.c_str(), chi2_per_bin_str.c_str() );

      chi2_param.Write();
    }
  }

} // namespace MicroBooNEHelper
