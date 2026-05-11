// MicroBooNEHelper.h
//
// assumes use of MicroBooNEBlockHandler and measurement structure like NumuCC1Pipm_2025_XSec_nu
//
// Helper functions used in different MicroBooNE measurements
//
// load_matrix:           load matrix from MicroBooNE data
//
// to_symmetric_matrix:   convert matrix to TMatrixDSym
//
// to_histogram:          convert matrix to TH1D 
//
// Write:                 take in standard information and write out



#ifndef MICROBOONE_HELPER_H_SEEN
#define MICROBOONE_HELPER_H_SEEN

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH1I.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "InteractionModes.h"
#include "StandardStacks.h"

#include "MicroBooNEBlockHandler.h"

class TH1D;
class SampleSettings;

namespace MicroBooNEHelper {

  // Load a matrix or a column vector from the MicroBooNE text-table format.
  //
  // Supported formats:
  //   numXbins <N> numYbins <M>  
  //   <binX> <binY> <value>
  //
  // and
  //   numXbins <N>
  //   <binX> <value>
  //
  // The returned matrix has dimensions N x M, with M = 1 for vectors.
  TMatrixD load_matrix( const std::string& input_file_name );

  // Convert a square TMatrixD into a TMatrixDSym.
  // The caller owns the returned matrix.
  TMatrixDSym* to_symmetric_matrix( const TMatrixD& mat );

  // Convert a single-column TMatrixD into a TH1D.
  // The caller owns the returned histogram.
  TH1D* to_histogram( const TMatrixD& vec );

  //take in results and write out standard information
  void WriteResults(const std::string& sampleName, TH1D* data_hist, TH1D* mc_hist_with_ac, TMatrixDSym* covar,
    std::shared_ptr<MicroBooNEBlockHandler> block_handler, std::vector<std::shared_ptr<TH1D>>& mc_slice_hists, const std::map<int,
    std::vector<std::shared_ptr<TH1D>>>& mc_mode_slice_hists, TrueModeStack* mc_modes);

} // namespace MicroBooNEHelper

#endif // MICROBOONE_HELPER_H_SEEN
