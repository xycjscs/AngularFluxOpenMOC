#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "Quadrature.h"


/**
 * @brief Constructor for VolumeCalculator calls the TraverseTracks
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information
 */
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseTracks(track_generator) {
}


/**
 * @brief FSR volumes are calculated and saved in the TrackGenerator's FSR
 *        volumes buffer
 * @details VolumeKernels are created and used to loop over all segments and
 *          tally each segments contribution to FSR volumes.
 */
void VolumeCalculator::execute() {
#pragma omp parallel
  {
    VolumeKernel kernel(_track_generator);
    loopOverTracks(&kernel);
  }
}


/**
 * @brief Constructor for TransportSweep calls the TraverseTracks
 *        constructor and allocates temporary memory for local scalar fluxes
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweep::TransportSweep(TrackGenerator* track_generator)
                              : TraverseTracks(track_generator) {
  _cpu_solver = NULL;

  /* Allocate temporary storage of FSR fluxes */
  int num_threads = omp_get_max_threads();
  int num_groups = track_generator->getGeometry()->getNumEnergyGroups();
  _thread_fsr_fluxes = new FP_PRECISION*[num_threads];

  /* Allocate fluxes for the number of groups plus some extra space to prevent
   * false sharing conflicts */
  int array_width = num_groups + 8;
  for (int i=0; i < num_threads; i++)
    _thread_fsr_fluxes[i] = new FP_PRECISION[array_width];
}


/**
 * @brief Destructor deletes temporary storage of local scalar fluxes
 */
TransportSweep::~TransportSweep() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++)
    delete [] _thread_fsr_fluxes[i];
  delete [] _thread_fsr_fluxes;
}


/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweep::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
  }
}


/**
 * @brief Sets the CPUSolver so that TransportSweep can apply MOC equations
 * @details This allows TransportSweep to transfer boundary fluxes from the
 *          CPUSolver and tally scalar fluxes
 * @param cpu_solver The CPUSolver which applies the MOC equations
 */
void TransportSweep::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Applies the MOC equations the Track and segments
 * @details The MOC equations are applied to each segment, attenuating the
 *          Track's angular flux and tallying FSR contributions. Finally,
 *          Track boundary fluxes are transferred.
 * @param track The Track for which the angular flux is attenuated and
 *        transferred
 * @param segments The segments over which the MOC equations are applied
 */
void TransportSweep::onTrack(Track* track, segment* segments) {

  /* Get array for temporary scalar flux storage */
  int tid = omp_get_thread_num();
  FP_PRECISION* thread_fsr_flux = _thread_fsr_fluxes[tid];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimAngleIndex();
  int num_segments = track->getNumSegments();
  FP_PRECISION* track_flux, partial_current;

  /* Correct azimuthal index to first octant */
  Quadrature* quad = _track_generator->getQuadrature();

  /* Get the forward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, true);

  /* Make a map of FSRs to cells */
  std::map<int,Cell*>& map_fsr_to_cell = _track_generator->getGeometry()->getMapFSRstoCells(); 

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {
    segment* curr_segment = &segments[s];
    segment* next_segment = &segments[std::min(s+1, num_segments-1)];
    
    log_printf(DEBUG, "Doing segment   %d in region %d", s, curr_segment->_region_id);
    log_printf(DEBUG, "Next segment is %d in region %d", s, curr_segment->_region_id);
    
    _cpu_solver->tallyScalarFlux(curr_segment, next_segment, azim_index, track_flux,
                                 thread_fsr_flux, map_fsr_to_cell);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, true, track_flux);

  /* Get the backward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, false);

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {
    segment* curr_segment = &segments[s];
    segment* next_segment = &segments[std::max(s-1, 0)];
    _cpu_solver->tallyScalarFlux(curr_segment, next_segment, azim_index, track_flux,
                                 thread_fsr_flux, map_fsr_to_cell);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, false);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, false, track_flux);
}
