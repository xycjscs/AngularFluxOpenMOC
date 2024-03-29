#include "CPUSolver.h"
#include <iostream>
using namespace std;					
/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param track_generator an optional pointer to the TrackGenerator
 */
CPUSolver::CPUSolver(TrackGenerator* track_generator)
  : Solver(track_generator) {

  setNumThreads(1);
  _FSR_locks = NULL;
  
  /* Initalizing to NULL to be able to check if they have been filled before */
  _reference_partial_currents = NULL;
  _reference_partial_currents_length = NULL;
  _current_start_row_index = NULL; // which index the currents in ref_partial_currents start for each cell_from 
  _current_column_index = NULL;    // which cell each element in ref_partial_currents is going to  
  
}


/**
 * @brief Fills an array of reference partial currents and two indexing arrays
 * @details This class method fills 
 *  - the array containing the reference partial currents,
 *  - the array containing the index of the row where currents from a given cell start
 *  - the array containing the cell_to of these currents
 * @param cell_from, cell from which the partial current comes from
 * @param cell_to, cell to which the partial current goes
 * @param p the polar angle index
 * @param group, the energy group of this partial current
 * @param ref_current, the value of the partial current
 */
void CPUSolver::setReferencePartialCurrents(int cell_from, int cell_to, int group, int p, double ref_current){

  /* Create row : cell_from indexes */
  int ii = 0;

  /* If the index of the beginning of the row hasn't been set yet */
  if(_current_start_row_index[cell_from] == 0){

    /* first time we fill in the currents is with cell_from = 0 */
    if(cell_from == 0){
      _current_start_row_index[cell_from] = 0;
    }
    /* not first time, let s look for when previous row ends */
    else{
      while(_current_column_index[_current_start_row_index[cell_from - 1] + ii] != -1){ 
          ii += 1;
      }
    _current_start_row_index[cell_from] = _current_start_row_index[cell_from - 1] + ii;
    }
  }

  /* to avoid a bug for other currents */
  if(ref_current == 0){
    ref_current = 1e-15; 
    log_printf(NORMAL, "Zero current replaced by 1e-15");
  }
  
  // Create column indexes too
  int row_start = _current_start_row_index[cell_from];
  int column_index;

  /* Loop through _reference_partial_currents vector and fill first zero positions */
  for(int ii = row_start; ii < 2 * _num_surfaces; ii++){

    if(_reference_partial_currents[ii][group * _num_polar_2 + p] == 0){
      _reference_partial_currents[ii][group * _num_polar_2 + p] = ref_current;
      _reference_partial_currents_length[ii][group * _num_polar_2 + p] = ref_current;
	  
      log_printf(DEBUG, "Set current for %d -> %d, group %d, polar %d", cell_from,
             cell_to, group, p);

      /* If first reference current to be set for this couple,
      save the cell_to for the column index */
      if(_current_column_index[ii] == -1){
        _current_column_index[ii] = cell_to;
      }
      /* If not first reference, and value in array is different, there's a bug */
      else if (_current_column_index[ii] != cell_to){
        log_printf(ERROR,  "wrong cell_to in column_index");
        abort();
      }
          
      column_index = ii - row_start;
	  std::cout << "cell from->to "<<cell_from<< " -> "<< cell_to<< "  index: "<< column_index +  row_start << std::endl;
      break;
    }
  }

  // Print all 3 arrays to check 
   //log_printf(DEBUG, "index of row start %d, and cell_to %d", row_start, column_index);
   //std::cout << "Start row index : ";
   //for(int ii = 0; ii < _num_FSRs; ii++){std::cout << _current_start_row_index[ii] << " ";}
   //std::cout << std::endl;
   //std::cout << "Cell to for each column : ";
   //for(int ii = 0; ii < 2*_num_surfaces; ii++){std::cout << _current_column_index[ii] << " ";}
   //std::cout << std::endl;

}


/**
 * @brief Access a vector of reference partial currents for a [cell_from, cell_to] couple
 * @details 
 * @param cell_from, cell from which the partial current comes from
 * @param cell_to, cell to which the partial current goes
 * @param index, pointer to the index of the partial current in the currents array
          saved to access the other partial current arrays without another search
 */
double* CPUSolver::getReferencePartialCurrents(int cell_from, int cell_to, int* index){
  
  int row_start = _current_start_row_index[cell_from];  //offset here if bad cell numbers
 
  /* Find column index */
  int column_index = 0;
  int row_end = _current_start_row_index[cell_from + 1] - 1;  //offset here if bad cell numbers
  for(int ii = row_start; ii < 2*_num_surfaces; ii++){
    //std::cout << _current_column_index[ii] << std::endl;
    if(_current_column_index[ii] == cell_to){ //FIXME  //offset here if bad cell numbers
      break;
    }
    column_index += 1;
  }
  *index = row_start + column_index;
   
  //log_printf(NORMAL, "cell from->to %d -> %d , index %d", cell_from, cell_to, row_start + column_index);
  /*log_printf(NORMAL, "cell from->to %d -> %d , index %d", cell_from, cell_to, *index);
  log_printf(NORMAL, "row_start %d  column_index %d", row_start, column_index);
  */
  return _reference_partial_currents[row_start + column_index];
}

double* CPUSolver::getReferencePartialCurrentsLength(int cell_from, int cell_to, int* index){
  
  int row_start = _current_start_row_index[cell_from];  //offset here if bad cell numbers
 
  /* Find column index */
  int column_index = 0;
  int row_end = _current_start_row_index[cell_from + 1] - 1;  //offset here if bad cell numbers
  for(int ii = row_start; ii < 2*_num_surfaces; ii++){
    //std::cout << _current_column_index[ii] << std::endl;
    if(_current_column_index[ii] == cell_to){ //FIXME  //offset here if bad cell numbers
      break;
    }
    column_index += 1;
  }
  *index = row_start + column_index;
   
  //log_printf(NORMAL, "cell from->to %d -> %d , index %d", cell_from, cell_to, row_start + column_index);
  /*log_printf(NORMAL, "cell from->to %d -> %d , index %d", cell_from, cell_to, *index);
  log_printf(NORMAL, "row_start %d  column_index %d", row_start, column_index);
  */
  return _reference_partial_currents_length[row_start + column_index];
}


/**
 * @brief Get a single value of the last iteration partial currents for a 
 * [surface_index, energy group, polar angle] combination
 */
double CPUSolver::getOngoingPartialCurrent(int index, int group, int azim, int p){
  return _ongoing_partial_currents[index][group * _num_polar_2 * _num_azim + azim * _num_polar_2 + p];
}

double CPUSolver::getOngoingPartialCurrentLength(int index, int group, int azim, int p){
  return _ongoing_partial_currents_length[index][group * _num_polar_2 * _num_azim + azim * _num_polar_2 + p];
}


/**
 * @brief Get a single value of the last iteration partial currents for a 
 * [cell_from, cell_to, energy group, polar angle] combination
 */
double CPUSolver::getAngularPartialCurrent(int cell_from, int cell_to, int group, int azim, int p){

  int index_partial_current;
  double current;
  if(cell_from == cell_to){
    log_printf(NORMAL, "current from a cell to same cell set to 0");
    current = 0.;
  }
  else{
    /* Get the surface index */
    double* ref_partial_current;
    ref_partial_current = getReferencePartialCurrents(cell_from, cell_to, 
                                                      &index_partial_current);

    /* Call the routine using the surface index */
    current = getOngoingPartialCurrent(index_partial_current, group, azim, p);
  }
  //std::cout << "cell from->to "<<cell_from<< " -> "<< cell_to<< "  index: "<< index_partial_current << std::endl;
  return current;
}

double CPUSolver::getAngularPartialCurrentLength(int cell_from, int cell_to, int group, int azim, int p){

  int index_partial_current;
  double current;
  if(cell_from == cell_to){
    log_printf(NORMAL, "current from a cell to same cell set to 0");
    current = 0.;
  }
  else{
    /* Get the surface index */
    double* ref_partial_current;
    ref_partial_current = getReferencePartialCurrentsLength(cell_from, cell_to, 
                                                      &index_partial_current);

    /* Call the routine using the surface index */
    current = getOngoingPartialCurrentLength(index_partial_current, group, azim, p);
  }
  std::cout << "cell from->to "<<cell_from<< " -> "<< cell_to<< "  index: "<< index_partial_current << std::endl;
  return current;
}

/**
 * @brief Replaces default value of num_surfaces 
 * @details This class method allows to use the exact number of surfaces output by another code
 * @param number_surfaces, new value of _num_surfaces
 */
void CPUSolver::setNumSurfaces(int number_surfaces){
    _num_surfaces = number_surfaces;
}


/**
 * @brief Resets array of partial currents
 */
void CPUSolver::resetOngoingPartialCurrentsArray(){
  for(int ii=0; ii < 2 * _num_surfaces; ii++){
      memset(_ongoing_partial_currents[ii], 0.0, sizeof(double) * _num_groups * _num_azim * _num_polar_2);
  }
}

void CPUSolver::resetOngoingPartialCurrentsLengthArray(){
  for(int ii=0; ii < 2 * _num_surfaces; ii++){
      memset(_ongoing_partial_currents_length[ii], 0.0, sizeof(double) * _num_groups * _num_azim * _num_polar_2);
  }
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
  return _num_threads;
}


/**
 * @brief Fills an array with the scalar fluxes.
 * @details This class method is a helper routine called by the OpenMOC
 *          Python "openmoc.krylov" module for Krylov subspace methods.
 *          Although this method appears to require two arguments, in
 *          reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_fluxes = num_groups * num_FSRs
 *          fluxes = solver.getFluxes(num_fluxes)
 * @endcode
 *
 * @param fluxes an array of FSR scalar fluxes in each energy group
 * @param num_fluxes the total number of FSR flux values
 */
void CPUSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {

  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since there are "
               "%d groups and %d FSRs which does not match the requested "
               "%d flux values", _num_groups, _num_FSRs, num_fluxes);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since they "
               "have not yet been allocated");

  /* If the user called setFluxes(...) they already have the flux */
  if (_user_fluxes && _scalar_flux == out_fluxes)
    return;

  /* Otherwise, copy the fluxes into the input array */
  else {
#pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        out_fluxes[r*_num_groups+e] = _scalar_flux(r,e);
    }
  }
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @param num_threads the number of threads
 */
void CPUSolver::setNumThreads(int num_threads) {

  if (num_threads <= 0)
    log_printf(ERROR, "Unable to set the number of threads to %d "
               "since it is less than or equal to 0", num_threads);

  /* Set the number of threads for OpenMP */
  _num_threads = num_threads;
  omp_set_num_threads(_num_threads);
}


/**
 * @brief Set the flux array for use in transport sweep source calculations.
 * @detail This is a helper method for the checkpoint restart capabilities,
 *         as well as the IRAMSolver in the openmoc.krylov submodule. This
 *         routine may be used as follows from within Python:
 *
 * @code
 *          fluxes = numpy.random.rand(num_FSRs * num_groups, dtype=np.float)
 *          solver.setFluxes(fluxes)
 * @endcode
 *
 *          NOTE: This routine stores a pointer to the fluxes for the Solver
 *          to use during transport sweeps and other calculations. Hence, the
 *          flux array pointer is shared between NumPy and the Solver.
 *
 * @param in_fluxes an array with the fluxes to use
 * @param num_fluxes the number of flux values (# groups x # FSRs)
 */
void CPUSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {
  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to set an array with %d flux values for %d "
               " groups and %d FSRs", num_fluxes, _num_groups, _num_FSRs);

  /* Allocate array if flux arrays have not yet been initialized */
  if (_scalar_flux == NULL)
    initializeFluxArrays();

  /* Set the scalar flux array pointer to the array passed in from NumPy */
  _scalar_flux = in_fluxes;
  _user_fluxes = true;
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method gets an array of OpenMP mutual exclusion locks
 *          for each FSR for use in the transport sweep algorithm.
 */
void CPUSolver::initializeFSRs() {
  Solver::initializeFSRs();
  _FSR_locks = _track_generator->getFSRLocks();
}


/**
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPUSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_start_flux != NULL)
    delete [] _start_flux;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  /* Allocate memory for the Track boundary flux arrays */
  try {
    int size = 2 * _tot_num_tracks * _polar_times_groups;
    _boundary_flux = new FP_PRECISION[size];
    _start_flux = new FP_PRECISION[size];

    /* Allocate an array for the FSR scalar flux */
    size = _num_FSRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
    _old_scalar_flux = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the fluxes");
  }
}


/**
 * @brief Allocates memory for reference, previous and ongoing Partial Currents
 * @details Deletes memory for old current arrays if they were allocated
 *          for a previous simulation.
 */
void CPUSolver::initializePartialCurrentArrays() {

  _num_FSRs = _geometry->getNumFSRs();
  _num_groups = _geometry->getNumEnergyGroups();
 
  if (_current_start_row_index != NULL)
    delete [] _current_start_row_index;
  
  if (_current_column_index != NULL)
    delete [] _current_column_index;
  
  /* Allocate arrays for the reference partial currents */
  log_printf(INFO, "Setting arrays for %d currents, %d fsrs, %d energy groups",
        2*_num_surfaces, _num_FSRs, _num_groups);

  _current_start_row_index = new int[_num_FSRs]();
  _current_column_index = new int[2*_num_surfaces]();
  
  
  /* Initialize row indexes as 0, to error if problem during initialization */
  _current_start_row_index = new int[_num_FSRs];
  for(int ii = 0; ii < _num_FSRs; ii++){_current_start_row_index[ii] = 0;}
  /* Initialize column indexes as -1, to error if problem during initialization */
  _current_column_index = new int[2*_num_surfaces];
  for(int ii = 0; ii < 2*_num_surfaces; ii++){_current_column_index[ii] = -1;}
  
  _ongoing_partial_currents = new double*[2*_num_surfaces];
  for(int ii=0; ii < 2*_num_surfaces; ii++){
      _ongoing_partial_currents[ii] = new double[_num_groups*_num_polar_2*_num_azim]();
  }
  _ongoing_partial_currents_length = new double*[2*_num_surfaces];
  for(int ii=0; ii < 2*_num_surfaces; ii++){
      _ongoing_partial_currents_length[ii] = new double[_num_groups*_num_polar_2*_num_azim]();
  }
  
  
  _reference_partial_currents = new double*[2*_num_surfaces];
  for(int ii=0; ii < 2*_num_surfaces; ii++){
      _reference_partial_currents[ii] = new double[_num_groups*_num_polar_2*_num_azim]();
  }
  _reference_partial_currents_length = new double*[2*_num_surfaces];
  for(int ii=0; ii < 2*_num_surfaces; ii++){
      _reference_partial_currents_length[ii] = new double[_num_groups*_num_polar_2*_num_azim]();
  }
}


/**
 * @brief Allocates memory for FSR source arrays.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void CPUSolver::initializeSourceArrays() {

  /* Delete old sources arrays if they exist */
  if (_reduced_sources != NULL)
    delete [] _reduced_sources;
  if (_fixed_sources != NULL)
    delete [] _fixed_sources;

  int size = _num_FSRs * _num_groups;

  /* Allocate memory for all source arrays */
  try {
    _reduced_sources = new FP_PRECISION[size];
    _fixed_sources = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR sources");
  }

  /* Initialize fixed sources to zero */
  memset(_fixed_sources, 0.0, sizeof(FP_PRECISION) * size);

  /* Populate fixed source array with any user-defined sources */
  initializeFixedSources();
}


/**
 * @brief Populates array of fixed sources assigned by FSR.
 */
void CPUSolver::initializeFixedSources() {

  Solver::initializeFixedSources();

  int fsr_id, group;
  std::pair<int, int> fsr_group_key;
  std::map< std::pair<int, int>, FP_PRECISION >::iterator fsr_iter;

  /* Populate fixed source array with any user-defined sources */
  for (fsr_iter = _fix_src_FSR_map.begin();
       fsr_iter != _fix_src_FSR_map.end(); ++fsr_iter) {

    /* Get the FSR with an assigned fixed source */
    fsr_group_key = fsr_iter->first;
    fsr_id = fsr_group_key.first;
    group = fsr_group_key.second;

    if (group <= 0 || group > _num_groups)
      log_printf(ERROR,"Unable to use fixed source for group %d in "
                 "a %d energy group problem", group, _num_groups);

    if (fsr_id < 0 || fsr_id >= _num_FSRs)
      log_printf(ERROR,"Unable to use fixed source for FSR %d with only "
                 "%d FSRs in the geometry", fsr_id, _num_FSRs);

    _fixed_sources(fsr_id, group-1) = _fix_src_FSR_map[fsr_group_key];
  }
}


/**
 * @brief Zero each Track's boundary fluxes for each energy group
 *        and polar angle in the "forward" and "reverse" directions.
 */
void CPUSolver::zeroTrackFluxes() {

#pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar_2; p++) {
        for (int e=0; e < _num_groups; e++) {
          _boundary_flux(t,d,p,e) = 0.0;
          _start_flux(t,d,p,e) = 0.0;
        }
      }
    }
  }
}


/**
 * @brief Copies values from the start flux into the boundary flux array
 *        for both the "forward" and "reverse" directions.
 */
void CPUSolver::copyBoundaryFluxes() {

  #pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar_2; p++) {
        for (int e=0; e < _num_groups; e++) {
          _boundary_flux(t, d, p, e) = _start_flux(t, d, p, e);
        }
      }
    }
  }
}


/**
 * @brief Set the scalar flux for each FSR and energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(r,e) = value;
  }
}


/**
 * @brief Stores the FSR scalar fluxes in the old scalar flux array.
 */
void CPUSolver::storeFSRFluxes() {

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _old_scalar_flux(r,e) = _scalar_flux(r,e);
  }
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void CPUSolver::normalizeFluxes() {

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;
  FP_PRECISION tot_fission_source;
  FP_PRECISION norm_factor;

  int size = _num_FSRs * _num_groups;
  FP_PRECISION* fission_sources = new FP_PRECISION[_num_FSRs * _num_groups];

  /* Compute total fission source for each FSR, energy group */
#pragma omp parallel for private(volume, nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    /* Get pointers to important data structures */
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
    volume = _FSR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
  }

  /* Compute the total fission source */
  tot_fission_source = pairwise_sum<FP_PRECISION>(fission_sources,size);

  /* Deallocate memory for fission source array */
  delete [] fission_sources;

  /* Normalize scalar fluxes in each FSR */
  norm_factor = 1.0 / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Norm. factor = %f",
             tot_fission_source, norm_factor);

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) *= norm_factor;
      _old_scalar_flux(r,e) *= norm_factor;
    }
  }

  /* Normalize angular boundary fluxes for each Track */
#pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar_2; p++) {
        for (int e=0; e < _num_groups; e++) {
          _boundary_flux(t,d,p,e) *= norm_factor;
          _start_flux(t,d,p,e) *= norm_factor;
        }
      }
    }
  }
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPUSolver::computeFSRSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s, fiss_mat;
    FP_PRECISION scatter_source, fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          scatter_sources[g_prime] = sigma_s * _scalar_flux(r,g_prime);
          fission_sources[g_prime] = fiss_mat * _scalar_flux(r,g_prime);
        }

        scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources,
                                                    _num_groups);
        fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                    _num_groups);
        fission_source /= _k_eff;

        /* Compute total (scatter+fission+fixed) reduced source */
        _reduced_sources(r,g) = _fixed_sources(r,g);
        _reduced_sources(r,g) += scatter_source + fission_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] fission_sources;
    delete [] scatter_sources;
  }
}

/**
 * @brief Computes the total fission source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRFissionSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION fiss_mat;
    FP_PRECISION fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          fission_sources[g_prime] = fiss_mat * _scalar_flux(r,g_prime);
        }

        fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                    _num_groups);

        /* Compute total (fission) reduced source */
        _reduced_sources(r,g) = fission_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] fission_sources;
  }
}

/**
 * @brief Computes the total scattering source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRScatterSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s;
    FP_PRECISION scatter_source;
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          scatter_sources[g_prime] = sigma_s * _scalar_flux(r,g_prime);
        }

        scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources,
                                                    _num_groups);

        /* Compute total (scatter) reduced source */
        _reduced_sources(r,g) = scatter_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] scatter_sources;
  }
}

/**
 * @brief Computes the residual between source/flux iterations.
 * @param res_type the type of residuals to compute
 *        (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)
 * @return the average residual in each FSR
 */
double CPUSolver::computeResidual(residualType res_type) {

  int norm;
  double residual;
  double* residuals = new double[_num_FSRs];
  memset(residuals, 0., _num_FSRs * sizeof(double));

  if (res_type == SCALAR_FLUX) {

    norm = _num_FSRs;

#pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        if (_old_scalar_flux(r,e) > 0.) {
          residuals[r] += pow((_scalar_flux(r,e) - _old_scalar_flux(r,e)) /
                               _old_scalar_flux(r,e), 2);
      }
    }
  }

  else if (res_type == FISSION_SOURCE) {

    if (_num_fissionable_FSRs == 0)
      log_printf(ERROR, "The Solver is unable to compute a "
                 "FISSION_SOURCE residual without fissionable FSRs");

    norm = _num_fissionable_FSRs;

#pragma omp parallel
    {

      double new_fission_source, old_fission_source;
      FP_PRECISION* nu_sigma_f;
      Material* material;

#pragma omp for schedule(guided)
      for (int r=0; r < _num_FSRs; r++) {
        new_fission_source = 0.;
        old_fission_source = 0.;
        material = _FSR_materials[r];

        if (material->isFissionable()) {
          nu_sigma_f = material->getNuSigmaF();

          for (int e=0; e < _num_groups; e++) {
            new_fission_source += _scalar_flux(r,e) * nu_sigma_f[e];
            old_fission_source += _old_scalar_flux(r,e) * nu_sigma_f[e];
          }

          if (old_fission_source > 0.)
            residuals[r] = pow((new_fission_source -  old_fission_source) /
                               old_fission_source, 2);
        }
      }
    }
  }

  else if (res_type == TOTAL_SOURCE) {

    norm = _num_FSRs;

#pragma omp parallel
    {

      double new_total_source, old_total_source;
      FP_PRECISION inverse_k_eff = 1.0 / _k_eff;
      FP_PRECISION* nu_sigma_f;
      Material* material;
      FP_PRECISION* sigma_s;

#pragma omp for schedule(guided)
      for (int r=0; r < _num_FSRs; r++) {
        new_total_source = 0.;
        old_total_source = 0.;
        material = _FSR_materials[r];
        sigma_s = material->getSigmaS();

        if (material->isFissionable()) {
          nu_sigma_f = material->getNuSigmaF();

          for (int e=0; e < _num_groups; e++) {
            new_total_source += _scalar_flux(r,e) * nu_sigma_f[e];
            old_total_source += _old_scalar_flux(r,e) * nu_sigma_f[e];
          }

          new_total_source *= inverse_k_eff;
          old_total_source *= inverse_k_eff;
        }

        /* Compute total scattering source for group G */
        for (int G=0; G < _num_groups; G++) {
          for (int g=0; g < _num_groups; g++) {
            new_total_source += sigma_s[G*_num_groups+g]
                * _scalar_flux(r,g);
            old_total_source += sigma_s[G*_num_groups+g]
                * _old_scalar_flux(r,g);
          }
        }

        if (old_total_source > 0.)
          residuals[r] = pow((new_total_source -  old_total_source) /
                             old_total_source, 2);
      }
    }
  }

  /* Sum up the residuals from each FSR and normalize */
  residual = pairwise_sum<double>(residuals, _num_FSRs);
  residual = sqrt(residual / norm);

  /* Deallocate memory for residuals array */
  delete [] residuals;

  return residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 */
void CPUSolver::computeKeff() {

  FP_PRECISION fission;
  FP_PRECISION* FSR_rates = new FP_PRECISION[_num_FSRs];
  FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

  /* Compute the old nu-fission rates in each FSR */
#pragma omp parallel
  {

    int tid = omp_get_thread_num() * _num_groups;
    Material* material;
    FP_PRECISION* sigma;
    FP_PRECISION volume;

#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      volume = _FSR_volumes[r];
      material = _FSR_materials[r];
      sigma = material->getNuSigmaF();

      for (int e=0; e < _num_groups; e++)
        group_rates[tid+e] = sigma[e] * _scalar_flux(r,e);

      FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
      FSR_rates[r] *= volume;
    }
  }

  /* Reduce new fission rates across FSRs */
  fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  _k_eff *= fission;

  delete [] FSR_rates;
  delete [] group_rates;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void CPUSolver::transportSweep() {

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->zeroCurrents();

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

  /* Copy starting flux to current flux */
  copyBoundaryFluxes();

  /* Tracks are traversed and the MOC equations from this CPUSolver are applied
     to all Tracks and corresponding segments */
  TransportSweep sweep_tracks(_track_generator);
  sweep_tracks.setCPUSolver(this);
  sweep_tracks.execute();
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the FSR
 *          scalar flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param next_segment a pointer to the next segment to know current to tally to
 * @param azim_index azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param fsr_flux a pointer to the temporary FSR flux buffer
 * @param map_fsr_to_cell a map to find which current to tally to
 */
void CPUSolver::tallyScalarFlux(segment* curr_segment, segment* next_segment,
                                int azim_index, FP_PRECISION* track_flux,
                                FP_PRECISION* fsr_flux,
                                std::map<int, Cell*>& map_fsr_to_cell) {
  
  int fsr_id = curr_segment->_region_id;
  int next_fsr_id = next_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exponential;
  
  // Temporary storage for partial currents (OTF allocation is slow, change!)
  double* track_current = new double [_num_groups * _num_polar_2]();
  double* track_current_length = new double [_num_groups * _num_polar_2]();

  Cell* Cell_from = map_fsr_to_cell.find(fsr_id)->second;
  Cell* Cell_to = map_fsr_to_cell.find(next_fsr_id)->second;
  int cell_from = Cell_from->getId() - 10000;
  int cell_to = Cell_to->getId() - 10000;

  /* Get index for storing partial currents */
  double* ref_partial_current;
  double* jump_condition;
  int index_partial_current;
  if(cell_from != cell_to){
    ref_partial_current = getReferencePartialCurrents(cell_from, cell_to, 
                                                      &index_partial_current);
  }

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this FSR and tally current*/
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar_2; p++) {
      exponential = _exp_evaluator->computeExponential(sigma_t[e] * length, p);
      delta_psi = (track_flux(p,e)-_reduced_sources(fsr_id,e)) * exponential;
      track_flux(p,e) -= delta_psi;
      
      fsr_flux[e] += delta_psi * _quadrature->getWeightInline(azim_index, p);
      
      // MODIFY this to tally a surface angular flux rather than a current
      track_current[e*_num_polar_2 + p] += 
           double(_quadrature->getWeightInline(azim_index, p) * track_flux(p,e));
      track_current_length[e*_num_polar_2 + p] += 
           double(_quadrature->getAzimSpacing(azim_index));  		   
	  
    }
  }

  /* Atomically increment the FSR scalar flux from the temporary array */
  omp_set_lock(&_FSR_locks[fsr_id]);
  {
    for (int e=0; e < _num_groups; e++){
      _scalar_flux(fsr_id,e) += fsr_flux[e];

      /* Increment the angular partial currents */
      if (cell_from != cell_to){
        for (int p=0; p<_num_polar_2; p++){
          _ongoing_partial_currents[index_partial_current]
                  [e*_num_polar_2*_num_azim + azim_index*_num_polar_2 + p]
                  += track_current[e*_num_polar_2 + p];
		  _ongoing_partial_currents_length[index_partial_current]
                  [e*_num_polar_2*_num_azim + azim_index*_num_polar_2 + p]
                  += track_current_length[e*_num_polar_2 + p];
        }
      }
    }
  }
  omp_unset_lock(&_FSR_locks[fsr_id]);

  /* De-allocate track current array */
  delete track_current;
  delete track_current_length;
}




/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index the azimuthal index for this segmenbt
 * @param track_flux a pointer to the Track's angular flux
 * @param fwd boolean indicating direction of integration along segment
 */
void CPUSolver::tallyCurrent(segment* curr_segment, int azim_index,
                             FP_PRECISION* track_flux, bool fwd) {

  /* Tally surface currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->tallyCurrent(curr_segment, track_flux, azim_index, fwd);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective and periodic boundary conditions, the outgoing
 *          boundary flux for the Track is given to the corresponding reflecting
 *          or periodic Track. For vacuum boundary conditions, the outgoing flux
 *          is tallied as leakage.
 * @param track_id the ID number for the Track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param direction the Track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the Track's outgoing angular flux
 */
void CPUSolver::transferBoundaryFlux(int track_id,
                                     int azim_index,
                                     bool direction,
                                     FP_PRECISION* track_flux) {
  int start;
  bool transfer_flux;
  int track_out_id;

  /* For the "forward" direction */
  if (direction) {
    start = _tracks[track_id]->isNextOut() * _polar_times_groups;
    transfer_flux = _tracks[track_id]->getTransferFluxOut();
    track_out_id = _tracks[track_id]->getTrackOut()->getUid();
  }

  /* For the "reverse" direction */
  else {
    start = _tracks[track_id]->isNextIn() * _polar_times_groups;
    transfer_flux = _tracks[track_id]->getTransferFluxIn();
    track_out_id = _tracks[track_id]->getTrackIn()->getUid();
  }

  FP_PRECISION* track_out_flux = &_start_flux(track_out_id, 0, 0, start);

  /* Loop over polar angles and energy groups */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar_2; p++)
      track_out_flux(p,e) = track_flux(p,e) * transfer_flux;
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPUSolver::addSourceToScalarFlux() {

  FP_PRECISION volume;
  FP_PRECISION* sigma_t;

  /* Add in source term and normalize flux to volume for each FSR */
  /* Loop over FSRs, energy groups */
#pragma omp parallel for private(volume, sigma_t) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    volume = _FSR_volumes[r];
    sigma_t = _FSR_materials[r]->getSigmaT();

    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) /= (sigma_t[e] * volume);
      _scalar_flux(r,e) += (FOUR_PI * _reduced_sources(r,e));
    }
  }
}


/**
 * @brief Computes the volume-integrated, energy-integrated nu-fission rate in
 *        each FSR and stores them in an array indexed by FSR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          FSR nu-fission rates as a NumPy array. An example of how this
 *          method can be called from Python is as follows:
 *
 * @code
 *          num_FSRs = geometry.getNumFSRs()
 *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
 * @endcode
 *
 * @param fission_rates an array to store the nu-fission rates (implicitly
 *                      passed in as a NumPy array from Python)
 * @param num_FSRs the number of FSRs passed in from Python
 */
void CPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to compute FSR fission rates since the "
               "source distribution has not been calculated");

  log_printf(INFO, "Computing FSR fission rates...");

  FP_PRECISION* sigma_f;
  FP_PRECISION volume;

  /* Initialize fission rates to zero */
  for (int r=0; r < _num_FSRs; r++)
    fission_rates[r] = 0.0;

  /* Loop over all FSRs and compute the volume-averaged fission rate */
#pragma omp parallel for private (sigma_f, volume) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    sigma_f = _FSR_materials[r]->getSigmaF();
    volume = _FSR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      fission_rates[r] += sigma_f[e] * _scalar_flux(r,e) * volume;
  }
}
