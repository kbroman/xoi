/**
 * coi_um.h
 *
 * Estimate the coincidence as a function of micron distance, with
 * data on XO locations in microns plus SC length in microns.
 **/

/*
 * n = number of cells
 * xoloc = XO locations in each cell (in microns)
 * n_xo = number of XO in each cell (length n)
 * sclength = SC lengths (in microns, length n)
 * group = vector of groups with common intensity function, {1, ..., n_group}
 * n_group = number of groups
 * intwindow = window for smoothing intensity function
 * coiwindow = window for smoothing coincidence
 * intloc = positions at which to calculate intensity
 * n_intloc = length of intloc
 * coiloc = values at which to calculate coincidence
 * n_coiloc = length of coiloc
 * intensity = vector of length n_intloc * n_group, to contain estimated intensities
 * coincidence = vector of length n_coiloc, to contain estimated coincidence

 */
void R_est_coi_um(int *n, double *xoloc, int *n_xo, double *sclength,
                  int *group, int *n_group, double *intwindow,
                  double *coiwindow, double *intloc, int *n_intloc,
                  double *coiloc, int *n_coiloc,
                  double *intensity, double *coincidence);

void est_coi_um(int n, double **XOLoc, int *n_xo, double *sclength,
                int *group, int n_group, double intwindow,
                double coiwindow, double *intloc, int n_intloc,
                double *coiloc, int n_coiloc,
                double **Intensity, double *coincidence);

