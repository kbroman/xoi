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
 * centromeres = positions of centromeres (in microns, length n)
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
void est_coi_um(int n, double **XOLoc, int *n_xo, double *sclength,
                double *centromeres, int *group, int n_group,
                double intwindow, double coiwindow,
                double *intloc, int n_intloc,
                double *coiloc, int n_coiloc,
                double **Intensity, double *coincidence);


/* to be called from R */
void R_est_coi_um(int *n, double *xoloc, int *n_xo, double *sclength,
                  double *centromeres, int *group, int *n_group,
                  double *intwindow, double *coiwindow,
                  double *intloc, int *n_intloc,
                  double *coiloc, int *n_coiloc,
                  double *intensity, double *coincidence);

/* estimate intensity function for one group */
void est_coi_um_intensity(int n, double **AdjustedXOPos, int *n_xo,
                          double *sclength, double *centromeres,
                          int *group, int which_group,
                          double intwindow,
                          double *intloc, int n_intloc,
                          double *intensity);

/* grab the intensities that correspond to each XOLoc position */
void grab_intensities(int n, double **XOLoc, int *n_xo,
                      int *group, double *intloc, int n_intloc,
                      double **Intensity, double **IntensityVal);

/* find index of element in vec that is closest to x */
/* with ties, we just pick the first one */
int find_index_of_closest_value(double x, int n, double *vec);

/* calculate the adjusted XO positions */
/* p-arm in (0,0.5); q-arm in (0.5, 1) */
void calc_adjusted_xo_pos(int n, double **XOLoc, int *n_xo,
                          double *sclength, double *centromeres,
                          double **AdjustedXOPos);

/* estimate coincidence */
void est_coi_um_coincidence(int n, double **XOLoc, double **IntensityVals,
                            int *n_xo, double *sclength, double *centromeres,
                            double intwindow, double coiwindow, double *coiloc,
                            int n_coiloc, double *coincidence);
