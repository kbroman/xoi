#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))

/* cleansing.c */
void R_cleansing(int *xovec, int *n_ind, int *n_xo, int *ob_ind, 
		 double *marker, double *criteria, int *cl);

void cleansing(int *xovec, int n_ind, int n_xo, int *ob_ind, 
	       double *marker, double criteria, int *cl);



/* coincidence.c */
void R_get_coincidence(int *xovec, double *int_dat, double *window, 
		       double *center, int *n_xo, int *n_pos, 
		       int *n_center, int *start_d, double *marker, 
		       double *coincidence);

void get_coincidence(int *xovec, double *int_dat, double window, 
		     double *center, int n_xo, int n_pos, int n_center, 
		     int start_d, double *marker, double *coincidence) ;



/* get_n_xo.c */
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno);

void R_get_N_xo(int *n_ind, int *n_pos, int *dat, int *n_xo);

int get_N_xo(int n_ind, int n_pos, int **Geno);



/* identify.c */
void R_identify_xo(int *sdat, int *n_ind, int *n_pos, int *n_xo, 
		   int *left, int *right, int *ind_id, int *ob_ind);

void identify_xo(int *sdat, int n_ind, int n_pos, int n_xo, int *left, 
		 int *right, int *ind_id, int *ob_ind);



/* intensity.c */
void R_get_intensity(int *xovec, double *window, double *center, 
		     int *n_pos, int *n_xo,  int *n_center, 
		     double *marker, double *intensity, int *n_ind);

void get_intensity(int *xovec, double window, double *center, 
		   int n_pos, int n_xo, int n_center, double *marker, 
		   double *intensity, int n_ind);
