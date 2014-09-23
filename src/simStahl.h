void simStahl(int *n_sim, double *nu, double *p, double *L,
              int *nxo, double *loc, int *max_nxo, 
              int *n_bins4start);

void simStahl_int(int n_sim, int m, double p, double L,
                  double Lstar, int *nxo, double **Loc, 
                  int max_nxo, int obligate_chiasma);

void R_simStahl_int(int *n_sim, int *m, double *p, double *L,
                    double *Lstar, int *nxo, double *loc, 
                    int *max_nxo, int *obligate_chiasma);

int random_int(int low, int high);

