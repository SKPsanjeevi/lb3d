double uncorr_error_quot(const double a, const double s_a, const double b, const double s_b);
double uncorr_error_product(const double a, const double s_a, const double b, const double s_b);

int find_extentminz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
int find_extentmaxz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
double gammadot(float ***vel, const int dx, const int j, const int k);

