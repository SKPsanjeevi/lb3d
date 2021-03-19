
typedef struct
{
      int           nx, ny, nz;
      char          datapath[160];
      char          outpath[160];
      char          inputfile[160];
      char          datafile[160];
      char          outfile[160];
      char          checkfile[160];
      unsigned long n_iteration, n_sci;
      int           ntask;   /* Single run by default, non taskfarming
			      */
      char          *folder;
      char          *gr_out_file;
      int           direction_z, direction_x, direction_y;

} Params;

typedef struct
{
      double min, max;

} Stats;


/* #define PATH         "../output/" */
#define PI           (3.14159265358979323846264338327950288)
#define DPI          (2 * PI)
#define L            (4.0)  /* Physical system legth, constant whatever
			       number of nodes params->nx in use, arbitrary units.
			       for now we use same length and no. of
			       nodes for each dimension. */ 
#define MAXLINE      (80)   /* Max length of 
			       input-file lines */
#define NAMELEN      (40)   /* Max length of
			       file of folder names (no path) */
#define PATHNAMELEN  (160)  /* Max length of
			       file of folder names (path included) */
#define VERBOSE      1
