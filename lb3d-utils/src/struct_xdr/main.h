typedef struct
{
      double         ***data;
      double         ***data2;
      double         ***data3;
} Dat;

typedef struct
{
      /* unsigned long nx, ny, nz; */
      int           nx, ny, nz;
      int           owd,arr;
      unsigned long kradius;
      char          datapath[160];
      char          outpath[160];
      char          inpfile[160];
      char          datafile[160];
      char          filetype[50];
      char          folder[160];
      char          studyname[160];
      char          field[10], outfile[240];
      unsigned long n_iteration, n_sci, sci_start;
      int           ntask;   /* Single run by default, non taskfarming
			      */
      
} Params;

typedef struct
{
      double min, max, avg;
      double min2, max2, avg2;
      double min3, max3, avg3;
} Stats;

typedef struct
{
      double         *S, *errS, *sum_sqS;
      unsigned long  *sum_shell;
      double         k1, k2, R1, R2;
      double         errk1, errk2, errR1, errR2;

} CASF;        /* Circularly Avg SF */

/* Define PATH where output data is, with trailing slash '/'
   inclusive!! 
   The program is designed to look out for PATH and then
   it will look for 'folder' 
   */

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
#define TIMECHARSTR  "8"
#define VERBOSE      0
