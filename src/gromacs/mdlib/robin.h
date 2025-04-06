/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
 

                               &@@@@@@@@@@@@@@@@@@@(                  
                       @@@@@@,..........................@@            
                 ,@@@@......................................@         
             &@@*................*%@@@@@@(....................@       
          @@....../@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@      
       @@.&@@@@@/ @@@@@@@@@@@@......@@@@@@@@@@@@@@@............@      
    (@@        @@@@@@@@@@@..........@@@@@@@@@@@@@@@............@      
             @@@@@@@@@@@@@..........@@@@@@@@@@@@@@%...........@@@     
           .@@@@@@@@@@@@@@.........@@@@@@@@@@@@@@............@@@@@@   
          (@@@@@@@@@@@@@@..........@@@@@@@@@@@,............@@@@@@@@@  
          @@@@@@@@@@@@@@@..........@@@@@(...............@@@@@@@@@@@@  
          @@@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@& 
          @@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@@  
          @@@@@@@@@@@@@@............................#@@@@@@@@@@@@@@@  
           @@@@@@@@@@@@(...........#@@@@..............@@@@@@@@@@@@@   
             @@@@@@@@@@............@@@@@@@@@@@..........@@@@@@@@@     
               @@@@@@@@............@@@@@@@@@@@@@@@#......@@@@@@       
                  @@@@............,@@@@@@@@@@@@@@@@@@@/....@@         
                     @.........@@@@@@@@@@@@@@@@@@@@@@@@@%@...@        
                    @,....&@@@@@@@@@@@@@@@@@@@@@@@.         @*@%      
                    *(             Protocol 5.0                .@   



         ^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^

             Written by Guy "Wayyne" Dayhoff II
	           @ the University of South Florida, 2020-2023

		           <version 5.0, March 2023>

         ^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^


\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


#include <tensorflow/c/c_api.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				  prototypes
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


void rp_og_vfg(rvec * gmx_restrict,const rvec * gmx_restrict,int);
void rp_v0_vfg(rvec * gmx_restrict,const rvec * gmx_restrict,int);
void rp_v1_vfg(rvec * gmx_restrict,const rvec * gmx_restrict,int);
void rp_v2_vfg(rvec * gmx_restrict,const rvec * gmx_restrict,int);
void rp_v3_vfg(rvec * gmx_restrict,const rvec * gmx_restrict,int);

void rp_describe_injection(int,rvec * gmx_restrict,bool);


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				     macros
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


#define STEPS 0
#define FEATS 1
#define DELAY 0
#define STRIDE 1
#define PAUSE 2
#define ROUNDS 3
#define VALI_DELAY 4
#define TEST_DELAY 5
#define VALI_STRIDE 6
#define TEST_STRIDE 7
#define TRUE_STRIDE 8

#define STUPID_BIG_NUMBER 999999

#define GMX_rvec rvec *gmx_restrict

#define DBL TF_DOUBLE
#define FLT TF_FLOAT

#define RP_HOLDING_DUMP 0
#define RP_TRAINING_DUMP 1
#define RP_VALIDATION_DUMP 2
#define RP_TESTING_DUMP 3
#define RP_FINISHED_DUMP 4

#define RP_PI 3.1415926535897932384626433832795028841
#define EQ_COM_DIST_FROM_OX 0.006555

//scrubbed 3/17/2023
#define do_nothing() do{;}while(false)

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				  structures
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


 /****************************************************************************\
 * Robin Protocol: Feature Set Structure
 * ----------------------------------------------------------------------------
 *
 * name: feature set name
 * dump_wdith: # of features dumped during non-active dump
 * input_width: # of features expected by grayson
 * output_width: # of features returned by grayson
 * sample_bloat: # of additional steps to add to the historical record
 * capture_pos: toggle the inclusion of positions in the historical record
 * whole_molecule_model: toggle injection_circuit injection method
 * map: forward mapping function (historical record ==> grayson input)
 * pam: reverse mapping function (grayson output ==> historical record)
 * vfg: velocity_from_grayson function used in the injection_circuit
 * tensor_type: precision of the grayson model (i.e. double vs float)
 *
 *****************************************************************************/


 struct rp_fset{
   const char 	*name;
   int 		dump_width;					
   int 		input_width;
   int 		output_width;
   int 		sample_bloat;
   bool 	capture_pos;
   bool 	whole_molecule_model;
   void 	*(*map)(int);
   void 	(*pam)(int,const GMX_rvec,double *);
   void 	(*vfg)(GMX_rvec,const GMX_rvec, int); 
   TF_DataType 	tensor_type;
 };


 /****************************************************************************\
 * Robin Protocol Structure
 * ----------------------------------------------------------------------------
 * inpur_rec: gmx input record (i.e. mdp options etc.)
 * natoms: number of atoms
 * Ngr: number of MLMD groups (i.e. number of models)
 * grpsz: number of reference indices per group
 * inj: counter for the number of MLMD injections performed during active run
 * final_rp_step: timestep to terminate robin protocol
 * dt_per_ps: conversion factor for nm/ps ==> nm/timestep 
 * active: true means performing injections, false means dumping data
 * dump_vels: toggle to dump velocities for both active and non-active dumps
 * dump_rvels: toggle to dump velocities that are in reverse time
 * balanced_diet: try to get a uniform non-active dump
 * dilation_factor: time dilation conversion number
 * dumped: used for non-active dump record 
 * histr_shape: dimensions of the historical record
 * input_shape: dimensions of the grayson input
 * fset: the feature set index corresponding to the employed feature set
 * reporter: the atom index used to enable verbosity
 * domain: specific molecule index to apply MLMD to when active
 * pattern: delay, stride, pause, rounds, vali_delay, test_delay, true_stride
 * aternate_patterns: switch between two patterns?
 * alternative_stride: if alternate_patterns, this is the alternative stride
 * alternative_pause: if alternate_patterns, this is the alterantive pause
 * pattern2: pattern placeholders...but only used stride and pause
 * threshold_type: type of threshold to apply (i.e. magnitude or component)
 * threshold: the threshold
 * map: forward map function (historical record ==> grayson input)
 * pam: reverse map function (grayson output ==> historical record)
 * rec: recording function used to construct the historical record
 * history: the historical record
 * forecast: storage for grayson predictions
 * rotation: rotation record for feature sets that rotate velocities
 * t0_magnitude: the t0 magnitude used to normalize features
 * ml_tog: per-atom MLMD toggle
 * models: the paths to grayson models
 * groups: the dominion of each grayson model
 * modulus: number of atoms per molecule (e.g. H2O = 3, O2 = 2)
 * sesh: tensorflow sessions (1 per model)
 * graph: tensorflow graphs (1 per model)
 * status: tensorflow status (1 per model)
 * io_opts: tensorflow i/o options (1 per model)
 * grayson_ndims: number of grayson input dimensions
 * grayson_inputsz: size of grayson input in bytes
 * grayson_dims: grayson input dimensions (i.e. batchs,steps,features)
 * dump_state: the status of a nonactive dump 
 * box: the cell dims for md
 * whole_molecule_model: true for whole molecule descriptors
 * inv_mass: pointer to the inv_mass array that gmx maintains
 * vfg: the velocity from grayson fn to use based on fset
 *
 *****************************************************************************/


typedef struct {
  t_inputrec *input_rec;
  int natoms;
  int Ngr;
  int non_active_strech;
  int *grpsz;
  int inj;
  int final_rp_step;
  double dt_per_ps;
  bool active;
  bool dump_vels;
  bool dump_rvels;
  bool balanced_diet;
  int dilation_factor;
  int **dumped;
  int histr_shape[2];
  int input_shape[2];
  int fset;
  int reporter;
  int domain;
  int pattern[8];
  bool alternate_patterns;
  int alternative_stride;
  int alternative_pause;
  int pattern2[8];
  int threshold_type;
  double threshold;
  int threshold_domain;
  void *(*map)(int);
  void (*pam)(int,const GMX_rvec,double *);
  void (*rec)(int,int,const GMX_rvec,GMX_rvec);
  real ***history;
  double **forecast;
  double **rotation;
  double **t0_magnitude;
  int *ml_tog;
  const char **models;
  int **groups;
  int modulus;
  TF_Session **sesh;
  TF_Graph **graph;
  TF_Status **status;
  TF_Output **io_opts;
  int grayson_ndims;
  int grayson_inputsz;
  int64_t grayson_dims[3];
  int dump_state;
  real box[3];
  bool whole_molecule_model;
  const real * gmx_restrict inv_mass;
  void (*vfg)(rvec * gmx_restrict,const rvec * gmx_restrict,int);
} robin_protocol;


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				   globals
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


robin_protocol *rp = NULL; 

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				  ancillary
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/17/2023
double *cross(double *vect_A, double *vect_B){
    double *cp = (double *)malloc(sizeof*cp*3);

    cp[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cp[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cp[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

   return cp;
}


//scrubbed 3/17/2023
double dot(double *v, double *u, int n){
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}


//scrubbed 3/17/2023
double mag(double *v, int n){
    double result = 0.0;
    for (int i = 0; i < n; i++)
      result += v[i]*v[i];
    return sqrt(result);
}


//scrubbed 3/17/2023
real dx_periodic(real x1, real x2, real xbox){
    real dx;
    dx=x2-x1;
    if(fabs(dx)>xbox/2.0)
        dx=xbox -(real)fabs(dx);
    return(dx);
}


//scrubbed 3/17/2023
int rp_get_group(int n){

    int mod = n%rp->modulus;

    for(int i=0;i<rp->Ngr;i++)
      for(int j=0;j<rp->grpsz[i];j++)
        if (rp->groups[i][j] == mod)
		return i;

    printf("RP_WARNING: Unable to get group for atom %d. Goodbye.\n",n);
    exit(1);
}


//scrubbed 3/16/2023
void rp_get_eqcom_from_gmx(int n, const rvec * gmx_restrict x, double com[3]){

    double h1[3], h2[3], m1, m2;

    h1[XX] = x[n+1][XX]-x[n][XX];
    h1[YY] = x[n+1][YY]-x[n][YY];
    h1[ZZ] = x[n+1][ZZ]-x[n][ZZ];

    h2[XX] = x[n+2][XX]-x[n][XX];
    h2[YY] = x[n+2][YY]-x[n][YY];
    h2[ZZ] = x[n+2][ZZ]-x[n][ZZ];

    m1 = mag(h1,3);
    m2 = mag(h2,3);

    h1[XX] /= m1;
    h1[YY] /= m1;
    h1[ZZ] /= m1;

    h2[XX] /= m2;
    h2[YY] /= m2;
    h2[ZZ] /= m2;

    h1[XX] += x[n][XX];
    h1[YY] += x[n][YY];
    h1[ZZ] += x[n][ZZ];

    h2[XX] += x[n][XX];
    h2[YY] += x[n][YY];
    h2[ZZ] += x[n][ZZ];

    com[XX] = ((x[n][XX]*15.999)+(h1[XX]*1.0078)+(h2[XX]*1.0078))/18.0146;
    com[YY] = ((x[n][YY]*15.999)+(h1[YY]*1.0078)+(h2[YY]*1.0078))/18.0146;
    com[ZZ] = ((x[n][ZZ]*15.999)+(h1[ZZ]*1.0078)+(h2[ZZ]*1.0078))/18.0146;

    com[XX] -= x[n][XX];
    com[YY] -= x[n][YY];
    com[ZZ] -= x[n][ZZ];

    m1 = mag(com,3);

    com[XX] /= m1;
    com[YY] /= m1;
    com[ZZ] /= m1;

    com[XX] *= EQ_COM_DIST_FROM_OX;
    com[YY] *= EQ_COM_DIST_FROM_OX;
    com[ZZ] *= EQ_COM_DIST_FROM_OX;

    com[XX] += x[n][XX];
    com[YY] += x[n][YY];
    com[ZZ] += x[n][ZZ];
}


//deprecated 3/17/2023
bool threshold_satisfied(int n){

    if (!strcmp(rp_threshold_types[rp->threshold_type],"magnitude")){
      double mag = rp->history[n][17][0]*rp->history[n][17][0];
      mag += rp->history[n][17][1]*rp->history[n][17][1];
      mag += rp->history[n][17][2]*rp->history[n][17][2];
     
      if (sqrt(mag) < rp->threshold){
        printf("skipping ML on %d due to |V| threshold\n",n);
        rp->ml_tog[n]--;
        return FALSE;
      }
    }
    else if (!strcmp(rp_threshold_types[rp->threshold_type],"xyz")){
      if (fabs(rp->history[n][17][0]-rp->history[n][16][0])<rp->threshold){
        printf("skipping ML on %d due to dVx threshold\n",n);
        return FALSE;
      }
      if (fabs(rp->history[n][17][1]-rp->history[n][16][1])<rp->threshold){
        printf("skipping ML on %d due to dVy threshold\n",n);
        return FALSE;
      }
      if (fabs(rp->history[n][17][2]-rp->history[n][16][2])<rp->threshold){
        printf("skipping ML on %d due to dVz threshold\n",n);
        return FALSE;
      }
    }
    return TRUE;
}

void *rp_normalized_map(double *map, int atom_index){

        int steps=rp->input_shape[STEPS];
        int feats=rp->input_shape[FEATS];
        
        float mu, sigma;

        //for each feature, compute the average (mu) and stddev (sigma)
        for (int f=0;f<feats;f++){
          mu = sigma = 0.f;
          for (int i=f;i<steps*feats;i+=feats)
                //track the sum to compute mu
                mu += map[i];

          //compute the average
          rp->rotation[atom_index][(f*2)+4] = mu /= steps;
            
          //compute the standard deviation
          for (int i=f;i<steps*feats;i+=feats)
                sigma += pow(map[i]-mu,2);
          sigma /= steps;
          rp->rotation[atom_index][(f*2)+3] = sigma = sqrt(sigma);

          //use mu and sigma to normalization & standardize the input
          //s.t. the mean is 0 and the stddev is 1
          for (int i=f;i<steps*feats;i+=feats)
                map[i] = (map[i] - mu) / sigma;
        }

        return (void *)map;
}

//scrubbed 3/17/2023
void *rp_shifted_map(double *map, int atom_index){

	int steps=rp->active ? rp->input_shape[STEPS] : rp->input_shape[STEPS] - 1;
	int feats=rp->input_shape[FEATS];
	
	double shift;

	for (int f=0;f<feats;f++){
	  for (int i=f;i<steps*feats;i+=feats)
		shift = map[i] = rp->t0_magnitude[atom_index][f+1];

	  for (int i=f;i<steps*feats;i+=feats)
		map[i] -= shift;
	}

	return (void *)map;
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				  rotations
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/17/2023
//computes counter-clockwise euler angle (postive theta) 
//has a range of 0 to 359.99999 with 360->0
double euler_angle(double opp, double adj){
  double e;

  if (opp >= 0 && adj >= 0)
      e =  fabs(atan(opp/adj));
  else if (opp >= 0 && adj < 0)
      e =  RP_PI-fabs(atan(opp/adj));
  else if (opp < 0 && adj >= 0)
      e = (2.0*RP_PI)-fabs(atan(opp/adj));
  else if (opp < 0 && adj < 0)
      e =  RP_PI+fabs(atan(opp/adj)); 

  return fmod(e,2.0*RP_PI);
}


//scrubbed 3/17/2023
//aka nutation
void bank(double *v, double theta){
     double xprime, yprime, zprime;

     xprime  = *(v+0);

     yprime  = *(v+1) * cos(theta);
     yprime += *(v+2) * -sin(theta);

     zprime  = *(v+1) * sin(theta);
     zprime += *(v+2) * cos(theta);

     //TODO: is the truncation hurting the models?
//    *(v+0) = fabs(xprime) < 1e-9 ? 0.0 : xprime;
//    *(v+1) = fabs(yprime) < 1e-9 ? 0.0 : yprime;
//    *(v+2) = fabs(zprime) < 1e-9 ? 0.0 : zprime;
      *(v+0) =  xprime;
      *(v+1) =  yprime;
      *(v+2) =  zprime;
}

//scrubbed 3/17/2023
void pitch(double *v, double alpha){
     double xprime, yprime, zprime;

     xprime  = *(v+0) * cos(alpha);
     xprime += *(v+2) * sin(alpha);
                     
     yprime  = *(v+1);
                     
     zprime  = *(v+0) *-sin(alpha);
     zprime += *(v+2) * cos(alpha);

     //TODO: is the truncation hurting the models?
//    *(v+0) = fabs(xprime) < 1e-9 ? 0.0 : xprime;
//    *(v+1) = fabs(yprime) < 1e-9 ? 0.0 : yprime;
//    *(v+2) = fabs(zprime) < 1e-9 ? 0.0 : zprime;
      *(v+0) =  xprime;
      *(v+1) =  yprime;
      *(v+2) =  zprime;

}


//scrubbed 3/17/2023
//aka precession
void yaw(double *v, double phi){
     double xprime, yprime, zprime;

     xprime  = *(v+0) * cos(phi);
     xprime += *(v+1) *-sin(phi);
                     
     yprime  = *(v+0) * sin(phi);
     yprime += *(v+1) * cos(phi);
                     
     zprime  = *(v+2);

     //TODO: is the truncation hurting the models?
//    *(v+0) = fabs(xprime) < 1e-9 ? 0.0 : xprime;
//    *(v+1) = fabs(yprime) < 1e-9 ? 0.0 : yprime;
//    *(v+2) = fabs(zprime) < 1e-9 ? 0.0 : zprime;
     *(v+0) =  xprime;
     *(v+1) =  yprime;
     *(v+2) =  zprime;
}

//scrubbed 3/17/2023
//rotate [counter-clockwise] the 3D vector a [inplace]
//about the 3D vector b (must be a unit vector) by theta [positive] radians
//negative theta will yield clockwise rotation
//via the orthogonal compoenent method (Rodrigues rotation)
void rodrot(double *a, double *b, double theta){
   double a_parr[3], a_perp[3], x0, x1, x2, *w_norm, b_hat[3];
 
   double mb = mag(b,3);
   b_hat[0] = b[0]/mb;
   b_hat[1] = b[1]/mb;
   b_hat[2] = b[2]/mb;

   x0 = dot(a,b_hat,3)/dot(b_hat,b_hat,3);

   a_parr[0] = x0*b_hat[0];
   a_parr[1] = x0*b_hat[1];
   a_parr[2] = x0*b_hat[2];

   a_perp[0] = a[0]-a_parr[0];
   a_perp[1] = a[1]-a_parr[1];
   a_perp[2] = a[2]-a_parr[2];

   w_norm = cross(b_hat,a_perp);

   x1 = cos(theta)/mag(a_perp,3);
   x2 = sin(theta)/mag(w_norm,3);

   a[0] = a_parr[0] + mag(a_perp,3)*(x1*a_perp[0] + x2*w_norm[0]);
   a[1] = a_parr[1] + mag(a_perp,3)*(x1*a_perp[1] + x2*w_norm[1]);
   a[2] = a_parr[2] + mag(a_perp,3)*(x1*a_perp[2] + x2*w_norm[2]);

   for (int i=0;i<3;i++)
     if (fabs(a[i]) < 1e-9)
       a[i] = 0.0;

   //housekeeping
   free(w_norm);
}


//scrubbed 3/17/2023
double *xz_angles(double *a, bool verbose){

  double *eas = (double*)malloc(sizeof*eas*2);

  double *A = (double*)malloc(sizeof*A*3);

  for (int i=0;i<3;i++){
	A[i] = a[i];
  }

  eas[0] = -euler_angle(A[2],A[1]);
  bank(A,eas[0]);

  //rotate ibx about z onto the x axis 
  eas[1] = -euler_angle(A[1],A[0]);
  yaw(A,eas[1]);

  if (verbose){
    printf("EXT: (% lf,% lf,% lf) = %lf\n",a[0],a[1],a[2],mag(A,3));
    printf("INT: (% lf,% lf,% lf) = %lf\n",A[0],A[1],A[2],mag(A,3));
    printf("ANG: (% lf,% lf,% lf)\n",eas[0],eas[1]);
  }

  //housekeeping
  free(A);

  return eas;
}


//scrubbed 3/17/2023
double *xzx_angles(double *a, double *v, bool verbose, int record_magnitudes){

  double *eas = (double*)malloc(sizeof*eas*3);

  double *A = (double*)malloc(sizeof*A*3);
  double *V = (double*)malloc(sizeof*V*3);

  for (int i=0;i<3;i++){
	A[i] = a[i];
	V[i] = v[i];
  }

  eas[0] = -euler_angle(A[2],A[1]);//+0.785;
  bank(A,eas[0]);
  bank(V,eas[0]);

  //rotate ibx about z onto the x axis 
  eas[1] = -euler_angle(A[1],A[0]);//+0.785;
  yaw(A,eas[1]);
  yaw(V,eas[1]);

  //rotate iby about x (now ibx) into the x-y plane
  eas[2] = -euler_angle(V[2],V[1]);//+0.785;
  bank(A,eas[2]);
  bank(V,eas[2]);

  if (record_magnitudes != -1){
    rp->t0_magnitude[record_magnitudes][0] = fabs(V[0]);
    rp->t0_magnitude[record_magnitudes][1] = fabs(V[1]);
    printf("recorded: %f and %f for atom %d\n",
	rp->t0_magnitude[record_magnitudes][0],
	rp->t0_magnitude[record_magnitudes][1],
	record_magnitudes);
    if (V[2] != 0.f)
	printf("Potential bug detected in xzx call, the Z component of the velocity is non-zero after rotations (%e)\n",V[2]);
  }

  if (verbose){
    printf("AOR: (% lf,% lf,% lf) = %lf\n",a[0],a[1],a[2],mag(A,3));
    printf("ANG: (% lf,% lf,% lf)\n",eas[0],eas[1],eas[2]);
    printf("Vst: (% lf,% lf,% lf)\n", v[0],v[1],v[2]);
    printf("Ast: (% lf,% lf,% lf)\n", a[0],a[1],a[2]);
    printf("Vfn: (% lf,% lf,% lf)\n", V[0],V[1],V[2]);
    printf("Afn: (% lf,% lf,% lf)\n", A[0],A[1],A[2]);

    double *A2 = (double*)malloc(sizeof*A*3);
    double *V2 = (double*)malloc(sizeof*V*3);

    for (int i=0;i<3;i++){
	A2[i] = a[i];
	V2[i] = v[i];
    }

    bank(A2,eas[0]);                               
     yaw(A2,eas[1]);                               
    bank(A2,eas[2]);                               

    bank(V2,eas[0]);                               
     yaw(V2,eas[1]);                               
    bank(V2,eas[2]);                               

    printf("Vsc: (% lf,% lf,% lf)\n", V2[0],V2[1],V2[2]);
    printf("Asc: (% lf,% lf,% lf)\n", A2[0],A2[1],A2[2]);
   
    //housekeeping
    free(A2);
    free(V2);
  }

  //housekeeping
  free(A);
  free(V);

  return eas;
}


//scrubbed 3/17/2023
double xzx_delta(double *a, double *v, double *prev_angles, bool verbose){

	double *A = (double*)malloc(sizeof*A*3);
	double *V = (double*)malloc(sizeof*V*3);
	
	for (int i=0;i<3;i++){
	      A[i] = a[i];
	      V[i] = v[i];
	}

	bank(A,prev_angles[0]); 
	 yaw(A,prev_angles[1]); 
	bank(A,prev_angles[2]); 
   
	if (verbose){ 
          printf("aor: (% lf,% lf,% lf) = %lf\n",a[0],a[1],a[2],mag(a,3));
          printf("xax: (% lf,% lf,% lf) = %lf\n\n",A[0],A[1],A[2],mag(A,3));
        }

	bank(V,prev_angles[0]); 
	 yaw(V,prev_angles[1]); 
	bank(V,prev_angles[2]); 

        double delta = atan2(V[2],V[1]);

	//housekeeping
	free(A);
	free(V);

        return delta;
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
		                neighbor searches
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/17/2023
float **rp_1d_neighbor_search(int atom_index, int vacancies,
	 const rvec * gmx_restrict x, const matrix box){

    float px,py,pz,d2;
    float **tenants = (float **)malloc(sizeof*tenants*vacancies);

    //initialize the [sorted] tenant list (index scales with distance)
    for (int v=0;v<vacancies;v++){
      tenants[v] = (float *)malloc(sizeof**tenants*2);
      tenants[v][0] = STUPID_BIG_NUMBER-v;
      tenants[v][1] = -v;
    }

    //loop over all [other] atoms in the system
    for (int i=0;i<rp->natoms;i++){
      if (i==atom_index) continue;

      //get the minimume image distances 
      px = dx_periodic(x[atom_index][XX],x[i][XX],box[XX][XX]);
      py = dx_periodic(x[atom_index][YY],x[i][YY],box[YY][YY]);
      pz = dx_periodic(x[atom_index][ZZ],x[i][ZZ],box[ZZ][ZZ]);

      //compute the squared distance
      d2 = ((px*px)+(py*py)+(pz*pz));

      //try to insert the tenant into the sorted list of tenants
      for (int v=0;v<vacancies;v++)
        if (d2 <= tenants[v][0]){
          for (int j=vacancies-1;j>v;j--){
            tenants[j][0] = tenants[j-1][0];
            tenants[j][1] = tenants[j-1][1];
          }
          tenants[v][0] = d2;
          tenants[v][1] = i;
          break;
        }
    }

    return tenants;
}


//scrubbed 3/17/2023
float **rp_3d_neighbor_search(int atom_index, int vacancies,
	const rvec * gmx_restrict x, const matrix box){

    //initialize tenant list
    float **tenants = (float**)malloc(sizeof*tenants*vacancies);

    for (int i=0;i<vacancies;i++){
      tenants[i] = (float*)malloc(sizeof**tenants*3);
      tenants[i][0] = STUPID_BIG_NUMBER;
      tenants[i][1] = STUPID_BIG_NUMBER;
      tenants[i][2] = STUPID_BIG_NUMBER;
    }

    //loop over all [other] atoms in the system
    for (int i=0;i<rp->natoms;i++){
      if (i==atom_index) continue;

      //get the minimum distance image
      real px = dx_periodic(x[atom_index][XX],x[i][XX],box[XX][XX]);
      real py = dx_periodic(x[atom_index][YY],x[i][YY],box[YY][YY]);
      real pz = dx_periodic(x[atom_index][ZZ],x[i][ZZ],box[ZZ][ZZ]);

      //compute the squared distance
      real d2 = ((px*px)+(py*py)+(pz*pz));

      //try to insert a new tenant
      for (int i=0;i<vacancies;i++){
        if (d2 < (tenants[i][0]*tenants[i][0]) +
		(tenants[i][1]*tenants[i][1]) + (tenants[i][2]*tenants[i][2])){
          for (int j=i;j<vacancies-1;j++){
            tenants[j+1][0] = tenants[j][0];
            tenants[j+1][1] = tenants[j][1];
            tenants[j+1][2] = tenants[j][2];
          }
          tenants[i][0] = px;
          tenants[i][1] = py;
          tenants[i][2] = pz;
          break;
        }
      }
    } 

    return tenants;
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
		       historical record construction
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/17/2023
void rp_rec3v(int n,int i,const GMX_rvec x, GMX_rvec v){

    //capture the velocities
    for (int k=0;k<3;k++)
      rp->history[n][i][k] = v[n][k];
}


//scrubbed 3/17/2023
void rp_make_molecules_whole(int n, real ***hr, bool verbose){

      //check i=0 positions for wrapping and correct p.r.n.
      int mod = n%rp->modulus;
      int hp = (mod == 0) ? n+1: n;
      int hs = (mod == 0) ? n+2: ((mod == 1) ? n+1 : n-1);
      int ox = (mod == 0) ?   n: ((mod == 1) ? n-1 : n-2);

      //make sure the molecule is not wrapped 
      double h1_bond[3], h2_bond[3], hh_bond[3], m1, m2, mh; 
     
      for (int k=0;k<3;k++){
        h1_bond[k] = hr[hp][0][k+3]-hr[ox][0][k+3];
        h2_bond[k] = hr[hs][0][k+3]-hr[ox][0][k+3];
        hh_bond[k] = hr[hp][0][k+3]-hr[hs][0][k+3];
      }
  
      m1 = mag(h1_bond,3);
      m2 = mag(h2_bond,3);
      mh = mag(hh_bond,3);
     
      //check if the oxygen is wrapped(m1 & m2 are wrapped but mh is not)
      if ((m1 > rp->box[XX]/2.0)
      &&  (m2 > rp->box[XX]/2.0)
      &&  (mh < rp->box[XX]/2.0)){

	if (verbose)
          printf("WRAPPED OX DETECTED: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);

	for (int k=0;k<3;k++)
	   if (hr[hp][0][k+3]-hr[ox][0][k+3] > rp->box[XX]/2.0)
	    hr[ox][0][k+3] = hr[ox][0][k+3] + rp->box[XX];
	   else if (hr[hp][0][k+3]-hr[ox][0][k+3] < -rp->box[XX]/2.0)
	    hr[ox][0][k+3] = hr[ox][0][k+3] - rp->box[XX];

 
        for (int k=0;k<3;k++){
          h1_bond[k] = hr[hp][0][k+3]-hr[ox][0][k+3];
          h2_bond[k] = hr[hs][0][k+3]-hr[ox][0][k+3];
          hh_bond[k] = hr[hp][0][k+3]-hr[hs][0][k+3];
        }
    
        m1 = mag(h1_bond,3);
        m2 = mag(h2_bond,3);
        mh = mag(hh_bond,3);

	if (verbose)
          printf("AFTER OX UNWRAPPING: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);
      }
      //check if h1 is wrapped
      else if ((m1 > rp->box[XX]/2.0)
           &&  (m2 < rp->box[XX]/2.0)
           &&  (mh > rp->box[XX]/2.0)){

	if (verbose)
          printf("WRAPPED H1 DETECTED: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);

	for (int k=0;k<3;k++)
	   if (hr[ox][0][k+3]-hr[hp][0][k+3] > rp->box[XX]/2.0)
	    hr[hp][0][k+3] = hr[hp][0][k+3] + rp->box[XX];
	   else if (hr[ox][0][k+3]-hr[hp][0][k+3] < -rp->box[XX]/2.0)
	    hr[hp][0][k+3] = hr[hp][0][k+3] - rp->box[XX];

 
        for (int k=0;k<3;k++){
          h1_bond[k] = hr[hp][0][k+3]-hr[ox][0][k+3];
          h2_bond[k] = hr[hs][0][k+3]-hr[ox][0][k+3];
          hh_bond[k] = hr[hp][0][k+3]-hr[hs][0][k+3];
        }
    
        m1 = mag(h1_bond,3);
        m2 = mag(h2_bond,3);
        mh = mag(hh_bond,3);

	if (verbose)
          printf("AFTER H1 UNWRAPPING: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);
      }
      //check if h2 is wrapped
      else if ((m1 < rp->box[XX]/2.0)
           &&  (m2 > rp->box[XX]/2.0)
           &&  (mh > rp->box[XX]/2.0)){


	if (verbose)
          printf("WRAPPED H2 DETECTED: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);

	for (int k=0;k<3;k++)
	   if (hr[ox][0][k+3]-hr[hs][0][k+3] > rp->box[XX]/2.0)
	    hr[hs][0][k+3] = hr[hs][0][k+3] + rp->box[XX];
	   else if (hr[ox][0][k+3]-hr[hs][0][k+3] < -rp->box[XX]/2.0)
	    hr[hs][0][k+3] = hr[hs][0][k+3] - rp->box[XX];

 
        for (int k=0;k<3;k++){
          h1_bond[k] = hr[hp][0][k+3]-hr[ox][0][k+3];
          h2_bond[k] = hr[hs][0][k+3]-hr[ox][0][k+3];
          hh_bond[k] = hr[hp][0][k+3]-hr[hs][0][k+3];
        }
    
        m1 = mag(h1_bond,3);
        m2 = mag(h2_bond,3);
        mh = mag(hh_bond,3);

	if (verbose)
          printf("AFTER H2 UNWRAPPING: h1: %f h2: %f hh: %f\n", 
	    m1, m2, mh);
      }
}



//scrubbed 3/17/2023
void rp_rec3v3r(int n,int i,const GMX_rvec x, GMX_rvec v){

    real ***hr = rp->history;

    //capture the velocities
    rp_rec3v(n,i,x,v);

    //capture the positions
    if (i == 0){

      if (n == rp->reporter)
        printf("\n\rRP: initializing historical record\n");

      for (int k=3;k<6;k++)
        rp->history[n][i][k] = x[n][k-3]; 

      return;
    }
    
    if (i == 1){
      if (n == rp->reporter)
        printf("RP: making sure molecules are whole\n");

      if (rp->whole_molecule_model)
        rp_make_molecules_whole(n,hr,false);
    }

    //use the new velocity to compute new position
    for (int k=3;k<6;k++)
      rp->history[n][i][k] = rp->history[n][i-1][k] 
			   + (v[n][k-3]/rp->dt_per_ps); 

 //   printf("hr[%d][%d] = %f %f %f\n",n,i,
//	rp->history[n][i][3],
//	rp->history[n][i][4],
//	rp->history[n][i][5]);
}


//scrubbed 5-20-2022 -Wayyne
void rp_round_robin(int n, const rvec * gmx_restrict x, rvec * gmx_restrict v, 
	gmx_int64_t step){

    //make sure we need a historical record
    if (rp->active && step > rp->final_rp_step) 
      return;

    //determine the working historical record index
    int i = step - (rp->pattern[DELAY] - rp->histr_shape[STEPS]);

        //shift historical record p.r.n. 
    if (i > rp->histr_shape[STEPS] - 1){
      i = rp->histr_shape[STEPS] - 1;
      for (int k=0;k<(rp->histr_shape[STEPS]-1);k++)
        for (int d=0;d<rp->histr_shape[FEATS];d++)
          rp->history[n][k][d] = rp->history[n][k+1][d];
    }
    //record the current state
    (*rp->rec)(n,i,x,v); 

    //bail out now if this is an inactive run
    if (!rp->active)
      return;
 
    //bail out now if we're not making prediction for this atom's group
    if (rp->models[rp_get_group(n)] == NULL && rp->active)
      return;

    //bail out now if this atom is outside of the specified domain
    if (rp->domain != -1 && floor(n/rp->modulus)+1.0 != rp->domain)
      return;

    //otherwise, determine the next effective step
    int nxt_estep = step + 1 - rp->pattern[DELAY];

    //and check if we should toggle on ML for the next effective step
    if (!rp->alternate_patterns || rp->alternative_stride == -1){

      double round = floor(nxt_estep/(rp->pattern[PAUSE] + rp->pattern[STRIDE]));

      if (rp->ml_tog[n] == 0 && nxt_estep >= 0 && rp->pattern[ROUNDS] > round)
        if (nxt_estep % (rp->pattern[PAUSE] + rp->pattern[STRIDE]) == 0)
          rp->ml_tog[n] = rp->pattern[STRIDE];

    }
    else{

      //a round consists of executing both patterns once
      double round = floor(nxt_estep/(rp->pattern[PAUSE] + rp->pattern[STRIDE] 
			+ rp->pattern2[PAUSE] + rp->pattern2[STRIDE]));

      //compute position in the round
      int pos = nxt_estep%(rp->pattern[PAUSE] + rp->pattern[STRIDE] 
		+ rp->pattern2[PAUSE] + rp->pattern2[STRIDE]);

      if (rp->ml_tog[n] == 0 && nxt_estep >= 0 && rp->pattern[ROUNDS] > round)
        if (pos == 0){
          rp->ml_tog[n] = rp->pattern[STRIDE];
//	  if (n == rp->reporter) printf("p1 active: effevtive step = %d\n",nxt_estep);
	}
        else if ( pos == rp->pattern[PAUSE] + rp->pattern[STRIDE]){
          rp->ml_tog[n] = rp->pattern2[STRIDE];
//	  if (n == rp->reporter) printf("p2 active: effevtive step = %d\n",nxt_estep);
	}
    }

}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
		              mapping macros
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/17/2023 - this isn't REALLY the equilibrium COM...
#define FIND_EQCOM(_X_,_Y_)						\
  tmp_b1[0] = hr[hp][_Y_][3]-hr[ox][_Y_][3];				\
  tmp_b1[1] = hr[hp][_Y_][4]-hr[ox][_Y_][4];				\
  tmp_b1[2] = hr[hp][_Y_][5]-hr[ox][_Y_][5];				\
  tmp_b2[0] = hr[hs][_Y_][3]-hr[ox][_Y_][3];				\
  tmp_b2[1] = hr[hs][_Y_][4]-hr[ox][_Y_][4];				\
  tmp_b2[2] = hr[hs][_Y_][5]-hr[ox][_Y_][5];				\
  mb1 = mag(tmp_b1,3);							\
  mb2 = mag(tmp_b2,3);							\
  for (int d=0;d<3;d++){						\
	tmp_b1[d] /= mb1;						\
	tmp_b2[d] /= mb2;						\
	tmp_b1[d] += hr[ox][_Y_][d+3];					\
	tmp_b2[d] += hr[ox][_Y_][d+3];					\
  }									\
  _X_[0] =  (double)(hr[ox][_Y_][3])*15.999;				\
  _X_[1] =  (double)(hr[ox][_Y_][4])*15.999;				\
  _X_[2] =  (double)(hr[ox][_Y_][5])*15.999;				\
  _X_[0] += tmp_b1[0]*1.0078;						\
  _X_[1] += tmp_b1[1]*1.0078;						\
  _X_[2] += tmp_b1[2]*1.0078;						\
  _X_[0] += tmp_b2[0]*1.0078;						\
  _X_[1] += tmp_b2[1]*1.0078;						\
  _X_[2] += tmp_b2[2]*1.0078;						\
  _X_[0] /= (15.999+1.0078+1.0078);					\
  _X_[1] /= (15.999+1.0078+1.0078);					\
  _X_[2] /= (15.999+1.0078+1.0078);					\


//scrubbed 3/17/2023 - this IS really the equilibirum COM 
#define FIND_EQCOM_v2(_X_,_Y_)						\
  FIND_EQCOM(_X_,_Y_)							\
  _X_[0] -= hr[ox][_Y_][3];						\
  _X_[1] -= hr[ox][_Y_][4];						\
  _X_[2] -= hr[ox][_Y_][5];						\
  mb1 = mag(_X_,3);							\
  _X_[0] = (_X_[0]/mb1)*EQ_COM_DIST_FROM_OX;				\
  _X_[1] = (_X_[1]/mb1)*EQ_COM_DIST_FROM_OX;				\
  _X_[2] = (_X_[1]/mb1)*EQ_COM_DIST_FROM_OX;				\
  _X_[0] += hr[ox][_Y_][3];						\
  _X_[1] += hr[ox][_Y_][4];						\
  _X_[2] += hr[ox][_Y_][5];						\


//scrubbed 3/17/2023 
#define FIND_O2COM(_X_,_Y_)\
  _X_[0] =  (double)(hr[0][_Y_][3]);					\
  _X_[1] =  (double)(hr[0][_Y_][4]);					\
  _X_[2] =  (double)(hr[0][_Y_][5]);					\
  _X_[0] += (double)(hr[1][_Y_][3]);					\
  _X_[1] += (double)(hr[1][_Y_][4]);					\
  _X_[2] += (double)(hr[1][_Y_][5]);					\
  _X_[0] /= 2.0;							\
  _X_[1] /= 2.0;							\
  _X_[2] /= 2.0;


//scrubbed 3/17/2023 
#define FIND_COM(_X_,_TIME_){						\
  int _start=pa-(pa%rp->modulus);					\
  int _stop=_start+rp->modulus;						\
									\
  float _total_mass = 1.f/rp->inv_mass[_start];				\
									\
  _X_[0] =  (double)(hr[_start][_TIME_][3])*1.f/rp->inv_mass[_start];	\
  _X_[1] =  (double)(hr[_start][_TIME_][4])*1.f/rp->inv_mass[_start];	\
  _X_[2] =  (double)(hr[_start][_TIME_][5])*1.f/rp->inv_mass[_start];	\
									\
  for (int I=_start+1;I<_stop;I++){					\
    _total_mass += 1.f/rp->inv_mass[I];					\
    _X_[0] +=  (double)(hr[I][_TIME_][3])*1.f/rp->inv_mass[I];		\
    _X_[1] +=  (double)(hr[I][_TIME_][4])*1.f/rp->inv_mass[I];		\
    _X_[2] +=  (double)(hr[I][_TIME_][5])*1.f/rp->inv_mass[I];		\
  }									\
									\
  _X_[0] /= _total_mass;						\
  _X_[1] /= _total_mass;						\
  _X_[2] /= _total_mass;}						\


//scrubbed 3/17/2023
#define PRIME_UELOCITIES_FOR(_ATOMINDEX_) 				\
    arr1[0] = (double)hr[_ATOMINDEX_][0][0];				\
    arr1[1] = (double)hr[_ATOMINDEX_][0][1];				\
    arr1[2] = (double)hr[_ATOMINDEX_][0][2];				\
    rp->t0_magnitude[_ATOMINDEX_][0] = mag(arr1,3);			\


//scrubbed 3/17/2023
#define RENDER_UELOCITIES_FOR(_X_,_Y_,_Z_)                              \
      mp[j+_Y_+0] = (double)(hr[_X_][_Z_][0]/rp->t0_magnitude[_X_][0]); \
      mp[j+_Y_+1] = (double)(hr[_X_][_Z_][1]/rp->t0_magnitude[_X_][0]); \
      mp[j+_Y_+2] = (double)(hr[_X_][_Z_][2]/rp->t0_magnitude[_X_][0]);	


//scrubbed 3/17/2023
#define RENDER_VELOCITIES_FOR(_X_,_Y_,_Z_)				\
      mp[j+_Y_+0] = (double)hr[_X_][_Z_][0];				\
      mp[j+_Y_+1] = (double)hr[_X_][_Z_][1];				\
      mp[j+_Y_+2] = (double)hr[_X_][_Z_][2];


//scrubbed 3/17/2023
#define BEGIN_MAP(_FSETNAME_)						\
  void *rp_##_FSETNAME_##_map(int n){					\
									\
    double *data, omega[3];                                             \
    double arr1[3], arr2[3], arr3[3], com[3][3], euler[3] = {0};        \
    double tmp_b1[3], tmp_b2[3], mb1, mb2;				\
									\
    double *mp = (double*)malloc(sizeof*data *				\
		(rp->input_shape[STEPS]*rp->input_shape[FEATS]));	\
									\
    real ***hr = rp->history;						\
									\
    int n_steps = rp->input_shape[STEPS];				\
    int n_feats = rp->input_shape[FEATS];				\
    int mod = n%rp->modulus;						\
    int hp = (mod == 0) ? n+1: n;					\
	    int hs = (mod == 0) ? n+2: ((mod == 1) ? n+1 : n-1);		\
	    int ox = (mod == 0) ? n  : ((mod == 1) ? n-1 : n-2);		\
	    int pa, sa, ta = hs;						\
										\
	    if (mod == 0){							\
	      pa = ox;								\
	      sa = hp;								\
	    }									\
	    else {								\
	      pa = hp;								\
	      sa = ox;								\
	    }									\
										\
	    double *c0_bp = (double *)malloc(sizeof*c0_bp*3);			\
	    double *c0_bs = (double *)malloc(sizeof*c0_bs*3);			\
	    double *c0_bo = (double *)malloc(sizeof*c0_bs*3);			\
										\
	    double *c1_bp = (double *)malloc(sizeof*c1_bp*3);			\
	    double *c1_bs = (double *)malloc(sizeof*c1_bs*3);			\
	    double *c1_bo = (double *)malloc(sizeof*c1_bs*3);			\


	//scrubbed 3/17/2023
	//i is our position in the input array
	//when the historical record is bloated, i+1 is the corresponding pos'n
	//when the historical record is NOT bloated, i is the corresponding pos'n
	#define LOOP_THROUGH_HISTORICAL_RECORD_AND				\
	    for(int i=0,j=0;i<n_steps;j=(++i)*n_feats){


	//scrubbed 3/17/2023
	#define END_iew_MAP							\
	    }									\
	    free(e1);								\
	    free(e2);								\
	    free(e3);								\
	    free(iby);								\
	    free(ibz);								\
	    free(c0_bp);							\
	    free(c0_bs);							\
	    free(c1_bp);							\
	    free(c1_bs);							\
	    return(void *)mp;							\
	 }


	//scrubbed 3/17/2023
	#define END_MAP								\
	    }									\
	    free(c0_bp);							\
	    free(c0_bs);							\
	    free(c1_bp);							\
	    free(c1_bs);							\
	    return(void *)mp;							\
	  }

	//scrubbed 3/17/2023
	#define END_MAP_AND(_X_)						\
	    }									\
	    _X_;								\
	    free(c0_bp);							\
	    free(c0_bs);							\
	    free(c1_bp);							\
	    free(c1_bs);							\
	    return(void *)mp;							\
	  }

	//scrubbed 3/17/2023
	#define RETURN_NORMALIZED_MAP						\
	     }									\
	     free(c0_bp);							\
	     free(c0_bs);							\
	     free(c1_bp);							\
	     free(c1_bs);							\
	     return rp_normalized_map(mp,pa);					\
	  }

	//scrubbed 3/17/2023
	#define RETURN_SHIFTED_MAP						\
	     }									\
	     free(c0_bp);							\
	     free(c0_bs);							\
	     free(c1_bp);							\
	     free(c1_bs);							\
	     return rp_shifted_map(mp,pa);					\
	  }

	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
					reverse maps
	\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


	//scrubbed 5-25-2022 - Wayyne
	void rp_NULL_pam(int n, const GMX_rvec x, double *forecast){
		do_nothing();
	}


	//scrubbed 5-25-2022 - Wayyne
	void rp_v_pam(int n, const GMX_rvec x, double *forecast){

	    int ref_index = n%rp->modulus; 

	    for (int i=0;i<3;i++)
		rp->forecast[n][i] = (real)(forecast[i]);
	}


	//scrubbed 5-25-2022 - Wayyne
	void rp_vC_pam(int n, const GMX_rvec x, double *forecast){
		
	    for (int i=0;i<3;i++)
		rp->forecast[n][i] = (real)(forecast[i]+forecast[i+3]);
	}

	//scrubbed 12-25-2022 - Wayyne
	void rp_u_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = forecast[i]*rp->t0_magnitude[n][0];
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_uR_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = forecast[i]*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_uRC_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = (forecast[i]+forecast[i+3])*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_DuR0_pam(int n, const GMX_rvec x, double *forecast){

	    rp->forecast[n][0] = (forecast[0]+forecast[3])*rp->t0_magnitude[n][0];
	    rp->forecast[n][1] = (forecast[1]+forecast[4])*rp->t0_magnitude[n][0];
	    rp->forecast[n][2] = (forecast[2]+forecast[5])*rp->t0_magnitude[n][0];

	     //bank(rp->forecast[n],RP_PI/4.f);
	    pitch(rp->forecast[n],RP_PI/4.f);
	      yaw(rp->forecast[n],RP_PI/4.f);
	     bank(rp->forecast[n],-rp->rotation[n][2]);
	      yaw(rp->forecast[n],-rp->rotation[n][1]);
	     bank(rp->forecast[n],-rp->rotation[n][0]);
	}

	//scrubbed 12-25-2022 - Wayyne
	void rp_9uRC_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      //rp->forecast[n][i] = (forecast[i]+forecast[i+3])*rp->t0_magnitude[n][0];
	      rp->forecast[n][i] = (forecast[i]+forecast[i+3]+forecast[i+6])*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}

	//scrubbed 12-25-2022 - Wayyne
	void rp_9uRCS_pam(int n, const GMX_rvec x, double *forecast){

	    //internal velocities (g=2 && angular)
	    rp->forecast[n][0] = (forecast[0]+forecast[3])*rp->t0_magnitude[n][0];
	    rp->forecast[n][1] = (forecast[1]+forecast[4])*rp->t0_magnitude[n][0];
	    rp->forecast[n][2] = (forecast[2]+forecast[5])*rp->t0_magnitude[n][0];

	    //rotate back to the external frame
	    bank(rp->forecast[n],-rp->rotation[n][2]);
	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);

	    //external velocities (CoM)
	    double Vcom[3];
	    Vcom[0] = forecast[3]*rp->t0_magnitude[n][1];
	    Vcom[1] = forecast[4]*rp->t0_magnitude[n][1];
	    Vcom[2] = forecast[5]*rp->t0_magnitude[n][1];

	    //rotate back to the external frame
	     yaw(rp->forecast[n],-rp->rotation[n][4]);
	    bank(rp->forecast[n],-rp->rotation[n][3]);
		
	    //combine internal and external for final forecast
	    rp->forecast[n][0] += Vcom[0];
	    rp->forecast[n][1] += Vcom[1];
	    rp->forecast[n][2] += Vcom[2];
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_9uRtC_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = (forecast[i]+forecast[i+3]+forecast[i+6])*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}



	//scrubbed 12-25-2022 - Wayyne
	void rp_9uRC2_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = ((forecast[i]*rp->t0_magnitude[n][1])+forecast[i+3])*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}

	//scrubbed 12-25-2022 - Wayyne
	void rp_6uRC_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = (forecast[i]*rp->t0_magnitude[n][0]) + (forecast[i+3]*rp->t0_magnitude[n][1]);

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_uRCS_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = ((forecast[i]+rp->t0_magnitude[n][i+1])+(forecast[i+3]+rp->t0_magnitude[n][i+4]))*rp->t0_magnitude[n][0];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}

	//scrubbed 12-25-2022 - Wayyne
	void rp_NSuR_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++){
	      rp->forecast[n][i] = forecast[i]*rp->rotation[n][(i*2)+3];
	      rp->forecast[n][i] += rp->rotation[n][(i*2)+4];
	      rp->forecast[n][i] *= rp->rotation[n][2];
	    }

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_vR_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = forecast[i];

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_vRs_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++)
	      rp->forecast[n][i] = forecast[i];

	    bank(rp->forecast[n],-rp->rotation[n][2]);
	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	//scrubbed 12-25-2022 - Wayyne
	void rp_NSvR_pam(int n, const GMX_rvec x, double *forecast){

	    for (int i=0;i<3;i++){
	      rp->forecast[n][i] = forecast[i]*rp->rotation[n][(i*2)+3];
	      rp->forecast[n][i] += rp->rotation[n][(i*2)+4];
	    }

	     yaw(rp->forecast[n],-rp->rotation[n][1]);
	    bank(rp->forecast[n],-rp->rotation[n][0]);
	}


	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
					forward maps
	\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


	/*raw gmx velocities*/
	BEGIN_MAP(3v)
	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_VELOCITIES_FOR(pa,0,i)
	    END_MAP


	/*gmx velocities rotated s.t. V_t0 is aligned with the x-axis*/
	BEGIN_MAP(3vR0)

	    //prime the rotation w.r.t. t=0
	    arr1[0] = (double)hr[pa][0][0];					
	    arr1[1] = (double)hr[pa][0][1];					
	    arr1[2] = (double)hr[pa][0][2];					


	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);
	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    fprintf(stdout,"md: % -18.6lf  % -18.6lf  % -18.6lf ang: % -18.6f  % -18.6f\n",
                (double)hr[pa][n_steps-1][0],
                (double)hr[pa][n_steps-1][1],
                (double)hr[pa][n_steps-1][2],
		rp->rotation[pa][1],
		rp->rotation[pa][0]);
       
	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_VELOCITIES_FOR(pa,0,i) //was i+i....but why?

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j]),rp->rotation[pa][0]);
	       yaw(&(mp[j]),rp->rotation[pa][1]);
	
	      mp[j+3] = rp->rotation[pa][0];
	      mp[j+4] = rp->rotation[pa][1];

	    END_MAP


	/*3vRv0 but normalized and standardized s.t.
	  for each feature the mean is 0 and the stddev is 1*/
	BEGIN_MAP(ns3vRv0)

	    //prime the rotation w.r.t. t=0
	    arr1[0] = (double)hr[pa][0][0];					
	    arr1[1] = (double)hr[pa][0][1];					
	    arr1[2] = (double)hr[pa][0][2];					

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);
	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_VELOCITIES_FOR(pa,0,i) //was i+i....but why?

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j]),rp->rotation[pa][0]);
	       yaw(&(mp[j]),rp->rotation[pa][1]);

	    RETURN_NORMALIZED_MAP


	BEGIN_MAP(3u)
	    PRIME_UELOCITIES_FOR(pa);
	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_UELOCITIES_FOR(pa,0,i)
	    END_MAP


	BEGIN_MAP(3uR0)

	    PRIME_UELOCITIES_FOR(pa);

	    arr1[0] = (double)hr[pa][0][0];					
	    arr1[1] = (double)hr[pa][0][1];					
	    arr1[2] = (double)hr[pa][0][2];					

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	
	    fprintf(stdout,"md: % -18.6lf  % -18.6lf  % -18.6lf ang: % -18.6f  % -18.6f\n",
	        (double)hr[pa][n_steps-1][0],
	        (double)hr[pa][n_steps-1][1],
	        (double)hr[pa][n_steps-1][2],
		rp->rotation[pa][1],
		rp->rotation[pa][0]);

	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_UELOCITIES_FOR(pa,0,i)

	       //rotate the velocity into the internal, inertial frame of reference
	       bank(&(mp[j]),rp->rotation[pa][0]);
		yaw(&(mp[j]),rp->rotation[pa][1]);

	      mp[j+3] = rp->rotation[pa][0];
              mp[j+4] = rp->rotation[pa][1];
              mp[j+5] = rp->t0_magnitude[pa][0];

	    END_MAP


	//scrubbed 3/18/2023 
	#define CAPTURE_COM_Vdt(_X_,_Y_)					\
	    FIND_COM(com[0],_Y_-1);						\
	    FIND_COM(com[1],_Y_);						\
	    _X_[0] = com[1][0] - com[0][0];					\
	    _X_[1] = com[1][1] - com[0][1];					\
	    _X_[2] = com[1][2] - com[0][2];					\
		

	#define CAPTURE_COM_Vgmx(_X_,_Y_)					\
	    CAPTURE_COM_Vdt(_X_,_Y_)						\
	    _X_[0] *= rp->dt_per_ps;						\
	    _X_[1] *= rp->dt_per_ps;						\
	    _X_[2] *= rp->dt_per_ps;						\


	BEGIN_MAP(3uCoMRv0)

	    //get magnitude at t=0 for normalization
	    CAPTURE_COM_Vgmx(com[2],1); 

	    arr1[0] = (double)hr[pa][1][0]-com[2][0];					
	    arr1[1] = (double)hr[pa][1][1]-com[2][1];					
	    arr1[2] = (double)hr[pa][1][2]-com[2][2];					

	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	 
	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      mp[j+0] = (double)((hr[pa][i+1][0]-com[2][0])/rp->t0_magnitude[pa][0]);
	      mp[j+1] = (double)((hr[pa][i+1][1]-com[2][1])/rp->t0_magnitude[pa][0]);
	      mp[j+2] = (double)((hr[pa][i+1][2]-com[2][2])/rp->t0_magnitude[pa][0]);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j]),rp->rotation[pa][0]);
	       yaw(&(mp[j]),rp->rotation[pa][1]);

	    END_MAP_AND(free(eas))


	BEGIN_MAP(D3v)

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //inbond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	    END_MAP

	BEGIN_MAP(D3uRv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);

	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];

	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];

	    END_MAP


	BEGIN_MAP(D6uRv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    //get magnitude at t=0 for normalization
	    CAPTURE_COM_Vgmx(com[2],1);

	    //A = Vtotal - Vcom
	    arr1[0] = (double)hr[pa][1][0]-com[2][0];
	    arr1[1] = (double)hr[pa][1][1]-com[2][1];
	    arr1[2] = (double)hr[pa][1][2]-com[2][2];

	    //B = R - CoM
	    *(c1_bp+0) = (double)(hr[pa][1][3]-com[1][0]);
	    *(c1_bp+1) = (double)(hr[pa][1][4]-com[1][1]);
	    *(c1_bp+2) = (double)(hr[pa][1][5]-com[1][2]);

	    //Bhat = B/mag(B)
	    double m_bp = mag(c1_bp,3);
	    *(c1_bp+0) /= m_bp;
	    *(c1_bp+1) /= m_bp;
	    *(c1_bp+2) /= m_bp;

	    //Pab = proj of A on B (in-bond velocity)
	    double d = dot(c1_bp,arr1,3);
	    *(c1_bp+0) *= d;
	    *(c1_bp+1) *= d;
	    *(c1_bp+2) *= d;

	    //angular velocity (A - Pab)
	    arr2[0] = arr1[0] - *(c1_bp+0);
	    arr2[1] = arr1[1] - *(c1_bp+1);
	    arr2[2] = arr1[2] - *(c1_bp+2);

	    rp->t0_magnitude[pa][0] = mag(arr2,3);
	    rp->t0_magnitude[pa][1] = mag(c1_bp,3);

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);

	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];

	      mp[j+3] /= rp->t0_magnitude[pa][1];
	      mp[j+4] /= rp->t0_magnitude[pa][1];
	      mp[j+5] /= rp->t0_magnitude[pa][1];

	    END_MAP

	BEGIN_MAP(D12Rv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);
	      bank(&(mp[j+6]),rp->rotation[pa][0]);
	       yaw(&(mp[j+6]),rp->rotation[pa][1]);

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];
	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];
	      mp[j+6] /= rp->t0_magnitude[pa][0];
	      mp[j+7] /= rp->t0_magnitude[pa][0];
	      mp[j+8] /= rp->t0_magnitude[pa][0];

	    END_MAP

	BEGIN_MAP(D9Rv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = (double)(arr1[0] - *(c1_bp+0));
	      mp[j+1] = (double)(arr1[1] - *(c1_bp+1));
	      mp[j+2] = (double)(arr1[2] - *(c1_bp+2));

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);
	      bank(&(mp[j+6]),rp->rotation[pa][0]);
	       yaw(&(mp[j+6]),rp->rotation[pa][1]);

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];

	      rp->t0_magnitude[pa][1] = mag(&(mp[j]),3);

	      mp[j+0] /= rp->t0_magnitude[pa][1];
	      mp[j+1] /= rp->t0_magnitude[pa][1];
	      mp[j+2] /= rp->t0_magnitude[pa][1];

	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];
	      mp[j+6] /= rp->t0_magnitude[pa][0];
	      mp[j+7] /= rp->t0_magnitude[pa][0];
	      mp[j+8] /= rp->t0_magnitude[pa][0];
	      mp[j+9]  = rp->t0_magnitude[pa][1];

	    END_MAP


	BEGIN_MAP(C9dRv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);
	      bank(&(mp[j+6]),rp->rotation[pa][0]);
	       yaw(&(mp[j+6]),rp->rotation[pa][1]);

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];
	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];
	      mp[j+6] /= rp->t0_magnitude[pa][0];
	      mp[j+7] /= rp->t0_magnitude[pa][0];
	      mp[j+8] /= rp->t0_magnitude[pa][0];
	   }

	   //the cumulative bit
	   LOOP_THROUGH_HISTORICAL_RECORD_AND
	       for(int k=i+1,J=(i+1)*n_feats;k<n_steps;J=(++k)*n_feats)
		 for (int h=0;h<9;h++)
		   mp[j+h] += mp[J+h]; 
	    END_MAP

	BEGIN_MAP(C9uRt0)

	    //get Vcom in GMX internal units (ps/nm)
	    CAPTURE_COM_Vgmx(com[2],1); 

	    //decompose Vtotal: A = Vtotal - Vcom
	    arr1[0] = (double)hr[pa][1][0]-com[2][0];					
	    arr1[1] = (double)hr[pa][1][1]-com[2][1];					
	    arr1[2] = (double)hr[pa][1][2]-com[2][2];					

	    //decompose Vtotal: B = R - CoM
	    *(c1_bp+0) = (double)(hr[pa][1][3]-com[1][0]);
	    *(c1_bp+1) = (double)(hr[pa][1][4]-com[1][1]);
	    *(c1_bp+2) = (double)(hr[pa][1][5]-com[1][2]);

	    //decompose Vtotal: Bhat = B/mag(B)
	    double m_bp = mag(c1_bp,3);
	    *(c1_bp+0) /= m_bp;
	    *(c1_bp+1) /= m_bp;
	    *(c1_bp+2) /= m_bp;
	 
	    //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	    double d = dot(c1_bp,arr1,3);
	    *(c1_bp+0) *= d;
	    *(c1_bp+1) *= d;
	    *(c1_bp+2) *= d;

	    //decompose Vtotal: W = A - Pab (angular/tangential velocity)
	    arr2[0] = arr1[0] - *(c1_bp+0);
	    arr2[1] = arr1[1] - *(c1_bp+1);
	    arr2[2] = arr1[2] - *(c1_bp+2);

	    //get angles to rotate t=0 decomposed velocities into internal frame
	    double *eas = xzx_angles(c1_bp,arr2,false,false);
	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	    rp->rotation[pa][2] = eas[2];

	    //capture Vtotal's t=0 magitude
	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get Vcom in GMX internal units (ps/nm)
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //decompose Vtotal: A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //decompose Vtotal: B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //decompose Vtotal: Bhat = B/mag(B)
	      m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	      d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //render angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //render in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //render CoM velocity
	      mp[j+6] = com[2][0] == 0.f ? 0.f : 1.0/com[2][0];
	      mp[j+7] = com[2][1] == 0.f ? 0.f : 1.0/com[2][1];
	      mp[j+8] = com[2][2] == 0.f ? 0.f : 1.0/com[2][2];

	      //rotate s.t. the in-bond velocity is on x_hat and angular velocity on y_hat at t=0
	      for (int k=0;k<3;k++){
		bank(&(mp[j+(3*k)]),rp->rotation[pa][0]);
		 yaw(&(mp[j+(3*k)]),rp->rotation[pa][1]);
		bank(&(mp[j+(3*k)]),rp->rotation[pa][2]);
	      }

	      //rescale s.t. t=0 is unit velocity
	      for (int k=0;k<9;k++)
		mp[j+k] /= rp->t0_magnitude[pa][0];
	   }

	   //the cumulative bit
	   LOOP_THROUGH_HISTORICAL_RECORD_AND
	       for(int k=i+1,J=(i+1)*n_feats;k<n_steps;J=(++k)*n_feats)
		 for (int h=0;h<9;h++)
		   //only accumulate CoM and angular velocities
		   if (h>5)
		     mp[j+h] += mp[J+h]; 
	       //take the inverse of the cumulative CoM velocity
	       for (int h=6;h<9;h++)
		 mp[j+h] = 1.f/mp[j+h];
	    END_MAP


	BEGIN_MAP(D9uRt0)
	    //get Vcom in GMX internal units (ps/nm)
	    CAPTURE_COM_Vgmx(com[2],1); 

	    //decompose Vtotal: A = Vtotal - Vcom
	    arr1[0] = (double)hr[pa][1][0]-com[2][0];					
	    arr1[1] = (double)hr[pa][1][1]-com[2][1];					
	    arr1[2] = (double)hr[pa][1][2]-com[2][2];					

	    //capture magnitude of A and Vcom at t=0
	    rp->t0_magnitude[pa][0] = mag(arr1,3);
	    rp->t0_magnitude[pa][1] = mag(com[2],3);

	    //decompose Vtotal: B = R - CoM
	    *(c1_bp+0) = (double)(hr[pa][1][3]-com[1][0]);
	    *(c1_bp+1) = (double)(hr[pa][1][4]-com[1][1]);
	    *(c1_bp+2) = (double)(hr[pa][1][5]-com[1][2]);

	    //decompose Vtotal: Bhat = B/mag(B)
	    double m_bp = mag(c1_bp,3);
	    *(c1_bp+0) /= m_bp;
	    *(c1_bp+1) /= m_bp;
	    *(c1_bp+2) /= m_bp;
	 
	    //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	    double d = dot(c1_bp,arr1,3);
	    arr3[0] = *(c1_bp+0) *= d;
	    arr3[1] = *(c1_bp+1) *= d;
	    arr3[2] = *(c1_bp+2) *= d;

	    //decompose Vtotal: W = A - Pab (angular/tangential velocity)
	    arr2[0] = arr1[0] - *(c1_bp+0);
	    arr2[1] = arr1[1] - *(c1_bp+1);
	    arr2[2] = arr1[2] - *(c1_bp+2);

	    //get angles to rotate t=0 decomposed velocities into internal frame
	    double *eas = xzx_angles(arr3,arr2,false,false);
	    rp->rotation[pa][0] = eas[0];//0.785;
	    rp->rotation[pa][1] = eas[1];//0.785;
	    rp->rotation[pa][2] = eas[2];//0.785;

	    //get angles to rotate t=0 Vcom onto x-axis
	    eas = xz_angles(com[2],false);
	    rp->rotation[pa][3] = eas[0];
	    rp->rotation[pa][4] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get Vcom in GMX internal units (ps/nm)
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //decompose Vtotal: A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //decompose Vtotal: B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //decompose Vtotal: Bhat = B/mag(B)
	      m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	      d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //render angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //render in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //render CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate s.t. the in-bond velocity is on x_hat and angular velocity on y_hat at t=0
	      for (int k=0;k<2;k++){
		bank(&(mp[j+(3*k)]),rp->rotation[pa][0]);
		 yaw(&(mp[j+(3*k)]),rp->rotation[pa][1]);
		bank(&(mp[j+(3*k)]),rp->rotation[pa][2]);
	      }

	      //rotate the Vcom s.t. t=0 Vcom is aligned with the x-ais
	      bank(&(mp[j+6]),rp->rotation[pa][3]);
	       yaw(&(mp[j+6]),rp->rotation[pa][4]);

	      //rescale s.t. t=0 is unit velocity
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];
	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];

	      mp[j+6] /= rp->t0_magnitude[pa][1];
	      mp[j+7] /= rp->t0_magnitude[pa][1];
	      mp[j+8] /= rp->t0_magnitude[pa][1];

	    END_MAP

	BEGIN_MAP(D9uRb0)
	    //get Vcom in GMX internal units (ps/nm)
	    CAPTURE_COM_Vgmx(com[2],1); 

	    //decompose Vtotal: A = Vtotal - Vcom
	    arr1[0] = (double)hr[pa][1][0]-com[2][0];					
	    arr1[1] = (double)hr[pa][1][1]-com[2][1];					
	    arr1[2] = (double)hr[pa][1][2]-com[2][2];					

	    //capture magnitude of A and Vcom at t=0
	    rp->t0_magnitude[pa][0] = mag(arr1,3);
	    rp->t0_magnitude[pa][1] = mag(com[2],3);

	    //decompose Vtotal: B = R - CoM
	    arr3[0] = *(c1_bp+0) = (double)(hr[pa][1][3]-com[1][0]);
	    arr3[1] = *(c1_bp+1) = (double)(hr[pa][1][4]-com[1][1]);
	    arr3[2] = *(c1_bp+2) = (double)(hr[pa][1][5]-com[1][2]);

	    //decompose Vtotal: Bhat = B/mag(B)
	    double m_bp = mag(c1_bp,3);
	    *(c1_bp+0) /= m_bp;
	    *(c1_bp+1) /= m_bp;
	    *(c1_bp+2) /= m_bp;
	 
	    //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	    double d = dot(c1_bp,arr1,3);
	    *(c1_bp+0) *= d;
	    *(c1_bp+1) *= d;
	    *(c1_bp+2) *= d;

	    //decompose Vtotal: W = A - Pab (angular/tangential velocity)
	    arr2[0] = arr1[0] - *(c1_bp+0);
	    arr2[1] = arr1[1] - *(c1_bp+1);
	    arr2[2] = arr1[2] - *(c1_bp+2);

	    //get angles to rotate t=0 decomposed velocities into internal frame
	    double *eas = xzx_angles(arr3,arr2,false,false);
	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	    rp->rotation[pa][2] = eas[2];

	    //get angles to rotate t=0 Vcom onto x-axis
	    eas = xz_angles(com[2],false);
	    rp->rotation[pa][3] = eas[0];
	    rp->rotation[pa][4] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get Vcom in GMX internal units (ps/nm)
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //decompose Vtotal: A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //decompose Vtotal: B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //decompose Vtotal: Bhat = B/mag(B)
	      m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //decompose Vtotal: Pab = proj of A on B (in-bond velocity)
	      d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //render angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //render in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //render CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate s.t. the in-bond velocity is on x_hat and angular velocity on y_hat at t=0
	      for (int k=0;k<2;k++){
		bank(&(mp[j+(3*k)]),rp->rotation[pa][0]);
		 yaw(&(mp[j+(3*k)]),rp->rotation[pa][1]);
		bank(&(mp[j+(3*k)]),rp->rotation[pa][2]);
	      }

	      //rotate the Vcom s.t. t=0 Vcom is aligned with the x-ais
	      bank(&(mp[j+6]),rp->rotation[pa][3]);
	       yaw(&(mp[j+6]),rp->rotation[pa][4]);

	      //rescale s.t. t=0 is unit velocity
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];
	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];

	      mp[j+6] /= rp->t0_magnitude[pa][1];
	      mp[j+7] /= rp->t0_magnitude[pa][1];
	      mp[j+8] /= rp->t0_magnitude[pa][1];

	    END_MAP

	BEGIN_MAP(DuR0_BACKUP_WORKING_3_15_2024)

	    //we must start with t=1 because we must compute v_com
	    //from positions in the historical record
	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	    
	   // bank(&(arr1[0]),rp->rotation[pa][0]);
	    // yaw(&(arr1[0]),rp->rotation[pa][1]);

	    arr1[0] = (double)(hr[pa][1][3] - hr[sa][1][3]);					
	    arr1[1] = (double)(hr[pa][1][4] - hr[sa][1][4]);					
	    arr1[2] = (double)(hr[pa][1][5] - hr[sa][1][5]);					

	    //what's happening to the bond? TODO: is an additional rotation required?
	   // printf("bond befr: %f %f %f\n",arr1[0],arr1[1],arr1[2]);
	   // bank(&(arr1[0]),rp->rotation[pa][0]);
	    // yaw(&(arr1[0]),rp->rotation[pa][1]);
	    //printf("bond aftr: %f %f %f\n",arr1[0],arr1[1],arr1[2]);

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get Vcom
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = V_internal = V_total - V_com
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM: fix is here, we're taking the bond from the prev state
	      *(c1_bp+0) = (double)(hr[pa][i][3]-com[0][0]);
	      *(c1_bp+1) = (double)(hr[pa][i][4]-com[0][1]);
	      *(c1_bp+2) = (double)(hr[pa][i][5]-com[0][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity) assuming no rotation
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d; 
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM velocity
	//      mp[j+6] = com[2][0];
	//      mp[j+7] = com[2][1];
	//      mp[j+8] = com[2][2];

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);
	      
	      //try to rotate the t0 velocit so its split between x y and z
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);

	//      bank(&(mp[j+6]),rp->rotation[pa][0]);
	//       yaw(&(mp[j+6]),rp->rotation[pa][1]);

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0]; //rotX
	      mp[j+1] /= rp->t0_magnitude[pa][0]*pow(-1.f,pa); //rotY - pow() is a hack to reflect atm2 only
	      mp[j+2] /= rp->t0_magnitude[pa][0]; //rotZ

	      mp[j+3] /= rp->t0_magnitude[pa][0]; //vibX
	      mp[j+4] /= rp->t0_magnitude[pa][0]*pow(-1.f,pa); //vibY
	      mp[j+5] /= rp->t0_magnitude[pa][0]; //vibZ
	      
	//      mp[j+6] /= rp->t0_magnitude[pa][0]; //comX
	//      mp[j+7] /= rp->t0_magnitude[pa][0]; //comY
	//      mp[j+8] /= rp->t0_magnitude[pa][0]; //comZ

	    END_MAP


	BEGIN_MAP(DuR0)

	    //++++++++++ NOTE ++++++++++++
	    // We must start with t=1 because we must compute v_com
	    //++++++++++ NOTE ++++++++++++
	  
            //get velocity at t=1 
	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					

	    //compute the velocities magnitude at t=1
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get CoM at t=0 (the first bond we project Vvib onto)
	    FIND_COM(arr2,0);

	    //get CoM->atom vector and make the CoM the origin
	    arr3[0] = (double)hr[pa][0][3] - arr2[0];
	    arr3[1] = (double)hr[pa][0][4] - arr2[1];
	    arr3[2] = (double)hr[pa][0][5] - arr2[2];

	    //get the X-Z-X euler angles at t=1 to:
	    // (1) place the bond vector at t=0 on the X-axis
	    // (2) place the velocity vector at t=1 in the XY-plane
	    double *eas = xzx_angles(arr3,arr1,false,-1);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	    rp->rotation[pa][2] = eas[2];
	
	    fprintf(stdout,"md: % -18.6lf  % -18.6lf  % -18.6lf ang: % -18.6f  % -18.6f\n",
                (double)hr[pa][n_steps-1][0],
                (double)hr[pa][n_steps-1][1],
                (double)hr[pa][n_steps-1][2],
		rp->rotation[pa][1],
		rp->rotation[pa][0]);

	    //we're ready to roll, construct the model input   
	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get Vcom
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = V_internal = V_total - V_com
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM: fix is here, we're taking the bond from the prev state
	      *(c1_bp+0) = (double)(hr[pa][i][3]-com[0][0]);
	      *(c1_bp+1) = (double)(hr[pa][i][4]-com[0][1]);
	      *(c1_bp+2) = (double)(hr[pa][i][5]-com[0][2]);

	      //current com->atm vector, needed for pam
	     //*(c0_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	     //*(c0_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	     //*(c0_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);
	     //rp->t0_magnitude[pa][1] = mag(c0_bp,3);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity) assuming no rotation
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d; 
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular/residual velocity (A - Pab)
	      mp[j+0] = (arr1[0] - *(c1_bp+0));
	      mp[j+1] = (arr1[1] - *(c1_bp+1));
	      mp[j+2] = (arr1[2] - *(c1_bp+2));
	      //mp[j+2] = (arr1[2]); //TODO: testing using the composed internal Z component

	      //in-bond/vibrational velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM/translational velocity
	      //mp[j+6] = com[2][0];
	      //mp[j+7] = com[2][1];
	      //mp[j+8] = com[2][2];

	      //rotate into the internal, inertial frame of reference
	       bank(&(mp[j+0]),rp->rotation[pa][0]);
		yaw(&(mp[j+0]),rp->rotation[pa][1]);
	       bank(&(mp[j+0]),rp->rotation[pa][2]);
	        yaw(&(mp[j+0]),-RP_PI/4.f);
	      pitch(&(mp[j+0]),-RP_PI/4.f);
	       //bank(&(mp[j+0]),-RP_PI/4.f);//TODO: testing
      
	       bank(&(mp[j+3]),rp->rotation[pa][0]);
		yaw(&(mp[j+3]),rp->rotation[pa][1]);
	       bank(&(mp[j+3]),rp->rotation[pa][2]);
	        yaw(&(mp[j+3]),-RP_PI/4.f);
	      pitch(&(mp[j+3]),-RP_PI/4.f);
	       //bank(&(mp[j+3]),-RP_PI/4.f); //TODO: testing

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0]; //rotX
	      mp[j+1] /= rp->t0_magnitude[pa][0]; //rotY 
	      mp[j+2] /= rp->t0_magnitude[pa][0]; //rotZ

	      mp[j+3] /= rp->t0_magnitude[pa][0]; //vibX
	      mp[j+4] /= rp->t0_magnitude[pa][0]; //vibY
	      mp[j+5] /= rp->t0_magnitude[pa][0]; //vibZ
	      
	      //mp[j+2] += mp[j+5];
	      //mp[j+6] /= rp->t0_magnitude[pa][0]; //comX
	      //mp[j+7] /= rp->t0_magnitude[pa][0]; //comY
	      //mp[j+8] /= rp->t0_magnitude[pa][0]; //comZ

	      mp[j+6] = rp->rotation[pa][0];
              mp[j+7] = rp->rotation[pa][1];
              mp[j+8] = rp->t0_magnitude[pa][0];
              mp[j+9] = rp->rotation[pa][2];
              mp[j+10] = -RP_PI/4.f;
              //mp[j+8] = (double)hr[pa][1][0];
              //mp[j+9] = (double)hr[pa][1][1];
              //mp[j+10] = (double)hr[pa][1][2];

	    END_MAP


	BEGIN_MAP(D9uRv0)

	    //we must start with t=1 because we must compute v_com
	    //from positions in the historical record
	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //CoM velocity
	      mp[j+6] = com[2][0];
	      mp[j+7] = com[2][1];
	      mp[j+8] = com[2][2];

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);
	      bank(&(mp[j+6]),rp->rotation[pa][0]);
	       yaw(&(mp[j+6]),rp->rotation[pa][1]);

	      //rescale
	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];
	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];
	      mp[j+6] /= rp->t0_magnitude[pa][0];
	      mp[j+7] /= rp->t0_magnitude[pa][0];
	      mp[j+8] /= rp->t0_magnitude[pa][0];

	    END_MAP

	BEGIN_MAP(D3uSRv0)

	    arr1[0] = (double)hr[pa][1][0];					
	    arr1[1] = (double)hr[pa][1][1];					
	    arr1[2] = (double)hr[pa][1][2];					
	    
	    rp->t0_magnitude[pa][0] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];

	    LOOP_THROUGH_HISTORICAL_RECORD_AND

	      //get magnitude at t=0 for normalization
	      CAPTURE_COM_Vgmx(com[2],i+1); 

	      //A = Vtotal - Vcom
	      arr1[0] = (double)hr[pa][i+1][0]-com[2][0];					
	      arr1[1] = (double)hr[pa][i+1][1]-com[2][1];					
	      arr1[2] = (double)hr[pa][i+1][2]-com[2][2];					

	      //B = R - CoM
	      *(c1_bp+0) = (double)(hr[pa][i+1][3]-com[1][0]);
	      *(c1_bp+1) = (double)(hr[pa][i+1][4]-com[1][1]);
	      *(c1_bp+2) = (double)(hr[pa][i+1][5]-com[1][2]);

	      //Bhat = B/mag(B)
	      double m_bp = mag(c1_bp,3);
	      *(c1_bp+0) /= m_bp;
	      *(c1_bp+1) /= m_bp;
	      *(c1_bp+2) /= m_bp;
	 
	      //Pab = proj of A on B (in-bond velocity)
	      double d = dot(c1_bp,arr1,3);
	      *(c1_bp+0) *= d;
	      *(c1_bp+1) *= d;
	      *(c1_bp+2) *= d;

	      //angular velocity (A - Pab)
	      mp[j+0] = arr1[0] - *(c1_bp+0);
	      mp[j+1] = arr1[1] - *(c1_bp+1);
	      mp[j+2] = arr1[2] - *(c1_bp+2);

	      //in-bond velocity
	      mp[j+3] = *(c1_bp+0);
	      mp[j+4] = *(c1_bp+1);
	      mp[j+5] = *(c1_bp+2);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+0]),rp->rotation[pa][0]);
	       yaw(&(mp[j+0]),rp->rotation[pa][1]);

	      //rotate the velocity into the internal, inertial frame of reference
	      bank(&(mp[j+3]),rp->rotation[pa][0]);
	       yaw(&(mp[j+3]),rp->rotation[pa][1]);

	      mp[j+0] /= rp->t0_magnitude[pa][0];
	      mp[j+1] /= rp->t0_magnitude[pa][0];
	      mp[j+2] /= rp->t0_magnitude[pa][0];

	      mp[j+3] /= rp->t0_magnitude[pa][0];
	      mp[j+4] /= rp->t0_magnitude[pa][0];
	      mp[j+5] /= rp->t0_magnitude[pa][0];

	    RETURN_SHIFTED_MAP

	BEGIN_MAP(ns3uRv0)

	    //prime the rotation w.r.t. t=0
	    arr1[0] = (double)hr[pa][0][0];					
	    arr1[1] = (double)hr[pa][0][1];					
	    arr1[2] = (double)hr[pa][0][2];					

	    //get magnitude at t=0 for normalization
	    double mt0 = rp->rotation[pa][2] = mag(arr1,3);

	    //get euler angles at t=0 
	    double *eas = xz_angles(arr1,FALSE);

	    rp->rotation[pa][0] = eas[0];
	    rp->rotation[pa][1] = eas[1];
	 
	    LOOP_THROUGH_HISTORICAL_RECORD_AND
	      RENDER_UELOCITIES_FOR(pa,0,i)

	       //rotate the velocity into the internal, inertial frame of reference
	       bank(&(mp[j]),rp->rotation[pa][0]);
		yaw(&(mp[j]),rp->rotation[pa][1]);

	    RETURN_NORMALIZED_MAP


	BEGIN_MAP(iew_v0)

	  double ibx[3], *iby, *ibz;
	  double *e1, *e2, *e3;
	  double  mx ,my, mz;

	  FIND_COM(com[0],0);

	  /*determine initial structure for the i's*/
	  {
	    *(c0_bo+0) = (double)(hr[ox][0][3]-com[0][0]);
	    *(c0_bo+1) = (double)(hr[ox][0][4]-com[0][1]);
	    *(c0_bo+2) = (double)(hr[ox][0][5]-com[0][2]);
	    
	    *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	    *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	    *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);
	    
	    *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	    *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	    *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);
	    
	    arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));
	  }

	  FIND_EQCOM(com[1],0);

	  /*determine initial rodrigues' angles for the w's*/
	  {
	    ibx[0] = (double)(hr[ox][0][3]-com[1][0]);
	    ibx[1] = (double)(hr[ox][0][4]-com[1][1]);
	    ibx[2] = (double)(hr[ox][0][5]-com[1][2]);

	    ibz = cross(ibx,c0_bs);
	    iby = cross(ibx,ibz);
	    
	    mx = mag(ibx,3);
	    my = mag(iby,3);
	    mz = mag(ibz,3);
	    
	    for (int i=0;i<3;i++){
	      ibx[i] /= mx;
	      iby[i] /= my;
	      ibz[i] /= -mz;
	    }
	    
	    e1 = xzx_angles(iby,ibz,false,false); //pitch
	    e2 = xzx_angles(ibz,iby,false,false); //yaw
	    e3 = xzx_angles(ibx,iby,false,false); //bank
	  }

	  /*determine initial euler angles for the e's*/
	  {
	    omega[0] = e3[0];
	    omega[1] = e3[1];
	    omega[2] = e3[2];
	  }

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
		FIND_COM(com[0],i+1);
		FIND_EQCOM(com[2],i+1);

		/*determine current e's (features 0-2)*/
		{
		  mp[j+0] = com[2][0]-com[1][0];
		  mp[j+1] = com[2][1]-com[1][1];
		  mp[j+2] = com[2][2]-com[1][2];
	   
		   bank(&(mp[j]),omega[0]);
		    yaw(&(mp[j]),omega[1]);
		   bank(&(mp[j]),omega[2]);
		}

		/*determine current i's (features 6-9)*/
		{
		  *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[0][0]);
		  *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[0][1]);
		  *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[0][2]);
		  
		  *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[0][0]);
		  *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[0][1]);
		  *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[0][2]);
		  
		  *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[0][0]);
		  *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[0][1]);
		  *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[0][2]);
		  
		  arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

		  mp[j+6] = mag(c1_bp,3)-mag(c0_bp,3);
		  mp[j+7] = arr3[1] - arr3[0];
		  mp[j+8] = mag(c1_bs,3)-mag(c0_bs,3);
		  mp[j+9] = mag(c1_bo,3)-mag(c0_bo,3);
		}

		/*determine current w's (features 3-5)*/
		{
		  ibx[0] = (double)(hr[ox][i+1][3]-com[2][0]);
		  ibx[1] = (double)(hr[ox][i+1][4]-com[2][1]);
		  ibx[2] = (double)(hr[ox][i+1][5]-com[2][2]);

		  free(ibz);
		  free(iby);

		  ibz = cross(ibx,c1_bs);
		  iby = cross(ibx,ibz);

		  mx = mag(ibx,3);
		  my = mag(iby,3);
		  mz = mag(ibz,3);

		  for (int d=0;d<3;d++){
		    ibx[d] /= mx;
		    iby[d] /= my;
		    ibz[d] /= -mz;
		  }

		  mp[j+3] = xzx_delta(iby,ibz,e1,false);
		  mp[j+4] = xzx_delta(ibz,iby,e2,false);
		  mp[j+5] = xzx_delta(ibx,iby,e3,false);

		  free(e1);
		  free(e2);
		  free(e3);

		  e1 = xzx_angles(iby,ibz,false,false); //pitch
		  e2 = xzx_angles(ibz,iby,false,false); //yaw
		  e3 = xzx_angles(ibx,iby,false,false); //bank
		}

	       /*copy data for the next iteration*/
	       {
		 arr3[0] = arr3[1]; 			//for i's
		 for (int q=0;q<3;q++){
		   *(c0_bp+q) = *(c1_bp+q);		//for i's
		   *(c0_bs+q) = *(c1_bs+q);		//for i's
		   *(c0_bo+q) = *(c1_bo+q);		//for i's

		   com[1][q] = com[2][q];		//for e's
		 }
	       }

	   END_iew_MAP


	BEGIN_MAP(iew_v1)

	  double ibx[3], *iby, *ibz;
	  double *e1, *e2, *e3;
	  double  mx ,my, mz;

	  /*determine initial structure for the i's*/
	  {
	    FIND_COM(com[0],0);

	    *(c0_bo+0) = (double)(hr[ox][0][3]-com[0][0]);
	    *(c0_bo+1) = (double)(hr[ox][0][4]-com[0][1]);
	    *(c0_bo+2) = (double)(hr[ox][0][5]-com[0][2]);
	    
	    *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	    *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	    *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);
	    
	    *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	    *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	    *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);
	    
	    arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));
	  }

	  /*determine initial rodrigues' angles for the w's*/
	  {
	    FIND_EQCOM(com[1],0);

	    ibx[0] = (double)(hr[ox][0][3]-com[1][0]);
	    ibx[1] = (double)(hr[ox][0][4]-com[1][1]);
	    ibx[2] = (double)(hr[ox][0][5]-com[1][2]);

	    ibz = cross(ibx,c0_bs);
	    iby = cross(ibx,ibz);
	    
	    mx = mag(ibx,3);
	    my = mag(iby,3);
	    mz = mag(ibz,3);
	    
	    for (int i=0;i<3;i++){
	      ibx[i] /= mx;
	      iby[i] /= my;
	      ibz[i] /= -mz;
	    }
	    
	    e1 = xzx_angles(iby,ibz,false,false); //pitch
	    e2 = xzx_angles(ibz,iby,false,false); //yaw
	    e3 = xzx_angles(ibx,iby,false,false); //bank
	  }

	  /*get the t=0 magnitude of e for normalization*/
	  {
	    FIND_EQCOM_v2(com[2],0);
	    FIND_EQCOM_v2(com[3],1);

	    arr1[0] = com[3][0]-com[2][0];
	    arr1[1] = com[3][1]-com[2][1];
	    arr1[2] = com[3][2]-com[2][2];

	    /*determine initial euler angles for the e's*/
	    double *e0 = xz_angles(arr1,false);

	    rp->rotation[ox][0] = omega[0] = e0[0];
	    rp->rotation[ox][1] = omega[1] = e0[1];

	    rp->t0_magnitude[ox][0] = mag(arr1,3);

	    /*housekeeping*/
	    free(e0);
	  }

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
		FIND_COM(com[0],i+1);
		FIND_EQCOM(com[1],i+1);
		FIND_EQCOM_v2(com[3],i+1);

		/*determine current e's (features 0-2)*/
		{
		  mp[j+0] = (com[3][0]-com[2][0])/rp->t0_magnitude[ox][0];
		  mp[j+1] = (com[3][1]-com[2][1])/rp->t0_magnitude[ox][0];
		  mp[j+2] = (com[3][2]-com[2][2])/rp->t0_magnitude[ox][0];
	   
		   bank(&(mp[j]),omega[0]);
		    yaw(&(mp[j]),omega[1]);
		}

		/*determine current i's (features 6-9)*/
		{
		  *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[0][0]);
		  *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[0][1]);
		  *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[0][2]);
		  
		  *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[0][0]);
		  *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[0][1]);
		  *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[0][2]);
		  
		  *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[0][0]);
		  *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[0][1]);
		  *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[0][2]);
		  
		  arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

		  mp[j+6] = mag(c1_bp,3)-mag(c0_bp,3);
		  mp[j+7] = arr3[1] - arr3[0];
		  mp[j+8] = mag(c1_bs,3)-mag(c0_bs,3);
		  mp[j+9] = mag(c1_bo,3)-mag(c0_bo,3);
		}

		/*determine current w's (features 3-5)*/
		{
		  ibx[0] = (double)(hr[ox][i+1][3]-com[1][0]);
		  ibx[1] = (double)(hr[ox][i+1][4]-com[1][1]);
		  ibx[2] = (double)(hr[ox][i+1][5]-com[1][2]);

		  free(ibz);
		  free(iby);

		  ibz = cross(ibx,c1_bs);
		  iby = cross(ibx,ibz);

		  mx = mag(ibx,3);
		  my = mag(iby,3);
		  mz = mag(ibz,3);

		  for (int d=0;d<3;d++){
		    ibx[d] /= mx;
		    iby[d] /= my;
		    ibz[d] /= -mz;
		  }

		  mp[j+3] = xzx_delta(iby,ibz,e1,false);
		  mp[j+4] = xzx_delta(ibz,iby,e2,false);
		  mp[j+5] = xzx_delta(ibx,iby,e3,false);

		  free(e1);
		  free(e2);
		  free(e3);

		  e1 = xzx_angles(iby,ibz,false,false); //pitch
		  e2 = xzx_angles(ibz,iby,false,false); //yaw
		  e3 = xzx_angles(ibx,iby,false,false); //bank
		}

	       /*copy data for the next iteration*/
	       {
		 arr3[0] = arr3[1]; 			//for i's
		 for (int q=0;q<3;q++){
		   *(c0_bp+q) = *(c1_bp+q);		//for i's
		   *(c0_bs+q) = *(c1_bs+q);		//for i's
		   *(c0_bo+q) = *(c1_bo+q);		//for i's

		   com[2][q] = com[3][q];		//for e's
		 }
	       }

	   END_iew_MAP


	BEGIN_MAP(iew_v2)

	  double ibx[3], *iby, *ibz;
	  double *e1, *e2, *e3;
	  double  mx ,my, mz;

	  FIND_EQCOM_v2(com[0],0);
	  FIND_EQCOM_v2(com[1],1);

	  /*determine initial structure for the i's*/
	  {

	    *(c0_bo+0) = (double)(hr[ox][0][3]-com[0][0]);
	    *(c0_bo+1) = (double)(hr[ox][0][4]-com[0][1]);
	    *(c0_bo+2) = (double)(hr[ox][0][5]-com[0][2]);
	    
	    *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	    *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	    *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);
	    
	    *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	    *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	    *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);
	    
	    arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));
	  }

	  /*determine initial rodrigues' angles for the w's*/
	  {
	    ibx[0] = (double)(hr[ox][0][3]-com[0][0]);
	    ibx[1] = (double)(hr[ox][0][4]-com[0][1]);
	    ibx[2] = (double)(hr[ox][0][5]-com[0][2]);

	    ibz = cross(ibx,c0_bs);
	    iby = cross(ibx,ibz);
	    
	    mx = mag(ibx,3);
	    my = mag(iby,3);
	    mz = mag(ibz,3);
	    
	    for (int i=0;i<3;i++){
	      ibx[i] /= mx;
	      iby[i] /= my;
	      ibz[i] /= -mz;
	    }
	    
	    e1 = xzx_angles(iby,ibz,false,false); //pitch
	    e2 = xzx_angles(ibz,iby,false,false); //yaw
	    e3 = xzx_angles(ibx,iby,false,false); //bank
	  }

	  /*get the t=0 magnitude of e for normalization*/
	  {

	    arr1[0] = com[1][0]-com[0][0];
	    arr1[1] = com[1][1]-com[0][1];
	    arr1[2] = com[1][2]-com[0][2];

	    /*determine initial euler angles for the e's*/
	    double *e0 = xz_angles(arr1,false);

	    rp->rotation[ox][0] = omega[0] = e0[0];
	    rp->rotation[ox][1] = omega[1] = e0[1];

	    rp->t0_magnitude[ox][0] = mag(arr1,3);

	    /*housekeeping*/
	    free(e0);
	  }

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
		FIND_EQCOM_v2(com[1],i+1);

		/*determine current e's (features 0-2)*/
		{
		  mp[j+0] = (com[1][0]-com[0][0])/rp->t0_magnitude[ox][0];
		  mp[j+1] = (com[1][1]-com[0][1])/rp->t0_magnitude[ox][0];
		  mp[j+2] = (com[1][2]-com[0][2])/rp->t0_magnitude[ox][0];
	   
		   bank(&(mp[j]),omega[0]);
		    yaw(&(mp[j]),omega[1]);
		}

		/*determine current i's (features 6-9)*/
		{
		  *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[1][0]);
		  *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[1][1]);
		  *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[1][2]);
		  
		  *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[1][0]);
		  *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[1][1]);
		  *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[1][2]);
		  
		  *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[1][0]);
		  *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[1][1]);
		  *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[1][2]);
		  
		  arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

		  mp[j+6] = mag(c1_bp,3)-mag(c0_bp,3);
		  mp[j+7] = arr3[1] - arr3[0];
		  mp[j+8] = mag(c1_bs,3)-mag(c0_bs,3);
		  mp[j+9] = mag(c1_bo,3)-mag(c0_bo,3);
		}

		/*determine current w's (features 3-5)*/
		{
		  ibx[0] = (double)(hr[ox][i+1][3]-com[1][0]);
		  ibx[1] = (double)(hr[ox][i+1][4]-com[1][1]);
		  ibx[2] = (double)(hr[ox][i+1][5]-com[1][2]);

		  free(ibz);
		  free(iby);

		  ibz = cross(ibx,c1_bs);
		  iby = cross(ibx,ibz);

		  mx = mag(ibx,3);
		  my = mag(iby,3);
		  mz = mag(ibz,3);

		  for (int d=0;d<3;d++){
		    ibx[d] /= mx;
		    iby[d] /= my;
		    ibz[d] /= -mz;
		  }

		  mp[j+3] = xzx_delta(iby,ibz,e1,false);
		  mp[j+4] = xzx_delta(ibz,iby,e2,false);
		  mp[j+5] = xzx_delta(ibx,iby,e3,false);

		  free(e1);
		  free(e2);
		  free(e3);

		  e1 = xzx_angles(iby,ibz,false,false); //pitch
		  e2 = xzx_angles(ibz,iby,false,false); //yaw
		  e3 = xzx_angles(ibx,iby,false,false); //bank
		}

	       /*copy data for the next iteration*/
	       {
		 arr3[0] = arr3[1]; 			//for i's
		 for (int q=0;q<3;q++){
		   *(c0_bp+q) = *(c1_bp+q);		//for i's
		   *(c0_bs+q) = *(c1_bs+q);		//for i's
		   *(c0_bo+q) = *(c1_bo+q);		//for i's

		   com[0][q] = com[1][q];		//for e's
		 }
	       }

	   END_iew_MAP


	BEGIN_MAP(iew_v3)

	  double ibx[3], *iby, *ibz;
	  double *e1, *e2, *e3, *e0;
	  double  mx ,my, mz;

	  FIND_EQCOM_v2(com[0],0);
	  FIND_EQCOM_v2(com[1],1);

	  /*determine initial structure for the i's*/
	  {
	    *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	    *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	    *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);
	    
	    *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	    *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	    *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);
	       
	    arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));
		  
	    *(c1_bs+0) = (double)(hr[hs][1][3]-com[1][0]);
	    *(c1_bs+1) = (double)(hr[hs][1][4]-com[1][1]);
	    *(c1_bs+2) = (double)(hr[hs][1][5]-com[1][2]);
	  }

	  /*determine initial rodrigues' angles for the w's*/
	  {
	    ibx[0] = (double)(hr[ox][0][3]-com[0][0]);
	    ibx[1] = (double)(hr[ox][0][4]-com[0][1]);
	    ibx[2] = (double)(hr[ox][0][5]-com[0][2]);

	    ibz = cross(ibx,c0_bs);
	    iby = cross(ibx,ibz);
	    
	    mx = mag(ibx,3);
	    my = mag(iby,3);
	    mz = mag(ibz,3);
	    
	    for (int i=0;i<3;i++){
	      ibx[i] /= mx;
	      iby[i] /= my;
	      ibz[i] /= -mz;
	    }
	    
	    e1 = xzx_angles(iby,ibz,false,false); //pitch
	    e2 = xzx_angles(ibz,iby,false,false); //yaw
	    e3 = xzx_angles(ibx,iby,false,false); //bank

	    ibx[0] = (double)(hr[ox][1][3]-com[1][0]);
	    ibx[1] = (double)(hr[ox][1][4]-com[1][1]);
	    ibx[2] = (double)(hr[ox][1][5]-com[1][2]);

	    free(ibz);
	    free(iby);

	    ibz = cross(ibx,c1_bs);
	    iby = cross(ibx,ibz);

	    mx = mag(ibx,3);
	    my = mag(iby,3);
	    mz = mag(ibz,3);

	    for (int d=0;d<3;d++){
	      ibx[d] /= mx;
	      iby[d] /= my;
	      ibz[d] /= -mz;
	    }

	    arr1[0] = xzx_delta(iby,ibz,e1,false);
	    arr1[1] = xzx_delta(ibz,iby,e2,false);
	    arr1[2] = xzx_delta(ibx,iby,e3,false);

	    rp->t0_magnitude[ox][1] = mag(arr1,3);
	  }

	  /*get the t=0 magnitude of e for normalization*/
	  {

	    arr1[0] = com[1][0]-com[0][0];
	    arr1[1] = com[1][1]-com[0][1];
	    arr1[2] = com[1][2]-com[0][2];

	    /*determine initial euler angles for the e's*/
	    e0 = xz_angles(arr1,false);

	    rp->rotation[ox][0] = e0[0];
	    rp->rotation[ox][1] = e0[1];

	    rp->t0_magnitude[ox][0] = mag(arr1,3);

	    /*housekeeping*/
	    free(e0);
	  }

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
		FIND_EQCOM_v2(com[1],i+1);

		/*determine current e's (features 0-2)*/
		{
		  mp[j+0] = (com[1][0]-com[0][0])/rp->t0_magnitude[ox][0];
		  mp[j+1] = (com[1][1]-com[0][1])/rp->t0_magnitude[ox][0];
		  mp[j+2] = (com[1][2]-com[0][2])/rp->t0_magnitude[ox][0];
	   
		   bank(&(mp[j]),rp->rotation[ox][0]);
		    yaw(&(mp[j]),rp->rotation[ox][1]);
		}

		/*determine current i's (features 6-9)*/
		{
		  *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[1][0]);
		  *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[1][1]);
		  *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[1][2]);
		  
		  *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[1][0]);
		  *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[1][1]);
		  *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[1][2]);
		  
		  arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

		  mp[j+6] = (mag(c1_bp,3)-mag(c0_bp,3));
		  mp[j+7] = (arr3[1] - arr3[0]);
		  mp[j+8] = (mag(c1_bs,3)-mag(c0_bs,3));
		}

		/*determine current w's (features 3-5)*/
		{
		  ibx[0] = (double)(hr[ox][i+1][3]-com[1][0]);
		  ibx[1] = (double)(hr[ox][i+1][4]-com[1][1]);
		  ibx[2] = (double)(hr[ox][i+1][5]-com[1][2]);

		  free(ibz);
		  free(iby);

		  ibz = cross(ibx,c1_bs);
		  iby = cross(ibx,ibz);

		  mx = mag(ibx,3);
		  my = mag(iby,3);
		  mz = mag(ibz,3);

		  for (int d=0;d<3;d++){
		    ibx[d] /= mx;
		    iby[d] /= my;
		    ibz[d] /= -mz;
		  }

		  mp[j+3] = xzx_delta(iby,ibz,e1,false)/rp->t0_magnitude[ox][1];
		  mp[j+4] = xzx_delta(ibz,iby,e2,false)/rp->t0_magnitude[ox][1];
		  mp[j+5] = xzx_delta(ibx,iby,e3,false)/rp->t0_magnitude[ox][1];
	 
		   bank(&(mp[j+3]),rp->rotation[ox][2]);
		    yaw(&(mp[j+3]),rp->rotation[ox][3]);

		  free(e1);
		  free(e2);
		  free(e3);

		  e1 = xzx_angles(iby,ibz,false,false); //pitch
		  e2 = xzx_angles(ibz,iby,false,false); //yaw
		  e3 = xzx_angles(ibx,iby,false,false); //bank
		}

	       /*copy data for the next iteration*/
	       {
		 arr3[0] = arr3[1]; 			//for i's
		 for (int q=0;q<3;q++){
		   *(c0_bp+q) = *(c1_bp+q);		//for i's
		   *(c0_bs+q) = *(c1_bs+q);		//for i's
		   *(c0_bo+q) = *(c1_bo+q);		//for i's

		   com[0][q] = com[1][q];		//for e's
		 }
	       }

	   END_iew_MAP


	BEGIN_MAP(eiew_v0)

	  /*compute com at t=0*/
	  FIND_EQCOM(com[0],0);

	  double ibx[3], *iby, *ibz, h2[3];

	  ibx[XX] = hr[ox][0][3]-com[0][0];
	  ibx[YY] = hr[ox][0][4]-com[0][1];
	  ibx[ZZ] = hr[ox][0][5]-com[0][2];

	   h2[XX] = hr[hs][0][3]-com[0][0];
	   h2[YY] = hr[hs][0][4]-com[0][1];
	   h2[ZZ] = hr[hs][0][5]-com[0][2];

	  ibz = cross(ibx,h2);
	  iby = cross(ibx,ibz);

	  double mx,my,mz;
	  mx = mag(ibx,3);
	  my = mag(iby,3);
	  mz = mag(ibz,3);

	  for (int i=0;i<3;i++){
	    ibx[i] /= mx;
	    iby[i] /= my;
	    ibz[i] /= -mz;
	  }

	  rp->rotation[pa] = xzx_angles(ibx,iby,false,false); //bank

	  free(ibz);
	  free(iby);

	  LOOP_THROUGH_HISTORICAL_RECORD_AND

		FIND_EQCOM(com[1],i+1);

		//map the external (i.e. com) velocity
		mp[j+0] = com[1][0]-com[0][0];
		mp[j+1] = com[1][1]-com[0][1];
		mp[j+2] = com[1][2]-com[0][2];

	       //rotate the external velocities into the internal frame using t=0 euler angles
		 bank(&(mp[j]),rp->rotation[pa][0]);
		  yaw(&(mp[j]),rp->rotation[pa][1]);
		 bank(&(mp[j]),rp->rotation[pa][2]);

		/*copy com for next iteration*/
		com[0][0] = com[1][0];
		com[0][1] = com[1][1];
		com[0][2] = com[1][2];

	  END_MAP


	BEGIN_MAP(eiew_v1)

	  /*get the t=0 magnitude of e for normalization*/
	  FIND_EQCOM_v2(com[0],0);
	  FIND_EQCOM_v2(com[1],1);

	  arr1[0] = com[1][0]-com[0][0];
	  arr1[1] = com[1][1]-com[0][1];
	  arr1[2] = com[1][2]-com[0][2];

	  /*determine initial euler angles for the e's*/
	  double *e0 = xz_angles(arr1,false);

	  rp->rotation[ox][0] = e0[0];
	  rp->rotation[ox][1] = e0[1];

	  rp->t0_magnitude[ox][0] = mag(arr1,3);

	  /*housekeeping*/
	  free(e0);

	  LOOP_THROUGH_HISTORICAL_RECORD_AND

		FIND_EQCOM_v2(com[1],i+1);

		mp[j+0] = (com[1][0]-com[0][0])/rp->t0_magnitude[ox][0];
		mp[j+1] = (com[1][1]-com[0][1])/rp->t0_magnitude[ox][0];
		mp[j+2] = (com[1][2]-com[0][2])/rp->t0_magnitude[ox][0];
	   
		   bank(&(mp[j]),rp->rotation[ox][0]);
		    yaw(&(mp[j]),rp->rotation[ox][1]);

		/*copy com for next iteration*/
		com[0][0] = com[1][0];
		com[0][1] = com[1][1];
		com[0][2] = com[1][2];

	  END_MAP


	BEGIN_MAP(iiew_v0)

	  /*correct for the number of features here*/
	  n_feats = 4;
	  mp = (double*)realloc(mp,sizeof*mp*(rp->input_shape[STEPS]*n_feats));

	  FIND_COM(com[0],0);
	  
	  /*capture ox->com at t=0*/
	  *(c0_bo+0) = (double)(hr[ox][0][3]-com[0][0]);
	  *(c0_bo+1) = (double)(hr[ox][0][4]-com[0][1]);
	  *(c0_bo+2) = (double)(hr[ox][0][5]-com[0][2]);

	  /*capture com->h1 at t=0*/
	  *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	  *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	  *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);

	  /*capture com->h2 at t=0*/
	  *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	  *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	  *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);

	  /*get h-com-h angle at t=0*/
	  arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
	       FIND_COM(com[1],i+1);
	       
	       *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[1][0]);
	       *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[1][1]);
	       *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[1][2]);
	       
	       *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[1][0]);
	       *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[1][1]);
	       *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[1][2]);
	       
	       *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[1][0]);
	       *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[1][1]);
	       *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[1][2]);
	       
	       arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

	       //map the internal velocities
	       mp[j+0] = mag(c1_bp,3)-mag(c0_bp,3);
	       mp[j+1] = arr3[1] - arr3[0];
	       mp[j+2] = mag(c1_bs,3)-mag(c0_bs,3);
	       mp[j+3] = mag(c1_bo,3)-mag(c0_bo,3);

	       /*copy the h-com-h angle for next iteration*/
	       arr3[0] = arr3[1];

	       /*copy c1s into c0s for next iteration*/
	       for (int q=0;q<3;q++){
		 *(c0_bp+q) = *(c1_bp+q);
		 *(c0_bs+q) = *(c1_bs+q);
		 *(c0_bo+q) = *(c1_bo+q);
	       }

	  END_MAP


	BEGIN_MAP(iiew_v2)

	  /*correct for the number of features here*/
	  n_feats = 4;
	  mp = (double*)realloc(mp,sizeof*mp*(rp->input_shape[STEPS]*n_feats));

	  FIND_EQCOM_v2(com[0],0);
	  
	  /*capture ox->com at t=0*/
	  *(c0_bo+0) = (double)(hr[ox][0][3]-com[0][0]);
	  *(c0_bo+1) = (double)(hr[ox][0][4]-com[0][1]);
	  *(c0_bo+2) = (double)(hr[ox][0][5]-com[0][2]);

	  /*capture com->h1 at t=0*/
	  *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	  *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	  *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);

	  /*capture com->h2 at t=0*/
	  *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	  *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	  *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);

	  /*get h-com-h angle at t=0*/
	  arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
	       FIND_EQCOM_v2(com[1],i+1);
	       
	       *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[1][0]);
	       *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[1][1]);
	       *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[1][2]);
	       
	       *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[1][0]);
	       *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[1][1]);
	       *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[1][2]);
	       
	       *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[1][0]);
	       *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[1][1]);
	       *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[1][2]);
	       
	       arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

	       //map the internal velocities
	       mp[j+0] = mag(c1_bp,3)-mag(c0_bp,3);
	       mp[j+1] = arr3[1] - arr3[0];
	       mp[j+2] = mag(c1_bs,3)-mag(c0_bs,3);
	       mp[j+3] = mag(c1_bo,3)-mag(c0_bo,3);

	       /*copy the h-com-h angle for next iteration*/
	       arr3[0] = arr3[1];

	       /*copy c1s into c0s for next iteration*/
	       for (int q=0;q<3;q++){
		 *(c0_bp+q) = *(c1_bp+q);
		 *(c0_bs+q) = *(c1_bs+q);
		 *(c0_bo+q) = *(c1_bo+q);
	       }

	  END_MAP


	BEGIN_MAP(iiew_v3)

	  /*get fixed EQ COM at t=0 and t=1*/
	  FIND_EQCOM_v2(com[0],0);
	  //FIND_EQCOM_v2(com[1],1);

	  /*get the bond vectors at t=0 and t=1*/
	  *(c0_bp+0) = (double)(hr[hp][0][3]-com[0][0]);
	  *(c0_bp+1) = (double)(hr[hp][0][4]-com[0][1]);
	  *(c0_bp+2) = (double)(hr[hp][0][5]-com[0][2]);
	  
	  *(c0_bs+0) = (double)(hr[hs][0][3]-com[0][0]);
	  *(c0_bs+1) = (double)(hr[hs][0][4]-com[0][1]);
	  *(c0_bs+2) = (double)(hr[hs][0][5]-com[0][2]);
	  
	  //*(c1_bp+0) = (double)(hr[hp][1][3]-com[1][0]);
	  //*(c1_bp+1) = (double)(hr[hp][1][4]-com[1][1]);
	  //*(c1_bp+2) = (double)(hr[hp][1][5]-com[1][2]);
	  
	  //*(c1_bs+0) = (double)(hr[hs][1][3]-com[1][0]);
	  //*(c1_bs+1) = (double)(hr[hs][1][4]-com[1][1]);
	  //*(c1_bs+2) = (double)(hr[hs][1][5]-com[1][2]);
	  
	  arr3[0] = acos(dot(c0_bp,c0_bs,3)/(mag(c0_bp,3)*mag(c0_bs,3)));
	  //arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

	  //arr1[0] = mag(c1_bp,3)-mag(c0_bp,3);
	  //arr1[1] = arr3[1] - arr3[0];
	  //arr1[2] = mag(c1_bs,3)-mag(c0_bs,3);

	  //rp->t0_magnitude[ox][2] = mag(arr1,3);

	  LOOP_THROUGH_HISTORICAL_RECORD_AND
		
	       FIND_EQCOM_v2(com[1],i+1);
	       
	       *(c1_bo+0) = (double)(hr[ox][i+1][3]-com[1][0]);
	       *(c1_bo+1) = (double)(hr[ox][i+1][4]-com[1][1]);
	       *(c1_bo+2) = (double)(hr[ox][i+1][5]-com[1][2]);
	       
	       *(c1_bp+0) = (double)(hr[hp][i+1][3]-com[1][0]);
	       *(c1_bp+1) = (double)(hr[hp][i+1][4]-com[1][1]);
	       *(c1_bp+2) = (double)(hr[hp][i+1][5]-com[1][2]);
	       
	       *(c1_bs+0) = (double)(hr[hs][i+1][3]-com[1][0]);
	       *(c1_bs+1) = (double)(hr[hs][i+1][4]-com[1][1]);
	       *(c1_bs+2) = (double)(hr[hs][i+1][5]-com[1][2]);
	       
	       arr3[1] = acos(dot(c1_bp,c1_bs,3)/(mag(c1_bp,3)*mag(c1_bs,3)));

	       //map the internal velocities
	       mp[j+0] = mag(c1_bp,3)-mag(c0_bp,3);
	       mp[j+1] = arr3[1] - arr3[0];
	       mp[j+2] = mag(c1_bs,3)-mag(c0_bs,3);

	       /*copy the h-com-h angle for next iteration*/
	       arr3[0] = arr3[1];

	       /*copy c1s into c0s for next iteration*/
	       for (int q=0;q<3;q++){
		 *(c0_bp+q) = *(c1_bp+q);
		 *(c0_bs+q) = *(c1_bs+q);
	       }

	  END_MAP


	BEGIN_MAP(wiew_v0)

	  FIND_EQCOM(com[0],0);

	  double ibx[3], *iby, *ibz, h2[3];

	  ibx[XX] = hr[ox][0][3]-com[0][0];
	  ibx[YY] = hr[ox][0][4]-com[0][1];
	  ibx[ZZ] = hr[ox][0][5]-com[0][2];

	   h2[XX] = hr[hs][0][3]-com[0][0];
	   h2[YY] = hr[hs][0][4]-com[0][1];
	   h2[ZZ] = hr[hs][0][5]-com[0][2];

	  ibz = cross(ibx,h2);
	  iby = cross(ibx,ibz);

	  double mx,my,mz;
	  mx = mag(ibx,3);
	  my = mag(iby,3);
	  mz = mag(ibz,3);

	  for (int i=0;i<3;i++){
	    ibx[i] /= mx;
	    iby[i] /= my;
	    ibz[i] /= -mz;
	  }

	  double *e1 = xzx_angles(iby,ibz,false,false); //pitch
	  double *e2 = xzx_angles(ibz,iby,false,false); //yaw
	  double *e3 = xzx_angles(ibx,iby,false,false); //bank

	  LOOP_THROUGH_HISTORICAL_RECORD_AND

		FIND_EQCOM(com[1],i+1);

		ibx[XX] = hr[ox][i+1][3]-com[1][0];
		ibx[YY] = hr[ox][i+1][4]-com[1][1];
		ibx[ZZ] = hr[ox][i+1][5]-com[1][2];
		
		 h2[XX] = hr[hs][i+1][3]-com[1][0];
		 h2[YY] = hr[hs][i+1][4]-com[1][1];
		 h2[ZZ] = hr[hs][i+1][5]-com[1][2];
		
		ibz = cross(ibx,h2);
		iby = cross(ibx,ibz);
		
		mx = mag(ibx,3);
		my = mag(iby,3);
		mz = mag(ibz,3);
		
		for (int d=0;d<3;d++){
		  ibx[d] /= mx;
		  iby[d] /= my;
		  ibz[d] /= -mz;
		}

		mp[j+0] = xzx_delta(iby,ibz,e1,false);
		mp[j+1] = xzx_delta(ibz,iby,e2,false);
		mp[j+2] = xzx_delta(ibx,iby,e3,false);

		free(e1);
		free(e2);
		free(e3);

		e1 = xzx_angles(iby,ibz,false,false); //pitch
		e2 = xzx_angles(ibz,iby,false,false); //yaw
		e3 = xzx_angles(ibx,iby,false,false); //bank

	  END_MAP


	BEGIN_MAP(wiew_v2)

	  FIND_EQCOM_v2(com[0],0);

	  double ibx[3], *iby, *ibz, h2[3];

	  ibx[XX] = hr[ox][0][3]-com[0][0];
	  ibx[YY] = hr[ox][0][4]-com[0][1];
	  ibx[ZZ] = hr[ox][0][5]-com[0][2];

	   h2[XX] = hr[hs][0][3]-com[0][0];
	   h2[YY] = hr[hs][0][4]-com[0][1];
	   h2[ZZ] = hr[hs][0][5]-com[0][2];

	  ibz = cross(ibx,h2);
	  iby = cross(ibx,ibz);

	  double mx,my,mz;
	  mx = mag(ibx,3);
	  my = mag(iby,3);
	  mz = mag(ibz,3);

	  for (int i=0;i<3;i++){
	    ibx[i] /= mx;
	    iby[i] /= my;
	    ibz[i] /= -mz;
	  }

	  double *e1 = xzx_angles(iby,ibz,false,false); //pitch
	  double *e2 = xzx_angles(ibz,iby,false,false); //yaw
	  double *e3 = xzx_angles(ibx,iby,false,false); //bank

	  LOOP_THROUGH_HISTORICAL_RECORD_AND

		FIND_EQCOM_v2(com[1],i+1);

		ibx[XX] = hr[ox][i+1][3]-com[1][0];
		ibx[YY] = hr[ox][i+1][4]-com[1][1];
		ibx[ZZ] = hr[ox][i+1][5]-com[1][2];
		
		 h2[XX] = hr[hs][i+1][3]-com[1][0];
		 h2[YY] = hr[hs][i+1][4]-com[1][1];
		 h2[ZZ] = hr[hs][i+1][5]-com[1][2];
		
		ibz = cross(ibx,h2);
		iby = cross(ibx,ibz);
		
		mx = mag(ibx,3);
		my = mag(iby,3);
		mz = mag(ibz,3);
		
		for (int d=0;d<3;d++){
		  ibx[d] /= mx;
		  iby[d] /= my;
		  ibz[d] /= -mz;
		}

		mp[j+0] = xzx_delta(iby,ibz,e1,false);
		mp[j+1] = xzx_delta(ibz,iby,e2,false);
		mp[j+2] = xzx_delta(ibx,iby,e3,false);

		free(e1);
		free(e2);
		free(e3);

		e1 = xzx_angles(iby,ibz,false,false); //pitch
		e2 = xzx_angles(ibz,iby,false,false); //yaw
		e3 = xzx_angles(ibx,iby,false,false); //bank

	  END_MAP


	BEGIN_MAP(wiew_v3)

	  FIND_EQCOM_v2(com[0],0);

	  double ibx[3], *iby, *ibz, h2[3];

	  ibx[XX] = hr[ox][0][3]-com[0][0];
	  ibx[YY] = hr[ox][0][4]-com[0][1];
	  ibx[ZZ] = hr[ox][0][5]-com[0][2];

	   h2[XX] = hr[hs][0][3]-com[0][0];
	   h2[YY] = hr[hs][0][4]-com[0][1];
	   h2[ZZ] = hr[hs][0][5]-com[0][2];

	  ibz = cross(ibx,h2);
	  iby = cross(ibx,ibz);

	  double mx,my,mz;
	  mx = mag(ibx,3);
	  my = mag(iby,3);
	  mz = mag(ibz,3);

	  for (int i=0;i<3;i++){
	    ibx[i] /= mx;
	    iby[i] /= my;
	    ibz[i] /= -mz;
	  }

	  double *e1 = xzx_angles(iby,ibz,false,false); //pitch
	  double *e2 = xzx_angles(ibz,iby,false,false); //yaw
	  double *e3 = xzx_angles(ibx,iby,false,false); //bank

	  FIND_EQCOM_v2(com[1],1);
	  
	  ibx[XX] = hr[ox][1][3]-com[1][0];
	  ibx[YY] = hr[ox][1][4]-com[1][1];
	  ibx[ZZ] = hr[ox][1][5]-com[1][2];
	  
	   h2[XX] = hr[hs][1][3]-com[1][0];
	   h2[YY] = hr[hs][1][4]-com[1][1];
	   h2[ZZ] = hr[hs][1][5]-com[1][2];
	  
	  ibz = cross(ibx,h2);
	  iby = cross(ibx,ibz);
	  
	  mx = mag(ibx,3);
	  my = mag(iby,3);
	  mz = mag(ibz,3);
	  
	  for (int d=0;d<3;d++){
	    ibx[d] /= mx;
	    iby[d] /= my;
	    ibz[d] /= -mz;
	  }
	  
	  arr1[0] = xzx_delta(iby,ibz,e1,false);
	  arr1[1] = xzx_delta(ibz,iby,e2,false);
	  arr1[2] = xzx_delta(ibx,iby,e3,false);

	  /*determine initial euler angles for the e's*/
	  double *e0 = xz_angles(arr1,false);

	  rp->rotation[ox][2] = e0[0];
	  rp->rotation[ox][3] = e0[1];

	  free(e0);

	  rp->t0_magnitude[ox][1] = mag(arr1,3);

	  LOOP_THROUGH_HISTORICAL_RECORD_AND

		FIND_EQCOM_v2(com[1],i+1);

		ibx[XX] = hr[ox][i+1][3]-com[1][0];
		ibx[YY] = hr[ox][i+1][4]-com[1][1];
		ibx[ZZ] = hr[ox][i+1][5]-com[1][2];
		
		 h2[XX] = hr[hs][i+1][3]-com[1][0];
		 h2[YY] = hr[hs][i+1][4]-com[1][1];
		 h2[ZZ] = hr[hs][i+1][5]-com[1][2];
		
		ibz = cross(ibx,h2);
		iby = cross(ibx,ibz);
		
		mx = mag(ibx,3);
		my = mag(iby,3);
		mz = mag(ibz,3);
		
		for (int d=0;d<3;d++){
		  ibx[d] /= mx;
		  iby[d] /= my;
		  ibz[d] /= -mz;
		}

		mp[j+0] = xzx_delta(iby,ibz,e1,false)/rp->t0_magnitude[ox][1];
		mp[j+1] = xzx_delta(ibz,iby,e2,false)/rp->t0_magnitude[ox][1];
		mp[j+2] = xzx_delta(ibx,iby,e3,false)/rp->t0_magnitude[ox][1];
	 
		bank(&(mp[j]),rp->rotation[ox][2]);
		 yaw(&(mp[j]),rp->rotation[ox][3]);

		free(e1);
		free(e2);
		free(e3);

		e1 = xzx_angles(iby,ibz,false,false); //pitch
		e2 = xzx_angles(ibz,iby,false,false); //yaw
		e3 = xzx_angles(ibx,iby,false,false); //bank

	  END_MAP


	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				    end of mapping macros
	\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


	#undef FIND_EQCOM
	#undef FIND_EQCOM_v2
	#undef FIND_COM
	#undef PRIME_UELOVITIES_FOR
	#undef RENDER_UELOCITIES_FOR
	#undef RENDER_VELOCITIES_FOR
	#undef BEGIN_MAP
	#undef LOOP_THROUGH_HISTORICAL_RECORD_AND
	#undef END_iew_MAP
	#undef END_MAP


	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
				      feature set table
	\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


	#define FSET(_A_,_B_,_C_,_D_,_E_,_F_,_G_,_H_,_I_)\
	 { #_A_,_B_,_C_,_D_,_E_,_F_,_G_,&rp_##_A_##_map,&rp_##_H_##_pam,_I_,DBL}


	#define NULL_FSET(_A_)\
	 { #_A_,0,0,0,0,0,0,NULL,NULL,NULL,DBL }


	const struct rp_fset rp_fset_table[] = {

	//FSET(name,#nad,#in,#out,bloat,cap_pos?,whol_mol?,pam,vfg)

	  /*general purpose, single-atom descriptors*/ 
	  FSET(3v,	 3, 3, 3, 0, false, false,    v, &rp_og_vfg ),
	  //FSET(3vR0,	 3, 3, 3, 0, false, false,   vR, &rp_og_vfg ), //before we included the back mapping constants for test set eval
	  FSET(3vR0,	 5, 3, 3, 0, false, false,   vR, &rp_og_vfg ),
	  FSET(ns3vRv0,	 3, 3, 3, 0, false, false, NSvR, &rp_og_vfg ),

	  FSET(3u,	 3, 3, 3, 0, false, false,    u, &rp_og_vfg ),
	  //FSET(3uR0,	 3, 3, 3, 0, false, false,   uR, &rp_og_vfg ), //before we included the back mapping constants for test set eval
	  FSET(3uR0,	 6, 3, 3, 0, false, false,   uR, &rp_og_vfg ), //before we included the back mapping constants for test set eval
	  FSET(3uCoMRv0, 3, 3, 3, 1,  true, false,   uR, &rp_og_vfg ),
	  FSET(D3v, 	 9, 6, 6, 1,  true, false,   vC, &rp_og_vfg ),
	  FSET(D3uRv0, 	 9, 6, 6, 1,  true, false,  uRC, &rp_og_vfg ),
	  FSET(D9uRv0, 	 9, 9, 9, 1,  true, false, 9uRC, &rp_og_vfg ),
	  //FSET(DuR0, 	 9, 6, 6, 1,  true, false, DuR0, &rp_og_vfg ), //before we included the back mapping constants for test set eval
	  FSET(DuR0, 	 11, 6, 6, 1,  true, false, DuR0, &rp_og_vfg ), //after ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	  FSET(D9uRt0,   9, 9, 9, 1,  true, false,9uRCS, &rp_og_vfg ),  
	  FSET(D9uRb0,   9, 9, 9, 1,  true, false,9uRCS, &rp_og_vfg ),  
	  FSET(C9uRt0, 	 9, 9, 9, 1,  true, false, 9uRC, &rp_og_vfg ),
	  FSET(C9dRv0, 	 9, 9, 9, 1,  true, false, 9uRC, &rp_og_vfg ),
	  FSET(D9Rv0, 	 10, 10, 10, 1,  true, false, 9uRC2,&rp_og_vfg ),
	  FSET(D12Rv0, 	 12,12,9, 1,  true, false, 9uRC, &rp_og_vfg ),
	  FSET(D6uRv0, 	 9, 6, 6, 1,  true, false, 6uRC, &rp_og_vfg ),
	  FSET(D3uSRv0,	 9, 6, 6, 1,  true, false, uRCS, &rp_og_vfg ),
	  FSET(ns3uRv0,	 3, 3, 3, 0, false, false, NSuR, &rp_og_vfg ),

	  /*water specific, whole-molecule descriptors*/
	  FSET(iew_v0,	10, 3, 3, 1,  true,  true, NULL, &rp_v0_vfg ),
	  FSET(iew_v1,	10, 3, 4, 1,  true,  true, NULL, &rp_v1_vfg ),
	  FSET(iew_v2,	10, 3, 4, 1,  true,  true, NULL, &rp_v2_vfg ),
	  FSET(iew_v3,	 9, 9, 3, 1,  true,  true, NULL, &rp_v3_vfg ),

	  /*end of the line*/
	  NULL_FSET(nullptr)
	};


	#undef NULLFSET
	#undef FSET


	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
					  dumping
	\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


	//scrubbed 5-19-2022 -Wayyne
	FILE **rp_press_check(gmx_int64_t *step){

	    int *p = rp->pattern, cond = 0;

	    //check the current condition of the nonactive data dump
	    cond += !p[VALI_DELAY] && !p[TEST_DELAY];
	    cond += *step >= p[VALI_DELAY] && p[VALI_DELAY] && !p[TEST_DELAY];
	    cond += *step >= p[TEST_DELAY] && p[TEST_DELAY];

	    if (cond){
	      //create a file pointer for each group
	      char buf[10000]; 
	      FILE **fp = (FILE **)malloc(sizeof*fp*rp->Ngr); 
	      const char *type = (p[VALI_DELAY] == p[TEST_DELAY]) ? 
		  "training": ((p[TEST_DELAY] == 0) ? "validation":"testing");
	      const char *mode = *step == p[DELAY] || *step == p[VALI_DELAY] || 
		  *step == p[TEST_DELAY] ? "w":"a";

	      for (int i=0;i<rp->Ngr;i++){
		sprintf(buf,"rp_nonactive_dump/%d.%s.%s.dat",
		  i,rp_fset_table[rp->fset].name,type);
		fp[i] = fopen(buf,mode);
	      }
	      return fp;
	    }
	    else
	      return NULL;
	}


	//scrubbed 5-19-2022 -Wayyne
	void rp_wipe_n_check(FILE **fp, gmx_int64_t *step){
	    
	    int *p = rp->pattern;

	    //loop over groups 
	    for (int g=0;g<rp->Ngr;g++){

	      //do some housekeeping
	      fclose(fp[g]);

	      //make sure we get the number of samples we want
	      if ((*step+1) == (p[DELAY] + p[STRIDE]) 
		&& (rp->dumped[g][0] < p[TRUE_STRIDE])){
		  if (rp->threshold)
		    p[STRIDE] += 1;
		  rp->input_rec->nsteps++;
		  break;
	      }
	      else if (p[VALI_DELAY] && (*step+1) == (p[VALI_DELAY] + p[VALI_STRIDE])
		&& (rp->dumped[g][0] < (int)ceil(p[TRUE_STRIDE]*.5))){
		  if (rp->threshold)
		    p[VALI_STRIDE] += 1;
		  rp->input_rec->nsteps++;
		  break;
	      }
	      else if (p[TEST_DELAY] && (*step+1) == (p[TEST_DELAY]+p[TEST_STRIDE])
		&& (rp->dumped[g][0] < (int)ceil(p[TRUE_STRIDE]*.5))){
		  if (rp->threshold) 
		    p[TEST_STRIDE] += 1;
		  rp->input_rec->nsteps++;
		  break;
	      }
	    }    

	    //check the status and update vali/test prn
	    if ((*step+1) == (p[DELAY] + p[STRIDE])){
		printf("\nRobin_Protocol: training dump complete.\n"); 
		for (int g=0;g<rp->Ngr;g++){
		  printf("\tgroup %d: %d/%d\n",g,rp->dumped[g][0],p[TRUE_STRIDE]);
		  rp->dumped[g][0] = 0;
		}
		fflush(stdout);
		p[VALI_DELAY] = *step + p[PAUSE] + 1;
		rp->dump_state = RP_VALIDATION_DUMP;
	    }
	    else if (p[VALI_DELAY] && (*step+1) == (p[VALI_DELAY] + p[VALI_STRIDE])){
		printf("\nRobin_Protocol: validation dump complete.\n"); 
		for (int g=0;g<rp->Ngr;g++){
		  printf("\tgroup %d: %d/%d\n",g,rp->dumped[g][0],(int)ceil(p[TRUE_STRIDE]*.5));
		  rp->dumped[g][0] = 0;
		}
		fflush(stdout);
		p[TEST_DELAY] = *step + p[PAUSE] + 1;
		rp->dump_state = RP_TESTING_DUMP;

		//ensure dump threshold is never applied to testing data
		rp->threshold = 0;
	    }
	    else if (p[TEST_DELAY] && (*step+1) == (p[TEST_DELAY]+p[TEST_STRIDE])){
		printf("\nRobin_Protocol: testing dump complete.\n"); 
		for (int g=0;g<rp->Ngr;g++)
		  printf("\tgroup %d: %d/%d\n",g,rp->dumped[0][0],(int)ceil(p[TRUE_STRIDE]*.5));
		fflush(stdout);
		rp->dump_vels = rp->dump_rvels = FALSE;
		rp->dump_state = RP_FINISHED_DUMP;
		rp->input_rec->nsteps = *step+1;
	    }
	 
	    //housekeeping
	    if (fp)
	      free(fp);
	}


	//there's [not] always room for seconds
	bool rp_room_for_seconds(double **sample, int grp){

	    int i = (rp->input_shape[STEPS]-1)*rp->input_shape[FEATS];

	    double target_x = (*sample)[i+0];
	    double target_y = (*sample)[i+1];
	    double target_z = (*sample)[i+2];

	    target_x *= target_x;
	    target_y *= target_y;
	    target_z *= target_z;
	  
	    double mag = sqrt(target_x+target_y+target_z);

	    int bin;

	    if (mag < 1e-7)
	      bin = 0;
	    else if (mag < 1e-6)
	      bin = 1;
	    else if (mag < 1e-5)
	      bin = 2;
	    else if (mag < 1e-4)
	      bin = 3;
	    else if (mag < 1e-3)
	      bin = 4;
	    else if (mag < 1e-2)
	      bin = 5;
	    else if (mag < 1e-1)
	      bin = 6;
	    else
	      bin = 7;
	      
	    if (rp->dumped[grp][bin] >= rp->threshold){
	      if (rp->pattern[VALI_DELAY] == 0)
		rp->pattern[STRIDE] += 1;
	      else if (rp->pattern[TEST_DELAY] == 0)
		rp->pattern[VALI_STRIDE] += 1;
	      else
		rp->pattern[TEST_STRIDE] += 1;
	      return false;
	    }

	    rp->dumped[grp][bin]++;

	    printf("\nPut one in grp %d's bin %d: mag=% .7e\n",grp,bin,mag);

	    return true;
	}


	//dump data in reverse time order
	void rev_dump(double *data, FILE *fp){

	  int steps = rp->input_shape[STEPS];

	  if (rp->dilation_factor != -1)
	    steps /= rp->dilation_factor;
	    
	  for (int i=steps;i>0;i--)
	    for (int j=0;j<rp->input_shape[FEATS];j++)
		fprintf(fp,"% .7e%s", 
		  j < 9 ? data[((i-1)*rp->input_shape[FEATS])+j]*-1.0 : 
		    data[((i-1)*rp->input_shape[FEATS])+j],
		  (j == rp->input_shape[FEATS]-1) ? "\n" : ",");
	}


	//reverse the historical record 
	void reverse_record(int n, bool recursive){

	    if (recursive){
	      //reverse the historical record for the entire molecule atom n belongs to
	      int mod = n%rp->modulus;
	      int hp = (mod == 0) ? n+1: n;
	      int hs = (mod == 0) ? n+2: ((mod == 1) ? n+1 : n-1);           
	      int ox = (mod == 0) ? n  : ((mod == 1) ? n-1 : n-2);           

	      reverse_record(ox,false);
	      reverse_record(hp,false);
	      reverse_record(hs,false);

	    }
	    else{
	      //alloc and init a temporary duplicate of the record
	      real **tmp = (real**)malloc(sizeof*tmp*rp->histr_shape[STEPS]);
	      for (int j=0;j<rp->histr_shape[STEPS];j++){
		tmp[j] = (real*)malloc(sizeof**tmp*rp->histr_shape[FEATS]);
		for (int f=0;f<rp->histr_shape[FEATS];f++)
		  tmp[j][f] = rp->history[n][j][f];
      }

      //"flip it and reverse it" -- Missy Elliot
      for (int i=rp->histr_shape[STEPS]-1,j=0;i>=0;i--,j++)
          for (int f=0;f<rp->histr_shape[FEATS];f++)
            rp->history[n][j][f] = f < 3 ? -tmp[i][f]:tmp[i][f];
      
      //housekeeping
      for (int j=0;j<rp->histr_shape[STEPS];j++)
        free(tmp[j]);
      free(tmp);
    }
}


void rp_nonactive_dump(gmx_int64_t step){

    FILE **fp;

    //set the dump state and make sure we should be dumping
    if (step < rp->pattern[DELAY]){
      rp->dump_state = RP_HOLDING_DUMP;
      return;
    }
    else if (step == rp->pattern[DELAY])
      rp->dump_state == RP_TRAINING_DUMP;
 
    if ((fp = rp_press_check(&step)) == NULL)
      return;

    //make sure we have a valid set of file pointers
    for (int i=0;i<rp->Ngr;i++)
      if (fp[i] == NULL)
        return;

    int max_data = rp->input_shape[STEPS]*rp->input_shape[FEATS];
    double *data = NULL;

    //write file header p.r.n.
    if (step == rp->pattern[DELAY] 
     || step == rp->pattern[VALI_DELAY] 
     || step == rp->pattern[TEST_DELAY])
        for (int i=0;i<rp->Ngr;i++)
          for (int j=1;j<=rp->input_shape[FEATS];j++)
            fprintf(fp[i],"\"f%d\"%s",j,j==rp->input_shape[FEATS] ? "\n":",");

    if (rp->balanced_diet){
      for (int g,n=0;n<rp->natoms;n++){
        data = (double*)(*rp->map)(n);
        g = rp_get_group(n);

        if (rp_room_for_seconds(&data,g)){
          for (int i=0;i<max_data;i++)
            fprintf(fp[g],"% .7e%s", data[i],
              (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ? 
                "\n" : ",");
        }
      }
    }
    //dilated dumping
    else if (rp->dilation_factor != -1){
      int compressed_max_data = 
        (rp->input_shape[STEPS]/rp->dilation_factor)*rp->input_shape[FEATS];
      double *compressed_data = 
        (double*)malloc(sizeof*compressed_data*compressed_max_data);

      for (int g,n=0;n<rp->natoms;n++){
        
            data = (double*)(*rp->map)(n);
            g = rp_get_group(n);

            //loop over dilated time steps
            int incr = rp->input_shape[FEATS]*rp->dilation_factor;
            for (int i=0,h=0;i<max_data;i+=incr,h++){
              double feats[12] = {0};

              //dilate the first 9 features 
              //TODO: this works for 9d3b,12d,and 9a3d fsets only
              for (int j=0;j<rp->dilation_factor;j++){
                for (int k=0;k<9;k++)
                  feats[k] += data[i+(j*rp->input_shape[FEATS])+k];
              }
        
              int reported_index = 
                i+rp->input_shape[FEATS]*(rp->dilation_factor-1);
              feats[9] = data[reported_index+9];
              feats[10] = data[reported_index+10];
              feats[11] = data[reported_index+11];

              for (int j=0;j<rp->input_shape[FEATS];j++)
                compressed_data[(h*rp->input_shape[FEATS])+j] = feats[j];
            }
  
        if (rp->dump_vels)
          for (int i=0;i<compressed_max_data;i++)
            fprintf(fp[g],"% .7e%s", compressed_data[i],
              (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ? 
                "\n" : ",");

        //reverse dump vels p.r.n.
        if (rp->dump_rvels)
          rev_dump(compressed_data,fp[g]);

            //TODO: the following code is useful to validate time dilation
            //fprintf(fp[g],"\nUNDILATED BELOW\n");
 
            //for (int i=0;i<max_data;i++)
              //fprintf(fp[g],"% .7e%s", data[i],
                //(i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ?
                  // "\n" : ",");

       }
       free(compressed_data); data = NULL;
    }
    //normal unbiased dumping
    else if (rp->dump_state != RP_TESTING_DUMP && rp->threshold == 0)
      for (int g,m,n=0;n<rp->natoms;n+=1){

        g = rp_get_group(n);

	//don't worry about atoms that have already dumped a suffient # of samples
	m = rp->dump_state == RP_VALIDATION_DUMP ? 
		rp->pattern[TRUE_STRIDE]*.5 : rp->pattern[TRUE_STRIDE];	

	//TODO: this is a tweak for dumping diatomic molecules
	m *= 2.0;

        if (rp->dumped[g][0] >= m){
	  continue;
        }

        //TODO: if we dump a sample from each atom at each timestep, then...
        //TODO: did I put this in because I was doing condensed phase?
	//if (n%rp->modulus == 1 && g == 0)
	//if (n%rp->modulus != step%rp->modulus)
	 // continue;

	//forward dump vels p.r.n.
        if (rp->dump_vels){

	  //map the historical record to grayson input
          data = (double*)(*rp->map)(n);

          //write the mapped input to disk
          for (int i=0;i<max_data;i++)
            fprintf(fp[g],"% .7e%s", data[i],
              (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ? 
                "\n" : ",");

	  rp->dumped[g][0]++;
	  //rp->dumped[1][0]++; //hack for dumping hydrogens without changing elsewhere

          //housekeeping
	  if (data != NULL){
		free(data);
		data = NULL;
	  }
	}

        if (rp->dumped[g][0] >= m)
          continue;

	//reverse dump vels p.r.n.
        if (rp->dump_rvels){
	  //recursively reverse the historical record
	  reverse_record(n,false);

	  //map the historical record to grayson input
	  data = (double*)(*rp->map)(n);

          //write the mapped input to disk
	  for (int i=0;i<max_data;i++)
	    fprintf(fp[g],"% .7e%s", data[i],
	      (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ?
		"\n" : ",");

	  rp->dumped[g][0]++;
	  //rp->dumped[1][0]++; //hack for dumping hydrogens without changing elsewhere

	  //recursively restore the historical record
	  reverse_record(n,false);

	  //housekeeping
	  if (data != NULL){
		free(data);
		data = NULL;
	  }
	}
      }
    //threshold dumping
    else if (rp->dump_state != RP_TESTING_DUMP && rp->threshold != 0){
      bool below_threshold; int m, n;
      for (int g,n=0;n<rp->natoms;n++){
        g = rp_get_group(n);

	//don't worry about atoms that have already dumped a suffient # of samples
	m = rp->dump_state == RP_VALIDATION_DUMP ? 
		rp->pattern[TRUE_STRIDE]*.5 : rp->pattern[TRUE_STRIDE];	

        if (rp->dumped[g][0] >= m){
	  continue;
        }

	//dump vels p.r.n.
        if (rp->dump_vels){
          data = (double*)(*rp->map)(n);
          below_threshold = false;

          for (int M = m = max_data-rp->input_shape[FEATS];m<M+rp->threshold_domain;m++)
            if (fabs(data[m]) > rp->threshold){
              below_threshold = true;
              break;
            }

          if (!below_threshold){
            for (int i=0;i<max_data;i++)
              fprintf(fp[g],"% .7e%s", data[i],
                (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ?
                  "\n" : ",");
	    rp->dumped[g][0]++;
	  }

	  if (data != NULL){free(data);data = NULL;}
	}

        if (rp->dumped[g][0] >= m){
	  continue;
        }

        //reverse dump vels p.r.n.
        if (rp->dump_rvels){

	  reverse_record(n,true);

	  data = (double*)(*rp->map)(n);
	  below_threshold = false;

          for (int M = m =max_data-rp->input_shape[FEATS];m<M+rp->threshold_domain;m++)
            if (fabs(data[m]) > rp->threshold){
              below_threshold = true;
              break;
            } 

          if (!below_threshold){
            for (int i=0;i<max_data;i++)
              fprintf(fp[g],"% .7e%s", data[i],
                (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ?
                  "\n" : ",");
	    rp->dumped[g][0]++;
	  }

	  reverse_record(n,true);

	  if (data != NULL){
		free(data);
		data = NULL;
	  }

	}
      }
    }
    //testing set dump (always f only)
    else{ 
      FILE *fp2; char buf[1000]; int m = (int)ceil(rp->pattern[TRUE_STRIDE]*0.5);

      //TODO: hack for diatomic molecules
      m*=2.0;

      for (int g,n=0;n<rp->natoms;n+=1){
	g = rp_get_group(n); 

        if (rp->dumped[g][0] >= m){
	  continue;
        }

//TODO: was this only here for condensed phase?
	//if (n%rp->modulus == 1 && g == 0)
	//if (n%rp->modulus != step%rp->modulus)
	  //continue;

        if (rp->dump_vels || rp->dump_rvels){
	  data = (double*)(*rp->map)(n);
          for (int i=0;i<max_data;i++)
            fprintf(fp[g],"% .7e%s", data[i],
              (i % rp->input_shape[FEATS] == (rp->input_shape[FEATS]-1)) ? 
                "\n" : ",");

	  rp->dumped[g][0]++;
	  //rp->dumped[1][0]++; //hack for dumping hydrogens without changing elsewhere

	  if (data != NULL){
	    free(data);
	    data = NULL;
	  }
	}
      }
    }

    //housekeeping 
    if (data != NULL){
	free(data);
	data = NULL;
    }

    rp_wipe_n_check(fp,&step);    
}

void rp_active_dump(int n,rvec * gmx_restrict v){
      char rp_tmp[10000] = {0};
      sprintf(rp_tmp,"rp_active_dump/%d.inj",n);
      FILE *rp_fp = fopen(rp_tmp,"a");
      fprintf(rp_fp,"inj%d,% .7e,% .7e,% .7e,target,% .7e,% .7e,% .7e\n",
          rp->inj,
          rp->forecast[n][0],
          rp->forecast[n][1],
          rp->forecast[n][2],
          v[n][0],
          v[n][1],
          v[n][2]);
      fclose(rp_fp);
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
			         forecasting
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/16/2023
void rp_deallocator(void *data, size_t sz, void *x){

	//this fn required for TensorFlow API, we don't use it though
	//it would be used by tensorflow to help with memory management
	//for other packagers like numpy whenever TF releases a tensor
	
	do_nothing(); 
}


//scrubbed 3/16/2023
TF_Tensor *rp_init_grayson_ivals(void *data, int *model_index){

    TF_Tensor *ivals;

    switch(*model_index){

	case 0:
	case 1:
	    rp->grayson_inputsz  = sizeof(double);
	    rp->grayson_inputsz *= rp->input_shape[STEPS]*rp->input_shape[FEATS];
	    rp->grayson_ndims    = 3,
	    rp->grayson_dims[0]  = 1;
	    rp->grayson_dims[1]  = rp->input_shape[STEPS];
	    rp->grayson_dims[2]  = rp->input_shape[FEATS];

	    ivals = TF_NewTensor(
		rp_fset_table[rp->fset].tensor_type,
		rp->grayson_dims,
		rp->grayson_ndims,
		data,
		rp->grayson_inputsz,
		&rp_deallocator,
		0
	    );
	break;

	case 2:
	    rp->grayson_inputsz  = sizeof(double);
	    rp->grayson_inputsz *= rp->input_shape[STEPS]*4;
	    rp->grayson_ndims    = 3,
	    rp->grayson_dims[0]  = 1;
	    rp->grayson_dims[1]  = rp->input_shape[STEPS];
	    rp->grayson_dims[2]  = 4;

	    ivals = TF_NewTensor(
		rp_fset_table[rp->fset].tensor_type,
		rp->grayson_dims,
		rp->grayson_ndims,
		data,
		rp->grayson_inputsz,
		&rp_deallocator,
		0
	    );
	   break;

	 case 3:
	    rp->grayson_inputsz  = sizeof(double);
	    rp->grayson_inputsz *= rp->input_shape[STEPS]*rp->input_shape[FEATS];
	    rp->grayson_ndims    = 3,
	    rp->grayson_dims[0]  = 1;
	    rp->grayson_dims[1]  = rp->input_shape[STEPS];
	    rp->grayson_dims[2]  = rp->input_shape[FEATS];

	    ivals = TF_NewTensor(
		rp_fset_table[rp->fset].tensor_type,
		rp->grayson_dims,
		rp->grayson_ndims,
		data,
		rp->grayson_inputsz,
		&rp_deallocator,
		0
	    );
	    model_index--;
	break;

	default:
	    ivals = TF_NewTensor(
		rp_fset_table[rp->fset].tensor_type,
		rp->grayson_dims,
		rp->grayson_ndims,
		data,
		rp->grayson_inputsz,
		&rp_deallocator,
		0
	    );
	  break;
    }

    return ivals;
}


//scrubbed 3/16/2023
double *rp_call_grayson(int n, void *data, int model_index){

    //setup the device-side input tensor
    TF_Tensor *ivals = rp_init_grayson_ivals(data,&model_index);

    //setup the device-side output tensor
    int64_t output_dims[2] = {1,rp_fset_table[rp->fset].output_width};
    TF_Tensor *ovals = TF_AllocateTensor(
	rp_fset_table[rp->fset].tensor_type, 			//TF_DataType
	output_dims, 		    	     			//dims
	2, 		     		     			//num_dims
	rp_fset_table[rp->fset].output_width*sizeof(double)	//len
    );

    //sanity check
    if (ivals == NULL || ovals == NULL){
       printf("RP_ERROR: Robin Protocol unable to allocate tensor(s) "
    	      "for Grayson call\n"); 
       exit(EXIT_FAILURE);
    }

    //call grayson, i.e. pass in ivals and populating ovals
    TF_SessionRun(
	rp->sesh[model_index],
	NULL,&rp->io_opts[model_index][0],&ivals,
	1,&rp->io_opts[model_index][1],&ovals,
	1,NULL,0,NULL,rp->status[n]
    );
   
    //sanity check
    if(TF_GetCode(rp->status[n]) != TF_OK){
        printf("%s",TF_Message(rp->status[n]));
        exit(EXIT_FAILURE);
    }

    //copy the predictions back to the host
    double *h_ovals = (double*)malloc(
	sizeof*h_ovals*rp_fset_table[rp->fset].output_width);

    memcpy(h_ovals,TF_TensorData(ovals),
	sizeof*h_ovals*rp_fset_table[rp->fset].output_width);

    //housekeeping
    TF_DeleteTensor(ivals);
    TF_DeleteTensor(ovals);

    free(data);

    return h_ovals;
}


#define OX n
#define H1 OX+1
#define H2 OX+2


//map historical record to grayson input, call grayson, compute forecasts
void rp_v1_vfg(rvec * gmx_restrict v,const rvec * gmx_restrict x,int n){

    //sanity check
    if (n%rp->modulus != 0){
    	rp->ml_tog[n]--;
        printf("-> RP_WARNING: vfg called for non-oxygen atom."
         " (this should not happen)\n");
  	return;
    }

    //get predictions for exterval, angular, and internal velocities
    double *e_pred = rp_call_grayson(n,rp_eiew_v1_map(n),0);
    double *w_pred = rp_call_grayson(n,rp_wiew_v0_map(n),1);
    double *i_pred = rp_call_grayson(n,rp_iiew_v0_map(n),2);

    //map the predictions back to the external frame forecast
    double ox[3], h1[3], h2[3], ev[3], com[3], eq_com[3];

    //apply i_pred to update the geometry
    {
      double m1, m2, mo, *norm;

      com[XX] = ((x[H1][XX]*1.0078)+(x[H2][XX]*1.0078)+(x[OX][XX]*15.999))/18.0146;
      com[YY] = ((x[H1][YY]*1.0078)+(x[H2][YY]*1.0078)+(x[OX][YY]*15.999))/18.0146;
      com[ZZ] = ((x[H1][ZZ]*1.0078)+(x[H2][ZZ]*1.0078)+(x[OX][ZZ]*15.999))/18.0146;

      //find the COM->OX vector
      ox[XX] = x[OX][XX]-com[XX];
      ox[YY] = x[OX][YY]-com[YY];
      ox[ZZ] = x[OX][ZZ]-com[ZZ];

      mo = mag(ox,3);

      //find the COM->H1 vector
      h1[XX] = x[H1][XX]-com[XX];
      h1[YY] = x[H1][YY]-com[YY];
      h1[ZZ] = x[H1][ZZ]-com[ZZ];

      m1 = mag(h1,3);

      //find the COM->H1 vector
      h2[XX] = x[H2][XX]-com[XX];
      h2[YY] = x[H2][YY]-com[YY];
      h2[ZZ] = x[H2][ZZ]-com[ZZ];

      m2 = mag(h2,3);

      //scale the "bonds" according to the predictions
      h1[XX] += (h1[XX]/m1)*(i_pred[0]);
      h1[YY] += (h1[YY]/m1)*(i_pred[0]);
      h1[ZZ] += (h1[ZZ]/m1)*(i_pred[0]);
 
      h2[XX] += (h2[XX]/m2)*(i_pred[2]);
      h2[YY] += (h2[YY]/m2)*(i_pred[2]);
      h2[ZZ] += (h2[ZZ]/m2)*(i_pred[2]);

      ox[XX] += (ox[XX]/mo)*(i_pred[3]);
      ox[YY] += (ox[YY]/mo)*(i_pred[3]);
      ox[ZZ] += (ox[ZZ]/mo)*(i_pred[3]);

      //3. update the H-COM/O-H angle
      norm = cross(h1,h2);

      rodrot(h1, norm,-i_pred[1]/2.0);
      rodrot(h2, norm, i_pred[1]/2.0);

      //housekeeping
      free(norm);
    }
        
   //apply w_pred to update the rotation
   {
     double *ibz, *iby, ibx[3];
     double tmp1[3],tmp2[3],m1,m2;
         
     //get the equilibrium com
     tmp1[XX] = x[H1][XX]-x[OX][XX];
     tmp1[YY] = x[H1][YY]-x[OX][YY];
     tmp1[ZZ] = x[H1][ZZ]-x[OX][ZZ];

     tmp2[XX] = x[H2][XX]-x[OX][XX];
     tmp2[YY] = x[H2][YY]-x[OX][YY];
     tmp2[ZZ] = x[H2][ZZ]-x[OX][ZZ];

     m1 = mag(tmp1,3);
     m2 = mag(tmp2,3);

     tmp1[XX] /= m1;
     tmp1[YY] /= m1;
     tmp1[ZZ] /= m1;

     tmp2[XX] /= m2;
     tmp2[YY] /= m2;
     tmp2[ZZ] /= m2;

     tmp1[XX] += x[OX][XX];
     tmp1[YY] += x[OX][YY];
     tmp1[ZZ] += x[OX][ZZ];

     tmp2[XX] += x[OX][XX];
     tmp2[YY] += x[OX][YY];
     tmp2[ZZ] += x[OX][ZZ];

     eq_com[XX] = ((x[OX][XX]*15.999)+(tmp1[XX]*1.0078)+(tmp2[XX]*1.0078))/18.0146;
     eq_com[YY] = ((x[OX][YY]*15.999)+(tmp1[YY]*1.0078)+(tmp2[YY]*1.0078))/18.0146;
     eq_com[ZZ] = ((x[OX][ZZ]*15.999)+(tmp1[ZZ]*1.0078)+(tmp2[ZZ]*1.0078))/18.0146;

     //get the bank axis
     ibx[XX] = x[OX][XX]-eq_com[XX];
     ibx[YY] = x[OX][YY]-eq_com[YY];
     ibx[ZZ] = x[OX][ZZ]-eq_com[ZZ];

     ibz = cross(ibx,h2);
     iby = cross(ibx,ibz);

     double mx,my,mz;
     mx = mag(ibx,3);
     my = mag(iby,3);
     mz = mag(ibz,3);

     for (int i=0;i<3;i++){
       ibx[i] /= mx;
       iby[i] /= my;
       ibz[i] /= -mz;
     }

     //apply the bank
     rodrot(h2 ,ibx,w_pred[2]);
     rodrot(h1 ,ibx,w_pred[2]);
     rodrot(ibz,ibx,w_pred[2]);
     rodrot(iby,ibx,w_pred[2]);

     //2. get the yaw axis and apply rotation about it
     rodrot(h1 ,ibz,w_pred[1]);
     rodrot(h2 ,ibz,w_pred[1]);
     rodrot(ox ,ibz,w_pred[1]);
     rodrot(ibx,ibz,w_pred[1]);
     rodrot(iby,ibz,w_pred[1]);

     //3. get the pitch axis and apply rotation about it
     rodrot(h1 ,iby,w_pred[0]);
     rodrot(h2 ,iby,w_pred[0]);
     rodrot(ox ,iby,w_pred[0]);
     rodrot(ibx,iby,w_pred[0]);
     rodrot(ibz,iby,w_pred[0]);
   }
	
   //apply e_pred to update the position of the com
   {
     ev[XX] = e_pred[XX]*rp->t0_magnitude[OX][0];
     ev[YY] = e_pred[YY]*rp->t0_magnitude[OX][0];
     ev[ZZ] = e_pred[ZZ]*rp->t0_magnitude[OX][0];
     
     //undo the rotations on e_pred - stored by a call to eiew_map
       yaw(ev,-rp->rotation[n][1]);
      bank(ev,-rp->rotation[n][0]);
     
     //apply the predicted translation to the com
     com[XX] += ev[XX];
     com[YY] += ev[YY];
     com[ZZ] += ev[ZZ];
   }

   //attach the com->X "bonds" back to the com to get external-frame positions
   ox[XX] += com[XX];
   ox[YY] += com[YY];
   ox[ZZ] += com[ZZ];

   h1[XX] += com[XX];
   h1[YY] += com[YY];
   h1[ZZ] += com[ZZ];

   h2[XX] += com[XX];
   h2[YY] += com[YY];
   h2[ZZ] += com[ZZ];

   //store the predicted velocities for subsequent application
   for (int d=0;d<3;d++){
       rp->forecast[OX][d] = (ox[d]-x[OX][d])*rp->dt_per_ps;
       rp->forecast[H1][d] = (h1[d]-x[H1][d])*rp->dt_per_ps;
       rp->forecast[H2][d] = (h2[d]-x[H2][d])*rp->dt_per_ps;
   }

   //housekeeping
   free(e_pred);
   free(w_pred);
   free(i_pred);
}


//map historical record to grayson input, call grayson, compute forecasts
void rp_v2_vfg(rvec * gmx_restrict v, const rvec * gmx_restrict x, int n){

    //sanity check
    if (n%rp->modulus != 0){
    	rp->ml_tog[n]--;
        printf("-> RP_WARNING: vfg called for non-oxygen atom."
         " (this should not happen)\n");
  	return;
    }

    //get predictions for exterval, angular, and internal velocities
    double *e_pred = rp_call_grayson(n,rp_eiew_v1_map(n),0);
    double *w_pred = rp_call_grayson(n,rp_wiew_v2_map(n),1);
    double *i_pred = rp_call_grayson(n,rp_iiew_v2_map(n),2);

    //map the predictions back to the external frame forecast
    double ox[3], h1[3], h2[3], ev[3], eq_com[3];

    //compute the equilibrium com
    rp_get_eqcom_from_gmx(n,x,eq_com);

    //apply i_pred to update the geometry
    {
      double m1, m2, mo, *norm;
     
      //find the COM->OX vector
      ox[XX] = x[OX][XX]-eq_com[XX];
      ox[YY] = x[OX][YY]-eq_com[YY];
      ox[ZZ] = x[OX][ZZ]-eq_com[ZZ];

      mo = mag(ox,3);

      //find the COM->H1 vector
      h1[XX] = x[H1][XX]-eq_com[XX];
      h1[YY] = x[H1][YY]-eq_com[YY];
      h1[ZZ] = x[H1][ZZ]-eq_com[ZZ];

      m1 = mag(h1,3);

      //find the COM->H1 vector
      h2[XX] = x[H2][XX]-eq_com[XX];
      h2[YY] = x[H2][YY]-eq_com[YY];
      h2[ZZ] = x[H2][ZZ]-eq_com[ZZ];

      m2 = mag(h2,3);

      //scale the "bonds" according to the predictions
      h1[XX] += (h1[XX]/m1)*(i_pred[0]);
      h1[YY] += (h1[YY]/m1)*(i_pred[0]);
      h1[ZZ] += (h1[ZZ]/m1)*(i_pred[0]);
 
      h2[XX] += (h2[XX]/m2)*(i_pred[2]);
      h2[YY] += (h2[YY]/m2)*(i_pred[2]);
      h2[ZZ] += (h2[ZZ]/m2)*(i_pred[2]);

      ox[XX] += (ox[XX]/mo)*(i_pred[3]);
      ox[YY] += (ox[YY]/mo)*(i_pred[3]);
      ox[ZZ] += (ox[ZZ]/mo)*(i_pred[3]);

      //3. update the H-COM/O-H angle
      norm = cross(h1,h2);

      rodrot(h1, norm,-i_pred[1]/2.0);
      rodrot(h2, norm, i_pred[1]/2.0);

      //housekeeping
      free(norm);
    }
        
   //apply w_pred to update the rotation
   {
     double *ibz, *iby, ibx[3];
         
     //get the bank axis
     ibx[XX] = ox[XX];
     ibx[YY] = ox[YY];
     ibx[ZZ] = ox[ZZ];

     ibz = cross(ibx,h2);
     iby = cross(ibx,ibz);

     double mx,my,mz;
     mx = mag(ibx,3);
     my = mag(iby,3);
     mz = mag(ibz,3);

     for (int i=0;i<3;i++){
       ibx[i] /= mx;
       iby[i] /= my;
       ibz[i] /= -mz;
     }

     //apply the bank
     rodrot(h2 ,ibx,w_pred[2]);
     rodrot(h1 ,ibx,w_pred[2]);
     rodrot(ibz,ibx,w_pred[2]);
     rodrot(iby,ibx,w_pred[2]);

     //2. get the yaw axis and apply rotation about it
     rodrot(h1 ,ibz,w_pred[1]);
     rodrot(h2 ,ibz,w_pred[1]);
     rodrot(ox ,ibz,w_pred[1]);
     rodrot(ibx,ibz,w_pred[1]);
     rodrot(iby,ibz,w_pred[1]);

     //3. get the pitch axis and apply rotation about it
     rodrot(h1 ,iby,w_pred[0]);
     rodrot(h2 ,iby,w_pred[0]);
     rodrot(ox ,iby,w_pred[0]);
     rodrot(ibx,iby,w_pred[0]);
     rodrot(ibz,iby,w_pred[0]);
   }
	
   //apply e_pred to update the position of the com
   {
      yaw(e_pred,-rp->rotation[n][1]);
     bank(e_pred,-rp->rotation[n][0]);
     
     //apply the predicted translation to the com
     eq_com[XX] += e_pred[XX]*rp->t0_magnitude[OX][0];
     eq_com[YY] += e_pred[YY]*rp->t0_magnitude[OX][0];
     eq_com[ZZ] += e_pred[ZZ]*rp->t0_magnitude[OX][0];
   }

   //attach the com->X "bonds" back to the com to get external-frame positions
   ox[XX] += eq_com[XX];
   ox[YY] += eq_com[YY];
   ox[ZZ] += eq_com[ZZ];

   h1[XX] += eq_com[XX];
   h1[YY] += eq_com[YY];
   h1[ZZ] += eq_com[ZZ];

   h2[XX] += eq_com[XX];
   h2[YY] += eq_com[YY];
   h2[ZZ] += eq_com[ZZ];

   //store the predicted velocities for subsequent application
   for (int d=0;d<3;d++){
       rp->forecast[OX][d] = (ox[d]-x[OX][d])*rp->dt_per_ps;
       rp->forecast[H1][d] = (h1[d]-x[H1][d])*rp->dt_per_ps;
       rp->forecast[H2][d] = (h2[d]-x[H2][d])*rp->dt_per_ps;
   }

   //housekeeping
   free(e_pred);
   free(w_pred);
   free(i_pred);
}


//map historical record to grayson input, call grayson, compute forecasts
void rp_v3_vfg(rvec * gmx_restrict v,const rvec * gmx_restrict x,int n){

    //sanity check
    if (n%rp->modulus != 0){
    	rp->ml_tog[n]--;
        printf("-> RP_WARNING: vfg called for non-oxygen atom."
         " (this should not happen)\n");
  	return;
    }

    //decrement the mlmd stride
    rp->ml_tog[n]--;

    //be verbose p.r.n.
    if (n==rp->reporter){
      printf("\n-> ML on for RP injection #%d\n",++rp->inj);
      fflush(stdout);
    }
  
    //get predictions for exterval, angular, and internal velocities
    double *e_pred = rp_call_grayson(n,rp_eiew_v1_map(n),0);
    double *w_pred = rp_call_grayson(n,rp_wiew_v3_map(n),1);
    double *i_pred = rp_call_grayson(n,rp_iiew_v3_map(n),3);

    //map the predictions back to the external frame forecast
    double ox[3], h1[3], h2[3], ev[3], eq_com[3], wv[3];

    //compute the equilibrium com
    rp_get_eqcom_from_gmx(n,x,eq_com);

    //apply i_pred to update the geometry
    {
      double m1, m2, mo, *norm;
      
      //find the COM->OX vector
      ox[XX] = x[OX][XX]-eq_com[XX];
      ox[YY] = x[OX][YY]-eq_com[YY];
      ox[ZZ] = x[OX][ZZ]-eq_com[ZZ];

      mo = mag(ox,3);

      //find the COM->H1 vector
      h1[XX] = x[H1][XX]-eq_com[XX];
      h1[YY] = x[H1][YY]-eq_com[YY];
      h1[ZZ] = x[H1][ZZ]-eq_com[ZZ];

      m1 = mag(h1,3);

      //find the COM->H1 vector
      h2[XX] = x[H2][XX]-eq_com[XX];
      h2[YY] = x[H2][YY]-eq_com[YY];
      h2[ZZ] = x[H2][ZZ]-eq_com[ZZ];

      m2 = mag(h2,3);

      //scale the "bonds" according to the predictions
      h1[XX] += (h1[XX]/m1)*(i_pred[0]);
      h1[YY] += (h1[YY]/m1)*(i_pred[0]);
      h1[ZZ] += (h1[ZZ]/m1)*(i_pred[0]);
 
      h2[XX] += (h2[XX]/m2)*(i_pred[2]);
      h2[YY] += (h2[YY]/m2)*(i_pred[2]);
      h2[ZZ] += (h2[ZZ]/m2)*(i_pred[2]);

      //3. update the H-COM/O-H angle
      norm = cross(h1,h2);

      rodrot(h1, norm,-(i_pred[1]/2.0));
      rodrot(h2, norm, (i_pred[1]/2.0));

      //housekeeping
      free(norm);
    }
        
   //apply w_pred to update the rotation
   {
     double *ibz, *iby, ibx[3];
         
     //get the bank axis
     ibx[XX] = x[OX][XX]-eq_com[XX];
     ibx[YY] = x[OX][YY]-eq_com[YY];
     ibx[ZZ] = x[OX][ZZ]-eq_com[ZZ];

     ibz = cross(ibx,h2);
     iby = cross(ibx,ibz);

     double mx,my,mz;
     mx = mag(ibx,3);
     my = mag(iby,3);
     mz = mag(ibz,3);

     for (int i=0;i<3;i++){
       ibx[i] /= mx;
       iby[i] /= my;
       ibz[i] /= -mz;
     }
 
     wv[XX] = w_pred[XX]*rp->t0_magnitude[OX][1];
     wv[YY] = w_pred[YY]*rp->t0_magnitude[OX][1];
     wv[ZZ] = w_pred[ZZ]*rp->t0_magnitude[OX][1];

     //apply the bank
     rodrot(h2 ,ibx,wv[2]);
     rodrot(h1 ,ibx,wv[2]);
     rodrot(ibz,ibx,wv[2]);
     rodrot(iby,ibx,wv[2]);

     //apply the yaw
     rodrot(h1 ,ibz,wv[1]);
     rodrot(h2 ,ibz,wv[1]);
     rodrot(ox ,ibz,wv[1]);
     rodrot(ibx,ibz,wv[1]);
     rodrot(iby,ibz,wv[1]);

     //apply the pitch
     rodrot(h1 ,iby,wv[0]);
     rodrot(h2 ,iby,wv[0]);
     rodrot(ox ,iby,wv[0]);
     rodrot(ibx,iby,wv[0]);
     rodrot(ibz,iby,wv[0]);
   }
	
   //apply e_pred to update the position of the com
   {
     ev[XX] = e_pred[XX]*rp->t0_magnitude[OX][0];
     ev[YY] = e_pred[YY]*rp->t0_magnitude[OX][0];
     ev[ZZ] = e_pred[ZZ]*rp->t0_magnitude[OX][0];
     
     //undo the rotations on e_pred - stored by a call to eiew_map
       yaw(ev,-rp->rotation[n][1]);
      bank(ev,-rp->rotation[n][0]);
     
     //apply the predicted translation to the com
     eq_com[XX] += ev[XX];
     eq_com[YY] += ev[YY];
     eq_com[ZZ] += ev[ZZ];
   }

   //attach the com->X "bonds" back to the com to get external-frame positions
   ox[XX] += eq_com[XX];
   ox[YY] += eq_com[YY];
   ox[ZZ] += eq_com[ZZ];

   h1[XX] += eq_com[XX];
   h1[YY] += eq_com[YY];
   h1[ZZ] += eq_com[ZZ];

   h2[XX] += eq_com[XX];
   h2[YY] += eq_com[YY];
   h2[ZZ] += eq_com[ZZ];

   //store the predicted velocities for subsequent application
   for (int d=0;d<3;d++){
       rp->forecast[OX][d] = (ox[d]-x[OX][d])*rp->dt_per_ps;
       rp->forecast[H1][d] = (h1[d]-x[H1][d])*rp->dt_per_ps;
       rp->forecast[H2][d] = (h2[d]-x[H2][d])*rp->dt_per_ps;
   }

   //housekeeping
   free(e_pred);
   free(w_pred);
   free(i_pred);
}


//map historical record to grayson input, call grayson, compute forecasts
void rp_v0_vfg(rvec * gmx_restrict v,const rvec * gmx_restrict x,int n){

    //sanity check
    if (n%rp->modulus != 0){
    	rp->ml_tog[n]--;
        printf("-> RP_WARNING: vfg called for non-oxygen atom."
         " (this should not happen)\n");
  	return;
    }

    //decrement the mlmd stride
    rp->ml_tog[n]--;

    //be verbose p.r.n.
    if (n==rp->reporter){
      printf("\n-> ML on for RP injection #%d\n",++rp->inj);
      fflush(stdout);
    }
  
    double tmp1[3], tmp2[3];

    //get predictions for exterval, angular, and internal velocities
    double *e_pred = rp_call_grayson(n,rp_eiew_v0_map(n),0);
    double *w_pred = rp_call_grayson(n,rp_wiew_v0_map(n),1);
    double *i_pred = rp_call_grayson(n,rp_iiew_v0_map(n),2);

    //map the predictions back to the external frame forecast
    double ox[3], h1[3], h2[3], ev[3], com[3];

    //find the current COM
    {
      com[XX] = ((x[H1][XX]*1.0078)+(x[H2][XX]*1.0078)+(x[OX][XX]*15.999))/18.0146;
      com[YY] = ((x[H1][YY]*1.0078)+(x[H2][YY]*1.0078)+(x[OX][YY]*15.999))/18.0146;
      com[ZZ] = ((x[H1][ZZ]*1.0078)+(x[H2][ZZ]*1.0078)+(x[OX][ZZ]*15.999))/18.0146;
    }

    //apply i_pred to update the geometry
    {
      double m1, m2, mo, *norm;

      //find the COM->OX vector
      ox[XX] = x[OX][XX]-com[XX];
      ox[YY] = x[OX][YY]-com[YY];
      ox[ZZ] = x[OX][ZZ]-com[ZZ];

      mo = mag(ox,3);

      //find the COM->H1 vector
      h1[XX] = x[H1][XX]-com[XX];
      h1[YY] = x[H1][YY]-com[YY];
      h1[ZZ] = x[H1][ZZ]-com[ZZ];

      m1 = mag(h1,3);

      //find the COM->H1 vector
      h2[XX] = x[H2][XX]-com[XX];
      h2[YY] = x[H2][YY]-com[YY];
      h2[ZZ] = x[H2][ZZ]-com[ZZ];

      m2 = mag(h2,3);

      //scale the "bonds" according to the predictions
      h1[XX] += (h1[XX]/m1)*(i_pred[0]);
      h1[YY] += (h1[YY]/m1)*(i_pred[0]);
      h1[ZZ] += (h1[ZZ]/m1)*(i_pred[0]);
 
      h2[XX] += (h2[XX]/m2)*(i_pred[2]);
      h2[YY] += (h2[YY]/m2)*(i_pred[2]);
      h2[ZZ] += (h2[ZZ]/m2)*(i_pred[2]);

      ox[XX] += (ox[XX]/mo)*(i_pred[3]);
      ox[YY] += (ox[YY]/mo)*(i_pred[3]);
      ox[ZZ] += (ox[ZZ]/mo)*(i_pred[3]);

      //3. update the H-COM/O-H angle
      norm = cross(h1,h2);

      rodrot(h1, norm,-i_pred[1]/2.0);
      rodrot(h2, norm, i_pred[1]/2.0);

      //housekeeping
      free(norm);
    }
	
   //apply e_pred to update the position of the com
   {
     //this prediction is in nm/dt, but is rotated s.t. ox->com 
     //  is on the x-axis and h1->com is in xy plane
     ev[XX] = e_pred[XX];
     ev[YY] = e_pred[YY];
     ev[ZZ] = e_pred[ZZ];
     
     //undo the rotations on e_pred - stored by a call to eiew_map
      bank(ev,-rp->rotation[n][2]);
       yaw(ev,-rp->rotation[n][1]);
      bank(ev,-rp->rotation[n][0]);
     
     //apply the predicted translation to the com
     com[XX] += ev[XX];
     com[YY] += ev[YY];
     com[ZZ] += ev[ZZ];
   }
        
   //apply w_pred to update the rotation
   {
     double *ibz, *iby, ibx[3];
     double tmp1[3],tmp2[3],m1,m2,eq_com[3];
         
     //get the equilibrium com
     tmp1[XX] = x[H1][XX]-x[OX][XX];
     tmp1[YY] = x[H1][YY]-x[OX][YY];
     tmp1[ZZ] = x[H1][ZZ]-x[OX][ZZ];

     tmp2[XX] = x[H2][XX]-x[OX][XX];
     tmp2[YY] = x[H2][YY]-x[OX][YY];
     tmp2[ZZ] = x[H2][ZZ]-x[OX][ZZ];

     m1 = mag(tmp1,3);
     m2 = mag(tmp2,3);

     tmp1[XX] /= m1;
     tmp1[YY] /= m1;
     tmp1[ZZ] /= m1;

     tmp2[XX] /= m2;
     tmp2[YY] /= m2;
     tmp2[ZZ] /= m2;

     tmp1[XX] += x[OX][XX];
     tmp1[YY] += x[OX][YY];
     tmp1[ZZ] += x[OX][ZZ];

     tmp2[XX] += x[OX][XX];
     tmp2[YY] += x[OX][YY];
     tmp2[ZZ] += x[OX][ZZ];

     eq_com[XX] = ((x[OX][XX]*15.999)+(tmp1[XX]*1.0078)+(tmp2[XX]*1.0078))/18.0146;
     eq_com[YY] = ((x[OX][YY]*15.999)+(tmp1[YY]*1.0078)+(tmp2[YY]*1.0078))/18.0146;
     eq_com[ZZ] = ((x[OX][ZZ]*15.999)+(tmp1[ZZ]*1.0078)+(tmp2[ZZ]*1.0078))/18.0146;

     //get the bank axis
     ibx[XX] = x[OX][XX]-eq_com[XX];
     ibx[YY] = x[OX][YY]-eq_com[YY];
     ibx[ZZ] = x[OX][ZZ]-eq_com[ZZ];

     ibz = cross(ibx,h2);
     iby = cross(ibx,ibz);

     double mx,my,mz;
     mx = mag(ibx,3);
     my = mag(iby,3);
     mz = mag(ibz,3);

     for (int i=0;i<3;i++){
       ibx[i] /= mx;
       iby[i] /= my;
       ibz[i] /= -mz;
     }

     //apply the bank
    rodrot(h2 ,ibx,w_pred[2]);
    rodrot(h1 ,ibx,w_pred[2]);
    rodrot(ibz,ibx,w_pred[2]);
    rodrot(iby,ibx,w_pred[2]);

    //2. get the yaw axis and apply rotation about it
    rodrot(h1 ,ibz,w_pred[1]);
    rodrot(h2 ,ibz,w_pred[1]);
    rodrot(ox ,ibz,w_pred[1]);
    rodrot(ibx,ibz,w_pred[1]);
    rodrot(iby,ibz,w_pred[1]);

    //3. get the pitch axis and apply rotation about it
    rodrot(h1 ,iby,w_pred[0]);
    rodrot(h2 ,iby,w_pred[0]);
    rodrot(ox ,iby,w_pred[0]);
    rodrot(ibx,iby,w_pred[0]);
    rodrot(ibz,iby,w_pred[0]);
   }

   //attach the com->X "bonds" back to the com to get external-frame positions
   ox[XX] += com[XX];
   ox[YY] += com[YY];
   ox[ZZ] += com[ZZ];

   h1[XX] += com[XX];
   h1[YY] += com[YY];
   h1[ZZ] += com[ZZ];

   h2[XX] += com[XX];
   h2[YY] += com[YY];
   h2[ZZ] += com[ZZ];

   //store the predicted velocities for subsequent application
   for (int d=0;d<3;d++){
       rp->forecast[OX][d] = (ox[d]-x[OX][d])*rp->dt_per_ps;
       rp->forecast[H1][d] = (h1[d]-x[H1][d])*rp->dt_per_ps;
       rp->forecast[H2][d] = (h2[d]-x[H2][d])*rp->dt_per_ps;
   }

   //housekeeping
   free(e_pred);
   free(w_pred);
   free(i_pred);
}


//scrubbed 3/16/2023
void rp_og_vfg(rvec * gmx_restrict v,const rvec * gmx_restrict x,int n){

    //get the atom's group index, which is also its model index
    int g = rp_get_group(n);

    //flatten atom history and map to Grayson input
    void *data = (*rp->map)(n); 

    //call grayson and apply the reverse map to the prediction
    double *pred = rp_call_grayson(n,data,g); 
    (*rp->pam)(n,x,pred); 

    //housekeeping
    free(pred);
}


//scrubbed 3/15/2023
bool rp_injection_circuit( const rvec * gmx_restrict f,
  rvec * gmx_restrict v,const rvec * gmx_restrict x, int a, int nrend){

    if (rp == NULL)
      return false;

    //make sure the robin protocol is active
    else if (!rp->active)
      return false;

    //handle whole molecule models
    else if (rp->whole_molecule_model == true)
      for (int n=a;n<a+GMX_SIMD_REAL_WIDTH;n++){
        if (n >= nrend)
          break;
        
	//make sure this is an ML step and only work on n%mod 
        if (rp->ml_tog[n] > 0 && n%rp->modulus==0 ){
      
   	 //decrement the mlmd stride
    	 rp->ml_tog[n]--;

	 //be verbose p.r.n.
    	 if (n==rp->reporter){
      	   printf("\n-> ML on for (whole-mol) RP injection #%d\n",++rp->inj);
           fflush(stdout);
    	 }

	 //get velocities from grayson
         (*rp->vfg)(v,x,n); 
 
   	 //reporter vebosity
   	 if (n == rp->reporter)
           rp_describe_injection(n,v,false);

	 //active dump 
         if (rp->dump_vels)
	   rp_active_dump(n,v);

	 //inject the forecasted velocities
         v[n+0][XX] = rp->forecast[n+0][XX];
         v[n+0][YY] = rp->forecast[n+0][YY];
         v[n+0][ZZ] = rp->forecast[n+0][ZZ];
      
         v[n+1][XX] = rp->forecast[n+1][XX];
         v[n+1][YY] = rp->forecast[n+1][YY];
         v[n+1][ZZ] = rp->forecast[n+1][ZZ];
      
         v[n+2][XX] = rp->forecast[n+2][XX];
         v[n+2][YY] = rp->forecast[n+2][YY];
         v[n+2][ZZ] = rp->forecast[n+2][ZZ];
        }
      }

    //handle single atom models
    else
      for (int n=a;n<a+GMX_SIMD_REAL_WIDTH;n++){
        if (n >= nrend)
          break;
 
        //make sure this is an ML step 
        if (rp->ml_tog[n] > 0){
     
   	 //decrement the mlmd stride
    	 rp->ml_tog[n]--;

	 //be verbose p.r.n.
    	 if (n==rp->reporter){
      	   printf("\n-> ML on for (single-atom) RP injection #%d\n",++rp->inj);
           fflush(stdout);
    	 }
      
	//get a velocity from grayson
        (*rp->vfg)(v,x,n); 
  
   	//reporter vebosity
   	if (n == rp->reporter)
          rp_describe_injection(n,v,true);
 
	//active dump 
        if (rp->dump_vels)
	  rp_active_dump(n,v);
     
	//inject the forecasted velocity
        v[n][XX] = rp->forecast[n][XX];
        v[n][YY] = rp->forecast[n][YY];
        v[n][ZZ] = rp->forecast[n][ZZ];

	//check if we're blowing up and try to exit gracefully
	if ( fabs(v[n][XX]) > rp->box[XX] 
	  || fabs(v[n][YY]) > rp->box[YY] 
	  || fabs(v[n][ZZ]) > rp->box[ZZ]){
	    printf("\n-> RP injected velocity exceeds box dimensions.\n"
		   "-> ROBIN PROTOCOL terminated. Goodbye.\n");	
	    rp->active = false;
	}
      }
    }

    return true;
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
			        miscellaneous				  
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/15/2023
void rp_splash(void){
    puts("\n\n\n"
    "                                &@@@@@@@@@@@@@@@@@@@(                  \n"
    "                        @@@@@@,..........................@@            \n"
    "                  ,@@@@......................................@         \n"
    "              &@@*................*%@@@@@@(....................@       \n"
    "           @@....../@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@      \n"
    "        @@.&@@@@@/ @@@@@@@@@@@@......@@@@@@@@@@@@@@@............@      \n"
    "     (@@        @@@@@@@@@@@..........@@@@@@@@@@@@@@@............@      \n"
    "              @@@@@@@@@@@@@..........@@@@@@@@@@@@@@%...........@@@     \n"
    "            .@@@@@@@@@@@@@@.........@@@@@@@@@@@@@@............@@@@@@   \n"
    "           (@@@@@@@@@@@@@@..........@@@@@@@@@@@,............@@@@@@@@@  \n"
    "           @@@@@@@@@@@@@@@..........@@@@@(...............@@@@@@@@@@@@  \n"
    "           @@@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@& \n"
    "           @@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@@  \n"
    "           @@@@@@@@@@@@@@............................#@@@@@@@@@@@@@@@  \n"
    "            @@@@@@@@@@@@(...........#@@@@..............@@@@@@@@@@@@@   \n"
    "              @@@@@@@@@@............@@@@@@@@@@@..........@@@@@@@@@     \n"
    "                @@@@@@@@............@@@@@@@@@@@@@@@#......@@@@@@       \n"
    "                   @@@@............,@@@@@@@@@@@@@@@@@@@/....@@         \n"
    "                      @.........@@@@@@@@@@@@@@@@@@@@@@@@@%@...@        \n"
    "                     @,....&@@@@@@@@@@@@@@@@@@@@@@@.         @*@%      \n"
    "                     *(             Protocol 5.0                .@   \n");
    fflush(stdout);
}


//scrubbed 3/16/2023
void rp_describe_injection(int n,rvec * gmx_restrict v,bool single_atom){

   for (int a=n;a<n+3;a++){
       fprintf(stdout,"atom %d -> inj%d: (% lf, % lf,% lf)"
		" target: (% lf, % lf, % lf) error: ( % lf, %lf, %lf)\n",
          a,
          rp->inj,
          rp->forecast[a][XX],
          rp->forecast[a][YY],
          rp->forecast[a][ZZ],
          v[a][XX],
          v[a][YY],
          v[a][ZZ],
	  v[a][XX]-rp->forecast[a][XX],
          v[a][YY]-rp->forecast[a][YY],
          v[a][ZZ]-rp->forecast[a][ZZ]);
	if (single_atom)
		break;
	}

   fflush(stdout);
}


//scrubbed 3/15/2023
//handles group string entries, be them ranges or stand-alone values
int rp_parse_group_substr(char *entry, int i, robin_protocol *rp){

    //for single indices
    char *p;
    if ((p = strchr(entry,'-')) == NULL){
      //adjust memory allocation
      rp->groups[i] = (int*)realloc(rp->groups[i],
                        sizeof**rp->groups*rp->grpsz[i]+1);
      //store entry
      rp->groups[i][rp->grpsz[i]] = atoi(entry);
      return 1;
    }

    //for ranges-of-indices 
    char *dup = strdup(entry); 

    *p = 0;
    dup = strchr(dup,'-');
    dup++;

    int b = atoi(entry), e = atoi(dup)+1;

    //adjust memory allocation
    rp->groups[i] = (int*)realloc(rp->groups[i],
                      sizeof**rp->groups*rp->grpsz[i]+(e-b));

    //store entire range of values
    for (int j=b;j<e;j++)
      rp->groups[i][rp->grpsz[i]+(j-b)] = j;

    return e-b;
}


//scrubbed 3/15/2023
//parses and processes the group string
int rp_parse_group_str(robin_protocol *rp, t_inputrec *ir){

    int modulus = 0; char *tok, *rest = ir->rp_groups;

    //loop over each detected group
    for (int i=0;i<rp->Ngr;i++){

      //pull a single group & strip the enclosing brackets
      tok = strtok_r(rest,":",&rest); 
      tok++;
      *strchr(tok,']') = 0;

      //process the group, entry by entry (the loop is leaky Â¯\_(ã)_/Â¯)
      if (strchr(tok,',') == NULL)
        rp->grpsz[i] += rp_parse_group_substr(tok,i,rp);
      else while (tok != NULL){
        char *dup = strdup(tok), *p; 

        if ((p = strchr(dup,',')) != NULL)
          *p = 0; 

        rp->grpsz[i] += rp_parse_group_substr(dup,i,rp);

        if ((tok = strchr(tok,',')) != NULL)
          tok++;
      }
      modulus += rp->grpsz[i];
    }

    printf("RP: modulus computed to be %d\n", modulus);
    return modulus;
}


//scrubbed 3/15/2023
void rp_parse_protocol_pattern(robin_protocol *p, char *pattern){
    char *tok, *rest = pattern;

    p->pattern[VALI_DELAY] = 0;
    p->pattern[TEST_DELAY] = 0;

    //first term in the delay
    tok = strtok_r(rest,",",&rest);
    p->pattern[DELAY] = atoi(tok);

    //second term is the stride
    tok = strtok_r(rest,",",&rest);

    int nmol = p->natoms/p->modulus;
    printf("RP: identified %d molecules comprising %d atoms\n",
	nmol,p->natoms);

    /* we interpret stride as 1 frame of data
       so we scale this by the # of waters so it works
       for both isolated water and bulk */

    //how many timesteps we need to take
    p->pattern[STRIDE] = atof(tok); 
   
    //TODO: for diatomic molecules, if both atoms contribute then we need to half the stride hack
    if (!p->active)
      p->pattern[STRIDE] /= 2.0;

    //timesteps * #water = how many samples we actually dump
    p->pattern[TRUE_STRIDE] = p->pattern[STRIDE]*nmol; 

    //set validation and testing stride to be half of the training stride
    p->pattern[VALI_STRIDE] = (int)ceil(p->pattern[STRIDE]*.5);
    p->pattern[TEST_STRIDE] = (int)ceil(p->pattern[STRIDE]*.5);
  
    //third term is the pause
    tok = strtok_r(rest,",",&rest);
    p->pattern[PAUSE] = atoi(tok);

    //fourth term is the number of rounds
    tok = strtok_r(rest,",",&rest);
    p->pattern[ROUNDS] = atoi(tok);

    //adjust the delay p.r.n. to accomodate the input requirment
    if (p->pattern[DELAY] < p->histr_shape[STEPS]){
        printf("RP: activation delay changed\n"
	       "\t-> user-specified: %d, updated: %d.\n",
	 	p->pattern[DELAY],p->histr_shape[STEPS]+1); //added +1 here and below at assignment for reproduction
	if (rp_fset_table[p->fset].sample_bloat != 0 && p->pattern[DELAY] >= p->histr_shape[STEPS]-rp_fset_table[p->fset].sample_bloat)
		printf("this adjustment is due to:\n\t1. sample bloat inherint to the feature set\n");
        if ( !p->active && p->non_active_strech)
		printf("\t2. non-active strech (=%d)for dumping targets\n",p->non_active_strech);
	
        p->pattern[DELAY] = p->histr_shape[STEPS] + 1;
    }

    //compute the final step for the active protocol to be enabled
    p->final_rp_step  = p->pattern[STRIDE]+p->pattern[PAUSE];
    p->final_rp_step *= p->pattern[ROUNDS];
    p->final_rp_step += p->pattern[DELAY];

    printf("\nRP_PATTERN:\n"
	   "\tdelay=%d\n"
	   "\tstride=%d\n"
	   "\ttrue stride=%d\n"
	   "\tpause=%d\n"
	   "\trounds=%d\n"
	   "\tvali-stride=%d\n"
	   "\ttest-stride=%d\n"
	   "\tfinal_step=%d\n\n",
	   p->pattern[DELAY], 
	   p->pattern[STRIDE],
	   p->pattern[TRUE_STRIDE],
	   p->pattern[PAUSE],
	   p->pattern[ROUNDS],
	   p->pattern[VALI_STRIDE],
	   p->pattern[TEST_STRIDE],
	   p->final_rp_step);

    //check for alternate pattern
    if (p->alternate_patterns && p->alternative_stride != -1){
	for (int j=0;j<8;j++)
	  p->pattern2[j] = p->pattern[j];

	p->pattern2[STRIDE] = p->alternative_stride;
	p->pattern2[PAUSE] = p->alternative_pause;

        printf("\nRP_ALTERNATE_PATTERN:\n"
	       "\tdelay=%d\n"
	       "\tstride=%d\n"
	       "\ttrue stride=%d\n"
	       "\tpause=%d\n"
	       "\trounds=%d\n"
	       "\tvali-stride=%d\n"
	       "\ttest-stride=%d\n"
	       "\tfinal_step=%d\n\n",
		p->pattern2[DELAY], 
		p->pattern2[STRIDE],
		p->pattern2[TRUE_STRIDE],
		p->pattern2[PAUSE],
		p->pattern2[ROUNDS],
		p->pattern2[VALI_STRIDE],
		p->pattern2[TEST_STRIDE],
		p->final_rp_step);
    }

    fflush(stdout);
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
			      memory management				  
\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


//scrubbed 3/15/2023
void init_rp_tensorflow_api(robin_protocol *p){

    struct stat buffer;
    const char *tags = "serve";
    TF_SessionOptions *sesh_opts = TF_NewSessionOptions();
    
    //for each user specified group in the input record (.mdp)
    for (int i=0;i<p->Ngr;i++){

      //make sure a model was specified and is accessible
      if (p->models[i] == NULL || stat(p->models[i],&buffer) == -1)
        continue;
     
      //allocate a new tensorflow graph 
      p->graph[i] = TF_NewGraph();
    
      //attempt to load the saved model
      p->sesh[i]  = TF_LoadSessionFromSavedModel(sesh_opts,NULL,
                    p->models[i],&tags,1,p->graph[i],NULL,p->status[i]);
    
      //initialize the models input options
      p->io_opts[i][0] = {
        TF_GraphOperationByName(p->graph[i],"serving_default_input_1"),0
      };
    
      //initialize the models output options
      p->io_opts[i][1] = {
        TF_GraphOperationByName(p->graph[i],"StatefulPartitionedCall"),0
      };
    
      //check to make sure things went according to plan here
      if(p->io_opts[i][0].oper == NULL || p->io_opts[i][1].oper == NULL){
            printf("ERROR: Failed GraphOperationByName\n");
            abort();
      }
    
      if (TF_GetCode(p->status[i]) != TF_OK){
          printf("Robin protocol failure.\nrp_model%d:%s\n",i,
          TF_Message(p->status[i]));
          abort();
      }

	printf("RP: loaded model for group %d: %s\n", i, p->models[i]);
    }
   
    //set the dimensions used for model calls
    //somewhat redundant bc of rp_call_grayson
    p->grayson_inputsz  = sizeof(double);
    p->grayson_inputsz  *= p->input_shape[STEPS]*p->input_shape[FEATS];
    p->grayson_ndims    = 3, 
    p->grayson_dims[0]  = 1;
    p->grayson_dims[1]  = p->input_shape[STEPS];
    p->grayson_dims[2]  = p->input_shape[FEATS];

    //housekeeping
    TF_DeleteSessionOptions(sesh_opts);
}


//scrubbed 3/15/2023
void init_rp_groups(robin_protocol *p, t_inputrec *ir){

      for (int i=0;i<strlen(ir->rp_groups);i++)
        if (ir->rp_groups[i] == '[')
          p->Ngr += 1;

      p->groups = (int **)malloc(sizeof*(p->groups)*p->Ngr);

      p->grpsz = (int *)calloc(p->Ngr,sizeof*(p->grpsz ));

      for (int i=0;i<p->Ngr;i++)
        p->groups[i] = (int *)malloc(0);

      p->modulus = rp_parse_group_str(p,ir);

      //init the protocol pattern (i.e. delay, stride, pause, rounds)  
      rp_parse_protocol_pattern(p, ir->rp_pattern);
    }


//scrubbed 3/15/2023
void init_rp_fset(robin_protocol *p, t_inputrec *ir){

    bool fMatch = false;

    for (int i=0;rp_fset_table[i].name != nullptr;i++){
      if (!strcmp(rp_fset_table[i].name,rp_fsets[ir->rp_fset])){
	 fMatch = true;
	 printf("RP: using %s feature set\n",rp_fsets[ir->rp_fset]);

         p->fset 		  = i;
         p->input_shape[FEATS] 	  = p->active ? rp_fset_table[i].input_width : rp_fset_table[i].dump_width;
         p->histr_shape[FEATS] 	  = rp_fset_table[i].capture_pos ? 6:3;
         p->rec 		  = rp_fset_table[i].capture_pos ? 
					&rp_rec3v3r: &rp_rec3v;
         p->map 		  = rp_fset_table[i].map;
         p->pam 		  = rp_fset_table[i].pam;
         p->vfg 		  = rp_fset_table[i].vfg;
         p->whole_molecule_model  = rp_fset_table[i].whole_molecule_model;
         p->histr_shape[STEPS] 	 += rp_fset_table[i].sample_bloat;

         if (!p->active){
           p->histr_shape[STEPS] += p->non_active_strech;
           p->input_shape[STEPS] += p->non_active_strech;
         }

         break;
      }
    }

    if (!fMatch){
	printf("\nRP_WARNING: The specified feature set '%s' does not exist.\n",
		rp_fsets[ir->rp_fset]);
        exit(EXIT_FAILURE);
    }
}


//scrubbed 3/15/2023
void init_rp_records(robin_protocol *p, t_state *state){

    p->history 	    = ( real***)malloc(sizeof*p->history*state->natoms);
    p->rotation     = (double**)malloc(sizeof*p->rotation*state->natoms);
    p->t0_magnitude = (double**)malloc(sizeof*p->t0_magnitude*state->natoms);

    //initialize malloc'd arrays p.r.n.
    for (int i=0;i<state->natoms;i++){
      p->history[i]	 = ( real**)malloc(
	sizeof*p->history[i]*p->histr_shape[STEPS]);

      p->rotation[i] 	 = (double*)malloc(sizeof*p->rotation[i]*9);
      p->t0_magnitude[i] = (double*)malloc(sizeof*p->t0_magnitude[i]*8);

      for (int j=0;j<p->histr_shape[STEPS];j++)
        p->history[i][j] = (real*)malloc(
          sizeof*p->history[i][j]*p->histr_shape[FEATS]);
    }
}


//scrubbed 3/15/2023
void init_rp_forecasting(robin_protocol *p, t_state *state, t_inputrec *ir){

  if (p->active){
    p->ml_tog 	= (         int*)malloc(sizeof*p->ml_tog*state->natoms);
    p->forecast = (     double**)malloc(sizeof*p->forecast*state->natoms);
    p->status 	= (  TF_Status**)malloc(sizeof*p->status*state->natoms);
    p->models 	= (const char **)malloc(sizeof*p->models*p->Ngr);
    p->graph 	= (  TF_Graph **)malloc(sizeof*p->graph*p->Ngr);
    p->io_opts 	= ( TF_Output **)malloc(sizeof*p->io_opts*p->Ngr);
    p->sesh 	= (TF_Session **)malloc(sizeof*p->sesh*p->Ngr);
    p->models 	= (const char **)malloc(sizeof*p->models*p->Ngr);
    
    for (int i=0;i<state->natoms;i++){
      p->ml_tog[i]   = 0;
      p->status[i]   = TF_NewStatus();
      p->forecast[i] = (double*)malloc(sizeof*p->forecast[i]*3);
    }
    
    for (int i=0;i<p->Ngr;i++){
      p->sesh[i] = NULL;

      if (ir->rp_models[i] != NULL){
        p->models[i] 	= ir->rp_models[i];
        p->io_opts[i] 	= (TF_Output *)malloc(sizeof**(p->io_opts)*2);
      }
      else
        p->models[i] = NULL;
    }
    
    //initialize the tensorflow api
    init_rp_tensorflow_api(p);
  }
}


//scrubbed 3/15/2023
void init_rp_nad(robin_protocol *p){

    //balanced diet 
    if (p->balanced_diet == true){
      p->dumped = (int**)malloc(sizeof*p->dumped*p->Ngr);
      for (int i=0;i<p->Ngr;i++)
        p->dumped[i] = (int*)calloc(8,sizeof**p->dumped);
    }
    else if (p->threshold != 0){
      printf("dumping threshold set to %e\n",p->threshold);
      fflush(stdout);

      p->threshold_domain = 3;
      char *name=strdup(rp_fset_table[p->fset].name);

      if (!strcmp((name+(strlen(name)-3)),"ROT") && name[5] != 'v')
	p->threshold_domain = 2;

      free(name);

      printf("will apply threshold to first %d features..\n",p->threshold_domain);
      fflush(stdout);
    }

    p->dumped = (int**)malloc(sizeof*p->dumped*p->Ngr);
    for (int i=0;i<p->Ngr;i++)
      p->dumped[i] = (int*)calloc(1,sizeof**p->dumped);
}
 

//scrubbed 3/15/2023
void announce_new_rp(robin_protocol *p){

    //anounce domain constraint p.r.n.
    if (p->domain > -1)
      printf("RP: constrained domain is enabled\n"
             "\t-> operating on residue %d\n",
             p->domain);

    //anounce time dilation p.r.n.
    if (p->non_active_strech != 1)
      printf("RP: time dilation is enabled (factor = %d)\n"
             "\t-> historical record length increased from %d to %d\n",
              p->dilation_factor,p->histr_shape[STEPS]/p->dilation_factor,
              p->histr_shape[STEPS]);

    //be verbose
    printf("\nRobin protocol initiated. "
           "\"Holy machine-learned molecular dynamics Batman!"
           "\" (Dick Grayson)\n\n");
}


//deprecated 3/15/2023
void init_rp_time_dilation(robin_protocol *p){
    if (p->dilation_factor != -1 && !p->active){
      p->histr_shape[STEPS] 	*= p->dilation_factor;
      p->input_shape[STEPS] 	*= p->dilation_factor;
      p->non_active_strech 	*= p->dilation_factor;
    }
}


//scrubbed 3/15/2023
robin_protocol *new_robin_protocol(t_state *state, t_inputrec *ir, const real * gmx_restrict invmass){
  
    //alloc structure 
    robin_protocol *p = (robin_protocol*)malloc(sizeof*p);

    //be flashy
    rp_splash();

    //init parameters specified by input record and system state
    p->input_rec 		= ir;
    p->reporter 		= ir->rp_reporter;
    p->alternative_stride 	= ir->rp_alternative_stride;
    p->alternative_pause 	= ir->rp_alternative_pause;
    p->alternate_patterns	= p->alternative_stride == -1 ? false : true;
    p->dump_vels 		= ir->rp_dump_vels;
    p->dump_rvels 		= ir->rp_dump_rvels;
    p->balanced_diet 		= ir->rp_balanced_diet;
    p->threshold_type 		= ir->rp_threshold_type;
    p->threshold 		= ir->rp_threshold;
    p->dilation_factor 		= ir->rp_dilation_factor;
    p->histr_shape[STEPS] 	= ir->rp_sample_sz;
    p->input_shape[STEPS] 	= ir->rp_sample_sz;
    p->active 			= ir->robin_protocol == 1 ? true : false;
    p->domain 			= ir->rp_domain;
    p->dt_per_ps 		= 1.0/ir->delta_t;
    p->box[0] 			= state->box[XX][XX];
    p->box[1] 			= state->box[YY][YY];
    p->box[2] 			= state->box[ZZ][ZZ];
    p->natoms 			= state->natoms;
    p->inj 			= 0;
    p->Ngr 			= 0;
    p->non_active_strech 	= 1;
    p->inv_mass			= invmass;

    //init the feature set
    init_rp_fset(p,ir);

    //init groups
    init_rp_groups(p,ir);

    //init record keeping arrays
    init_rp_records(p,state);

    //init forecasting whenever doing an active run
    init_rp_forecasting(p,state,ir);

    //init non-active dump histograms p.r.n.
    init_rp_nad(p);

    //be verbose
    announce_new_rp(p);

    return p;
}


//scrubbed 3/15/2023
void init_robin_protocol(t_state *state, t_inputrec *inputrec,t_mdatoms *md){

    //allocate and initialize a new robin protocol     
    rp = new_robin_protocol(state,inputrec,md->invmass);

    //create directories for active/non-active dump data p.r.n.
    struct stat st = {0};

    if (rp->active && rp->dump_vels){
      if (stat("rp_active_dump", &st) == -1)
        mkdir("rp_active_dump", 0777);
    }
    else if (!rp->active && (rp->dump_vels || rp->dump_rvels)){
      if (stat("rp_nonactive_dump", &st) == -1)
        mkdir("rp_nonactive_dump", 0777);
    }
}


//scrubbed 3/15/2023
void free_robin_protocol(robin_protocol *p){ 

    //release allocations only present during active runs
    if (p->active){
      for (int i=0;i<p->Ngr;i++) if (p->sesh[i] != NULL){
          TF_DeleteGraph(p->graph[i]);
          TF_DeleteSession(p->sesh[i],p->status[i]);
          TF_DeleteStatus(p->status[i]);
          free(p->io_opts[i]);
        }

      for (int i=0;i<p->natoms;i++)
        free(p->forecast[i]);

      free(p->graph);
      free(p->sesh);
      free(p->status);
      free(p->io_opts);
      free(p->models);
      free(p->ml_tog);
      free(p->forecast);
    }

    //release all other allocations
    for (int i=0;i<p->natoms;i++){
      for (int j=0;j<p->histr_shape[STEPS];j++)
        free(p->history[i][j]);

      free(p->history[i]);
      free(p->rotation[i]);
      free(p->t0_magnitude[i]);
    }

    free(p->rotation); 
    free(p->history); 

    free(p); 
    p = NULL;

    //be verbose
    printf("\nFree Bird! Robin protocol released."
           " \"Cause I'm as free as a bird now\" (Lynyrd Skynyrd)\n\n");
}


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
 

                               &@@@@@@@@@@@@@@@@@@@(                  
                       @@@@@@,..........................@@            
                 ,@@@@......................................@         
             &@@*................*%@@@@@@(....................@       
          @@....../@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@      
       @@.&@@@@@/ @@@@@@@@@@@@......@@@@@@@@@@@@@@@............@      
    (@@        @@@@@@@@@@@..........@@@@@@@@@@@@@@@............@      
             @@@@@@@@@@@@@..........@@@@@@@@@@@@@@%...........@@@     
           .@@@@@@@@@@@@@@.........@@@@@@@@@@@@@@............@@@@@@   
          (@@@@@@@@@@@@@@..........@@@@@@@@@@@,............@@@@@@@@@  
          @@@@@@@@@@@@@@@..........@@@@@(...............@@@@@@@@@@@@  
          @@@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@& 
          @@@@@@@@@@@@@@...........................@@@@@@@@@@@@@@@@@  
          @@@@@@@@@@@@@@............................#@@@@@@@@@@@@@@@  
           @@@@@@@@@@@@(...........#@@@@..............@@@@@@@@@@@@@   
             @@@@@@@@@@............@@@@@@@@@@@..........@@@@@@@@@     
               @@@@@@@@............@@@@@@@@@@@@@@@#......@@@@@@       
                  @@@@............,@@@@@@@@@@@@@@@@@@@/....@@         
                     @.........@@@@@@@@@@@@@@@@@@@@@@@@@%@...@        
                    @,....&@@@@@@@@@@@@@@@@@@@@@@@.         @*@%      
                    *(             Protocol 5.0                .@   



         ^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^

             Written by Guy "Wayyne" Dayhoff II
	           @ the University of South Florida, 2020-2023

		           <version 5.0, March 2023>

         ^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^%^


\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

