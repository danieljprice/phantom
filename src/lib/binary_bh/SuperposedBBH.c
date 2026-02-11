/*@@
   @file      SuperposedBBH.c
   @date      May 31, 2023
   @author    Luciano Combi
   @desc      This is the bare function that can be obtained from the Mathematica notebook.
   @enddesc 
@@*/

#include <math.h>
#include <stdlib.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define NDIM 4
#define X1 0
#define Y1 1
#define Z1 2
#define X2 3
#define Y2 4
#define Z2 5
#define VX1 6
#define VY1 7
#define VZ1 8
#define VX2 9
#define VY2 10
#define VZ2 11
#define AX1 12
#define AY1 13
#define AZ1 14
#define AX2 15
#define AY2 16
#define AZ2 17
#define M1T 18
#define M2T 19
#define AST_adjust_mass1 1.
#define AST_adjust_mass2 1.
#define AST_a1_buffer 1.e-4
#define AST_a2_buffer 1.e-4
#define AST_cutoff_floor 0.1

/* Function declarations */
extern void SuperposedBBH(const double *xx, double gcov[][NDIM], const double *traj_array);

/* Function definition */
/* Expressions obtained from SKS_NoAcc_3D_spinarb.nb */
extern void SuperposedBBH(const double *xx, double gcov[][NDIM], const double *traj_array)
{
  /* Mask */
  double    x =  xx[0];
  double    y =  xx[1];
  double    z =  xx[2];

  /* Superposition components*/
  double KS1[4][4];
  double KS2[4][4];
  double J1[4][4];
  double J2[4][4];

  /* Load trajectories */
  double xi1x = traj_array[X1]; 
  double xi1y = traj_array[Y1]; 
  double xi1z = traj_array[Z1];
  double xi2x = traj_array[X2];
  double xi2y = traj_array[Y2];
  double xi2z = traj_array[Z2];
  double v1x  = traj_array[VX1] + 1e-40;
  double v1y  = traj_array[VY1] + 1e-40;
  double v1z  = traj_array[VZ1] + 1e-40;
  double v2x =  traj_array[VX2] + 1e-40;
  double v2y =  traj_array[VY2] + 1e-40;
  double v2z =  traj_array[VZ2] + 1e-40;
  
  double v2  =  sqrt( v2x * v2x + v2y * v2y + v2z * v2z );
  double v1  =  sqrt( v1x * v1x + v1y * v1y + v1z * v1z ); 
  
  double a1x  = traj_array[AX1];
  double a1y  = traj_array[AY1];
  double a1z  = traj_array[AZ1];
  
  double a2x =  traj_array[AX2];
  double a2y =  traj_array[AY2];
  double a2z =  traj_array[AZ2];
  
  double m1_t = traj_array[M1T];
  double m2_t = traj_array[M2T];

  double a1_t = sqrt( a1x*a1x + a1y*a1y + a1z*a1z + 1e-40) ;
  double a2_t = sqrt( a2x*a2x + a2y*a2y + a2z*a2z + 1e-40) ;
 
  /* Load coordinates */  

   double oo1 = v1 * v1;
   double oo2 = oo1 * -1;
   double oo3 = 1 + oo2;
   double oo4 = sqrt(oo3);
   double oo5 = 1 / oo4;
   double oo6 = x * -1;
   double oo7 = oo6 + xi1x;
   double oo8 = v1x * oo7;
   double oo9 = y * -1;
   double oo10 = z * -1;
   double oo11 = v2 * v2;
   double oo12 = oo11 * -1;
   double oo13 = 1 + oo12;
   double oo14 = sqrt(oo13);
   double oo15 = 1 / oo14;
   double oo16 = oo6 + xi2x;
   double oo17 = v2x * oo16;
   double oo18 = xi1x * -1;
   double oo19 = 1 / oo1;
   double oo20 = -1 + oo4;
   double oo21 = xi1y * -1;
   double oo22 = xi1z * -1;
   double oo23 = xi2x * -1;
   double oo24 = 1 / oo11;
   double oo25 = -1 + oo14;
   double oo26 = xi2y * -1;
   double oo27 = xi2z * -1;
   double oo28 = xi1y * v1y;
   double oo29 = xi1z * v1z;
   double oo30 = v1y * (y * -1);
   double oo31 = v1z * (z * -1);
   double oo32 = oo28 + (oo29 + (oo30 + (oo31 + oo8)));
   double oo33 = xi2y * v2y;
   double oo34 = xi2z * v2z;
   double oo35 = v2y * (y * -1);
   double oo36 = v2z * (z * -1);
   double oo37 = oo17 + (oo33 + (oo34 + (oo35 + oo36)));
   double x0BH1 = (oo8 + ((oo9 + xi1y) * v1y + (oo10 + xi1z) * v1z)) * oo5;
   double x0BH2 = (oo17 + ((oo9 + xi2y) * v2y + (oo10 + xi2z) * v2z)) * oo15;
   double x1BH1 = (oo18 + x) - oo20 * (oo5 * (v1x * (((oo18 + x) * v1x + ((oo21 + y) * v1y + (oo22 + z) * v1z)) * oo19)));
   double x1BH2 = (oo23 + x) - oo24 * (oo25 * (v2x * (((oo23 + x) * v2x + ((oo26 + y) * v2y + (oo27 + z) * v2z)) * oo15)));
   double x2BH1 = oo21 + (oo20 * (oo32 * (oo5 * (v1y * oo19))) + y);
   double x2BH2 = oo26 + (oo24 * (oo25 * (oo37 * (v2y * oo15))) + y);
   double x3BH1 = oo22 + (oo20 * (oo32 * (oo5 * (v1z * oo19))) + z);
   double x3BH2 = oo27 + (oo24 * (oo25 * (oo37 * (v2z * oo15))) + z);

 
  /* Adjust mass */
  /* This is useful for reducing the effective mass of each BH */
  /* Adjust by hand to get the correct irreducible mass of the BH */
  double a1 = a1_t * AST_adjust_mass1;
  double m1 = m1_t * AST_adjust_mass1;
  double a2 = a2_t * AST_adjust_mass2;
  double m2 = m2_t * AST_adjust_mass2;
 
  //============================================// 
  // Regularize horizon and apply excision mask //
  //============================================//

  /* Define radius with respect to BH frame */
  double rBH1 = sqrt( x1BH1*x1BH1 + x2BH1*x2BH1 + x3BH1*x3BH1) ;
  double rBH2 = sqrt( x1BH2*x1BH2 + x2BH2*x2BH2 + x3BH2*x3BH2) ;

   /* Define radius cutoff */
  double rBH1_Cutoff = fabs(a1) * ( 1.0 + AST_a1_buffer) + AST_cutoff_floor ;
  double rBH2_Cutoff = fabs(a2) * ( 1.0 + AST_a2_buffer) + AST_cutoff_floor ;

  /* Apply excision */
  if ((rBH1) < rBH1_Cutoff) { if(x3BH1>0) {x3BH1 = rBH1_Cutoff;} else {x3BH1 = -1.0*rBH1_Cutoff;}}
  if ((rBH2) < rBH2_Cutoff) { if(x3BH2>0) {x3BH2 = rBH2_Cutoff;} else {x3BH2 = -1.0*rBH2_Cutoff;}}
 

  //=================//
  //     Metric      //
  //=================//
  double o1 = 1.4142135623730951;
  double o2 = 1 / o1;
  double o3 = a1x * a1x;
  double o4 = o3 * -1;
  double o5 = a1z * a1z;
  double o6 = o5 * -1;
  double o7 = a2x * a2x;
  double o8 = o7 * -1;
  double o9 = x1BH1 * x1BH1;
  double o10 = x2BH1 * x2BH1;
  double o11 = x3BH1 * x3BH1;
  double o12 = x1BH1 * a1x;
  double o13 = x2BH1 * a2x;
  double o14 = x3BH1 * a1z;
  double o15 = o12 + (o13 + o14);
  double o16 = o15 * o15;
  double o17 = o16 * 4;
  double o18 = o10 + (o11 + (o4 + (o6 + (o8 + o9))));
  double o19 = o18 * o18;
  double o20 = o17 + o19;
  double o21 = sqrt(o20);
  double o22 = o10 + (o11 + (o21 + (o4 + (o6 + (o8 + o9)))));
  double o23 = pow(o22, 1.5);
  double o24 = o22 * o22;
  double o25 = o24 * 0.25;
  double o26 = o16 + o25;
  double o27 = 1 / o26;
  double o28 = x2BH1 * a1z;
  double o29 = a2x * (x3BH1 * -1);
  double o30 = sqrt(o22);
  double o31 = 1 / o30;
  double o32 = o1 * (o15 * (o31 * a1x));
  double o33 = o30 * (x1BH1 * o2);
  double o34 = o28 + (o29 + (o32 + o33));
  double o35 = o22 * 0.5;
  double o36 = o3 + (o35 + (o5 + o7));
  double o37 = 1 / o36;
  double o38 = o2 * (o23 * (o27 * (o34 * (o37 * m1))));
  double o39 = a1z * (x1BH1 * -1);
  double o40 = x3BH1 * a1x;
  double o41 = o1 * (o15 * (o31 * a2x));
  double o42 = o30 * (x2BH1 * o2);
  double o43 = o39 + (o40 + (o41 + o42));
  double o44 = o2 * (o23 * (o27 * (o37 * (o43 * m1))));
  double o45 = x1BH1 * a2x;
  double o46 = a1x * (x2BH1 * -1);
  double o47 = o1 * (o15 * (o31 * a1z));
  double o48 = o30 * (x3BH1 * o2);
  double o49 = o45 + (o46 + (o47 + o48));
  double o50 = o2 * (o23 * (o27 * (o37 * (o49 * m1))));
  double o51 = o36 * o36;
  double o52 = 1 / o51;
  double o53 = o2 * (o23 * (o27 * (o34 * (o43 * (o52 * m1)))));
  double o54 = o2 * (o23 * (o27 * (o34 * (o49 * (o52 * m1)))));
  double o55 = o2 * (o23 * (o27 * (o43 * (o49 * (o52 * m1)))));
  double o56 = a2y * a2y;
  double o57 = o56 * -1;
  double o58 = a2z * a2z;
  double o59 = o58 * -1;
  double o60 = x1BH2 * x1BH2;
  double o61 = x2BH2 * x2BH2;
  double o62 = x3BH2 * x3BH2;
  double o63 = x1BH2 * a2x;
  double o64 = x2BH2 * a2y;
  double o65 = x3BH2 * a2z;
  double o66 = o63 + (o64 + o65);
  double o67 = o66 * o66;
  double o68 = o67 * 4;
  double o69 = o57 + (o59 + (o60 + (o61 + (o62 + o8))));
  double o70 = o69 * o69;
  double o71 = o68 + o70;
  double o72 = sqrt(o71);
  double o73 = o57 + (o59 + (o60 + (o61 + (o62 + (o72 + o8)))));
  double o74 = pow(o73, 1.5);
  double o75 = o73 * o73;
  double o76 = o75 * 0.25;
  double o77 = o67 + o76;
  double o78 = 1 / o77;
  double o79 = x2BH2 * a2z;
  double o80 = a2y * (x3BH2 * -1);
  double o81 = sqrt(o73);
  double o82 = 1 / o81;
  double o83 = o1 * (o66 * (o82 * a2x));
  double o84 = o81 * (x1BH2 * o2);
  double o85 = o79 + (o80 + (o83 + o84));
  double o86 = o73 * 0.5;
  double o87 = o56 + (o58 + (o7 + o86));
  double o88 = 1 / o87;
  double o89 = o2 * (o74 * (o78 * (o85 * (o88 * m2))));
  double o90 = a2z * (x1BH2 * -1);
  double o91 = x3BH2 * a2x;
  double o92 = o1 * (o66 * (o82 * a2y));
  double o93 = o81 * (x2BH2 * o2);
  double o94 = o90 + (o91 + (o92 + o93));
  double o95 = o2 * (o74 * (o78 * (o88 * (o94 * m2))));
  double o96 = x1BH2 * a2y;
  double o97 = a2x * (x2BH2 * -1);
  double o98 = o1 * (o66 * (o82 * a2z));
  double o99 = o81 * (x3BH2 * o2);
  double o100 = o96 + (o97 + (o98 + o99));
  double o101 = o100 * (o2 * (o74 * (o78 * (o88 * m2))));
  double o102 = o87 * o87;
  double o103 = 1 / o102;
  double o104 = o103 * (o2 * (o74 * (o78 * (o85 * (o94 * m2)))));
  double o105 = o100 * (o103 * (o2 * (o74 * (o78 * (o85 * m2)))));
  double o106 = o100 * (o103 * (o2 * (o74 * (o78 * (o94 * m2)))));
  double o107 = v1 * v1;
  double o108 = o107 * -1;
  double o109 = 1 + o108;
  double o110 = sqrt(o109);
  double o111 = 1 / o110;
  double o112 = o111 * (v1x * -1);
  double o113 = o111 * (v1y * -1);
  double o114 = o111 * (v1z * -1);
  double o115 = 1 / o107;
  double o116 = -1 + o111;
  double o117 = o116 * (v1x * (v1y * o115));
  double o118 = o116 * (v1x * (v1z * o115));
  double o119 = o116 * (v1y * (v1z * o115));
  double o120 = v2 * v2;
  double o121 = o120 * -1;
  double o122 = 1 + o121;
  double o123 = sqrt(o122);
  double o124 = 1 / o123;
  double o125 = o124 * (v2x * -1);
  double o126 = o124 * (v2y * -1);
  double o127 = o124 * (v2z * -1);
  double o128 = 1 / o120;
  double o129 = -1 + o124;
  double o130 = o129 * (v2x * (v2y * o128));
  double o131 = o129 * (v2x * (v2z * o128));
  double o132 = o129 * (v2y * (v2z * o128));
  KS1[0][0] = o2 * (o23 * (o27 * m1));
  KS1[0][1] = o38;
  KS1[0][2] = o44;
  KS1[0][3] = o50;
  KS1[1][0] = o38;
  KS1[1][1] = o2 * (o23 * (o27 * ((o34 * o34) * (o52 * m1))));
  KS1[1][2] = o53;
  KS1[1][3] = o54;
  KS1[2][0] = o44;
  KS1[2][1] = o53;
  KS1[2][2] = o2 * (o23 * (o27 * ((o43 * o43) * (o52 * m1))));
  KS1[2][3] = o55;
  KS1[3][0] = o50;
  KS1[3][1] = o54;
  KS1[3][2] = o55;
  KS1[3][3] = o2 * (o23 * (o27 * ((o49 * o49) * (o52 * m1))));
  KS2[0][0] = o2 * (o74 * (o78 * m2));
  KS2[0][1] = o89;
  KS2[0][2] = o95;
  KS2[0][3] = o101;
  KS2[1][0] = o89;
  KS2[1][1] = o103 * (o2 * (o74 * (o78 * ((o85 * o85) * m2))));
  KS2[1][2] = o104;
  KS2[1][3] = o105;
  KS2[2][0] = o95;
  KS2[2][1] = o104;
  KS2[2][2] = o103 * (o2 * (o74 * (o78 * ((o94 * o94) * m2))));
  KS2[2][3] = o106;
  KS2[3][0] = o101;
  KS2[3][1] = o105;
  KS2[3][2] = o106;
  KS2[3][3] = (o100 * o100) * (o103 * (o2 * (o74 * (o78 * m2))));
  J1[0][0] = o111;
  J1[0][1] = o112;
  J1[0][2] = o113;
  J1[0][3] = o114;
  J1[1][0] = o112;
  J1[1][1] = 1 + o116 * ((v1x * v1x) * o115);
  J1[1][2] = o117;
  J1[1][3] = o118;
  J1[2][0] = o113;
  J1[2][1] = o117;
  J1[2][2] = 1 + o116 * ((v1y * v1y) * o115);
  J1[2][3] = o119;
  J1[3][0] = o114;
  J1[3][1] = o118;
  J1[3][2] = o119;
  J1[3][3] = 1 + o116 * ((v1z * v1z) * o115);
  J2[0][0] = o124;
  J2[0][1] = o125;
  J2[0][2] = o126;
  J2[0][3] = o127;
  J2[1][0] = o125;
  J2[1][1] = 1 + o129 * ((v2x * v2x) * o128);
  J2[1][2] = o130;
  J2[1][3] = o131;
  J2[2][0] = o126;
  J2[2][1] = o130;
  J2[2][2] = 1 + o129 * ((v2y * v2y) * o128);
  J2[2][3] = o132;
  J2[3][0] = o127;
  J2[3][1] = o131;
  J2[3][2] = o132;
  J2[3][3] = 1 + o129 * ((v2z * v2z) * o128);
  /* Initialize the flat part */
  double eta[4][4] = {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  for (int i=0; i < 4; i++ ){ for (int j=0; j < 4; j++ ){ gcov[i][j] = eta[i][j]; }}

  /* Load symmetric gcov (from chatGPT3)*/
  for (int i = 0; i < 4; ++i) {
      for (int j = i; j < 4; ++j) {
  
          double sum = 0.0;
          for (int m = 0; m < 4; ++m) {
              double term1 = J2[m][i];
              double term2 = J1[m][i];
  
              for (int n = 0; n < 4; ++n) {
                  double term3 = J2[n][j];
                  double term4 = J1[n][j];
  
                  sum += (term1 * term3 * KS2[m][n] + term2 * term4 * KS1[m][n]);
              }
          }
  
          gcov[i][j] += sum;
          gcov[j][i] = gcov[i][j];
      }
  }

  return;
}

