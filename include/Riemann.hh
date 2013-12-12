#ifndef HEADER_HH
#define HEADER_HH

typedef int INTEGER;
typedef double REAL;

#include <iostream>
#include <fstream>  
#include <iomanip>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

const REAL PI=4*atan(1.0);
const REAL GAMMA=1.4;           /**< specific heat ratio of air */

//constants relating to GAMMA
const REAL G1=(GAMMA-1.0)/(2.0*GAMMA);
const REAL G2=(GAMMA+1.0)/(2.0*GAMMA);
const REAL G3=2.0*GAMMA/(GAMMA-1.0);
const REAL G4=2.0/(GAMMA-1.0);
const REAL G5=2.0/(GAMMA+1.0);
const REAL G6=(GAMMA-1.0)/(GAMMA+1.0);
const REAL G7=0.5*(GAMMA-1.0);
const REAL G8=1.0/GAMMA;
const REAL G9=GAMMA-1.0;

class SHOCK
{
protected:

  REAL xL,xC,xR;                /**< geometry */
  REAL rhoL,uL,pL,eL;           /**< left state */
  REAL rhoR,uR,pR,eR;           /**< right state */
  REAL rho_star,u_star,p_star;  /**< middle state */
  REAL rho_star_L,rho_star_R;   /*!<  */
  REAL cL,cR,c;                 /**< speed of sound */
  REAL c_star_L,c_star_R;       /*!< speed of sound */
  REAL sL,sR;                   /*!< shock speed */
  REAL sHL,sTL;                 /*!< rarefaction head and tail, left */
  REAL sHR,sTR;                 /*!< rafefaction head and tail, right*/
  REAL du;                      /**< uR-uL */
  void Pressure_U();
  void Display();
  void PREFUN(REAL&,REAL&,REAL,REAL,REAL,REAL);
  void STARTE(REAL&);
  void Density();
  void Sample(REAL,REAL,REAL&,REAL&,REAL&,REAL&);
public:
  SHOCK();
  SHOCK(REAL, REAL, REAL, REAL, REAL, REAL);
  ~SHOCK();
  void Set(REAL, REAL, REAL, REAL, REAL, REAL);
  void ExactRiemann(REAL,REAL,REAL);
  void LocalRiemann(REAL,REAL,REAL,REAL,REAL,REAL,REAL&,REAL&,REAL&,REAL&);

};

class NumericalSHOCK: public SHOCK
{
protected:
  INTEGER Scheme;               /*!< 0:first order; 1:MUSCL-third order; 
                                     2:WENO 3; 3:WENO 5*/
  INTEGER N;                    /*!< number of mesh cells */
  REAL dx;                      /*!< mesh step size */
  REAL dt;                      /*!< time step size */
  REAL *Mx;                     /*!< x coordinate of the mesh cell center */
  REAL *Mq1;                    /*!< density at the mesh cell center */
  REAL *Mq2;                    /*!< moment at the mesh cell center */
  REAL *Mq3;                    /*!< Total energy at the mesh cell center */
  REAL *Mq4;                    /*!< Internal enery at the mesh cell center */

  REAL *MF1;                    /*!< flux for density */
  REAL *MF2;                    /*!< flux for momentum */
  REAL *MF3;                    /*!< flux for enery */

  void MeshSet();               /*!< set the mesh */
  void Init();                  /*!< Initialise the flow */
  void Advance(REAL);           /*!< advance the solution */
  void Flux();                  /*!< compute the flux first-order*/
  void Flux_SecondOrderBackup();      /*!< compute the flux second-order */
  void Flux_MUSCL();            /*!< compute the flux second-order & third-order */
  void Flux_WENO3();            /*!< third order WENO scheme */
  void Flux_WENO5();            /*!< fifth order WENO scheme */
  void FaceFluxGudunov(REAL,REAL,REAL,REAL,REAL,REAL,REAL&,REAL&,REAL&); /*!< Gudunov's scheme */
  void FaceFluxHLL(REAL,REAL,REAL,REAL,REAL,REAL,REAL&,REAL&,REAL&); /*!< HLL scheme */
  void FaceFluxHLLC(REAL,REAL,REAL,REAL,REAL,REAL,REAL&,REAL&,REAL&); /*!< HLLC scheme */
  void FaceFluxHLLCBACKUP(REAL,REAL,REAL,REAL,REAL,REAL,REAL&,REAL&,REAL&); /*!< HLLC scheme */
public:
  NumericalSHOCK();             /*!< set the condition */
  NumericalSHOCK(REAL, REAL, REAL, REAL, REAL, REAL, INTEGER);
  ~NumericalSHOCK();
  void Out();
  void Solution(REAL);          /*!< get the solution */
};

//HLLC riemann 2d and 3d solver
void riemann_hllc_2d_(REAL wL[4], REAL wR[4], 
                      REAL &kas, REAL &kbs, REAL &Lab, 
                      REAL &ddt, REAL GI[4]);
void riemann_hllc_3d_(REAL wL[5], REAL wR[5], 
                      REAL n[3], REAL &Lab, 
                      REAL &ddt, REAL GI[5]);


#endif
