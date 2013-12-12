#include "Riemann.hh"

using namespace std;

NumericalSHOCK::NumericalSHOCK()
{
  //clear the variables
  xL=xC=xR=0.0;
  rhoL=uL=pL=rhoR=uR=pR=0.0;
  rho_star=u_star=p_star=0.0;
  rho_star_L=rho_star_R=0.0;
  cL=cR=c=0.0;
  c_star_L=c_star_R=0.0;
  sL=sR=sHL=sTL=sHR=sTR=0.0;


  //clear the arrays
  Mx=0;
  Mq1=0;  Mq2=0; Mq3=0; Mq4=0;
  MF1=0;  MF2=0; MF3=0; 
  
  Scheme=0;
}



NumericalSHOCK::NumericalSHOCK(REAL qL1, REAL qL2, REAL qL3, REAL qR1, REAL qR2, REAL qR3,
                               INTEGER scheme)
{
  //clear the variables
  xL=xC=xR=0.0;
  rhoL=uL=pL=rhoR=uR=pR=0.0;
  rho_star=u_star=p_star=0.0;
  rho_star_L=rho_star_R=0.0;
  cL=cR=c=0.0;
  c_star_L=c_star_R=0.0;
  sL=sR=sHL=sTL=sHR=sTR=0.0;

  //clear the arrays
  Mx=0;
  Mq1=0;  Mq2=0; Mq3=0; Mq4=0;
  MF1=0;  MF2=0; MF3=0; 

  //set the left state
  rhoL=qL1;
  uL=qL2;
  pL=qL3;
  eL=pL/G9/rhoL;

  //set the right state
  rhoR=qR1;
  uR=qR2;
  pR=qR3;
  eR=pR/G9/rhoR;

  //compute speed of sound
  cL=sqrt(GAMMA*pL/rhoL);
  cR=sqrt(GAMMA*pR/rhoR);

  //compute du=uR-uL;
  du=uR-uL;

  //scheme
  Scheme=scheme;
}

//! the destructor
/*! 
  
 */
NumericalSHOCK::~NumericalSHOCK()
{
  if(Mx) delete[] Mx; Mx=0;

  if(Mq1) delete[] Mq1; Mq1=0; 
  if(Mq2) delete[] Mq2; Mq2=0; 
  if(Mq3) delete[] Mq3; Mq3=0; 
  if(Mq4) delete[] Mq4; Mq4=0; 

  if(MF1) delete[] MF1; MF1=0; 
  if(MF2) delete[] MF2; MF2=0; 
  if(MF3) delete[] MF3; MF3=0; 

  cout<<"Done in NumericalSHOCK"<<endl;
} 

void NumericalSHOCK::MeshSet()
{
  N=1001;
  dx=1.0/(N-1);

  if(Mx==NULL) Mx =new REAL[N]; 
  
  for(INTEGER i=0;i<N;i++)
    Mx[i]=i*dx-0.5;


  if(Mq1==NULL) Mq1=new REAL[N]; 
  if(Mq2==NULL) Mq2=new REAL[N]; 
  if(Mq3==NULL) Mq3=new REAL[N]; 
  if(Mq4==NULL) Mq4=new REAL[N]; 

  if(MF1==NULL) MF1=new REAL[N-1]; 
  if(MF2==NULL) MF2=new REAL[N-1]; 
  if(MF3==NULL) MF3=new REAL[N-1]; 

  for(INTEGER i=0;i<N-1;i++)
    {
      MF1[i]=0;
      MF2[i]=0;
      MF3[i]=0;
    }


}

void NumericalSHOCK::Init()
{
  for(INTEGER i=0;i<N;i++){
    if(Mx[i]<=0.0){
      Mq1[i]=rhoL;
      Mq2[i]=rhoL*uL;
      Mq3[i]=pL/G9+0.5*rhoL*uL*uL;
      Mq4[i]=eL;
    }
    else{
      Mq1[i]=rhoR;
      Mq2[i]=rhoR*uR;
      Mq3[i]=pR/G9+0.5*rhoR*uR*uR;
      Mq4[i]=eR;     
    }
  }
}


void NumericalSHOCK::Out()
{
  ofstream ft;
  if(Scheme==0)  ft.open("Sod_X_Numerical-FirstOrder.dat");
  if(Scheme==1)  ft.open("Sod_X_Numerical-MUSCL.dat");
  if(Scheme==2)  ft.open("Sod_X_Numerical-WENO3.dat");
  if(Scheme==3)  ft.open("Sod_X_Numerical-WENO5.dat");
  for(INTEGER i=0;i<N;i++){
    ft<<Mx[i]<<" "
      <<Mq1[i]<<" "
      <<Mq2[i]/Mq1[i]<<" "
      <<G9*(Mq3[i]-0.5*Mq2[i]*Mq2[i]/Mq1[i])<<" "
      <<Mq4[i]<<endl;
  }
  ft.close();

}


void NumericalSHOCK::Flux()
{
  //left boundary condition - transmissive
  Mq1[0]=Mq1[1];
  Mq2[0]=Mq2[1];
  Mq3[0]=Mq3[1];
  //right boundary condition - transmissive
  Mq1[N-1]=Mq1[N-2];
  Mq2[N-1]=Mq2[N-2];
  Mq3[N-1]=Mq3[N-2];


  REAL qL1,qL2,qL3;             /*!< Left side, density,velocity,pressure */
  REAL qR1,qR2,qR3;             /*!< Right side, density,velocity,pressure */

  //loop the faces 
  for(INTEGER i=0;i<N-1;i++)
    {
      qL1=Mq1[i];   
      qL2=Mq2[i]/qL1;  
      qL3=G9*(Mq3[i]-0.5*qL1*qL2*qL2);

      qR1=Mq1[i+1];  
      qR2=Mq2[i+1]/qR1;  
      qR3=G9*(Mq3[i+1]-0.5*qR1*qR2*qR2);

      //FaceFluxGudunov(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
      FaceFluxHLLC(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
    }
}


void NumericalSHOCK::Flux_SecondOrderBackup()
{
  REAL const epsi=1.0e-10;
  //left boundary condition - transmissive
  Mq1[0]=Mq1[1];
  Mq2[0]=Mq2[1];
  Mq3[0]=Mq3[1];
  //right boundary condition - transmissive
  Mq1[N-1]=Mq1[N-2];
  Mq2[N-1]=Mq2[N-2];
  Mq3[N-1]=Mq3[N-2];


  REAL rhoA,uA,pA;
  REAL rhoB,uB,pB;
  REAL rhoC,uC,pC;
  REAL rhoD,uD,pD;

  REAL qL1,qL2,qL3;             /*!< Left side, density,velocity,pressure */
  REAL qR1,qR2,qR3;             /*!< Right side, density,velocity,pressure */

  //loop the faces 
  for(INTEGER i=0;i<N-1;i++)
    {

      rhoB=Mq1[i];
      uB=Mq2[i]/rhoB;
      pB=G9*(Mq3[i]-0.5*rhoB*uB*uB);

      rhoC=Mq1[i+1];
      uC=Mq2[i+1]/rhoC;
      pC=G9*(Mq3[i+1]-0.5*rhoC*uC*uC);

      if(i==0){
        //left bound
        rhoA=rhoB;
        uA=uB;
        pA=pB;

        rhoD=Mq1[i+2];
        uD=Mq2[i+2]/rhoD;
        pD=G9*(Mq3[i+2]-0.5*rhoD*uD*uD);
      }
      else if(i==N-2){
        //right bound
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]*0.5*rhoA*uA*uA);

        rhoD=rhoC;
        uD=uC;
        pD=pC;

      }
      else{
        //inner cell
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]*0.5*rhoA*uA*uA);
      }

      REAL aR1=rhoD-rhoC;
      REAL aR2=uD-uC;
      REAL aR3=pD-pC;

      REAL bR1=rhoC-rhoB;
      REAL bR2=uC-uB;
      REAL bR3=pC-pB;

      REAL aL1=bR1;
      REAL aL2=bR2;
      REAL aL3=bR3;

      REAL bL1=rhoB-rhoA;
      REAL bL2=uB-uA;
      REAL bL3=pB-pA;

      //k=0 
      // REAL dR1=(aR1*(bR1*bR1+epsi)+bR1*(aR1*aR1+epsi))/(aR1*aR1+bR1*bR1+2*epsi);
      // REAL dR2=(aR2*(bR2*bR2+epsi)+bR2*(aR2*aR2+epsi))/(aR2*aR2+bR2*bR2+2*epsi);      
      // REAL dR3=(aR3*(bR3*bR3+epsi)+bR3*(aR3*aR3+epsi))/(aR3*aR3+bR3*bR3+2*epsi);

      // REAL dL1=(aL1*(bL1*bL1+epsi)+bL1*(aL1*aL1+epsi))/(aL1*aL1+bL1*bL1+2*epsi);
      // REAL dL2=(aL2*(bL2*bL2+epsi)+bL2*(aL2*aL2+epsi))/(aL2*aL2+bL2*bL2+2*epsi);      
      // REAL dL3=(aL3*(bL3*bL3+epsi)+bL3*(aL3*aL3+epsi))/(aL3*aL3+bL3*bL3+2*epsi);

      //k=1/3
      REAL dR1=(bR1*(2*aR1*aR1+epsi)+aR1*(bR1*bR1+2*epsi))/(2*aR1*aR1+2*bR1*bR1-aR1*bR1+3*epsi);
      REAL dR2=(bR2*(2*aR2*aR2+epsi)+aR2*(bR2*bR2+2*epsi))/(2*aR2*aR2+2*bR2*bR2-aR2*bR2+3*epsi);
      REAL dR3=(bR3*(2*aR3*aR3+epsi)+aR3*(bR3*bR3+2*epsi))/(2*aR3*aR3+2*bR3*bR3-aR3*bR3+3*epsi);

      REAL dL1=(bL1*(2*aL1*aL1+epsi)+aL1*(bL1*bL1+2*epsi))/(2*aL1*aL1+2*bL1*bL1-aL1*bL1+3*epsi);
      REAL dL2=(bL2*(2*aL2*aL2+epsi)+aL2*(bL2*bL2+2*epsi))/(2*aL2*aL2+2*bL2*bL2-aL2*bL2+3*epsi);
      REAL dL3=(bL3*(2*aL3*aL3+epsi)+aL3*(bL3*bL3+2*epsi))/(2*aL3*aL3+2*bL3*bL3-aL3*bL3+3*epsi);

      

      qL1=Mq1[i];   
      qL2=Mq2[i]/qL1;  
      qL3=G9*(Mq3[i]-0.5*qL1*qL2*qL2);

        
      qR1=Mq1[i+1];  
      qR2=Mq2[i+1]/qR1;  
      qR3=G9*(Mq3[i+1]-0.5*qR1*qR2*qR2);

      
      qL1+=0.5*dL1;
      qL2+=0.5*dL2;
      qL3+=0.5*dL3;

      qR1-=0.5*dR1;
      qR2-=0.5*dR2;
      qR3-=0.5*dR3;
      
      //FaceFluxGudunov(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
      FaceFluxHLLC(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);

    }
}


void NumericalSHOCK::Flux_MUSCL()
{
  REAL const epsi=1.0e-12;
  //left boundary condition - transmissive
  Mq1[0]=Mq1[1];
  Mq2[0]=Mq2[1];
  Mq3[0]=Mq3[1];

  //right boundary condition - transmissive
  Mq1[N-1]=Mq1[N-2];
  Mq2[N-1]=Mq2[N-2];
  Mq3[N-1]=Mq3[N-2];

  REAL rhoA,uA,pA;
  REAL rhoB,uB,pB;
  REAL rhoC,uC,pC;
  REAL rhoD,uD,pD;

  REAL qL1,qL2,qL3;             /*!< Left side, density,velocity,pressure */
  REAL qR1,qR2,qR3;             /*!< Right side, density,velocity,pressure */

  //loop the faces 
  for(INTEGER i=0;i<N-1;i++)
    {

      rhoB=Mq1[i];
      uB=Mq2[i]/rhoB;
      pB=G9*(Mq3[i]-0.5*rhoB*uB*uB);

      rhoC=Mq1[i+1];
      uC=Mq2[i+1]/rhoC;
      pC=G9*(Mq3[i+1]-0.5*rhoC*uC*uC);

      if(i==0){
        //left bound
        rhoA=rhoB;
        uA=uB;
        pA=pB;

        rhoD=Mq1[i+2];
        uD=Mq2[i+2]/rhoD;
        pD=G9*(Mq3[i+2]-0.5*rhoD*uD*uD);
      }
      else if(i==N-2){
        //right bound
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]*0.5*rhoA*uA*uA);

        rhoD=rhoC;
        uD=uC;
        pD=pC;

      }
      else{
        //inner cell
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]-0.5*rhoA*uA*uA);

        rhoD=Mq1[i+2];
        uD=Mq2[i+2]/rhoD;
        pD=G9*(Mq3[i+2]-0.5*rhoD*uD*uD);
      }

      REAL aR1=rhoD-rhoC;
      REAL aR2=uD-uC;
      REAL aR3=pD-pC;

      REAL bR1=rhoC-rhoB;
      REAL bR2=uC-uB;
      REAL bR3=pC-pB;

      REAL aL1=bR1;
      REAL aL2=bR2;
      REAL aL3=bR3;

      REAL bL1=rhoB-rhoA;
      REAL bL2=uB-uA;
      REAL bL3=pB-pA;

      //Limiter
      REAL sL1=max(0.0,(2*bL1*aL1+epsi)/(bL1*bL1+aL1*aL1+epsi));
      REAL sL2=max(0.0,(2*bL2*aL2+epsi)/(bL2*bL2+aL2*aL2+epsi));
      REAL sL3=max(0.0,(2*bL3*aL3+epsi)/(bL3*bL3+aL3*aL3+epsi));

      REAL sR1=max(0.0,(2*aR1*aL1+epsi)/(aR1*aR1+aL1*aL1+epsi));
      REAL sR2=max(0.0,(2*aR2*aL2+epsi)/(aR2*aR2+aL2*aL2+epsi));
      REAL sR3=max(0.0,(2*aR3*aL3+epsi)/(aR3*aR3+aL3*aL3+epsi));

      REAL k=1.0/3.0;

      //REAL dL1=((1-k*sL1)*bL1+(1+k*sL1)*aL1)*sL1/2.0;
      //REAL dL2=((1-k*sL2)*bL2+(1+k*sL2)*aL2)*sL2/2.0;      
      //REAL dL3=((1-k*sL3)*bL3+(1+k*sL3)*aL3)*sL3/2.0;      

      //REAL dR1=((1-k*sR1)*aR1+(1+k*sR1)*aL1)*sR1/2.0;
      //REAL dR2=((1-k*sR2)*aR2+(1+k*sR2)*aL2)*sR1/2.0;
      //REAL dR3=((1-k*sR3)*aR3+(1+k*sR3)*aL3)*sR1/2.0;

       REAL dR1=(bR1*(2*aR1*aR1+epsi)+aR1*(bR1*bR1+2*epsi))/(2*aR1*aR1+2*bR1*bR1-aR1*bR1+3*epsi);
       REAL dR2=(bR2*(2*aR2*aR2+epsi)+aR2*(bR2*bR2+2*epsi))/(2*aR2*aR2+2*bR2*bR2-aR2*bR2+3*epsi);
       REAL dR3=(bR3*(2*aR3*aR3+epsi)+aR3*(bR3*bR3+2*epsi))/(2*aR3*aR3+2*bR3*bR3-aR3*bR3+3*epsi);

       REAL dL1=(bL1*(2*aL1*aL1+epsi)+aL1*(bL1*bL1+2*epsi))/(2*aL1*aL1+2*bL1*bL1-aL1*bL1+3*epsi);
       REAL dL2=(bL2*(2*aL2*aL2+epsi)+aL2*(bL2*bL2+2*epsi))/(2*aL2*aL2+2*bL2*bL2-aL2*bL2+3*epsi);
       REAL dL3=(bL3*(2*aL3*aL3+epsi)+aL3*(bL3*bL3+2*epsi))/(2*aL3*aL3+2*bL3*bL3-aL3*bL3+3*epsi);


      // qL1=Mq1[i];   
      // qL2=Mq2[i]/qL1;  
      // qL3=G9*(Mq3[i]-0.5*qL1*qL2*qL2);

        
      // qR1=Mq1[i+1];  
      // qR2=Mq2[i+1]/qR1;  
      // qR3=G9*(Mq3[i+1]-0.5*qR1*qR2*qR2);

      qL1=rhoB;
      qL2=uB;
      qL3=pB;

      qR1=rhoC;
      qR2=uC;
      qR3=pC;

      
      qL1+=0.5*dL1;
      qL2+=0.5*dL2;
      qL3+=0.5*dL3;

      qR1-=0.5*dR1;
      qR2-=0.5*dR2;
      qR3-=0.5*dR3;
      
      //FaceFluxGudunov(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
      FaceFluxHLLC(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);

    }
}

void NumericalSHOCK::Flux_WENO3()
{
  REAL const epsi=1.0e-6;

  //left boundary condition - transmissive
  Mq1[0]=Mq1[1];
  Mq2[0]=Mq2[1];
  Mq3[0]=Mq3[1];
  //right boundary condition - transmissive
  Mq1[N-1]=Mq1[N-2];
  Mq2[N-1]=Mq2[N-2];
  Mq3[N-1]=Mq3[N-2];

  REAL rhoA,uA,pA;
  REAL rhoB,uB,pB;
  REAL rhoC,uC,pC;
  REAL rhoD,uD,pD;

  REAL qL1,qL2,qL3;             /*!< Left side, density,velocity,pressure */
  REAL qR1,qR2,qR3;             /*!< Right side, density,velocity,pressure */

  REAL ftmp;
  
  //loop the faces 
  for(INTEGER i=0;i<N-1;i++)
    {

      rhoB=Mq1[i];
      uB=Mq2[i]/rhoB;
      pB=G9*(Mq3[i]-0.5*rhoB*uB*uB);

      rhoC=Mq1[i+1];
      uC=Mq2[i+1]/rhoC;
      pC=G9*(Mq3[i+1]-0.5*rhoC*uC*uC);

      if(i==0){
        //left bound
        rhoA=rhoB;
        uA=uB;
        pA=pB;

        rhoD=Mq1[i+2];
        uD=Mq2[i+2]/rhoD;
        pD=G9*(Mq3[i+2]-0.5*rhoD*uD*uD);
      }
      else if(i==N-2){
        //right bound
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]-0.5*rhoA*uA*uA);

        rhoD=rhoC;
        uD=uC;
        pD=pC;

      }
      else{
        //inner cell
        rhoA=Mq1[i-1];
        uA=Mq2[i-1]/rhoA;
        pA=G9*(Mq3[i-1]-0.5*rhoA*uA*uA);

        rhoD=Mq1[i+2];
        uD=Mq2[i+2]/rhoD;
        pD=G9*(Mq3[i+2]-0.5*rhoD*uD*uD);
      }

      //left side qL
      /*
      {
        REAL s1rho=(rhoB-rhoA)*(rhoB-rhoA);
        REAL s1u  =(uB-uA)*(uB-uA);
        REAL s1p  =(pB-pA)*(pB-pA);

        REAL s2rho=(rhoB-rhoC)*(rhoB-rhoC);
        REAL s2u  =(uB-uC)*(uB-uC);
        REAL s2p  =(pB-pC)*(pB-pC);

        REAL const C1=1.0/3.0;
        REAL const C2=2.0/3.0;

        REAL alpha1rho=C1/(s1rho+epsi)/(s1rho+epsi);
        REAL alpha1u  =C1/(s1u+epsi)/(s1u+epsi);
        REAL alpha1p  =C1/(s1p+epsi)/(s2p+epsi);

        REAL alpha2rho=C2/(s2rho+epsi)/(s2rho+epsi);
        REAL alpha2u  =C2/(s2u+epsi)/(s2u+epsi);
        REAL alpha2p  =C2/(s2p+epsi)/(s2p+epsi);

        REAL w1rho = alpha1rho/(alpha1rho+alpha2rho);
        REAL w1u   = alpha1u/(alpha1u+alpha2u);
        REAL w1p   = alpha1p/(alpha1p+alpha2p);

        REAL w2rho = alpha2rho/(alpha1rho+alpha2rho);
        REAL w2u   = alpha2u/(alpha1u+alpha2u);
        REAL w2p   = alpha2p/(alpha1p+alpha2p);

        qL1= w1rho*(-0.5*rhoA+1.5*rhoB)+w2rho*(0.5*rhoB+0.5*rhoC);  //rho
        qL2= w1u*(-0.5*uA+1.5*uB)+w2u*(0.5*uB+0.5*uC);              //u
        qL3= w1p*(-0.5*pA+1.5*pB)+w2p*(0.5*pB+0.5*pC);              //p
      }
      */

      REAL C1=1.0/3.0;
      REAL C2=2.0/3.0;

      ftmp=rhoB-rhoA; REAL s1rho=ftmp*ftmp;
      ftmp=uB-uA;     REAL s1u  =ftmp*ftmp;
      ftmp=pB-pA;     REAL s1p  =ftmp*ftmp;

      ftmp=rhoB-rhoC; REAL s2rho=ftmp*ftmp;
      ftmp=uB-uC;     REAL s2u  =ftmp*ftmp;
      ftmp=pB-pC;     REAL s2p  =ftmp*ftmp;

      ftmp=s1rho+epsi; REAL alpha1rho=C1/ftmp/ftmp;
      ftmp=s1u+epsi;   REAL alpha1u  =C1/ftmp/ftmp;
      ftmp=s1p+epsi;   REAL alpha1p  =C1/ftmp/ftmp;

      ftmp=s2rho+epsi;  REAL alpha2rho=C2/ftmp/ftmp;
      ftmp=s2u+epsi;    REAL alpha2u  =C2/ftmp/ftmp;
      ftmp=s2p+epsi;    REAL alpha2p  =C2/ftmp/ftmp;

      REAL w1rho = alpha1rho/(alpha1rho+alpha2rho);
      REAL w1u   = alpha1u/(alpha1u+alpha2u);
      REAL w1p   = alpha1p/(alpha1p+alpha2p);

      REAL w2rho = 1.0-w1rho;
      REAL w2u   = 1.0-w1u;
      REAL w2p   = 1.0-w1p;

      qL1= (w1rho*(3*rhoB-rhoA)+w2rho*(rhoB+rhoC))*0.5;  //rho
      qL2= (w1u*(3*uB-uA)+w2u*(uB+uC))*0.5;              //u
      qL3= (w1p*(3*pB-pA)+w2p*(pB+pC))*0.5;              //p
      
      //right side qR
      /*
      {
        REAL s1rho = (rhoB-rhoC)*(rhoB-rhoC);
        REAL s1u   = (uB-uC)*(uB-uC);
        REAL s1p   = (pB-pC)*(pB-pC);

        REAL s2rho = (rhoC-rhoD)*(rhoC-rhoD);
        REAL s2u   = (uC-uD)*(uC-uD);
        REAL s2p   = (pC-pD)*(pC-pD);

        REAL const C1=2.0/3.0;
        REAL const C2=1.0/3.0;

        REAL alpha1rho=C1/(s1rho+epsi)/(s1rho+epsi);
        REAL alpha1u  =C1/(s1u+epsi)/(s1u+epsi);
        REAL alpha1p  =C1/(s1p+epsi)/(s2p+epsi);

        REAL alpha2rho=C2/(s2rho+epsi)/(s2rho+epsi);
        REAL alpha2u  =C2/(s2u+epsi)/(s2u+epsi);
        REAL alpha2p  =C2/(s2p+epsi)/(s2p+epsi);

        REAL w1rho = alpha1rho/(alpha1rho+alpha2rho);
        REAL w1u   = alpha1u/(alpha1u+alpha2u);
        REAL w1p   = alpha1p/(alpha1p+alpha2p);

        REAL w2rho = alpha2rho/(alpha1rho+alpha2rho);
        REAL w2u   = alpha2u/(alpha1u+alpha2u);
        REAL w2p   = alpha2p/(alpha1p+alpha2p);

        qR1= w1rho*(0.5*rhoB+0.5*rhoC)+w2rho*(1.5*rhoC-0.5*rhoD);
        qR2= w1u*(0.5*uB+0.5*uC)+w2u*(1.5*uC-0.5*uD);
        qR3= w1p*(0.5*pB+0.5*pC)+w2p*(1.5*pC-0.5*pD);
      }      
      */
     
      C1=2.0/3.0;
      C2=1.0/3.0;

      s1rho = s2rho;
      s1u   = s2u;
      s1p   = s2p;

      ftmp=rhoC-rhoD;  s2rho = ftmp*ftmp;
      ftmp=uC-uD;      s2u   = ftmp*ftmp;
      ftmp=pC-pD;      s2p   = ftmp*ftmp;

      ftmp=s1rho+epsi; alpha1rho=C1/ftmp/ftmp;
      ftmp=s1u+epsi;   alpha1u  =C1/ftmp/ftmp;
      ftmp=s1p+epsi;   alpha1p  =C1/ftmp/ftmp;

      ftmp=s2rho+epsi;  alpha2rho=C2/ftmp/ftmp;
      ftmp=s2u+epsi;    alpha2u  =C2/ftmp/ftmp;
      ftmp=s2p+epsi;    alpha2p  =C2/ftmp/ftmp;

      w1rho = alpha1rho/(alpha1rho+alpha2rho);
      w1u   = alpha1u/(alpha1u+alpha2u);
      w1p   = alpha1p/(alpha1p+alpha2p);

      w2rho = 1.0-w1rho;
      w2u   = 1.0-w1u;
      w2p   = 1.0-w1p;

      qR1= (w1rho*(rhoB+rhoC)+w2rho*(3*rhoC-rhoD))*0.5;
      qR2= (w1u*(uB+uC)+w2u*(3*uC-uD))*0.5;
      qR3= (w1p*(pB+pC)+w2p*(3*pC-pD))*0.5;
    
      //FaceFluxGudunov(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
      FaceFluxHLLC(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
    }

 }

void NumericalSHOCK::Flux_WENO5()
{
  REAL const epsi=1.0e-6;
  //left boundary condition - transmissive
  Mq1[0]=Mq1[1];
  Mq2[0]=Mq2[1];
  Mq3[0]=Mq3[1];
  //right boundary condition - transmissive
  Mq1[N-1]=Mq1[N-2];
  Mq2[N-1]=Mq2[N-2];
  Mq3[N-1]=Mq3[N-2];

  REAL rhoA,uA,pA;   //point: i-2
  REAL rhoB,uB,pB;   //point: i-1
  REAL rhoC,uC,pC;   //point: i
  REAL rhoD,uD,pD;   //point: i+1
  REAL rhoE,uE,pE;   //point: i+2
  REAL rhoF,uF,pF;   //point: i+3

  REAL qL1,qL2,qL3;             /*!< Left side, density,velocity,pressure */
  REAL qR1,qR2,qR3;             /*!< Right side, density,velocity,pressure */

  REAL ftmp;

  //loop the faces
  for(INTEGER i=0;i<N-1;i++)
    {
      INTEGER iA=i-2;
      INTEGER iB=i-1;
      INTEGER iC=i;
      INTEGER iD=i+1;
      INTEGER iE=i+2;
      INTEGER iF=i+3;

      //point A
      if(iA>=0)
        {//in the domain
          rhoA = Mq1[iA];
          uA   = Mq2[iA]/rhoA;
          pA   = G9*(Mq3[iA]-0.5*rhoA*uA*uA);
        }
      else
        {//outside of the domain
          rhoA = Mq1[0];
          uA   = Mq2[0]/rhoA;
          pA   = G9*(Mq3[0]-0.5*rhoA*uA*uA);
        }
       
      //point B
      if(iB>=0)
        {//in the domain
          rhoB = Mq1[iB];
          uB   = Mq2[iB]/rhoB;
          pB   = G9*(Mq3[iB]-0.5*rhoB*uB*uB);
        }
      else
        {//outside of the domain
          rhoB = Mq1[0];
          uB   = Mq2[0]/rhoB;
          pB   = G9*(Mq3[0]-0.5*rhoB*uB*uB);
        }

      //point C

          rhoC = Mq1[iC];
          uC   = Mq2[iC]/rhoC;
          pC   = G9*(Mq3[iC]-0.5*rhoC*uC*uC);


      //point D

          rhoD = Mq1[iD];
          uD   = Mq2[iD]/rhoD;
          pD   = G9*(Mq3[iD]-0.5*rhoD*uD*uD);



      //point E
      if(iE<=N-1)
        {//in the domain
          rhoE = Mq1[iE];
          uE   = Mq2[iE]/rhoE;
          pE   = G9*(Mq3[iE]-0.5*rhoE*uE*uE);
        }
      else
        {//outside of the domain
          rhoE = Mq1[N-1];
          uE   = Mq2[N-1];
          pE   = G9*(Mq3[N-1]-0.5*Mq2[N-1]*Mq2[N-1]/Mq1[N-1]);
        }

      //point F
      if(iF<=N-1)
        {//in the domain
          rhoF = Mq1[iF];
          uF   = Mq2[iF]/rhoF;
          pF   = G9*(Mq3[iF]-0.5*rhoF*uF*uF);
        }
      else
        {//outside of the domain
          rhoF = Mq1[N-1];
          uF   = Mq2[N-1];
          pF   = G9*(Mq3[N-1]-0.5*Mq2[N-1]*Mq2[N-1]/Mq1[N-1]);
        }

      //compute the smooth coefficient, beta0,beta1,beta2,beta3,beta4,beta5
      REAL fArho=rhoA-rhoB;  REAL fAu=uA-uB;  REAL fAp=pA-pB;
      REAL fBrho=rhoB-rhoC;  REAL fBu=uB-uC;  REAL fBp=pB-pC;
      REAL fCrho=rhoC-rhoD;  REAL fCu=uC-uD;  REAL fCp=pC-pD;
      REAL fDrho=rhoD-rhoE;  REAL fDu=uD-uE;  REAL fDp=pD-pE;
      REAL fErho=rhoE-rhoF;  REAL fEu=uE-uF;  REAL fEp=pE-pF;

      ftmp=fArho-fBrho; REAL ftmprhoAB=ftmp*ftmp*13/12.0; ftmp=fArho-3*fBrho; REAL beta0rho=ftmprhoAB+0.25*ftmp*ftmp;
      ftmp=fBrho-fCrho; REAL ftmprhoBC=ftmp*ftmp*13/12.0; ftmp=fBrho+fCrho;   REAL beta1rho=ftmprhoBC+0.25*ftmp*ftmp;
      ftmp=fCrho-fDrho; REAL ftmprhoCD=ftmp*ftmp*13/12.0; ftmp=3*fCrho-fDrho; REAL beta2rho=ftmprhoCD+0.25*ftmp*ftmp;

      ftmp=fAu-fBu; REAL ftmpuAB=ftmp*ftmp*13/12.0; ftmp=fAu-3*fBu; REAL beta0u=ftmpuAB+0.25*ftmp*ftmp;
      ftmp=fBu-fCu; REAL ftmpuBC=ftmp*ftmp*13/12.0; ftmp=fBu+fCu;   REAL beta1u=ftmpuBC+0.25*ftmp*ftmp;
      ftmp=fCu-fDu; REAL ftmpuCD=ftmp*ftmp*13/12.0; ftmp=3*fCu-fDu; REAL beta2u=ftmpuCD+0.25*ftmp*ftmp;

      ftmp=fAp-fBp; REAL ftmppAB=ftmp*ftmp*13/12.0; ftmp=fAp-3*fBp; REAL beta0p=ftmppAB+0.25*ftmp*ftmp;
      ftmp=fBp-fCp; REAL ftmppBC=ftmp*ftmp*13/12.0; ftmp=fBp+fCp;   REAL beta1p=ftmppBC+0.25*ftmp*ftmp;
      ftmp=fCp-fDp; REAL ftmppCD=ftmp*ftmp*13/12.0; ftmp=3*fCp-fDp; REAL beta2p=ftmppCD+0.25*ftmp*ftmp;

      // REAL beta0u=13/12.0*(fAu-fBu)*(fAu-fBu)+0.25*(fAu-3*fBu)*(fAu-3*fBu);
      // REAL beta1u=13/12.0*(fBu-fCu)*(fBu-fCu)+0.25*(fBu+fCu)*(fBu+fCu);
      // REAL beta2u=13/12.0*(fCu-fDu)*(fCu-fDu)+0.25*(3*fCu-fDu)*(3*fCu-fDu);

      // REAL beta0p=13/12.0*(fAp-fBp)*(fAp-fBp)+0.25*(fAp-3*fBp)*(fAp-3*fBp);
      // REAL beta1p=13/12.0*(fBp-fCp)*(fBp-fCp)+0.25*(fBp+fCp)*(fBp+fCp);
      // REAL beta2p=13/12.0*(fCp-fDp)*(fCp-fDp)+0.25*(3*fCp-fDp)*(3*fCp-fDp);

      ftmp=fBrho-3*fCrho; REAL beta3rho=ftmprhoBC+0.25*ftmp*ftmp;
      ftmp=fCrho+fDrho;   REAL beta4rho=ftmprhoCD+0.25*ftmp*ftmp;
      ftmp=fDrho-fErho;   REAL ftmprhoDE=ftmp*ftmp*13/12.0; ftmp=3*fDrho-fErho; REAL beta5rho=ftmprhoDE+0.25*ftmp*ftmp;

      ftmp=fBu-3*fCu; REAL beta3u=ftmpuBC+0.25*ftmp*ftmp;
      ftmp=fCu+fDu;   REAL beta4u=ftmpuCD+0.25*ftmp*ftmp;
      ftmp=fDu-fEu;   REAL ftmpuDE=ftmp*ftmp*13/12.0; ftmp=3*fDu-fEu; REAL beta5u=ftmpuDE+0.25*ftmp*ftmp;

      ftmp=fBp-3*fCp; REAL beta3p=ftmppBC+0.25*ftmp*ftmp;
      ftmp=fCp+fDp;   REAL beta4p=ftmppCD+0.25*ftmp*ftmp;
      ftmp=fDp-fEp;   REAL ftmppDE=ftmp*ftmp*13/12.0; ftmp=3*fDp-fEp; REAL beta5p=ftmppDE+0.25*ftmp*ftmp;

      // REAL beta3rho=13/12.0*(fBrho-fCrho)*(fBrho-fCrho)+0.25*(fBrho-3*fCrho)*(fBrho-3*fCrho);
      // REAL beta4rho=13/12.0*(fCrho-fDrho)*(fCrho-fDrho)+0.25*(fCrho+fDrho)*(fCrho+fDrho);
      // REAL beta5rho=13/12.0*(fDrho-fErho)*(fDrho-fErho)+0.25*(3*fDrho-fErho)*(3*fDrho-fErho);

      // REAL beta3u=13/12.0*(fBu-fCu)*(fBu-fCu)+0.25*(fBu-3*fCu)*(fBu-3*fCu);
      // REAL beta4u=13/12.0*(fCu-fDu)*(fCu-fDu)+0.25*(fCu+fDu)*(fCu+fDu);
      // REAL beta5u=13/12.0*(fDu-fEu)*(fDu-fEu)+0.25*(3*fDu-fEu)*(3*fDu-fEu);

      // REAL beta3p=13/12.0*(fBp-fCp)*(fBp-fCp)+0.25*(fBp-3*fCp)*(fBp-3*fCp);
      // REAL beta4p=13/12.0*(fCp-fDp)*(fCp-fDp)+0.25*(fCp+fDp)*(fCp+fDp);
      // REAL beta5p=13/12.0*(fDp-fEp)*(fDp-fEp)+0.25*(3*fDp-fEp)*(3*fDp-fEp);

      //Compute alpha0,alpha1,alpha2,alph3,alpha4,alpha5
      REAL alpha0rho=0.1/(beta0rho+epsi)/(beta0rho+epsi);
      REAL alpha1rho=0.6/(beta1rho+epsi)/(beta1rho+epsi);
      REAL alpha2rho=0.3/(beta2rho+epsi)/(beta2rho+epsi);

      REAL alpha0u=0.1/(beta0u+epsi)/(beta0u+epsi);
      REAL alpha1u=0.6/(beta1u+epsi)/(beta1u+epsi);
      REAL alpha2u=0.3/(beta2u+epsi)/(beta2u+epsi);

      REAL alpha0p=0.1/(beta0p+epsi)/(beta0p+epsi);
      REAL alpha1p=0.6/(beta1p+epsi)/(beta1p+epsi);
      REAL alpha2p=0.3/(beta2p+epsi)/(beta2p+epsi);

      REAL alpha3rho=0.3/(beta3rho+epsi)/(beta3rho+epsi);
      REAL alpha4rho=0.6/(beta4rho+epsi)/(beta4rho+epsi);
      REAL alpha5rho=0.1/(beta5rho+epsi)/(beta5rho+epsi);

      REAL alpha3u=0.3/(beta3u+epsi)/(beta3u+epsi);
      REAL alpha4u=0.6/(beta4u+epsi)/(beta4u+epsi);
      REAL alpha5u=0.1/(beta5u+epsi)/(beta5u+epsi);

      REAL alpha3p=0.3/(beta3p+epsi)/(beta3p+epsi);
      REAL alpha4p=0.6/(beta4p+epsi)/(beta4p+epsi);
      REAL alpha5p=0.1/(beta5p+epsi)/(beta5p+epsi);

      //compute w0,w1,w2,w3,w4,w5
      ftmp=alpha0rho+alpha1rho+alpha2rho;
      REAL w0rho=alpha0rho/ftmp;
      REAL w1rho=alpha1rho/ftmp;
      REAL w2rho=alpha2rho/ftmp;

      ftmp=alpha0u+alpha1u+alpha2u;
      REAL w0u=alpha0u/ftmp;
      REAL w1u=alpha1u/ftmp;
      REAL w2u=alpha2u/ftmp;

      ftmp=alpha0p+alpha1p+alpha2p;
      REAL w0p=alpha0p/ftmp;
      REAL w1p=alpha1p/ftmp;
      REAL w2p=alpha2p/ftmp;

      ftmp=alpha3rho+alpha4rho+alpha5rho;
      REAL w3rho=alpha3rho/ftmp;
      REAL w4rho=alpha4rho/ftmp;
      REAL w5rho=alpha5rho/ftmp;

      ftmp=alpha3u+alpha4u+alpha5u;
      REAL w3u=alpha3u/ftmp;
      REAL w4u=alpha4u/ftmp;
      REAL w5u=alpha5u/ftmp;

      ftmp=alpha3p+alpha4p+alpha5p;
      REAL w3p=alpha3p/ftmp;
      REAL w4p=alpha4p/ftmp;
      REAL w5p=alpha5p/ftmp;        

      //compute q0,q1,q2,q3,q4,q5
      REAL q0rho=( 2*rhoA-7*rhoB+11*rhoC)/6.0;
      REAL q1rho=(  -rhoB+5*rhoC+ 2*rhoD)/6.0;
      REAL q2rho=( 2*rhoC+5*rhoD-   rhoE)/6.0;

      REAL q0u=( 2*uA-7*uB+11*uC)/6.0;
      REAL q1u=(  -uB+5*uC+ 2*uD)/6.0;
      REAL q2u=( 2*uC+5*uD-   uE)/6.0;

      REAL q0p=( 2*pA-7*pB+11*pC)/6.0;
      REAL q1p=(  -pB+5*pC+ 2*pD)/6.0;
      REAL q2p=( 2*pC+5*pD-   pE)/6.0;

      // REAL q0u= 2.0/6*uA-7.0/6*uB+11.0/6*uC;
      // REAL q1u=-1.0/6*uB+5.0/6*uC+ 2.0/6*uD;
      // REAL q2u= 2.0/6*uC+5.0/6*uD- 1.0/6*uE;

      // REAL q0p= 2.0/6*pA-7.0/6*pB+11.0/6*pC;
      // REAL q1p=-1.0/6*pB+5.0/6*pC+ 2.0/6*pD;
      // REAL q2p= 2.0/6*pC+5.0/6*pD- 1.0/6*pE;

      REAL q3rho= q1rho;
      REAL q4rho= q2rho;
      REAL q5rho= (11*rhoD-7*rhoE+2*rhoF)/6.0;

      REAL q3u= q1u;
      REAL q4u= q2u;
      REAL q5u= (11*uD-7*uE+2*uF)/6.0;

      REAL q3p= q1p;
      REAL q4p= q2p;
      REAL q5p= (11*pD-7*pE+2*pF)/6.0;

      //compute qL,qR
      // qL1=w0rho*q0rho+w1rho*q1rho+w2rho*q2rho;
      // qL2=w0rho*q0u  +w1rho*q1u  +w2rho*q2u;
      // qL3=w0rho*q0p  +w1rho*q1p  +w2rho*q2p;

      // qR1=w3rho*q3rho+w4rho*q4rho+w5rho*q5rho;
      // qR2=w3rho*q3u  +w4rho*q4u  +w5rho*q5u;
      // qR3=w3rho*q3p  +w4rho*q4p  +w5rho*q5p;

      qL1=w0rho*q0rho+w1rho*q1rho+w2rho*q2rho;
      qL2=w0u*q0u+w1u*q1u+w2u*q2u;
      qL3=w0p*q0p+w1p*q1p+w2p*q2p;

      qR1=w3rho*q3rho+w4rho*q4rho+w5rho*q5rho;
      qR2=w3u*q3u+w4u*q4u+w5u*q5u;
      qR3=w3p*q3p+w4p*q4p+w5p*q5p;

      //Riemann solver
      FaceFluxHLLC(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
      //FaceFluxGudunov(qL1,qL2,qL3,qR1,qR2,qR3,MF1[i],MF2[i],MF3[i]);
    }

}

void NumericalSHOCK::FaceFluxGudunov(REAL qL1, REAL qL2, REAL qL3, 
                                     REAL qR1, REAL qR2, REAL qR3,
                                     REAL& MF1, REAL& MF2, REAL& MF3)
{
  REAL q1,q2,q3,q4;             /*!< density,velocity,pressure,internal energy */
  LocalRiemann(qL1,qL2,qL3,qR1,qR2,qR3,q1,q2,q3,q4);
  MF1=q1*q2;
  MF2=q1*q2*q2+q3;
  MF3=(q3*GAMMA/G9+0.5*q1*q2*q2)*q2;
}

//! compute the interface flux F
/*! 
  
  \param qL1 left state density
  \param qL2 left state velocity
  \param qL3 left state pressure
  \param qR1 right state density
  \param qR2 right state velocity
  \param qR3 right state pressure
  \param MF1 flux for density
  \param MF2 flux for moment
  \param MF3 flux for total energy
*/
void NumericalSHOCK::FaceFluxHLLCBACKUP(REAL rhoL, REAL uL, REAL pL, 
                                  REAL rhoR, REAL uR, REAL pR,
                                  REAL& MF1, REAL& MF2, REAL& MF3)
{
  //0->basic variables for left and right state
  REAL cL=sqrt(GAMMA*pL/rhoL);  /*!< sound speed */
  REAL HL=GAMMA/(GAMMA-1)*pL/rhoL+0.5*uL*uL; /*!< enthapy */

  REAL cR=sqrt(GAMMA*pR/rhoR);   /*!< soud speed */
  REAL HR=GAMMA/(GAMMA-1)*pR/rhoR+0.5*uR*uR; /*!< enthapy */

  //1->Roe's average
  REAL f=sqrt(rhoR/rhoL);
  REAL Au=(uL+uR*f)/(1+f);
  REAL AH=(HL+HR*f)/(1+f);
  REAL Ac=sqrt((GAMMA-1)*(AH-0.5*uL*uL));

  //2->Normal velocity
  REAL AU=Au;

  //3->Left state wave speed
  REAL UL=uL;
  REAL SL=min(UL-cL,AU-Ac);

  //4->Right state wave speed
  REAL UR=uR;
  REAL SR=max(UR+cR,AU+Ac);

  //5->Middle state wave speed
  REAL g1=rhoR*UR*(SR-UR)-rhoL*UL*(SL-UL)+(pL-pR);
  REAL g2=rhoR*(SR-UR)-rhoL*(SL-UL);
  REAL SM=g1/g2;

  //6->The interface flux
  if(SL>0){
    //F=FL
    MF1=rhoL*UL;
    MF2=rhoL*uL*UL+pL;
    MF3=rhoL*HL*UL;    
    return;
  }

  if(SL<=0&&0<SM){
    //F=FL_star
    REAL w1L=rhoL;              /*!< left state density */
    REAL w2L=rhoL*uL;           /*!< left state moment rho*u */
    REAL w3L=pL/(GAMMA-1)+0.5*rhoL*uL*uL; /*!< left state energy */

    REAL F1L=rhoL*UL;           /*!< flux for density */
    REAL F2L=rhoL*uL*UL+pL;     /*!< flux for moment */
    REAL F3L=rhoL*HL*UL;        /*!< flux for energy */

    REAL pM=rhoL*(UL-SL)*(UL-SM)+pL; /*!< middle state pressure p* */

    REAL f1=SL-UL;
    REAL f2=SL-SM;
    REAL f3=f1/f2;
    REAL f4=(pM-pL);
    REAL f5=f4/f2;

    REAL w1M=rhoL*f3;           /*!< middle left state density */
    REAL w2M=rhoL*uL*f3+f5;     /*!< middle left state moment rho*u */
    REAL w3M=w3L*f3+(-pL*UL+pM*SM)/f2; /*!< middle left state energy E */
    
    MF1=F1L+SL*(w1M-w1L);
    MF2=F2L+SL*(w2M-w2L);
    MF3=F3L+SL*(w3M-w3L);

    return;
  }

  if(SM<=0&&0<=SR){
    //F=FR_star
    REAL w1R=rhoR;              /*!< right state density */
    REAL w2R=rhoR*uR;           /*!< right state moment rho*u */
    REAL w3R=pR/(GAMMA-1)+0.5*rhoR*uR*uR; /*!< right state energy */

    REAL F1R=rhoR*UR;           /*!< flux for density */
    REAL F2R=rhoR*uR*UR+pR;     /*!< flux for moment */
    REAL F3R=rhoR*HR*UR;        /*!< flux for energy */

    REAL pM=rhoL*(UL-SL)*(UL-SM)+pL; /*!< middle state pressure p* */

    REAL f1=SR-UR;
    REAL f2=SR-SM;
    REAL f3=f1/f2;
    REAL f4=(pM-pR);
    REAL f5=f4/f2;

    REAL w1M=rhoR*f3;           /*!< middle right state density */
    REAL w2M=rhoR*uR*f3+f5;     /*!< middle right state moment rho*u */
    REAL w3M=w3R*f3+(-pR*UR+pM*SM)/f2; /*!< middle right state energy E */
    
    MF1=F1R+SR*(w1M-w1R);
    MF2=F2R+SR*(w2M-w2R);
    MF3=F3R+SR*(w3M-w3R);

    return;
  }

  if(SR<0){
    //F=FR
    MF1=rhoR*UR;
    MF2=rhoR*uR*UR+pR;
    MF3=rhoR*HR*UR;
    return;
  }
}


void NumericalSHOCK::FaceFluxHLLC(REAL rhoL, REAL uL, REAL pL, 
                                  REAL rhoR, REAL uR, REAL pR,
                                  REAL& MF1, REAL& MF2, REAL& MF3)
{
  REAL wL[4],wR[4];
  REAL kas,kbs,Lab;
  REAL ddt, GI[4];
  wL[0]=rhoL;
  wL[1]=rhoL*uL;
  wL[2]=0;
  wL[3]=pL/(GAMMA-1)+0.5*rhoL*uL*uL;

  wR[0]=rhoR;
  wR[1]=rhoR*uR;
  wR[2]=0;
  wR[3]=pR/(GAMMA-1)+0.5*rhoR*uR*uR;

  kas=1.0;
  kbs=0.0;
  Lab=1.0;

  riemann_hllc_2d_(wL,wR,kas,kbs,Lab,ddt,GI);

  MF1=GI[0];
  MF2=GI[1];
  MF3=GI[3];
  
}

void NumericalSHOCK::Advance(REAL T)
{

  //set time step
  dt=0.00001;

  INTEGER ITmax=T/dt;
  //for(INTEGER j=0;j<25000;j++)
  for(INTEGER j=0;j<ITmax;j++)
    {

      //compute the flux
      if(Scheme==0)
        Flux();
      else if(Scheme==1)
        Flux_MUSCL();
      else if(Scheme==2)
        Flux_WENO3();
      else if(Scheme==3)
        Flux_WENO5();

      //update
      for(INTEGER i=1;i<N-1;i++)
        {
          Mq1[i]=Mq1[i]-dt/dx*(MF1[i]-MF1[i-1]);
          Mq2[i]=Mq2[i]-dt/dx*(MF2[i]-MF2[i-1]);
          Mq3[i]=Mq3[i]-dt/dx*(MF3[i]-MF3[i-1]);
          Mq4[i]=(Mq3[i]-0.5*Mq2[i]*Mq2[i]/Mq1[i])/Mq1[i];
        }
      
      cout<<"it="<<j+1<<endl;
    }
}

void NumericalSHOCK::Solution(REAL T)
{
  //set up the mesh
  MeshSet();

  //set up the initial flow field
  Init();

  //time advancing
  Advance(T);
}
