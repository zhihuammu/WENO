#include "Riemann.hh"

using namespace std;


//! Sod problem
/*! 
  
 */
void runCase0()
{
  //exact Riemann solver
  SHOCK SOD(1.0,0.0,1.0,0.125,0.0,0.1);   //case 1 - Toro's shock tube
  SOD.ExactRiemann(-0.5,0.5,0.2);

  //numeircal Riemann solver
  NumericalSHOCK NSOD(1.0,0.0,1.0,0.125,0.0,0.1,0);  //case 1
  NSOD.Solution(0.2);
  NSOD.Out();

  NumericalSHOCK NSOD1(1.0,0.0,1.0,0.125,0.0,0.1,1);  //case 1
  NSOD1.Solution(0.2);
  NSOD1.Out();

  NumericalSHOCK NSOD2(1.0,0.0,1.0,0.125,0.0,0.1,2);  //case 1
  NSOD2.Solution(0.2);
  NSOD2.Out();

  NumericalSHOCK NSOD3(1.0,0.0,1.0,0.125,0.0,0.1,3);  //case 1
  NSOD3.Solution(0.2);
  NSOD3.Out();
}

//! Toro's shock tube
/*! 
  
  A modified Sod problem
 */
void runCase1()
{
  //exact Riemann solver
  SHOCK SOD(1.0,0.75,1.0,0.125,0.0,0.1);   //case 1 - Toro's shock tube
  SOD.ExactRiemann(-0.5,0.5,0.2);

  //numeircal Riemann solver
  NumericalSHOCK NSOD(1.0,0.75,1.0,0.125,0.0,0.1,0);  //case 1
  NSOD.Solution(0.2);
  NSOD.Out();

  NumericalSHOCK NSOD1(1.0,0.75,1.0,0.125,0.0,0.1,1);  //case 1
  NSOD1.Solution(0.2);
  NSOD1.Out();

  NumericalSHOCK NSOD2(1.0,0.75,1.0,0.125,0.0,0.1,2);  //case 1
  NSOD2.Solution(0.2);
  NSOD2.Out();

  NumericalSHOCK NSOD3(1.0,0.75,1.0,0.125,0.0,0.1,3);  //case 1
  NSOD3.Solution(0.2);
  NSOD3.Out();


}

void runCase2()
{
  SHOCK SOD2(1.0,-2.0,0.4,1.0,2.0,0.4);  //case 2
  SOD2.ExactRiemann(-0.5,0.5,0.15);

  NumericalSHOCK NSOD(1.0,-2.0,0.4,1.0,2.0,0.4,0);   //case 2
  NSOD.Solution(0.15);
  NSOD.Out();

  NumericalSHOCK NSOD1(1.0,-2.0,0.4,1.0,2.0,0.4,1);   //case 2
  NSOD1.Solution(0.15);
  NSOD1.Out();

  NumericalSHOCK NSOD2(1.0,-2.0,0.4,1.0,2.0,0.4,2);   //case 2
  NSOD2.Solution(0.15);
  NSOD2.Out();

  NumericalSHOCK NSOD3(1.0,-2.0,0.4,1.0,2.0,0.4,3);   //case 2
  NSOD3.Solution(0.15);
  NSOD3.Out();
}
void runCase3()
{
  SHOCK SOD3(1.0,0.0,1000.0,1.0,0.0,0.01);  //case 3
  SOD3.ExactRiemann(-0.5,0.5,0.012);

  NumericalSHOCK NSOD(1.0,0.0,1000.0,1.0,0.0,0.01,0);   //case 2
  NSOD.Solution(0.012);
  NSOD.Out();

  NumericalSHOCK NSOD1(1.0,0.0,1000.0,1.0,0.0,0.01,1);   //case 2
  NSOD1.Solution(0.012);
  NSOD1.Out();

  NumericalSHOCK NSOD2(1.0,0.0,1000.0,1.0,0.0,0.01,2);   //case 2
  NSOD2.Solution(0.012);
  NSOD2.Out();

  NumericalSHOCK NSOD3(1.0,0.0,1000.0,1.0,0.0,0.01,3);   //case 2
  NSOD3.Solution(0.012);
  NSOD3.Out();
}

void runCase4()
{
  SHOCK SOD4(1.0,0.0,0.01,1.0,0.0,100.0);  //case 4
  SOD4.ExactRiemann(-0.5,0.5,0.035);

  NumericalSHOCK NSOD(1.0,0.0,0.01,1.0,0.0,100.0,0);   //case 2
  NSOD.Solution(0.035);
  NSOD.Out();

  NumericalSHOCK NSOD1(1.0,0.0,0.01,1.0,0.0,100.0,1);   //case 2
  NSOD1.Solution(0.035);
  NSOD1.Out();

  NumericalSHOCK NSOD2(1.0,0.0,0.01,1.0,0.0,100.0,2);   //case 2
  NSOD2.Solution(0.035);
  NSOD2.Out();

  NumericalSHOCK NSOD3(1.0,0.0,0.01,1.0,0.0,100.0,3);   //case 2
  NSOD3.Solution(0.035);
  NSOD3.Out();
}

void runCase5()
{
  SHOCK SOD5(5.99924,19.5975,460.894,5.99242,-6.19633,46.0950);  //case 5
  SOD5.ExactRiemann(-0.5,0.5,0.035);

  NumericalSHOCK NSOD(5.99924,19.5975,460.894,5.99242,-6.19633,46.0950,0);   //case 2
  NSOD.Solution(0.035);
  NSOD.Out();

  NumericalSHOCK NSOD1(5.99924,19.5975,460.894,5.99242,-6.19633,46.0950,1);   //case 2
  NSOD1.Solution(0.035);
  NSOD1.Out();

  NumericalSHOCK NSOD2(5.99924,19.5975,460.894,5.99242,-6.19633,46.0950,2);   //case 2
  NSOD2.Solution(0.035);
  NSOD2.Out();

  NumericalSHOCK NSOD3(5.99924,19.5975,460.894,5.99242,-6.19633,46.0950,3);   //case 2
  NSOD3.Solution(0.035);
  NSOD3.Out();
}
void runCase6()
{
  // SHOCK SOD6(0.445,0.698876404,3.52773,0.5,0.0,0.571);
  // SOD6.ExactRiemann(-0.5,0.5,0.15);

  // NumericalSHOCK NSOD(0.445,0.698876404,3.52773,0.5,0.0,0.571,0);   //case 2
  // NSOD.Solution(0.15);
  // NSOD.Out();

  // NumericalSHOCK NSOD1(0.445,0.698876404,3.52773,0.5,0.0,0.571,1);   //case 2
  // NSOD1.Solution(0.15);
  // NSOD1.Out();

  // NumericalSHOCK NSOD2(0.445,0.698876404,3.52773,0.5,0.0,0.571,2);   //case 2
  // NSOD2.Solution(0.15);
  // NSOD2.Out();

  // NumericalSHOCK NSOD3(0.445,0.698876404,3.52773,0.5,0.0,0.571,3);   //case 2
  // NSOD3.Solution(0.15);
  // NSOD3.Out();


  SHOCK SOD6(0.5,0.0,0.571,0.445,-0.698876404,3.52773);
  SOD6.ExactRiemann(-0.5,0.5,0.15);

  NumericalSHOCK NSOD(0.5,0.0,0.571,0.445,-0.698876404,3.52773,0);   //case 2
  NSOD.Solution(0.15);
  NSOD.Out();

  NumericalSHOCK NSOD1(0.5,0.0,0.571,0.445,-0.698876404,3.52773,1);   //case 2
  NSOD1.Solution(0.15);
  NSOD1.Out();

  NumericalSHOCK NSOD2(0.5,0.0,0.571,0.445,-0.698876404,3.52773,2);   //case 2
  NSOD2.Solution(0.15);
  NSOD2.Out();

  NumericalSHOCK NSOD3(0.5,0.0,0.571,0.445,-0.698876404,3.52773,3);   //case 2
  NSOD3.Solution(0.15);
  NSOD3.Out();


}


int main()
{

  cout<<"Hello Riemann"<<endl;

  
  runCase6();
  
  return 0;
}
