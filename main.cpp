#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>

//we include our libraries

using namespace std;
double rand_uniform ( long rand_int );
double randU ();
int iseed = time(0);
double rand_exp ( double mfp );
double rand_gauss ( double u, double v );
double estar_xenon(double Ei_xenon);
double estar_argon(double Ei_argon);
double estar_water(double Ei_water);

//we include our functions

double Energy_datapoint[49]={1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 12.5, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 900.0, 1000.0};
// KeV

double LET_datapoint_xenon[49]={28540, 28630, 27810, 26700, 25540, 23360, 21470, 19860, 18490, 17310, 16290, 15400, 14610, 13280, 12200, 11310, 10550, 9097, 8041, 7238, 6604, 5664, 4997, 4497, 4108, 3796, 3539, 3324, 3142, 2849, 2623, 2443, 2298, 2030, 1848, 1717, 1619, 1483, 1396, 1336, 1294, 1264, 1243, 1227, 1216, 1203, 1199, 1200, 1204 };
// KeV cm^2/g

double LET_datapoint_argon[49]={64650, 58080, 52750, 48370, 44720, 39000, 34710, 31360, 28670, 26450, 24590, 23010, 21640, 19390, 17620, 16170, 14980, 12710, 11110, 9916, 8983, 7619, 6666, 5959, 5414, 4978, 4623, 4327, 4076, 3675, 3368, 3125, 2927, 2567, 2323, 2147, 2015, 1832, 1713, 1631, 1572, 1529, 1497, 1472, 1453, 1428, 1414, 1406, 1404};
// KeV cm^2/g

double LET_datapoint_water[49]={11980, 10370, 91770, 82570, 75220, 64170, 56210, 50180, 45430, 41590, 38410, 35730, 33440, 29720, 26820, 24480,  22560, 18980, 16470, 14610, 13180, 11100, 9657, 8596, 7781, 7134, 6607, 6170, 5801, 5211, 4761, 4407, 4119, 3596, 3242, 2988, 2798, 2533, 2360, 2241, 2154, 2090, 2041, 2003, 1972, 1926, 1896, 1876, 1862};
// KeV cm^2/g

//we include our data

double xenon_dx=1;
double argon_dx=1;
double water_dx=1; 

//and initialize our tracklengths


int main ( int argc, char** argv ) {

  std::cout << std::fixed;
  std::cout<< std::setprecision(8);
  double initialEnergy;
  double Ei_xenon, Ef_xenon, dE_xenon, LET_xenon, stoppingPower_xenon;
  double Ei_argon, Ef_argon, dE_argon, LET_argon, stoppingPower_argon;
  double Ei_water, Ef_water, dE_water, LET_water, stoppingPower_water;
  double xenon_density, argon_density, water_density;
  
  //we initialize our variables

  ofstream xenon_tracklengths;
  ofstream argon_tracklengths;
  ofstream water_tracklengths;
  std::remove("xenon_tracklengths.txt");
  std::remove("argon_tracklengths.txt");
  std::remove("water_tracklengths.txt");
  
  //remove old data files

  double xenon_trackLength = 0.;
  double argon_trackLength = 0.;
  double water_trackLength = 0.;
  
  long numElectrons; 
  cout << "#electrons: ";
  cin >> numElectrons;

  cout << "Xenon Density [g/cm^3] = ";
  cin >> xenon_density;
  
  cout << "Argon Density [g/cm^3] = ";
  cin >> argon_density;
  
  cout << "Water Density [g/cm^3] = ";
  cin >> water_density;
  
  cout << "Energy [keV] = ";
  cin >> initialEnergy;

  //gather our information about the system from the user


  if (initialEnergy>=1 && initialEnergy<2) {
     xenon_dx =xenon_dx * 0.0000001; //cm
     argon_dx =argon_dx * 0.0000001; //cm
     water_dx =water_dx * 0.0000001; //cm
  }
  else if (initialEnergy>=2 && initialEnergy<3) {
     xenon_dx =xenon_dx * 0.0000003; //cm
     argon_dx =argon_dx * 0.0000003; //cm
     water_dx =water_dx * 0.0000003; //cm
  }
  else if (initialEnergy>=3 && initialEnergy<5) {
     xenon_dx =xenon_dx * 0.0000005; //cm
     argon_dx =argon_dx * 0.0000005; //cm
     water_dx =water_dx * 0.0000005; //cm
  }
  else if (initialEnergy>=5 && initialEnergy<15) {
     xenon_dx =xenon_dx * 0.000002; //cm
     argon_dx =argon_dx * 0.000002; //cm
     water_dx =water_dx * 0.000002; //cm
  }
  else if (initialEnergy>=15 && initialEnergy<25) {
     xenon_dx =xenon_dx * 0.000005; //cm
     argon_dx =argon_dx * 0.000005; //cm
     water_dx =water_dx * 0.000005; //cm
  }
  else if (initialEnergy>=25 && initialEnergy<50) {
     xenon_dx =xenon_dx * 0.00001; //cm
     argon_dx =argon_dx * 0.00001; //cm
     water_dx =water_dx * 0.00001; //cm
  }
  else if (initialEnergy>=50 && initialEnergy<100) {
     xenon_dx =xenon_dx * 0.00002; //cm
     argon_dx =argon_dx * 0.00002; //cm
     water_dx =water_dx * 0.00002; //cm
  }
  else if (initialEnergy>=100 && initialEnergy<250) {
     xenon_dx =xenon_dx * 0.00005; //cm
     argon_dx =argon_dx * 0.00005; //cm
     water_dx =water_dx * 0.00005; //cm
  }
  else if (initialEnergy>=250 && initialEnergy<=1000) {
     xenon_dx =xenon_dx * 0.0001; //cm
     argon_dx =argon_dx * 0.0001; //cm
     water_dx =water_dx * 0.0001; //cm
  }
  
  //find the steplength from initial energy
  
  for ( long i = 0; i < numElectrons; i++ ) {
    
    Ei_xenon = initialEnergy; xenon_trackLength = 0.; Ef_xenon = Ei_xenon; 
    
    while ( Ef_xenon > 0 ) {
      LET_xenon = pow ( 10., estar_xenon ( Ei_xenon ) );
      stoppingPower_xenon = LET_xenon*xenon_density * (1. + 0.15 * rand_gauss(rand_uniform(rand()),rand_uniform(rand())));
      Ef_xenon = Ei_xenon - stoppingPower_xenon*xenon_dx; 
      xenon_trackLength += xenon_dx;
      Ei_xenon = Ef_xenon; 
    }
    xenon_tracklengths.open("xenon_tracklengths.txt",ios::app);
    xenon_tracklengths << xenon_trackLength <<endl;
    xenon_tracklengths.close();
   
    
  }
  
  //find our tracklengths in xenon for the given amount of electrons and record them

  for ( long n = 0; n < numElectrons; n++ ) {
    
    Ei_argon = initialEnergy; argon_trackLength = 0.; Ef_argon = Ei_argon; 
    
    while ( Ef_argon > 0 ) {
      LET_argon = pow ( 10., estar_argon ( Ei_argon ) );
      stoppingPower_argon = LET_argon*argon_density * (1. + 0.15 * rand_gauss(rand_uniform(rand()),rand_uniform(rand())));
      Ef_argon = Ei_argon - stoppingPower_argon*argon_dx; 
      argon_trackLength += argon_dx;
      Ei_argon = Ef_argon; 
    }
    argon_tracklengths.open("argon_tracklengths.txt" , ios::app);
    argon_tracklengths << argon_trackLength << endl;
    argon_tracklengths.close();
    
    
  }
  //and again in argon
  
   for ( long m = 0; m < numElectrons; m++ ) {
    
    Ei_water = initialEnergy; water_trackLength = 0.; Ef_water = Ei_water; 
    
    while ( Ef_water > 0 ) {
      LET_water = pow ( 10., estar_water ( Ei_water ) );
      stoppingPower_water = LET_water*water_density * (1. + 0.15 * rand_gauss(rand_uniform(rand()),rand_uniform(rand())));
      Ef_water = Ei_water - stoppingPower_water*water_dx; 
      water_trackLength += water_dx;
      Ei_water = Ef_water; 
    }
    
    water_tracklengths.open("water_tracklengths.txt" , ios::app);
    water_tracklengths << water_trackLength << endl;
    water_tracklengths.close();
    
  }
  
  //and again in water

  return 0;

}


double rand_uniform ( long rand_int ) {
  return double(rand_int) / (double)RAND_MAX;
}
double rand_exp ( double mfp ) {
  return -mfp*log(rand_uniform(rand()));
}
double rand_gauss ( double u, double v ) {
  return sqrt(-2.*log(u))*cos(2*M_PI*v);
}

//we express different probability distributions

double estar_xenon ( double Ei_xenon ) {

 
  for (int k=0; k<49; k++) {
    if ( log10(Energy_datapoint[k]) <= log10(Ei_xenon) && log10(Ei_xenon) < log10(Energy_datapoint[k+1] ) ){
      return log10(LET_datapoint_xenon[k])+(log10(LET_datapoint_xenon[k+1])-log10(LET_datapoint_xenon[k]))*((log10(Ei_xenon)-log10(Energy_datapoint[k]))/(log10(Energy_datapoint[k])-log10(Energy_datapoint[k+1])));
    }
    else if (Ei_xenon=1000){
      return log10(1204);
    }
  }
}

double estar_argon ( double Ei_argon ) {

 
  for (int j=0; j<49; j++) {
    if ( log10(Energy_datapoint[j]) <= log10(Ei_argon) && log10(Ei_argon) < log10(Energy_datapoint[j+1] ) ){
      return log10(LET_datapoint_argon[j])+(log10(LET_datapoint_argon[j+1])-log10(LET_datapoint_argon[j]))*((log10(Ei_argon)-log10(Energy_datapoint[j]))/(log10(Energy_datapoint[j])-log10(Energy_datapoint[j+1])));
    }
    else if (Ei_argon=1000){
      return log10(1204);
    }
  }
}

double estar_water ( double Ei_water ) {

 
  for (int l=0; l<49; l++) {
    if ( log10(Energy_datapoint[l]) <= log10(Ei_water) && log10(Ei_water) < log10(Energy_datapoint[l+1] ) ){
      return log10(LET_datapoint_water[l])+(log10(LET_datapoint_water[l+1])-log10(LET_datapoint_water[l]))*((log10(Ei_water)-log10(Energy_datapoint[l]))/(log10(Energy_datapoint[l])-log10(Energy_datapoint[l+1])));
    }
    else if (Ei_water=1000){
      return log10(1204);
    }
  }
}

//We lograithmically look at the stopping powers of the different media

