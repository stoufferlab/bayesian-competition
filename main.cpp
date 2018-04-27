//      ----------
//      LIBRAIRIES
//      ----------
#include <cstdlib>                              // librairie standard
#include <iostream>                             // librairie flux i/o avec ecran
#include <fstream>                              // librairie flux i/o avec fichiers
#include <cmath>                                // librairie maths
#include <ctime>                                // librairie temps
#include <cstring>                              // librairie chaine de characteres
#include <vector>				// librairie conteneur vector
#include <algorithm>				// librairie algorithmes STL
#include "AnnPlant_model.hpp"

using namespace std;

int main()
{
  // set seed
  srand(time(0));

  // define oject model and the number of species
  class Model model(2);

  // define parameters for the simulation (length, thresholds, maximum attractor period checked)
  model.set_simu(1000,1e-6,1e-8,10);

  // define the growth-competition model function 
  model.set_model(5);
  // read parameter values from input file
  model.set_growthrate("lambda.txt");
  model.set_compet("alpha.txt");

  // run simulations for a quasi Monte-Carlo sample of initial conditions
  model.run_qMCsample_initcond(0,5,10,"./RES/tests/");


  return 0;
}
