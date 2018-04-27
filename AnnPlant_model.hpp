//      ----------
//      LIBRAIRIES
//      ----------
#include <cstdlib>                              // librairie standard
#include <iostream>                             // librairie flux i/o avec ecran
#include <fstream>                              // librairie flux i/o avec fichiers
#include <cmath>                                // librairie maths
#include <ctime>                                // librairie temps
#include <cstring>                              // librairie chaine de characteres
#include <vector>                               // librairie conteneur vector
#include <algorithm>

using namespace std;


#ifndef _DEF_MODEL_
#define _DEF_MODEL_

//      -----------------------------------------------------------------------------------
//	CLASS THAT DESCRIBES AND SIMULATES THE MODEL OF ANNUAL PLANT GROWTH AND COMPETITION
//      -----------------------------------------------------------------------------------
class Model
{
//-----------------
// member functions
//-----------------
  public:

  // constructor
  Model(int const n):m_n(n),m_nSquare(n*n)
  {
    m_lambda.resize(m_n);
    m_alpha.resize(m_nSquare);
  }
  // destructor
  ~Model()
  {
    for (int i=0;i<=m_Tmax;++i) delete m_Nsimu[i];
  }

  // set model
  void set_model(int const mod)
  {
    m_model = mod;
    switch(mod)
    {
      case 1: cout << "Model: Holling-like function" << endl; break;
      case 3: cout << "Model: exponential function" << endl; break;
      case 4: cout << "Model: Holling-like with coefficients as exponents" << endl; break;
      case 5: cout << "Model: linear function" << endl; break;
    }
  }

  // set simulation
  void set_simu(int const Tmax, double const thrExtinct, double const thrCV, int const max_attract_period)
  {
    m_Nsimu.resize(Tmax+1);
    for (int i=0;i<=Tmax;++i)
    {
      m_Nsimu[i] = new vector<double>;
      m_Nsimu[i]->resize(m_n);
    }
    m_Tmax = Tmax;
    cout << "Simulation length: " << Tmax << endl;
    m_thrExtinct = thrExtinct;
    cout << "Extinction threshold: " << thrExtinct << endl;
    m_thrCV = thrCV;
    cout << "Convergence threshold: " << thrCV << endl;
    m_max_attract_period = max_attract_period;
    cout << "Maximum attractor period checked: " << max_attract_period << endl;
  }

  // set initial condition (same for all)
  void set_init_allthesame(int const N0)
  {
    for (int i=0;i<m_n;++i) m_Nsimu[0]->at(i) = N0;
    cout << "All species biomass initialized to: " << N0 << endl;
  }

  // set initial condition (random between two numbers)
  void set_init_random(double const Nmin, double const Nmax)
  {
    double const range(Nmax-Nmin);
    for (int i=0;i<m_n;++i) m_Nsimu[0]->at(i) = Nmin + (rand()/(double)RAND_MAX) * range;
    cout << "Species biomass initialized to random values between: " << Nmin << " and " << Nmax << endl;
  }

  // read growth rate in a file
  void set_growthrate(const char* growth)
  {
    ifstream myfile (growth);
    if(myfile.is_open())
    {
      int i(0);
      cout << "Annual growth rates:" << endl;
      while ( i < m_n && myfile >> m_lambda[i] )
      {
	cout << m_lambda[i] << endl;
	++i;
      }
      myfile.close();
    }
    else cerr << "Unable to open file" << endl;
  }

  // read competition coefficients in a file
  void set_compet(const char* compet)
  {
    ifstream myfile (compet);
    if(myfile.is_open())
    {
      int i(0),j(0);
      cout << "Competition coefficients:" << endl;
      while ( i < m_nSquare && myfile >> m_alpha[i] )
      {
	cout << m_alpha[i] << " ";
	++i;++j;
	if (j==m_n)
	{
	  cout << endl;
	  j = 0;
	}
      }
      myfile.close();
    }
    else cerr << "Unable to open file" << endl;
  }

  // compute N(t+1) for all species with model 1
  void timeincrement_model1(int const t)
  {
    double temp;
    int k(0),tM1(t-1);
    for (int i=0;i<m_n;++i)
    {
      temp = 1;
      for (int j=0;j<m_n;++j)
      {
	temp += m_Nsimu[tM1]->at(j) * m_alpha[k];
	++k;
      }
      m_Nsimu[t]->at(i) = m_Nsimu[tM1]->at(i) * m_lambda[i] / temp;
    }
  }

  // compute N(t+1) for all species with model 3
  void timeincrement_model3(int const t)
  {
    double temp;
    int k(0),tM1(t-1);
    for (int i=0;i<m_n;++i)
    {
      temp = 0;
      for (int j=0;j<m_n;++j)
      {
	temp -= m_Nsimu[tM1]->at(j) * m_alpha[k];
	++k;
      }
      m_Nsimu[t]->at(i) = m_Nsimu[tM1]->at(i) * m_lambda[i] * exp(temp);
    }
  }

  // compute N(t+1) for all species with model 4
  void timeincrement_model4(int const t)
  {
    double temp;
    int k(0),tM1(t-1);
    for (int i=0;i<m_n;++i)
    {
      temp = 1;
      for (int j=0;j<m_n;++j)
      {
	temp += pow(m_Nsimu[tM1]->at(j),m_alpha[k]);
	++k;
      }
      m_Nsimu[t]->at(i) = m_Nsimu[tM1]->at(i) * m_lambda[i] / temp;
    }
  }

  // compute N(t+1) for all species with model 5
  void timeincrement_model5(int const t)
  {
    double temp;
    int k(0),tM1(t-1);
    for (int i=0;i<m_n;++i)
    {
      temp = m_lambda[i];
      for (int j=0;j<m_n;++j)
      {
	temp -= m_Nsimu[tM1]->at(j) * m_alpha[k];
	++k;
      }
      m_Nsimu[t]->at(i) = m_Nsimu[tM1]->at(i) * temp;
    }
  }

  // run simulation until Tmax
  void run_simu()
  {
    m_attract_period = -1;			// no attractor reached

    cout << endl << "Run simulation:" << endl << "Time";
    for (int i=0;i<m_n;++i) cout << "\tSp. " << i;
    cout << endl;
    cout << 0;
    for (int i=0;i<m_n;++i) cout << "\t" << m_Nsimu[0]->at(i);
    cout << endl;

    // optimized: using a pointer to member function a each time step would be time consuming, so as testing the condition "m_model" at each time step
    switch (m_model)
    {
      case 1:
	for (int t=1;t<=m_Tmax;++t)
	{
	  timeincrement_model1(t);
	  check_extinct_simu(t);	// check extinctions
	  print_simu(t);
	  if (after_timeincrement(t))
	  {
	    m_Tend = t;
	    break;
	  }
	}
	break;
      case 3:
	for (int t=1;t<=m_Tmax;++t)
	{
	  timeincrement_model3(t);
	  check_extinct_simu(t);	// check extinctions
	  print_simu(t);
	  if (after_timeincrement(t))
	  {
	    m_Tend = t;
	    break;
	  }
	}
	break;
      case 4:
	for (int t=1;t<=m_Tmax;++t)
	{
	  timeincrement_model4(t);
	  check_extinct_simu(t);	// check extinctions
	  print_simu(t);
	  if (after_timeincrement(t))
	  {
	    m_Tend = t;
	    break;
	  }
	}
	break;
      case 5:
	for (int t=1;t<=m_Tmax;++t)
	{
	  timeincrement_model5(t);
	  check_extinct_simu(t);	// check extinctions
	  print_simu(t);
	  if (after_timeincrement(t))
	  {
	    m_Tend = t;
	    break;
	  }
	}
	break;
    }

    switch (m_attract_period)
    {
      case -1: cout << "No attractor reached" << endl; break;
      case 1: cout << "Equilibrium" << endl; break;
      default: cout << "Attractor of period " << m_attract_period << endl; break;
    }
  }

  // check species biomass lower than extinction threshold
  void check_extinct_simu(int const t){for (int i=0;i<m_n;++i) if (m_Nsimu[t]->at(i) < m_thrExtinct) m_Nsimu[t]->at(i) = 0;}

  // check equilibrium is reached
  bool reach_attractor(const int t, const int p)
  {
    double temp(0);
    int const tMp(t-p);
    for (int i=0;i<m_n;++i)
    {
      if (m_Nsimu[t]->at(i)>0) temp += fabs(m_Nsimu[t]->at(i)-m_Nsimu[tMp]->at(i)) / m_Nsimu[tMp]->at(i);
    }
    temp /= m_n;
    return (temp<m_thrCV);
  }

  // after each time step: check extinctions and if attractor is reached
  bool after_timeincrement(int const t)
  {
    for (int p=1;p<m_max_attract_period;++p)
    {
      if (t>=p)
      {
	if (reach_attractor(t,p))
	{
	  m_attract_period = p;	// save attractor period
	  return true;		// signal to stop simulation
	}
      }
    }
    return false;	// continue simulation
  }

  // print a time step of simulation
  void const print_simu(int const t)
  {
    cout << t;
    for (int i=0;i<m_n;++i) cout << "\t" << m_Nsimu[t]->at(i);
    cout << endl;
  }

  // save simulation results
  void const write_simu(const char* simu)
  {
    int Tini;
    ofstream myfile (simu);
    if(myfile.is_open())
    {
      if (m_attract_period == -1)
      {
	Tini = 0;
	m_Tend = m_Tmax;
      }
      else Tini = m_Tend - m_attract_period;
      for (int t=Tini;t<=m_Tend;++t)
      {
	for (int i=0;i<m_n;++i) myfile << m_Nsimu[t]->at(i) << " ";
	myfile << endl;
      }
      myfile.close();
    }
    else cerr << "Unable to open file" << endl;
  }

  // big function that run simulations for a quasi Monte-Carlo sampling of the initial condition space (min and max the same for each species)
  void run_qMCsample_initcond(double const min_all, double const max_all, int const n_dim, char const * Rep)
  {
    double const deltaIC((max_all-min_all)/(n_dim));	// size of the subspace for each random draw (quasi Monte-Carlo sampling)
    double const init_max(min_all+deltaIC);
    vector<double> min_sp,max_sp;
    min_sp.assign(n_dim,min_all);
    max_sp.resize(n_dim,min_all+deltaIC);
    int const n_simu(pow(n_dim,m_n));
    char numfich[10],namefich[50];

    for (int i=0;i<n_simu;++i)
    {
      // print
/*      for (int k=0;k<m_n;++k) cout << min_sp[k] << " ";
      cout << endl;
      for (int k=0;k<m_n;++k) cout << max_sp[k] << " ";
      cout << endl << endl;
*/

      // draw initial condition and run simulation
      for (int k=0;k<m_n;++k) m_Nsimu[0]->at(k) = min_sp[k] + (rand()/(double)RAND_MAX) * deltaIC;
      run_simu();

      // write in file
      sprintf(numfich,"%d",i);
      strcpy(namefich,Rep);
      strcat(namefich,"IC_");
      strcat(namefich,numfich);
      strcat(namefich,".txt");
      write_simu(namefich);

      // change subspace for next draw
      for (int k=0;k<m_n;++k)
      {
	if (max_sp[k]<max_all)
	{
	  min_sp[k] += deltaIC;
	  max_sp[k] += deltaIC;
	  break;
	}
	else
	{
	  min_sp[k] = min_all;
	  max_sp[k] = init_max;
	}
      }
    }
  }


//-----------------
// member variables
//-----------------
  private:

  // model definition
  int const m_n;			// number of species (n)
  int const m_nSquare;			// square number of species (number of intra- and inter-competition coefficients and interactions)
  vector<double> m_lambda;		// maximum annual growth rate (size: n)
  vector<double> m_alpha;		// intra and inter-specific competition coefficients (size: n^2)
  double m_b;				// b exponent in model 2
  int m_model;				// code to select model (1: Holling-like function, 3: exponential function, 4: Holling-like with coefficients as exponents, 5: linear function)

  // simulation and analysis
  int m_Tmax;				// maximum simulation time
  int m_Tend;				// last simulation time (if stopped earlier)
  double m_thrExtinct;			// threshold to define extinction
  double m_thrCV;			// threshold to define convergence on an attractor
  int m_max_attract_period;		// maximum period of attractor checked
  int m_attract_period;			// period of the attractor reached (0 if equilibrium, -1 if no attractor reached)
  vector<vector<double>*> m_Nsimu;	// plant population biomass (size: Tmax, each element is a point to a vector of size n)
};



#endif
