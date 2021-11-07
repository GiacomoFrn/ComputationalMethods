#ifndef System_h
#define System_h

#include<vector>
using namespace std;

class System {

 public:

  // constructor
  System(int Num, double kf, double mass, double dT, double Amp, double al=0, double be=0 );
  // deleted copy constructor and assignment to prevent unadvertent copy
  System           ( const System& x ) = delete;
  System& operator=( const System& x ) = delete;

  // destructor
  virtual ~System();

  void Update();
  void Compute();
  void Initialize();
  
  void   computeEnergy(  int n );
  double getEnergy( int n );
  
  double FunctionPr( int i );
  double FunctionNe( int i );

  double getPrVel( int i ) const;
  double getNeVel( int i ) const;
  double getPrPos( int i ) const;
  double getNePos( int i ) const;
  
  double getNperiod( int n );
  double getN() const;

 private:
 
	int         N; // number of oscillators without the two fixed extremals
	double      k; // hooke constant
	double      m; // mass of oscillators
	double deltaT; // delta time
	double      A; // amplitude
	double      a; // quadratic coefficient
	double 		b; // 3^ coefficient
	double 		w; // oscillation frequency

	struct Position{
		double pr; // present position
		double ne; // next position
	};
	
	struct Velocity{
		double pr; // present position
		double ne; // next position
	};
	
	// position and velocity vectors of pointers 
	vector<Position> Npos;
	vector<Velocity> Nvel;
	
	// energy coefficient vectors
	vector<double> nC;
	vector<double> nDC;
	
	// energy vectors;
	vector<double> nE;
};

#endif

