#ifndef IsingLattice_h
#define IsingLattice_h

#include <vector>
#include <fstream>

class IsingLattice {

 public:

  // constructor
  IsingLattice( int nx, int ny, double j);
  // deleted copy constructor and assignment to prevent unadvertent copy
  IsingLattice           ( const IsingLattice& x ) = delete;
  IsingLattice& operator=( const IsingLattice& x ) = delete;
  
  // destructor
  virtual ~IsingLattice();
  
  // shuffle spin configuration (total or singol spin)
  void singleShuffle( double r );
  void totalShuffle();
  
  // print spinMat
  void printSpinMat();
  
  // write spinMat
  void writeSpinMat(std::ofstream &file);
  
  // compute energy 
  double getEn() const;
  
  // compute magnetization
  double getMagn() const;
  
  // get energy difference
  double getDiffEn() const;
  
 private:
  
  double Nx;
  double Ny;
  double J;
  
  std::vector<int> spinMat;
  
  double diffEn;
  
};

#endif

