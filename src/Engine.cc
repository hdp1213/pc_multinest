//--------------------------------------------------------------------------
//
// Description:
// 	class Engine : see header file (Engine.hh) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "Engine.hh"
// ------------------------
// Collaborating classes --
//-------------------------
//--------------------
// C++
//--------------------
#include <numeric>
#include <iostream>
#include <stdexcept>
//--------------------
// C
//----------------

//---------------
// Constructors --
//----------------
Engine::Engine():m_lmax(-1) {
  // constructor
}

Engine::Engine(int lmax):m_lmax(lmax) {
  // constructor
}
//--------------
// Destructor --
//--------------

//-----------------
// Member functions --
//-----------------

void
Engine::write_Cls(std::ostream &out){
  // Create vector of l values
  std::vector<unsigned> lvec(m_lmax-1, 1);
  lvec[0] = 2;
  std::partial_sum(lvec.begin(), lvec.end(), lvec.begin());

  // Get all spectra values
  std::vector<double> cltt, clte, clee, clbb, clpp, cltp, clep;
  bool hasLensing = false;
  try {
    get_Cls(lvec, cltt, clte, clee, clbb);
    hasLensing = get_lensing_Cls(lvec, clpp, cltp, clep);
  }
  catch (std::exception &e) {
    std::cerr << "[ERROR] write_Cls(): " << e.what() << std::endl;
  }

  std::cout.precision(16);
  for (std::size_t i = 0; i < lvec.size(); i++) {
    out << lvec[i] << "\t"
        << cltt[i] << "\t"
        << clte[i] << "\t"
        << clee[i] << "\t"
        << clbb[i];

    if (hasLensing) {
      out << "\t"
          << clpp[i] << "\t"
          << cltp[i] << "\t"
          << clep[i];
    }

    out << std::endl;
  }
}
