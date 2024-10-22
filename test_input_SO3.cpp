// test_input.cpp - a test driver for UNtoSO3 class.
//
// License: BSD 2-Clause (https://opensource.org/licenses/BSD-2-Clause)
//
// Copyright (c) 2019, Daniel Langr
// All rights reserved.
//
// Program implements the U(2l+1) to SO(3) the input irrep [f] specified by the orbital angular momentum l,
// and its number of fours, threes, twos, ones, and zeros read from the standard input.
// For instance, for the input U(21) irrep [f] = [2,2,2,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
// the user should provide the following numbers: 5 0 0 6 1 14.
//
// The program performs the U(2l+1) to SO(3) reduction and calculates the sum of the dimensions
// of resulting SO(3) irrpes multiplied by their level dimensionalities, and print it to the
// standard output. For instance, for the input irrep specified above, the output should read:
// SO(3) irreps total dim = 2168999910
//
// This sum should be equal to dim[f], which can be calculated analytically with the support
// of rational numbers. The program performs this calculation as well if the Boost library
// is available and uses its Boost.Rational sublibrary. Availability of Boost is indicated by
// users by definition of HAVE_BOOST preprocessor symbol.
// For the input irrep [f] specified above, the program should first print out:
// U(N) irrep dim = 2168999910

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#ifdef HAVE_BOOST
#include <boost/rational.hpp>
#endif

// #define UNTOSO3_DISABLE_TCE
//#define UNTOSO3_DISABLE_UNORDERED
#define UNTOSO3_DISABLE_PRECALC
#define UNTOSO3_ENABLE_OPENMP
#include "UNtoSO3.h"

#ifdef HAVE_BOOST
// Impelments analytical formula for calculation of a dimension of a generic U(N) irrep [f].
// [f] is specified by its labels passed as an array (generally any indexed data structure)
// argument irrep.
template <typename T>
unsigned long dim(const T & irrep) {
   const auto N = irrep.size();
   boost::rational<unsigned long> result{1};
   for (uint32_t l = 2; l <= N; l++)
      for (uint32_t k = 1; k <= l - 1; k++)
         result *= { irrep[k - 1] - irrep[l - 1] + l - k, l - k };

   assert(result.denominator() == 1);
   return result.numerator();
}
#endif

// Implements analytical formula for calculcation of a dimension of an input SO(3) irrep.
// (Does not require rational arithmetics.)
unsigned long dim(const UNtoSO3<>::SO3Weight & l) {
   return 2*l+1;
}

int main() {
   // angular momentum level
   unsigned long l;
   // specification of intput U(N) irrep
   unsigned short n4, n3, n2, n1, n0;
   std::cin >> l >> n4 >> n3 >> n2 >> n1 >> n0;

   if (n4 + n3 + n2 + n1 + n0 != 2*l+1)
      throw std::invalid_argument("Arguments mismatch!");

#ifdef HAVE_BOOST
   // analytical calculation of dim([f])
   std::vector<unsigned long> f(n4, 4);
   std::fill_n(std::back_inserter(f), n3, 3);
   std::fill_n(std::back_inserter(f), n2, 2);
   std::fill_n(std::back_inserter(f), n1, 1);
   std::fill_n(std::back_inserter(f), n0, 0);
   std::cout << "U(N) irrep dim = " << dim(f) << std::endl;
#endif

   UNtoSO3<> gen;
   // generate weight vector for a given l
   gen.generateM(l);
   // generation of SO(3) irreps in the input U(N) irrep [f]
   gen.generateSO3Weights(n4, n3, n2, n1, n0);
   // calculated sum
   unsigned long sum = 0;
   // iteration over generated SO(3) weights
   for (const auto & pair : gen.multMap()) {
      // get SO(3) weight lables
      const auto & weight = pair.first;
      // get its level dimensionality if its nonzero and the SO(3) weight is a SO(3) irrep
      if (auto D_l = gen.getLevelDimensionality(weight)) {
         // add contribution of this SO(3) irrep to the sum
         sum += D_l * dim(weight);
         std::cout << " [" << weight << "]"
                   << " : " << D_l
                   << std::endl;
     }
   }
   std::cout << "SO(3) irreps total dim = " << sum << std::endl;
}
