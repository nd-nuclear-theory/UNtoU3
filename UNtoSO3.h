// UNtoSO3.h - a header file that implements the algorithm for reduction of an U(N) irreducible representation (irrep) into SO(3) irreps.
//
// License: BSD 2-Clause (https://opensource.org/licenses/BSD-2-Clause)
//
// Copyright (c) 2019, Daniel Langr
// All rights reserved.
//
// 03/01/21 (pjf): Modify algorithm to allow U(N) irreps with up to 4 particles per state.

#ifndef UNTOSO3_H
#define UNTOSO3_H

// Optionally, define before inclsion of this header file:
//
// #define UNTOSO3_ENABLE_OPENMP     : enable parallelization of the algorithm based on OpenMP
// #define UNTOSO3_DISABLE_TCE       : disable tail call elimination in recursive calls
// #define UNTOSO3_DISABLE_PRECALC   : disable precalculation of low Gelfand pattern rows to SO(3) weights

#ifndef UNTOSO3_DISABLE_PRECALC
#error Pre-calc not yet implemented for U(4).
#endif

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#include <unordered_map>

namespace UNtoU3 {

// Generates SO(3) weights and their multiplicities in an input U(N) irrep and allows to evaluate their level dimensionalities.
// Labels of U(N) are limited to {4,3,2,1,0}.
//
// T - type of the computer representation of SO(3) weight labels
// U - type used for storing multiplicities of SO(3) weights
//
// Example:
//    UNtoSO3<> gen;
//    gen.generateXYZ(2);             // n=2, N=(n+1)*(n+2)/2
//    gen.generateSO3Weights(0, 0, 1, 4, 1); // n4=0, n3=0, n2=1, n1=4, n0=1 represent number of twos, ones, and zeros in the input U(N) irrep, respectively
//    for (const auto& pair : gen.multMap()) { // iterate over SO(3) weights
//       const auto & w = pair.first; // get weight labels
//       if (auto D_l = gen.getLevelDimensionality(w)) // get weight level dimensionality
//          std::cout << "[" << w[0] << "," << w[1] << "," << w[2] << "] : " << D_l << std::endl;
//    }
template <typename T = int32_t, typename U = uint32_t>
class UNtoSO3 {
   public:
      UNtoSO3();
      ~UNtoSO3();

      // type for storing labels of SO(3) weights
      using SO3Weight = T;

      using SO3MultMap = std::unordered_map<SO3Weight, U>;

      // type of the representation of a single Gelfand pattern row, which contain its number of fours, threes, twos, ones, and zeros
      using GelfandRow = std::array<uint16_t, 5>;
      // type of a Galfand pattern labels
      using GRT = GelfandRow::value_type;

      // Generates HO quanta vectors for given nth HO level.
      // Need to be used befor generateSO3Weights member function is called.
      void generateM(int l);

      // Generates SO(3) weights and their multiplicities for an input U(N) irrep [f].
      // [f] is specified by the number of twos n2, ones n1, and zeros n0.
      // N=n2+n1+n0 must be equal to (n+1)*(n+2)/2, where n was used as an argument of generateXYZ.
      void generateSO3Weights(GelfandRow gpr);
      inline void generateSO3Weights(uint16_t n4, uint16_t n3, uint16_t n2, uint16_t n1, uint16_t n0) { generateSO3Weights({n4,n3,n2,n1,n0}); }

      // Provides an access to the table of SO(3) weights and their multiplicities generated by generateSO3Weights.
      // Returns a constant reference to the computer representation of this table of type SO3MultMap.
      const SO3MultMap& multMap() const { return mult_; }

      // Get level dimensionality for a given SO(3) weight passed as an argument.
      U getLevelDimensionality(const SO3Weight&) const;

   private:
      // weight vectors generated by generateM
      std::vector<int8_t> m_;
      int8_t l_;
      // table of resulting SO(3) irreps and their multiplicities
      SO3MultMap mult_;

#ifndef UNTOSO3_DISABLE_PRECALC
      // arrays used to store precalculated contributions of low-level Gelfand patterns into resulting SO(3) weights
      std::array<std::array<std::array<std::array<std::array<uint8_t, 4>, 4>, 4>, 4>, 4> cnt_;
      std::array<std::array<std::array<std::array<std::array<int16_t*, 4>, 4>, 4>, 4>, 4> ptr_;
      std::array<int16_t, 45> conts_;
#endif

#ifdef UNTOSO3_ENABLE_OPENMP
      // thread-local tables for generated SO(3) weights and their multiplicites, which are finally merged into mult_
      static SO3MultMap* mult_tl_;
#pragma omp threadprivate(mult_tl_)
#endif

      // Recursive function for generation of Gelfand patterns.
      // It calls itself for all possible lower Gelfand pattern rows generated by the input Gelfand pattern row.
      // At the end of recursion, it increments the multiplicity of resulting SO(3) weight ith mult_ (or mult_tl_).
      // gpr - representation of an input gelfand pattern row
      // pp - partial contribution of higher Gelfand pattern rows to the generated SO(3) weights
      void generateSO3WeightsRec(GelfandRow gpr, SO3Weight pp);

#ifndef UNTOSO3_DISABLE_PRECALC
      // initialize arrays cnt_ and ptr_
      void init_cnt_ptr();
      // initialize array conts_
      void init_conts();
#endif
};

#ifdef UNTOSO3_ENABLE_OPENMP
template <typename T, typename U>
typename UNtoSO3<T, U>::SO3MultMap* UNtoSO3<T, U>::mult_tl_;
#endif

template <typename T, typename U>
UNtoSO3<T, U>::UNtoSO3()
{
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp parallel
   mult_tl_ = new SO3MultMap{};
#endif

#ifndef UNTOSO3_DISABLE_PRECALC
   init_cnt_ptr();
#endif
}

template <typename T, typename U>
UNtoSO3<T, U>::~UNtoSO3()
{
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp parallel
   delete mult_tl_;
#endif
}

template <typename T, typename U>
void UNtoSO3<T, U>::generateM(int l)
{
   m_.clear();

   l_ = l;
   for (int m = -l; m <= l; m++) {
      m_.push_back(m);
   }

#ifndef UNTOSO3_DISABLE_PRECALC
   init_conts();
#endif
}

template <typename T, typename U>
U UNtoSO3<T, U>::getLevelDimensionality(const SO3Weight& l) const
{
   if (l < 0) return 0;
   auto so3_mult = mult_.find(l);
   assert(so3_mult != mult_.end());
   auto mult = so3_mult->second;

   so3_mult = mult_.find(l + 1);
   mult -= (so3_mult == mult_.end()) ? 0 : so3_mult->second;

   return mult;
}

template <typename T, typename U>
void UNtoSO3<T, U>::generateSO3Weights(GelfandRow gpr)
{
   mult_.clear(); // ???
// for (auto & e : mult_) e.second = 0; // ???

#ifdef UNTOSO3_ENABLE_OPENMP

#pragma omp parallel
   {
      mult_tl_->clear(); // ???
   // for (auto & e : *mult_tl_) e.second = 0; // ???
#pragma omp single
      generateSO3WeightsRec(gpr, 0);
#pragma omp critical
      for (const auto temp : *mult_tl_)
         mult_[temp.first] += temp.second;
   }

#else  /* UNTOSO3_ENABLE_OPENMP */

   generateSO3WeightsRec(gpr, 0);

#endif /* UNTOSO3_ENABLE_OPENMP */
}

template <typename T, typename U>
void UNtoSO3<T, U>::generateSO3WeightsRec(GelfandRow gpr, SO3Weight pp)
{
   auto& [n4, n3, n2, n1, n0] = gpr;
   std::size_t N = n4 + n3 + n2 + n1 + n0 - 1;

#ifndef UNTOSO3_DISABLE_TCE
   while
#else
   if
#endif
   (N >
#ifndef UNTOSO3_DISABLE_PRECALC
   2
#else
   0
#endif
   ) {
       if (n4 > 0) {
           if ((n3 > 0) || (n2 > 0) || (n1 > 0) || (n0 > 0)) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { (GRT)(n4 - 1), n3, n2, n1, n0 }, pp + 4 * m_[N]);
               if (n2 > 0) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), (GRT)(n2 - 1), n1, n0 }, pp + 3 * m_[N]);
                   if (n0 > 0)
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                       generateSO3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), (GRT)(n2 - 1), (GRT)(n1 + 1), (GRT)(n0 - 1) }, pp + 2 * m_[N]);
               }
               if (n1 > 0) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), n2, (GRT)(n1 - 1), n0 }, pp + 2 * m_[N]);
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), n3, (GRT)(n2 + 1), (GRT)(n1 - 1), n0 }, pp + 3 * m_[N]);
               }
               if (n0 > 0) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), n2, n1, (GRT)(n0 - 1) }, pp + m_[N]);
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), n3, (GRT)(n2 + 1), n1, (GRT)(n0 - 1) }, pp + 2 * m_[N]);
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { (GRT)(n4 - 1), n3, n2, (GRT)(n1 + 1), (GRT)(n0 - 1) }, pp + 3 * m_[N]);
               }
           }
           else {  // if ((n3 == 0) && (n2 == 0) && (n1 == 0) && (n0 == 0))
#ifndef UNTOSO3_DISABLE_TCE
               n4--;
               pp += 4 * m_[N];
#else /* UNTOSO3_DISABLE_TCE */
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { (GRT)(n4 - 1), n3, n2, n1, n0 }, pp + 4 * m_[N]);
#endif /* UNTOSO3_DISABLE_TCE */
           }
       }

       if (n3 > 0) {
           if ((n2 > 0) || (n1 > 0) || (n0 > 0)) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, (GRT)(n3 - 1), n2, n1, n0 }, pp + 3 * m_[N]);
               if (n1 > 0)
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { n4, (GRT)(n3 - 1), (GRT)(n2 + 1), (GRT)(n1 - 1), n0 }, pp + 2 * m_[N]);
               if (n0 > 0) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { n4, (GRT)(n3 - 1), (GRT)(n2 + 1), n1, (GRT)(n0 - 1) }, pp + m_[N]);
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { n4, (GRT)(n3 - 1), n2, (GRT)(n1 + 1), (GRT)(n0 - 1) }, pp + 2 * m_[N]);
               }
           }
           else {  // if ((n2 == 0) && (n1 == 0) && (n0 == 0))
#ifndef UNTOSO3_DISABLE_TCE
               n3--;
               pp += 3 * m_[N];
#else /* UNTOSO3_DISABLE_TCE */
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, (GRT)(n3 - 1), n2, n1, n0 }, pp + 3 * m_[N]);
#endif /* UNTOSO3_DISABLE_TCE */
           }
       }

       if (n2 > 0) {
           if ((n1 > 0) || (n0 > 0)) {
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, n3, (GRT)(n2 - 1), n1, n0 }, pp + 2 * m_[N]);
               if (n0 > 0)
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateSO3WeightsRec( { n4, n3, (GRT)(n2 - 1), (GRT)(n1 + 1), (GRT)(n0 - 1) }, pp + m_[N]);
           }
           else {  // if ((n1 == 0) && (n0 == 0))
#ifndef UNTOSO3_DISABLE_TCE
               n2--;
               pp += 2 * m_[N];
#else /* UNTOSO3_DISABLE_TCE */
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, n3, (GRT)(n2 - 1), n1, n0 }, pp + 2 * m_[N]);
#endif /* UNTOSO3_DISABLE_TCE */
           }
       }

       if (n1 > 0)
           if (n0 > 0)
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, n3, n2, (GRT)(n1 - 1), n0 }, pp + m_[N]);
           else {  // if (n0 == 0)
#ifndef UNTOSO3_DISABLE_TCE
              n1--;
              pp += m_[N];
#else /* UNTOSO3_DISABLE_TCE */
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, n3, n2, (GRT)(n1 - 1), n0 }, pp + m_[N]);
#endif /* UNTOSO3_DISABLE_TCE */
           }

       if (n0 > 0) {
#ifndef UNTOSO3_DISABLE_TCE
           n0--;
#else /* UNTOSO3_DISABLE_TCE */
#ifdef UNTOSO3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateSO3WeightsRec( { n4, n3, n2, n1, (GRT)(n0 - 1) }, pp);
#endif /* UNTOSO3_DISABLE_TCE */
       }

#ifndef UNTOSO3_DISABLE_TCE
       N--;
#endif
   }
#ifdef UNTOSO3_DISABLE_TCE
   else {
#endif

#ifndef UNTOSO3_DISABLE_PRECALC

   auto c = cnt_[n4][n3][n2][n1][n0];
   if (c > 0) {
      auto p = ptr_[n4][n3][n2][n1][n0];
      for (int i = 0; i < c; i++) {
         SO3Weight pp_;
         pp_ = pp + *p++;
#ifdef UNTOSO3_ENABLE_OPENMP
         (*mult_tl_)[pp_] += 1;
#else
         mult_[pp_] += 1;
#endif
      }
   }
   else  {
#ifdef UNTOSO3_ENABLE_OPENMP
      (*mult_tl_)[pp] += 1;
#else
      mult_[pp] += 1;
#endif
   }

#else /* UNTOSO3_DISABLE_PRECALC */

   auto temp = 4 * n4 + 3 * n3 + 2 * n2 + n1;
   pp += temp * m_[0];

#ifdef UNTOSO3_ENABLE_OPENMP
   (*mult_tl_)[pp] += 1;
#else
   mult_[pp] += 1;
#endif

#endif /* UNTOSO3_DISABLE_PRECALC */

#ifdef UNTOSO3_DISABLE_TCE
   }
#endif
}

#ifndef UNTOSO3_DISABLE_PRECALC

template <typename T, typename U>
void UNtoSO3<T, U>::init_cnt_ptr()
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
               for (int m = 0; m < 4; m++)
                  cnt_[i][j][k][l][m] = 0;

   cnt_[0][0][3][0][0] = 1; cnt_[0][0][0][3][0] = 1;                     //  2
   cnt_[0][0][2][1][0] = 3; cnt_[0][0][2][0][1] = 6; cnt_[0][0][1][2][0] = 3;  // 12
   cnt_[0][0][0][2][1] = 3; cnt_[0][0][1][0][2] = 6; cnt_[0][0][0][1][2] = 3;  // 12
   cnt_[0][0][1][1][1] = 8;                                        //  8
   cnt_[0][0][2][0][0] = 1; cnt_[0][0][0][2][0] = 1;                     //  2
   cnt_[0][0][1][1][0] = 2; cnt_[0][0][1][0][1] = 3; cnt_[0][0][0][1][1] = 2;  //  7
   cnt_[0][0][1][0][0] = 1; cnt_[0][0][0][1][0] = 1;                     //  2

   auto p = conts_.data();
   ptr_[0][0][3][0][0] = p; p += 1; ptr_[0][0][0][3][0] = p; p += 1;
   ptr_[0][0][2][1][0] = p; p += 3; ptr_[0][0][2][0][1] = p; p += 6; ptr_[0][0][1][2][0] = p; p += 3;
   ptr_[0][0][0][2][1] = p; p += 3; ptr_[0][0][1][0][2] = p; p += 6; ptr_[0][0][0][1][2] = p; p += 3;
   ptr_[0][0][1][1][1] = p; p += 8;
   ptr_[0][0][2][0][0] = p; p += 1; ptr_[0][0][0][2][0] = p; p += 1;
   ptr_[0][0][1][1][0] = p; p += 2; ptr_[0][0][1][0][1] = p; p += 3; ptr_[0][0][0][1][1] = p; p += 2;
   ptr_[0][0][1][0][0] = p; p += 1; ptr_[0][0][0][1][0] = p;
}

template <typename T, typename U>
void UNtoSO3<T, U>::init_conts()
{
   auto p = conts_.data();
   // (3,0,0)
   *p++ = 2 * m_[2] + 2 * m_[1] + 2 * m_[0];

   // (0,3,0)
   *p++ = 1 * m_[2] + 1 * m_[1] + 1 * m_[0];

   // (2,1,0)
   *p++ = 1 * m_[2] + 2 * m_[1] + 2 * m_[0];

   *p++ = 2 * m_[2] + 1 * m_[1] + 2 * m_[0];

   *p++ = 2 * m_[2] + 2 * m_[1] + 1 * m_[0];

   // (2,0,1)
   *p++ = 0 * m_[2] + 2 * m_[1] + 2 * m_[0];

   *p++ = 1 * m_[2] + 1 * m_[1] + 2 * m_[0];

   *p++ = 1 * m_[2] + 2 * m_[1] + 1 * m_[0];

   *p++ = 2 * m_[2] + 0 * m_[1] + 2 * m_[0];

   *p++ = 2 * m_[2] + 1 * m_[1] + 1 * m_[0];

   *p++ = 2 * m_[2] + 2 * m_[1] + 0 * m_[0];

   // (1,2,0)
   *p++ = 1 * m_[2] + 1 * m_[1] + 2 * m_[0];

   *p++ = 1 * m_[2] + 2 * m_[1] + 1 * m_[0];

   *p++ = 2 * m_[2] + 1 * m_[1] + 1 * m_[0];

   // (0,2,1)
   *p++ = 0 * m_[2] + 1 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[2] + 0 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[2] + 1 * m_[1] + 0 * m_[0];

   // (1,0,2)
   *p++ = 0 * m_[2] + 0 * m_[1] + 2 * m_[0];

   *p++ = 0 * m_[2] + 1 * m_[1] + 1 * m_[0];

   *p++ = 0 * m_[2] + 2 * m_[1] + 0 * m_[0];

   *p++ = 1 * m_[2] + 0 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[2] + 1 * m_[1] + 0 * m_[0];

   *p++ = 2 * m_[2] + 0 * m_[1] + 0 * m_[0];

   // (0,1,2)
   *p++ = 0 * m_[2] + 0 * m_[1] + 1 * m_[0];

   *p++ = 0 * m_[2] + 1 * m_[1] + 0 * m_[0];

   *p++ = 1 * m_[2] + 0 * m_[1] + 0 * m_[0];

   // (1,1,1)
   *p++ = 0 * m_[2] + 1 * m_[1] + 2 * m_[0];

   *p++ = 0 * m_[2] + 2 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[2] + 0 * m_[1] + 2 * m_[0];

   *p++ = 1 * m_[2] + 1 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[2] + 2 * m_[1] + 0 * m_[0];

   *p++ = 1 * m_[2] + 1 * m_[1] + 1 * m_[0];

   *p++ = 2 * m_[2] + 0 * m_[1] + 1 * m_[0];

   *p++ = 2 * m_[2] + 1 * m_[1] + 0 * m_[0];

   // (2,0,0)
   *p++ = 2 * m_[1] + 2 * m_[0];

   // (0,2,0)
   *p++ = 1 * m_[1] + 1 * m_[0];

   // (1,1,0)
   *p++ = 1 * m_[1] + 2 * m_[0];

   *p++ = 2 * m_[1] + 1 * m_[0];

   // (1,0,1)
   *p++ = 0 * m_[1] + 2 * m_[0];

   *p++ = 2 * m_[1] + 0 * m_[0];

   *p++ = 1 * m_[1] + 1 * m_[0];

   // (0,1,1)
   *p++ = 0 * m_[1] + 1 * m_[0];

   *p++ = 1 * m_[1] + 0 * m_[0];

   // (1,0,0)
   *p++ = 2 * m_[0];

   // (0,1,0)
   *p++ = 1 * m_[0];
}

#endif /* UNTOSO3_DISABLE_PRECALC */

}  // namespace UNtoU3

#ifdef UNTOU3_GLOBAL_NAMESPACE
using namespace UNtoU3;
#endif /* UNTOU3_GLOBAL_NAMESPACE*/

#endif /* UNTOSO3_H */
