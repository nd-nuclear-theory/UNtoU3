// UNtoU3.h - a header file that implements the algorithm for reduction of an U(N) irreducible representation (irrep) into U(3) irreps.
//
// License: BSD 2-Clause (https://opensource.org/licenses/BSD-2-Clause)
//
// Copyright (c) 2019, Daniel Langr
// All rights reserved.
//
// 03/01/21 (pjf): Modify algorithm to allow U(N) irreps with up to 4 particles per state.

#ifndef UNTOU3_H
#define UNTOU3_H

// Optionally, define before inclsion of this header file:
//
// #define UNTOU3_ENABLE_OPENMP     : enable parallelization of the algorithm based on OpenMP
// #define UNTOU3_DISABLE_TCE       : disable tail call elimination in recursive calls
// #define UNTOU3_DISABLE_UNORDERED : disable the use of a hash table for U(3) weights (binary search tree is used instead)
// #define UNTOU3_DISABLE_PRECALC   : disable precalculation of low Gelfand pattern rows to U(3) weights

#ifndef UNTOU3_DISABLE_PRECALC
#error Pre-calc not yet implemented for U(4).
#endif

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#ifndef UNTOU3_DISABLE_UNORDERED
#include <unordered_map>
#else
#include <map>
#endif

#ifndef UNTOU3_DISABLE_UNORDERED
// An auxiliary struct that implements a hasher for an array of 3 numbers
template <typename T>
struct array_3_hasher
{
   std::size_t operator()(const std::array<T, 3>& key) const
   {
      return key[0] + (key[1] << 8) + (key[2] << 16);
   }
};
#endif /* UNTOU3_DISABLE_UNORDERED */

// Generates U(3) weights and their multiplicities in an input U(N) irrep and allows to evaluate their level dimensionalities.
// Labels of U(N) are limited to {4,3,2,1,0}.
//
// T - type of the computer representation of U(3) weight labels
// U - type used for storing multiplicities of U(3) weights
//
// Example:
//    UNtoU3<> gen;
//    gen.generateXYZ(2);             // n=2, N=(n+1)*(n+2)/2
//    gen.generateU3Weights(0, 0, 1, 4, 1); // n4=0, n3=0, n2=1, n1=4, n0=1 represent number of twos, ones, and zeros in the input U(N) irrep, respectively
//    for (const auto& pair : gen.multMap()) { // iterate over U(3) weights
//       const auto & w = pair.first; // get weight labels
//       if (auto D_l = gen.getLevelDimensionality(w)) // get weight level dimensionality
//          std::cout << "[" << w[0] << "," << w[1] << "," << w[2] << "] : " << D_l << std::endl;
//    }
template <typename T = uint32_t, typename U = uint32_t>
class UNtoU3 {
   public:
      UNtoU3();
      ~UNtoU3();

      // type for storing labels of U(3) weights
      using U3Weight = std::array<T, 3>;

      // type of the data structure used for storing U(3) weights and their multiplicities
#ifndef UNTOU3_DISABLE_UNORDERED
      using U3MultMap = std::unordered_map<U3Weight, U, array_3_hasher<T>>;
#else
      using U3MultMap = std::map<U3Weight, U>;
#endif /* U3MultMap */

      // type of the representation of a single Gelfand pattern row, which contain its number of fours, threes, twos, ones, and zeros
      using GelfandRow = std::array<uint16_t, 5>;
      // type of a Galfand pattern labels
      using GRT = GelfandRow::value_type;

      // definition of the order of axis for weight vectors
      enum { NZ, NX, NY };

      // Generates HO quanta vectors for given nth HO level.
      // Need to be used befor generateU3Weights member function is called.
      void generateXYZ(int n);

      // Generates U(3) weights and their multiplicities for an input U(N) irrep [f].
      // [f] is specified by the number of twos n2, ones n1, and zeros n0.
      // N=n2+n1+n0 must be equal to (n+1)*(n+2)/2, where n was used as an argument of generateXYZ.
      void generateU3Weights(GelfandRow gpr);
      inline void generateU3Weights(uint16_t n4, uint16_t n3, uint16_t n2, uint16_t n1, uint16_t n0) { generateU3Weights({n4,n3,n2,n1,n0}); }

      // Provides an access to the table of U(3) weights and their multiplicities generated by generateU3Weights.
      // Returns a constant reference to the computer representation of this table of type U3MultMap.
      const U3MultMap& multMap() const { return mult_; }

      // Get level dimensionality for a given U(3) weight passed as an argument.
      U getLevelDimensionality(const U3Weight&) const;

   private:
      // HO quanta vectors generated by generateXYZ
      std::array<std::vector<uint32_t>, 3> xyz_;
      // table of resulting U(3) irreps and their multiplicities
      U3MultMap mult_;

#ifndef UNTOU3_DISABLE_PRECALC
      // arrays used to store precalculated contributions of low-level Gelfand patterns into resulting U(3) weights
      std::array<std::array<std::array<std::array<std::array<uint8_t, 4>, 4>, 4>, 4>, 4> cnt_;
      std::array<std::array<std::array<std::array<std::array<uint16_t*, 4>, 4>, 4>, 4>, 4> ptr_;
      std::array<uint16_t, 3 * 45> conts_;
#endif

#ifdef UNTOU3_ENABLE_OPENMP
      // thread-local tables for generated U(3) weights and their multiplicites, which are finally merged into mult_
      static U3MultMap* mult_tl_;
#pragma omp threadprivate(mult_tl_)
#endif

      // Recursive function for generation of Gelfand patterns.
      // It calls itself for all possible lower Gelfand pattern rows generated by the input Gelfand pattern row.
      // At the end of recursion, it increments the multiplicity of resulting U(3) weight ith mult_ (or mult_tl_).
      // gpr - representation of an input gelfand pattern row
      // pp - partial contribution of higher Gelfand pattern rows to the generated U(3) weights
      void generateU3WeightsRec(GelfandRow gpr, U3Weight pp);

#ifndef UNTOU3_DISABLE_PRECALC
      // initialize arrays cnt_ and ptr_
      void init_cnt_ptr();
      // initialize array conts_
      void init_conts();
#endif
};

#ifdef UNTOU3_ENABLE_OPENMP
template <typename T, typename U>
typename UNtoU3<T, U>::U3MultMap* UNtoU3<T, U>::mult_tl_;
#endif

template <typename T, typename U>
UNtoU3<T, U>::UNtoU3()
{
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp parallel
   mult_tl_ = new U3MultMap{};
#endif

#ifndef UNTOU3_DISABLE_PRECALC
   init_cnt_ptr();
#endif
}

template <typename T, typename U>
UNtoU3<T, U>::~UNtoU3()
{
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp parallel
   delete mult_tl_;
#endif
}

template <typename T, typename U>
void UNtoU3<T, U>::generateXYZ(int n)
{
   for (auto & v : xyz_) v.clear();

   for (int k = 0; k <= n; k++) {
      uint32_t nz = n - k;
      for (int nx = k; nx >= 0; nx--) {
         xyz_[NX].push_back(nx);
         xyz_[NY].push_back(n - nz - nx);
         xyz_[NZ].push_back(nz);
      }
   }

#ifndef UNTOU3_DISABLE_PRECALC
   init_conts();
#endif
}

template <typename T, typename U>
U UNtoU3<T, U>::getLevelDimensionality(const U3Weight& labels) const
{
   T f1 = labels[0], f2 = labels[1], f3 = labels[2];
   if ((f1 < f2) || (f2 < f3)) return 0;

   auto u3_mult = mult_.find({f1, f2, f3});
   assert(u3_mult != mult_.end());
   auto mult = u3_mult->second;

   u3_mult = mult_.find({f1 + 1, f2 + 1, f3 - 2});
   mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 2, f2 - 1, f3 - 1});
   mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 2, f2, f3 - 2});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 1, f2 - 1, f3});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1, f2 + 1, f3 - 1});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   return mult;
}

template <typename T, typename U>
void UNtoU3<T, U>::generateU3Weights(GelfandRow gpr)
{
   mult_.clear(); // ???
// for (auto & e : mult_) e.second = 0; // ???

#ifdef UNTOU3_ENABLE_OPENMP

#pragma omp parallel
   {
      mult_tl_->clear(); // ???
   // for (auto & e : *mult_tl_) e.second = 0; // ???
#pragma omp single
      generateU3WeightsRec(gpr, {0, 0, 0});
#pragma omp critical
      for (const auto temp : *mult_tl_)
         mult_[temp.first] += temp.second;
   }

#else  /* UNTOU3_ENABLE_OPENMP */

   generateU3WeightsRec(gpr, {0, 0, 0});

#endif /* UNTOU3_ENABLE_OPENMP */
}

template <typename T, typename U>
void UNtoU3<T, U>::generateU3WeightsRec(GelfandRow gpr, U3Weight pp)
{
   auto& [n4, n3, n2, n1, n0] = gpr;
   size_t N = n4 + n3 + n2 + n1 + n0 - 1;

#ifndef UNTOU3_DISABLE_TCE
   while
#else
   if
#endif
   (N >
#ifndef UNTOU3_DISABLE_PRECALC
   2
#else
   0
#endif
   ) {
       if (n4 > 0) {
           if ((n3 > 0) || (n2 > 0) || (n1 > 0) || (n0 > 0)) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { (GRT)(n4 - 1), n3, n2, n1, n0 },
                       { pp[0] + 4 * xyz_[0][N], pp[1] + 4 * xyz_[1][N], pp[2] + 4 * xyz_[2][N] });
               if (n2 > 0) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), (GRT)(n2 - 1), n1, n0 },
                           { pp[0] + 3 * xyz_[0][N], pp[1] + 3 * xyz_[1][N], pp[2] + 3 * xyz_[2][N] });
                   if (n0 > 0)
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                       generateU3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), (GRT)(n2 - 1), (GRT)(n1 + 1), (GRT)(n0 - 1) },
                               { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
               }
               if (n1 > 0) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), n2, (GRT)(n1 - 1), n0 },
                           { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), n3, (GRT)(n2 + 1), (GRT)(n1 - 1), n0 },
                           { pp[0] + 3 * xyz_[0][N], pp[1] + 3 * xyz_[1][N], pp[2] + 3 * xyz_[2][N] });
               }
               if (n0 > 0) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), (GRT)(n3 + 1), n2, n1, (GRT)(n0 - 1) },
                           { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), n3, (GRT)(n2 + 1), n1, (GRT)(n0 - 1) },
                           { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(n4 - 1), n3, n2, (GRT)(n1 + 1), (GRT)(n0 - 1) },
                           { pp[0] + 3 * xyz_[0][N], pp[1] + 3 * xyz_[1][N], pp[2] + 3 * xyz_[2][N] });
               }
           }
           else {  // if ((n3 == 0) && (n2 == 0) && (n1 == 0) && (n0 == 0))
#ifndef UNTOU3_DISABLE_TCE
               n4--;
               pp[0] += 4 * xyz_[0][N]; pp[1] += 4 * xyz_[1][N]; pp[2] += 4 * xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { (GRT)(n4 - 1), n3, n2, n1, n0 },
                       { pp[0] + 4 * xyz_[0][N], pp[1] + 4 * xyz_[1][N], pp[2] + 4 * xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }
       }

       if (n3 > 0) {
           if ((n2 > 0) || (n1 > 0) || (n0 > 0)) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, (GRT)(n3 - 1), n2, n1, n0 },
                       { pp[0] + 3 * xyz_[0][N], pp[1] + 3 * xyz_[1][N], pp[2] + 3 * xyz_[2][N] });
               if (n1 > 0)
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { n4, (GRT)(n3 - 1), (GRT)(n2 + 1), (GRT)(n1 - 1), n0 },
                           { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
               if (n0 > 0) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { n4, (GRT)(n3 - 1), (GRT)(n2 + 1), n1, (GRT)(n0 - 1) },
                           { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { n4, (GRT)(n3 - 1), n2, (GRT)(n1 + 1), (GRT)(n0 - 1) },
                           { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
               }
           }
           else {  // if ((n2 == 0) && (n1 == 0) && (n0 == 0))
#ifndef UNTOU3_DISABLE_TCE
               n3--;
               pp[0] += 3 * xyz_[0][N]; pp[1] += 3 * xyz_[1][N]; pp[2] += 3 * xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, (GRT)(n3 - 1), n2, n1, n0 },
                       { pp[0] + 3 * xyz_[0][N], pp[1] + 3 * xyz_[1][N], pp[2] + 3 * xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }
       }

       if (n2 > 0) {
           if ((n1 > 0) || (n0 > 0)) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, n3, (GRT)(n2 - 1), n1, n0 },
                       { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
               if (n0 > 0)
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { n4, n3, (GRT)(n2 - 1), (GRT)(n1 + 1), (GRT)(n0 - 1) },
                           { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
           }
           else {  // if ((n1 == 0) && (n0 == 0))
#ifndef UNTOU3_DISABLE_TCE
               n2--;
               pp[0] += 2 * xyz_[0][N]; pp[1] += 2 * xyz_[1][N]; pp[2] += 2 * xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, n3, (GRT)(n2 - 1), n1, n0 },
                       { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }
       }

       if (n1 > 0)
           if (n0 > 0)
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, n3, n2, (GRT)(n1 - 1), n0 },
                       { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
           else {  // if (n0 == 0)
#ifndef UNTOU3_DISABLE_TCE
              n1--;
              pp[0] += xyz_[0][N]; pp[1] += xyz_[1][N]; pp[2] += xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, n3, n2, (GRT)(n1 - 1), n0 },
                       { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }

       if (n0 > 0) {
#ifndef UNTOU3_DISABLE_TCE
           n0--;
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { n4, n3, n2, n1, (GRT)(n0 - 1) }, { pp[0], pp[1], pp[2] });
#endif /* UNTOU3_DISABLE_TCE */
       }

#ifndef UNTOU3_DISABLE_TCE
       N--;
#endif
   }
#ifdef UNTOU3_DISABLE_TCE
   else {
#endif

#ifndef UNTOU3_DISABLE_PRECALC

   auto c = cnt_[n4][n3][n2][n1][n0];
   if (c > 0) {
      auto p = ptr_[n4][n3][n2][n1][n0];
      for (int i = 0; i < c; i++) {
         U3Weight pp_;
         pp_[0] = pp[0] + *p++;
         pp_[1] = pp[1] + *p++;
         pp_[2] = pp[2] + *p++;
#ifdef UNTOU3_ENABLE_OPENMP
         (*mult_tl_)[pp_] += 1;
#else
         mult_[pp_] += 1;
#endif
      }
   }
   else  {
#ifdef UNTOU3_ENABLE_OPENMP
      (*mult_tl_)[pp] += 1;
#else
      mult_[pp] += 1;
#endif
   }

#else /* UNTOU3_DISABLE_PRECALC */

   auto temp = 4 * n4 + 3 * n3 + 2 * n2 + n1;
   pp[0] += temp * xyz_[0][0];
   pp[1] += temp * xyz_[1][0];
   pp[2] += temp * xyz_[2][0];

#ifdef UNTOU3_ENABLE_OPENMP
   (*mult_tl_)[pp] += 1;
#else
   mult_[pp] += 1;
#endif

#endif /* UNTOU3_DISABLE_PRECALC */

#ifdef UNTOU3_DISABLE_TCE
   }
#endif
}

#ifndef UNTOU3_DISABLE_PRECALC

template <typename T, typename U>
void UNtoU3<T, U>::init_cnt_ptr()
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
   ptr_[0][0][3][0][0] = p; p += 1 * 3; ptr_[0][0][0][3][0] = p; p += 1 * 3;
   ptr_[0][0][2][1][0] = p; p += 3 * 3; ptr_[0][0][2][0][1] = p; p += 6 * 3; ptr_[0][0][1][2][0] = p; p += 3 * 3;
   ptr_[0][0][0][2][1] = p; p += 3 * 3; ptr_[0][0][1][0][2] = p; p += 6 * 3; ptr_[0][0][0][1][2] = p; p += 3 * 3;
   ptr_[0][0][1][1][1] = p; p += 8 * 3;
   ptr_[0][0][2][0][0] = p; p += 1 * 3; ptr_[0][0][0][2][0] = p; p += 1 * 3;
   ptr_[0][0][1][1][0] = p; p += 2 * 3; ptr_[0][0][1][0][1] = p; p += 3 * 3; ptr_[0][0][0][1][1] = p; p += 2 * 3;
   ptr_[0][0][1][0][0] = p; p += 1 * 3; ptr_[0][0][0][1][0] = p;
}

template <typename T, typename U>
void UNtoU3<T, U>::init_conts()
{
   auto p = conts_.data();
   // (3,0,0)
   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   // (0,3,0)
   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (2,1,0)
   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   // (2,0,1)
   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,2,0)
   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (0,2,1)
   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,0,2)
   *p++ = 0 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 0 * xyz_[2][0];

   // (0,1,2)
   *p++ = 0 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,1,1)
   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (2,0,0)
   *p++ = 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 2 * xyz_[2][0];

   // (0,2,0)
   *p++ = 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (1,1,0)
   *p++ = 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 1 * xyz_[2][0];

   // (1,0,1)
   *p++ = 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (0,1,1)
   *p++ = 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,0,0)
   *p++ = 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][0];

   // (0,1,0)
   *p++ = 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][0];
}

#endif /* UNTOU3_DISABLE_PRECALC */

#endif /* UNTOU3_H */
