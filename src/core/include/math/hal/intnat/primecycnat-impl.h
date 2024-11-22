#ifndef PRIME_CYC_IMPL
#define PRIME_CYC_IMPL

#include "math/math-hal.h"
#include "math/nbtheory.h"

#include "utils/debug.h"
#include "utils/exception.h"
#include "utils/utilities.h"

#include <map>
#include <vector>
#include <iostream>

#include <immintrin.h>

#include "math/hal/intnat/primecycnat.h"

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_bitReverseTableBase2n3;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_base2n3RootPreconTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_base2n3RootTableByModulusRoot;

template <typename VecType>
std::map<usint, std::array<usint, 4>> primecyc::RaderFFTNat<VecType>::m_Base2n3Info;

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_forwardPermutation;

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_inversePermutation;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_rootTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_rootPreconTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_inverseRootTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, std::vector<typename VecType::Integer>> primecyc::RaderFFTNat<VecType>::m_inverseRootPreconTableByModulusRoot;

template <typename VecType>
std::map<usint, bool> primecyc::RaderFFTNat<VecType>::m_enabled;

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::PreComputeIsomorphism(usint cycloOrder) {
    usint g = primecycutil::findPrimitiveRoot(cycloOrder);

    std::vector<usint> forward(cycloOrder);
    std::vector<usint> inverse(cycloOrder);

    uint64_t t = 1;

    for (usint i = 0; i < cycloOrder - 1; i++) {
        if (t == cycloOrder - 1) {
            forward[i] = 0;
            inverse[0] = i;
        } else {
            forward[i] = t;
        }
        inverse[t] = i;
        t = t * g % cycloOrder;
    }
    forward[cycloOrder - 1] = 1;

    m_forwardPermutation[cycloOrder] = forward;
    m_inversePermutation[cycloOrder] = inverse;
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::PreComputeRootTable(usint cycloOrder, const ModulusRoot<IntType>& nttModulusRoot) {
    usint tot = cycloOrder - 1;
    auto rootOfUnityOrder = nttModulusRoot.second.ModExp(tot, nttModulusRoot.first);
    auto rootOfUnityTot = nttModulusRoot.second.ModExp(cycloOrder, nttModulusRoot.first);
    auto &indices = m_inversePermutation[cycloOrder];

    IntType z = rootOfUnityOrder.ModMulFast(IntType(tot).ModInverse(nttModulusRoot.first), nttModulusRoot.first);

    std::vector<typename VecType::Integer> rootTable(tot);
    std::vector<typename VecType::Integer> rootPreconTable(tot);
    std::vector<typename VecType::Integer> inverseRootTable(tot);
    std::vector<typename VecType::Integer> inverseRootPreconTable(tot);

    for (usint i = 1; i <= tot; i++) {
        if (indices[i] == 0) {
            rootTable[0] = z;
        } else {
            rootTable[tot - indices[i]] = z;
        }
        z.ModMulFastEq(rootOfUnityOrder, nttModulusRoot.first);
    }

    std::vector<typename VecType::Integer> rootTableT(tot);
    ForwardFFTBase2n3(rootTable, nttModulusRoot.first, rootOfUnityTot, rootTableT);

    IntType w = IntType(tot).ModInverse(nttModulusRoot.first).ModExp(2, nttModulusRoot.first);

    for (usint i = 0; i < tot; i++) {
        inverseRootTable[i] = rootTableT[i].ModInverse(nttModulusRoot.first).ModMulFast(w, nttModulusRoot.first);
        rootPreconTable[i] = rootTableT[i].PrepModMulConst(nttModulusRoot.first);
        inverseRootPreconTable[i] = inverseRootTable[i].PrepModMulConst(nttModulusRoot.first);
    }

    m_rootTableByModulusRoot[nttModulusRoot] = rootTableT;
    m_inverseRootTableByModulusRoot[nttModulusRoot] = inverseRootTable;
    m_rootPreconTableByModulusRoot[nttModulusRoot] = rootPreconTable;
    m_inverseRootPreconTableByModulusRoot[nttModulusRoot] = inverseRootPreconTable;
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::PreComputeBitReverseTableBase2n3(usint order) {
    usint u = 0, U = 1, v = 0, V = order;
    for (; V % 2 == 0; V >>= 1) {
        u++;
        U <<= 1;
    }
    for (usint t = 1; t < V; t *= 3) {
        v++;
    }

    m_Base2n3Info[order] = {u, U, v, V};

    std::vector<usint> bitReverseTable(order);
    for (usint i = 0; i < V; i++) {
        for (usint j = 0; j < U; j++) {
            bitReverseTable[i * U + j] = primecycutil::reverseBits(j, u, 2) * V + primecycutil::reverseBits(i, v, 3);
        }
    }

    m_bitReverseTableBase2n3[order] = bitReverseTable;
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::PreComputeBase2n3RootTable(usint order, const ModulusRoot<IntType>& nttModulusRoot) {
    auto rootOfUnity = nttModulusRoot.second;
    typename VecType::Integer z = 1;

    std::vector<typename VecType::Integer> rootTable(order);
    std::vector<typename VecType::Integer> rootPreconTable(order);

    for (usint i = 0; i < order; i++) {
        rootTable[i] = z;
        rootPreconTable[i] = z.PrepModMulConst(nttModulusRoot.first);
        z.ModMulFastEq(rootOfUnity, nttModulusRoot.first);
    }

    m_base2n3RootTableByModulusRoot[nttModulusRoot] = rootTable;
    m_base2n3RootPreconTableByModulusRoot[nttModulusRoot] = rootPreconTable;
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::ForwardFFTBase2n3(const VecType& element, const IntType& rootOfUnity, VecType* result) {
    using IntType = primecyc::RaderFFTNat<VecType>::IntType;

    usint n = element.GetLength();

    auto modulus = element.GetModulus();
    result->SetModulus(modulus);

    if (m_bitReverseTableBase2n3.find(n) == m_bitReverseTableBase2n3.end()) {
        PreComputeBitReverseTableBase2n3(n);
    }

    if (m_base2n3RootTableByModulusRoot.find({modulus, rootOfUnity}) == m_base2n3RootTableByModulusRoot.end()) {
        PreComputeBase2n3RootTable(n, {modulus, rootOfUnity});
    }

    const auto &indices = m_bitReverseTableBase2n3[n];
    for (usint i = 0; i < n; i++) {
        (*result)[i] = element[indices[i]];
    }

    auto [u, U, v, V] = m_Base2n3Info[n];

    const auto &rootTable = m_base2n3RootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &rootTablePrecon = m_base2n3RootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    usint l0 = 1, l1 = 1, d = n;

    for (usint i = 0; i < u; i++) {
        l0 = 1 << i;
        l1 = 1 << (i + 1);
        d = n >> (i + 1);
        for (usint j = 0; j != n; j += l1) {
            usint ind_jk = j;
            usint ind_jkl0 = j + l0;
            for (usint k = 0; k < l0; k++) {
                IntType y1 = (*result)[ind_jkl0].ModMulFastConst(rootTable[k * d], modulus, rootTablePrecon[k * d]);
                (*result)[ind_jkl0] = (*result)[ind_jk].ModSubFast(y1, modulus);
                (*result)[ind_jk].ModAddFastEq(y1, modulus);
                ind_jk++;
                ind_jkl0++;
            }
        }
    }

    IntType z3 = rootTable[n / 3], z32 = rootTable[2 * n / 3];
    IntType z3precon = rootTablePrecon[n / 3], z32precon = rootTablePrecon[2 * n / 3];

    for (usint i = 0; i < v; i++) {
        l0 = U;
        l1 = U * 3;
        d = n / (3 * U);
        for (usint t = 0; t < i; t++) {
            l0 *= 3;
            l1 *= 3;
            d /= 3;
        }
        for (usint j = 0; j != n; j += l1) {
            usint ind_kd = 0;
            usint ind_jk = j;
            usint ind_jkl0 = j + l0;
            usint ind_jk2l0 = j + 2 * l0;
            for (usint k = 0; k < l0; k++) {
                IntType y1 = (*result)[ind_jkl0].ModMulFastConst(rootTable[ind_kd], modulus, rootTablePrecon[ind_kd]);
                IntType y2 = (*result)[ind_jk2l0].ModMulFastConst(rootTable[ind_kd * 2], modulus, rootTablePrecon[ind_kd * 2]);
                IntType y0 = y1.ModAddFast(y2, modulus);
                IntType w = y1.ModMulFastConst(z3, modulus, z3precon).ModAddFast(y2.ModMulFastConst(z32, modulus, z32precon), modulus);
                (*result)[ind_jkl0] = (*result)[ind_jk].ModAddFast(w, modulus);
                (*result)[ind_jk2l0] = (*result)[ind_jk].ModSubFast(y0.ModAddFast(w, modulus), modulus);
                (*result)[ind_jk].ModAddFastEq(y0, modulus);
                ind_kd += d;
                ind_jk++;
                ind_jkl0++;
                ind_jk2l0++;
            }
        }
    }
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::ForwardFFTBase2n3(const std::vector<IntType> &element, const IntType &modulus, const IntType &rootOfUnity, std::vector<IntType> &result) {
    using IntType = primecyc::RaderFFTNat<VecType>::IntType;

    usint n = element.size();

    if (m_bitReverseTableBase2n3.find(n) == m_bitReverseTableBase2n3.end()) {
        PreComputeBitReverseTableBase2n3(n);
    }

    if (m_base2n3RootTableByModulusRoot.find({modulus, rootOfUnity}) == m_base2n3RootTableByModulusRoot.end()) {
        PreComputeBase2n3RootTable(n, {modulus, rootOfUnity});
    }

    const auto &indices = m_bitReverseTableBase2n3[n];
    for (usint i = 0; i < n; i++) {
        result[i] = element[indices[i]];
    }

    auto [u, U, v, V] = m_Base2n3Info[n];

    const auto &rootTable = m_base2n3RootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &rootTablePrecon = m_base2n3RootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    usint l0 = 1, l1 = 1, d = n;

    for (usint i = 0; i < u; i++) {
        l0 = 1 << i;
        l1 = 1 << (i + 1);
        d = n >> (i + 1);
        for (usint j = 0; j != n; j += l1) {
            usint ind_jk = j;
            for (usint k = 0; k < l0; k++) {
                IntType y1 = result[ind_jk + l0].ModMulFastConst(rootTable[k * d], modulus, rootTablePrecon[k * d]);
                result[ind_jk + l0] = result[ind_jk].ModSubFast(y1, modulus);
                result[ind_jk].ModAddFastEq(y1, modulus);
                ind_jk++;
            }
        }
    }

    IntType z3 = rootTable[n / 3], z32 = rootTable[2 * n / 3];
    IntType z3precon = rootTablePrecon[n / 3], z32precon = rootTablePrecon[2 * n / 3];

    for (usint i = 0; i < v; i++) {
        l0 = U;
        l1 = U * 3;
        d = n / (3 * U);
        for (usint t = 0; t < i; t++) {
            l0 *= 3;
            l1 *= 3;
            d /= 3;
        }
        for (usint j = 0; j != n; j += l1) {
            usint ind_kd = 0;
            usint ind_jk = j;
            usint ind_jkl0 = j + l0;
            usint ind_jk2l0 = j + 2 * l0;
            for (usint k = 0; k < l0; k++) {
                IntType y1 = result[ind_jkl0].ModMulFastConst(rootTable[ind_kd], modulus, rootTablePrecon[ind_kd]);
                IntType y2 = result[ind_jk2l0].ModMulFastConst(rootTable[ind_kd * 2], modulus, rootTablePrecon[ind_kd * 2]);
                IntType y0 = y1.ModAddFast(y2, modulus);
                IntType w = y1.ModMulFastConst(z3, modulus, z3precon).ModAddFast(y2.ModMulFastConst(z32, modulus, z32precon), modulus);
                result[ind_jkl0] = result[ind_jk].ModAddFast(w, modulus);
                result[ind_jk2l0] = result[ind_jk].ModSubFast(y0.ModAddFast(w, modulus), modulus);
                result[ind_jk].ModAddFastEq(y0, modulus);
                ind_kd += d;
                ind_jk++;
                ind_jkl0++;
                ind_jk2l0++;
            }
        }
    }
}

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m, uint64_t b_inv) {
    uint64_t q = (uint64_t)(((unsigned __int128)a * b_inv) >> 64);
    uint64_t y = a * b - q * m;
    return y >= m ? y - m : y;
}

template <typename VecType>
void primecyc::RaderFFTNat<VecType>::ForwardFFTBase2n3AVX(const std::vector<uint64_t> &element, uint64_t modulus, uint64_t rootOfUnity, std::vector<uint64_t> &result) {

    usint n = element.size();
    uint64_t Q = modulus;

    if (m_bitReverseTableBase2n3.find(n) == m_bitReverseTableBase2n3.end()) {
        PreComputeBitReverseTableBase2n3(n);
    }

    if (m_base2n3RootTableByModulusRoot.find({modulus, rootOfUnity}) == m_base2n3RootTableByModulusRoot.end()) {
        PreComputeBase2n3RootTable(n, {modulus, rootOfUnity});
    }

    const auto &indices = m_bitReverseTableBase2n3[n];
    for (usint i = 0; i < n; i++) {
        result[i] = element[indices[i]];
    }

    auto [u, U, v, V] = m_Base2n3Info[n];

    const auto &rootTable = m_base2n3RootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &rootTablePrecon = m_base2n3RootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    static std::vector<uint64_t> rootTableAVX;
    rootTableAVX.resize(n);
    for (usint i = 0; i < n; i++) {
        rootTableAVX[i] = rootTable[i].ConvertToInt();
    }
    static std::vector<uint64_t> rootTablePreconAVX;
    rootTablePreconAVX.resize(n);
    for (usint i = 0; i < n; i++) {
        rootTablePreconAVX[i] = rootTablePrecon[i].ConvertToInt();
    }

    usint l0 = 1, l1 = 1, d = n;

    __m128i Q_vec = _mm_set1_epi64x(Q);
    __m128i Q_minus_one_vec = _mm_set1_epi64x(Q - 1);

    __m256i Q_vec_256 = _mm256_set1_epi64x(Q);
    __m256i Q_minus_one_vec_256 = _mm256_set1_epi64x(Q - 1);

    // i = 0
    {
        size_t j = 0;

        // Process pairs of elements using AVX2
        for (j = 0; j + 3 < n; j += 4) {
            // Load elements
            __m128i Vj_lo = _mm_loadu_si128((__m128i*)&result[j]);     // [result[j], result[j+1]]
            __m128i Vj_hi = _mm_loadu_si128((__m128i*)&result[j + 2]); // [result[j+2], result[j+3]]

            // Unpack to get V0 and V1
            __m128i V0 = _mm_unpacklo_epi64(Vj_lo, Vj_hi); 
            __m128i V1 = _mm_unpackhi_epi64(Vj_lo, Vj_hi);

            // Compute sums and differences
            __m128i Vsum = _mm_add_epi64(V0, V1);
            __m128i Vdiff = _mm_sub_epi64(V0, V1);

            // Compute mask for sum >= Q
            __m128i mask_sum = _mm_cmpgt_epi64(Vsum, Q_minus_one_vec);

            // Adjust sum where sum >= Q: adj_sum = sum - Q where mask is true
            __m128i adj_sum = _mm_sub_epi64(Vsum, _mm_and_si128(mask_sum, Q_vec));

            __m128i mask_diff = _mm_cmpgt_epi64(V1, V0);
            __m128i adj_diff = _mm_add_epi64(Vdiff, _mm_and_si128(mask_diff, Q_vec));

            // Store results back
            _mm_storel_epi64((__m128i*)&result[j], adj_sum);                 // result[j]
            _mm_storeh_pd((double*)&result[j + 2], _mm_castsi128_pd(adj_sum)); // result[j+2]
            _mm_storel_epi64((__m128i*)&result[j + 1], adj_diff);                   // result[j+1]
            _mm_storeh_pd((double*)&result[j + 3], _mm_castsi128_pd(adj_diff));     // result[j+3]
        }

        // Process remaining elements
        for (; j < n; j += 2) {
            auto t = result[j + 1];
            if (result[j] < t) {
                result[j + 1] = result[j] + Q - t;
            } else {
                result[j + 1] = result[j] - t;
            }
            result[j] += t;
            if (result[j] >= Q) {
                result[j] -= Q;
            }
        }
    }

    for (usint i = 1; i < u; i++) {
        l0 = 1 << i;
        l1 = 1 << (i + 1);
        d = n >> (i + 1);
        for (usint j = 0; j < n; j += l1) {
            usint k = 0;
            for (; k + 3 < l0; k += 4) {
                uint64_t y0 = mulmod(result[j + k + l0], rootTableAVX[k * d], Q, rootTablePreconAVX[k * d]);
                uint64_t y1 = mulmod(result[j + k + l0 + 1], rootTableAVX[(k + 1) * d], Q, rootTablePreconAVX[(k + 1) * d]);
                uint64_t y2 = mulmod(result[j + k + l0 + 2], rootTableAVX[(k + 2) * d], Q, rootTablePreconAVX[(k + 2) * d]);
                uint64_t y3 = mulmod(result[j + k + l0 + 3], rootTableAVX[(k + 3) * d], Q, rootTablePreconAVX[(k + 3) * d]);

                // Load result[j + k] to result[j + k + 3]
                __m256i rvec = _mm256_loadu_si256((__m256i*)&result[j + k]);

                // Set yvec = [y3, y2, y1, y0]
                __m256i yvec = _mm256_set_epi64x(y3, y2, y1, y0);

                // Compute sum = rvec + yvec
                __m256i sum = _mm256_add_epi64(rvec, yvec);

                // Compute mask for sum >= Q
                __m256i mask_sum = _mm256_cmpgt_epi64(sum, Q_minus_one_vec_256);

                // Adjust sum where sum >= Q: adj_sum = sum - Q where mask is true
                __m256i adj_sum = _mm256_sub_epi64(sum, _mm256_and_si256(mask_sum, Q_vec_256));

                // Store adjusted sum back to result[j + k] to result[j + k + 3]
                _mm256_storeu_si256((__m256i*)&result[j + k], adj_sum);

                // Compute diff = rvec - yvec
                __m256i diff = _mm256_sub_epi64(rvec, yvec);

                // Compute mask for rvec < yvec
                __m256i mask_diff = _mm256_cmpgt_epi64(yvec, rvec);

                // Adjust diff where rvec < yvec: adj_diff = diff + Q where mask is true
                __m256i adj_diff = _mm256_add_epi64(diff, _mm256_and_si256(mask_diff, Q_vec_256));

                // Store adjusted diff to result[j + k + l0] to result[j + k + l0 + 3]
                _mm256_storeu_si256((__m256i*)&result[j + k + l0], adj_diff);
            }

            for (; k < l0; k += 2) {
                uint64_t y0 = mulmod(result[j + k + l0], rootTableAVX[k * d], Q, rootTablePreconAVX[k * d]);
                uint64_t y1 = mulmod(result[j + k + l0 + 1], rootTableAVX[(k + 1) * d], Q, rootTablePreconAVX[(k + 1) * d]);
                // Load result[j + k] and result[j + k + 1]
                __m128i rvec = _mm_loadu_si128((__m128i*)&result[j + k]);

                // Set yvec = [y1, y0]
                __m128i yvec = _mm_set_epi64x(y1, y0);

                // Compute sum = rvec + yvec
                __m128i sum = _mm_add_epi64(rvec, yvec);

                // Compute mask for sum >= Q
                __m128i mask_sum = _mm_cmpgt_epi64(sum, Q_minus_one_vec);

                // Adjust sum where sum >= Q: adj_sum = sum - Q where mask is true
                __m128i adj_sum = _mm_sub_epi64(sum, _mm_and_si128(mask_sum, Q_vec));

                // Store adjusted sum back to result[j + k] and result[j + k + 1]
                _mm_storeu_si128((__m128i*)&result[j + k], adj_sum);

                // Compute diff = rvec - yvec
                __m128i diff = _mm_sub_epi64(rvec, yvec);

                // Compute mask for rvec < yvec
                __m128i mask_diff = _mm_cmpgt_epi64(yvec, rvec);

                // Adjust diff where rvec < yvec: adj_diff = diff + Q where mask is true
                __m128i adj_diff = _mm_add_epi64(diff, _mm_and_si128(mask_diff, Q_vec));

                // Store adjusted diff to result[j + k + l0] and result[j + k + l0 + 1]
                _mm_storeu_si128((__m128i*)&result[j + k + l0], adj_diff);
            }
        }
    }

    uint64_t z3 = rootTableAVX[n / 3], z32 = rootTableAVX[2 * n / 3];
    uint64_t z3precon = rootTablePreconAVX[n / 3], z32precon = rootTablePreconAVX[2 * n / 3];

    for (usint i = 0; i < v; i++) {
        l0 = U;
        l1 = U * 3;
        d = n / (3 * U);
        for (usint t = 0; t < i; t++) {
            l0 *= 3;
            l1 *= 3;
            d /= 3;
        }
        for (usint j = 0; j != n; j += l1) {
            usint k = 0;

            for (; k + 3 < l0; k += 4) {
                uint64_t y01 = mulmod(result[j + k + l0], rootTableAVX[k * d], Q, rootTablePreconAVX[k * d]);
                uint64_t y02 = mulmod(result[j + k + l0 * 2], rootTableAVX[k * d * 2], Q, rootTablePreconAVX[k * d * 2]);

                uint64_t y11 = mulmod(result[j + k + l0 + 1], rootTableAVX[(k + 1) * d], Q, rootTablePreconAVX[(k + 1) * d]);
                uint64_t y12 = mulmod(result[j + k + l0 * 2 + 1], rootTableAVX[(k + 1) * d * 2], Q, rootTablePreconAVX[(k + 1) * d * 2]);

                uint64_t y21 = mulmod(result[j + k + l0 + 2], rootTableAVX[(k + 2) * d], Q, rootTablePreconAVX[(k + 2) * d]);
                uint64_t y22 = mulmod(result[j + k + l0 * 2 + 2], rootTableAVX[(k + 2) * d * 2], Q, rootTablePreconAVX[(k + 2) * d * 2]);

                uint64_t y31 = mulmod(result[j + k + l0 + 3], rootTableAVX[(k + 3) * d], Q, rootTablePreconAVX[(k + 3) * d]);
                uint64_t y32 = mulmod(result[j + k + l0 * 2 + 3], rootTableAVX[(k + 3) * d * 2], Q, rootTablePreconAVX[(k + 3) * d * 2]);

                uint64_t y00 = y01 + y02;
                uint64_t y10 = y11 + y12;
                uint64_t y20 = y21 + y22;
                uint64_t y30 = y31 + y32;

                uint64_t w0 = mulmod(y01, z3, Q, z3precon) + mulmod(y02, z32, Q, z32precon);
                uint64_t w1 = mulmod(y11, z3, Q, z3precon) + mulmod(y12, z32, Q, z32precon);
                uint64_t w2 = mulmod(y21, z3, Q, z3precon) + mulmod(y22, z32, Q, z32precon);
                uint64_t w3 = mulmod(y31, z3, Q, z3precon) + mulmod(y32, z32, Q, z32precon);

                __m256i w_vec = _mm256_set_epi64x(w3, w2, w1, w0);
                __m256i w_mask = _mm256_cmpgt_epi64(w_vec, Q_minus_one_vec_256);
                w_vec = _mm256_sub_epi64(w_vec, _mm256_and_si256(w_mask, Q_vec_256));

                __m256i rvec = _mm256_loadu_si256((__m256i*)&result[j + k]);
                __m256i reg1 = _mm256_add_epi64(rvec, w_vec);
                __m256i mask_1 = _mm256_cmpgt_epi64(reg1, Q_minus_one_vec_256);
                reg1 = _mm256_sub_epi64(reg1, _mm256_and_si256(mask_1, Q_vec_256));
                _mm256_storeu_si256((__m256i*)&result[j + k + l0], reg1);

                __m256i y0vec = _mm256_set_epi64x(y30, y20, y10, y00);
                __m256i y0_mask = _mm256_cmpgt_epi64(y0vec, Q_minus_one_vec_256);
                y0vec = _mm256_sub_epi64(y0vec, _mm256_and_si256(y0_mask, Q_vec_256));

                __m256i y0wvec = _mm256_add_epi64(y0vec, w_vec);
                __m256i y0w_mask = _mm256_cmpgt_epi64(y0wvec, Q_minus_one_vec_256);
                y0wvec = _mm256_sub_epi64(y0wvec, _mm256_and_si256(y0w_mask, Q_vec_256));

                __m256i reg2 = _mm256_sub_epi64(rvec, y0wvec);
                __m256i mask_2 = _mm256_cmpgt_epi64(y0wvec, rvec);
                reg2 = _mm256_add_epi64(reg2, _mm256_and_si256(mask_2, Q_vec_256));
                _mm256_storeu_si256((__m256i*)&result[j + k + l0 + l0], reg2);

                __m256i reg0 = _mm256_add_epi64(rvec, y0vec);
                __m256i mask_0 = _mm256_cmpgt_epi64(reg0, Q_minus_one_vec_256);
                reg0 = _mm256_sub_epi64(reg0, _mm256_and_si256(mask_0, Q_vec_256 ));
                _mm256_storeu_si256((__m256i*)&result[j + k], reg0);
            }

            for (; k < l0; k += 2) {
                uint64_t y01 = mulmod(result[j + k + l0], rootTableAVX[k * d], Q, rootTablePreconAVX[k * d]);
                uint64_t y02 = mulmod(result[j + k + l0 * 2], rootTableAVX[k * d * 2], Q, rootTablePreconAVX[k * d * 2]);

                uint64_t y11 = mulmod(result[j + k + l0 + 1], rootTableAVX[(k + 1) * d], Q, rootTablePreconAVX[(k + 1) * d]);
                uint64_t y12 = mulmod(result[j + k + l0 * 2 + 1], rootTableAVX[(k + 1) * d * 2], Q, rootTablePreconAVX[(k + 1) * d * 2]);

                uint64_t y00 = y01 + y02;
                uint64_t y10 = y11 + y12;

                uint64_t w0 = mulmod(y01, z3, Q, z3precon) + mulmod(y02, z32, Q, z32precon);
                uint64_t w1 = mulmod(y11, z3, Q, z3precon) + mulmod(y12, z32, Q, z32precon);

                __m128i w_vec = _mm_set_epi64x(w1, w0);
                __m128i w_mask = _mm_cmpgt_epi64(w_vec, Q_minus_one_vec);
                w_vec = _mm_sub_epi64(w_vec, _mm_and_si128(w_mask, Q_vec));

                __m128i rvec = _mm_loadu_si128((__m128i*)&result[j + k]);
                __m128i reg1 = _mm_add_epi64(rvec, w_vec);
                __m128i mask_1 = _mm_cmpgt_epi64(reg1, Q_minus_one_vec);
                reg1 = _mm_sub_epi64(reg1, _mm_and_si128(mask_1, Q_vec));
                _mm_storeu_si128((__m128i*)&result[j + k + l0], reg1);

                __m128i y0vec = _mm_set_epi64x(y10, y00);
                __m128i y0_mask = _mm_cmpgt_epi64(y0vec, Q_minus_one_vec);
                y0vec = _mm_sub_epi64(y0vec, _mm_and_si128(y0_mask, Q_vec));

                __m128i y0wvec = _mm_add_epi64(y0vec, w_vec);
                __m128i y0w_mask = _mm_cmpgt_epi64(y0wvec, Q_minus_one_vec);
                y0wvec = _mm_sub_epi64(y0wvec, _mm_and_si128(y0w_mask, Q_vec));

                __m128i reg2 = _mm_sub_epi64(rvec, y0wvec);
                __m128i mask_2 = _mm_cmpgt_epi64(y0wvec, rvec);
                reg2 = _mm_add_epi64(reg2, _mm_and_si128(mask_2, Q_vec));
                _mm_storeu_si128((__m128i*)&result[j + k + l0 + l0], reg2);

                __m128i reg0 = _mm_add_epi64(rvec, y0vec);
                __m128i mask_0 = _mm_cmpgt_epi64(reg0, Q_minus_one_vec);
                reg0 = _mm_sub_epi64(reg0, _mm_and_si128(mask_0, Q_vec));
                _mm_storeu_si128((__m128i*)&result[j + k], reg0);
            }
        }
    }
}

template <typename VecType>
VecType primecyc::RaderFFTNat<VecType>::ForwardRader(const VecType& element, const IntType& rootOfUnity) {
    usint tot = element.GetLength();
    
    auto modulus = element.GetModulus();
    auto order = tot + 1;

    if (m_forwardPermutation.find(order) == m_forwardPermutation.end()) {
        PreComputeIsomorphism(order);
    }

    if (m_rootTableByModulusRoot.find({modulus, rootOfUnity}) == m_rootTableByModulusRoot.end()) {
        PreComputeRootTable(order, {modulus, rootOfUnity});
    }

    const auto &indices = m_inversePermutation[order];
    const auto &rootsT = m_rootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &rootsTPrecon = m_rootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    std::vector<IntType> temp(tot);
    std::vector<IntType> out(tot);

    temp[0] = IntType(0).ModSub(element[0], modulus);

#pragma GCC ivdep
    for (usint i = 1; i < tot; i++) {
        temp[i] = temp[0] + element[i];
        if (temp[i] >= modulus) {
            temp[i] -= modulus;
        }
    }

    auto rootOfUnityTot = rootOfUnity.ModExp(order, modulus);

#pragma GCC ivdep
    for (usint i = 0; i < tot; i++) {
        out[indices[i]] = temp[i];
    }

    ForwardFFTBase2n3(out, modulus, rootOfUnityTot, temp);

#pragma GCC ivdep
    for (usint i = 0; i < tot; i++) {
        temp[i].ModMulFastConstEq(rootsT[i], modulus, rootsTPrecon[i]);
    }

    ForwardFFTBase2n3(temp, modulus, rootOfUnityTot.ModExp(tot - 1, modulus), out);

    const auto &forward = m_forwardPermutation[order];
#pragma GCC ivdep
    for (usint i = 0; i < tot; i++) {
        if (forward[i] == 0) {
            temp[tot - 1] = out[i];
        } else {
            temp[forward[i] - 1] = out[(tot - i) % tot];
        }
    }

    VecType result(tot, modulus);
    for (usint i = 0; i < tot; i++) {
        result[i] = temp[i];
    }

    return result;
}

template <typename VecType>
VecType primecyc::RaderFFTNat<VecType>::ForwardRaderPermute(const VecType& element, const IntType& rootOfUnity) {
    usint tot = element.GetLength();
    
    auto modulus = element.GetModulus();
    auto order = tot + 1;

    if (m_forwardPermutation.find(order) == m_forwardPermutation.end()) {
        PreComputeIsomorphism(order);
    }

    if (m_rootTableByModulusRoot.find({modulus, rootOfUnity}) == m_rootTableByModulusRoot.end()) {
        PreComputeRootTable(order, {modulus, rootOfUnity});
    }

    const auto &indices = m_inversePermutation[order];
    const auto &rootsT = m_rootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &rootsTPrecon = m_rootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    auto temp = VecType(tot, modulus);
    auto out = VecType(tot, modulus);

    out[0] = element[tot - 1];
    for (usint i = 1; i < tot; i++) {
        out[i] = element[i - 1];
    }

    auto rootOfUnityTot = rootOfUnity.ModExp(order, modulus);

    for (usint i = 0; i < tot; i++) {
        temp[indices[i]] = out[i];
    }

    ForwardFFTBase2n3(temp, rootOfUnityTot, &out);

    for (usint i = 0; i < tot; i++) {
        out[i].ModMulFastConstEq(rootsT[i], modulus, rootsTPrecon[i]);
    }

    ForwardFFTBase2n3(out, rootOfUnityTot.ModExp(tot - 1, modulus), &temp);

    const auto &forward = m_forwardPermutation[order];

    for (usint i = 0; i < tot; i++) {
        if (forward[i] == 0) {
            out[tot - 1] = temp[i];
        } else {
            out[forward[i] - 1] = temp[(tot - i) % tot];
        }
    }

    return out;
}

template <typename VecType>
VecType primecyc::RaderFFTNat<VecType>::InverseRader(const VecType& element, const IntType& rootOfUnity) {
    usint tot = element.GetLength();

    auto modulus = element.GetModulus();
    auto order = tot + 1;

    if (m_forwardPermutation.find(order) == m_forwardPermutation.end()) {
        PreComputeIsomorphism(order);
    }

    if (m_rootTableByModulusRoot.find({modulus, rootOfUnity}) == m_rootTableByModulusRoot.end()) {
        PreComputeRootTable(order, {modulus, rootOfUnity});
    }

    const auto &indices = m_inversePermutation[order];
    const auto &invRootsT = m_inverseRootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &invRootsTPrecon = m_inverseRootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    const auto &forward = m_forwardPermutation[order];

    auto temp = VecType(tot, modulus);
    auto out = VecType(tot, modulus);

    auto rootOfUnityTot = rootOfUnity.ModExp(order, modulus);

    for (usint i = 0; i < tot; i++) {
        if (forward[i] == 0) {
            temp[i] = element[tot - 1];
        } else {
            temp[(tot - i) % tot] = element[forward[i] - 1];
        }
    }

    ForwardFFTBase2n3(temp, rootOfUnityTot, &out);

    for (usint i = 0; i < tot; i++) {
        out[i].ModMulFastConstEq(invRootsT[i], modulus, invRootsTPrecon[i]);
    }

    ForwardFFTBase2n3(out, rootOfUnityTot.ModExp(tot - 1, modulus), &temp);

    for (usint i = 0; i < tot; i++) {
        out[i] = temp[indices[i]];
    }

    out[0] = IntType(0).ModSubFast(out[0], modulus);
    for (usint i = 1; i < tot; i++) {
        out[i] = out[0].ModAddFast(out[i], modulus);
    }

    return out;
}

template <typename VecType>
VecType primecyc::RaderFFTNat<VecType>::InverseRaderPermute(const VecType& element, const IntType& rootOfUnity) {
    usint tot = element.GetLength();

    auto modulus = element.GetModulus();
    auto order = tot + 1;

    if (m_forwardPermutation.find(order) == m_forwardPermutation.end()) {
        PreComputeIsomorphism(order);
    }

    if (m_rootTableByModulusRoot.find({modulus, rootOfUnity}) == m_rootTableByModulusRoot.end()) {
        PreComputeRootTable(order, {modulus, rootOfUnity});
    }

    const auto &indices = m_inversePermutation[order];
    const auto &invRootsT = m_inverseRootTableByModulusRoot[{modulus, rootOfUnity}];
    const auto &invRootsTPrecon = m_inverseRootPreconTableByModulusRoot[{modulus, rootOfUnity}];

    auto temp = VecType(tot, modulus);
    auto out = VecType(tot, modulus);

    const auto &forward = m_forwardPermutation[order];

    auto rootOfUnityTot = rootOfUnity.ModExp(order, modulus);

    for (usint i = 0; i < tot; i++) {
        if (forward[i] == 0) {
            temp[i] = element[tot - 1];
        } else {
            temp[(tot - i) % tot] = element[forward[i] - 1];
        }
    }

    ForwardFFTBase2n3(temp, rootOfUnityTot, &out);

    for (usint i = 0; i < tot; i++) {
        out[i].ModMulFastConstEq(invRootsT[i], modulus, invRootsTPrecon[i]);
    }

    ForwardFFTBase2n3(out, rootOfUnityTot.ModExp(tot - 1, modulus), &temp);

    for (usint i = 1; i < tot; i++) {
        out[i - 1] = temp[indices[i]];
    }
    out[tot - 1] = temp[indices[0]];

    return out;
}

#endif // PRIME_CYC_IMPL
