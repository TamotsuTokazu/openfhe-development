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

#include "math/hal/intnat/primecycnat.h"

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_bitReverseTableBase2n3;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_base2n3RootPreconTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_base2n3RootTableByModulusRoot;

template <typename VecType>
std::map<usint, std::array<usint, 4>> primecyc::RaderFFTNat<VecType>::m_Base2n3Info;

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_forwardPermutation;

template <typename VecType>
std::map<usint, std::vector<usint>> primecyc::RaderFFTNat<VecType>::m_inversePermutation;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_rootTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_rootPreconTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_inverseRootTableByModulusRoot;

template <typename VecType>
std::map<primecyc::ModulusRoot<typename VecType::Integer>, VecType> primecyc::RaderFFTNat<VecType>::m_inverseRootPreconTableByModulusRoot;

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

    VecType rootTable(tot, nttModulusRoot.first);
    VecType rootPreconTable(tot, nttModulusRoot.first);
    VecType inverseRootTable(tot, nttModulusRoot.first);
    VecType inverseRootPreconTable(tot, nttModulusRoot.first);

    for (usint i = 1; i <= tot; i++) {
        if (indices[i] == 0) {
            rootTable[0] = z;
        } else {
            rootTable[tot - indices[i]] = z;
        }
        z.ModMulFastEq(rootOfUnityOrder, nttModulusRoot.first);
    }

    VecType rootTableT(tot, nttModulusRoot.first);
    ForwardFFTBase2n3(rootTable, rootOfUnityTot, &rootTableT);

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

    VecType rootTable(order, nttModulusRoot.first);
    VecType rootPreconTable(order, nttModulusRoot.first);

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
    if (result->GetLength() != n) {
        OPENFHE_THROW(lbcrypto::math_error, "size of input element and size of output element not of same size");
    }

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
        l1 *= 2;
        d /= 2;
        for (usint j = 0; j != n; j += l1) {
            for (usint k = 0; k < l0; k++) {
                const auto &o = rootTable[k * d], &oPrecon = rootTablePrecon[k * d];
                auto &a0 = (*result)[j + k], &a1 = (*result)[j + k + l0];
                IntType y1 = a1.ModMulFastConst(o, modulus, oPrecon);
                a1 = a0.ModSubFast(y1, modulus);
                a0.ModAddFastEq(y1, modulus);
            }
        }
        l0 *= 2;
    }

    IntType z3 = rootOfUnity.ModExp(n / 3, modulus), z32 = z3.ModMulFast(z3, modulus);
    IntType z3precon = z3.PrepModMulConst(modulus), z32precon = z32.PrepModMulConst(modulus);

    for (usint i = 0; i < v; i++) {
        l1 *= 3;
        d /= 3;
        for (usint j = 0; j != n; j += l1) {
            for (usint k = 0; k < l0; k++) {
                const auto &o = rootTable[k * d], &o2 = rootTable[k * d * 2];
                const auto &oPrecon = rootTablePrecon[k * d], &o2Precon = rootTablePrecon[k * d * 2];
                IntType &a0 = (*result)[j + k], &a1 = (*result)[j + k + l0], &a2 = (*result)[j + k + 2 * l0];
                IntType y1 = a1.ModMulFastConst(o, modulus, oPrecon), y2 = a2.ModMulFastConst(o2, modulus, o2Precon), y0 = y1.ModAddFast(y2, modulus);
                IntType w = y1.ModMulFastConst(z3, modulus, z3precon).ModAddFast(y2.ModMulFastConst(z32, modulus, z32precon), modulus);
                a1 = a0.ModAddFast(w, modulus);
                a2 = a0.ModSubFast(y0.ModAddFast(w, modulus), modulus);
                a0.ModAddFastEq(y0, modulus);
            }
        }
        l0 *= 3;
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

    auto temp = VecType(tot, modulus);
    auto out = VecType(tot, modulus);

    temp[0] = IntType(0).ModSub(element[0], modulus);
    for (usint i = 1; i < tot; i++) {
        temp[i] = temp[0].ModAddFast(element[i], modulus);
    }

    auto rootOfUnityTot = rootOfUnity.ModExp(order, modulus);

    for (usint i = 0; i < tot; i++) {
        out[indices[i]] = temp[i];
    }

    ForwardFFTBase2n3(out, rootOfUnityTot, &temp);

    for (usint i = 0; i < tot; i++) {
        temp[i].ModMulFastConstEq(rootsT[i], modulus, rootsTPrecon[i]);
    }

    ForwardFFTBase2n3(temp, rootOfUnityTot.ModExp(tot - 1, modulus), &out);

    const auto &forward = m_forwardPermutation[order];

    for (usint i = 0; i < tot; i++) {
        if (forward[i] == 0) {
            temp[tot - 1] = out[i];
        } else {
            temp[forward[i] - 1] = out[(tot - i) % tot];
        }
    }

    return temp;
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

#endif
