
#ifndef PRIME_CYC
#define PRIME_CYC

#include "utils/debug.h"
#include "utils/exception.h"
#include "utils/utilities.h"

#include <map>
#include <vector>

namespace primecycutil {

inline usint reverseBits(usint num, usint bits, usint base) {
    usint ans = 0;
    for (uint i = 0; i < bits; i++) {
        ans = ans * base + num % base;
        num /= base;
    }
    return ans;
}

// assert modulus is 32-bit
inline usint powMod(usint base, usint exponent, usint modulus) {
    uint64_t t = base, result = 1;
    for (; exponent; exponent >>= 1) {
        if (exponent & 1) {
            result = (result * t) % modulus;
        }
        t = t * t % modulus;
    }
    return (usint)result;
}

inline usint findPrimitiveRoot(usint modulus) {
    std::vector<usint> primeFactors;
    usint groupOrder = modulus - 1;
    for (usint t = groupOrder, i = 2; i <= t; i++) {
        if (t % i == 0) {
            primeFactors.push_back(i);
            while (t % i == 0) {
                t /= i;
            }
        }
    }
    for (auto &i : primeFactors) {
        i = groupOrder / i;
    }
    for (usint i = 2; i < modulus; i++) {
        bool flag = true;
        for (auto j : primeFactors) {
            if (powMod(i, j, modulus) == 1) {
                flag = false;
                break;
            }
        }
        if (flag) {
            return i;
        }
    }
    return 0;
}

}

namespace primecyc {

template <typename IntType>
using ModulusRoot = std::pair<IntType, IntType>;

template <typename VecType>
class RaderFFTNat {
    using IntType = typename VecType::Integer;

public:

    static std::map<usint, std::vector<usint>> m_bitReverseTableBase2n3;
    static std::map<ModulusRoot<IntType>, VecType> m_base2n3RootTableByModulusRoot;
    static std::map<ModulusRoot<IntType>, VecType> m_base2n3RootPreconTableByModulusRoot;
    static std::map<usint, std::array<usint, 4>> m_Base2n3Info;

    static std::map<usint, std::vector<usint>> m_forwardPermutation;
    static std::map<usint, std::vector<usint>> m_inversePermutation;

    static std::map<ModulusRoot<IntType>, VecType> m_rootTableByModulusRoot;
    static std::map<ModulusRoot<IntType>, VecType> m_rootPreconTableByModulusRoot;

    static std::map<ModulusRoot<IntType>, VecType> m_inverseRootTableByModulusRoot;
    static std::map<ModulusRoot<IntType>, VecType> m_inverseRootPreconTableByModulusRoot;

    static std::map<usint, bool> m_enabled;

    void PreComputeIsomorphism(usint cycloOrder);

    void PreComputeRootTable(usint cycloOrder, const ModulusRoot<IntType>& nttModulusRoot);

    void PreComputeBitReverseTableBase2n3(usint cycloOrder);

    void PreComputeBase2n3RootTable(usint cycloOrder, const ModulusRoot<IntType>& nttModulusRoot);

    void ForwardFFTBase2n3(const VecType& element, const IntType& rootOfUnity, VecType* result);

    VecType ForwardRader(const VecType& element, const IntType& rootOfUnity);
    VecType ForwardRaderPermute(const VecType& element, const IntType& rootOfUnity);

    VecType InverseRader(const VecType& element, const IntType& rootOfUnity);
    VecType InverseRaderPermute(const VecType& element, const IntType& rootOfUnity);
};

}

#include "math/hal/intnat/primecycnat-impl.h"

#endif