#ifndef _TENSORNAT_H_
#define _TENSORNAT_H_

#include "utils/utilities.h"

namespace primecyc {

template <typename IntType>
using ModulusRoot = std::pair<IntType, IntType>;

template <typename VecType>
class TensorFFTNat {
    using IntType = typename VecType::Integer;

public:
    static std::map<usint, std::pair<usint, usint>> m_factor;

    VecType ForwardTransform(const VecType& element, const IntType& root, usint cycloOrder);

    VecType InverseTransform(const VecType& element, const IntType& root, usint cycloOrder);
};

}

#include "math/hal/intnat/tensornat-impl.h"

#endif // _TENSORNAT_H_