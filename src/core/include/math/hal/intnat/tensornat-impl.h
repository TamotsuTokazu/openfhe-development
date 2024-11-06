#ifndef _TENSORNAT_IMPL_H_
#define _TENSORNAT_IMPL_H_

#include "math/math-hal.h"
#include "math/nbtheory.h"

#include "math/hal/intnat/tensornat.h"

template <typename VecType>
std::map<usint, std::pair<usint, usint>> primecyc::TensorFFTNat<VecType>::m_factor;

template <typename VecType>
VecType primecyc::TensorFFTNat<VecType>::ForwardTransform(const VecType& element, const IntType& root, usint cycloOrder) {
    usint tot = element.GetLength();

    const auto& modulus = element.GetModulus();
    auto [p, q] = m_factor[cycloOrder];

    IntType rootp = root.ModExp(q * (q - 1), modulus);
    IntType rootq = root.ModExp(p * (p - 1), modulus);

    VecType result(tot, modulus);

    VecType tempp(p - 1, modulus);
    VecType tempq(q - 1, modulus);

    for (usint i = 0; i < p - 1; i++) {
        for (usint j = 0; j < q - 1; j++) {
            tempq[j] = element[(q - 1) * i + j];
        }
        tempq = RaderFFTNat<VecType>().ForwardRaderPermute(tempq, rootq);
        for (usint j = 0; j < q - 1; j++) {
            result[(q - 1) * i + j] = tempq[j];
        }
    }

    for (usint j = 0; j < q - 1; j++) {
        for (usint i = 0; i < p - 1; i++) {
            tempp[i] = result[(q - 1) * i + j];
        }
        tempp = RaderFFTNat<VecType>().ForwardRaderPermute(tempp, rootp);
        for (usint i = 0; i < p - 1; i++) {
            result[(q - 1) * i + j] = tempp[i];
        }
    }

    return result;
}

template <typename VecType>
VecType primecyc::TensorFFTNat<VecType>::InverseTransform(const VecType& element, const IntType& root, usint cycloOrder) {
    usint tot = element.GetLength();

    const auto& modulus = element.GetModulus();
    auto [p, q] = m_factor[cycloOrder];

    IntType rootp = root.ModExp(q * (q - 1), modulus);
    IntType rootq = root.ModExp(p * (p - 1), modulus);

    VecType result(tot, modulus);

    VecType tempp(p - 1, modulus);
    VecType tempq(q - 1, modulus);

    for (usint j = 0; j < q - 1; j++) {
        for (usint i = 0; i < p - 1; i++) {
            tempp[i] = element[(q - 1) * i + j];
        }
        tempp = RaderFFTNat<VecType>().InverseRaderPermute(tempp, rootp);
        for (usint i = 0; i < p - 1; i++) {
            result[(q - 1) * i + j] = tempp[i];
        }
    }

    for (usint i = 0; i < p - 1; i++) {
        for (usint j = 0; j < q - 1; j++) {
            tempq[j] = result[(q - 1) * i + j];
        }
        tempq = RaderFFTNat<VecType>().InverseRaderPermute(tempq, rootq);
        for (usint j = 0; j < q - 1; j++) {
            result[(q - 1) * i + j] = tempq[j];
        }
    }

    return result;
}





#endif // _TENSORNAT_IMPL_H_