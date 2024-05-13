// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"

// System include(s)
#include <cstddef>

namespace Acts::TracccPlugin::Utils {

template <std::size_t N, typename dvector_t>
Acts::ActsVector<N> newVector(const dvector_t& dvec){
    Acts::ActsVector<N> res;
    for(std::size_t i = 0; i < N; i++){
        res(i) = static_cast<Acts::ActsScalar>(dvec[0][i]);
    }
    return res;
}

template <std::size_t N, typename matrixNxN_t>
Acts::ActsSquareMatrix<N> newSqaureMatrix(const matrixNxN_t& mat){
    Acts::ActsSquareMatrix<N> res;
    for(std::size_t x = 0; x < N; x++){
        for(std::size_t y = 0; y < N; y++){
            res(x,y) = static_cast<Acts::ActsScalar>(mat[x][y]);
        }
    }
    return res;
}

}