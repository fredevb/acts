// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Project include(s).
#include "traccc/utils/algorithm.hpp"

// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

// traccc plugin
#include "Acts/Plugins/Traccc/TracccConversion.hpp"

// System include(s).
#include <string>

class DoubleInt : public traccc::algorithm<int(const int &)> {
public:
    virtual int operator()(const int &i) const override { return 2 * i; }
};

BOOST_AUTO_TEST_CASE(algorithm_basic) {
    DoubleInt i;

    BOOST_TEST(i(1) == 2);
    BOOST_TEST(i(-1) == -2);
    BOOST_TEST(i(5) == 10);
}

BOOST_AUTO_TEST_CASE(Traccc_Convert){
    BOOST_TEST(Acts::TracccConversion::test());
}
