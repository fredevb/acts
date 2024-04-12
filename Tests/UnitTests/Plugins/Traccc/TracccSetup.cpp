// Project include(s).
#include "traccc/utils/algorithm.hpp"

// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

// traccc plugin
#include "Acts/Plugins/Traccc/TracccConverter.hpp"

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
    BOOST_TEST(acts::traccc_conversion::test());
}