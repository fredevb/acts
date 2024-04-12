// covfie core
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// covfie plugin
#include "Acts/Plugins/Covfie/CovfieConverter.hpp"

// Boost.Test include(s)
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Covfie_ConstantField1D) {
    using field_t =
        covfie::field<covfie::backend::constant<covfie::vector::float1,
                                                covfie::vector::float1>>;

    field_t f(
        covfie::make_parameter_pack(field_t::backend_t::configuration_t{2.f}));
    field_t::view_t v(f);

    for (float x = -100.f; x <= 100.f; x += 1.f) {
        BOOST_TEST(v.at(x)[0] == 2.f);
    }
}

BOOST_AUTO_TEST_CASE(Covfie_Convert){
    BOOST_TEST(acts::covfie_conversion::test());
}