// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Traccc/TracccChainAlgorithm.hpp"
#include "ActsExamples/Traccc/Chain.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
namespace py = pybind11;

template class ActsExamples::TracccChainAlgorithm<ActsExamples::Chain::Host>;

namespace Acts::Python {

template <typename platform_t>
void declareTracccAlgorithm(py::module &m, const std::string &platformName) {
  using config_type = typename ActsExamples::Chain::Chain<platform_t>::Config;
  py::class_<config_type, std::shared_ptr<config_type>>(m, (std::string("TracccChainConfig") + platformName).c_str())
  .def(py::init<>());
  
  using algorithm_type = typename ActsExamples::TracccChainAlgorithm<platform_t>;
  ACTS_PYTHON_DECLARE_ALGORITHM(
    algorithm_type, m,
    (std::string("TracccChainAlgorithm") + platformName).c_str(), inputCells,
    outputTracks, trackingGeometry, field, digitizationConfigs, chainConfig);
}

void addTracccChain(Context& ctx) {
  auto m = ctx.get("examples");
  declareTracccAlgorithm<ActsExamples::Chain::Host>(m, "Host");
}

}  // namespace Acts::Python
