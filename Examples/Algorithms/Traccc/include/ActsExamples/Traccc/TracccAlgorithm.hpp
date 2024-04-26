// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Construct a traccc algorithm
class TracccAlgorithm final : public IAlgorithm {
public:
struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// (optional) Input clusters for each measurement
    std::string inputClusters;
    /// Output fitted tracks collection.
    std::string outputTracks;
    /// Type erased fitter function.
    std::shared_ptr<TrackFitterFunction> fit;
    /// Pick a single track for debugging (-1 process all tracks)
    int pickTrack = -1;
    // Type erased calibrator for the measurements
    std::shared_ptr<MeasurementCalibrator> calibrator;
};

/// Construct the traccc algorithm.
///
/// @param cfg is the algorithm configuration
/// @param lvl is the logging level
TracccAlgorithm(Config cfg, Acts::Logging::Level lvl);

/// Run the algorithm.
///
/// @param ctx is the algorithm context with event information
/// @return a process code indication success or failure
ProcessCode execute(const AlgorithmContext& ctx) const final override;

/// Const access to the config
const Config& config() const { return m_cfg; }

private:
Config m_cfg;

ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                        "InputMeasurements"};
ReadDataHandle<IndexSourceLinkContainer> m_inputSourceLinks{
    this, "InputSourceLinks"};
ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                        "InputProtoTracks"};
ReadDataHandle<TrackParametersContainer> m_inputInitialTrackParameters{
    this, "InputInitialTrackParameters"};

ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};

WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples