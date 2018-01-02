

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "SameEventSkimmer/SameEventSkimmer.h"
#include "SameEventSkimmer/McSkimmer.h"
#include "SameEventSkimmer/SingleTrackSkimmer.h"
#include "SameEventSkimmer/PidQAMaker.h"
#include "MixedEventSkimmer/MixedEventAnalyzer.h"
#include "SingleTrackAna/SingleTrackAna.h"
#include "SingleTrackAna/WeightedTrackMeans.h"
#include "SameEventSkimmer/SignalRebinner.h"
#include "SameEventSkimmer/EfficiencyTableMaker.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);


	TaskFactory::registerTaskRunner<FemtoDstSkimmer>( "FemtoDstSkimmer" );
	TaskFactory::registerTaskRunner<MixedEventAnalyzer>( "MixedEventAnalyzer" );
	TaskFactory::registerTaskRunner<SameEventSkimmer>( "SameEventSkimmer" );
	TaskFactory::registerTaskRunner<SingleTrackSkimmer>( "SingleTrackSkimmer" );
	TaskFactory::registerTaskRunner<SingleTrackAna>( "SingleTrackAna" );
	TaskFactory::registerTaskRunner<WeightedTrackMeans>( "WeightedTrackMeans" );
	TaskFactory::registerTaskRunner<SignalRebinner>( "SignalRebinner" );
	TaskFactory::registerTaskRunner<PidQAMaker>( "PidQAMaker" );

	TaskFactory::registerTaskRunner<McSkimmer>( "McSkimmer" );
	TaskFactory::registerTaskRunner<EfficiencyTableMaker>( "EfficiencyTableMaker" );

	TaskEngine engine( argc, argv );

	return 0;
}
