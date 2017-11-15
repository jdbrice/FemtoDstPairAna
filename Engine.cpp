

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "SameEventSkimmer/SameEventSkimmer.h"
#include "SameEventSkimmer/PidQAMaker.h"
#include "MixedEventSkimmer/MixedEventAnalyzer.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);


	TaskFactory::registerTaskRunner<FemtoDstSkimmer>( "FemtoDstSkimmer" );
	TaskFactory::registerTaskRunner<MixedEventAnalyzer>( "MixedEventAnalyzer" );
	TaskFactory::registerTaskRunner<SameEventSkimmer>( "SameEventSkimmer" );
	TaskFactory::registerTaskRunner<PidQAMaker>( "PidQAMaker" );

	TaskEngine engine( argc, argv );

	return 0;
}
