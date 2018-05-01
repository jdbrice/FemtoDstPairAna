#ifndef K0S_ANA_H
#define K0S_ANA_H


#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "vendor/loguru.h"

class K0SAna : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;
	FemtoTrackProxy _proxy2;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;
public:
	K0SAna() {}
	~K0SAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

	}

protected:

	virtual void analyzeEvent(){

		_event = _rEvent.get();


		size_t nTracks = _rTracks.N();
		for (size_t i = 0; i < nTracks; i++ ){
		}
	}
	
};


#endif