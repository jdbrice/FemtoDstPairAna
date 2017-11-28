#ifndef SINGLE_TRACK_SKIMMER_H
#define SINGLE_TRACK_SKIMMER_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "vendor/loguru.h"

#include "Filters/MtdTrackFilter.h"
#include "Filters/MuonMLPFilter.h"


class SingleTrackSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	MtdTrackFilter _trackFilter;
	MuonMLPFilter _mlp;

	int nPos, znPos;
	int nNeg, znNeg;
	float nPosMLP, nNegMLP;



public:
	virtual const char* classname() const {return "SingleTrackSkimmer";}
	SingleTrackSkimmer() {}
	~SingleTrackSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		_trackFilter.load( config, nodePath + ".TrackFilter" );
		_mlp.load( config, nodePath + ".MuonMLPFilter" );
	}

protected:

	virtual void analyzeEvent(){

		_event = _rEvent.get();
		
		size_t nTracks = _rTracks.N();

		nPos = 0;
		nNeg = 0;

		nPosMLP = 0.0;
		nNegMLP = 0.0;

		znPos = 0;
		znNeg = 0;

		float dvz = _event->mPrimaryVertex_mX3 - _event->mWeight;

		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			
			if ( false == _trackFilter.pass( _proxy ) )
				continue;
			float r = _mlp.evaluate( _proxy );

			FemtoTrack * t = _proxy._track;
			book->fill( "mPt", t->mPt * t->charge() );
			book->fill( "mEta", t->mEta );
			book->fill( "mPhi", t->mPhi );
			book->fill( "mlp", r );
			book->fill( "mlp_vs_mPt", t->mPt * t->charge(), r );

			if ( _proxy._track->charge() > 0 ){
				nPos++;
				nPosMLP += r;
				book->fill( "mlp_pos_vs_vz", _event->mPrimaryVertex_mX3, r );
				book->fill( "mlp_pos_vs_grefmult", _event->mGRefMult, r );
				book->fill( "mlp_pos_vs_dvz", dvz, r );
			}
			if ( _proxy._track->charge() < 0 ){
				nNeg++;
				nNegMLP += r;
				book->fill( "mlp_neg_vs_vz", _event->mPrimaryVertex_mX3, r );
				book->fill( "mlp_neg_vs_grefmult", _event->mGRefMult, r );
				book->fill( "mlp_neg_vs_dvz", dvz, r );
			}
		}


		book->fill( "nPos", nPos );
		book->fill( "nNeg", nNeg );
		book->fill( "nPos_nNeg", nNeg, nPos );

		book->fill( "nPos_vs_grefmult", _event->mGRefMult, nPos );
		book->fill( "nNeg_vs_grefmult", _event->mGRefMult, nNeg );

		book->fill( "nPos_vs_vz", _event->mPrimaryVertex_mX3, nPos );
		book->fill( "nNeg_vs_vz", _event->mPrimaryVertex_mX3, nNeg );

	} // analyse Event

};

#endif