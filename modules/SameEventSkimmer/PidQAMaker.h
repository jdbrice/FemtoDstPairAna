#ifndef PID_QA_MAKER_H
#define PID_QA_MAKER_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "vendor/loguru.h"

#include "Filters/TrackFilter.h"
#include "Filters/MuonMlpFilter.h"


class PidQAMaker : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;
	FemtoTrackProxy _proxy2;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	TrackFilter _trackFilter;
	MuonMLPFilter _mlp;


public:
	virtual const char* classname() const {return "PidQAMaker";}
	PidQAMaker() {}
	~PidQAMaker() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		_trackFilter.load( config, nodePath + ".TrackFilter" );
		_mlp.load( config, nodePath + ".MuonMLPFilter" );
	}

protected:

	virtual void fill( 	FemtoTrackProxy &_ltp, FemtoTrackProxy &_sltp, string prefix = "" ){

		TLorentzVector lv1 = _ltp._track->lv( 0.105 );
		TLorentzVector lv2 = _sltp._track->lv( 0.105 );
		TLorentzVector lv = lv1 + lv2;

		book->fill( prefix + "mlpl_vs_mass", lv.M(), _ltp._pid );
		book->fill( prefix + "mlpsl_vs_mass", lv.M(), _sltp._pid );
		book->fill( prefix + "mlpsum_vs_mass", lv.M(), _sltp._pid + _ltp._pid );
		book->fill( prefix + "mlp_vs_mass", lv.M(), _sltp._pid, 1 );
		book->fill( prefix + "mlp_vs_mass", lv.M(), _ltp._pid, 1 );
	}

	virtual void analyze_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		int chargeSum = _tp1._track->charge() + _tp2._track->charge();
		FemtoTrackProxy *ltp = &_tp1; 		// leading pt 
		FemtoTrackProxy *sltp = &_tp2; 		// sub-leading pt

		if ( fabs(_tp2._track->mPt) > fabs( _tp1._track->mPt ) ){
			ltp  = &_tp2;
			sltp = &_tp1;
		}

		if ( chargeSum == 2 ){
			fill( *ltp, *sltp, "lsp_" );
		}
		else if ( chargeSum == -2 ){
			fill( *ltp, *sltp, "lsn_" );
		}
		else if ( chargeSum == 0 ){
			fill( *ltp, *sltp, "uls_" );
		}
	}

	virtual void analyzeEvent(){

		_event = _rEvent.get();
		
		size_t nTracks = _rTracks.N();
		if ( nTracks < 2 )
			return;

		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			if ( false == _trackFilter.pass( _proxy ) ){
				// LOG_F( INFO, "FAIL TRACK" );
				continue;
			}

			_proxy._pid =  _mlp.evaluate( _proxy );

			string prefix = "pos_";
			if ( _proxy._track->charge() < 0 )
				prefix = "neg_";

			book->fill( prefix + "mlp_vs_pT", _proxy._track->mPt, _proxy._pid );



			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rTracks, _rMtdPid );
				_proxy2._pid = _mlp.evaluate( _proxy2 );
				if ( false == _trackFilter.pass( _proxy2 ) )
					continue;
				analyze_pair( _proxy, _proxy2 );
			} // loop on tracks j

		} // loop on tracks i

	} // analyse Event

};

#endif