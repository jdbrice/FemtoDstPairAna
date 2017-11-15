#ifndef SAME_EVENT_SKIMMER_H
#define SAME_EVENT_SKIMMER_H

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
#include "Filters/MuonMlpFilter.h"


class SameEventSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;
	FemtoTrackProxy _proxy2;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	MtdTrackFilter _trackFilter;
	MuonMLPFilter _mlp;

	int sig_nPos;
	int sig_nNeg;
	int nSig_pairs;

	int nPos;
	int nNeg;
	int nPairs;



public:
	virtual const char* classname() const {return "SameEventSkimmer";}
	SameEventSkimmer() {}
	~SameEventSkimmer() {}

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
		
		nSig_pairs = 0;
		sig_nPos = 0;
		sig_nNeg = 0;

		

		_event = _rEvent.get();
		
		size_t nTracks = _rTracks.N();
		// if ( nTracks < 2 )
		// 	return;

		nPos = 0;
		nNeg = 0;
		nPairs = 0;
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			if ( false == _trackFilter.pass( _proxy ) )
				continue;
			if ( _proxy._track->charge() > 0 )
				nPos++;
			if ( _proxy._track->charge() < 0 )
				nNeg++;
		}



		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			if ( false == _trackFilter.pass( _proxy ) )
				continue;
			_proxy._pid = _mlp.evaluate( _proxy );

			if ( _proxy._pid > 0.6 ){
				if ( _proxy._track->charge() > 0 )
					sig_nPos++;
				if ( _proxy._track->charge() < 0 )
					sig_nNeg++;
			}


			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rTracks, _rMtdPid );
				_proxy2._pid = _mlp.evaluate( _proxy2 );
				if ( false == _trackFilter.pass( _proxy2 ) )
					continue;
				analyze_pair( _proxy, _proxy2 );
			}

		}


		// Fill Aggregate quantities
		book->fill( "nSig", nSig_pairs );
		book->fill( "sig_nPos", sig_nPos );
		book->fill( "sig_nNeg", sig_nNeg );
		book->fill( "sig_nPos_nNeg", sig_nNeg, sig_nPos );

		book->fill( "nPairs", nPairs );
		book->fill( "nPos", nPos );
		book->fill( "nNeg", nNeg );
		book->fill( "nPos_nNeg", nNeg, nPos );



	} // analyse Event
	

	virtual void fill_prefix( 	FemtoTrackProxy &_ltp, FemtoTrackProxy &_sltp, string prefix ){

		TLorentzVector lv1 = _ltp._track->lv( 0.105 );
		TLorentzVector lv2 = _sltp._track->lv( 0.105 );
		TLorentzVector lv = lv1 + lv2;

		book->fill( prefix + "mass_vs_mlp", lv.M() );

		book->fill( prefix + "mlp_vs_mlp", _ltp._pid, _sltp._pid );

		float r_sig = sqrt( pow( _ltp._pid - 1.0, 2 ) + pow( _sltp._pid - 1.0, 2 ) );
		float r_bg = sqrt( pow( _ltp._pid, 2 ) + pow( _sltp._pid, 2 ) );
		book->fill( prefix + "mass_vs_bgr", lv.M(), r_bg );
		book->fill( prefix + "mass_vs_sigr", lv.M(), r_sig );

		if ( "uls_" == prefix ){
			book->fill( "np_mass", lv.M(), nPos );
			book->fill( "nn_mass", lv.M(), nNeg );
		}
	}

	virtual void fill( 	FemtoTrackProxy &_ltp, FemtoTrackProxy &_sltp ){
		TLorentzVector lv1 = _ltp._track->lv( 0.105 );
		TLorentzVector lv2 = _sltp._track->lv( 0.105 );
		TLorentzVector lv = lv1 + lv2;

		float r_sig = sqrt( pow( _ltp._pid - 1.0, 2 ) + pow( _sltp._pid - 1.0, 2 ) );
		float r_bg = sqrt( pow( _ltp._pid, 2 ) + pow( _sltp._pid, 2 ) );

		if ( _ltp._pid > 0.6 && _sltp._pid > 0.6 ){
			book->fill( "sig_mass", lv.M() );
			nSig_pairs++;

		} else if ( _ltp._pid > 0.6 && _sltp._pid > 0.2 && _sltp._pid < 0.5 ){
			book->fill( "bg_mass", lv.M() );
		} else if ( _sltp._pid > 0.6 && _ltp._pid > 0.2 && _ltp._pid < 0.5 ){
			book->fill( "bg_mass", lv.M() );
		} else if ( sqrt( pow( _ltp._pid, 2 ) + pow( _sltp._pid, 2 ) ) < 0.1 ){
			book->fill( "purebg_mass", lv.M() );
		}

		




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
			fill_prefix( *ltp, *sltp, "lsp_" );
		}
		else if ( chargeSum == -2 ){
			fill_prefix( *ltp, *sltp, "lsn_" );
		}
		else if ( chargeSum == 0 ){
			fill_prefix( *ltp, *sltp, "uls_" );
			fill( *ltp, *sltp );
			nPairs++;
		}


	}

};

#endif