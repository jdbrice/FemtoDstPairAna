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
#include "Filters/MuonMLPFilter.h"


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

	float sMin;



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

		sMin = 0.8;
	}

protected:

	virtual void analyzeEvent(){

		_event = _rEvent.get();
		
		size_t nTracks = _rTracks.N();

		nPos = 0;
		nNeg = 0;
		nPairs = 0;


		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			
			if ( nullptr == _proxy._mtdPid  ) continue;
			if ( _proxy._track->mPt < 0.01 ) continue;

			if ( false == _trackFilter.pass( _proxy ) )
				continue;

			
			if ( _proxy._track->charge() > 0 )
				nPos++;
			if ( _proxy._track->charge() < 0 )
				nNeg++;


			_proxy._pid = _mlp.evaluate( _proxy );

			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rTracks, _rMtdPid );

				if ( nullptr == _proxy2._mtdPid  ) continue;
				if ( _proxy2._track->mPt < 0.01 ) continue;

				if ( false == _trackFilter.pass( _proxy2 ) )
					continue;
				
				_proxy2._pid = _mlp.evaluate( _proxy2 );
				analyze_pair( _proxy, _proxy2 );
			}
		}


		// Fill Aggregate quantities
		book->fill( "nPairs", nPairs );
		book->fill( "nPos", nPos );
		book->fill( "nNeg", nNeg );
		book->fill( "nPos_nNeg", nNeg, nPos );

	} // analyse Event
	

	bool signal_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		if ( _tp1._pid > sMin && _tp2._pid > sMin )
			return true;
		return false;
	}

	bool mixed_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		if ( ( _tp1._pid > 0.6 && _tp2._pid < 0.4 ) || ( _tp2._pid > 0.6 && _tp1._pid < 0.4 ) )
			return true;
		return false;
	}

	virtual void fill_pid( TLorentzVector &lv, FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2, string prefix ){

		book->fill( prefix + "mass_vs_mlp", lv.M() );
		book->fill( prefix + "mlp_vs_mlp", _tp1._pid, _tp2._pid );

		float r_sig = sqrt( pow( _tp1._pid - 1.0, 2 ) + pow( _tp2._pid - 1.0, 2 ) );
		float r_bg = sqrt( pow( _tp1._pid, 2 ) + pow( _tp2._pid, 2 ) );

		book->fill( prefix + "mass_vs_bgr", lv.M(), r_bg );
		book->fill( prefix + "mass_vs_sigr", lv.M(), r_sig );

		// if ( "lsp_" == prefix || "lsn_" == prefix ){
		// 	if ( signal_pair( _tp1, _tp2 ) ){
		// 		book->fill( "ls_mass", lv.M() );
		// 		book->fill( "ls_pt_vs_mass", lv.M(), lv.Pt() );
		// 	}
		// }
	}

	virtual void fill( 	FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		TLorentzVector lv1 = _tp1._track->lv( 0.105 );
		TLorentzVector lv2 = _tp2._track->lv( 0.105 );
		TLorentzVector lv = lv1 + lv2;

		float r_sig = sqrt( pow( _tp1._pid - 1.0, 2 ) + pow( _tp2._pid - 1.0, 2 ) );
		float r_bg = sqrt( pow( _tp1._pid, 2 ) + pow( _tp2._pid, 2 ) );

		float dca_r = sqrt( pow( _tp1._track->gDCA(), 2 ) + pow(_tp2._track->gDCA(), 2) );

		if ( signal_pair( _tp1, _tp2 ) ){
			book->fill( "sig_mass", lv.M() );
			book->fill( "sig_pt_vs_mass", lv.M(), lv.Pt() );
			nSig_pairs++;
		} 
		
		
	}

	virtual void analyze_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2 ){
		if ( _tp1._track->mId == _tp2._track->mId) return;

		int chargeSum = _tp1._track->charge() + _tp2._track->charge();

		const float muon_mass = 0.1056583745; // GeV/c^2
		TLorentzVector lv = _tp1._track->lv( muon_mass ) + _tp2._track->lv( muon_mass );

		if ( lv.Pt() < 2.0 )
			return;


		if ( chargeSum == 2 ){
			// fill_pid( lv, _tp1, _tp2, "lsp_" );
			// fill_pid( lv, _tp1, _tp2, "ls_" );
			if ( signal_pair( _tp1, _tp2 ) ){
				book->fill( "ls_mass", lv.M() );
				book->fill( "ls_pt_vs_mass", lv.M(), lv.Pt() );
			}
		}
		else if ( chargeSum == -2 ){
			// fill_pid( lv, _tp1, _tp2, "lsn_" );
			// fill_pid( lv, _tp1, _tp2, "ls_" );
			if ( signal_pair( _tp1, _tp2 ) ){
				book->fill( "ls_mass", lv.M() );
				book->fill( "ls_pt_vs_mass", lv.M(), lv.Pt() );
			}
		}
		else if ( chargeSum == 0 ){
			// fill_pid( lv, _tp1, _tp2, "uls_" );
			
			book->fill( "ruls_mass", lv.M() );
			book->fill( "ruls_pt_vs_mass", lv.M(), lv.Pt() );

			if ( signal_pair( _tp1, _tp2 ) ){
				book->fill( "uls_mass", lv.M() );
				book->fill( "uls_pt_vs_mass", lv.M(), lv.Pt() );
			} else if ( _tp1._track->gDCA() > 1.0 || _tp2._track->gDCA() > 1.0 ){
				if ( mixed_pair( _tp1, _tp2 ) ){
					book->fill( "bg_mass", lv.M() );
					book->fill( "bg_pt_vs_mass", lv.M(), lv.Pt() );
				} else {
					book->fill( "purebg_mass", lv.M() );
					book->fill( "purebg_pt_vs_mass", lv.M(), lv.Pt() );
				}
			}

			nPairs++;
		}


	}

};

#endif