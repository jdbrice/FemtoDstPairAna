#ifndef MC_SKIMMER_H
#define MC_SKIMMER_H

#include "TreeAnalyzer.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "vendor/loguru.h"

#include "Filters/MtdTrackFilter.h"
#include "Filters/MuonMLPFilter.h"


class McSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;
	FemtoTrackProxy _proxy;
	FemtoTrackProxy _proxy2;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMcTrack> _rMcTracks;
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
	virtual const char* classname() const {return "McSkimmer";}
	McSkimmer() {}
	~McSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();


		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMcTracks.setup( chain, "McTracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		_trackFilter.load( config, nodePath + ".TrackFilter" );
		_mlp.load( config, nodePath + ".MuonMLPFilter" );

		sMin = 0.8;
	}

protected:

	virtual void analyzeEvent(){

		_event = _rEvent.get();
		
		size_t nTracks = _rTracks.N();
		size_t nMcTracks = _rMcTracks.N();

		nPos = 0;
		nNeg = 0;
		nPairs = 0;

		////////////////////////////////////////////////////////
		/// MC
		/// 
		for (size_t i = 0; i < nMcTracks; i++ ){
			_proxy.assemble( i, _rMcTracks, _rTracks, _rMtdPid );
			
			if ( _proxy._mcTrack->mPt < 0.1 ) continue;

			for (size_t j = i; j < nMcTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rMcTracks, _rTracks, _rMtdPid );

				if ( _proxy2._mcTrack->mPt < 0.1 ) continue;
				_proxy2._pid = 0;
				analyze_pair( _proxy, _proxy2, "mc_" );
			}
		}

		////////////////////////////////////////////////////////
		/// TPC
		/// 
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			if ( _proxy._track->mPt < 0.1 ) continue;
			if ( _proxy._track->mMcIndex < 0 ) continue;
			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rTracks, _rMtdPid );
				if ( _proxy2._track->mPt < 0.1 ) continue;
				if ( _proxy2._track->mMcIndex < 0 ) continue;
				analyze_pair( _proxy, _proxy2, "tpc_" );
			}
		}


		////////////////////////////////////////////////////////
		/// MTD
		/// 
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			
			if ( nullptr == _proxy._mtdPid  ) continue;
			if ( _proxy._track->mPt < 0.01 ) continue;

			if ( false == _trackFilter.pass( _proxy ) )
				continue;
			if ( _proxy._track->mMcIndex < 0 ) continue;

			
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

				if ( _proxy2._track->mMcIndex < 0 ) continue;
				
				_proxy2._pid = _mlp.evaluate( _proxy2 );
				analyze_pair( _proxy, _proxy2, "mtd_" );
			}
		}


		// Fill Aggregate quantities
		book->fill( "nPairs", nPairs );
		book->fill( "nPos", nPos );
		book->fill( "nNeg", nNeg );
		book->fill( "nPos_nNeg", nNeg, nPos );

	} // analyse Event

	virtual void analyze_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2, string _prefix = "mtd_" ){
		
		int chargeSum = -999; 

		const float muon_mass = 0.1056583745; // muon mass in GeV/c^2
		TLorentzVector lv1, lv2, lv; 

		if ( nullptr != _tp1._track && nullptr != _tp2._track && "mc_" != _prefix ){
			if ( _tp1._track->mId == _tp2._track->mId) return;
			chargeSum = _tp1._track->charge() + _tp2._track->charge();
			
			lv1 = _tp1._track->lv( muon_mass );
			lv2 = _tp2._track->lv( muon_mass );
		} else {
			if ( _tp1._mcTrack->mId == _tp2._mcTrack->mId) return;
			chargeSum = _tp1._mcTrack->mCharge + _tp2._mcTrack->mCharge;
			
			lv1 = _tp1._mcTrack->lv( muon_mass );
			lv2 = _tp2._mcTrack->lv( muon_mass );
		}

		lv = lv1 + lv2;

		if ( 0 == chargeSum ){
			// pair mass vs. d1 pt vs. d2 pt
			book->fill( _prefix + "mass_pt1_pt2", lv.M(), lv1.Pt(), lv2.Pt() );
			
			// pair pt vs d1 pt vs. d2 pt
			book->fill( _prefix + "pt_pt1_pt2", lv.Pt(), lv1.Pt(), lv2.Pt() );

			// pair pt vs. pair mass
			book->fill( _prefix + "pt_vs_mass", lv.M(), lv.Pt() );
			nPairs++;
		}


	}

};

#endif