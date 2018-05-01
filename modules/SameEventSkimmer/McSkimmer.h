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

	string cloop = "mc";


public:
	virtual const char* classname() const {return "McSkimmer";}
	McSkimmer() {}
	~McSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		loguru::add_file( TString::Format("McSkimmer_%d.log", config.getInt("jobIndex", -1) ).Data(), loguru::Truncate, loguru::Verbosity_INFO );

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMcTracks.setup( chain, "McTracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		_trackFilter.load( config, nodePath + ".TrackFilter" );
		// _mlp.load( config, nodePath + ".MuonMLPFilter" );

		sMin = 0.8;
	}

	bool keepMcTrack( FemtoTrackProxy &_tp ){
		if ( "mc" == cloop ) book->fill( "mc_pass", "All" );
		if ( nullptr == _tp._mcTrack )
			return false;
		if ( "mc" == cloop ) book->fill("mc_pass", "nnull");
		if ( _tp._mcTrack->mParentIndex >= 0 )
			return false;
		if ( "mc" == cloop ) book->fill( "mc_pass", "primary" );

		if ( fabs(_tp._mcTrack->mEta) > 0.5 )
			return false;
		if ( "mc" == cloop ) book->fill( "mc_pass", "#eta" );
		
		return true;
	} 


	bool keepTrack( FemtoTrackProxy &_tp ){
		book->fill( cloop + "_pass", "All" );
		if ( nullptr == _tp._track )
			return false;
		book->fill( cloop + "_pass", "nnull" );
		
		if ( nullptr == _tp._mcTrack )
			return false;
		book->fill( cloop + "_pass", "nomc" );

		if ( fabs(_tp._track->mNHitsFit) < 20 )
			return false;
		book->fill( cloop + "_pass", "nhfit" );

		if (_tp._track->mNHitsDedx < 15 )
			return false;
		book->fill( cloop + "_pass", "nhdedx" );

		if ( fabs(_tp._track->mEta) > 0.5 )
			return false;
		book->fill( cloop + "_pass", "#eta" );

		return true;
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
		cloop = "mc";
		LOG_F( 1, "nMcTracks=%zu", nMcTracks );
		for (size_t i = 0; i < nMcTracks; i++ ){
			_proxy.assemble( i, _rMcTracks, _rTracks, _rMtdPid );
			
			if ( _proxy._mcTrack->mPt < 0.1 ) continue;
			if ( false == keepMcTrack( _proxy ) ) continue;

			analyze_track( _proxy, "mc_" );

			for (size_t j = i; j < nMcTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rMcTracks, _rTracks, _rMtdPid );

				if ( _proxy2._mcTrack->mPt < 0.1 ) continue;
				if ( false == keepMcTrack( _proxy2 ) ) continue;
				_proxy2._pid = 0;

				if ( _proxy._mcTrack->mCharge == 1 && _proxy2._mcTrack->mCharge == -1 )
					analyze_pair( _proxy, _proxy2, "mc_" );
				else if ( _proxy._mcTrack->mCharge == -1 && _proxy2._mcTrack->mCharge == 1 )
					analyze_pair( _proxy2, _proxy, "mc_" );
				
			}
		}

		////////////////////////////////////////////////////////
		/// TPC
		/// 
		cloop = "tpc";
		LOG_F( 1, "TPC nTracks=%zu", nTracks );
		for (size_t i = 0; i < nTracks; i++ ){
			
			_proxy.assemble( i, _rTracks, _rMtdPid );
			if ( _proxy._track->mPt < 0.1 ) continue;
			if ( _proxy._track->mMcIndex < 0 ) continue;

			_proxy._mcTrack = _rMcTracks.get( _proxy._track->mMcIndex );

			if ( false == keepMcTrack( _proxy ) ) continue;
			if ( false == keepTrack(_proxy) ) continue;

			analyze_track( _proxy, "tpc_" );

			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				
				_proxy2.assemble( j, _rTracks, _rMtdPid );
				
				if ( _proxy2._track->mPt < 0.1 ) continue;
				if ( _proxy2._track->mMcIndex < 0 ) continue;

				_proxy2._mcTrack = _rMcTracks.get( _proxy2._track->mMcIndex );

				if ( false == keepMcTrack( _proxy2 ) ) continue;
				if ( false == keepTrack( _proxy2 ) ) continue;

				if ( _proxy._track->charge() == 1 && _proxy2._track->charge() == -1 )
					analyze_pair( _proxy, _proxy2, "tpc_" );
				else if ( _proxy._track->charge() == -1 && _proxy2._track->charge() == 1 )
					analyze_pair( _proxy2, _proxy, "tpc_" );
			}
		}


		////////////////////////////////////////////////////////
		/// MTD
		/// 
		cloop = "mtd";
		LOG_F( 1, "MTD nTracks=%zu", nTracks );
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rMtdPid );
			
			if ( nullptr == _proxy._mtdPid  ) continue;
			if ( _proxy._track->mPt < 0.01 ) continue;

			if ( false == _trackFilter.pass( _proxy ) )
				continue;
			if ( _proxy._track->mMcIndex < 0 ) continue;
			_proxy._mcTrack = _rMcTracks.get( _proxy._track->mMcIndex );
			if ( false == keepMcTrack( _proxy ) ) continue;
			if (false == keepTrack(_proxy)) continue;

			if ( _proxy._track->charge() > 0 )
				nPos++;
			if ( _proxy._track->charge() < 0 )
				nNeg++;

			analyze_track( _proxy, "mtd_" );

			// _proxy._pid = _mlp.evaluate( _proxy );

			for (size_t j = i; j < nTracks; j++ ){
				if ( i == j ) continue;
				_proxy2.assemble( j, _rTracks, _rMtdPid );

				if ( nullptr == _proxy2._mtdPid  ) continue;
				if ( _proxy2._track->mPt < 0.01 ) continue;

				if ( false == _trackFilter.pass( _proxy2 ) )
					continue;

				if ( _proxy2._track->mMcIndex < 0 ) continue;
				_proxy2._mcTrack = _rMcTracks.get( _proxy2._track->mMcIndex );
				if ( false == keepMcTrack( _proxy2 ) ) continue;
				if ( false == keepTrack(_proxy2)) continue;

				// _proxy2._pid = _mlp.evaluate( _proxy2 );
				if ( _proxy._track->charge() == 1 && _proxy2._track->charge() == -1 )
					analyze_pair( _proxy, _proxy2, "mtd_" );
				else if ( _proxy._track->charge() == -1 && _proxy2._track->charge() == 1 )
					analyze_pair( _proxy2, _proxy, "mtd_" );
			}
		}


		// Fill Aggregate quantities
		book->fill( "nPairs", nPairs );
		book->fill( "nPos", nPos );
		book->fill( "nNeg", nNeg );
		book->fill( "nPos_nNeg", nNeg, nPos );

	} // analyse Event

	string chargeString( int charge ){
		if ( charge > 0 )
			return "p";
		if ( charge < 0 )
			return "m";
		return "";
	}

	string chargeStringLong( int charge ){
		if ( charge > 0 )
			return "pos";
		if ( charge < 0 )
			return "neg";
		return "";
	}


	virtual void analyze_track( FemtoTrackProxy &_tp, string _prefix = "mc_" ){
		book->fill( _prefix + chargeString(_tp._mcTrack->mCharge) + "_mPt", _tp._mcTrack->mPt );
		book->fill( _prefix + chargeString(_tp._mcTrack->mCharge) + "_mEta", _tp._mcTrack->mEta );
		book->fill( _prefix + chargeString(_tp._mcTrack->mCharge) + "_mPhi", _tp._mcTrack->mPhi );


		if ( "mc_" == _prefix )
			book->fill( chargeStringLong( _tp._mcTrack->mCharge ) + "_mc", _tp._mcTrack->mPhi, _tp._mcTrack->mEta, _tp._mcTrack->mPt );
		if ( "tpc_" == _prefix )
			book->fill( chargeStringLong( _tp._mcTrack->mCharge ) + "_rc", _tp._mcTrack->mPhi, _tp._mcTrack->mEta, _tp._mcTrack->mPt );
		if ( "mtd_" == _prefix )
			book->fill( chargeStringLong( _tp._mcTrack->mCharge ) + "_mtd", _tp._mcTrack->mPhi, _tp._mcTrack->mEta, _tp._mcTrack->mPt );
		
	}

	virtual void analyze_pair( FemtoTrackProxy &_tp1, FemtoTrackProxy &_tp2, string _prefix = "mtd_" ){
		
		int chargeSum = -999; 

		const float muon_mass = 0.1056583745; // muon mass in GeV/c^2
		TLorentzVector lv1, lv2, lv; 
		int c1 = 0, c2 = 0;

		if ( nullptr != _tp1._track && nullptr != _tp2._track && "mc_" != _prefix ){
			if ( _tp1._track->mId == _tp2._track->mId) return;
		}

		if ( _tp1._mcTrack->mId == _tp2._mcTrack->mId) return;
		chargeSum = _tp1._mcTrack->mCharge + _tp2._mcTrack->mCharge;
		c1 = _tp1._mcTrack->mCharge;
		c2 = _tp2._mcTrack->mCharge;

		lv1 = _tp1._mcTrack->lv( muon_mass );
		lv2 = _tp2._mcTrack->lv( muon_mass );
		

		lv = lv1 + lv2;

		LOG_F( 1, "prefix=%s, Y=%f, chargeSum=%d", _prefix.c_str(), lv.Rapidity(), chargeSum );
		if ( fabs( lv.Rapidity() ) > 0.5 ) return;

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