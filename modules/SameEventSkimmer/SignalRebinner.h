#ifndef SIGNAL_REBIN_H
#define SIGNAL_REBIN_H

#include "HistoAnalyzer.h"

using namespace jdb;


#include "vendor/loguru.h"

class SignalRebinner : public HistoAnalyzer
{
public:
	SignalRebinner() {}
	~SignalRebinner() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();
		LOG_F( INFO, "test" );
	}

	TH1* do_rebin( string in, string out ){
		TH1 * hsr = getH1D( in );

		HistoBins signal_mass_bins;
		signal_mass_bins.load( config, "bins.signalMass" );

		TH1 * hsb = hsr->Rebin( signal_mass_bins.nBins(), out.c_str(), signal_mass_bins.bins.data() );
		hsb->Scale( 1.0, "width" );
		return hsb;
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);

		book->cd();
		// book->makeAll( nodePath + ".histograms" );

		TH1 * hbgr = getH1D( "bg_mass" );

		// rough background scaler
		float pid_eff_pos_0p8 = 0.202710;
		float pid_eff_neg_0p8 = 0.219116;
		float NRawBG = 472717.116543;
		float nBG = NRawBG  * 0.665655 * pid_eff_neg_0p8 * pid_eff_pos_0p8;// (1.0/(9.71860690858818210e-01*9.69813530022919257e-01))
		
		hbgr->Scale( nBG / hbgr->Integral() );
		// TH1 * hbg = (TH1*)hbgr->Clone( "bg_mass" );


		// getH1D("uls_mass")->Clone( "uls_mass" );
		// TH1 * hls = (TH1*)getH1D("ls_mass")->Clone( "ls_mass" );

		TH1 * huls = do_rebin( "uls_mass", "uls_mass" );
		TH1 * hls = do_rebin( "ls_mass", "ls_mass" );
		TH1 * hbg = do_rebin( "bg_mass", "bg_mass" );
		do_rebin( "purebg_mass", "purebg_mass" );

		hls->SetLineColor( kRed );
		hbg->SetLineColor( kBlack );

		TH1 * hs1 = (TH1*)huls->Clone( "sig_1" );
		TH1 * hs2 = (TH1*)huls->Clone( "sig_2" );


		hs1->Add( hbg, -1 );
		hs2->Add( hls, -1 );
		


	}
	
};



#endif