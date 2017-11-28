#ifndef SINGLE_TRACK_ANA_H
#define SINGLE_TRACK_ANA_H

#include "HistoAnalyzer.h"

using namespace jdb;


#include "vendor/loguru.h"

class SingleTrackAna : public HistoAnalyzer
{
public:
	SingleTrackAna() {}
	~SingleTrackAna() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();
		LOG_F( INFO, "test" );
	}

	double weighted_mean( TH1 * h, double &error ){
		if ( nullptr == h )
			return -999;
		double v = 0, w = 0, ve = 0, we = 0;
		for ( int i = 1; i <= h->GetXaxis()->GetNbins(); i++  ){
			v += (i-1) * h->GetBinContent(i);
			ve += (i-1) * h->GetBinError(i);
			w += h->GetBinContent(i);
			we += h->GetBinError(i);
		}
		if ( w == 0 )
			return 0;

		error = (v/w) * sqrt( pow(ve/v, 2) + pow(we/w, 2) );
		LOG_F( 2, "v=%f, w=%f, ve=%f, we=%f, error=%f", v, w, ve, we, error );
		return v / w;
	}

	void make_profile( string name, int nRebin = 1 ){
		LOG_SCOPE_FUNCTION(INFO);
		LOG_F( INFO, "Making profile for %s", name.c_str() );
		TH2 * h2 = (TH2*)( getH2D( name )->Clone(name.c_str()) );
		TH1 * hmean = book->get( "mean_" + name );
		hmean->Rebin( nRebin );
		LOG_F( INFO, "h2=%p, h=%p", h2, hmean );
		int iBin = 1;
		for ( int i = 1; i < h2->GetXaxis()->GetNbins(); i+=nRebin ){
			TH1 * h = h2->ProjectionY("temp", i, i+nRebin-1);
			double we = 0;
			double wm = weighted_mean( h, we );
			hmean->SetBinContent( iBin, wm );
			hmean->SetBinError( iBin, we );
			LOG_F( 2, "error=%f", we );
			
			delete h;
			iBin++;
		}
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);

		book->cd();
		book->makeAll( nodePath + ".histograms" );

		vector<string> names = {
			"nPos_vs_vz", "nNeg_vs_vz", "nPos_vs_grefmult", "nNeg_vs_grefmult"
		};
		vector<int> nbins = {70, 70, 1, 1};
		
		for ( int i = 0; i < names.size(); i++  ){
			string n = names[i];
			int bs = nbins[i];
			make_profile( n, bs );
		}



	}
	
};



#endif