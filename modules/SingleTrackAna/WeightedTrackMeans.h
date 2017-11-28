#ifndef WEIGHTED_TRACK_MEANS
#define WEIGHTED_TRACK_MEANS


#include "HistoAnalyzer.h"

using namespace jdb;


class WeightedTrackMeans : public HistoAnalyzer
{
public:
	WeightedTrackMeans() {}
	~WeightedTrackMeans() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();
	}

	double compute_weighted_value( TH1 * functional, TH1 * source, double &error ){
		double v = 0, w =0, ve = 0, we = 0;
		
		for ( int i = 1; i <= source->GetXaxis()->GetNbins(); i++ ){
			float sv = source->GetBinContent(i);
			float se = source->GetBinError(i);
			float fv = functional->GetBinContent(i);
			float fe = functional->GetBinError(i);

			if ( 0 == fv ) fe = 0;
			LOG_F( INFO, "sv=%f+-%f, fv=%f+/-%f, v=%f, w=%f", sv,se, fv,fe, v, w );
			v += fv * sv;
			if ( fv != 0 && sv != 0 )
				ve += (fv * sv) * sqrt( pow(se/sv,2) + pow( fe/fv,2 ) );
			w += sv;
			we += se;
		}

		if ( w == 0 )
			return  -1;
		error = (v/w) * sqrt( pow(ve/v,2)+pow(we/w,2) );
		return v / w;
	}

	virtual void make(){

		book->cd();
		TH2 * mb_grm2 = (TH2*)(getH2D("nPos_vs_grefmult", "mb")->Clone("mb_nPos_vs_grefmult"));
		TH2 * dimuon_grm2 = (TH2*)(getH2D("nPos_vs_grefmult", "dimuon")->Clone("dimuon_nPos_vs_grefmult"));

		TH2 * mb_vz2 = (TH2*)(getH2D("nPos_vs_vz", "mb")->Clone("mb_nPos_vs_vz"));
		TH2 * dimuon_vz2 = (TH2*)(getH2D("nPos_vs_vz", "dimuon")->Clone("dimuon_nPos_vs_vz"));


		TH1 * mb_grm = mb_grm2->ProjectionX("mb_grefmult");
		TH1 * dimuon_grm = dimuon_grm2->ProjectionX("dimuon_grefmult");

		TH1 * mb_vz = mb_vz2->ProjectionX("mb_vz");
		TH1 * dimuon_vz = dimuon_vz2->ProjectionX("dimuon_vz");

		mb_grm->SetLineColor(kRed);
		mb_vz->SetLineColor(kRed);


		double ep = 0;
		double mp = compute_weighted_value( getH1D("mean_nPos_vs_grefmult", "mb"), dimuon_grm, ep );
		LOG_F( INFO, "m = %f +/- %f", mp, ep );

		double en = 0;
		double mn = compute_weighted_value( getH1D("mean_nNeg_vs_grefmult", "mb"), dimuon_grm, en );
		LOG_F( INFO, "m = %f +/- %f", mn, en );

		LOG_F( INFO, "nbg = %f", mn*mp * (9.71860690858818210e-01*9.69813530022919257e-01) * 3.54e11 );

	}
	
};



#endif

