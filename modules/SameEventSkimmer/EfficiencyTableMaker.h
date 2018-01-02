#ifndef EFFICIENCY_TABLE_MAKER_H
#define EFFICIENCY_TABLE_MAKER_H


#include "HistoAnalyzer.h"
using namespace jdb;

#include "vendor/loguru.h"

class EfficiencyTableMaker : public HistoAnalyzer {

public:

	virtual void initialize(){
		LOG_F( INFO, "" );
	}

	virtual void make(){
		LOG_F(INFO, "");
		book->cd();

		TH1 * hmc = getH1D( "mc_mass_pt1_pt2" );
		TH1 * hmtd = getH1D( "mtd_mass_pt1_pt2" );

		TH3 * heff = (TH3*)hmtd->Clone( "mtd_eff_mass_pt1_pt2" );

		heff->Divide( hmc );


	}

};


#endif