/*
* @Author: Daniel
* @Date:   2017-05-05 14:22:49
* @Last Modified by:   Daniel
* @Last Modified time: 2017-05-26 10:16:12
*/
#ifndef MTD_TRACK_FILTER_H
#define MTD_TRACK_FILTER_H

#include "XmlConfig.h"
#include "XmlRange.h"

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// Filters
#include "Filters/TrackFilter.h"

// STL
#include <map>

#include "vendor/loguru.h"

class MtdTrackFilter : public TrackFilter{
public:
	MtdTrackFilter() {}
	MtdTrackFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~MtdTrackFilter() {}

	bool pass( FemtoTrackProxy &_proxy ){
		bool res = TrackFilter::pass( _proxy );
		if ( false == res)
			return false;

		if ( nullptr == _proxy._mtdPid )
			return false;

		if ( range.count( "mtdMatchFlag" ) && !range[ "mtdMatchFlag" ].inInclusiveRange( _proxy._mtdPid->mMatchFlag ) )
			return false;

		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif