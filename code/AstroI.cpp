
#include <AstroI.h>
#include <sstream>
#include <iostream>
#include <vector>

#include <cstring>

#include <GenmapParams.h>
#include "MathUtils.h"

#include "DBAstro.h"
#include "Astro.h"
#include <Freeze/Freeze.h>

using namespace std;

::Astro::Matrix
Astro::AstroInterfaceI::calculateMapKey(const ::Astro::SeqEvtKey& keys,
                                        const Ice::Current& current)
{
	CtsGenParams params;
	params.Load(0,NULL);

	vector<double> raVector, decVector;

	int numcol = 0;
	long nrows = 0;
	int status = 0;
	double l = 0, b = 0;
	double x = 0, y = 0;
	int i = 0, ii = 0;
	double the = 0;
	long mxdim= params.mxdim; // dimension (in pixels) of the map
	//unsigned short A[mxdim][mxdim]; //counts map
	Astro::Matrix A(mxdim, SimpleSeq(mxdim));
	int selectedEvents = 0;

	for (int var = 0; var < keys.size(); ++var) {
		raVector.push_back(keys[var].ra);
		decVector.push_back(keys[var].dec);
	}

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			A[i][ii] = 0;
		}
	}

	double baa = params.ba * DEG2RAD;
	double laa = params.la * DEG2RAD;

	int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
	long naxis    =   2;  /* 2-dimensional image                            */
	long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */

	const double obtlimit = 104407200.;
	if (params.tmin<obtlimit)
		status = 1005;
	else
		//status = addfile(evtFits, params);
	std::cout << "AG_ctsmapgen0...................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;

	nrows = raVector.size();

	cout << nrows << endl;
	//double ra, dec;
	switch (params.projection) {
	case ARC:
		for (long k = 0; k<nrows; ++k) {
			//sostituire tutto questo con l'accesso al singolo elemento k del vector<double>
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);

			Euler(raVector[k], decVector[k], &l, &b, 1);
				//Euler(ra, dec, &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
			y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

			i = (int)floor(((-x+(params.mdim/2.))/params.mres));
			ii = (int)floor(((y+(params.mdim/2.))/params.mres));

			if (params.inmap(i,ii)) {
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
		break;
	case AIT:
		for (long k = 0; k<nrows; ++k) {
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			//Euler(ra, dec, &l, &b, 1);

		Euler(raVector[k], decVector[k], &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			l=l-laa;

			if ( l < M_PI  )
				l=-l;
			else
				l=2*M_PI -l;

			x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;
			y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

			i=(int)floor(((x+(params.mdim/2.))/params.mres));
			ii=(int)floor(((y+(params.mdim/2.))/params.mres));

			if (params.inmap(i,ii)) {
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
		break;
	}

	return A;
}

::Astro::Matrix
Astro::AstroInterfaceI::calculateMapVector(const ::Astro::Ra& raVector,
                                           const ::Astro::Dec& decVector,
                                           const Ice::Current& current)
{
	CtsGenParams params;
	params.Load(0,NULL);
	int numcol = 0;
	long nrows = 0;
	int status = 0;
	double l = 0, b = 0;
	double x = 0, y = 0;
	int i = 0, ii = 0;
	double the = 0;
	long mxdim= params.mxdim; // dimension (in pixels) of the map
	//unsigned short A[mxdim][mxdim]; //counts map
	Astro::Matrix A(mxdim, SimpleSeq(mxdim));
	int selectedEvents = 0;

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			A[i][ii] = 0;
		}
	}

	double baa = params.ba * DEG2RAD;
	double laa = params.la * DEG2RAD;

	int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
	long naxis    =   2;  /* 2-dimensional image                            */
	long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */

	const double obtlimit = 104407200.;
	if (params.tmin<obtlimit)
		status = 1005;
	else
		//status = addfile(evtFits, params);
	std::cout << "AG_ctsmapgen0...................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;

	nrows = raVector.size();

	cout << nrows << endl;
	//double ra, dec;
	switch (params.projection) {
	case ARC:
		for (long k = 0; k<nrows; ++k) {
			//sostituire tutto questo con l'accesso al singolo elemento k del vector<double>
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);

			Euler(raVector[k], decVector[k], &l, &b, 1);
				//Euler(ra, dec, &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
			y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

			i = (int)floor(((-x+(params.mdim/2.))/params.mres));
			ii = (int)floor(((y+(params.mdim/2.))/params.mres));

			if (params.inmap(i,ii)) {
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
		break;
	case AIT:
		for (long k = 0; k<nrows; ++k) {
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			//Euler(ra, dec, &l, &b, 1);

		Euler(raVector[k], decVector[k], &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			l=l-laa;

			if ( l < M_PI  )
				l=-l;
			else
				l=2*M_PI -l;

			x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;
			y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

			i=(int)floor(((x+(params.mdim/2.))/params.mres));
			ii=(int)floor(((y+(params.mdim/2.))/params.mres));

			if (params.inmap(i,ii)) {
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
		break;
	}

    return A;
}

void
Astro::AstroInterfaceI::shutdown(const Ice::Current& current)
{
}
