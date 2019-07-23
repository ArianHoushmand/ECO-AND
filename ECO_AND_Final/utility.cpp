#include "utility.h"
#include <stdio.h>
#include <math.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double fittoca_case2(double t1, void *parameters)
{
	struct fittoca_case2_params *p = (struct fittoca_case2_params *) parameters;

	double spd = p->spd;
	double final_time = p->final_time;
	double vmax = p->vmax;
	double amax = p->amax;

	return (final_time - (2 * (vmax - spd) / amax - t1))*vmax + spd * t1 + 0.5 * amax * pow(t1, 2) + (spd + amax * t1)*(2 * (vmax - spd) / amax - 2 * t1) + pow(amax, 2) / 6.0 * pow(2 * (vmax - spd) / amax - 2 * t1, 3) / (spd - vmax + (2 * (vmax - spd) / amax - t1) *amax) - 1;
}

double fittoca_case3(double t1, void *parameters)
{
	struct fittoca_case3_params *p = (struct fittoca_case3_params *) parameters;

	double spd = p->spd;
	double final_time = p->final_time;
	double amax = p->amax;

	return spd * t1 + 0.5 * amax * pow(t1, 2) + (spd + amax * t1) * (final_time - t1) + amax * pow(final_time - t1, 2) / 3.0 - 1;
}

double fittocd_case1(double t1, void *parameters)
{
	struct fittocd_case1_params *p = (struct fittocd_case1_params *) parameters;

	double spd = p->spd;
	double final_time = p->final_time;
	double vmin = p->vmin;
	double amin = p->amin;

	return (final_time - (2 * (vmin - spd) / amin - t1))*vmin + spd * t1 + 0.5*amin* pow(t1, 2) + (spd + amin * t1)*(2 * (vmin - spd) / amin - 2 * t1) + pow(amin, 2) / 6.0 * pow(2 * (vmin - spd) / amin - 2 * t1, 3) / (spd - vmin + (2 * (vmin - spd) / amin - t1)*amin) - 1;
}

double fittocd_case2(double t1, void *parameters)
{
	struct fittocd_case2_params *p = (struct fittocd_case2_params *) parameters;

	double spd = p->spd;
	double final_time = p->final_time;
	double amin = p->amin;

	return spd * t1 + 0.5*amin*pow(t1, 2) + (spd + amin * t1)*(final_time - t1) + amin / 3.0*pow(final_time - t1, 2) - 1;
}

double fttoc_solve(double tf, void *parameters)
{
	struct fttoc_params *p = (struct fttoc_params *) parameters;

	double spd = p->spd;
	double rho_t = p->rho_t;
	double rho_u = p->rho_u;

	return spd * tf + rho_t * pow(tf, 3) / (3 * rho_u * (spd + sqrt(pow(spd, 2) + rho_t * pow(tf, 2) / rho_u))) - 1;
}

 double haversine(double lat1, double lon1, double lat2, double lon2)
{
	// distance between latitudes
	// and longitudes
	double dLat = (lat2 - lat1) *
		M_PI / 180.0;
	double dLon = (lon2 - lon1) *
		M_PI / 180.0;

	// convert to radians
	lat1 = (lat1)* M_PI / 180.0;
	lat2 = (lat2)* M_PI / 180.0;

	// apply formulae
	double a = pow(sin(dLat / 2), 2) +
		pow(sin(dLon / 2), 2) *
		cos(lat1) * cos(lat2);
	double rad = 6371;
	double c = 2 * asin(sqrt(a));
	return rad * c * 1000; // distance in meters
}
