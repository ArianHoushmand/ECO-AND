#pragma once
struct fittoca_case2_params
{
	double spd, final_time, vmax, amax;
};

struct fittoca_case3_params
{
	double spd, final_time, amax;
};

struct fittocd_case1_params
{
	double spd, final_time, vmin, amin;
};

struct fittocd_case2_params
{
	double spd, final_time, amin;
};

struct fttoc_params
{
	double spd, rho_t, rho_u;
};

double fittoca_case2(double t1, void *parameters);
double fittoca_case3(double t1, void *parameters);
double fittocd_case1(double t1, void *parameters);
double fittocd_case2(double t1, void *parameters);
double fttoc_solve(double tf, void *parameters);
double haversine(double lat1, double lon1, double lat2, double lon2);
