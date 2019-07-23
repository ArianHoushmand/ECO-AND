#include "ecoand.h"
#include "utility.h"

ecoand::ecoand()
{
}


ecoand::~ecoand()
{
}

#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */

void ecoand::init(long id, long lat, long lon, double tm, double pos, double spd, double acc)
{
	veh_id = id;
	current_lat = lat;
	current_lon = lon;
	cur_time = tm;
	cur_pos = pos;
	cur_spd = spd;
	cur_acc = acc;
	in_ecoand = 0;
	lead_id = -1;

	lead_acc = 0;
	lead_dist = 0;
	lead_spd = 0;

	rho_u = 0.0100;
	rho_t = 0.2500;

	mode = 0;
}

void ecoand::init(long id, long lat, long lon, double tm, double pos, double spd, double acc, long f_id, double f_acc, double f_spd_diff, double f_dist)
{
	veh_id = id;
	current_lat = lat;
	current_lon = lon;
	cur_time = tm;
	cur_pos = pos;
	cur_spd = spd;
	cur_acc = acc;


	rho_u = 0.0100;
	rho_t = 0.2500;

	lead_id = f_id;
	lead_acc = f_acc;
	lead_dist = f_dist;
	lead_spd = spd - f_spd_diff;

	in_ecoand = 0;
	mode = 0;
}


void ecoand::set_limits(double max_spd, double min_spd, double max_acc, double min_acc)
{
	vmax = max_spd;
	vmin = min_spd;
	amax = max_acc;
	amin = min_acc;
}

void ecoand::update_status(long lk, long lane, double tm, double spd, double acc, double pos)
{
	cur_lk = lk;
	cur_ln = lane;
	cur_time = tm;
	cur_spd = spd;
	cur_acc = acc;
	cur_pos = pos;
}

void ecoand::update_status(long lk, long lane, double tm, double spd, double acc, double pos, long f_id, double f_acc, double f_spd_diff, double f_dist)
{
	cur_lk = lk;
	cur_ln = lane;
	cur_time = tm;
	cur_spd = spd;
	cur_acc = acc;
	cur_pos = pos;

	lead_id = f_id;
	lead_acc = f_acc;
	lead_dist = f_dist;
	lead_spd = spd - f_spd_diff;
}

vector<double> ecoand::fittoca(double spd, double dist_to_sig, double final_time)
{
	map<double, vector<double>> obj_values;
	double h = 0, obj = 0, acc = 0, t1 = 0, t2 = 0;
	vector<double> results;

	//case 1
	if (final_time <= (vmax - spd) / amax)
	{
		h = spd * final_time + 0.5 * amax * pow(final_time, 2) - 1;
		obj = pow(amax, 2) * final_time;
		acc = amax;
		t1 = final_time;
		t2 = final_time;
	}
	else
	{
		h = vmax * final_time - 0.5 * pow(vmax - spd, 2) / amax - 1;
		obj = amax * (vmax - spd);
		acc = amax;
		t1 = (vmax - spd) / amax;
		t2 = final_time;
	}

	if (h == 0)
	{
		obj_values[obj].push_back(obj);
		obj_values[obj].push_back(t1);
		obj_values[obj].push_back(t2);
		obj_values[obj].push_back(acc);
		obj_values[obj].push_back(1);
	}
	else
	{
		obj = 10000;
		obj_values[obj].push_back(obj);
		obj_values[obj].push_back(0);
		obj_values[obj].push_back(0);
		obj_values[obj].push_back(0);
		obj_values[obj].push_back(1);
	}

	//case 2
	int status = GSL_CONTINUE;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;
	double r = 0;
	double x_lo = 0.0, x_hi = final_time;
	gsl_function F;

	struct fittoca_case2_params fa2params = { spd, final_time, vmax, amax };
	F.function = &fittoca_case2;
	F.params = &fa2params;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_set_error_handler_off();
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	while (status == GSL_CONTINUE && iter < max_iter)
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.01);
	}
	gsl_root_fsolver_free(s);

	if (status == GSL_SUCCESS)
	{
		t1 = r;
		t2 = 2.0 * (vmax - spd) / amax - t1;

		if (t1 > 0 && t1 < t2 && t2 < final_time)
		{
			obj = t1 * pow(amax, 2) + pow(amax, 4) / 12.0 * pow(t2 - t1, 3) / pow(spd - vmax + t2 * amax, 2);
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(t1);
				obj_values[obj].push_back(t2);
				obj_values[obj].push_back(amax);
				obj_values[obj].push_back(2);
			}
		}
		else
		{
			obj = 10000;
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(2);
			}
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(2);
		}
	}

	//case 3
	r = 0;
	status = GSL_CONTINUE;
	iter = 0;
	struct fittoca_case3_params fa3params = { spd, final_time, amax };
	x_lo = 0.0;
	x_hi = final_time;
	gsl_function F2;
	const gsl_root_fsolver_type* T2;
	gsl_root_fsolver* s2;

	F2.function = &fittoca_case3;
	F2.params = &fa3params;
	T2 = gsl_root_fsolver_brent;
	s2 = gsl_root_fsolver_alloc(T2);
	gsl_set_error_handler_off();
	gsl_root_fsolver_set(s2, &F2, x_lo, x_hi);

	while (status == GSL_CONTINUE && iter < max_iter)
	{
		iter++;
		status = gsl_root_fsolver_iterate(s2);
		r = gsl_root_fsolver_root(s2);
		x_lo = gsl_root_fsolver_x_lower(s2);
		x_hi = gsl_root_fsolver_x_upper(s2);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.01);
	}
	gsl_root_fsolver_free(s2);

	if (status == GSL_SUCCESS)
	{
		t1 = r;
		if (t1 > 0 && t1 < final_time)
		{
			obj = pow(amax, 2) * (final_time + 2.0 * t1) / 3.0;
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(t1);
				obj_values[obj].push_back(final_time);
				obj_values[obj].push_back(amax);
				obj_values[obj].push_back(3);
			}
		}
		else
		{
			obj = 10000;
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(3);
			}
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(3);
		}
	}

	//case 4
	t1 = 0;
	t2 = (dist_to_sig - final_time * vmax) / (spd + 2.0 * (vmax - spd) / 3.0 - vmax);
	acc = 2.0 * (vmax - spd) / t2;

	if (t2 > 0 && t2 < final_time && acc < amax)
	{
		obj = 4.0 * pow(vmax - spd, 2) / (3.0 * t2);
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(t1);
			obj_values[obj].push_back(t2);
			obj_values[obj].push_back(acc);
			obj_values[obj].push_back(4);
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(4);
		}
	}

	//case 5
	t1 = 0;
	t2 = final_time;
	double final_spd = (3.0 / 2.0) * (dist_to_sig - spd * final_time) / final_time + spd;
	acc = 2.0 * (final_spd - spd) / final_time;

	if (final_spd<vmax && final_spd>spd && acc > 0 && acc < amax)
	{
		obj = 3.0 * pow(dist_to_sig - spd * final_spd, 2) / pow(final_spd, 3);
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(t1);
			obj_values[obj].push_back(t2);
			obj_values[obj].push_back(acc);
			obj_values[obj].push_back(5);
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(5);
		}
	}
	results = obj_values.begin()->second;

	return results;
}

vector<double> ecoand::fittocd(double spd, double dist_to_sig, double final_time)
{
	map<double, vector<double>> obj_values;
	double h = 0, obj = 0, acc = 0, t1 = 0, t2 = 0, final_spd = 0;
	vector<double> results;

	//case 1
	int status = GSL_CONTINUE;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;
	double r = 0;
	double x_lo = 0.0, x_hi = final_time;
	gsl_function F;

	struct fittocd_case1_params fd1params = { spd, final_time, vmin, amin };
	F.function = &fittocd_case1;
	F.params = &fd1params;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_set_error_handler_off();
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	while (status == GSL_CONTINUE && iter < max_iter)
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.01);
	}
	gsl_root_fsolver_free(s);

	if (status == GSL_SUCCESS)
	{
		t1 = r;
		t2 = 2.0 * (vmin - spd) / amin - t1;
		final_spd = vmin;

		if (t1 >= 0 && t1 < t2 && t2 <= final_time)
		{
			obj = t1 * pow(amin, 2) + pow(amin, 4) / 12.0 * pow(t2 - t1, 3) / pow(spd - vmin + t2 * amin, 2);
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(t1);
				obj_values[obj].push_back(t2);
				obj_values[obj].push_back(amin);
				obj_values[obj].push_back(1);
			}
		}
		else
		{
			obj = 10000;

			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(1);
			}
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(1);
		}
	}

	//case 2
	r = 0;
	status = GSL_CONTINUE;
	iter = 0;
	struct fittocd_case2_params fd2params = { spd, final_time, amin };
	const gsl_root_fsolver_type* T2;
	gsl_root_fsolver* s2;
	x_lo = 0.0;
	x_hi = final_time;
	gsl_function F2;

	F2.function = &fittocd_case2;
	F2.params = &fd2params;
	T2 = gsl_root_fsolver_brent;
	s2 = gsl_root_fsolver_alloc(T2);
	gsl_set_error_handler_off();
	gsl_root_fsolver_set(s2, &F2, x_lo, x_hi);

	while (status == GSL_CONTINUE && iter < max_iter)
	{
		iter++;
		status = gsl_root_fsolver_iterate(s2);
		r = gsl_root_fsolver_root(s2);
		x_lo = gsl_root_fsolver_x_lower(s2);
		x_hi = gsl_root_fsolver_x_upper(s2);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.01);
	}
	gsl_root_fsolver_free(s2);

	if (status == GSL_SUCCESS)
	{
		t1 = r;
		final_spd = spd + (t1 + final_time) / 2.0 * amin;
		if (t1 > 0 && t1 < final_time)
		{
			obj = pow(amin, 2) * (final_time + 2.0 * t1) / 3.0;
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(t1);
				obj_values[obj].push_back(final_time);
				obj_values[obj].push_back(amin);
				obj_values[obj].push_back(2);
			}
		}
		else
		{
			obj = 10000;
			if (obj_values.count(obj) == 0)
			{
				obj_values[obj].push_back(obj);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(0);
				obj_values[obj].push_back(2);
			}
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(2);
		}
	}

	//case 3
	t1 = 0;
	t2 = (3.0 * dist_to_sig - 3.0 * final_time * vmin) / (spd - vmin);
	acc = 2.0 * (vmin - spd) / t2;
	final_spd = vmin;
	if (t2 > 0 && t2 < final_time && acc > amin)
	{
		obj = 4.0 * pow(vmin - spd, 2) / (3.0 * t2);
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(t1);
			obj_values[obj].push_back(t2);
			obj_values[obj].push_back(acc);
			obj_values[obj].push_back(3);
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(3);
		}
	}

	//case 4
	t1 = 0;
	t2 = final_time;
	final_spd = (3.0 / 2.0) * (dist_to_sig - spd * final_time) / final_time + spd;
	acc = 2.0 * (final_spd - spd) / final_time;

	if (final_spd <= vmax && final_spd >= vmin)
	{
		obj = 3.0 * pow(dist_to_sig - spd * final_time, 2) / pow(final_time, 3);
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(t1);
			obj_values[obj].push_back(t2);
			obj_values[obj].push_back(acc);
			obj_values[obj].push_back(4);
		}
	}
	else
	{
		obj = 10000;
		if (obj_values.count(obj) == 0)
		{
			obj_values[obj].push_back(obj);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(0);
			obj_values[obj].push_back(4);
		}
	}

	results = obj_values.begin()->second;

	return results;
}

vector<double> ecoand::rorg(double cur_time, double cur_pos, double spd, long next_sig)
{
	vector<double> results;

	double est_time = 0, index_to_go = 0, time_use_acc = 0, time_use_dec = 0;
	double sig_cyc_time = sig_cyc_times[next_sig];
	string sig_state = sig_states[next_sig];
	double time_to_next_red = sig_tm_next_red[next_sig];
	double time_to_next_green = sig_tm_next_green[next_sig];

	if (sig_ini_times.count(next_sig) != 0) //if sig_ini_time already calculated
	{
		est_time = sig_ini_times[next_sig];
	}
	else
	{
		double pos_sig = sig_lengths[next_sig];
		double dist_to_sig = pos_sig - cur_pos;

		if (spd > 0)
			est_time = cur_time + dist_to_sig / spd;

		sig_ini_times[next_sig] = est_time;
	}

	double time_to_sig = est_time - cur_time;
	double over_cyc_time = remainder(time_to_sig, sig_cyc_time);

	if (strcmp(sig_state.c_str(), "GREEN")) // green
	{
		if ((over_cyc_time <= time_to_next_red) || (over_cyc_time >= time_to_next_green))
		{
			index_to_go = 1;
			time_use_acc = time_to_sig;
			time_use_dec = time_to_sig;
		}
		else
		{
			index_to_go = 0;
			time_use_acc = time_to_next_red + (time_to_sig - over_cyc_time);
			time_use_dec = time_to_next_green + (time_to_sig - over_cyc_time);
		}
	}
	else //amber or red
	{
		if ((over_cyc_time >= time_to_next_green) && (over_cyc_time <= time_to_next_red))
		{
			index_to_go = 1;
			time_use_acc = time_to_sig;
			time_use_dec = time_to_sig;
		}
		else
		{
			index_to_go = 0;
			time_use_acc = time_to_next_green + (time_to_sig - over_cyc_time) - (sig_cyc_time - time_to_next_red + time_to_next_green);
			time_use_dec = time_to_next_red + (time_to_sig - over_cyc_time);
		}
	}
	results.push_back(index_to_go);
	results.push_back(time_use_acc);
	results.push_back(time_use_dec);

	return results;
}

void ecoand::set_ctrls(double final_time, double time_step, double acc, double spd)
{
	ctrls.clear();
	double speed = acc * time_step + spd;

	ctrls.push_back(final_time);
	ctrls.push_back(acc);
	ctrls.push_back(speed);
}


void ecoand::set_ctrls(string mode, double cur_time, double t1, double t2, double tf, double acc, double case_number, double obj)
{
	ctrls.clear();
	double t_1, t_2, t_3, a1, a2, b1, b2, a3, b3;

	if (mode == "fittoca")
	{
		if (case_number == 1)
		{
			if (obj != 10000)
			{
				t_1 = t1 + cur_time;
				t_2 = t2 + cur_time;
				a1 = 0;
				b1 = acc;
				a2 = 0;
				b2 = 0;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
			}
		}
		else if (case_number == 2)
		{
			if (obj != 10000)
			{
				t_1 = t1 + cur_time;
				t_2 = t2 + cur_time;
				t_3 = tf;
				a1 = 0;
				b1 = acc;
				a2 = -acc / (t2 - t1);
				b2 = acc;
				a3 = 0;
				b3 = 0;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(t_3);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
				ctrls.push_back(a3);
				ctrls.push_back(b3);
			}
		}
		else if (case_number == 3)
		{
			if (obj != 10000)
			{
				t_1 = t1 + cur_time;
				t_2 = tf;
				a1 = 0;
				b1 = acc;
				a2 = -acc / (tf - t_1);
				b2 = acc;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
			}
		}
		else if (case_number == 4)
		{
			if (obj != 10000)
			{
				t_1 = t2 + cur_time;
				t_2 = tf;
				a1 = -acc / t2;
				b1 = acc;
				a2 = 0;
				b2 = 0;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
			}
		}
		else if (case_number == 5)
		{
			if (obj != 10000)
			{
				t_1 = tf;
				a1 = -acc / (tf - cur_time);
				b1 = acc;

				ctrls.push_back(t_1);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
			}
		}
	}
	else if (mode == "fittocd")
	{
		if (case_number == 1)
		{
			if (obj != 10000)
			{
				t_1 = t1 + cur_time;
				t_2 = t2 + cur_time;
				t_3 = tf;
				a1 = 0;
				b1 = acc;
				a2 = -2 * acc / (t2 - t1);
				b2 = acc;
				a3 = 0;
				b3 = 0;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(t_3);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
				ctrls.push_back(a3);
				ctrls.push_back(b3);
			}
		}
		else if (case_number == 2)
		{
			if (obj != 10000)
			{
				t_1 = t1 + cur_time;
				t_2 = tf;
				a1 = 0;
				b1 = acc;
				a2 = -2 * acc / (tf - t_1);
				b2 = acc;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
			}
		}
		else if (case_number == 3)
		{
			if (obj != 10000)
			{
				t_1 = t2 + cur_time;
				t_2 = tf;
				a1 = -acc / t2;
				b1 = acc;
				a2 = 0;
				b2 = 0;

				ctrls.push_back(t_1);
				ctrls.push_back(t_2);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
				ctrls.push_back(a2);
				ctrls.push_back(b2);
			}
		}
		else if (case_number == 4)
		{
			if (obj != 10000)
			{
				t_1 = tf;
				a1 = -acc / (tf - cur_time);
				b1 = acc;

				ctrls.push_back(t_1);
				ctrls.push_back(a1);
				ctrls.push_back(b1);
			}
		}
	}
	else if (mode == "cruise")
	{
		t_1 = tf;
		a1 = 0;
		b1 = 0;

		ctrls.push_back(t_1);
		ctrls.push_back(a1);
		ctrls.push_back(b1);
	}
}

void ecoand::set_ctrl_modes(vector<double> modes)
{
	ctrl_modes.clear();

	ctrl_modes = modes;
}


//void ecoand::cal_ini_time(double cur_time, double cur_pos, double spd, long next_sig)
//{
//	double est_time = 0;
//	double pos_sig = sig_lengths[next_sig];
//	double dist_to_sig = pos_sig - cur_pos;
//
//	if (spd > 0)
//		est_time = cur_time + dist_to_sig / spd;
//
//	sig_ini_times[next_sig] = est_time;
//}

void ecoand::cal_ini_time(double cur_time, double dist_to_nxt_sig, double spd, long next_sig_id)
{
	double est_time = 0;
	if (spd > 0)
		est_time = cur_time + dist_to_nxt_sig / spd;

	sig_ini_times[next_sig_id] = est_time;
}

vector<double> ecoand::fttoc(double spd, double rho_t, double rho_u, double dist_to_sig) //determine optimal times for free time optimization problem
{
	double v0 = spd;
	double vM = vmax;
	double aM = amax;
	double l = dist_to_sig;
	vector<double> results;
	double t1 = 0, t2 = 0, ud = 0, tf = -1, cas = 0, v;

	if (v0 / vM < 1 - pow(aM, 2) * rho_u / rho_t)
	{
		if (l >= (pow(vM, 2) - pow(v0, 2)) / (2 * aM) + aM * pow(vM, 2) * rho_u / rho_t - pow(aM, 3) * pow(vM, 2) * pow(rho_u, 2) / (6 * pow(rho_t, 2)))
		{
			cas = 1;
			t1 = ((1 - pow(aM, 2) * rho_u / rho_t) * vM - v0) / aM;
			t2 = 2 * aM * vM * rho_u / rho_t;
			ud = rho_t / (2 * rho_u * vM);
			tf = t1 + t2 + (l - (pow(vM, 2) - pow(v0, 2)) / (2 * aM) - aM * pow(vM, 2) * rho_u / rho_t + pow(aM, 3) * pow(vM, 2) * pow(rho_u, 2) / (6 * pow(rho_t, 2))) / vM;
		}
		else
		{
			cas = 2;
			v = sqrt((2 * aM * l + pow(v0, 2)) / (1 + 4 * pow(aM, 2) * rho_u / (rho_t * (1 - rho_u * pow(aM, 2) / rho_t)) + 8 * pow(rho_u, 2) * pow(aM, 4) / (3 * pow(rho_t, 2) * pow(1 - rho_u * pow(aM, 2) / rho_t, 2))));
			t1 = (v - v0) / aM;
			t2 = 2 * aM * (rho_u / rho_t) * v / (1 - rho_u * pow(aM, 2) / rho_t);
			ud = rho_t * (1 - rho_u * pow(aM, 2) / rho_t) / (2 * rho_u * v);
			tf = t1 + t2;
		}
	}
	else
	{
		if (l >= 2 * v0 * sqrt((vM - v0) * vM * rho_u / rho_t) + 4 * (vM - v0) * sqrt((vM - v0) * vM * rho_u / rho_t) / 3)
		{
			cas = 3;
			t1 = 0;
			t2 = 2 * sqrt((vM - v0) * vM * rho_u / rho_t);
			ud = rho_t / (2 * rho_u * vM);
			tf = t2 + (l - 2 * v0 * sqrt((vM - v0) * vM * rho_u / rho_t) - 4 * (vM - v0) * sqrt((vM - v0) * vM * rho_u / rho_t) / 3) / vM;
		}
		else
		{
			cas = 4;

			int status = GSL_CONTINUE;
			int iter = 0, max_iter = 100;
			const gsl_root_fsolver_type* T;
			gsl_root_fsolver* s;
			double r = 0;
			double x_lo = 0.0, x_hi = 1000;
			gsl_function F;

			struct fttoc_params ftparams = { spd, rho_t, rho_u };
			F.function = &fttoc_solve;
			F.params = &ftparams;
			T = gsl_root_fsolver_brent;
			s = gsl_root_fsolver_alloc(T);
			gsl_set_error_handler_off();
			gsl_root_fsolver_set(s, &F, x_lo, x_hi);

			while (status == GSL_CONTINUE && iter < max_iter)
			{
				iter++;
				status = gsl_root_fsolver_iterate(s);
				r = gsl_root_fsolver_root(s);
				x_lo = gsl_root_fsolver_x_lower(s);
				x_hi = gsl_root_fsolver_x_upper(s);
				status = gsl_root_test_interval(x_lo, x_hi,
					0, 0.01);
			}
			gsl_root_fsolver_free(s);

			if (status == GSL_SUCCESS)
			{
				tf = r;
			}
			else
			{
				tf = -1;
			}
			if (tf != -1)
			{
				t1 = 0;
				t2 = tf;
				ud = rho_t / (rho_u * (v0 + sqrt(pow(v0, 2) + rho_t * pow(tf, 2) / rho_u)));
			}
		}
	}
	results.push_back(t1);
	results.push_back(t2);
	results.push_back(ud);
	results.push_back(tf);
	results.push_back(cas);

	return results;
}


void ecoand::get_sig_links(string sig_link_str, string sig_group_str, map<long, double> link_signals, string sig_cyc_time_str)
{
	sig_ctrl_lks.clear();
	sig_ctrl_groups.clear();

	vector<string> ctrl_links_temp, ctrl_groups_temp, cycle_temp;
	boost::split(ctrl_links_temp, sig_link_str, boost::is_any_of(","), boost::token_compress_on);
	boost::split(ctrl_groups_temp, sig_group_str, boost::is_any_of(","), boost::token_compress_on);
	boost::split(cycle_temp, sig_cyc_time_str, boost::is_any_of(","), boost::token_compress_on);

	for (size_t i = 0; i < ctrl_links_temp.size(); i++)
	{
		if (ctrl_links_temp[i] != "")
		{
			long lk = stoi(ctrl_links_temp[i]);
			if (find(sigs.begin(), sigs.end(), lk) == sigs.end())
				sigs.push_back(lk);
			sig_cyc_times[lk] = stol(cycle_temp[i]);

			if (link_signals.count(lk) != 0)
			{
				string sig_group = to_string(link_signals[lk]);
				size_t pos = sig_group.find('.');

				long sig = stol(sig_group.substr(0, pos));
				long group = stol(sig_group.substr(pos + 1, 1));

				sig_ctrl_groups[lk] = group;
				sig_ctrl_lks[lk] = sig;
			}
		}
	}
}



void ecoand::get_sig_states(string sig_state_str, double tm_nxt_green, double tm_nxt_red, double cyc_time, long sig_ID)
{

	sig_states[sig_ID] = sig_state_str;
	sig_tm_next_green[sig_ID] = tm_nxt_green;
	sig_tm_next_red[sig_ID] = tm_nxt_red;
	sig_cyc_times[sig_ID] = cyc_time;
}



void ecoand::get_sig_pos(long cur_lat, long cur_lon, long nxt_sig_lat, long nxt_sig_lon)
{
	// call SPaT message for current position and heading and get the location of traffic light
	dist_to_sig = haversine(cur_lat, cur_lon, nxt_sig_lat, nxt_sig_lon);

}



vector<double> ecoand::cal_final_time_ecoand( map<long, vector<long>> sigs_vehs,
	map<long, vector<double>> sigs_vehs_times, double des_headway, long sig_ID)
{
	double vf = -1, acc_dec = 0;
	vector<double> results;
	long last_veh_id = 0;
	double last_veh_time = 0;
	if (sigs_vehs[sig_ID].size() > 0)
	{
		last_veh_id = sigs_vehs[sig_ID].back();
		last_veh_time = sigs_vehs_times[sig_ID].back();
	}

	if (sigs_vehs[sig_ID].size() == 0 || (last_veh_time > 0 && last_veh_time <= cur_time)) //first vehicle approaching the intersection from one direction (link of the signal head -- unique for each directions)
	{
		if (sig_states[sig_ID] == "GREEN") //if current state is green
		{
			double left_green_time = sig_tm_next_red[sig_ID]; //time until red
			double final_time_green = left_green_time + cur_time; //maximum final time

			double max_time = (sig_lengths[sig_ID] - cur_pos) / vmin + cur_time;
			double min_time = (sig_lengths[sig_ID] - cur_pos) / vmax + cur_time;

			if (sig_ini_times[sig_ID] <= final_time_green) //if initial time is smaller than green end time
			{
				vf = sig_ini_times[sig_ID];
			}
			else if (min_time <= final_time_green) // if with maximum speed, green end time is smaller than minimum time
			{
				vf = final_time_green;
				acc_dec = 1; //acceleration
			}
			else
			{
				double final_time_new = sig_tm_next_green[sig_ID] + cur_time; //time until next green
				vf = final_time_new;
				acc_dec = -1;
			}
		}
		else //amber or red
		{
			double left_red_time = sig_tm_next_green[sig_ID]; //time until green
			double final_time_green_start = left_red_time + cur_time; // earliest available green, minimum final time
			double final_time_green_end = sig_tm_next_red[sig_ID] + cur_time; //time until next red, maximum final time

			double max_time = (sig_lengths[sig_ID] - cur_pos) / vmin + cur_time; //deceleration time
			double min_time = (sig_lengths[sig_ID] - cur_pos) / vmax + cur_time; //acceleration time

			if (sig_ini_times[sig_ID] >= final_time_green_start && sig_ini_times[sig_ID] <= final_time_green_end) //if initial time is smaller than green end time
			{
				vf = sig_ini_times[sig_ID];
			}
			else if (sig_ini_times[sig_ID] < final_time_green_start) // if initial time is earlier than red end time -- deceleration
			{
				if (max_time >= final_time_green_start)
				{
					vf = final_time_green_start;
					acc_dec = -1; //deceleration
				}
				else
				{
					vf = -1; //need to stop for red
				}
			}
			else if (sig_ini_times[sig_ID] > final_time_green_end)
			{
				double final_time_new = final_time_green_start + sig_cyc_times[sig_ID]; //next green time

				if (min_time <= final_time_green_end) // if with maximum speed, green end time is smaller than minimum time
				{
					vf = final_time_green_end;
					acc_dec = 1; //acceleration
				}
				else if (max_time >= final_time_new)
				{
					vf = final_time_new;
					acc_dec = -1;
				}
				else
				{
					vf = -1;//need to stop for red
				}
			}
		}
	}
	else if (last_veh_id != veh_id)
	{
		if (sig_states[sig_ID] == "GREEN") //if current state is green, then check last vehicle entry time
		{
			double left_green_time = sig_tm_next_red[sig_ID]; //time until red
			double final_time_green = left_green_time + cur_time;
			double max_time = (sig_lengths[sig_ID] - cur_pos) / vmin + cur_time;
			double min_time = (sig_lengths[sig_ID] - cur_pos) / vmax + cur_time; //useful

			if (last_veh_time + des_headway <= final_time_green)
			{
				if (sig_ini_times[sig_ID] <= final_time_green)
				{
					vf = max(sig_ini_times[sig_ID], last_veh_time + des_headway);
					acc_dec = sig_ini_times[sig_ID] < last_veh_time + des_headway ? -1 : 0;
				}
				else
				{
					if (min_time <= final_time_green)
					{
						vf = max(last_veh_time + des_headway, min_time);
						acc_dec = 1;
					}
					else
					{
						double final_time_new = sig_tm_next_green[sig_ID] + cur_time; //time until next green
						if (max_time > final_time_new)
						{
							vf = final_time_new;
							acc_dec = -1; //deceleration
						}
						else
							vf = -1;
					}
				}
			}
			else //last_veh_time + des_headway > final_time_green
			{
				double final_time_new = sig_tm_next_green[sig_ID] + cur_time; //time until next green
				double final_time_new_exit = sig_tm_next_green[sig_ID] + sig_cyc_times[sig_ID] + cur_time;//time when next green ends

				if (last_veh_time <= final_time_green)
				{
					if (sig_ini_times[sig_ID] <= final_time_new)
					{
						if (max_time >= final_time_new) // if with maximum speed, green end time is smaller than minimum time
						{
							vf = final_time_new;
							acc_dec = -1; //deceleration
						}
						else
							vf = -1;
					}
					else if (sig_ini_times[sig_ID] > final_time_new && sig_ini_times[sig_ID] <= final_time_new_exit)
					{
						vf = sig_ini_times[sig_ID]; //constant speed
					}
					else if (sig_ini_times[sig_ID] >= final_time_new_exit)
					{
						if (min_time <= final_time_new_exit)
						{
							vf = final_time_new_exit;
							acc_dec = 1;
						}
						else
						{
							vf = -1;
						}
					}
				}
				else //last_veh_time > final_time_green -- then last_veh_time should be at least final_time_new
				{
					if (last_veh_time + des_headway <= final_time_new_exit)
					{
						if (sig_ini_times[sig_ID] > final_time_new && sig_ini_times[sig_ID] <= final_time_new_exit)
						{
							vf = max(sig_ini_times[sig_ID], last_veh_time + des_headway);
							acc_dec = sig_ini_times[sig_ID] < last_veh_time + des_headway ? -1 : 0;
						}
						else if (sig_ini_times[sig_ID] <= final_time_new)
						{
							if (max_time > final_time_new && max_time <= final_time_new_exit) // if with maximum speed, green end time is smaller than minimum time
							{
								vf = min(max_time, last_veh_time + des_headway);
								acc_dec = -1; //deceleration
							}
							else if (max_time > final_time_new_exit)
							{
								vf = last_veh_time + des_headway;
								acc_dec = -1; //deceleration
							}
							else
								vf = -1;
						}
						else // sig_ini_times[sig] > final_time_new_exit
						{
							if (min_time <= final_time_new_exit)
							{
								vf = max(min_time, last_veh_time + des_headway);
								acc_dec = 1;
							}
							else//(min_time > final_time_new_exit)
							{
								vf = -1;
							}
						}
					}
					else
					{
						vf = -1;
					}
				}
			}
		}
		else // red or amber
		{
			double left_red_time = sig_tm_next_green[sig_ID]; //time until green
			double final_time_enter = left_red_time + cur_time; // earliest available green
			double final_time_exit = sig_tm_next_red[sig_ID] + cur_time; //time until next red
			double max_time = (sig_lengths[sig_ID] - cur_pos) / vmin + cur_time;
			double min_time = (sig_lengths[sig_ID] - cur_pos) / vmax + cur_time; //useful

			if (last_veh_time + des_headway > final_time_enter && last_veh_time + des_headway <= final_time_exit) // if potential final time for this vehicle is green
			{
				if (sig_ini_times[sig_ID] == last_veh_time + des_headway)
					acc_dec = 0;
				else if (sig_ini_times[sig_ID] > final_time_exit)
				{
					if (min_time < final_time_exit)
					{
						vf = max(min_time, last_veh_time + des_headway);
						acc_dec = 1;
					}
					else
					{
						vf = -1;
					}
				}
				else if (sig_ini_times[sig_ID] <= final_time_enter)
				{
					if (max_time > final_time_exit)
					{
						vf = last_veh_time + des_headway;
						acc_dec = -1;
					}
					else if (max_time > final_time_enter)
					{
						vf = min(max_time, last_veh_time + des_headway);
						acc_dec = -1;
					}
					else
					{
						vf = -1;
					}
				}
				else
				{
					vf = max(sig_ini_times[sig_ID], last_veh_time + des_headway);
					acc_dec = sig_ini_times[sig_ID] < vf ? -1 : 1;
				}
			}
			else if (last_veh_time + des_headway > final_time_exit) // if potential final time for this vehicle is next red
			{
				vf = -1;
			}
		}
	}
	results.push_back(vf);
	results.push_back(acc_dec);

	return results;
}



