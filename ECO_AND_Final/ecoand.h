#pragma once

#include <stdio.h>
#include <vector>
#include <map>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <boost\algorithm\string\split.hpp>
#include <boost\algorithm\string\classification.hpp>

using namespace std;

class ecoand
{
public:
	ecoand();
	~ecoand();

public:
	long veh_id, cur_ln, cur_lk;
	double cur_time, cur_pos, cur_spd, cur_acc;
	long current_lat, current_lon; // lat long position of the car (from GPS)
	double vmin, vmax, amin, amax;
	long nxt_sig_lat, nxt_sig_lon;
	//long next_sig, next_sig_group;
	//string next_sig_state;
	//double next_sig_red_len, next_sig_green_len, time_to_next_green, time_to_next_red, next_sig_dist;
	vector<double> ctrls; //next_sig_final_time, next_sig_acc, next_sig_speed;
	vector<double> ctrl_modes;

	long lead_id;
	double lead_spd, lead_dist, lead_acc;
	long in_ecoand, mode;
	long in_ctrl;
	double rho_u, rho_t; //weights of acceleration and time
	double dist_to_sig; // distance to the first signal ahead in meters

	map<long, vector<long>> conflict_links; //to store conflict links
	vector<long> rt_links; // to store all the links en route
	vector<long> sigs; // to store the signal controller link in current order

	map<long, long> sig_cyc_times;
	map<long, long> sig_ctrl_groups; //to store signal controllers
	map<long, long> sig_ctrl_lks; //to store the link information of the signal heads
	map<long, double> sig_lengths; //to store segment lengths for vehicle route (devided by signals)

	map<long, double> sig_tm_next_green; //to store dynamic (time to next green) information
	map<long, double> sig_tm_next_red; // to store dynamic (time to next red) information
	map<long, string> sig_states; // to store current signal state information

	map<long, double> sig_times; //to store assigned final times at intersections
	map<long, double> sig_ini_times; //to store initial final times at intersections
	map<long, double> sig_in_ctrl_times; //to store entry times at intersections

public:
	void init(long id, double tm, double pos, double spd);
	void init(long id, double tm, double pos, double spd, long f_id, double f_spd_diff, double f_dist);
	void update_status(long lk, long lane, double tm, double spd, double acc, double pos);
	void update_status(long lk, long lane, double tm, double spd, double acc, double pos, long f_id, double f_acc, double f_spd_diff, double f_dist);

	vector<double> fittoca(double spd, double dist_to_sig, double final_time);
	vector<double> fittocd(double spd, double dist_to_sig, double final_time);
	vector<double> fttoc(double spd, double rho_t, double rho_u, double dist_to_sig);
	vector<double> rorg(double cur_time, double cur_pos, double spd, long next_sig);
	void set_limits(double max_spd, double min_spd, double max_acc, double min_acc);
	void cal_ini_time(double cur_time, double dist_to_sig, double spd, long next_sig);
	void get_sig_links(string sig_link_str, string sig_group_str, map<long, double> link_signals, string sig_cyc_time_str);
	void get_sig_states(string sig_state_str, double tm_nxt_green, double tm_nxt_red, double cyc_tm, long sig_ID);
	void get_sig_pos(long cur_lat, long cur_lon, long sig_lat, long sig_lon);

	void set_ctrls(double final_time, double time_step, double acc, double spd);
	void set_ctrls(string mode, double cur_time, double t1, double t2, double tf, double acc, double case_number, double obj);

	void set_ctrl_modes(vector<double> modes);

	vector<double> cal_final_time_ecoand(map<long, vector<long>> sigs_vehs,
		map<long, vector<double>> sigs_vehs_times, double des_headway, long sig_ID);


	//long double solve_cubic_ts(double tc, double tm);
	//long double solve_cubic_tm(double tc, double ts, double l, double sigma, double delta);
};
