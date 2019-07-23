// ECO_AND_Arian_Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ecoand.h"
#include "utility.h"
#include <list>
#include <map>
#include <vector>
//#include <wingdi.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <iostream>
#include <ctime>
#include <functional>
#include <utility>
#include <set>
#include <boost\bind.hpp>
#include <cmath>
#include <math.h>


#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */

//function declaration
//static double haversine(double lat1, double lon1, double lat2, double lon2);

map<long, vector<long>> sigs_vehs; //intersection vehicles -- intersection, list of vehicle ids
map<long, vector<double>> sigs_vehs_times; //intersection vehicle final times -- vehicle id, vehicle final time

double des_headway = 1.2;
double current_time = 20;

double max_spd_highway = 23; // 50mph for highway
double max_spd_urban = 16; // 40mph for urban intersection
double min_spd = 2; // 5mph (queue start speed in vissim)
double min_acc = -2.0; // (desired deceleration for normal vehicles)
double max_acc = 2.0; // (desired acceleration for normal vehicles)

double safe_dist = 1.5;
double safe_headway = 1.2;
double  desired_acceleration = 0.0;
vector<double> cal_final_time_ecoand(ecoand& e, long sig);

int main()
{
	ecoand c;
	int leadId = -1; //there is no car in the front
	double dist_traveled = 10; // distance travelled since the beginning of ECOAND mode
	double current_time = 20;
	long current_veh_id = 1;
	double current_spd = 12;
	double current_acc = 0;
	long driving_state = 0;
	double self_veh_length = 2;
	double lead_veh_length = 2;
	double ctrl_len = 150; // min distance to traffic light for activating ecoAND
	long current_veh_lat = 42.348535;
	long current_veh_lon = -71.116543;
	double current_veh_head = 0;
	double dist_to_sig = 145;

	long sig_lat = 42.348732;
	long sig_lon = -71.118099;

	double lead_spd_diff = 0;
	double lead_dist = 200;
	double lead_acc = 0;

	double current_position = 0;
	double max_acceleration = 0;
	long switch_to_car_following = 0;

	double max_spd_highway = 23; // 50mph for highway
	double max_spd_urban = 16; // 40mph for urban intersection
	double min_spd = 2; // 5mph (queue start speed in vissim)
	double min_acc = -2.0; // (desired deceleration for normal vehicles)
	double max_acc = 2.0; // (desired acceleration for normal vehicles)


	string sig_state = "GREEN";
	double sig_tm_nxt_green = 12;
	double sig_tm_nxt_red = 0.5;
	double sig_cyc_time = 10;
	long sig_id = 1;


	double des_headway = 1.2;
	double safe_dist = 1.5;
	double safe_headway = 1.2;

	if (leadId == -1)
		c.init(current_veh_id, current_veh_lat, current_veh_lon, current_time, dist_traveled, current_spd, current_acc);
	else
		c.init(current_veh_id, current_veh_lat, current_veh_lon, current_time, dist_traveled, current_spd, current_acc, leadId, lead_acc, lead_spd_diff, lead_dist);
	c.mode = 3;

	c.set_limits(max_spd_urban, min_spd, max_acc, min_acc);
	//get_sig_pos(current_veh_lat, current_veh_lon, sig_lat, sig_lon);

	c.dist_to_sig = dist_to_sig;

	if (c.dist_to_sig < ctrl_len)
	{
		c.cal_ini_time(current_time, c.dist_to_sig, current_spd, sig_id); //calculate initial time, and then add to intersection record
		c.in_ecoand = 1;

		c.get_sig_states(sig_state, sig_tm_nxt_green, sig_tm_nxt_red, sig_cyc_time, sig_id);
		vector<double> final_results = cal_final_time_ecoand(c, sig_id);
		c.set_ctrl_modes(final_results);
		c.sig_times[sig_id] = final_results[0];
		c.sig_in_ctrl_times[sig_id] = current_time;

		std::cout << final_results[0] << "\n";
		std::cout << final_results[1] << "\n";

		double acc_com;
		vector<double> ctrl_results;

		if (final_results[0] == -1)
		{
			c.in_ctrl = 0;
		}
		else
		{
			c.in_ctrl = 1;

			if (final_results[1] == 0)
			{
				acc_com = 0;
				c.set_ctrls("cruise", current_time, 0, 0, final_results[0], 0, 0, 0);
			}
			else if (final_results[1] == -1) //deceleration
			{								//0: objective value; 1: t1; 2: t2; 3: acc; 4: case_number
				ctrl_results = c.fittocd(current_spd, c.dist_to_sig, final_results[0] - current_time);
				acc_com = ctrl_results[3];
				if (acc_com > min_acc && acc_com <= max_acc)
				{
					c.set_ctrls("fittocd", current_time, ctrl_results[1], ctrl_results[2], final_results[0], ctrl_results[3], ctrl_results[4], ctrl_results[0]);
				}
				else
				{
					c.in_ctrl = 0;
				}
			}
			else if (final_results[1] == 1)
			{
				ctrl_results = c.fittoca(current_spd, c.dist_to_sig, final_results[0] - current_time);
				acc_com = ctrl_results[3];
				if (acc_com > min_acc && acc_com <= max_acc)
				{
					c.set_ctrls("fittoca", current_time, ctrl_results[1], ctrl_results[2], final_results[0], ctrl_results[3], ctrl_results[4], ctrl_results[0]);
				}
				else
				{
					c.in_ctrl = 0;
				}
			}
		}

		sigs_vehs[sig_id].push_back(current_veh_id);
		sigs_vehs_times[sig_id].push_back(final_results[0]);
	}
	else
	{
		c.in_ecoand = 0;
	}

	//std::cout << "Hello World!\n";
	std::cout << c.in_ctrl << "\n";
	std::cout << c.in_ecoand << "\n";
	std::cout << "Hi" << "\n";
//////////////////////////////////////////////////
	double acc = 0;
	if (c.in_ctrl == 1)
	{
		double t1 = 0, t2 = 0, t3 = 0, a1 = 0, a2 = 0, b1 = 0, b2 = 0, a3 = 0, b3 = 0;

		long count = c.ctrls.size();
		vector<double> ego_ctrls = c.ctrls;
		double in_ctrl_time = c.sig_in_ctrl_times[sig_id];

		if (count == 3)
		{//first element is the final time
			t1 = ego_ctrls[0];
			a1 = ego_ctrls[1];
			b1 = ego_ctrls[2];
			if (current_time <= t1)
				acc = a1 * (current_time - in_ctrl_time) + b1;

		}
		else if (count == 6)
		{
			t1 = ego_ctrls[0];
			t2 = ego_ctrls[1];
			a1 = ego_ctrls[2];
			b1 = ego_ctrls[3];
			a2 = ego_ctrls[4];
			b2 = ego_ctrls[5];

			if (current_time <= t1)
				acc = a1 * (current_time - in_ctrl_time) + b1;
			else if (current_time <= t2)
				acc = a2 * (current_time - t1) + b2;
		}
		else if (count == 9)
		{
			t1 = ego_ctrls[0];
			t2 = ego_ctrls[1];
			t3 = ego_ctrls[2];
			a1 = ego_ctrls[3];
			b1 = ego_ctrls[4];
			a2 = ego_ctrls[5];
			b2 = ego_ctrls[6];
			a3 = ego_ctrls[7];
			b3 = ego_ctrls[8];

			if (current_time <= t1)
				acc = a1 * (current_time - in_ctrl_time) + b1;
			else if (current_time <= t2)
				acc = a2 * (current_time - t1) + b2;
			else if (current_time <= t3)
				acc = a3 * (current_time - t2) + b3;
		}
		//safety check
		if (leadId == -1 || (leadId != -1 && lead_dist - lead_veh_length >= current_spd * safe_headway + safe_dist))
			desired_acceleration = acc;
		else
		{
			c.in_ctrl = 0;
		}
	}
	std::cout << desired_acceleration << "\n";
	std::cout << c.ctrl_modes[0] << "\n";

//////////////////////////////////////////////


}


//////////////////////==Global functions==////////////////////////
vector<double> cal_final_time_ecoand(ecoand& e, long sig)
{
	double vf = -1, acc_dec = 0;
	vector<double> results;
	long last_veh_id = 0;
	double last_veh_time = 0;
	if (sigs_vehs[sig].size() > 0)
	{
		last_veh_id = sigs_vehs[sig].back();
		last_veh_time = sigs_vehs_times[sig].back();
	}

	if (sigs_vehs[sig].size() == 0 || (last_veh_time > 0 && last_veh_time <= current_time)) //first vehicle approaching the intersection from one direction (link of the signal head -- unique for each directions)
	{
		if (e.sig_states[sig] == "GREEN") //if current state is green
		{
			double left_green_time = e.sig_tm_next_red[sig]; //time until red
			double final_time_green = left_green_time + current_time; //maximum final time

			double max_time = (e.dist_to_sig) / min_spd + current_time;
			double min_time = (e.dist_to_sig) / max_spd_urban + current_time;

			if (e.sig_ini_times[sig] <= final_time_green) //if initial time is smaller than green end time
			{
				vf = e.sig_ini_times[sig];
			}
			else if (min_time <= final_time_green) // if with maximum speed, green end time is smaller than minimum time
			{
				vf = final_time_green;
				acc_dec = 1; //acceleration
			}
			else
			{
				double final_time_new = e.sig_tm_next_green[sig] + current_time; //time until next green
				vf = final_time_new;
				acc_dec = -1;
			}
		}
		else //amber or red
		{
			double left_red_time = e.sig_tm_next_green[sig]; //time until green
			double final_time_green_start = left_red_time + current_time; // earliest available green, minimum final time
			double final_time_green_end = e.sig_tm_next_red[sig] + current_time; //time until next red, maximum final time

			double max_time = (e.dist_to_sig) / min_spd + current_time; //deceleration time
			double min_time = (e.dist_to_sig) / max_spd_urban + current_time; //acceleration time

			if (e.sig_ini_times[sig] >= final_time_green_start && e.sig_ini_times[sig] <= final_time_green_end) //if initial time is smaller than green end time
			{
				vf = e.sig_ini_times[sig];
			}
			else if (e.sig_ini_times[sig] < final_time_green_start) // if initial time is earlier than red end time -- deceleration
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
			else if (e.sig_ini_times[sig] > final_time_green_end)
			{
				double final_time_new = final_time_green_start + e.sig_cyc_times[sig]; //next green time

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
	else if (last_veh_id != e.veh_id)
	{
		if (e.sig_states[sig] == "GREEN") //if current state is green, then check last vehicle entry time
		{
			double left_green_time = e.sig_tm_next_red[sig]; //time until red
			double final_time_green = left_green_time + current_time;
			double max_time = (e.dist_to_sig) / min_spd + current_time;
			double min_time = (e.dist_to_sig) / max_spd_urban + current_time; //useful

			if (last_veh_time + des_headway <= final_time_green)
			{
				if (e.sig_ini_times[sig] <= final_time_green)
				{
					vf = max(e.sig_ini_times[sig], last_veh_time + des_headway);
					acc_dec = e.sig_ini_times[sig] < last_veh_time + des_headway ? -1 : 0;
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
						double final_time_new = e.sig_tm_next_green[sig] + current_time; //time until next green
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
				double final_time_new = e.sig_tm_next_green[sig] + current_time; //time until next green
				double final_time_new_exit = e.sig_tm_next_green[sig] + e.sig_cyc_times[sig] + current_time;//time when next green ends

				if (last_veh_time <= final_time_green)
				{
					if (e.sig_ini_times[sig] <= final_time_new)
					{
						if (max_time >= final_time_new) // if with maximum speed, green end time is smaller than minimum time
						{
							vf = final_time_new;
							acc_dec = -1; //deceleration
						}
						else
							vf = -1;
					}
					else if (e.sig_ini_times[sig] > final_time_new && e.sig_ini_times[sig] <= final_time_new_exit)
					{
						vf = e.sig_ini_times[sig]; //constant speed
					}
					else if (e.sig_ini_times[sig] >= final_time_new_exit)
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
						if (e.sig_ini_times[sig] > final_time_new && e.sig_ini_times[sig] <= final_time_new_exit)
						{
							vf = max(e.sig_ini_times[sig], last_veh_time + des_headway);
							acc_dec = e.sig_ini_times[sig] < last_veh_time + des_headway ? -1 : 0;
						}
						else if (e.sig_ini_times[sig] <= final_time_new)
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
						else //e.sig_ini_times[sig] > final_time_new_exit
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
			double left_red_time = e.sig_tm_next_green[sig]; //time until green
			double final_time_enter = left_red_time + current_time; // earliest available green
			double final_time_exit = e.sig_tm_next_red[sig] + current_time; //time until next red
			double max_time = (e.dist_to_sig) / min_spd + current_time;
			double min_time = (e.dist_to_sig) / max_spd_urban + current_time; //useful

			if (last_veh_time + des_headway > final_time_enter && last_veh_time + des_headway <= final_time_exit) // if potential final time for this vehicle is green
			{
				if (e.sig_ini_times[sig] == last_veh_time + des_headway)
					acc_dec = 0;
				else if (e.sig_ini_times[sig] > final_time_exit)
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
				else if (e.sig_ini_times[sig] <= final_time_enter)
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
					vf = max(e.sig_ini_times[sig], last_veh_time + des_headway);
					acc_dec = e.sig_ini_times[sig] < vf ? -1 : 1;
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




