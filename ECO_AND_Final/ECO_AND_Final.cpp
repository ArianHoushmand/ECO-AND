// ECO_AND_Finall_Test.cpp
//

#include <iostream>
#include "ecoand.h"
#include "utility.h"
#include <list>
#include <map>
#include <vector>
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


vector<double> cal_ecoand(long current_veh_id, double current_spd, double dist_traveled,
	long lead_id, double lead_tf, double lead_dist, double lead_spd_diff, long sig_id, string sig_state,
	double sig_tm_nxt_green, double sig_tm_nxt_red, double sig_cyc_time, double dist_to_sig, double current_time);

int main()
{
	//Inputes:
	// Ego  vehicle info
	double dist_traveled = 10; // distance travelled since the beginning
	double current_time = 19.5;
	long current_veh_id = 1;
	double current_spd = 12;
	// Lead vehicle information
	long lead_id = 2; //-1 if there is no car in the front, otherwise pass the vehicle ID
	double lead_spd_diff = 0; // spd_lead - spd_ego
	double lead_dist = 40; // distance of the lead vehicle from the ego car
	double lead_tf = 28; // final arriavla time of the lead vehicle at the intersection (should be derived from VISSIM)
	// Traggic light info
	string sig_state = "GREEN"; // options: "GREEN", "RED", "YELLOW"
	double sig_tm_nxt_green = 10.5; // Remaining time to the next green light
	double sig_tm_nxt_red = 0.5;// Remaining time to the next red light
	double sig_cyc_time = 10; // total cycle time of the light: total red time + total greem time
	long sig_id = 1; // Signal ID (we can leave it as it since we only have one traffic light)
	double dist_to_sig = 140; // distance of the ego car to the traffic light

	// Output vector:
	// output[0]: desired acceleration
	// output[1]: terminal time (arriavle time) at the intersection
	// output[2]: mode: 3 if in eco_and, 0 otherwise
	vector<double> output;

	output = cal_ecoand(current_veh_id, current_spd, dist_traveled, lead_id, lead_tf,
		lead_dist, lead_spd_diff, sig_id, sig_state, sig_tm_nxt_green, sig_tm_nxt_red,
		sig_cyc_time, dist_to_sig, current_time);

std:cout << "Acceleration: " << output[0] << "\n";
	std::cout << "Tf: " << output[1] << "\n";
	std::cout << "Control Mode " << output[2] << "\n";



	//////////////////////////////////////////////


}


vector<double> cal_ecoand(long current_veh_id, double current_spd, double dist_traveled,
	long lead_id, double lead_tf, double lead_dist, double lead_spd_diff, long sig_id, string sig_state,
	double sig_tm_nxt_green, double sig_tm_nxt_red, double sig_cyc_time, double dist_to_sig, double current_time)
{
	ecoand c;
	vector<double> results;
	double tf = -1;

	double max_spd = 16; // 36mph for urban intersection
	double min_spd = 2; // 5mph (queue start speed in vissim)
	double min_acc = -2.0; // (desired deceleration for normal vehicles)
	double max_acc = 2.0; // (desired acceleration for normal vehicles)

	double ctrl_len = 150; // min distance to traffic light for activating ecoAND

	double des_headway = 1.2;
	double safe_dist = 1.5;
	double safe_headway = 1.2;
	double desired_acceleration = 0.0;
	double lead_veh_length = 2;


	long mode = 0;

	map<long, vector<long>> sigs_vehs; //intersection vehicles -- intersection, list of vehicle ids
	map<long, vector<double>> sigs_vehs_times; //intersection vehicle final times -- vehicle id, vehicle final time



	if (lead_id == -1)
	{
		c.init(current_veh_id, current_time, dist_traveled, current_spd);
	}
	else
	{
		c.init(current_veh_id, current_time, dist_traveled, current_spd, lead_id, lead_spd_diff, lead_dist);
		sigs_vehs[sig_id].push_back(lead_id);
		sigs_vehs_times[sig_id].push_back(lead_tf);
	}
	c.mode = 3;
	c.set_limits(max_spd, min_spd, max_acc, min_acc);

	c.dist_to_sig = dist_to_sig;

	if (c.dist_to_sig < ctrl_len)
	{
		c.cal_ini_time(current_time, c.dist_to_sig, current_spd, sig_id); //calculate initial time, and then add to intersection record
		c.in_ecoand = 1;

		c.get_sig_states(sig_state, sig_tm_nxt_green, sig_tm_nxt_red, sig_cyc_time, sig_id);
		vector<double> final_results = c.cal_final_time_ecoand(sigs_vehs, sigs_vehs_times, des_headway, sig_id);
		c.set_ctrl_modes(final_results);
		c.sig_times[sig_id] = final_results[0];
		c.sig_in_ctrl_times[sig_id] = current_time;

		//std::cout << final_results[0] << "\n";
		//std::cout << final_results[1] << "\n";

		double acc_com;
		vector<double> ctrl_results;

		if (final_results[0] == -1)
		{
			c.in_ctrl = 0;
			c.mode = 0;
			mode = 0;
		}
		else
		{
			c.in_ctrl = 1;
			c.mode = 3;
			mode = 3;

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
					c.mode = 0;
					mode = 0;
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
					c.mode = 0;
					mode = 0;
				}
			}
		}

		sigs_vehs[sig_id].push_back(current_veh_id);
		sigs_vehs_times[sig_id].push_back(final_results[0]);
		tf = final_results[0];
	}
	else
	{
		c.in_ecoand = 0;
		c.mode = 0;
		mode = 0;
	}

	//std::cout << c.in_ctrl << "\n";
	//std::cout << c.in_ecoand << "\n";
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
		if (lead_id == -1 || (lead_id != -1 && lead_dist - lead_veh_length >= current_spd * safe_headway + safe_dist))
			desired_acceleration = acc;
		else
		{
			c.in_ctrl = 0;
			c.mode = 0;
		}
	}

	results.push_back(acc);
	results.push_back(tf);
	results.push_back(c.mode);

	return results;


}


