// ECO_AND_Arian_Test.cpp : This file contains the 'main' function. Program execution begins and ends there.
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


int main()
{
	ecoand c;
	int leadId = -1; //there is no car in the front
	double dist_traveled = 10; // distance travelled since the beginning of ECOAND mode
	double current_time = 19.5;
	long current_veh_id = 1;
	double current_spd = 12;
	double current_acc = 0;
	long driving_state = 0;
	double lead_veh_length = 2;
	double ctrl_len = 150; // min distance to traffic light for activating ecoAND
	double current_veh_head = 0;
	double dist_to_sig = 140;


	double lead_spd_diff = 0;
	double lead_dist = 200;
	double lead_acc = 0;

	double current_position = 0;
	double max_acceleration = 0;
	long switch_to_car_following = 0;
	double  desired_acceleration = 0.0;

	double max_spd_highway = 23; // 50mph for highway
	double max_spd_urban = 16; // 40mph for urban intersection
	double min_spd = 2; // 5mph (queue start speed in vissim)
	double min_acc = -2.0; // (desired deceleration for normal vehicles)
	double max_acc = 2.0; // (desired acceleration for normal vehicles)


	string sig_state = "GREEN";
	double sig_tm_nxt_green = 10.5;
	double sig_tm_nxt_red = 0.5;
	double sig_cyc_time = 10;
	long sig_id = 1;

	double des_headway = 1.2;
	double safe_dist = 1.5;
	double safe_headway = 1.2;

	map<long, vector<long>> sigs_vehs; //intersection vehicles -- intersection, list of vehicle ids
	map<long, vector<double>> sigs_vehs_times; //intersection vehicle final times -- vehicle id, vehicle final time



	if (leadId == -1)
		c.init(current_veh_id, current_time, dist_traveled, current_spd, current_acc);
	else
		c.init(current_veh_id, current_time, dist_traveled, current_spd, current_acc, leadId, lead_acc, lead_spd_diff, lead_dist);
	c.mode = 3;
	c.set_limits(max_spd_urban, min_spd, max_acc, min_acc);

	c.dist_to_sig = dist_to_sig;

	if (c.dist_to_sig < ctrl_len)
	{
		c.cal_ini_time(current_time, c.dist_to_sig, current_spd, sig_id); //calculate initial time, and then add to intersection record
		c.in_ecoand = 1;

		c.get_sig_states(sig_state, sig_tm_nxt_green, sig_tm_nxt_red, sig_cyc_time, sig_id);
		//vector<double> final_results = cal_final_time_ecoand(c, sig_id);
		vector<double> final_results = c.cal_final_time_ecoand(sigs_vehs, sigs_vehs_times, des_headway, sig_id);
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


