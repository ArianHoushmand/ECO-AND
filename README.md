# Economical Arrival and Departure (ECO-AND)

This repository contains C++ codes for calculating the optimal speed and acceleration profiles for an autonomous vehicle approaching a traffic light. We assume the the vehicle has V2V communication capabilities, and through SPaT (Signal Phase and Timing) messages we can receive the traffic light information (i.e. cycle time, remaining time to red, and green). The objective function in the optimal control problem is to minimize the overall fuel consumption and travel time for the car approaching the traffic light without ever fully stopping at the intersection. 

 



## Code Architecture

### 1- ECO_AND_Final.cpp
 This is the main code for the simulation.
### 2- ecoand.cpp
Includes the ecoand class in which we calculate the optimal control problem
### 3- ecoand.h
Includes the header file for the ecoand class.

## Related Publications
* Meng, Xiangyu, and Christos G. Cassandras. "Optimal control of autonomous vehicles for non-stop signalized intersection crossing." In 2018 IEEE Conference on Decision and Control (CDC), pp. 6988-6993. IEEE, 2018.


## Contributers
Arian Houshmand and Liuhui Zhao

