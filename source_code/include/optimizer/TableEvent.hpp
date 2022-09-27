#ifndef TABLEEVENT_H_
#define TABLEEVENT_H_


#include <iostream>

#include "CarIterator_header.hpp"


/// <summary>
/// Possible actions of the vehicle.
/// </summary>
enum class Actions
{
	// Car generation.
	START,
	// Car movement (movement across one edge).
	RUNNING,
	// Car starts charging.
	CHARGING,
	// Car just ended charging
	END_CHARGING,
	// Car waiting in the line for charging.
	WAITING_LINE,
	// Car waiting for the return or arrived to the end.
	WAITING_RETURN
};


/// <summary>
/// Class to store info about the simulation event (for discrete simulation).
/// </summary>
class TableEvent
{
public:
	// Time of the start of the event
	double ActionTime_;
	// Iterator to `CarMap` contaioner (pair of car ID and the Car object itself).
	int64_t CarID_;
	// To identify action.
	Actions Action_;


	TableEvent(int64_t carID, double actionTime, Actions action);
};

#endif