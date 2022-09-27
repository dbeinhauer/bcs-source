#ifndef CAR_H_
#define CAR_H_


#include "Vehicle.hpp"


class Car : public Vehicle
{
public:
	Car(double startTime, double waiting, MapPosition start, MapPosition end,
		double consumption, double batteryLevel, double relativeSpeed);

	// Car actions:
	void StartReturning(double actualTime);
	
	// Get information about the car:
	bool IsReturning();
	double GetWaitingTime();

private:
	// How long will the car wait until it starts returning to start.
	double waitingTime_;
	// Flag if car is returning to the start.
	bool returning_;
};

#endif