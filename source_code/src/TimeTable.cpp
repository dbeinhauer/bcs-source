#include "optimizer/TimeTable.hpp"


/// <summary>
/// Initializes the timetable.
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
TimeTable::TimeTable(SimulationParameters& simulationParameters)
	:trafficSimulator_(simulationParameters.NumClosestStations, 
		simulationParameters.SegmentLength)
{
	this->ActualTime_ = 0;
}


/// <summary>
/// Loads map of the simulator from the given source and 
/// sets proper simulator parameters.
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <returns>Returns `true` if loading was sucessfull, else `false`.</returns>
bool TimeTable::LoadMap(SimulationParameters& simulationParameters)
{
	return this->trafficSimulator_.LoadMap(simulationParameters);
}


/// <summary>
/// Randomly generates charging stations.
/// </summary>
/// <param name="numStations">Number of station to be generated.</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
void TimeTable::RandomlyGenerateChargingStations(int32_t numStations,
	SimulationParameters& simulationPararemetrs)
{
	this->trafficSimulator_
		.GenerateChargingStations(numStations, simulationPararemetrs);
}


/// <summary>
/// Resets all simulation information.
/// </summary>
void TimeTable::ResetSimulation()
{
	// Reset time and number of run down vehicles.
	this->ActualTime_ = 0;
	this->CarsBatteryRunDown_ = 0;
	this->RunDownCharging = 0;
	this->NumCharging_ = 0;
	this->NumNotFinishedCars = 0;

	// Delete all scheduled actions and clear all info about the simulation.
	this->allEvents_ = 
		std::priority_queue<TableEvent, std::vector<TableEvent>, MyGreater>();
	this->travelDurations_.clear();
	this->batteryDifferences_.clear();
	this->chargingLevels_.clear();
	this->waitingTimes_.clear();
	this->runDownPositions_.clear();
	this->numReturnedVehicles_ = 0;

	// Reset all more low level simulation properties.
	this->trafficSimulator_.ResetSimulation();
}


/// <summary>
/// Executes the simulation.
/// </summary>
/// <param name="numStations">Number of charging stations.</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
void TimeTable::RunSimulation(
	int32_t numStations, SimulationParameters& simulationParameters, bool logs)
{
	// Initialize charging station information containers.
	this->chargingLevels_ = std::vector<std::vector<double>>(numStations);
	this->waitingTimes_ = std::vector<std::pair<int32_t, double>>(numStations);

	// Generate first two cars.
	this->AddEvent(TableEvent(this->trafficSimulator_.GenerateCar(
		this->ActualTime_, simulationParameters), 0, Actions::RUNNING));
	
	this->AddEvent(TableEvent(this->trafficSimulator_.GenerateCar(
		this->ActualTime_, simulationParameters), 
		this->trafficSimulator_.GenerateNextDepartureTime(),
		Actions::START));


	// Run the rest of the simulation.
	this->RunRestSimulation(simulationParameters, logs);
}



/// <summary>
/// Adds event to the timetable.
/// </summary>
/// <param name="tableEvent">Event to be added.</param>
void TimeTable::AddEvent(TableEvent tableEvent)
{
	this->allEvents_.push(tableEvent);
}


/// <summary>
/// </summary>
/// <returns>Returns const reference to next event to be done.</returns>
const TableEvent& TimeTable::GetNextEvent()
{
	return this->allEvents_.top();
}


/// <summary>
/// Executes simulation after generation of first 2 cars.
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
void TimeTable::RunRestSimulation(
	SimulationParameters& simulationParameters, bool logs)
{
	double lastTime = 0;

	while (this->allEvents_.top().ActionTime_ < simulationParameters.SimulationTime)
	// Run simulation until simulation time exceeds.
	{
		// Get next event.
		const TableEvent& nextEvent = this->allEvents_.top();
		this->ActualTime_ = nextEvent.ActionTime_;

		if (this->ActualTime_ > lastTime + 100)
		// After 10 simulated minutes write info to stdout.
		{
			lastTime += 100;
			std::cout << "Time " << lastTime  << " elapsed!" << std::endl;
		}

		switch (nextEvent.Action_)
		// Choose the appropriate action.
		{
			case Actions::START:
				// Vehicle is starting.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: START " << "Car ID : " << 
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_START(simulationParameters);
				break;

			case Actions::RUNNING:
				// Vehicle is running.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: RUNNING " << "Car ID : " << 
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_RUNNIG(simulationParameters, logs);
				break;

			case Actions::CHARGING:
				// Vehicle is charged in the charging staion.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: CHARGING " << "Car ID : " <<
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_CHARGING();
				break;

			case Actions::END_CHARGING:
				// Vehicle ended charging.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: END_CHARGING " << "Car ID : " << 
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_END_CHARGING(simulationParameters);
				break;

			case Actions::WAITING_LINE:
				// Vehicle is waiting in the line for charging.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: WAITING_LINE " << "Car ID : " <<
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_WAITING_LINE();
				break;

			case Actions::WAITING_RETURN:
				// Vehicle is waiting for the returning to start.
				if (logs)
				// Write info to stdout.
				{
					std::cout << "Action time: " << this->ActualTime_ << 
						", action type: WAITING_RETURN " << "Car ID : " <<
						nextEvent.CarID_ << std::endl;
				}
				this->reaction_WAITING_RETURN(logs);
				break;

			default:
				break;
		}

		// No next event (probably error)
		if (this->allEvents_.empty())
			return;
	}

	this->NumNotFinishedCars = this->trafficSimulator_.GetNumCars();
}


/// <summary>
/// </summary>
/// <returns>Returns reference to vector of all finished travel 
/// durations of the simulation.</returns>
std::vector<double>& TimeTable::GetTravelDurations()
{
	return this->travelDurations_;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to vector of all battery 
/// differences between start and end of the route.</returns>
std::vector<double>& TimeTable::GetBatteryDifferences()
{
	return this->batteryDifferences_;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to vector of all charging level for
/// each charging station.</returns>
std::vector<std::vector<double>>& TimeTable::GetChargingLevels()
{
	return this->chargingLevels_;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to vector of number of customers and 
/// sum of their waiting times for each charging station.</returns>
std::vector<std::pair<int32_t, double>>& TimeTable::GetWaitingTimes()
{
	return this->waitingTimes_;
}


/// <summary>
/// </summary>
/// <returns>Returns map with key edge identifier and value reference to vector of 
/// battery info for each segment of the corresponding edge.</returns>
std::map<edge_t, std::vector<BateryPair>&> TimeTable::GetAllMapSegmentsInfo()
{
	return this->trafficSimulator_.GetMap().GetAllSegments();
}


/// <summary>
/// </summary>
/// <returns>Returns reference to `TrafficSimulator` object.</returns>
TrafficSimulator& TimeTable::GetTrafficSimulator()
{
	return this->trafficSimulator_;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to all battery run down positions 
/// of each city.</returns>
std::map<int, std::vector<MapPosition>>& TimeTable::GetRunDownPositions()
{
	return this->runDownPositions_;
}


int32_t TimeTable::GetNumReturnedVehicles()
{
	return this->numReturnedVehicles_;
}



/// <summary>
/// Changes action of the car to RUNNING. 
/// </summary>
void TimeTable::reaction_START(SimulationParameters& simulationParameters)
{
	// Pop action before starting new (to avoid problems with priority queue).
	int64_t carID = this->allEvents_.top().CarID_;
	this->allEvents_.pop();

	// Schedule next car generation.
	this->AddEvent(
		TableEvent(this->trafficSimulator_
			.GenerateCar(this->ActualTime_, simulationParameters),
			this->ActualTime_ + this->trafficSimulator_.GenerateNextDepartureTime(),
			Actions::START));

	// Start running with the vehicle.
	this->AddEvent(TableEvent(carID, this->ActualTime_, Actions::RUNNING));
}



/// <summary>
/// Properly reacts to the action for crossing next edge (crosses the edge and 
/// moves to another or enters charging station or gets to the end).
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
void TimeTable::reaction_RUNNIG(SimulationParameters& simulationParameters, bool logs)
{
	// Get current car.
	int64_t carID = this->allEvents_.top().CarID_;
	CarIterator carIterator = this->trafficSimulator_.GetCarIterator(carID);

	Car& car = (*carIterator).second;

	// Remove current action from the actions container.
	this->allEvents_.pop();
	

	double transitionTime;
	if (car.IsInsideEdge())
	// Starting the path -> just get to the closest node.
	{
		if (logs)
		// Write info to the stdout.
		{
			std::cout << "Is STARTING" << std::endl;
		}

		// Get transition time.
		transitionTime = car.GetNextNodeTransitTime(
			this->trafficSimulator_.GetMap(),
			car.GetPosition().EdgeSegmentID,
			car.IsIncreasing(car.GetPosition().FurtherVertexID,
				car.GetPosition().CloserVertexID), 
			-1);


		if (!car.MoveFirstSegment(this->trafficSimulator_.GetMap(), logs))
		// Car battery run down -> store info and delete the car.
		{
			if (car.GoingToStation())
				this->RunDownCharging++;

			this->CarsBatteryRunDown_++;
			this->addRunDownCar(car);
			this->trafficSimulator_.DeleteCar((*carIterator).first);
			return;
		}

		// Add new event for crossing the next edge.
		this->AddEvent(
			TableEvent(carID, this->ActualTime_ + transitionTime, Actions::RUNNING));
	}

	else if (car.IsFinishing())
	// Car is either finishing the route or entering the charging station.
	{
		if (logs)
		// Write info to the stdout.
		{
			std::cout << "Is FINISHING" << std::endl;
		}


		// Find the ending position (end of the path or charging station).
		MapPosition endPosition = car.GetEndPosition();
		if (car.GoingToStation())
		// Car entering the charging station.
		{
			endPosition = car.GetStationPosition();
		}


		if (car.GetPosition().CloserVertexID != endPosition.CloserVertexID)
		// If car is in the further vertex of the final edge 
		// -> just change further to closer vertex to get to the end 
		// from the opposite side (to avoid moving to the other end of 
		// the edge and then returning to the middle of the edge).
		{
			endPosition.ChangeParameters(endPosition.FurtherVertexID,
				endPosition.CloserVertexID,
				endPosition.EdgeSegmentID,
				this->trafficSimulator_.GetMap().GetGraph());
		}


		// Determine whether moving inside the edge or 
		// starting from one end of the edge.
		int secondSegment = -1;
		const MapPosition& currentPosition = car.GetPosition();
		if (currentPosition.EdgeSegmentID != 0 &&
			currentPosition.EdgeSegmentID != this->trafficSimulator_.GetMap()
			.GetGraph()[currentPosition.EdgeID].GetNumSegments() - 1)
		// Movement only inside the edge (not on the end).
		{
			secondSegment = currentPosition.EdgeSegmentID;
		}


		// Get transition time.
		transitionTime = car.GetNextNodeTransitTime(
			this->trafficSimulator_.GetMap(),
			endPosition.EdgeSegmentID,
			// Moving to the middle of the edge -> change order 
			// of the start and end for the `IsIncreasing()` function.
			car.IsIncreasing(endPosition.FurtherVertexID,
				endPosition.CloserVertexID),
			secondSegment);


		if (!car.MoveToFinalSegment(this->trafficSimulator_.GetMap(), logs))
		// Car battery run down -> store info and delete the car.
		{
			if (car.GoingToStation())
				this->RunDownCharging++;
			
			this->CarsBatteryRunDown_++;
			this->addRunDownCar(car);
			this->trafficSimulator_.DeleteCar((*carIterator).first);
			return;
		}

		// Change the event of the car to either waiting for the charging or waiting in the end.
		if (car.GoingToStation())
		// Car waiting for charging.
		{
			car.SetStartLineWaitingTime(this->ActualTime_);
			this->AddEvent(
				TableEvent(carID, this->ActualTime_ + transitionTime,
					Actions::WAITING_LINE));
		}

		else
		// Car waiting in the end of the road.
		{
			this->AddEvent(
				TableEvent(carID, this->ActualTime_ + transitionTime,
					Actions::WAITING_RETURN));
		}
	}

	else
	// Car is just crossing the edge (doesn't end or entering the charging station).
	{
		if (logs)
		// Write info to the stdout.
		{
			std::cout << "Is CROSSING" << std::endl;
		}

		// We want to pass throught the whole edge (doesn't matter from which side 
		//of the edge we start to determine transition time).
		transitionTime = car.GetNextNodeTransitTime(
			this->trafficSimulator_.GetMap(), 0, true, -1);


		if (!car.MoveToNextNode(this->trafficSimulator_.GetMap(), logs))
		// Car battery run down -> store info and delete the car.
		{
			if (car.GoingToStation())
				this->RunDownCharging++;
			
			this->CarsBatteryRunDown_++;
			this->addRunDownCar(car);
			this->trafficSimulator_.DeleteCar((*carIterator).first);
			return;
		}


		if (!car.GoingToStation())
		// Car is not going to the charging station 
		// -> decide whether it should start going towards it.
		{
			if (!car.TriedCharging())
			// Car already hasn't tried to go to charging station 
			// -> it is possible to try to go to the station
			{

				// Check if the battery level is high enough to 
				// not consider charging at all.
				double batteryTreshold = 
					this->trafficSimulator_.GenerateBatteryTreshold();
				if (car.BatteryLevelToGoCharging(
					simulationParameters.ChargingTreshold, 
					simulationParameters.NotChargingTreshold, batteryTreshold))
				// Level of the battery is low enough to plan path 
				// considering visiting charging station -> update the path
				{
					// Vehicle tried to go charging (won't try it anymore).
					car.SetChargingTry(true);

					if (car.GoCharging(batteryTreshold, simulationParameters))
					// Battery level low enough to go charging 
					// -> update vehicle path which will visit the charging station
					{
						car.UpdateVehiclePath(
							car.GetPosition().FurtherVertexID,
							car.GetEndPosition().CloserVertexID,
							this->trafficSimulator_.GetMap(),
							batteryTreshold,
							true,
							simulationParameters);
					}
				}
			}
		}


		// Add new event for next edge crossing.
		this->AddEvent(TableEvent(carID, this->ActualTime_ + transitionTime, 
			Actions::RUNNING));
	}
}


/// <summary>
/// Applies charging on the vehicle and plans its end.
/// </summary>
void TimeTable::reaction_CHARGING()
{
	// Get car object.
	int64_t carID = this->allEvents_.top().CarID_;
	CarIterator carIterator = this->trafficSimulator_.GetCarIterator(carID);

	// Remove current action from the timetable.
	this->allEvents_.pop();

	// Already used percentage of the battery.
	double emptyBatteryLevel = 1 - (*carIterator).second.GetBatteryLevel();


	// Do proper actions in the charging station.
	std::unique_ptr<ChargingStation>& chargingStation = 
		this->trafficSimulator_.GetMap()
		.GetChargingStation((*carIterator).second.GetChargingStationID());

	// Store info about battery level.
	this->chargingLevels_[(*chargingStation).StationID_]
		.push_back(emptyBatteryLevel);

	// Compute waiting time for charging the battery to full capacity.
	double batteryChargingTime = 
		(*chargingStation).GetChargingTime(emptyBatteryLevel);


	double waitingTime = this->ActualTime_ - (*carIterator).second.GetStartLineWaitingTime();

	// Update info about waiting times in the station.
	this->waitingTimes_[(*chargingStation).StationID_].first++;
	this->waitingTimes_[(*chargingStation).StationID_].second += 
		waitingTime;//this->ActualTime_ - (*carIterator).second.GetStartLineWaitingTime();

	// if ((*carIterator).second.GetExpectedPathBatteryLevel(
	// 	(waitingTime * (*carIterator).second.GetVelocity())) <= 0)
	// // Waiting too long -> let the vehicle run down of the battery 
	// // (to ensure that all vehicles won't end in the same station).
	// {
	// 	batteryChargingTime = 0;
	// 	emptyBatteryLevel = -1;
	// 	// std::cout << "Long waiting" << std::endl;
	// }

	// By default charge to full capacity or null the battery level (if waiting too long).
	(*carIterator).second.ChargeBattery(emptyBatteryLevel);

	// Plan end of the charging.
	this->AddEvent(TableEvent(carID,
		this->ActualTime_ + batteryChargingTime,
		Actions::END_CHARGING));
}


/// <summary>
/// Ends charging of the car and starts moving on the rest of the path.
/// Also moves next car in the waiting line to the charging station.
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
void TimeTable::reaction_END_CHARGING(SimulationParameters& simulationParameters)
{
	// Get current car.
	int64_t carID = this->allEvents_.top().CarID_;
	CarIterator carIterator = this->trafficSimulator_.GetCarIterator(carID);
	Car& car = (*carIterator).second;

	// Removes action from the timetable.
	this->allEvents_.pop();

	
	// Current charging station.
	std::unique_ptr<ChargingStation>& chargingStation = this->trafficSimulator_.GetMap()
		.GetChargingStation(car.GetChargingStationID());

	// Leave the charging station.
	(*chargingStation).RemoveCustomer();


	// Remove customer from the waiting line and plan his start charging time.
	CarIterator nextChargedCarIterator =
		this->trafficSimulator_.GetCarIterator((*chargingStation).NextCustomer());

	if (nextChargedCarIterator != this->trafficSimulator_.GetCarIterator(-1))
	// If some next car in the charging station queue exists -> plan its charging time. 
	{
		// this->waitingTimes_[(*chargingStation).StationID_].first++;

		if ((*nextChargedCarIterator).first == carID)
		// Error car ID check.
		{
			std::cout << "Same next customer ID and current car" << std::endl;
			throw "Same next customer ID and current car";
		}
		
		// Plan charging of the next vehicle.
		this->AddEvent(TableEvent((*nextChargedCarIterator).first,
			this->ActualTime_,
			Actions::CHARGING));
	}

	// Update car path (no battery treshold 
	// (to not consider visiting charging station)).
	car.UpdateVehiclePath(car.GetPosition().CloserVertexID,
		car.GetEndPosition().CloserVertexID,
		this->trafficSimulator_.GetMap(),
		0, 
		false,
		simulationParameters);

	// Start running with the current vehicle.
	this->AddEvent(TableEvent(carID, this->ActualTime_, Actions::RUNNING));
}


/// <summary>
/// Do proper actions during waiting in the charging station line (either waits
/// for the free slot or starts charging (when capacity is not full)).
/// </summary>
void TimeTable::reaction_WAITING_LINE()
{
	// Get current car.
	int64_t carID = this->allEvents_.top().CarID_;
	CarIterator carIterator = this->trafficSimulator_.GetCarIterator(carID);

	// Remove event from the timetable.
	this->allEvents_.pop();


	// Get charging station.
	std::unique_ptr<ChargingStation>& chargingStation = 
		this->trafficSimulator_.GetMap()
		.GetChargingStation((*carIterator).second.GetChargingStationID());

	(*chargingStation).AddCustomer((*carIterator).first);

	if ((*chargingStation).NumCustomers_ <= (*chargingStation).Capacity_)
	// Go charging immediately (less customers than maximal capacity).
	{
		// Just add customer (waiting time was 0).
		// this->waitingTimes_[(*chargingStation).StationID_].first++;
		
		// Remove current car from the queue (current car is the only car in the queue,
		// because the capacity is higher than current number of customers).
		int next = (*chargingStation).NextCustomer();
		if (next != carID)
		// Validity check.
		{
			std::cout << "Different customer ID than current car" << std::endl;
			throw "Different customer ID than current car";
		}

		// Immediately plan charging.
		this->AddEvent(TableEvent(carID, this->ActualTime_, Actions::CHARGING));
	}
}


/// <summary>
/// Either wait for the returning or ends route of the vehicle.
/// </summary>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
void TimeTable::reaction_WAITING_RETURN(bool logs)
{
	this->numReturnedVehicles_++;

	// Get current car.
	int64_t carID = this->allEvents_.top().CarID_;
	CarIterator carIterator = this->trafficSimulator_.GetCarIterator(carID);
	Car& car = (*carIterator).second;

	// Remove event from the timetable.
	this->allEvents_.pop();


	// Save info about path duration.
	this->travelDurations_.push_back(this->ActualTime_ - car.GetStartTime());
	this->batteryDifferences_
		.push_back(car.GetStartBatteryLevel() - car.GetBatteryLevel());

	// std::cout << car.GetStartBatteryLevel() << std::endl;

	if (car.IsReturning())
	// Car returned from the finish to original start 
	// -> delete car (not further actions).
	{
		if (logs)
		// Write info to stdout.
		{
			std::cout << "Car RETURNED to START" << std::endl;
		}
		
		this->trafficSimulator_.DeleteCar((*carIterator).first);
		return;
	}

	// Car arrived to the end -> wait there randomly generated amount of time 
	// -> start returning to the start.
	car.StartReturning(this->ActualTime_ + car.GetWaitingTime());


	// Plan running.
	this->AddEvent(TableEvent(carID, this->ActualTime_ + car.GetWaitingTime(),
		Actions::RUNNING));
}


/// <summary>
/// Adds `car` position to the container of all positions where some car ran 
/// down with battery (information for optimization).
/// </summary>
/// <param name="car">Reference to car object which ran down.</param>
void TimeTable::addRunDownCar(Car& car)
{
	// Get city ID.
	int32_t runDownCityID = this->trafficSimulator_.GetMap()
		.GetGraph()[car.GetPosition().CloserVertexID].CityID_;

	// Save run down position.
	auto cityIterator = this->runDownPositions_.find(runDownCityID);
	if (cityIterator == this->runDownPositions_.end())
	// First car run down in the city -> create proper container
	{
		std::vector<MapPosition> runDownVector{ car.GetPosition() };
		this->runDownPositions_.insert(
			std::make_pair(runDownCityID, std::move(runDownVector)));
	}
	else
	// There are already run down cars in the city -> just add the car there
	{
		cityIterator->second.push_back(car.GetPosition());
	}
}