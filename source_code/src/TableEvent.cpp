#include "optimizer/TableEvent.hpp"


/// <summary>
/// Initializes the event object.
/// </summary>
/// <param name="carID">ID of the vehicle of the event.</param>
/// <param name="actionTime">Time of the event.</param>
/// <param name="action">Type of the action.</param>
TableEvent::TableEvent(int64_t carID, double actionTime, Actions action)
	: CarID_(std::move(carID)), ActionTime_(std::move(actionTime)),
	Action_(std::move(action))
{ }