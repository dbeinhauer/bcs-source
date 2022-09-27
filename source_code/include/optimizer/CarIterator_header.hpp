#ifndef CARITERATOR_HEADER_H_
#define CARITERATOR_HEADER_H_

#include <map>

#include "Car.hpp"

// Map to find car object by its ID.
typedef std::map<uint64_t, Car> CarMap;
typedef CarMap::iterator CarIterator;

#endif // !CARITERATOR_HEADER_H_

