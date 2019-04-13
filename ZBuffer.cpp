//============================================================================
// @name        : ZBuffer.cpp
// @author      : Ward Gauderis
// @date        : 19/03/2019
// @version     : 
// @copyright   : Project Software Engineering - BA1 Informatica - Ward Gauderis - University of Antwerp
// @description : 
//============================================================================
#include "ZBuffer.h"
#include <limits>

ZBuffer::ZBuffer(const unsigned int width, const unsigned int height) {
	resize(height);
	for (auto &row: *this) {
		row.resize(width, std::numeric_limits<double>::infinity());
	}
}
