#include "FrameConverter.h"

//----------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------
FrameConverter::FrameConverter(const Logger& logger)
:_logger(logger)
{
	blankingPercentage = 0.1;
	blankingLeadingPercentage = 0.95;
}
