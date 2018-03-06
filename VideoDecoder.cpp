#include "VideoDecoder.h"

//----------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------
VideoDecoder::VideoDecoder(const SyncDetector& syncDetector, const FrameBuilder& frameBuilder, const Logger& logger)
:_syncDetector(syncDetector), _frameBuilder(frameBuilder), _logger(logger)
{
	blankingPercentage = 0.1;
	blankingLeadingPercentage = 0.95;

	lineWidthInPixels = 930;
	maxChunkSizeInBytes = 1024*1024*1024;
}
