#include "FrameConverter.h"
#include <cmath>
#include <fstream>

//----------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------
FrameConverter::FrameConverter()
{
	//##FIX## We've seen as low as 0.125 for the Mega Drive here. Perhaps we should be using the front porch, and
	//relying on averaging to cancel out the colour burst? This doesn't work with a quantum signal though, as unlike an
	//analog signal we've got sampling bias potentially causing us to have a much higher or lower average sample value.
	//We'd have to fit the sample points to a curve and take a difference of the volumes below and above in order to
	//avoid this issue. Then again, higher sample rates avoid this problem, and we'll probably get a better result from
	//a longer averaged run than an often noisy and short trailing run from the last line.
	//double blankingSampleRatioToSyncWidth = 0.2;

	//##TODO## Comment these constants
	blankingPercentage = 0.1;
	blankingLeadingPercentage = 0.95;
	syncSearchPosIncrement = 0.0001;
	syncAmplitudeMinTolerance = 0.125;
	slopeValueFlatToleranceAsPercentage = 0.3;
}

//----------------------------------------------------------------------------------------
// Line analysis methods
//----------------------------------------------------------------------------------------
//##TODO## This function does a very good job, but not as good as the colour burst phase offset detection approach used
//in ld-decode. We should consider that approach the primary method, and use our approach here as a secondary method.
//##TODO## Then again, it's questionable locking to colour burst is actually the right approach. A real monitor locks to
//the sync signal, and with a monochrome signal there might be no colour burst at all to lock to. Perform more
//comparisons which focus on the actual image region, and ignore the synchronization of the colour burst between lines,
//to determine if the colour burst appears to be a more or less reliable sync point.
void FrameConverter::FindPreciseLineStartEndPos(const SyncDetector::SyncInfo& leadingSyncInfo, const SyncDetector::SyncInfo& followingSyncInfo, const tk::spline& lineSpline, double& preciseLineStartPos, double& preciseLineEndPos) const
{
	//##FIX## Make constants like these configurable
	//##DEBUG## Pretty good for sonic 2, not as good when large min/max range
	//double syncSearchPosIncrement = 0.0001;
	//double syncAmplitudeMinTolerance = 0.2;
	//double slopeValueFlatToleranceAsPercentage = 0.3;

	size_t lineSplineSampleRange = followingSyncInfo.endSampleNo - leadingSyncInfo.startSampleNo;

	// Find the point at which we cross the minimum threshold for ending the leading sync run
	double leadingSyncEndPosInSpline = (double)(leadingSyncInfo.endSampleNo - leadingSyncInfo.startSampleNo) / (double)lineSplineSampleRange;
	double leadingSyncEndMinimumThreshold = leadingSyncInfo.averageSyncLevel + (leadingSyncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance);
	double leadingSyncEndSearchPosInSpline = leadingSyncEndPosInSpline;
	while (lineSpline(leadingSyncEndSearchPosInSpline) < leadingSyncEndMinimumThreshold)
	{
		leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
	}

	// Find the point at which the leading sync run levels out to an acceptably flat slope, and mark it as the precise
	// line start pos.
	double leadingSyncLastSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
	double leadingSyncCurrentSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	double leadingSyncSlopeValueFlatTolerance = leadingSyncInfo.approxMinMaxSampleRange * slopeValueFlatToleranceAsPercentage;
	while (std::fabs(leadingSyncCurrentSlopeValue - leadingSyncLastSlopeValue) > leadingSyncSlopeValueFlatTolerance)
	{
		leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
		leadingSyncLastSlopeValue = leadingSyncCurrentSlopeValue;
		leadingSyncCurrentSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	}
	preciseLineStartPos = leadingSyncEndSearchPosInSpline;

	// Find the point at which we cross the minimum threshold for starting the following sync run
	double followingSyncStartPosInSpline = (double)(followingSyncInfo.startSampleNo - leadingSyncInfo.startSampleNo) / (double)lineSplineSampleRange;
	double followingSyncStartMinimumThreshold = followingSyncInfo.averageSyncLevel + (followingSyncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance);
	double followingSyncStartSearchPosInSpline = followingSyncStartPosInSpline;
	while (lineSpline(followingSyncStartSearchPosInSpline) < followingSyncStartMinimumThreshold)
	{
		followingSyncStartSearchPosInSpline -= syncSearchPosIncrement;
	}

	// Find the point at which the following sync run levels out to an acceptably flat slope, and mark it as the precise
	// line end pos.
	double followingSyncLastSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	followingSyncStartSearchPosInSpline += syncSearchPosIncrement;
	double followingSyncCurrentSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	double followingSyncSlopeValueFlatTolerance = followingSyncInfo.approxMinMaxSampleRange * slopeValueFlatToleranceAsPercentage;
	while (std::fabs(followingSyncCurrentSlopeValue - followingSyncLastSlopeValue) > followingSyncSlopeValueFlatTolerance)
	{
		followingSyncStartSearchPosInSpline -= syncSearchPosIncrement;
		followingSyncLastSlopeValue = followingSyncCurrentSlopeValue;
		followingSyncCurrentSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	}
	preciseLineEndPos = followingSyncStartSearchPosInSpline;
}
