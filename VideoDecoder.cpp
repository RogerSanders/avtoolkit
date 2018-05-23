#include "VideoDecoder.h"
#include "StatisticsHelpers.h"

//----------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------
VideoDecoder::VideoDecoder(const SyncDetector& syncDetector, const FrameBuilder& frameBuilder, const Logger& log)
:_syncDetector(syncDetector), _frameBuilder(frameBuilder), _log(log)
{
	//##TODO## Comment these
	blankingPercentage = 0.1;
	blankingLeadingPercentage = 0.95;
	burstWaveSkipCount = 2;
	useAverageFieldBlankingLevel = true;
	useAverageFieldBurstWaveFrequency = true;
	matchColorBurstPhaseBetweenLines = true;
	rawOutputOnly = false;
	decodeColor = true;
	useIRE7Point5 = true;
	decodeAsYUV = true;
	forceMonoOutputForColorDecoding = false;
	forceColorBurstWaveFrequency = false;
	waveFrequencyTolerance = 0.001;

	lineWidthInPixels = 930;
	maxChunkSizeInBytes = 1024*1024*1024;
}

//----------------------------------------------------------------------------------------
// Video conversion methods
//----------------------------------------------------------------------------------------
bool VideoDecoder::GetWavePeakPositionsForLine(const FrameBuilder::LineInfo& lineInfo, std::vector<double>& wavePeakPositions, bool& firstWavePeakIsPositive) const
{
	// Initialize our output values
	firstWavePeakIsPositive = false;
	wavePeakPositions.clear();

	// If there aren't enough colour burst wave entries in the target line to calculate any peak positions, return
	// false.
	if (lineInfo.colorBurstWaves.size() <= (burstWaveSkipCount * 2))
	{
		return false;
	}

	// Calculate the positions of each peak of the color burst wave. Note that it doesn't matter here if we have zero
	// level wrong, and our positive and negative wave oscillations are different lengths as a result. The peak
	// positions of each oscillation will be calculated correctly despite this issue, since the waves will still be
	// symmetrical. Also note that we potentially skip some wave oscillations at the beginning and end of the detected
	// colour burst wave set. We do this as it has been observed that the last few waves at either end are often
	// distorted, being abnormally shorter or longer than the correct period.
	wavePeakPositions.reserve(lineInfo.colorBurstWaves.size() - (burstWaveSkipCount * 2));
	for (size_t waveIndex = burstWaveSkipCount; waveIndex < (lineInfo.colorBurstWaves.size() - burstWaveSkipCount); ++waveIndex)
	{
		const FrameBuilder::ColorBurstWaveInfo& waveInfo = lineInfo.colorBurstWaves[waveIndex];
		if (wavePeakPositions.empty())
		{
			firstWavePeakIsPositive = waveInfo.isPositive;
		}
		double wavePeakPos = waveInfo.startPos + ((waveInfo.endPos - waveInfo.startPos) / 2);
		wavePeakPositions.push_back(wavePeakPos);
	}
	return true;
}

//----------------------------------------------------------------------------------------
void VideoDecoder::CalculateIRELevelsForField(const FrameBuilder::FieldInfo& fieldInfo, double& ireLevel0ForField, double& ireLevel100ForField) const
{
	// Calculate the average sync and blanking levels across the field
	double averageSyncLevelForField = 0.0;
	double averageBlankingLevelForField = 0.0;
	for (const FrameBuilder::LineInfo& lineInfo : fieldInfo.lines)
	{
		averageSyncLevelForField += lineInfo.averageSyncLevel;
		averageBlankingLevelForField += lineInfo.averageBlankingLevel;
	}
	averageSyncLevelForField /= (double)fieldInfo.lines.size();
	averageBlankingLevelForField /= (double)fieldInfo.lines.size();

	// Calculate the IRE levels for this line from the average sync and blanking levels. Sync occurs at IRE -40, and
	// blanking occurs at IRE 0. From these two reference levels, we can establish our IRE scale.
	//##FIX## Update this comment, and find a way to de-duplicate this calculation with FrameBuilder. It should
	//probably calculate the average levels and IRE scale for the fields for us.
	//##FIX## In a PAL signal, sync occurs at -43 IRE, and the colour burst has a 21.5 IRE amplitude with 10 cycles.
	ireLevel0ForField = averageBlankingLevelForField;
	ireLevel100ForField = averageBlankingLevelForField + ((averageBlankingLevelForField - averageSyncLevelForField) * (1.0 / 0.4));
}

//----------------------------------------------------------------------------------------
double VideoDecoder::CalculateColorBurstWaveFrequencyForLine(const FrameBuilder::LineInfo& lineInfo, std::vector<double>& wavePeakPositions) const
{
	// Calculate the frequency of the colour burst wave as the interval between the first and last wave peak
	// positions, divided by the number of wave oscillations between them.
	return (wavePeakPositions.back() - wavePeakPositions.front()) / (wavePeakPositions.size() - 1);
}

//----------------------------------------------------------------------------------------
double VideoDecoder::CalculateColorBurstWaveFrequencyForField(const FrameBuilder::FieldInfo& fieldInfo) const
{
	// Build a set of all calculated burst wave frequencies for each line in the field
	std::vector<double> burstWaveFrequenciesForField;
	burstWaveFrequenciesForField.reserve(fieldInfo.lines.size());
	for (const FrameBuilder::LineInfo& lineInfo : fieldInfo.lines)
	{
		std::vector<double> wavePeakPositions;
		bool firstWavePeakIsPositive;
		if ((lineInfo.leadingSyncInfo.type == SyncDetector::SyncType::Horizontal) && GetWavePeakPositionsForLine(lineInfo, wavePeakPositions, firstWavePeakIsPositive))
		{
			burstWaveFrequenciesForField.push_back(CalculateColorBurstWaveFrequencyForLine(lineInfo, wavePeakPositions));
		}
	}

	// Set the initial field burst wave frequency as the average value across all lines in the field
	double burstWaveFrequencyForField;
	burstWaveFrequencyForField = std::accumulate(burstWaveFrequenciesForField.cbegin(), burstWaveFrequenciesForField.cend(), 0.0) / (double)burstWaveFrequenciesForField.size();

	// Use the initial burst wave frequency we calculated from a simple average, and use it to build a set of wave
	// frequencies for lines in this frame which are close to the estimated frequency.
	std::vector<double> burstWaveFrequenciesAveraged;
	double waveFrequencyEffectiveTolerance = burstWaveFrequencyForField * waveFrequencyTolerance;
	for (size_t lineNo = 0; lineNo < (fieldInfo.lines.size() - 1); ++lineNo)
	{
		// Calculate the positions of each peak of the colour burst wave for this line
		const FrameBuilder::LineInfo& lineInfo = fieldInfo.lines[lineNo];
		std::vector<double> wavePeakPositions;
		bool firstWavePeakIsPositive;
		if ((lineInfo.leadingSyncInfo.type != SyncDetector::SyncType::Horizontal) || !GetWavePeakPositionsForLine(lineInfo, wavePeakPositions, firstWavePeakIsPositive))
		{
			continue;
		}

		// Calculate the positions of each peak of the colour burst wave in the following line
		const FrameBuilder::LineInfo& nextLineInfo = fieldInfo.lines[lineNo+1];
		std::vector<double> nextLineWavePeakPositions;
		bool nextLineFirstWavePeakIsPositive;
		if (!GetWavePeakPositionsForLine(nextLineInfo, nextLineWavePeakPositions, nextLineFirstWavePeakIsPositive))
		{
			continue;
		}

		// Calculate the best phase lock between the current line and the following line which most closely matches the
		// estimated colour burst wave frequency
		size_t bestStartingPeakPosition;
		bool bestWavePeakPositionIsPositive;
		double bestDisplacementValue;
		PhaseLockColorBurstSamples(wavePeakPositions, firstWavePeakIsPositive, nextLineWavePeakPositions, nextLineFirstWavePeakIsPositive, burstWaveFrequencyForField, bestStartingPeakPosition, bestWavePeakPositionIsPositive, bestDisplacementValue);
		double burstWaveFrequency = burstWaveFrequencyForField + bestDisplacementValue;

		// If the phase lock is outside our tolerance, reject it, and advance to the next line. We do this because we
		// don't want samples which are wildly off from biasing our frequency calculations. We exclude these cases here,
		// and focus on the ones which are close to our expectations to refine the colour burst frequency.
		if (bestDisplacementValue > waveFrequencyEffectiveTolerance)
		{
			continue;
		}

		// Add the calculated frequency to the set of burst wave frequencies for this line
		burstWaveFrequenciesAveraged.push_back(burstWaveFrequency);
	}

	// If we rejected too many line samples, fall back to the initial simple average as the burst wave frequency for the
	// field, and return it to the caller.
	if (burstWaveFrequenciesAveraged.size() < (fieldInfo.lines.size() / 10))
	{
		//##DEBUG##
		std::wcout << "Line tolerance check failed!\n";
		return burstWaveFrequencyForField;
	}

	// Set the field burst wave frequency as the median value of the built set of line burst wave frequencies, and
	// return it to the caller.
	//##TODO## Decode between these two approaches
	burstWaveFrequencyForField = FindMedianValue(burstWaveFrequenciesAveraged.begin(), burstWaveFrequenciesAveraged.end());
	//burstWaveFrequencyForField = std::accumulate(burstWaveFrequenciesAveraged.cbegin(), burstWaveFrequenciesAveraged.cend(), 0.0) / (double)burstWaveFrequenciesAveraged.size();
	return burstWaveFrequencyForField;
}

//----------------------------------------------------------------------------------------
void VideoDecoder::PhaseLockColorBurstSamples(std::vector<double> wavePeakPositions, bool firstWavePeakIsPositive, std::vector<double> targetWavePeakPositions, bool firstTargetWavePeakIsPositive, double burstWaveFrequency, size_t& bestMatchPeakPositionIndex, bool& bestMatchPeakPositionIsPositive, double& phaseLockWaveFrequencyDisplacement) const
{
	// Initialize our running values for tracking the best match found so far
	size_t bestSourceWavePeakPosition = 0;
	bool bestWavePeakPositionIsPositive = firstWavePeakIsPositive;
	double bestDisplacementValue = burstWaveFrequency;
	double bestEffectiveFrequency = burstWaveFrequency;

	// Compare all source wave peak positions with all target wave peak positions, and identify the source wave peak
	// position which achieves the best phase lock with a target wave peak, based on the relative positions of each wave
	// and the provided wave frequency.
	size_t sourceWavePeakPositionIndex = 0;
	while (sourceWavePeakPositionIndex < wavePeakPositions.size())
	{
		// Retrieve the next source wave peak position
		double sourceWavePeak = wavePeakPositions[sourceWavePeakPositionIndex];

		// Compare the current source wave peak position with each target wave peak position
		size_t targetWavePeakPositionIndex = 0;
		while (targetWavePeakPositionIndex < targetWavePeakPositions.size())
		{
			// Retrieve the next target wave peak position
			double targetWavePeak = targetWavePeakPositions[targetWavePeakPositionIndex];

			// Calculate the wave oscillation count between the source and target peak positions, and the displacement
			// from the ideal phase lock between them.
			double distanceBetweenPeaks = std::abs(targetWavePeak - sourceWavePeak);
			double waveBurstsBetweenPeaks = distanceBetweenPeaks / burstWaveFrequency;
			unsigned int waveBurstsBetweenPeaksRounded = (unsigned int)(waveBurstsBetweenPeaks + 0.5);
			double displacementFromIdealIntervalSigned = waveBurstsBetweenPeaks - (double)waveBurstsBetweenPeaksRounded;
			double displacementFromIdealInterval = std::abs(displacementFromIdealIntervalSigned);
			double effectiveFrequency = distanceBetweenPeaks / (double)waveBurstsBetweenPeaksRounded;

			// If the polarity of the matched wave oscillation is inverted, adjust the displacement to be relative to
			// the adjacent oscillation with the correct polarity.
			bool expectedMatchPeakIsPositive = ((sourceWavePeakPositionIndex + waveBurstsBetweenPeaksRounded) & 1) == (firstWavePeakIsPositive ? 0 : 1);
			bool matchPeakIsPositive = (targetWavePeakPositionIndex & 1) == (firstTargetWavePeakIsPositive ? 0 : 1);
			if (expectedMatchPeakIsPositive != matchPeakIsPositive)
			{
				displacementFromIdealInterval = burstWaveFrequency - displacementFromIdealInterval;
				effectiveFrequency = (distanceBetweenPeaks - burstWaveFrequency) / (double)waveBurstsBetweenPeaksRounded;
			}

			// If this match is the best one so far, save information on it.
			if (displacementFromIdealInterval <= bestDisplacementValue)
			{
				bestSourceWavePeakPosition = sourceWavePeakPositionIndex;
				bestDisplacementValue = displacementFromIdealInterval;
				bestEffectiveFrequency = effectiveFrequency;
				bestWavePeakPositionIsPositive = (sourceWavePeakPositionIndex & 1) == (firstWavePeakIsPositive ? 0 : 1);
			}

			// Advance to the next wave peak in the target set of peak positions
			++targetWavePeakPositionIndex;
		}

		// Advance to the next wave peak in the source set of peak positions
		++sourceWavePeakPositionIndex;
	}

	// Return the best identified match to the caller
	bestMatchPeakPositionIndex = bestSourceWavePeakPosition;
	bestMatchPeakPositionIsPositive = bestWavePeakPositionIsPositive;
	phaseLockWaveFrequencyDisplacement = bestEffectiveFrequency - burstWaveFrequency;
}

//----------------------------------------------------------------------------------------
void VideoDecoder::ConvertYIQToRGB(double sampleY, double sampleI, double sampleQ, double& red, double& green, double& blue)
{
	// YIQ to RGB
	// http://dystopiancode.blogspot.com.au/2012/06/yiq-rgb-conversion-algorithms.html
	// https://en.wikipedia.org/wiki/YIQ#From_YIQ_to_RGB
	//##TODO## Confirm the correct values here. We've seen different formulas online.
	double sampleIScaled = sampleI * 0.5957;
	double sampleQScaled = sampleQ * 0.5226;
	red   = sampleY + ( 0.9563*sampleIScaled) + ( 0.6210*sampleQScaled);
	green = sampleY + (-0.2721*sampleIScaled) + (-0.6474*sampleQScaled);
	blue  = sampleY + (-1.1070*sampleIScaled) + ( 1.7046*sampleQScaled);
}

//----------------------------------------------------------------------------------------
void VideoDecoder::ConvertYUVToRGB(double sampleY, double sampleU, double sampleV, double& red, double& green, double& blue)
{
	// YUV to RGB
	// http://softpixel.com/~cwright/programming/colorspace/yuv/
	// https://en.wikipedia.org/wiki/YUV#Conversion_to/from_RGB
	//##TODO## Confirm the correct values here. We've seen different formulas online.
	double sampleUScaled = sampleU * 0.5;
	double sampleVScaled = sampleV * 0.5;
	//red   = sampleY + ( 1.4075 * sampleVScaled);
	//green = sampleY + (-0.3455 * sampleUScaled) + (-0.7169 * sampleVScaled);
	//blue  = sampleY + ( 1.7790 * sampleUScaled);
	red   = sampleY + ( 1.13983 * sampleVScaled);
	green = sampleY + (-0.39465 * sampleUScaled) + (-0.58060 * sampleVScaled);
	blue  = sampleY + ( 2.03211 * sampleUScaled);
}
