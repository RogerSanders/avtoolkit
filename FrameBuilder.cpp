#include "FrameBuilder.h"
//##DEBUG##
#include <iostream>

//----------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------
FrameBuilder::FrameBuilder(const Logger& logger)
:_logger(logger)
{
	// By default, we want to combine interlaced fields together into a single frame.
	combineInterlacedFields = true;

	// This should be true for NTSC, false for PAL and SECAM.
	interlaceHalfLineInFirstField = true;

	//##FIX## We've seen as low as 0.125 for the Mega Drive here. Perhaps we should be using the front porch, and
	//relying on averaging to cancel out the colour burst? This doesn't work with a quantum signal though, as unlike an
	//analog signal we've got sampling bias potentially causing us to have a much higher or lower average sample value.
	//We'd have to fit the sample points to a curve and take a difference of the volumes below and above in order to
	//avoid this issue. Then again, higher sample rates avoid this problem, and we'll probably get a better result from
	//a longer averaged run than an often noisy and short trailing run from the last line.
	//double blankingSampleRatioToSyncWidth = 0.2;

	//##TODO## Comment these constants
	syncAmplitudeMinTolerance = 0.125;

	slopeValueFlatTolerance = 0.003;

	// We need a value for this, as there appears to be no way to safely determine the length of back porch without it,
	// which we need to establish IRE levels and detect the colour burst. NTSC signal specs state 0.075H +/- 0.005H for
	// hsync from IRE 0 to 0, with 0.02H front porch length, and 0.145H min (y) from front porch end to back porch end.
	// This gives us a 0.8 min ratio outlined in the NTSC specs. We go a bit lower to allow for non-conformant signals.
	//##TODO## Cross-check with PAL specs
	syncLengthToBackPorchMinRatio = 0.7;

	blankingUpsampleRatio = 10.0;

	// Exactly 9 cycles are supposed to be required for the color burst according to NTSC specs, meaning 18 half
	// oscillations should be present. In practice, the number of valid oscillations which meet the minimum amplitude
	// requirements varies. We require a lower number by default to give some tolerance for non-conforming signals.
	colorBurstMinimumHalfOscillationCount = 12;

	// The colour burst amplitude is supposed to be IRE 40 according to NTSC specs. We use this value to search for
	// peaks from the blanking level, so a value of 0.2 would represent a 40 IRE peak. Real world signals have shown
	// considerably burst amplitude levels than this however, so we only look for a burst at half that level to give
	// some tolerance.
	colorBurstIREDetectionThreshold = 0.06;

	colorBurstLineSyncTolerance = 0.48;

	// When enabled, this option allows the colour burst to be used to error-correct horizontal sync, to provide better
	// alignment between lines in a field. This option is recommended, as it provides more accurate synchronization than
	// the rising edge of the horizontal sync pulse can provide alone.
	useColorBurstForLineSyncCorrection = true;

	//##DEBUG## Pretty good for sonic 2, not as good when large min/max range
	//double syncSearchPosIncrement = 0.0001;
	//double syncAmplitudeMinTolerance = 0.2;
	//double slopeValueFlatToleranceAsPercentage = 0.3;

	// This value provides the tolerance to use when determining whether a colour burst wave is valid during the repair
	// process. Wave entries are only considered valid if they fit within this tolerance of the average length based on
	// other wave entries detected on that line.
	colorBurstRepairWaveLengthTolerance = 0.1;

	// When searching for trigger locations on rising and falling edges, an interative search is necessary to locate the
	// approximate location where the threshold is crossed. This tolerance value affects how precisely we lock to these
	// sample positions.
	slopeDetectionTolerance = 0.001;
}

//----------------------------------------------------------------------------------------
// Frame detection methods
//----------------------------------------------------------------------------------------
std::vector<FrameBuilder::FrameInfo> FrameBuilder::DetectFrames(const std::list<FieldInfo>& fields) const
{
	// Combine the fields into frames
	std::vector<FrameInfo> frames;
	frames.reserve(fields.size());
	FrameInfo currentFrameInfo;
	int lastFieldLineCount;
	for (auto fieldIterator = fields.cbegin(); fieldIterator != fields.cend(); ++fieldIterator)
	{
		// Retrieve information about the current field
		const FieldInfo& fieldInfo = *fieldIterator;
		int fieldLineCount = fieldInfo.lineCount;

		// If the current frame already has a field, and that field isn't interlaced with the current field, add the
		// current frame to the list of complete frames, and start a new frame. Note that we assume here that two
		// adjacent fields, with the second field having one more line than the preceding field form an interlaced pair
		// of fields within a frame. This correctly detects the basic even/odd interlaced field pairs for traditional
		// analog video streams, however we don't validate correct sync timing by searching for half-length lines at the
		// start and end of the field with the extra line. Assuming this relationship between adjacent fields should be
		// safe for non-contrived input data.
		//##FIX## NTSC and PAL reverse these concepts. In NTSC, the half-line occurs halfway through the frame. In PAL
		//and SECAM, the half-line occurs at the end of the frame.
		if (!currentFrameInfo.fields.empty() && (!combineInterlacedFields || (currentFrameInfo.fields.size() >= 2) || (!interlaceHalfLineInFirstField && ((lastFieldLineCount + 1) != fieldLineCount)) || (interlaceHalfLineInFirstField && ((lastFieldLineCount - 1) != fieldLineCount))))
		{
			frames.push_back(std::move(currentFrameInfo));
			currentFrameInfo = FrameInfo();
		}

		// Update information about the last processed field
		lastFieldLineCount = fieldLineCount;

		// Add this field to the current frame
		currentFrameInfo.fields.push_back(fieldInfo);
	}
	if (!currentFrameInfo.fields.empty())
	{
		frames.push_back(std::move(currentFrameInfo));
	}

	// If we have an incomplete interlaced frame at the start or end of the frame list, trim them now.
	if (combineInterlacedFields && (frames.size() > 1))
	{
		// If only half of an interlaced field was found in the first frame, drop it.
		auto firstFrameIterator = frames.begin();
		const FrameInfo& firstFrameInfo = *firstFrameIterator;
		const FrameInfo& secondFrameInfo = *std::next(firstFrameIterator);
		if ((secondFrameInfo.fields.size() > firstFrameInfo.fields.size()) && (secondFrameInfo.fields.back().lineCount == firstFrameInfo.fields.back().lineCount))
		{
			frames.erase(firstFrameIterator);
		}

		// If only half of an interlaced field was found in the last frame, drop it.
		auto lastFrameIterator = frames.rbegin();
		const FrameInfo& lastFrameInfo = *lastFrameIterator;
		const FrameInfo& secondLastFrameInfo = *std::next(lastFrameIterator);
		if ((secondLastFrameInfo.fields.size() > lastFrameInfo.fields.size()) && (secondLastFrameInfo.fields.back().lineCount == lastFrameInfo.fields.back().lineCount))
		{
			frames.erase(std::next(lastFrameIterator).base());
		}
	}

	// Return the list of frames to the caller
	return frames;
}

//----------------------------------------------------------------------------------------
// Colour burst methods
//----------------------------------------------------------------------------------------
bool FrameBuilder::ValidateColorBurst(unsigned int minimumWaveCount, const std::vector<ColorBurstWaveInfo>& burstWaves) const
{
	// Ensure the detected burst waves appear to represent successive positive and negative oscillations, and calculate
	// the oscillation count.
	bool waveAppearsValid = true;
	unsigned int oscillationCount = 0;
	bool firstWave = true;
	bool lastWaveWasPositive;
	auto burstWaveIterator = burstWaves.cbegin();
	while (burstWaveIterator != burstWaves.cend())
	{
		// Ensure the current wave swings in the opposite direction from the previous one
		const auto& waveInfo = *burstWaveIterator;
		if (!firstWave && (lastWaveWasPositive == waveInfo.isPositive))
		{
			waveAppearsValid = false;
			break;
		}

		// Increment the oscillation count for the waveform
		++oscillationCount;

		// Advance to the next wave
		lastWaveWasPositive = waveInfo.isPositive;
		firstWave = false;
		++burstWaveIterator;
	}

	// If the wave doesn't represent valid successive positive and negative transitions of the required threshold,
	// return false. Note that we still return the detected waves in this case, so that the caller can error correct if
	// required.
	if (!waveAppearsValid || (burstWaves.size() < minimumWaveCount))
	{
		return false;
	}

	// If we appear to have found a valid colour burst, return true.
	return true;
}

//----------------------------------------------------------------------------------------
bool FrameBuilder::PerformColorBurstLineSyncCorrection(FieldInfo& fieldInfo) const
{
	// Calculate a reference position in the colour burst to use to perform fine line synchronization
	bool foundReferenceColorBurst = false;
	double referenceColorBurstStartPosLineRelative;
	double referenceColorBurstWaveStride;
	bool referenceColorBurstPositive;
	bool referenceColorBurstLineOdd;
	unsigned int referenceColorBurstIterator = 0;
	while (!foundReferenceColorBurst && (referenceColorBurstIterator < fieldInfo.lines.size()))
	{
		const LineInfo& lineInfo = fieldInfo.lines[referenceColorBurstIterator];
		if (lineInfo.colorBurstValid)
		{
			const ColorBurstWaveInfo& referenceWave = lineInfo.colorBurstWaves[lineInfo.colorBurstWaves.size() / 2];
			referenceColorBurstStartPosLineRelative = referenceWave.startPos - lineInfo.backPorchStartPos;
			referenceColorBurstPositive = referenceWave.isPositive;
			referenceColorBurstLineOdd = ((referenceColorBurstIterator % 2) != 0);
			referenceColorBurstWaveStride = referenceWave.endPos - referenceWave.startPos;
			foundReferenceColorBurst = true;
		}
		++referenceColorBurstIterator;
	}
	if (!foundReferenceColorBurst)
	{
		return false;
	}

	// If we found a reference colour burst sync point, use it to perform fine adjustment on the line position across
	// the field.
	double accumulatedSyncCorrection = 0.0;
	unsigned int syncCorrectionEntries = 0;
	double colorBurstLineSyncRange = referenceColorBurstWaveStride * colorBurstLineSyncTolerance;
	double undefinedSyncCorrectionValue = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> syncCorrectionValues(fieldInfo.lines.size(), undefinedSyncCorrectionValue);
	for (unsigned int i = 0; i < fieldInfo.lines.size(); ++i)
	{
		// Retrieve info for the current line
		LineInfo& lineInfo = fieldInfo.lines[i];
		if (!lineInfo.colorBurstValid)
		{
			continue;
		}
		bool lineIsOdd = ((i % 2) != 0);

		// Search for a colour burst wave to phase lock with the reference wave
		size_t colorBurstWaveSearchIndex = 0;
		bool foundSyncCorrection = false;
		double syncCorrection;
		while (colorBurstWaveSearchIndex < lineInfo.colorBurstWaves.size())
		{
			const ColorBurstWaveInfo& waveInfo = lineInfo.colorBurstWaves[colorBurstWaveSearchIndex];
			if ((waveInfo.isPositive == referenceColorBurstPositive) || (referenceColorBurstLineOdd != lineIsOdd))
			{
				double colorBurstStartPosLineRelative = waveInfo.startPos - lineInfo.backPorchStartPos;
				if ((colorBurstStartPosLineRelative >= (referenceColorBurstStartPosLineRelative - colorBurstLineSyncRange)) && (colorBurstStartPosLineRelative <= (referenceColorBurstStartPosLineRelative + colorBurstLineSyncRange)))
				{
					syncCorrection = colorBurstStartPosLineRelative - referenceColorBurstStartPosLineRelative;
					foundSyncCorrection = true;
					break;
				}
			}
			++colorBurstWaveSearchIndex;
		}
		if (!foundSyncCorrection)
		{
			_logger.Error("Color burst sync correction failure! Start sample: {0}, End sample: {1}", lineInfo.leadingSyncInfo.startSampleNo, lineInfo.followingSyncInfo.startSampleNo);
			continue;
		}

		// Calculate a sync correction offset for this line based on our phase lock with the reference colour burst wave
		syncCorrectionValues[i] = syncCorrection;
		accumulatedSyncCorrection += syncCorrection;
		++syncCorrectionEntries;
	}

	// We used a single line to establish our colour burst phase lock, but that line may deviate from the average start
	// position for the lines in this field. To compensate for that, adjust the calculated sync correction values by an
	// accumulated sync correction value across the field, to ensure we align all the lines together at the mean start
	// position for this field.
	double averageSyncCorrection = accumulatedSyncCorrection / (double)syncCorrectionEntries;

	// Apply sync correction for all lines in the field
	for (unsigned int i = 0; i < fieldInfo.lines.size(); ++i)
	{
		// Retrieve info for the current line
		LineInfo& lineInfo = fieldInfo.lines[i];
		double syncCorrectionForLine = syncCorrectionValues[i];
		if (std::isnan(syncCorrectionForLine))
		{
			continue;
		}

		// Offset the line sync correction value by the average correction value
		syncCorrectionForLine -= averageSyncCorrection;

		// Offset our line position by the sync correction
		lineInfo.colorBurstSyncCorrection = syncCorrectionForLine;
		lineInfo.backPorchStartPos += lineInfo.colorBurstSyncCorrection;
		lineInfo.frontPorchEndPos += lineInfo.colorBurstSyncCorrection;
	}

	// Since we've performed line sync correction using the colour burst, return true.
	return true;
}
