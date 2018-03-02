#include "SplineHelpers.h"
#include <sstream>

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
double FrameBuilder::FindSyncRisingEdgeSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo) const
{
	// Find the point at which we cross the minimum threshold for ending the sync run
	double syncEndMinimumThreshold = (SampleType)(syncInfo.averageSyncLevel + (syncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance));
	SampleType syncEndMinimumThresholdAsSampleType = (SampleType)syncEndMinimumThreshold;
	size_t syncEndSearchPos = syncInfo.endSampleNo;
	while (inputData[syncEndSearchPos] < syncEndMinimumThresholdAsSampleType)
	{
		++syncEndSearchPos;
	}
	double sampleOffset = CreateSplineCatmullRomUniform((double)inputData[syncEndSearchPos - 2], (double)inputData[syncEndSearchPos - 1], (double)inputData[syncEndSearchPos], (double)inputData[syncEndSearchPos + 1]).Reverse(syncEndMinimumThreshold, slopeDetectionTolerance);
	return (double)(syncEndSearchPos - 1) + sampleOffset;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
double FrameBuilder::FindSyncFallingEdgeSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo) const
{
	// Find the point at which we cross the minimum threshold for starting the sync run
	double syncStartMinimumThreshold = (syncInfo.averageSyncLevel + (syncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance));
	SampleType syncStartMinimumThresholdAsSampleType = (SampleType)syncStartMinimumThreshold;
	size_t syncStartSearchPos = syncInfo.startSampleNo;
	while (inputData[syncStartSearchPos] < syncStartMinimumThreshold)
	{
		--syncStartSearchPos;
	}
	double sampleOffset = CreateSplineCatmullRomUniform((double)inputData[syncStartSearchPos - 1], (double)inputData[syncStartSearchPos], (double)inputData[syncStartSearchPos + 1], (double)inputData[syncStartSearchPos + 2]).Reverse(syncStartMinimumThreshold, slopeDetectionTolerance);
	return (double)syncStartSearchPos + sampleOffset;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
double FrameBuilder::FindSyncRisingEdgeEndSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo, double risingEdgePos) const
{
	// Find the point at which the leading sync run levels out to an acceptably flat slope
	size_t searchPos = (size_t)risingEdgePos;
	while ((inputData[searchPos+1] > inputData[searchPos]) && (((double)(inputData[searchPos+1] - inputData[searchPos]) / syncInfo.approxMinMaxSampleRange) > slopeValueFlatTolerance))
	{
		++searchPos;
	}
	++searchPos;
	return (double)searchPos;
}

//----------------------------------------------------------------------------------------------------------------------
// Frame detection methods
//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
std::list<FrameBuilder::FieldInfo> FrameBuilder::DetectFields(const std::vector<SampleType>& sampleData, const std::list<SyncDetector::SyncInfo>& syncEvents) const
{
	// Build a list of complete fields from the list of sync events
	std::list<FieldInfo> fields;
	FieldInfo currentField;
	bool currentFieldHasKnownStartPos = false;
	bool foundHSyncInCurrentField = false;
	unsigned int lineCount = 0;
	size_t lastSyncEventCount = 0;
	for (auto syncInfo : syncEvents)
	{
		// If this is a vertical sync event and we've found at least one horizontal sync event in the current field,
		// we've found the start of a new field. In this case, we need to complete the last field and reset the current
		// field.
		//##FIX## NTSC signal specs state line 1 starts with the first equalizing pulse, while PAL and SECAM specs state
		//that line 1 starts with the first vsync pulse. We currently implement the NTSC approach here, but this should
		//be selectable by the user.
		if ((syncInfo.type != SyncDetector::SyncType::Horizontal) && foundHSyncInCurrentField)
		{
			// If the current field has a known start position, add it to the list of complete fields.
			if (currentFieldHasKnownStartPos)
			{
				lastSyncEventCount = currentField.syncEvents.size();
				currentField.lineCount = lineCount;
				currentField.endSampleNo = syncInfo.startSampleNo;
				currentField.followingSyncEvent = syncInfo;
				fields.push_back(std::move(currentField));
			}

			// Reset the current field
			currentField = FieldInfo();
			currentField.syncEvents.reserve(lastSyncEventCount);
			currentFieldHasKnownStartPos = true;
			foundHSyncInCurrentField = false;
			lineCount = 0;
		}

		// If this is a horizontal sync event, increment the line count for this field.
		if (syncInfo.type == SyncDetector::SyncType::Horizontal)
		{
			++lineCount;
		}

		// Flag whether we've found a hsync event in the current field
		foundHSyncInCurrentField |= (syncInfo.type == SyncDetector::SyncType::Horizontal);

		// If the current field doesn't have a known start position, advance to the next sync event.
		if (!currentFieldHasKnownStartPos)
		{
			continue;
		}

		// Add this sync event to the list of sync events for the current field
		currentField.syncEvents.push_back(syncInfo);
	}

	// Return the list of fields to the caller
	return fields;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void FrameBuilder::DetectLines(const std::vector<SampleType>& sampleData, std::vector<FrameInfo>& frames, unsigned int threadCount) const
{
	// Determine the number of threads to use for this operation
	size_t frameCount = frames.size();
	if (threadCount <= 0)
	{
		unsigned int coreCount = std::thread::hardware_concurrency();
		threadCount = (coreCount > 0) ? coreCount : 4;
	}
	if (threadCount > frameCount)
	{
		threadCount = 1;
	}

	// Spawn worker threads to perform line detection in parallel
	size_t chunkCount = threadCount;
	size_t framesPerChunk = frameCount / chunkCount;
	std::vector<std::thread> workerThreads;
	for (unsigned int i = 0; i < chunkCount; ++i)
	{
		workerThreads.emplace_back(std::thread([&, i]
		{
			size_t frameNo = i * framesPerChunk;
			size_t lastFrameNoForChunk = frameNo + framesPerChunk;
			while (frameNo < lastFrameNoForChunk)
			{
				for (FrameBuilder::FieldInfo& fieldInfo : frames[frameNo].fields)
				{
					DetectLines(sampleData, fieldInfo);
				}
				++frameNo;
			}
		}));
	}
	for (auto& entry : workerThreads)
	{
		entry.join();
	}
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void FrameBuilder::DetectLines(const std::vector<SampleType>& sampleData, FieldInfo& fieldInfo) const
{
	// Reserve enough space in the line buffer for up to the number of sync events we have in the field
	fieldInfo.lines.reserve(fieldInfo.syncEvents.size());

	// Build a set of lines from the sync events in this field
	auto syncEventIterator = fieldInfo.syncEvents.cbegin();
	auto syncEventStartLineIterator = syncEventIterator;
	while (syncEventIterator != fieldInfo.syncEvents.cend())
	{
		bool foundHalfLineEvent = false;
		while ((syncEventIterator != fieldInfo.syncEvents.cend()) && (syncEventIterator->type != SyncDetector::SyncType::Horizontal))
		{
			bool syncEventIsHalfLineEvent = (syncEventIterator->type == SyncDetector::SyncType::Equalizer) || (syncEventIterator->type == SyncDetector::SyncType::Vertical);
			if (foundHalfLineEvent && syncEventIsHalfLineEvent)
			{
				break;
			}
			foundHalfLineEvent |= syncEventIsHalfLineEvent;

			++syncEventIterator;
		}
		if (syncEventIterator != fieldInfo.syncEvents.cend())
		{
			++syncEventIterator;
			const SyncDetector::SyncInfo& nextSyncEvent = (syncEventIterator != fieldInfo.syncEvents.cend()) ? *syncEventIterator : fieldInfo.followingSyncEvent;
			LineInfo lineInfo = {};
			lineInfo.leadingSyncInfo = *syncEventStartLineIterator;
			lineInfo.followingSyncInfo = nextSyncEvent;
			fieldInfo.lines.emplace_back(std::move(lineInfo));
			syncEventStartLineIterator = syncEventIterator;
		}
	}

	// Process each detected line, building advanced information for each of them.
	for (LineInfo& lineInfo : fieldInfo.lines)
	{
		lineInfo.backPorchStartPos = FindSyncRisingEdgeSamplePos(sampleData, lineInfo.leadingSyncInfo);
		lineInfo.frontPorchEndPos = FindSyncFallingEdgeSamplePos(sampleData, lineInfo.followingSyncInfo);
		double backPorchFlatStartPos = FindSyncRisingEdgeEndSamplePos(sampleData, lineInfo.leadingSyncInfo, lineInfo.backPorchStartPos);

		// Calculate the average level of the back porch region for the line
		//##TODO## Determine the right approach here
		std::vector<SampleType> upsampledBackPorch;
		double upsampledBackPorchStartPos;
		double upsampledBackPorchEndPos;
		if (lineInfo.leadingSyncInfo.type == SyncDetector::SyncType::Horizontal)
		{
			size_t hsyncLength = lineInfo.leadingSyncInfo.endSampleNo - lineInfo.leadingSyncInfo.startSampleNo;
			double backPorchSafeLength = (double)hsyncLength * syncLengthToBackPorchMinRatio;
			upsampledBackPorchStartPos = backPorchFlatStartPos;
			upsampledBackPorchEndPos = upsampledBackPorchStartPos + backPorchSafeLength;
			upsampledBackPorch.resize((size_t)(backPorchSafeLength * blankingUpsampleRatio));
			CubicInterpolateCatmullRom(sampleData.data(), upsampledBackPorchStartPos, upsampledBackPorchEndPos, upsampledBackPorch);

			//##FIX## This gave very poor results. Run some tests to figure out why.
			//lineInfo.averageBlankingLevel = std::accumulate(upsampledBackPorch.begin(), upsampledBackPorch.end(), 0.0) / (double)upsampledBackPorch.size();
			//##DEBUG##
			double averageBlankingLevelOld = std::accumulate(upsampledBackPorch.begin(), upsampledBackPorch.end(), 0.0) / (double)upsampledBackPorch.size();
			auto backPorchMinMax = std::minmax_element(upsampledBackPorch.begin(), upsampledBackPorch.end());
			lineInfo.averageBlankingLevel = (double)*backPorchMinMax.first + (((double)*backPorchMinMax.second - (double)*backPorchMinMax.first) / 2);

			//size_t backPorchSafeLengthInSamples = (size_t)backPorchSafeLength;
			//size_t backPorchStartSampleNo = (size_t)(lineInfo.backPorchStartPos + 0.5);
			//size_t backPorchEndSampleNo = backPorchStartSampleNo + backPorchSafeLengthInSamples;
			//lineInfo.averageBlankingLevel = std::accumulate(sampleData.begin() + backPorchStartSampleNo, sampleData.begin() + backPorchEndSampleNo, 0.0) / (double)(backPorchEndSampleNo - backPorchStartSampleNo);
		}
		else
		{
			size_t backPorchStartSampleNo = (size_t)(backPorchFlatStartPos + 0.5);
			lineInfo.averageBlankingLevel = sampleData[backPorchStartSampleNo];
		}

		// Calculate the average level of the leading sync region for the line
		//##TODO## Compare this with the average calculated by the run
		std::vector<SampleType> upsampledSync((size_t)((lineInfo.leadingSyncInfo.endSampleNo - lineInfo.leadingSyncInfo.startSampleNo) * blankingUpsampleRatio));
		CubicInterpolateCatmullRom(sampleData.data(), (double)lineInfo.leadingSyncInfo.startSampleNo, (double)lineInfo.leadingSyncInfo.endSampleNo, upsampledSync);
		lineInfo.averageSyncLevel = std::accumulate(upsampledSync.begin(), upsampledSync.end(), 0.0) / (double)upsampledSync.size();

		// Calculate the IRE levels for this line from the average sync and blanking levels. Sync occurs at IRE -40, and
		// blanking occurs at IRE 0. From these two reference levels, we can establish our IRE scale.
		//##FIX## In a PAL signal, sync occurs at -43 IRE, and the colour burst has a 21.5 IRE amplitude with 10 cycles.
		lineInfo.ireLevel0 = lineInfo.averageBlankingLevel;
		lineInfo.ireLevel100 = lineInfo.averageBlankingLevel + ((lineInfo.averageBlankingLevel - lineInfo.averageSyncLevel) * (1.0 / 0.4));

		// Attempt to locate a colour burst signal in the back porch region of the line
		if (lineInfo.leadingSyncInfo.type == SyncDetector::SyncType::Horizontal)
		{
			// Attempt to detect the color burst for this line
			//##TODO## Handle this properly if the color burst is completely absent. We shouldn't consider monochrome
			//signals to be "broken".
			SampleType colorBurstMinAmplitude = (SampleType)((lineInfo.ireLevel100 - lineInfo.ireLevel0) * colorBurstIREDetectionThreshold);
			DetectColorBurst(upsampledBackPorch, (SampleType)lineInfo.ireLevel0, colorBurstMinAmplitude, lineInfo.colorBurstWaves);
			lineInfo.colorBurstValid = ValidateColorBurst(colorBurstMinimumHalfOscillationCount, lineInfo.colorBurstWaves);

			// If the colour burst isn't valid, attempt to repair it.
			if (!lineInfo.colorBurstValid)
			{
				lineInfo.colorBurstValid = RepairColorBurst(upsampledBackPorch, (SampleType)lineInfo.ireLevel0, colorBurstMinAmplitude, lineInfo.colorBurstWaves);
			}

			// If the colour burst repair failed, report an error.
			//##TODO## Add an extra processing step after line detection to predict missing colour bursts from
			//surrounding lines in the field. Do the same for fields to predict broken vsync signals, and for lines to
			//predict broken hsync pulses.
			if (!lineInfo.colorBurstValid)
			{
				//##DEBUG##
				std::cout << "Failed to detect colour burst on hsync!\t" << std::fixed << lineInfo.leadingSyncInfo.startSampleNo << "\t" << backPorchFlatStartPos << "\t" << lineInfo.backPorchStartPos << "\t" << lineInfo.frontPorchEndPos << "\n";
			}

			// Convert our colour burst wave positions back to absolute sample positions
			for (ColorBurstWaveInfo& waveInfo : lineInfo.colorBurstWaves)
			{
				waveInfo.startPos = upsampledBackPorchStartPos + (waveInfo.startPos / blankingUpsampleRatio);
				waveInfo.endPos = upsampledBackPorchStartPos + (waveInfo.endPos / blankingUpsampleRatio);
			}
		}
	}

	// Perform more precise sync correction using the phase of the colour burst if requested. Note that we do this here
	// rather than leaving it to caller to perform this extra step, as we handle automatic threading on this method, and
	// it's more efficient to perform this operation within the worker threads while line detection is being done in
	// parallel.
	if (useColorBurstForLineSyncCorrection)
	{
		PerformColorBurstLineSyncCorrection(fieldInfo);
	}
}

//----------------------------------------------------------------------------------------
// Colour burst methods
//----------------------------------------------------------------------------------------
template<class SampleType>
void FrameBuilder::DetectColorBurst(const std::vector<SampleType>& backPorchData, SampleType zeroLevel, SampleType burstAmplitude, std::vector<ColorBurstWaveInfo>& burstWaves) const
{
	// Detect all possible burst waves in the sample region
	burstWaves.reserve(30);
	size_t runStartPos = 0;
	SampleType runPeakSample = backPorchData[0];
	bool runIsPositive = (runPeakSample > zeroLevel);
	for (unsigned int i = 1; i < backPorchData.size(); ++i)
	{
		// If the current sample crosses the zero-level, accept the previous run as a burst wave if it meets the burst
		// amplitude, and start a new run from the current sample.
		SampleType currentSample = backPorchData[i];
		bool currentSampleIsPositive = (currentSample > zeroLevel);
		bool crossedZeroLevel = (runIsPositive != currentSampleIsPositive);
		if (crossedZeroLevel)
		{
			// If the previous run reached burst amplitude, add it to the list of detected burst waves.
			bool runReachedBurstAmplitude = (runIsPositive ? ((zeroLevel + burstAmplitude) <= runPeakSample) : ((zeroLevel - burstAmplitude) >= runPeakSample));
			if (runReachedBurstAmplitude)
			{
				ColorBurstWaveInfo waveInfo;
				waveInfo.isPositive = runIsPositive;
				waveInfo.peakLevel = (double)runPeakSample;
				waveInfo.startPos = (double)runStartPos;
				waveInfo.endPos = (double)i;
				waveInfo.isPredicted = false;
				burstWaves.push_back(std::move(waveInfo));
			}

			// Start a new run from the current sample, and advance to the next sample.
			runStartPos = i;
			runPeakSample = currentSample;
			runIsPositive = currentSampleIsPositive;
			continue;
		}

		// Update the peak sample level for this run
		runPeakSample = ((runIsPositive ? currentSample > runPeakSample : currentSample < runPeakSample) ? currentSample : runPeakSample);
	}
}

//----------------------------------------------------------------------------------------
template<class SampleType>
bool FrameBuilder::RepairColorBurst(const std::vector<SampleType>& backPorchData, SampleType zeroLevel, SampleType burstAmplitude, std::vector<ColorBurstWaveInfo>& burstWaves) const
{
	// Trim leading and trailing wave entries where successive waves have the same "polarity". This can happen
	// frequently due to noise and other artifacts.
	while ((burstWaves.size() > 1) && (burstWaves[0].isPositive == burstWaves[1].isPositive))
	{
		burstWaves.erase(burstWaves.begin());
	}
	while ((burstWaves.size() > 1) && (burstWaves[burstWaves.size()-2].isPositive == burstWaves[burstWaves.size()-1].isPositive))
	{
		burstWaves.resize(burstWaves.size()-1);
	}

	// If there aren't enough wave entries left to form a full oscillation, return false. We need at least two entries
	// in order to perform our repair steps below.
	if (burstWaves.size() <= 1)
	{
		return false;
	}

	// Build a map of burst wave lengths in samples, snapping close values together. We do this to help determine the
	// average wave length below. We don't want to take a simple arithmetic mean however, as we expect bad entries with
	// wildly different lengths to be common, and we don't want them skewing the results. We use these numbers to
	// essentially "vote" on the average wave length here. Also note that we track separate averages for positive and
	// negative burst wave entries. We do this because it's possible our blanking level average is off, possibly from
	// noise or other interference, in which case our average sample length may be different for positive and negative
	// pulses.
	//##TODO## Compare this approach with an arithmetic median calculation.
	const double colorBurstAverageCalculationTolerance = 0.05;
	std::map<double, unsigned int> waveLengthAveragesPositive;
	std::map<double, unsigned int> waveLengthAveragesNegative;
	for (size_t i = 0; i < burstWaves.size(); ++i)
	{
		// Start scanning from the middle of the burst wave entries. We do this, because we expect bad wave entries
		// to be most frequent on either end of the range, and we want to start building our averages from good
		// entries first, since we "snap" length averages to previous entries, and we don't want to seed our
		// averages with bad data.
		size_t waveIndex = (i + (burstWaves.size() / 2));
		waveIndex = (waveIndex >= burstWaves.size() ? waveIndex - burstWaves.size() : waveIndex);

		// Build our set of common burst wave lengths
		double waveLength = burstWaves[waveIndex].endPos - burstWaves[waveIndex].startPos;
		double waveTolerance = waveLength * colorBurstAverageCalculationTolerance;
		auto& waveLengthAverages = (burstWaves[waveIndex].isPositive ? waveLengthAveragesPositive : waveLengthAveragesNegative);
		auto waveLengthAveragesIterator = waveLengthAverages.begin();
		while (waveLengthAveragesIterator != waveLengthAverages.end())
		{
			if ((waveLength <= (waveLengthAveragesIterator->first + waveTolerance)) && (waveLength >= (waveLengthAveragesIterator->first - waveTolerance)))
			{
				++waveLengthAveragesIterator->second;
				break;
			}
			++waveLengthAveragesIterator;
		}
		if (waveLengthAveragesIterator == waveLengthAverages.end())
		{
			waveLengthAverages[waveLength] = 1;
		}
	}

	// Select the highest occurring wave lengths as the average wave lengths for our positive and negative waves
	double waveLengthAveragePositive = 0.0;
	unsigned int waveLengthAverageHitCountPositive = 0;
	for (auto waveLengthAveragesEntry : waveLengthAveragesPositive)
	{
		if (waveLengthAveragesEntry.second >= waveLengthAverageHitCountPositive)
		{
			waveLengthAveragePositive = waveLengthAveragesEntry.first;
			waveLengthAverageHitCountPositive = waveLengthAveragesEntry.second;
		}
	}
	double waveLengthAverageNegative = 0.0;
	unsigned int waveLengthAverageHitCountNegative = 0;
	for (auto waveLengthAveragesEntry : waveLengthAveragesNegative)
	{
		if (waveLengthAveragesEntry.second >= waveLengthAverageHitCountNegative)
		{
			waveLengthAverageNegative = waveLengthAveragesEntry.first;
			waveLengthAverageHitCountNegative = waveLengthAveragesEntry.second;
		}
	}

	// Perform error correction on the burst waves
	double waveLengthAverageTolerancePositive = waveLengthAveragePositive * colorBurstRepairWaveLengthTolerance;
	double waveLengthAverageToleranceNegative = waveLengthAverageNegative * colorBurstRepairWaveLengthTolerance;
	double waveLengthAverageCombined = (waveLengthAveragePositive + waveLengthAverageNegative) / 2.0;
	double waveLengthAverageToleranceCombined = (waveLengthAverageTolerancePositive + waveLengthAverageToleranceNegative) / 2.0;
	unsigned int burstWavesIndex = 0;
	while (burstWavesIndex < (burstWaves.size() - 1))
	{
		double waveLengthAverage;
		double waveLengthAverageOpposite;
		double waveLengthAverageTolerance;
		double waveLengthAverageToleranceOpposite;
		if (burstWaves[burstWavesIndex].isPositive)
		{
			waveLengthAverage = waveLengthAveragePositive;
			waveLengthAverageOpposite = waveLengthAverageNegative;
			waveLengthAverageTolerance = waveLengthAverageTolerancePositive;
			waveLengthAverageToleranceOpposite = waveLengthAverageToleranceNegative;
		}
		else
		{
			waveLengthAverage = waveLengthAverageNegative;
			waveLengthAverageOpposite = waveLengthAveragePositive;
			waveLengthAverageTolerance = waveLengthAverageToleranceNegative;
			waveLengthAverageToleranceOpposite = waveLengthAverageTolerancePositive;
		}

		// Calculate the length of the next burst wave, and the length of the gap to the next burst wave. In a valid
		// colour burst, there should be no gap between each entry, as we detect each wave start and end position at the
		// zero crossing point relative to the blanking level.
		double waveLength = burstWaves[burstWavesIndex].endPos - burstWaves[burstWavesIndex].startPos;
		double waveGapLength = burstWaves[burstWavesIndex+1].startPos - burstWaves[burstWavesIndex].endPos;

		// If a gap is present between this burst wave and the next burst wave, calculate various parameters we'll need
		// to assess if we can generate entries to fill the gap.
		double followingWaveLength;
		double followingWaveLengthAverage;
		double followingWaveLengthAverageTolerance;
		unsigned int averageWaveGapFillCount;
		double minGapLengthForFill;
		double maxGapLengthForFill;
		if ((waveGapLength > 0) && (burstWavesIndex > 0))
		{
			followingWaveLength = burstWaves[burstWavesIndex+1].endPos - burstWaves[burstWavesIndex+1].startPos;
			followingWaveLengthAverage = (burstWaves[burstWavesIndex+1].isPositive ? waveLengthAveragePositive : waveLengthAverageNegative);
			followingWaveLengthAverageTolerance = (burstWaves[burstWavesIndex+1].isPositive ? waveLengthAverageTolerancePositive : waveLengthAverageToleranceNegative);
			averageWaveGapFillCount = (unsigned int)((waveGapLength + (waveLengthAverageCombined / 2)) / waveLengthAverageCombined);
			minGapLengthForFill = averageWaveGapFillCount * (waveLengthAverageCombined - waveLengthAverageToleranceCombined);
			maxGapLengthForFill = averageWaveGapFillCount * (waveLengthAverageCombined + waveLengthAverageToleranceCombined);
		}

		// Perform error correction on this wave entry
		if (!burstWaves[burstWavesIndex].isPredicted && ((waveLength < (waveLengthAverage - waveLengthAverageTolerance)) || (waveLength > (waveLengthAverage + waveLengthAverageTolerance))))
		{
			// If we haven't accepted any wave entries as valid yet, and this wave entry doesn't fit within our
			// tolerance of the average wave length, remove it. Note that we skip this step for predicted wave entries.
			// Our predicted waves should generally be within tolerance, but due to round off error they may fall very
			// slightly outside the range when created. If we don't exclude them here, we can enter an infinite loop.
			burstWaves.erase(burstWaves.begin() + burstWavesIndex);

			// If the wave entry we erased wasn't the first one, step back to the previous wave entry so we can
			// compare the new next entry with the previous one. If we don't do this and the entry we just erased
			// was the second last one for example, we'll abort the loop without even examining the last entry.
			if (burstWavesIndex > 0)
			{
				--burstWavesIndex;
			}
			continue;
		}
		else if
		(
		    // A gap is present before the next wave, and we've already found at least two waves that match.
		    (waveGapLength > 0)
		    && (burstWavesIndex > 0)
		    // The following wave is within tolerance of its expected length
		    && (followingWaveLength >= (followingWaveLengthAverage - followingWaveLengthAverageTolerance))
		    && (followingWaveLength <= (followingWaveLengthAverage + followingWaveLengthAverageTolerance))
		    &&
		    (
		        // We can fit a single wave oscillation in the gap, and the polarity of the resulting oscillations will be consistent.
		        ((waveGapLength >= (waveLengthAverageOpposite - waveLengthAverageToleranceOpposite)) && (waveGapLength <= (waveLengthAverageOpposite + waveLengthAverageToleranceOpposite)) && (burstWaves[burstWavesIndex].isPositive == burstWaves[burstWavesIndex+1].isPositive))
		        // We can fit multiple wave oscillations in the gap, and the polarity of the resulting oscillations will be consistent.
		     || ((waveGapLength >= minGapLengthForFill) && (waveGapLength <= maxGapLengthForFill) && (((averageWaveGapFillCount % 2) != 0) == (burstWaves[burstWavesIndex].isPositive == burstWaves[burstWavesIndex+1].isPositive)))
		    )
		)
		{
			// If we need to insert a wave entry before the next wave entry in order to close a gap, do it now.
			//##TODO## Consider using the separate averages for positive and negative waves to fill the gap. We'll need
			// to scale each one by the overall gap length.
			double gapFillWaveLength = waveGapLength / (double)averageWaveGapFillCount;
			while (waveGapLength > waveLengthAverageToleranceCombined)
			{
				ColorBurstWaveInfo predictedWaveInfo;
				predictedWaveInfo.isPositive = !burstWaves[burstWavesIndex].isPositive;
				predictedWaveInfo.startPos = burstWaves[burstWavesIndex].endPos;
				predictedWaveInfo.endPos = predictedWaveInfo.startPos + gapFillWaveLength;
				predictedWaveInfo.endPos = (predictedWaveInfo.endPos > burstWaves[burstWavesIndex+1].startPos) ? burstWaves[burstWavesIndex+1].startPos : predictedWaveInfo.endPos;
				predictedWaveInfo.peakLevel = burstWaves[burstWavesIndex].peakLevel;
				predictedWaveInfo.isPredicted = true;
				burstWaves.insert(burstWaves.begin() + (burstWavesIndex + 1), predictedWaveInfo);

				++burstWavesIndex;
				waveGapLength = burstWaves[burstWavesIndex+1].startPos - burstWaves[burstWavesIndex].endPos;
			}

			// If we've left a small gap to the next wave entry, close it now. This could occur due to round off error.
			if (waveGapLength > 0)
			{
				burstWaves[burstWavesIndex].endPos = burstWaves[burstWavesIndex+1].startPos;
			}
			continue;
		}
		else if ((waveGapLength > 0) || (burstWaves[burstWavesIndex].isPositive == burstWaves[burstWavesIndex+1].isPositive))
		{
			// If there's a gap to the next wave entry but it's not large enough to fit a valid wave oscillation into,
			// or if we don't have differing polarity of two consecutive wave entries, we need to erase either the
			// current wave entry or following wave entry. If we haven't found at least two good consecutive burst waves
			// (one full oscillation) yet, we assume the leading sample is wrong and that's why we have a gap. In this
			// case, we erase the current wave entry, otherwise we assume the following sample is wrong, in which case
			// we erase the following wave sample.
			if (burstWavesIndex < 2)
			{
				// Erase the current burst wave entry
				burstWaves.erase(burstWaves.begin() + burstWavesIndex);

				// If the wave entry we erased wasn't the first one, step back to the previous wave entry so we can
				// compare the new next entry with the previous one. If we don't do this and the entry we just erased
				// was the second last one for example, we'll abort the loop without even examining the last entry.
				if (burstWavesIndex > 0)
				{
					--burstWavesIndex;
				}
			}
			else
			{
				// Erase the following burst wave entry
				burstWaves.erase(burstWaves.begin() + (burstWavesIndex + 1));
			}
			continue;
		}

		// Advance to the next wave entry
		++burstWavesIndex;
	}

	// Attempt to validate the color burst after the repair operation, and return the result the caller.
	bool colorBurstValid = ValidateColorBurst(colorBurstMinimumHalfOscillationCount, burstWaves);
	return colorBurstValid;
}

//----------------------------------------------------------------------------------------
// IRE conversion methods
//----------------------------------------------------------------------------------------
template<class SampleType>
float FrameBuilder::SampleToIRE(SampleType sampleValue, float ireLevel0, float ireLevel100) const
{
	return ((float)sampleValue - ireLevel0) * (100.0f / (ireLevel100 - ireLevel0));
}

//----------------------------------------------------------------------------------------
template<class SampleType>
SampleType FrameBuilder::IREToSample(float ire, float ireLevel0, float ireLevel100) const
{
	return (SampleType)(((ire * ((ireLevel100 - ireLevel0) / 100.0f)) + ireLevel0) + (std::numeric_limits<T>::is_integer ? 0.5f : 0.0f));
}
