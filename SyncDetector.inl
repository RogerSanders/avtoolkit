#include <numeric>
#include <cmath>
#include <algorithm>
#include <thread>

//----------------------------------------------------------------------------------------------------------------------
// Sync detection methods
//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
std::list<SyncDetector::SyncPulseInfo> SyncDetector::DetectSyncPulses(const std::vector<SampleType>& sampleData, size_t initialSampleNo, size_t sampleCount, size_t initialMinMaxSampleNo, size_t minMaxSampleCount, unsigned int threadCount) const
{
	// If the initial sample number isn't valid, abort any further processing.
	if (initialSampleNo >= sampleData.size())
	{
		return std::list<SyncPulseInfo>();
	}

	// If no sample count has been specified or it is out of range, default to all sample data from the initial sample.
	sampleCount = ((sampleCount <= 0) || ((sampleCount + initialSampleNo) > sampleData.size())) ? sampleData.size() - initialSampleNo: sampleCount;
	minMaxSampleCount = ((minMaxSampleCount <= 0) || ((minMaxSampleCount + initialMinMaxSampleNo) > sampleData.size())) ? sampleData.size() - initialMinMaxSampleNo: minMaxSampleCount;
	size_t lastSampleNo = initialSampleNo + sampleCount;
	size_t lastMinMaxSampleNo = initialMinMaxSampleNo + minMaxSampleCount;

	// Determine the number of threads to use for this operation
	if (threadCount <= 0)
	{
		unsigned int coreCount = std::thread::hardware_concurrency();
		threadCount = (coreCount > 0) ? coreCount : 4;
	}
	if (threadCount > sampleCount)
	{
		threadCount = 1;
	}

	// If we're not using a sliding min/max window, retrieve the min/max values for run detection by scanning over all
	// the sample data. Since we're using this same calculated info for the entire sample set, it's more efficient to
	// pre-calculate the min/max info where threading is being used, as this information can be calculated once and
	// shared between all threads.
	MinMaxWindowInfo<SampleType> minMaxWindowInfo;
	minMaxWindowInfo.slidingWindowEnabled = enableMinMaxSlidingWindow;
	minMaxWindowInfo.scanSampleStartNo = initialMinMaxSampleNo;
	minMaxWindowInfo.scanSampleEndNo = lastMinMaxSampleNo;
	minMaxWindowInfo.minMaxValuesPopulated = false;
	if (!enableMinMaxSlidingWindow)
	{
		ObtainMinMaxValues(sampleData, initialSampleNo, minMaxWindowInfo);
	}

	// Build our initial list of run entries for the sample data
	std::list<RunEntry<SampleType>> runEntries;
	if (threadCount < 2)
	{
		// Build a list of runs in the sample data using a single linear process on one thread
		CutRunEntry<SampleType> cutRunEntry;
		runEntries = ExtractRunsFromSampleData(sampleData, initialSampleNo, lastSampleNo, minMaxWindowInfo, cutRunEntry, false, false, false);
	}
	else
	{
		// Break the sample data into chunks, and build separate lists of runs from each chunk in parallel using
		// multiple threads.
		unsigned int chunkCount = threadCount;
		size_t samplesPerChunk = sampleCount / chunkCount;
		std::vector<std::thread> extractRunsThreads;
		std::vector<std::list<RunEntry<SampleType>>> extractRunsResults(chunkCount);
		std::vector<size_t> chunkStartPositions(chunkCount);
		std::vector<size_t> chunkEndPositions(chunkCount);
		std::vector<CutRunEntry<SampleType>> cutRunEntries(chunkCount);
		for (unsigned int i = 0; i < chunkCount; ++i)
		{
			extractRunsThreads.emplace_back(std::thread([&, i]
			{
				chunkStartPositions[i] = initialSampleNo + (i * samplesPerChunk);
				chunkEndPositions[i] = (i == (threadCount - 1)) ? lastSampleNo : chunkStartPositions[i] + samplesPerChunk;
				MinMaxWindowInfo<SampleType> minMaxWindowInfoForChunk = minMaxWindowInfo;
				extractRunsResults[i] = ExtractRunsFromSampleData(sampleData, chunkStartPositions[i], chunkEndPositions[i], minMaxWindowInfoForChunk, cutRunEntries[i], false, false, false);
			}));
		}
		for (auto& entry : extractRunsThreads)
		{
			entry.join();
		}

		// Merge the separate lists of runs from each chunk together
		for (unsigned int i = 0; i < chunkCount; ++i)
		{
			MergeChunkResults(sampleData, lastSampleNo, runEntries, chunkStartPositions[i], chunkEndPositions[i], extractRunsResults[i], cutRunEntries[i], minMaxWindowInfo);
		}
	}

	// We have a list of run entries, but we now need to filter down to a set of entries with the lowest values, which
	// show a reasonably stable level between successive occurrences. This should correctly isolate sync runs from runs
	// occurring during blanking or active scan. Note that the output from this isn't corrected for "broken" runs at
	// this stage caused by spikes or incorrectly embedded colour bursts during sync pulses.
	FilterRunEntriesToSyncCandidates(runEntries);

	// We've correctly detected sync runs by this point, but now we need to error correct. We do this by merging close
	// runs together, if the data separating them passes a series of tests around length and average value. The result
	// of this stage should be a series of runs representing sync pulses in the source data.
	ErrorCorrectRunEntrySyncCandidates(sampleData, runEntries);

	// At this point, we're in really good shape. The only anomaly is the possibility that merged runs with strong
	// swings at the join point have resulted in batches of samples inappropriately included or omitted at the ends of
	// the merged run. We compensate for that here by attempting to extend each run forwards and backwards, then
	// performing a final trim on both ends to clean it up.
	CleanRunEntryEdges(sampleData, runEntries, initialSampleNo, lastSampleNo);

	//##FIX## Our final cleaning stage needs extra work to more aggressively take in samples which are separated by a
	//trivial number of samples from our run, which are in range of the average or below it. Right now a minor bump
	//prevents edge samples being folded back in.
	//##FIX## See sample run 176087977 length 125 in FantasiaCompositeSigned.bin for the target case.

	// Build our sync pulse objects from the run entries
	std::list<SyncPulseInfo> syncEntries;
	for (const auto& runEntry : runEntries)
	{
		// Populate the sync entry
		SyncPulseInfo syncEntry;
		syncEntry.startSampleNo = runEntry.startSampleNo;
		syncEntry.endSampleNo = runEntry.endSampleNo;
		syncEntry.approxMinMaxSampleRange = runEntry.initialMinMaxSampleRange;

		// Calculate the average sample value of the sync region
		syncEntry.averageSyncLevel = std::accumulate(sampleData.begin() + syncEntry.startSampleNo, sampleData.begin() + syncEntry.endSampleNo, 0.0) / (double)(syncEntry.endSampleNo - syncEntry.startSampleNo);

		// Add this sync entry to the list of identified sync entries
		syncEntries.push_back(syncEntry);
	}

	// Return the list of sync pulse objects to the caller
	return syncEntries;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
std::list<SyncDetector::SyncInfo> SyncDetector::DetectSyncEvents(const std::vector<SampleType>& sampleData, const std::list<SyncPulseInfo>& syncPulses) const
{
	// At this point, we start assuming significant aspects of video encoding. The first assumption is that the sync
	// pulses we've identified are in fact hsync and vsync pulses, and possibly equalizer pulses. The second assumption
	// is that the basic target horizontal and vertical sync rate of the video stream doesn't change throughout the
	// sample data. We still allow for large variation in actual intervals between sync pulses, however barring
	// supported variations like changing between interlaced and progressive scan video, we make the assumption here
	// that the original intent of the data stream is to output video in a consistent mode, and variations in time are
	// an anomaly in need of correction.

	// At this point, we hopefully have a set of sync pulses falling into two or three main length categories. We do a
	// quick count of various sync count lengths here:
	std::map<size_t, size_t> syncPulseLengthCounts;
	for (SyncPulseInfo pulseInfo : syncPulses)
	{
		size_t pulseLength = pulseInfo.endSampleNo - pulseInfo.startSampleNo;
		++syncPulseLengthCounts[pulseLength];
	}

	// The most common sync pulse should be horizontal sync by far, so we start by detecting those pulses.
	size_t hsyncAveragePulseLength = GetMostCommonSyncPulseLength(syncPulseLengthCounts);
	size_t hsyncPulseLengthTolerance = (size_t)((double)hsyncAveragePulseLength * syncPulseLengthVariationTolerance);

	// Erase hsync pulse events from the sync pulse length counts
	auto syncPulseLengthCountIterator = syncPulseLengthCounts.begin();
	while (syncPulseLengthCountIterator != syncPulseLengthCounts.end())
	{
		size_t syncPulseLength = syncPulseLengthCountIterator->first;
		auto nextIterator = std::next(syncPulseLengthCountIterator, 1);
		if ((size_t)std::abs((long long)syncPulseLength - (long long)hsyncAveragePulseLength) <= hsyncPulseLengthTolerance)
		{
			syncPulseLengthCounts.erase(syncPulseLengthCountIterator);
		}
		syncPulseLengthCountIterator = nextIterator;
	}

	// Extract equalizer pulses if present. The second most common occurring sample should either be vsync pulses if
	// there are no equalizer pulses present, or equalizer pulses. Since equalizer pulses are shorter than hsync, and
	// vsync pulses are larger, we check the length of the next largest pulse length to determine whether equalizer
	// pulses are present, and extract them if they are.
	bool hasEqualizerPulses = false;
	size_t vsyncAveragePulseLength = 0;
	size_t vsyncPulseLengthTolerance = 0;
	size_t equalizerAveragePulseLength = 0;
	size_t equalizerPulseLengthTolerance = 0;
	size_t nextAveragePulseLength = GetMostCommonSyncPulseLength(syncPulseLengthCounts);
	if (nextAveragePulseLength < hsyncAveragePulseLength)
	{
		// Equalizer pulses are present, so record information on the equalizer pulses
		equalizerAveragePulseLength = nextAveragePulseLength;
		equalizerPulseLengthTolerance = (size_t)((double)equalizerAveragePulseLength * syncPulseLengthVariationTolerance);
		hasEqualizerPulses = true;

		// Erase equalizer pulse events from the sync pulse length counts
		syncPulseLengthCountIterator = syncPulseLengthCounts.begin();
		while (syncPulseLengthCountIterator != syncPulseLengthCounts.end())
		{
			size_t syncPulseLength = syncPulseLengthCountIterator->first;
			auto nextIterator = std::next(syncPulseLengthCountIterator, 1);
			if ((size_t)std::abs((long long)syncPulseLength - (long long)equalizerAveragePulseLength) <= equalizerPulseLengthTolerance)
			{
				syncPulseLengthCounts.erase(syncPulseLengthCountIterator);
			}
			syncPulseLengthCountIterator = nextIterator;
		}

		// Obtain the next average pulse length, which should now refer to vsync pulses.
		nextAveragePulseLength = GetMostCommonSyncPulseLength(syncPulseLengthCounts);
	}

	// Vsync pulses should be all that's left, so extract them from the set of sync pulses.
	vsyncAveragePulseLength = nextAveragePulseLength;
	vsyncPulseLengthTolerance = (size_t)((double)vsyncAveragePulseLength * syncPulseLengthVariationTolerance);

	// Erase vsync pulse events from the sync pulse length counts
	syncPulseLengthCountIterator = syncPulseLengthCounts.begin();
	while (syncPulseLengthCountIterator != syncPulseLengthCounts.end())
	{
		size_t syncPulseLength = syncPulseLengthCountIterator->first;
		auto nextIterator = std::next(syncPulseLengthCountIterator, 1);
		if ((size_t)std::abs((long long)syncPulseLength - (long long)vsyncAveragePulseLength) <= vsyncPulseLengthTolerance)
		{
			syncPulseLengthCounts.erase(syncPulseLengthCountIterator);
		}
		syncPulseLengthCountIterator = nextIterator;
	}

	// Build our list of sync info structures from the sync pulses
	std::list<SyncInfo> syncInfoList;
	for (SyncPulseInfo pulseInfo : syncPulses)
	{
		// Obtain the next sync pulse
		size_t pulseLength = pulseInfo.endSampleNo - pulseInfo.startSampleNo;

		// Attempt to detect the sync type of this sync pulse. If no type can be determined, skip it.
		SyncType syncType;
		if ((size_t)std::abs((long long)pulseLength - (long long)hsyncAveragePulseLength) <= hsyncPulseLengthTolerance)
		{
			syncType = SyncType::Horizontal;
		}
		else if ((size_t)std::abs((long long)pulseLength - (long long)equalizerAveragePulseLength) <= equalizerPulseLengthTolerance)
		{
			syncType = SyncType::Equalizer;
		}
		else if ((size_t)std::abs((long long)pulseLength - (long long)vsyncAveragePulseLength) <= vsyncPulseLengthTolerance)
		{
			syncType = SyncType::Vertical;
		}
		else
		{
			syncType = SyncType::Unknown;
		}

		// Add this sync pulse to the list of detected sync events
		SyncInfo syncInfo;
		syncInfo.type = syncType;
		syncInfo.startSampleNo = pulseInfo.startSampleNo;
		syncInfo.endSampleNo = pulseInfo.endSampleNo;
		syncInfo.averageSyncLevel = pulseInfo.averageSyncLevel;
		syncInfo.approxMinMaxSampleRange = pulseInfo.approxMinMaxSampleRange;
		syncInfoList.push_back(syncInfo);
	}

	// Return the set of detected sync events to the caller
	return syncInfoList;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::ObtainMinMaxValues(const std::vector<SampleType>& sampleData, size_t samplePos, MinMaxWindowInfo<SampleType>& minMaxInfo) const
{
	// Obtain the min/max values for the target sample position
	size_t minMaxSlidingWindowStartPos = minMaxInfo.scanSampleStartNo;
	size_t minMaxSlidingWindowEndPos = minMaxInfo.scanSampleEndNo;
	if (minMaxInfo.slidingWindowEnabled)
	{
		// Update the sliding min/max window for the current position
		size_t minMaxSlidingWindowQuarterLengthInSamples = (size_t)(minMaxSlidingWindowLength * (double)(minMaxInfo.scanSampleEndNo - minMaxInfo.scanSampleStartNo)) / 4;
		size_t minMaxSlidingWindowHalfLengthInSamples = minMaxSlidingWindowQuarterLengthInSamples * 2;
		size_t windowReferencePos = (samplePos / minMaxSlidingWindowHalfLengthInSamples) * minMaxSlidingWindowHalfLengthInSamples;
		minMaxSlidingWindowStartPos = ((minMaxInfo.scanSampleStartNo + minMaxSlidingWindowQuarterLengthInSamples) >= windowReferencePos) ? minMaxInfo.scanSampleStartNo : windowReferencePos - minMaxSlidingWindowQuarterLengthInSamples;
		minMaxSlidingWindowEndPos = windowReferencePos + minMaxSlidingWindowHalfLengthInSamples + minMaxSlidingWindowQuarterLengthInSamples;
		minMaxInfo.scanStartPos = windowReferencePos;
		minMaxInfo.scanResetPos = minMaxSlidingWindowEndPos - minMaxSlidingWindowQuarterLengthInSamples;
		if (minMaxSlidingWindowEndPos > minMaxInfo.scanSampleEndNo)
		{
			minMaxSlidingWindowEndPos = minMaxInfo.scanSampleEndNo;
			minMaxInfo.scanResetPos = windowReferencePos + minMaxSlidingWindowHalfLengthInSamples;
		}
	}

	// Retrieve the new min/max values from the window
	auto minMaxResult = std::minmax_element(sampleData.cbegin() + minMaxSlidingWindowStartPos, sampleData.begin() + minMaxSlidingWindowEndPos);
	size_t minValueIndex = minMaxResult.first - sampleData.cbegin();
	size_t maxValueIndex = minMaxResult.second - sampleData.cbegin();
	minMaxInfo.minSampleValue = sampleData[minValueIndex];
	minMaxInfo.maxSampleValue = sampleData[maxValueIndex];
	minMaxInfo.minMaxSampleRange = (double)minMaxInfo.maxSampleValue - (double)minMaxInfo.minSampleValue;
	minMaxInfo.minMaxValuesPopulated = true;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
std::list<SyncDetector::RunEntry<SampleType>> SyncDetector::ExtractRunsFromSampleData(const std::vector<SampleType>& sampleData, size_t sampleStartNo, size_t sampleEndNo, MinMaxWindowInfo<SampleType>& minMaxInfo, CutRunEntry<SampleType>& cutRunEntry, bool cutRunEntryPopulated, bool resolveCutRunEntryOnly, bool getNextRunEntryOnly) const
{
	// Populate an initial run entry to fill as we perform run detection
	size_t currentSampleNo = sampleStartNo;
	RunEntry<SampleType> runEntry;
	if (cutRunEntryPopulated)
	{
		runEntry = cutRunEntry.runEntry;
	}
	else
	{
		runEntry.startSampleNo = sampleStartNo;
		runEntry.minValue = sampleData[sampleStartNo];
		runEntry.maxValue = runEntry.minValue;
		++currentSampleNo;
	}

	// Build our list of run entries for the sample set of data
	MinMaxWindowInfo<SampleType> minMaxInfoAlternate = minMaxInfo;
	std::list<RunEntry<SampleType>> runEntries;
	bool minMaxSampleRangeUpdateRequired = true;
	double minMaxRunRangeDevianceInSamples = 0;
	while (currentSampleNo < sampleEndNo)
	{
		// If we're using a sliding min/max window and we need to update it for the new position, do so now.
		if (!minMaxInfo.minMaxValuesPopulated || (minMaxInfo.slidingWindowEnabled && ((currentSampleNo >= minMaxInfo.scanResetPos) || (currentSampleNo < minMaxInfo.scanStartPos))))
		{
			// We keep a pair of the most recent min/max structures here, so we can swap back between them efficiently.
			// Since we can end up scanning back over the same region numerous times when looking for a run, keeping a
			// pair of entries buffered here avoids repeated and expensive resolution of the min/max info when a sliding
			// window is in use.
			if (!minMaxInfoAlternate.minMaxValuesPopulated || (minMaxInfoAlternate.slidingWindowEnabled && ((currentSampleNo >= minMaxInfoAlternate.scanResetPos) || (currentSampleNo < minMaxInfoAlternate.scanStartPos))))
			{
				// Shift the current min/max info into the alternate buffer
				minMaxInfoAlternate = minMaxInfo;

				// Retrieve the new min/max values for the current position
				ObtainMinMaxValues(sampleData, currentSampleNo, minMaxInfo);

				// Since we've just updated the min/max window, flag that the min/max sample range needs to be updated.
				minMaxSampleRangeUpdateRequired = true;
			}
			else
			{
				// Switch the current min/max info with the alternate buffer
				std::swap(minMaxInfo, minMaxInfoAlternate);
			}
		}

		// If we updated the min/max values for the new sliding window position, update any dependent values now.
		if (minMaxSampleRangeUpdateRequired)
		{
			minMaxRunRangeDevianceInSamples = (minMaxInfo.minMaxSampleRange * runRangeValueDeviance);
			minMaxSampleRangeUpdateRequired = false;
		}

		// Calculate the effect on the sample deviance in the current run if we add this sample into it
		SampleType currentSample = sampleData[currentSampleNo];
		SampleType minElementInRun = std::min(runEntry.minValue, currentSample);
		SampleType maxElementInRun = std::max(runEntry.maxValue, currentSample);
		double minMaxRunRange = (double)maxElementInRun - (double)minElementInRun;

		// If the current sample can be folded into the existing run, combine it now, and advance to the next sample.
		if (minMaxRunRange <= minMaxRunRangeDevianceInSamples)
		{
			runEntry.minValue = minElementInRun;
			runEntry.maxValue = maxElementInRun;
			++currentSampleNo;
			continue;
		}

		// If there's only a single sample in the current run, reset it to the current sample, and advance to the next
		// sample.
		if ((runEntry.startSampleNo + 1) == currentSampleNo)
		{
			if (resolveCutRunEntryOnly)
			{
				break;
			}
			runEntry.startSampleNo = currentSampleNo;
			runEntry.minValue = currentSample;
			runEntry.maxValue = currentSample;
			++currentSampleNo;
			continue;
		}

		// Calculate the effect on the sample deviance in the current run if we add this sample, and drop the first
		// sample.
		SampleType minElementInTrimmedRun = *std::min_element(sampleData.cbegin() + runEntry.startSampleNo + 1, sampleData.cbegin() + currentSampleNo + 1);
		SampleType maxElementInTrimmedRun = *std::max_element(sampleData.cbegin() + runEntry.startSampleNo + 1, sampleData.cbegin() + currentSampleNo + 1);
		double minMaxTrimmedRunRange = (double)maxElementInTrimmedRun - (double)minElementInTrimmedRun;

		// If we can continue the run by dropping the leading sample, do it now, and advance to the next sample.
		if (minMaxTrimmedRunRange <= minMaxRunRangeDevianceInSamples)
		{
			runEntry.startSampleNo += 1;
			runEntry.minValue = minElementInTrimmedRun;
			runEntry.maxValue = maxElementInTrimmedRun;
			++currentSampleNo;
			continue;
		}

		// We've reached the end of the current run. If it doesn't meet the minimum number of samples for a run, discard
		// it, and return to one sample past the start of the current run to start scanning again.
		if ((currentSampleNo - runEntry.startSampleNo) < minRunRangeSampleLength)
		{
			if (resolveCutRunEntryOnly)
			{
				break;
			}

			// Reset the current run entry to start from one past the original start position for the discarded run
			runEntry.startSampleNo = runEntry.startSampleNo + 1;
			runEntry.minValue = sampleData[runEntry.startSampleNo];
			runEntry.maxValue = runEntry.minValue;

			// Reset the current search position to the next potential sample in the reset run. Note that this is two
			// samples past the original start position for the discarded run, as we've already pulled one sample in to
			// start the new run.
			currentSampleNo = runEntry.startSampleNo + 1;
			continue;
		}

		// We're terminating this run, so record this sample as the end sample number.
		runEntry.endSampleNo = currentSampleNo;

		// We've reached the end of the current run, and now we need to trim it. A run can only be terminated by a
		// positive or negative swing from the average, and we expect some samples on either end of the run may be the
		// start of this trend. As we're aiming specifically for sync pulses however, which are characterized by dips
		// rather than peaks, we don't trim any low samples here, only high samples. This has consistently been shown to
		// provide the best results, as it ensures only curves running upwards are clipped by the trimming operation.
		size_t runLengthBeforeTrim = runEntry.endSampleNo - runEntry.startSampleNo;
		size_t runLengthMaxTrimSampleCount = (size_t)((double)runLengthBeforeTrim * runRangeMaxTrimLength);
		double runRange = (double)runEntry.maxValue - (double)runEntry.minValue;
		SampleType trimThreshold = (SampleType)((runRange - (runRange * runRangeTrimTolerance)) / 2.0);
		SampleType trimSampleTooHighThreshold = runEntry.maxValue - trimThreshold;
		SampleType trimSampleTooLowThreshold = runEntry.minValue;
		TrimRunEntry(sampleData, runEntry, trimSampleTooHighThreshold, trimSampleTooLowThreshold, runLengthMaxTrimSampleCount);

		// Now that we've trimmed the run, calculate final values for the min, max, and average sample values in the
		// run.
		runEntry.minValue = *std::min_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
		runEntry.maxValue = *std::max_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
		runEntry.averageSample = std::accumulate(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo, 0.0) / (double)(runEntry.endSampleNo - runEntry.startSampleNo);

		// We're accepting this run entry as a valid run, so add it to the list of identified runs.
		runEntry.initialMinMaxSampleRange = minMaxInfo.minMaxSampleRange;
		runEntries.push_back(runEntry);

		// Reset the current run entry
		runEntry.startSampleNo = runEntry.endSampleNo;
		runEntry.minValue = currentSample;
		runEntry.maxValue = currentSample;

		// Advance to the next search point. Note that we need to rebase to the final run end position here, to account
		// for the fact that we trimmed the resulting run entry, and that run entry may now end prior to the current
		// position. We want to re-cover the trimmed samples in this case, so we reset the search position to the final
		// end position of the run.
		currentSampleNo = runEntry.endSampleNo;

		// If we've been requested to resolve the cut run entry or the first run entry only, abort further searching for
		// runs.
		if (resolveCutRunEntryOnly || getNextRunEntryOnly)
		{
			break;
		}
	}

	// Record information on the cut run at the end of this set of runs
	if (!resolveCutRunEntryOnly)
	{
		cutRunEntry.minMaxInfo = minMaxInfo;
		cutRunEntry.runEntry = runEntry;
		cutRunEntry.runEntry.endSampleNo = currentSampleNo;
	}

	// Return the list of run entries to the caller
	return runEntries;
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::MergeChunkResults(const std::vector<SampleType>& sampleData, size_t sampleEndNo, std::list<RunEntry<SampleType>>& runEntries, size_t chunkStartPos, size_t chunkEndPos, std::list<RunEntry<SampleType>>& chunkRunEntries, CutRunEntry<SampleType>& cutRunEntry, MinMaxWindowInfo<SampleType>& minMaxInfo) const
{
	// If we've already got combined sample data, it's possible there's some overlap between that data and the
	// data for this chunk, as we complete cut runs between the chunk boundaries, and part of that cut run may
	// have been detected as a run in this chunk. We correct any overlap here by modifying the run entries for
	// this chunk to eliminate redundancy before we merge them.
	if (!runEntries.empty() && !chunkRunEntries.empty())
	{
		auto chunkRunEntryIterator = chunkRunEntries.begin();
		while (chunkRunEntryIterator != chunkRunEntries.end())
		{
			// Retrieve the last completed run entry and the next unmerged chunk entry
			const auto& lastRunEntry = *runEntries.rbegin();
			const auto& firstChunkRunEntry = *chunkRunEntryIterator;

			// If the first run entry for this chunk entirely overlaps with the completed run set, remove it.
			// If we end up in sync with the completed entries afterwards, abort further join processing
			// otherwise we restart join processing with the next run entry.
			bool scanToEndOfChunk = false;
			if (lastRunEntry.endSampleNo >= firstChunkRunEntry.endSampleNo)
			{
				// If the run entry we're dropping has the same end position as the last completed run, we're in
				// sync with the completed run set now, so we flag that here and move on with erasing the
				// overlapping entry.
				bool runEntriesInSync = lastRunEntry.endSampleNo == firstChunkRunEntry.endSampleNo;

				// Since the next run entry entirely overlaps with the completed run entries, remove it now.
				chunkRunEntryIterator = chunkRunEntries.erase(chunkRunEntryIterator);

				// If there are more unresolved chunk entries remaining, either stop or restart the join
				// process.
				if (chunkRunEntryIterator != chunkRunEntries.end())
				{
					// Since we're in sync with the completed runs now, abort fixing the join location at the
					// start of this chunk and move on with the merge process.
					if (runEntriesInSync)
					{
						break;
					}

					// Since there's at least one more run entry remaining in this chunk, restart the loop and
					// compare the next run entry with the last resolved run entry.
					continue;
				}

				// Since there are no more entries in this chunk, flag that we need to scan to the end of the
				// chunk region to resolve the next run.
				scanToEndOfChunk = true;
			}

			// The first run entry for this chunk is possibly invalid, or we've erased all runs up to the end of
			// the chunk. In this case, we need to re-scan the region following the end of our completed run
			// entries up to the chunk end point. We obtain only the first located run if there are other runs
			// remaining in this chunk, otherwise we obtain all runs up to the end of the chunk sample region.
			// Note that we do in fact need to re-scan up to the following run in the case of a partially
			// overlapping run, as it's possible the partially overlapping run entry we're about to cut intruded
			// a significant margin into the next run entry start position. If we don't re-evaluate the run that
			// follows it, it's possible we'll miss valid leading samples. Also note that we always need to
			// re-evaluate the first entry in a chunk when joining. Even if the last run stops before the first
			// run in this chunk, it's possible that a run was completed running into the sample region for this
			// chunk, and it's possible that samples within that region were included in creating a run within
			// this chunk that later got trimmed when it was completed, but initially influenced the run average
			// and affected its final position and length.
			CutRunEntry<SampleType> cutRunEntryForRescan;
			MinMaxWindowInfo<SampleType> minMaxWindowInfoForRescan = cutRunEntry.minMaxInfo;
			std::list<RunEntry<SampleType>> newRunEntries = ExtractRunsFromSampleData(sampleData, lastRunEntry.endSampleNo, chunkEndPos, minMaxWindowInfoForRescan, cutRunEntryForRescan, false, false, !scanToEndOfChunk);

			// If we scanned up to the end of this chunk, replace the cut run entry for the chunk.
			if (scanToEndOfChunk || newRunEntries.empty())
			{
				cutRunEntry = cutRunEntryForRescan;
				chunkRunEntries.clear();
				break;
			}

			// Merge the list of new runs generated in the overlapped region into the combined set of runs
			runEntries.splice(runEntries.end(), std::move(newRunEntries));
		}
	}

	// Merge the list of runs for this chunk into the combined set of runs
	runEntries.splice(runEntries.end(), std::move(chunkRunEntries));

	// It's possible completing a cut run entry from a previous chunk has completely enveloped this chunk, in
	// which case the cut run entry here might be covering duplicate space. To handle this case and ensure we
	// get the same results as a single-threaded run, we range test the last completed run entry with the cut
	// run entry for this chunk. If the last entry has an end position after the cut run start position, we
	// re-evaluate the runs following that position up to the end of this chunk, and replace the existing cut
	// run entry with the result.
	if (!runEntries.empty())
	{
		const auto& lastRunEntry = *runEntries.rbegin();
		if ((lastRunEntry.endSampleNo < chunkEndPos) && (lastRunEntry.endSampleNo > cutRunEntry.runEntry.startSampleNo))
		{
			// The first run entry for this chunk partially overlaps with our completed run entries. In this
			// case, we need to re-scan the region following the end of our completed run entries. Here we
			// calculate the region in which to re-scan.
			size_t startScanIndex = lastRunEntry.endSampleNo;
			size_t endScanIndex = chunkEndPos;

			// Obtain the new list of runs in the scan region
			MinMaxWindowInfo<SampleType> minMaxWindowInfoForRescan = cutRunEntry.minMaxInfo;
			CutRunEntry<SampleType> cutRunEntryForRescan;
			std::list<RunEntry<SampleType>> newRunEntries = ExtractRunsFromSampleData(sampleData, startScanIndex, endScanIndex, minMaxWindowInfoForRescan, cutRunEntryForRescan, false, false, false);

			// Merge the list of runs for this chunk into the combined set of runs
			runEntries.splice(runEntries.end(), std::move(newRunEntries));

			// Replace the cut run entry for the chunk
			cutRunEntry = cutRunEntryForRescan;
		}
	}

	// If there's another chunk of sample data after this one, complete the cut run at the end of this chunk,
	// and append any generated runs until we reach the start position of the next chunk. Note that it's
	// possible the cut run entry will actually resolve to before the end of the chunk, as it may pull in
	// samples from the next chunk to complete, then be trimmed back within the bounds of the current chunk. To
	// handle this, we loop until we produce a run which passes the end of the current chunk, or we reach the
	// end of the sample data.
	size_t lastRunEntryEndPos = (runEntries.empty() ? chunkStartPos : runEntries.rbegin()->endSampleNo);
	while (lastRunEntryEndPos < chunkEndPos)
	{
		// Attempt to resolve the supplied cut run into a run entry
		MinMaxWindowInfo<SampleType> minMaxWindowInfoForRescan = cutRunEntry.minMaxInfo;
		std::list<RunEntry<SampleType>> newRunEntries = ExtractRunsFromSampleData(sampleData, cutRunEntry.runEntry.endSampleNo, sampleEndNo, minMaxWindowInfoForRescan, cutRunEntry, true, false, true);
		if (newRunEntries.empty())
		{
			break;
		}
		runEntries.splice(runEntries.end(), std::move(newRunEntries));
		lastRunEntryEndPos = runEntries.rbegin()->endSampleNo;
	}
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::FilterRunEntriesToSyncCandidates(std::list<RunEntry<SampleType>>& runEntries) const
{
	auto runEntriesIterator = runEntries.begin();
	while (runEntriesIterator != runEntries.end())
	{
		// Obtain information on the current run entry
		const RunEntry<SampleType>& runEntry = *runEntriesIterator;
		double adjacentRunMaxDevianceInSamples = runEntry.initialMinMaxSampleRange * successiveRunAverageDifferenceTolerance;

		// Compare with the previous run entry
		if (runEntriesIterator != runEntries.begin())
		{
			auto previousRunEntryIterator = std::prev(runEntriesIterator, 1);
			const RunEntry<SampleType>& previousRunEntry = *previousRunEntryIterator;
			if ((previousRunEntry.averageSample - runEntry.averageSample) > adjacentRunMaxDevianceInSamples)
			{
				runEntries.erase(previousRunEntryIterator);
				continue;
			}
			else if ((runEntry.averageSample - previousRunEntry.averageSample) > adjacentRunMaxDevianceInSamples)
			{
				runEntriesIterator = runEntries.erase(runEntriesIterator);
				continue;
			}
		}

		// Compare with the next run entry
		auto nextRunEntryIterator = std::next(runEntriesIterator, 1);
		if (nextRunEntryIterator != runEntries.end())
		{
			const RunEntry<SampleType>& nextRunEntry = *nextRunEntryIterator;
			if ((nextRunEntry.averageSample - runEntry.averageSample) > adjacentRunMaxDevianceInSamples)
			{
				runEntries.erase(nextRunEntryIterator);
				continue;
			}
			else if ((runEntry.averageSample - nextRunEntry.averageSample) > adjacentRunMaxDevianceInSamples)
			{
				runEntriesIterator = runEntries.erase(runEntriesIterator);
				continue;
			}
		}

		// This run entry has passed, so advance to the next element.
		++runEntriesIterator;
	}
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::ErrorCorrectRunEntrySyncCandidates(const std::vector<SampleType>& sampleData, std::list<RunEntry<SampleType>>& runEntries) const
{
	auto runEntriesIterator = runEntries.begin();
	while (runEntriesIterator != runEntries.end())
	{
		//##FIX## Do something better with this
		RunEntry<SampleType>& runEntry = *runEntriesIterator;
		double ignoreSpikeToleranceInSamples = runEntry.initialMinMaxSampleRange * errorMergeRunAverageDifferenceTolerance;
		double adjacentRunMaxDevianceInSamples = runEntry.initialMinMaxSampleRange * successiveRunAverageDifferenceTolerance;

		// Compare with the previous run entry
		size_t runEntryLength = runEntry.endSampleNo - runEntry.startSampleNo;
		if (runEntriesIterator != runEntries.begin())
		{
			auto previousRunEntryIterator = std::prev(runEntriesIterator, 1);
			const RunEntry<SampleType>& previousRunEntry = *previousRunEntryIterator;
			size_t previousRunEntryLength = previousRunEntry.endSampleNo - previousRunEntry.startSampleNo;
			size_t longestRunEntry = std::max(runEntryLength, previousRunEntryLength);

			if (((runEntry.startSampleNo - previousRunEntry.endSampleNo) < longestRunEntry)
			 && (std::fabs(previousRunEntry.averageSample - runEntry.averageSample) <= adjacentRunMaxDevianceInSamples))
			{
				double averageIfMerged = std::accumulate(sampleData.cbegin() + previousRunEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo, 0.0) / (double)(runEntry.endSampleNo - previousRunEntry.startSampleNo);
				if (((runEntry.startSampleNo - previousRunEntry.endSampleNo) < ((runEntryLength + previousRunEntryLength) * errorMergeRunLengthTolerance))
				 || (std::fabs(previousRunEntry.averageSample - averageIfMerged) <= ignoreSpikeToleranceInSamples) || (std::fabs(runEntry.averageSample - averageIfMerged) <= ignoreSpikeToleranceInSamples))
				{
					runEntry.startSampleNo = previousRunEntry.startSampleNo;
					runEntry.averageSample = averageIfMerged;
					//##TODO## We currently don't use min/max from here onwards. Should we omit this step, and leave the
					//data incorrect?
					//runEntry.minValue = *std::min_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
					//runEntry.maxValue = *std::max_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
					runEntries.erase(previousRunEntryIterator);
					continue;
				}
			}
		}

		// Compare with the next run entry
		auto nextRunEntryIterator = std::next(runEntriesIterator, 1);
		if (nextRunEntryIterator != runEntries.end())
		{
			const RunEntry<SampleType>& nextRunEntry = *nextRunEntryIterator;
			size_t nextRunEntryLength = nextRunEntry.endSampleNo - nextRunEntry.startSampleNo;
			size_t longestRunEntry = std::max(runEntryLength, nextRunEntryLength);

			if (((nextRunEntry.startSampleNo - runEntry.endSampleNo) < longestRunEntry)
			 && (std::fabs(nextRunEntry.averageSample - runEntry.averageSample) <= adjacentRunMaxDevianceInSamples))
			{
				double averageIfMerged = std::accumulate(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + nextRunEntry.endSampleNo, 0.0) / (double)(nextRunEntry.endSampleNo - runEntry.startSampleNo);
				if (((nextRunEntry.startSampleNo - runEntry.endSampleNo) < ((runEntryLength + nextRunEntryLength) * errorMergeRunLengthTolerance))
				 || (std::fabs(nextRunEntry.averageSample - averageIfMerged) <= ignoreSpikeToleranceInSamples) || (std::fabs(runEntry.averageSample - averageIfMerged) <= ignoreSpikeToleranceInSamples))
				{
					runEntry.endSampleNo = nextRunEntry.endSampleNo;
					runEntry.averageSample = averageIfMerged;
					//runEntry.minValue = *std::min_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
					//runEntry.maxValue = *std::max_element(sampleData.cbegin() + runEntry.startSampleNo, sampleData.cbegin() + runEntry.endSampleNo);
					runEntries.erase(nextRunEntryIterator);
					continue;
				}
			}
		}

		// This run entry has passed, so advance to the next element.
		++runEntriesIterator;
	}
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::CleanRunEntryEdges(const std::vector<SampleType>& sampleData, std::list<RunEntry<SampleType>>& runEntries, size_t sampleStartNo, size_t sampleEndNo) const
{
	for (auto& runEntry : runEntries)
	{
		// Extend the beginning and end of the run entry to include extra samples that are within tolerance of the
		// average sample, or below it.
		double minMaxRunRangeDevianceInSamples = runEntry.initialMinMaxSampleRange * errorMergeRunAverageDifferenceTolerance;
		size_t currentSampleIndex = runEntry.endSampleNo + 1;
		while (currentSampleIndex < sampleEndNo)
		{
			SampleType currentSample = sampleData[currentSampleIndex++];
			if (((double)currentSample - runEntry.averageSample) > minMaxRunRangeDevianceInSamples)
			{
				break;
			}
			runEntry.endSampleNo = currentSampleIndex;
		}
		currentSampleIndex = runEntry.startSampleNo;
		while (currentSampleIndex > sampleStartNo)
		{
			SampleType currentSample = sampleData[--currentSampleIndex];
			if (((double)currentSample - runEntry.averageSample) > minMaxRunRangeDevianceInSamples)
			{
				break;
			}
			runEntry.startSampleNo = currentSampleIndex;
		}

		// Trim the edges of this run entry. Note that we don't trim samples for being too low here. We've already
		// filtered our runs down to sync pulses, and since our sync pulses are supposed to be the lowest signal level,
		// anything low is something we want to keep.
		size_t runLengthBeforeTrim = runEntry.endSampleNo - runEntry.startSampleNo;
		size_t runLengthMaxTrimSampleCount = (size_t)((double)runLengthBeforeTrim * runRangeMaxTrimLength);
		SampleType trimSampleTooHighThreshold = (SampleType)(runEntry.averageSample + (runEntry.initialMinMaxSampleRange * finalTrimTolerance));
		SampleType trimSampleTooLowThreshold = std::numeric_limits<SampleType>::min();
		TrimRunEntry(sampleData, runEntry, trimSampleTooHighThreshold, trimSampleTooLowThreshold, runLengthMaxTrimSampleCount);
	}
}

//----------------------------------------------------------------------------------------------------------------------
template<class SampleType>
void SyncDetector::TrimRunEntry(const std::vector<SampleType>& sampleData, RunEntry<SampleType>& runEntry, SampleType trimSampleTooHighThreshold, SampleType trimSampleTooLowThreshold, size_t maxTrimSampleCount) const
{
	// Trim from the start of the run
	size_t trimRunStartSampleNo = runEntry.startSampleNo;
	size_t trimRunEndSampleNo = trimRunStartSampleNo + maxTrimSampleCount;
	while (trimRunStartSampleNo < trimRunEndSampleNo)
	{
		// Obtain the next trim candidate
		SampleType trimSample = sampleData[trimRunStartSampleNo++];

		// If the next sample is outside the trim tolerance, cut it.
		if ((trimSample > trimSampleTooHighThreshold) || (trimSample < trimSampleTooLowThreshold))
		{
			runEntry.startSampleNo += 1;
			continue;
		}

		// If we've hit a sample we don't want to cut, abort any further processing.
		break;
	}

	// Trim from the end of the run
	trimRunStartSampleNo = runEntry.endSampleNo - 1;
	trimRunEndSampleNo = trimRunStartSampleNo - maxTrimSampleCount;
	while (trimRunStartSampleNo > trimRunEndSampleNo)
	{
		// Obtain the next trim candidate
		SampleType trimSample = sampleData[trimRunStartSampleNo--];

		// If the next sample is outside the trim tolerance, cut it.
		if ((trimSample > trimSampleTooHighThreshold) || (trimSample < trimSampleTooLowThreshold))
		{
			runEntry.endSampleNo -= 1;
			continue;
		}

		// If we've hit a sample we don't want to cut, abort any further processing.
		break;
	}
}
