#include "SyncDetector.h"

//----------------------------------------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------------------------------------
SyncDetector::SyncDetector(const Logger& logger)
:_logger(logger)
{
	// Set our initial setting values
	RestoreDefaultSettings();
}

//----------------------------------------------------------------------------------------------------------------------
// Settings methods
//----------------------------------------------------------------------------------------------------------------------
void SyncDetector::RestoreDefaultSettings()
{
	// These values modify the tolerance of run detection when identifying sequences of similar values. These tolerance
	// values are relative to the detected scale of the data itself, and aren't specific to any particular signal
	// structure. If there is a large amount of noise in the data however, these values may need to be increased to
	// improve sync detection.
	runRangeValueDeviance = 0.1;
	runRangeTrimTolerance = 0.5;

	// This value determines how many leading or trailing samples are permitted to be trimmed from a run. The default
	// value should be sufficient for any signal.
	runRangeMaxTrimLength = 0.2;

	// This value controls the sample tolerance from the average sample value in our final merged runs when trimming is
	// being performed. Lower values more aggressively trim outlying samples from the edges of runs, while higher values
	// allow more samples to be maintained. The default value here should be sufficient for any signal.
	finalTrimTolerance = 0.03;

	// This value determines how many consecutive samples must meet the run value deviance test for a run to be formed.
	// Although it is based on sample counts, it shouldn't need adjustment with higher or lower resolution samples. High
	// resolution samples will generate more runs that pass this test, but they'll still be merged together or discarded
	// during the process correctly. We don't care about any runs less than 20 samples long for video detection, as none
	// of the triggers we're looking for can be shorter than this while still having a signal with enough resolution to
	// hold useful content.
	minRunRangeSampleLength = 20;

	// When folding successive runs together, this determines how similar the runs must be in average value in order to
	// combine them. The default value here should be sufficient for any signal.
	successiveRunAverageDifferenceTolerance = 0.1;

	// When blocks of non-conforming samples are located between two compatible runs, where at least one of those runs
	// is longer than the non-conforming sample block, this constant adjusts how much the non-conforming samples are
	// allowed to influence the average value of the resulting combined run. The default value here should be
	// sufficient, as it allows shorter spikes with large deviation through, while requiring longer blocks to be closer
	// to the trend of the run they're merged with. Note in particular that an incorrectly embedded colour burst during
	// sync is well handled by a small tolerance, as it is a sine wave that averages out. If a signal contains a large
	// amount of noise, due to laser rot or other source issues, a higher tolerance value here might improve sync
	// detection. Care should be taken not to make this value too large, as increasing the tolerance here makes it more
	// likely that equalization pulses within the vsync region will be folded into vsync.
	errorMergeRunAverageDifferenceTolerance = 0.007;

	// When blocks of non-conforming samples are located between two compatible runs of significantly longer length,
	// this constant adjusts how long a block of non-conforming samples needs to be in order to allow merging it into
	// the run, regardless of how it affects the average sample value of the resulting run. This is useful to allow
	// sharp spikes or other strong deviations in a localized space to be folded into a run, independent of the average
	// calculation. The total length of the runs on either side are summed, and multiplied by this value. If the length
	// of the non-conforming data is less than or equal to the result, it is allowed to be folded into the run. The
	// default here is carefully chosen for a specific property, which is that it doesn't allow equalizing pulses within
	// the vsync region to be folded in. A value of 0.07 reaches the point where this might start to occur. We still
	// allow a relatively long length of data to be combined however, to improve tolerance for unexpected signal
	// variance during the sync period.
	errorMergeRunLengthTolerance = 0.05;

	// These values control the min/max sliding window. Composite video can be at any base amplitude and scale, and the
	// sliding window allows us to scale our tolerance tests when searching for sync markers dynamically, according to
	// the base and peak amplitude of the video content nearby. If a sliding window isn't used, the entire input is
	// scanned to determine the min/max values from the outset. Calculating a single min/max for the entire input is
	// faster, and gives more consistent results if the video signal is output at a relatively stable level. If your
	// signal varies in amplitude however, or if it is captured in "a/c" mode and shifts in base amplitude when large
	// shifts in intensity take place (IE, going from all black to all white), a sliding window will generally give
	// better results in these cases. The sliding window length can be adjusted to constrain it to a more localized
	// region if required. Since we don't know or assume anything about the video sample rate (IE, it could be 30MSPS or
	// 30GSPS), this value is based on the total number of samples. It may be more appropriate to scale this value to be
	// relative to real time, such as using a window size of 1 second real time.
	enableMinMaxSlidingWindow = false;
	minMaxSlidingWindowLength = 0.01;

	// This value defines how much the length of a given sync pulse is allowed to vary from other sync events in the
	// data stream while still being considered valid.
	syncPulseLengthVariationTolerance = 0.25;
}

//----------------------------------------------------------------------------------------------------------------------
// Sync detection methods
//----------------------------------------------------------------------------------------------------------------------
size_t SyncDetector::GetMostCommonSyncPulseLength(const std::map<size_t, size_t>& syncPulseLengthCounts) const
{
	size_t mostCommonPulseLength = 0;
	size_t mostCommonPulseLengthSampleCount = 0;
	for (auto lengthCountEntry : syncPulseLengthCounts)
	{
		if (lengthCountEntry.second > mostCommonPulseLengthSampleCount)
		{
			mostCommonPulseLengthSampleCount = lengthCountEntry.second;
			mostCommonPulseLength = lengthCountEntry.first;
		}
	}
	return mostCommonPulseLength;
}
