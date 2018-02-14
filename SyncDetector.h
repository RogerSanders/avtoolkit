#ifndef __SYNCDETECTOR_H__
#define __SYNCDETECTOR_H__
#include <list>
#include <vector>
#include <map>

class SyncDetector
{
public:
	// Enumerations
	enum class SyncType
	{
		Horizontal,
		Vertical,
		Equalizer,
		Unknown,
	};

public:
	// Structures
	struct SyncPulseInfo
	{
		size_t startSampleNo;
		size_t endSampleNo;
		double averageSyncLevel;
		double approxMinMaxSampleRange;
	};
	struct SyncInfo
	{
		SyncType type;
		size_t startSampleNo;
		size_t endSampleNo;
		double averageSyncLevel;
		//double averageBlankingLevel;
		double approxMinMaxSampleRange;
	};

public:
	// Constructors
	SyncDetector();

	// Settings methods
	void RestoreDefaultSettings();

	// Sync detection methods
	template<class T>
	std::list<SyncPulseInfo> DetectSyncPulses(const std::vector<T>& sampleData, size_t initialSampleNo = 0, size_t sampleCount = 0, unsigned int threadCount = 0) const;
	template<class T>
	std::list<SyncInfo> DetectSyncEvents(const std::vector<T>& sampleData, const std::list<SyncPulseInfo>& syncPulses) const;

private:
	// Structures
	template <class T>
	struct MinMaxWindowInfo
	{
		bool minMaxValuesPopulated;
		double minMaxSampleRange;
		//##DEBUG##
		T minSampleValue;
		T maxSampleValue;

		bool slidingWindowEnabled;
		size_t scanSampleStartNo;
		size_t scanSampleEndNo;
		size_t scanStartPos;
		size_t scanResetPos;
	};
	template <class T>
	struct RunEntry
	{
		size_t startSampleNo;
		size_t endSampleNo;
		T minValue;
		T maxValue;
		double initialMinMaxSampleRange;
		double averageSample;
	};
	template <class T>
	struct CutRunEntry
	{
		RunEntry<T> runEntry;
		MinMaxWindowInfo<T> minMaxInfo;
	};

private:
	// Sync detection methods
	template<class T>
	void ObtainMinMaxValues(const std::vector<T>& sampleData, size_t samplePos, MinMaxWindowInfo<T>& minMaxInfo) const;
	template<class T>
	std::list<RunEntry<T>> ExtractRunsFromSampleData(const std::vector<T>& sampleData, size_t sampleStartNo, size_t sampleEndNo, MinMaxWindowInfo<T>& minMaxInfo, CutRunEntry<T>& cutRunEntry, bool cutRunEntryPopulated, bool resolveCutRunEntryOnly, bool getNextRunEntryOnly) const;
	template<class T>
	void MergeChunkResults(const std::vector<T>& sampleData, size_t sampleEndNo, std::list<RunEntry<T>>& runEntries, size_t chunkStartPos, size_t chunkEndPos, std::list<RunEntry<T>>& chunkRunEntries, CutRunEntry<T>& cutRunEntry, MinMaxWindowInfo<T>& minMaxInfo) const;
	template<class T>
	void FilterRunEntriesToSyncCandidates(std::list<RunEntry<T>>& runEntries) const;
	template<class T>
	void ErrorCorrectRunEntrySyncCandidates(const std::vector<T>& sampleData, std::list<RunEntry<T>>& runEntries) const;
	template<class T>
	void CleanRunEntryEdges(const std::vector<T>& sampleData, std::list<RunEntry<T>>& runEntries, size_t sampleStartNo, size_t sampleEndNo) const;
	template<class T>
	void TrimRunEntry(const std::vector<T>& sampleData, RunEntry<T>& runEntry, T trimSampleTooHighThreshold, T trimSampleTooLowThreshold, size_t maxTrimSampleCount) const;
	size_t GetMostCommonSyncPulseLength(const std::map<size_t, size_t>& syncPulseLengthCounts) const;

private:
	//##FIX## We've seen as low as 0.125 for the Mega Drive here. Perhaps we should be using the front porch, and
	//relying on averaging to cancel out the colour burst? This doesn't work with a quantum signal though, as unlike an
	//analog signal we've got sampling bias potentially causing us to have a much higher or lower average sample value.
	//We'd have to fit the sample points to a curve and take a difference of the volumes below and above in order to
	//avoid this issue. Then again, higher sample rates avoid this problem, and we'll probably get a better result from
	//a longer averaged run than an often noisy and short trailing run from the last line.
	const double blankingSampleRatioToSyncWidth = 0.2;

//##DEBUG##
public:
	double runRangeValueDeviance;
	double runRangeTrimTolerance;
	double runRangeMaxTrimLength;
	double finalTrimTolerance;
	unsigned long long minRunRangeSampleLength;
	double successiveRunAverageDifferenceTolerance;
	double errorMergeRunAverageDifferenceTolerance;
	double errorMergeRunLengthTolerance;
	bool enableMinMaxSlidingWindow;
	double minMaxSlidingWindowLength;
	//##TODO## Add parameters to force a min/max sample range as percentage (0.0 - 1.0) of sample range (std::numeric_limits),
	//and skip all min/max calculation in this case.
};

#include "SyncDetector.inl"
#endif
