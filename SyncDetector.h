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
		double approxMinMaxSampleRange;
	};

public:
	// Constructors
	SyncDetector();

	// Settings methods
	void RestoreDefaultSettings();

	// Sync detection methods
	template<class SampleType>
	std::list<SyncPulseInfo> DetectSyncPulses(const std::vector<SampleType>& sampleData, size_t initialSampleNo = 0, size_t sampleCount = 0, unsigned int threadCount = 0) const;
	template<class SampleType>
	std::list<SyncInfo> DetectSyncEvents(const std::vector<SampleType>& sampleData, const std::list<SyncPulseInfo>& syncPulses) const;

private:
	// Structures
	template <class SampleType>
	struct MinMaxWindowInfo
	{
		bool minMaxValuesPopulated;
		double minMaxSampleRange;
		//##DEBUG##
		SampleType minSampleValue;
		SampleType maxSampleValue;

		bool slidingWindowEnabled;
		size_t scanSampleStartNo;
		size_t scanSampleEndNo;
		size_t scanStartPos;
		size_t scanResetPos;
	};
	template <class SampleType>
	struct RunEntry
	{
		size_t startSampleNo;
		size_t endSampleNo;
		SampleType minValue;
		SampleType maxValue;
		double initialMinMaxSampleRange;
		double averageSample;
	};
	template <class SampleType>
	struct CutRunEntry
	{
		RunEntry<SampleType> runEntry;
		MinMaxWindowInfo<SampleType> minMaxInfo;
	};

private:
	// Sync detection methods
	template<class SampleType>
	void ObtainMinMaxValues(const std::vector<SampleType>& sampleData, size_t samplePos, MinMaxWindowInfo<SampleType>& minMaxInfo) const;
	template<class SampleType>
	std::list<RunEntry<SampleType>> ExtractRunsFromSampleData(const std::vector<SampleType>& sampleData, size_t sampleStartNo, size_t sampleEndNo, MinMaxWindowInfo<SampleType>& minMaxInfo, CutRunEntry<SampleType>& cutRunEntry, bool cutRunEntryPopulated, bool resolveCutRunEntryOnly, bool getNextRunEntryOnly) const;
	template<class SampleType>
	void MergeChunkResults(const std::vector<SampleType>& sampleData, size_t sampleEndNo, std::list<RunEntry<SampleType>>& runEntries, size_t chunkStartPos, size_t chunkEndPos, std::list<RunEntry<SampleType>>& chunkRunEntries, CutRunEntry<SampleType>& cutRunEntry, MinMaxWindowInfo<SampleType>& minMaxInfo) const;
	template<class SampleType>
	void FilterRunEntriesToSyncCandidates(std::list<RunEntry<SampleType>>& runEntries) const;
	template<class SampleType>
	void ErrorCorrectRunEntrySyncCandidates(const std::vector<SampleType>& sampleData, std::list<RunEntry<SampleType>>& runEntries) const;
	template<class SampleType>
	void CleanRunEntryEdges(const std::vector<SampleType>& sampleData, std::list<RunEntry<SampleType>>& runEntries, size_t sampleStartNo, size_t sampleEndNo) const;
	template<class SampleType>
	void TrimRunEntry(const std::vector<SampleType>& sampleData, RunEntry<SampleType>& runEntry, SampleType trimSampleTooHighThreshold, SampleType trimSampleTooLowThreshold, size_t maxTrimSampleCount) const;
	size_t GetMostCommonSyncPulseLength(const std::map<size_t, size_t>& syncPulseLengthCounts) const;

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
	double syncPulseLengthVariationTolerance;
	//##TODO## Add parameters to force a min/max sample range as percentage (0.0 - 1.0) of sample range (std::numeric_limits),
	//and skip all min/max calculation in this case.
};

#include "SyncDetector.inl"
#endif
