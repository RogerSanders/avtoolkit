#ifndef __FRAMEBUILDER_H__
#define __FRAMEBUILDER_H__
#include "Logger.h"
#include "SyncDetector.h"
#include <list>

class FrameBuilder
{
public:
	// Structures
	struct ColorBurstWaveInfo
	{
		double startPos;
		double endPos;
		double peakLevel;
		bool isPositive;
		bool isPredicted;
	};
	struct LineInfo
	{
		SyncDetector::SyncInfo leadingSyncInfo;
		SyncDetector::SyncInfo followingSyncInfo;
		double backPorchStartPos;
		double frontPorchEndPos;
		double averageSyncLevel;
		double averageBlankingLevel;
		//##TODO##
		bool halfLine;
		double ireLevel0;
		double ireLevel100;
		bool colorBurstValid;
		double colorBurstSyncCorrection;
		std::vector<ColorBurstWaveInfo> colorBurstWaves;
	};
	struct FieldInfo
	{
		std::vector<SyncDetector::SyncInfo> syncEvents;
		int lineCount;
		size_t endSampleNo;
		SyncDetector::SyncInfo followingSyncEvent;
		std::vector<LineInfo> lines;
	};
	struct FrameInfo
	{
		std::vector<FieldInfo> fields;
	};

public:
	// Constructors
	FrameBuilder(const Logger& logger);

	// Frame detection methods
	template<class SampleType>
	std::list<FieldInfo> DetectFields(const std::vector<SampleType>& sampleData, const std::list<SyncDetector::SyncInfo>& syncEvents) const;
	std::vector<FrameInfo> DetectFrames(const std::list<FieldInfo>& fields) const;
	template<class SampleType>
	void DetectLines(const std::vector<SampleType>& sampleData, FieldInfo& fieldInfo) const;
	template<class SampleType>
	void DetectLines(const std::vector<SampleType>& sampleData, std::vector<FrameInfo>& frames, unsigned int threadCount = 0) const;

private:
	// Trigger detection methods
	template<class SampleType>
	double FindSyncRisingEdgeSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo) const;
	template<class SampleType>
	double FindSyncFallingEdgeSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo) const;
	template<class SampleType>
	double FindSyncRisingEdgeEndSamplePos(const std::vector<SampleType>& inputData, const SyncDetector::SyncInfo& syncInfo, double risingEdgePos) const;

	// Colour burst methods
	template<class SampleType>
	void DetectColorBurst(const std::vector<SampleType>& backPorchData, SampleType zeroLevel, SampleType burstAmplitude, std::vector<ColorBurstWaveInfo>& burstWaves) const;
	bool ValidateColorBurst(unsigned int minimumWaveCount, const std::vector<ColorBurstWaveInfo>& burstWaves) const;
	template<class SampleType>
	bool RepairColorBurst(const std::vector<SampleType>& backPorchData, SampleType zeroLevel, SampleType burstAmplitude, std::vector<ColorBurstWaveInfo>& burstWaves) const;
	bool PerformColorBurstLineSyncCorrection(FieldInfo& fieldInfo) const;

	// IRE conversion methods
	template<class SampleType>
	float SampleToIRE(SampleType sampleValue, float ireLevel0, float ireLevel100) const;
	template<class SampleType>
	SampleType IREToSample(float ire, float ireLevel0, float ireLevel100) const;

	// Math helper methods
	template<class Iter>
	static double FindMedianValue(Iter first, Iter last);

private:
	const Logger& _logger;

//##DEBUG##
public:
	bool combineInterlacedFields;
	bool interlaceHalfLineInFirstField;
	double syncAmplitudeMinTolerance;
	double slopeDetectionTolerance;
	double slopeValueFlatTolerance;
	double syncLengthToBackPorchMinRatio;
	double blankingUpsampleRatio;
	double colorBurstIREDetectionThreshold;
	unsigned int colorBurstMinimumHalfOscillationCount;
	double colorBurstLineSyncTolerance;
	bool useColorBurstForLineSyncCorrection;
	double colorBurstRepairWaveLengthTolerance;
};

#include "FrameBuilder.inl"
#endif
