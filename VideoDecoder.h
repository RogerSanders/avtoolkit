#ifndef __VIDEODECODER_H__
#define __VIDEODECODER_H__
#include "Logger.h"
#include "FileSystemInterop.h"
#include "FrameBuilder.h"
#include "SyncDetector.h"
#include <string>
#include <vector>

class VideoDecoder
{
public:
	// Constructors
	VideoDecoder(const SyncDetector& syncDetector, const FrameBuilder& frameBuilder, const Logger& log);

	// Video conversion methods
	template<class SampleType>
	bool ConvertCompositeVideoToImages(const PathString& inputFilePath, const PathString& outputFolderPath, const PathString& outputFileNameBase) const;
	template<class SampleType>
	void WriteFramesToBMP(const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::vector<FrameBuilder::FrameInfo>& frames, size_t initialFrameNo, unsigned int threadCount = 0) const;
	template<class SampleType>
	bool WriteFrameToBMP(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo) const;

private:
	// Enumerations
	enum class BitmapCompressionType : unsigned int
	{
		RGB = 0,
		RLE8 = 1,
		RLE4 = 2,
		Bitfields = 3,
		JPEG = 4,
		PNG = 5,
	};

private:
	// Structures
	struct BitmapFileHeader
	{
		unsigned short bfType;
		unsigned int bfSize;
		short bfReserved1;
		short bfReserved2;
		unsigned int bfOffBits;
	};
	struct BitmapInfoHeader
	{
		unsigned int biSize;
		int biWidth;
		int biHeight;
		short biPlanes;
		short biBitCount;
		BitmapCompressionType biCompression;
		unsigned int biSizeImage;
		int biXPelsPerMeter;
		int biYPelsPerMeter;
		unsigned int biClrUsed;
		unsigned int biClrImportant;
	};
	struct DecodeLinePixelDataBuffer
	{
		std::vector<double> carrierWaveISamplePoints;
		std::vector<double> carrierWaveQSamplePoints;
		std::vector<double> samplePointsQ;
		std::vector<double> samplePointsI;
		std::vector<double> samplePointsY;
		std::vector<double> lineResolutionOutputR;
		std::vector<double> lineResolutionOutputG;
		std::vector<double> lineResolutionOutputB;
	};

private:
	// Video conversion methods
	bool GetWavePeakPositionsForLine(const FrameBuilder::LineInfo& lineInfo, std::vector<double>& wavePeakPositions, bool& firstWavePeakIsPositive) const;
	void CalculateIRELevelsForField(const FrameBuilder::FieldInfo& fieldInfo, double& ireLevel0ForField, double& ireLevel100ForField) const;
	double CalculateColorBurstWaveFrequencyForLine(const FrameBuilder::LineInfo& lineInfo, std::vector<double>& wavePeakPositions) const;
	double CalculateColorBurstWaveFrequencyForField(const FrameBuilder::FieldInfo& fieldInfo) const;
	void PhaseLockColorBurstSamples(std::vector<double> wavePeakPositions, bool firstWavePeakIsPositive, std::vector<double> targetWavePeakPositions, bool firstTargetWavePeakIsPositive, double burstWaveFrequency, size_t& bestMatchPeakPositionIndex, bool& bestMatchPeakPositionIsPositive, double& phaseLockDisplacement) const;
	static void ConvertYIQToRGB(double sampleY, double sampleI, double sampleQ, double& red, double& green, double& blue);
	static void ConvertYUVToRGB(double sampleY, double sampleU, double sampleV, double& red, double& green, double& blue);
	template<class SampleType>
	void DecodeMonochromeLinePixelDataForActiveScanRegion(const std::vector<SampleType>& sampleData, const FrameBuilder::FieldInfo& fieldInfo, const FrameBuilder::LineInfo& lineInfo, size_t lineNo, double preciseLineStartPos, double preciseLineEndPos, double ireLevel0ForField, double ireLevel100ForField, unsigned int outputPixelCount, std::vector<double>& outputDataR, std::vector<double>& outputDataG, std::vector<double>& outputDataB, DecodeLinePixelDataBuffer& buffer) const;
	template<class SampleType>
	bool DecodeColorLinePixelDataForActiveScanRegion(const std::vector<SampleType>& sampleData, const FrameBuilder::FieldInfo& fieldInfo, const FrameBuilder::LineInfo& lineInfo, size_t lineNo, double preciseLineStartPos, double preciseLineEndPos, double burstWaveFrequencyForField, double ireLevel0ForField, double ireLevel100ForField, unsigned int outputPixelCount, std::vector<double>& outputDataR, std::vector<double>& outputDataG, std::vector<double>& outputDataB, DecodeLinePixelDataBuffer& buffer) const;

private:
	const Logger& _log;
	const SyncDetector& _syncDetector;
	const FrameBuilder& _frameBuilder;

	//##DEBUG##
public:
	double blankingPercentage;
	double blankingLeadingPercentage;
	unsigned int lineWidthInPixels;
	size_t maxChunkSizeInBytes;
	size_t burstWaveSkipCount;
	bool useAverageFieldBlankingLevel;
	bool useAverageFieldBurstWaveFrequency;
	bool matchColorBurstPhaseBetweenLines;
	bool rawOutputOnly;
	bool decodeColor;
	bool useIRE7Point5;
	bool decodeAsYUV;
	double waveFrequencyTolerance;
	bool forceMonoOutputForColorDecoding;
	bool forceColorBurstWaveFrequency;
	double forcedColorBurstWaveFrequency;
};

#include "VideoDecoder.inl"
#endif
