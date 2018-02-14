#ifndef __FRAMECONVERTER_H__
#define __FRAMECONVERTER_H__
#include "FileSystemInterop.h"
#include "FrameBuilder.h"
#include "SyncDetector.h"
#include <string>
#include <vector>
#include "spline.h"

class FrameConverter
{
public:
	// Constructors
	FrameConverter();

	// Frame conversion methods
	template<class SampleType>
	void WriteFrameToFile(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo, unsigned int pixelsPerLine) const;

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

private:
	// Line analysis methods
	void FindPreciseLineStartEndPos(const SyncDetector::SyncInfo& leadingSyncInfo, const SyncDetector::SyncInfo& followingSyncInfo, const tk::spline& lineSpline, double& preciseLineStartPos, double& preciseLineEndPos) const;
	//##TODO##
	//template<class SampleType>
	//bool DetectSyncBurst();

	// IRE conversion methods
	template<class SampleType>
	float SampleToIRE(SampleType sampleValue, float ireLevel0, float ireLevel100) const;
	template<class SampleType>
	SampleType IREToSample(float ire, float ireLevel0, float ireLevel100) const;

	// Resampling methods
	template<class SampleType>
	void ResampleLinear(const std::vector<SampleType>& inputData, std::vector<SampleType>& outputData) const;
	template<class SampleType>
	void ResampleCubic(const std::vector<SampleType>& inputData, std::vector<SampleType>& outputData) const;

	//##DEBUG##
public:
	double blankingPercentage;
	double blankingLeadingPercentage;
	double syncSearchPosIncrement;
	double syncAmplitudeMinTolerance;
	double slopeValueFlatToleranceAsPercentage;
};

#include "FrameConverter.inl"
#endif
