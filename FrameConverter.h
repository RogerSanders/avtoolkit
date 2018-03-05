#ifndef __FRAMECONVERTER_H__
#define __FRAMECONVERTER_H__
#include "Logger.h"
#include "FileSystemInterop.h"
#include "FrameBuilder.h"
#include "SyncDetector.h"
#include <string>
#include <vector>

class FrameConverter
{
public:
	// Constructors
	FrameConverter(const Logger& logger);

	// Frame conversion methods
	template<class SampleType>
	void WriteFrameToBMP(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo, unsigned int pixelsPerLine) const;

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

//##DEBUG##
public:
	const Logger& _logger;
	double blankingPercentage;
	double blankingLeadingPercentage;
};

#include "FrameConverter.inl"
#endif
