#include "SplineHelpers.h"

//----------------------------------------------------------------------------------------
// Frame conversion methods
//----------------------------------------------------------------------------------------
template<class SampleType>
void FrameConverter::WriteFrameToBMP(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo, unsigned int pixelsPerLine) const
{
	// Determine the total number of lines in this frame
	unsigned int lineCount = 0;
	for (const auto& fieldInfo : frameInfo.fields)
	{
		lineCount += (unsigned int)fieldInfo.lines.size();
	}

	// Calculate the amount of padding on each line. Lines are padded out to DWORD boundaries.
	unsigned int bitsPerPixel = 24;
	unsigned int lineByteCount = pixelsPerLine * 3;
	unsigned int linePaddingByteCount = 0;
	if((lineByteCount % sizeof(unsigned int)) != 0)
	{
		linePaddingByteCount = sizeof(unsigned int) - (lineByteCount % sizeof(unsigned int));
	}

	// Calculate the sizes of our header structures. We need to do this as rather than rely on the structure size in
	// memory, as padding can interfere with the layout.
	BitmapFileHeader fileHeader = {};
	BitmapInfoHeader bitmapHeader = {};
	unsigned int bitmapFileHeaderSize = sizeof(fileHeader.bfType) + sizeof(fileHeader.bfSize) + sizeof(fileHeader.bfReserved1) + sizeof(fileHeader.bfReserved2) + sizeof(fileHeader.bfOffBits);
	unsigned int bitmapInfoHeaderSize = sizeof(bitmapHeader.biSize) + sizeof(bitmapHeader.biWidth) + sizeof(bitmapHeader.biHeight) + sizeof(bitmapHeader.biPlanes) + sizeof(bitmapHeader.biBitCount) + sizeof(bitmapHeader.biCompression) + sizeof(bitmapHeader.biSizeImage) + sizeof(bitmapHeader.biXPelsPerMeter) + sizeof(bitmapHeader.biYPelsPerMeter) + sizeof(bitmapHeader.biClrUsed) + sizeof(bitmapHeader.biClrImportant);

	// Calculate the totals and offsets we need to write to the file
	unsigned int pixelDataOffset = bitmapFileHeaderSize + bitmapInfoHeaderSize;
	unsigned int pixelDataSize = (unsigned int)((lineByteCount + linePaddingByteCount) * lineCount);
	unsigned int fileSize = pixelDataOffset + pixelDataSize;

	// Write the bitmap file header to the file
	char typeByte1 = 'B';
	char typeByte2 = 'M';
	fileHeader.bfType = ((unsigned short)typeByte2 << 8) | (unsigned short)typeByte1;
	fileHeader.bfSize = fileSize;
	fileHeader.bfReserved1 = 0;
	fileHeader.bfReserved2 = 0;
	fileHeader.bfOffBits = pixelDataOffset;

	// Write the bitmap info header to the file
	bitmapHeader.biSize = bitmapInfoHeaderSize;
	bitmapHeader.biWidth = (int)pixelsPerLine;
	bitmapHeader.biHeight = -((int)lineCount);
	bitmapHeader.biPlanes = 1;
	bitmapHeader.biBitCount = (short)bitsPerPixel;
	bitmapHeader.biCompression = BitmapCompressionType::RGB;
	bitmapHeader.biSizeImage = (int)fileSize;
	bitmapHeader.biXPelsPerMeter = 0;
	bitmapHeader.biYPelsPerMeter = 0;
	bitmapHeader.biClrUsed = 0;
	bitmapHeader.biClrImportant = 0;

	// Create an output file for this frame
	std::ofstream outfile(outputFilePath, std::ios_base::binary);
	outfile.write((char*)&fileHeader.bfType, sizeof(fileHeader.bfType));
	outfile.write((char*)&fileHeader.bfSize, sizeof(fileHeader.bfSize));
	outfile.write((char*)&fileHeader.bfReserved1, sizeof(fileHeader.bfReserved1));
	outfile.write((char*)&fileHeader.bfReserved2, sizeof(fileHeader.bfReserved2));
	outfile.write((char*)&fileHeader.bfOffBits, sizeof(fileHeader.bfOffBits));
	outfile.write((char*)&bitmapHeader.biSize, sizeof(bitmapHeader.biSize));
	outfile.write((char*)&bitmapHeader.biWidth, sizeof(bitmapHeader.biWidth));
	outfile.write((char*)&bitmapHeader.biHeight, sizeof(bitmapHeader.biHeight));
	outfile.write((char*)&bitmapHeader.biPlanes, sizeof(bitmapHeader.biPlanes));
	outfile.write((char*)&bitmapHeader.biBitCount, sizeof(bitmapHeader.biBitCount));
	outfile.write((char*)&bitmapHeader.biCompression, sizeof(bitmapHeader.biCompression));
	outfile.write((char*)&bitmapHeader.biSizeImage, sizeof(bitmapHeader.biSizeImage));
	outfile.write((char*)&bitmapHeader.biXPelsPerMeter, sizeof(bitmapHeader.biXPelsPerMeter));
	outfile.write((char*)&bitmapHeader.biYPelsPerMeter, sizeof(bitmapHeader.biYPelsPerMeter));
	outfile.write((char*)&bitmapHeader.biClrUsed, sizeof(bitmapHeader.biClrUsed));
	outfile.write((char*)&bitmapHeader.biClrImportant, sizeof(bitmapHeader.biClrImportant));

	// Write the pixel data to the file
	std::vector<SampleType> leadingData;
	std::vector<SampleType> activeScanData;
	std::vector<SampleType> followingData;
	std::vector<unsigned char> outputData(pixelsPerLine*3);
	float sampleConversionFactor = (float)((double)std::numeric_limits<unsigned char>::max() / ((double)std::numeric_limits<SampleType>::max() + -(double)std::numeric_limits<SampleType>::min()));
	for (const FrameBuilder::FieldInfo& fieldInfo : frameInfo.fields)
	{
		for (const FrameBuilder::LineInfo& lineInfo : fieldInfo.lines)
		{
			double preciseLineStartPos = lineInfo.backPorchStartPos;
			double preciseLineEndPos = lineInfo.frontPorchEndPos;
			double preciseLineWidth = (preciseLineEndPos - preciseLineStartPos);

			unsigned int outputSampleCountForLine = (unsigned int)(outputData.size() / 3);
			unsigned int blankingLeadingSampleCount = (unsigned int)(blankingPercentage * blankingLeadingPercentage * (double)outputSampleCountForLine);
			unsigned int blankingFollowingSampleCount = (unsigned int)(blankingPercentage * (double)outputSampleCountForLine) - blankingLeadingSampleCount;
			//##FIX## Make the FindPreciseLineStartEndPos method handle alignment for non-hsync events.
			if (lineInfo.leadingSyncInfo.type != SyncDetector::SyncType::Horizontal)
			{
				preciseLineStartPos = (double)lineInfo.leadingSyncInfo.startSampleNo;
				preciseLineEndPos = (double)lineInfo.followingSyncInfo.startSampleNo;
				preciseLineWidth = preciseLineEndPos - preciseLineStartPos;
				blankingLeadingSampleCount = 0;
			}
			unsigned int activeImageSampleCount = outputSampleCountForLine - (blankingLeadingSampleCount + blankingFollowingSampleCount);
			unsigned int outputDataPos = 0;

			//##FIX## Make CubicInterpolateCatmullRom accept offsets into an existing buffer, and use a combined buffer
			//here.
			//##FIX## Replace all our uses of raw pointers with vectors and indexes
			leadingData.resize(blankingLeadingSampleCount);
			CubicInterpolateCatmullRom(sampleData.data(), (double)lineInfo.leadingSyncInfo.startSampleNo, preciseLineStartPos, leadingData);
			activeScanData.resize(activeImageSampleCount);
			CubicInterpolateCatmullRom(sampleData.data(), preciseLineStartPos, preciseLineEndPos, activeScanData);
			followingData.resize(blankingFollowingSampleCount);
			CubicInterpolateCatmullRom(sampleData.data(), preciseLineEndPos, lineInfo.followingSyncInfo.startSampleNo + ((lineInfo.followingSyncInfo.endSampleNo - lineInfo.followingSyncInfo.startSampleNo) * (1.0 - blankingLeadingPercentage)), followingData);

			for (unsigned int i = 0; i < leadingData.size(); ++i)
			{
				SampleType sampleValue = leadingData[i];
				unsigned char sampleValueAsByte = (unsigned char)(((float)sampleValue + -(float)std::numeric_limits<SampleType>::min()) * sampleConversionFactor);
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
			}
			for (unsigned int i = 0; i < activeScanData.size(); ++i)
			{
				SampleType sampleValue = activeScanData[i];
				unsigned char sampleValueAsByte = (unsigned char)(((float)sampleValue + -(float)std::numeric_limits<SampleType>::min()) * sampleConversionFactor);
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
			}
			for (unsigned int i = 0; i < followingData.size(); ++i)
			{
				SampleType sampleValue = followingData[i];
				unsigned char sampleValueAsByte = (unsigned char)(((float)sampleValue + -(float)std::numeric_limits<SampleType>::min()) * sampleConversionFactor);
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
				outputData[outputDataPos++] = sampleValueAsByte;
			}

			outfile.write((char*)&outputData[0], outputData.size());
			for (unsigned int i = 0; i < linePaddingByteCount; ++i)
			{
				outfile.write("", 1);
			}
		}
	}
}
