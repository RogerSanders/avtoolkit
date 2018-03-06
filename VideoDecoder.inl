#include "FrameBuilder.h"
#include "SplineHelpers.h"
#include "SyncDetector.h"
#include <fstream>

//----------------------------------------------------------------------------------------
// Video conversion methods
//----------------------------------------------------------------------------------------
template<class SampleType>
bool VideoDecoder::ConvertCompositeVideoToImages(const PathString& inputFilePath, const PathString& outputFolderPath, const PathString& outputFileNameBase) const
{
	// Write a log message about the decode operation
	_logger.Info("Decoding \"{0}\"", inputFilePath);

	// Attempt to open the target file, and calculate the size of the file in bytes.
	std::ifstream infile(inputFilePath, std::ios_base::binary);
	if (infile.fail())
	{
		_logger.Error("Failed to open input file \"{0}\"", inputFilePath);
		return false;
	}
	infile.seekg(0, infile.end);
	size_t fileLengthInBytes = infile.tellg();
	infile.seekg(0, infile.beg);

	// Process the input file and write detected frames out to disk
	std::vector<SampleType> fileData;
	size_t currentFilePos = 0;
	bool firstPass = true;
	size_t sampleCountInChunk = 0;
	size_t continuePos;
	size_t startFrameNo = 0;
	while (currentFilePos < fileLengthInBytes)
	{
		// Calculate the size of the next chunk of file data to read
		size_t chunkSizeInBytes = (fileLengthInBytes - currentFilePos);
		chunkSizeInBytes = (chunkSizeInBytes > maxChunkSizeInBytes) ? maxChunkSizeInBytes : chunkSizeInBytes;
		sampleCountInChunk = (chunkSizeInBytes / sizeof(SampleType));

		// Prepare the memory buffer to receive file data
		size_t fileReadStartPos;
		size_t sampleScanningStartPos;
		if (firstPass)
		{
			fileData.resize(sampleCountInChunk);
			sampleScanningStartPos = 0;
			fileReadStartPos = 0;
			firstPass = false;
		}
		else
		{
			// Move the last half of the previous chunk data to the start of the data buffer
			size_t lastSampleCopySize = ((maxChunkSizeInBytes / sizeof(SampleType)) / 2);
			std::vector<SampleType> tempFileData(lastSampleCopySize + sampleCountInChunk);
			std::copy(fileData.end() - lastSampleCopySize, fileData.end(), tempFileData.begin());
			sampleScanningStartPos = lastSampleCopySize - (fileData.size() - continuePos);
			fileReadStartPos = lastSampleCopySize;
			fileData = std::move(tempFileData);
		}

		// Read in the next chunk of file data and store it in the buffer
		_logger.Info("Reading {0} bytes from file pos {1}", chunkSizeInBytes, currentFilePos);
		infile.read((char*)&fileData[fileReadStartPos], chunkSizeInBytes);
		currentFilePos += chunkSizeInBytes;

		// Detect sync pulses in the sample data
		_logger.Info("Detecting sync pulses");
		std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = _syncDetector.DetectSyncPulses(fileData, sampleScanningStartPos);

		// Build sync events from the raw sync pulses
		_logger.Info("Extracting sync events");
		std::list<SyncDetector::SyncInfo> syncInfo = _syncDetector.DetectSyncEvents(fileData, syncPulseInfo);

		// Collect the sync events into frames
		_logger.Info("Detecting frames");
		std::list<FrameBuilder::FieldInfo> fields = _frameBuilder.DetectFields(fileData, syncInfo);
		std::vector<FrameBuilder::FrameInfo> frames = _frameBuilder.DetectFrames(fields);

		// Decode line information within each detected frame
		_logger.Info("Detecting lines");
		_frameBuilder.DetectLines(fileData, frames);

		// Write each detected frame out to an image file
		_logger.Info("Writing frames");
		WriteFramesToBMP(outputFolderPath, outputFileNameBase, fileData, frames, startFrameNo);
		startFrameNo += frames.size();

		// Calculate a position to resume the search from when the next chunk of file data is read
		if (!frames.empty())
		{
			const FrameBuilder::FieldInfo& lastFieldInfo = frames.back().fields.back();
			const SyncDetector::SyncInfo& targetSyncInfo = lastFieldInfo.syncEvents[lastFieldInfo.syncEvents.size() / 2];
			continuePos = targetSyncInfo.startSampleNo;
		}
		else
		{
			continuePos = 0;
		}
	}

	// Return the result the caller
	return true;
}

//----------------------------------------------------------------------------------------
template<class SampleType>
void VideoDecoder::WriteFramesToBMP(const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::vector<FrameBuilder::FrameInfo>& frames, size_t initialFrameNo, unsigned int threadCount) const
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

	size_t chunkCount = threadCount;
	size_t framesPerChunk = frameCount / chunkCount;
	size_t extraFrames = frameCount - (framesPerChunk * chunkCount);
	std::vector<std::thread> workerThreads;
	for (unsigned int i = 0; i < chunkCount; ++i)
	{
		workerThreads.emplace_back(std::thread([&, i]
		{
			size_t frameNo = (i * framesPerChunk) + (i > 0 ? extraFrames : 0);
			size_t lastFrameNoForChunk = frameNo + (framesPerChunk + (i == 0 ? extraFrames : 0));
			while (frameNo < lastFrameNoForChunk)
			{
				const FrameBuilder::FrameInfo& frameEntry = frames[frameNo];
				PathString outputFilePath = outputFolderPath + PathSeparatorChar + outputFileNameBase + ToPathString(std::to_string(initialFrameNo + frameNo) + ".bmp");
				WriteFrameToBMP(outputFilePath, sampleData, frameEntry);
				++frameNo;
			}
		}));
	}
	for (auto& entry : workerThreads)
	{
		entry.join();
	}
}

//----------------------------------------------------------------------------------------
template<class SampleType>
bool VideoDecoder::WriteFrameToBMP(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo) const
{
	// Determine the total number of lines in this frame
	unsigned int lineCount = 0;
	for (const auto& fieldInfo : frameInfo.fields)
	{
		lineCount += (unsigned int)fieldInfo.lines.size();
	}

	// Calculate the amount of padding on each line. Lines are padded out to DWORD boundaries.
	unsigned int bitsPerPixel = 24;
	unsigned int lineByteCount = lineWidthInPixels * 3;
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
	bitmapHeader.biWidth = (int)lineWidthInPixels;
	bitmapHeader.biHeight = -((int)lineCount);
	bitmapHeader.biPlanes = 1;
	bitmapHeader.biBitCount = (short)bitsPerPixel;
	bitmapHeader.biCompression = BitmapCompressionType::RGB;
	bitmapHeader.biSizeImage = (int)fileSize;
	bitmapHeader.biXPelsPerMeter = 0;
	bitmapHeader.biYPelsPerMeter = 0;
	bitmapHeader.biClrUsed = 0;
	bitmapHeader.biClrImportant = 0;

	// Attempt to create an output file for this frame
	std::ofstream outfile(outputFilePath, std::ios_base::binary);
	if (outfile.fail())
	{
		_logger.Error("Failed to create output file \"{0}\"", outputFilePath);
		return false;
	}

	// Write the bitmap file headers to the output file
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
	std::vector<unsigned char> outputData(lineWidthInPixels*3);
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
	return true;
}
