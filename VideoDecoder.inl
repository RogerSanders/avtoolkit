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
	_log.Info("Decoding \"{0}\"", inputFilePath);

	// Attempt to open the target file, and calculate the size of the file in bytes.
	std::ifstream infile(inputFilePath, std::ios_base::binary);
	if (infile.fail())
	{
		_log.Error("Failed to open input file \"{0}\"", inputFilePath);
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
		_log.Info("Reading {0} bytes from file pos {1}", chunkSizeInBytes, currentFilePos);
		infile.read((char*)&fileData[fileReadStartPos], chunkSizeInBytes);
		currentFilePos += chunkSizeInBytes;

		// Detect sync pulses in the sample data
		_log.Info("Detecting sync pulses");
		std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = _syncDetector.DetectSyncPulses(fileData, sampleScanningStartPos);

		// Build sync events from the raw sync pulses
		_log.Info("Extracting sync events");
		std::list<SyncDetector::SyncInfo> syncInfo = _syncDetector.DetectSyncEvents(fileData, syncPulseInfo);

		// Collect the sync events into frames
		_log.Info("Detecting frames");
		std::list<FrameBuilder::FieldInfo> fields = _frameBuilder.DetectFields(fileData, syncInfo);
		std::vector<FrameBuilder::FrameInfo> frames = _frameBuilder.DetectFrames(fields);

		// Decode line information within each detected frame
		_log.Info("Detecting lines");
		_frameBuilder.DetectLines(fileData, frames);

		// Write each detected frame out to an image file
		_log.Info("Writing frames");
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
//##TODO## Separate the decoding of data from the writing of it to an image file
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

	// Calculate the sizes of our header structures. We need to do this rather than rely on the structure size in
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
		_log.Error("Failed to create output file \"{0}\"", outputFilePath);
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

	// Pre-calculate conversion factors we'll use during line decoding
	typedef unsigned char OutputSampleType;
	double normalizedSampleToOutputSampleConversionFactor = (double)std::numeric_limits<OutputSampleType>::max() - (double)std::numeric_limits<OutputSampleType>::min();
	double inputSampleToOutputSampleConversionFactor = normalizedSampleToOutputSampleConversionFactor / ((double)std::numeric_limits<SampleType>::max() - (double)std::numeric_limits<SampleType>::min());
	double inputSampleToNormalizedSampleConversionFactor = 1.0 / ((double)std::numeric_limits<SampleType>::max() - (double)std::numeric_limits<SampleType>::min());

	// Allocate scratch buffers we'll use during line decoding
	DecodeLinePixelDataBuffer scratchBuffer;
	std::vector<SampleType> leadingData;
	std::vector<SampleType> activeScanData;
	std::vector<double> activeScanDataR;
	std::vector<double> activeScanDataG;
	std::vector<double> activeScanDataB;
	std::vector<SampleType> followingData;
	std::vector<OutputSampleType> outputData(lineWidthInPixels*3);

	// Decode each field in this frame, and output it to the target image file.
	for (const FrameBuilder::FieldInfo& fieldInfo : frameInfo.fields)
	{
		// Calculate the average sync and blanking levels across the field if requested
		double ireLevel0ForField = 0.0;
		double ireLevel100ForField = 0.0;
		if (useAverageFieldBlankingLevel)
		{
			CalculateIRELevelsForField(fieldInfo, ireLevel0ForField, ireLevel100ForField);
		}

		// Calculate the average burst wave frequency across the field if requested
		double burstWaveFrequencyForField = 0.0;
		if (useAverageFieldBurstWaveFrequency)
		{
			burstWaveFrequencyForField = CalculateColorBurstWaveFrequencyForField(fieldInfo);
		}

		// Decode each line in this field, and output it to the target image file.
		for (size_t lineNo = 0; lineNo < fieldInfo.lines.size(); ++lineNo)
		{
			// Calculate the exact position of this line in the sample data
			const FrameBuilder::LineInfo& lineInfo = fieldInfo.lines[lineNo];
			double preciseLineStartPos = lineInfo.backPorchStartPos;
			double preciseLineEndPos = lineInfo.frontPorchEndPos;

			// Calculate the number of pixels in this line to output in the active scan and border regions
			unsigned int blankingLeadingPixelCount = (unsigned int)(blankingPercentage * blankingLeadingPercentage * (double)lineWidthInPixels);
			unsigned int blankingFollowingPixelCount = (unsigned int)(blankingPercentage * (double)lineWidthInPixels) - blankingLeadingPixelCount;
			unsigned int activeImagePixelCount = lineWidthInPixels - (blankingLeadingPixelCount + blankingFollowingPixelCount);

			// Build pixel data for the left border region of the line, up to the start of the active scan region.
			leadingData.resize(blankingLeadingPixelCount);
			CubicInterpolateCatmullRom(sampleData.data(), (double)lineInfo.leadingSyncInfo.startSampleNo, preciseLineStartPos, leadingData);

			// Build pixel data for the active scan region of the line
			activeScanData.resize(activeImagePixelCount);
			activeScanDataR.resize(activeImagePixelCount);
			activeScanDataG.resize(activeImagePixelCount);
			activeScanDataB.resize(activeImagePixelCount);
			if (!rawOutputOnly && (lineInfo.leadingSyncInfo.type == SyncDetector::SyncType::Horizontal))
			{
				// Decode the active scan region as visible pixels
				bool colorDecodingPerformed = false;
				if (decodeColor)
				{
					colorDecodingPerformed = DecodeColorLinePixelDataForActiveScanRegion(sampleData, fieldInfo, lineInfo, lineNo, preciseLineStartPos, preciseLineEndPos, burstWaveFrequencyForField, ireLevel0ForField, ireLevel100ForField, activeImagePixelCount, activeScanDataR, activeScanDataG, activeScanDataB, scratchBuffer);
				}
				if (!colorDecodingPerformed)
				{
					DecodeMonochromeLinePixelDataForActiveScanRegion(sampleData, fieldInfo, lineInfo, lineNo, preciseLineStartPos, preciseLineEndPos, ireLevel0ForField, ireLevel100ForField, activeImagePixelCount, activeScanDataR, activeScanDataG, activeScanDataB, scratchBuffer);
				}
			}
			else
			{
				// Rescale the raw sample data for the active scan region of the line to match the output sample count
				CubicInterpolateCatmullRom(sampleData.data(), preciseLineStartPos, preciseLineEndPos, activeScanData);

				// Convert the sample data into a normalized form, and load it into the output colour channel buffers.
				for (size_t i = 0; i < activeImagePixelCount; ++i)
				{
					double sampleValue = ((double)activeScanData[i] + -(double)std::numeric_limits<SampleType>::min()) * inputSampleToNormalizedSampleConversionFactor;
					activeScanDataR[i] = sampleValue;
					activeScanDataG[i] = sampleValue;
					activeScanDataB[i] = sampleValue;
				}
			}

			// Build pixel data for the right border region of the line, from the end of active scan up to the start of
			// the next line.
			followingData.resize(blankingFollowingPixelCount);
			CubicInterpolateCatmullRom(sampleData.data(), preciseLineEndPos, lineInfo.followingSyncInfo.startSampleNo + ((lineInfo.followingSyncInfo.endSampleNo - lineInfo.followingSyncInfo.startSampleNo) * (1.0 - blankingLeadingPercentage)), followingData);

			// Output the left border region of the line
			size_t outputDataPos = 0;
			for (size_t i = 0; i < blankingLeadingPixelCount; ++i)
			{
				SampleType sampleValue = leadingData[i];
				OutputSampleType outputSampleValue = (OutputSampleType)((double)std::numeric_limits<OutputSampleType>::min() + (((double)sampleValue + -(double)std::numeric_limits<SampleType>::min()) * inputSampleToOutputSampleConversionFactor));
				outputData[outputDataPos++] = outputSampleValue;
				outputData[outputDataPos++] = outputSampleValue;
				outputData[outputDataPos++] = outputSampleValue;
			}

			// Output the active scan region of the line
			for (size_t i = 0; i < activeImagePixelCount; ++i)
			{
				OutputSampleType sampleValueAsByteR = (OutputSampleType)((double)std::numeric_limits<OutputSampleType>::min() + (activeScanDataR[i] * normalizedSampleToOutputSampleConversionFactor));
				OutputSampleType sampleValueAsByteG = (OutputSampleType)((double)std::numeric_limits<OutputSampleType>::min() + (activeScanDataG[i] * normalizedSampleToOutputSampleConversionFactor));
				OutputSampleType sampleValueAsByteB = (OutputSampleType)((double)std::numeric_limits<OutputSampleType>::min() + (activeScanDataB[i] * normalizedSampleToOutputSampleConversionFactor));
				outputData[outputDataPos++] = sampleValueAsByteB;
				outputData[outputDataPos++] = sampleValueAsByteG;
				outputData[outputDataPos++] = sampleValueAsByteR;
			}

			// Output the right border region of the line
			for (size_t i = 0; i < blankingFollowingPixelCount; ++i)
			{
				SampleType sampleValue = followingData[i];
				OutputSampleType outputSampleValue = (OutputSampleType)((double)std::numeric_limits<OutputSampleType>::min() + (((double)sampleValue + -(double)std::numeric_limits<SampleType>::min()) * inputSampleToOutputSampleConversionFactor));
				outputData[outputDataPos++] = outputSampleValue;
				outputData[outputDataPos++] = outputSampleValue;
				outputData[outputDataPos++] = outputSampleValue;
			}

			// Write the pixel data for this line to the output file
			outfile.write((char*)&outputData[0], outputData.size());
			for (unsigned int i = 0; i < linePaddingByteCount; ++i)
			{
				outfile.write("", 1);
			}
		}
	}
	return true;
}

//----------------------------------------------------------------------------------------
template<class SampleType>
void VideoDecoder::DecodeMonochromeLinePixelDataForActiveScanRegion(const std::vector<SampleType>& sampleData, const FrameBuilder::FieldInfo& fieldInfo, const FrameBuilder::LineInfo& lineInfo, size_t lineNo, double preciseLineStartPos, double preciseLineEndPos, double ireLevel0ForField, double ireLevel100ForField, unsigned int outputPixelCount, std::vector<double>& outputDataR, std::vector<double>& outputDataG, std::vector<double>& outputDataB, DecodeLinePixelDataBuffer& buffer) const
{
	//##TODO##
	// If we're decoding a monochrome image, we pass the decoded "Y" channel (luma) back as the content for all
	// three colour channels, and advance to the next sample.
}

//----------------------------------------------------------------------------------------
template<class SampleType>
bool VideoDecoder::DecodeColorLinePixelDataForActiveScanRegion(const std::vector<SampleType>& sampleData, const FrameBuilder::FieldInfo& fieldInfo, const FrameBuilder::LineInfo& lineInfo, size_t lineNo, double preciseLineStartPos, double preciseLineEndPos, double burstWaveFrequencyForField, double ireLevel0ForField, double ireLevel100ForField, unsigned int outputPixelCount, std::vector<double>& outputDataR, std::vector<double>& outputDataG, std::vector<double>& outputDataB, DecodeLinePixelDataBuffer& buffer) const
{
	// Retrieve our temp storage buffers from the supplied buffer object. We use this approach to allow us to reuse
	// scratch buffers (which have already reserved appropriate amounts of memory) for successive lines. This avoids
	// expensive heap allocation and deallocation for our larger buffers.
	std::vector<double>& carrierWaveISamplePoints = buffer.carrierWaveISamplePoints;
	std::vector<double>& carrierWaveQSamplePoints = buffer.carrierWaveQSamplePoints;
	std::vector<double>& samplePointsQ = buffer.samplePointsQ;
	std::vector<double>& samplePointsI = buffer.samplePointsI;
	std::vector<double>& samplePointsY = buffer.samplePointsY;
	std::vector<double>& lineResolutionOutputR = buffer.lineResolutionOutputR;
	std::vector<double>& lineResolutionOutputG = buffer.lineResolutionOutputG;
	std::vector<double>& lineResolutionOutputB = buffer.lineResolutionOutputB;

	// Calculate the positions of each peak of the colour burst wave for this line
	std::vector<double> wavePeakPositions;
	bool firstWavePeakIsPositive;
	if (!GetWavePeakPositionsForLine(lineInfo, wavePeakPositions, firstWavePeakIsPositive))
	{
		return false;
	}

	// Determine the half-frequency in input samples of the colour burst wave to use for this line
	double burstWaveFrequency;
	if (forceColorBurstWaveFrequency)
	{
		burstWaveFrequency = forcedColorBurstWaveFrequency;
	}
	else if (useAverageFieldBurstWaveFrequency)
	{
		burstWaveFrequency = burstWaveFrequencyForField;
	}
	else
	{
		burstWaveFrequency = CalculateColorBurstWaveFrequencyForLine(lineInfo, wavePeakPositions);
	}

	// Calculate the region of the sample data to extract for the active scan region of this line, and load it into a
	// cubic spline for sampling.
	const unsigned int activeScanSampleWaveOverrun = 4;
	size_t activeScanFirstSampleNo = (size_t)preciseLineStartPos - ((size_t)burstWaveFrequency * 2 * activeScanSampleWaveOverrun);
	size_t activeScanLastSampleNo = (size_t)preciseLineEndPos + ((size_t)burstWaveFrequency * 2 * activeScanSampleWaveOverrun);
	CubicSpline<SampleType> activeScanAsSpline(sampleData.data() + activeScanFirstSampleNo, activeScanLastSampleNo - activeScanFirstSampleNo);

	// Select an initial search position to use to phase lock with the colour burst
	bool firstQWaveIsPositive;
	double carrierWaveSearchStartPos;
	if (!matchColorBurstPhaseBetweenLines || ((lineNo + 1) >= fieldInfo.lines.size()))
	{
		// Select the middle colour burst wave peak as the position we phase lock the colour burst to. We use the middle
		// wave entry as entries towards the ends are more prone to being incorrect.
		firstQWaveIsPositive = firstWavePeakIsPositive;
		carrierWaveSearchStartPos = wavePeakPositions[(wavePeakPositions.size() / 2) & ~1];
	}
	else
	{
		// Calculate the positions of each peak of the colour burst wave in the following line
		const FrameBuilder::LineInfo& nextLineInfo = fieldInfo.lines[lineNo+1];
		std::vector<double> nextLineWavePeakPositions;
		bool nextLineFirstWavePeakIsPositive;
		GetWavePeakPositionsForLine(nextLineInfo, nextLineWavePeakPositions, nextLineFirstWavePeakIsPositive);

		// Evaluate our phase lock with the colour burst on the following line, and determine the wave peak position
		// from the current line that gives the best phase lock with the calculated colour burst wave frequency.
		size_t bestStartingPeakPosition;
		bool bestWavePeakPositionIsPositive;
		double bestDisplacementValue;
		PhaseLockColorBurstSamples(wavePeakPositions, firstWavePeakIsPositive, nextLineWavePeakPositions, nextLineFirstWavePeakIsPositive, burstWaveFrequency, bestStartingPeakPosition, bestWavePeakPositionIsPositive, bestDisplacementValue);

		// Select the burst wave peak with the best phase lock to the following line as the position we phase lock the
		// colour burst to.
		firstQWaveIsPositive = bestWavePeakPositionIsPositive;
		carrierWaveSearchStartPos = wavePeakPositions[bestStartingPeakPosition];
	}

	// As per "Video Demystified", 3rd edition, chapter 8, page 240, YUV and YIQ are equivalent colour spaces, with Q
	// corresponding to U (B-Y), and I corresponding to V(R-Y), except the axes are rotated 33 degrees counterclockwise
	// under YIQ compared with YUV. The "cross-talk" between I and Q which occurs when you sample with this 33 degree
	// phase shift creates a hue variation which is equal to the difference between the colour spaces. Although the
	// original video signal may have been encoded colour either in YIQ at a 33 degree offset from the colour burst, or
	// in YUV synchronized with the colour burst, the relationship between these colour planes means it shouldn't
	// actually matter which encoding scheme was used. We can decode in either YUV or YIQ, and get an equivalent
	// result. In practice, slight hue variation is evident between the two colour planes, so where the true original
	// plane which was used to encode the colour signal is known, it would be preferable to use that. It isn't possible
	// to determine this information from the signal itself however, and it's unlikely to be a significant factor, as
	// the variation in hue between these planes is small in comparison to the natural variation in hue which would
	// occur between real analog displays in the consumer market. Signals are also likely to have been mastered in many
	// cases to take hue variations into account which occur on consumer hardware, which typically decodes as YUV. It's
	// also simpler to decode YUV, since we don't have to rotate 33 degrees back. Also note that PAL has only ever
	// encoded as YUV. For these reasons, YUV decoding is a better default. If we're decoding in true YIQ mode here
	// however, we perform a 33 degree phase shift from the reference colour burst, to get the correct sampling location
	// for the YIQ data. Note that we always refer to our sample buffers and variables using the terms "I" and "Q"
	// below, but where YUV decoding is being performed, Q corresponds to U, and I corresponds to V.
	//##TODO## Consider renaming the variables from "Q" and "I" to "U" and "V" respectively, since YUV is the preferred
	//colour space for decoding.
	if (!decodeAsYUV)
	{
		carrierWaveSearchStartPos -= (33.0 / 90.0);
	}

	// Calculate the reference start positions for the I and Q carrier waves
	while (carrierWaveSearchStartPos >= (preciseLineStartPos - burstWaveFrequency))
	{
		carrierWaveSearchStartPos -= burstWaveFrequency * 2;
	}
	double carrierWaveQStartPos = carrierWaveSearchStartPos;
	double carrierWaveIStartPos = carrierWaveSearchStartPos - (burstWaveFrequency / 2);
	size_t positiveOscillationOffset = (firstQWaveIsPositive ? 1 : 0);

	// Build up a set of sample points for the I carrier wave at each peak position only. This should consist of only
	// the I and Y signals.
	size_t carrierWaveReserveSize = (size_t)((double)(activeScanLastSampleNo - activeScanFirstSampleNo) / burstWaveFrequency) + 10;
	carrierWaveISamplePoints.clear();
	carrierWaveISamplePoints.reserve(carrierWaveReserveSize);
	double carrierWaveISearchPos = carrierWaveIStartPos;
	while (carrierWaveISearchPos < ((double)activeScanLastSampleNo + (burstWaveFrequency / 2)))
	{
		double samplePoint = activeScanAsSpline.Evaluate(carrierWaveISearchPos - (double)activeScanFirstSampleNo);
		carrierWaveISamplePoints.push_back(samplePoint);
		carrierWaveISearchPos += burstWaveFrequency;
	}
	double carrierWaveIEndPos = carrierWaveISearchPos;

	// Build up a set of sample points for the Q carrier wave at each peak position only. This should consist of only
	// the Q and Y signals.
	carrierWaveQSamplePoints.clear();
	carrierWaveQSamplePoints.reserve(carrierWaveReserveSize);
	double carrierWaveQSearchPos = carrierWaveQStartPos;
	while (carrierWaveQSearchPos < ((double)activeScanLastSampleNo + (burstWaveFrequency / 2)))
	{
		double samplePoint = activeScanAsSpline.Evaluate(carrierWaveQSearchPos - (double)activeScanFirstSampleNo);
		carrierWaveQSamplePoints.push_back(samplePoint);
		carrierWaveQSearchPos += burstWaveFrequency;
	}
	double carrierWaveQEndPos = carrierWaveQSearchPos;

	// Calculate sample points for the I signal, without Y.
	samplePointsI.resize(carrierWaveISamplePoints.size() - 1);
	for (size_t i = positiveOscillationOffset; i < (carrierWaveISamplePoints.size() - 1); i += 2)
	{
		samplePointsI[i] = (carrierWaveISamplePoints[i] - carrierWaveISamplePoints[i+1]) / 2.0;
	}
	for (size_t i = 1 - positiveOscillationOffset; i < (carrierWaveISamplePoints.size() - 1); i += 2)
	{
		samplePointsI[i] = (carrierWaveISamplePoints[i+1] - carrierWaveISamplePoints[i]) / 2.0;
	}

	// Calculate sample points for the Q signal, without Y.
	samplePointsQ.resize(carrierWaveQSamplePoints.size() - 1);
	for (size_t i = positiveOscillationOffset; i < (carrierWaveQSamplePoints.size() - 1); i += 2)
	{
		samplePointsQ[i] = (carrierWaveQSamplePoints[i] - carrierWaveQSamplePoints[i+1]) / 2.0;
	}
	for (size_t i = 1 - positiveOscillationOffset; i < (carrierWaveQSamplePoints.size() - 1); i += 2)
	{
		samplePointsQ[i] = (carrierWaveQSamplePoints[i+1] - carrierWaveQSamplePoints[i]) / 2.0;
	}

	// Calculate sample points for the Y signal, without I or Q.
	double samplePointsYStartPos = (carrierWaveQStartPos < carrierWaveIStartPos ? carrierWaveQStartPos : carrierWaveIStartPos);
	double samplePointsYEndPos = (carrierWaveQEndPos > carrierWaveIEndPos ? carrierWaveQEndPos : carrierWaveIEndPos);
	size_t samplePointsYCount = samplePointsQ.size() + samplePointsI.size();
	samplePointsY.resize(samplePointsYCount);
	for (size_t i = 0, sourcePos = 0; i < samplePointsYCount; i += 2, ++sourcePos)
	{
		if ((sourcePos % 2) == positiveOscillationOffset)
		{
			samplePointsY[i] = carrierWaveISamplePoints[sourcePos] - samplePointsI[sourcePos];
		}
		else
		{
			samplePointsY[i] = carrierWaveISamplePoints[sourcePos] + samplePointsI[sourcePos];
		}
	}
	for (size_t i = 1, sourcePos = 0; i < samplePointsYCount; i += 2, ++sourcePos)
	{
		if ((sourcePos % 2) == positiveOscillationOffset)
		{
			samplePointsY[i] = carrierWaveQSamplePoints[sourcePos] - samplePointsQ[sourcePos];
		}
		else
		{
			samplePointsY[i] = carrierWaveQSamplePoints[sourcePos] + samplePointsQ[sourcePos];
		}
	}

	// Calculate the IRE levels to use when converting the active scan signals
	//##FIX## These levels are based on NTSC signal specs. Allow for PAL specification here.
	double ireLevel0 = (useAverageFieldBlankingLevel ? ireLevel0ForField : lineInfo.ireLevel0);
	double ireLevel100 = (useAverageFieldBlankingLevel ? ireLevel100ForField : lineInfo.ireLevel100);
	double ireLevel7Point5 = ireLevel0 + ((ireLevel100ForField - lineInfo.ireLevel0) * 0.075);
	double baseIRELevel = (useIRE7Point5 ? ireLevel7Point5 : ireLevel0);
	const double ireLevelForMaxIQSignal = 20.0;
	double maxLevelForIQSignals = (double)(_frameBuilder.IREToSample<SampleType>((float)ireLevelForMaxIQSignal, (float)ireLevel0, (float)ireLevel100) - (SampleType)ireLevel0);
	double sampleIQScaleFactor = 1.0 / (maxLevelForIQSignals * 2);

	// Output pixel data to the colour buffers for the active scan region of this line
	lineResolutionOutputR.resize(samplePointsYCount);
	lineResolutionOutputG.resize(samplePointsYCount);
	lineResolutionOutputB.resize(samplePointsYCount);
	for (size_t samplePointYNo = 0; samplePointYNo < samplePointsYCount; ++samplePointYNo)
	{
		// Calculate the "Y" sample data (luma)
		double sampleY = samplePointsY[samplePointYNo];
		sampleY = (std::min(std::max(baseIRELevel, sampleY), ireLevel100) - baseIRELevel) / (ireLevel100 - baseIRELevel);

		// Calculate the "I" sample data at this location, either as the direct sample corresponding with this point, or
		// as the average of the surrounding two samples.
		size_t samplePointIQNo = samplePointYNo / 2;
		bool currentSampleIsI = (samplePointYNo - samplePointIQNo) == 0;
		double sampleI = ((currentSampleIsI || ((samplePointIQNo + 1) > samplePointsI.size())) ? samplePointsI[samplePointIQNo] : ((samplePointsI[samplePointIQNo] + samplePointsI[samplePointIQNo + 1]) / 2.0));
		sampleI *= sampleIQScaleFactor;
		sampleI = std::min(1.0, std::max(-1.0, sampleI));

		// Calculate the "Q" sample data at this location, either as the direct sample corresponding with this point, or
		// as the average of the surrounding two samples.
		double sampleQ = ((!currentSampleIsI || ((samplePointIQNo + 1) > samplePointsQ.size())) ? samplePointsQ[samplePointIQNo] : ((samplePointsQ[samplePointIQNo] + samplePointsQ[samplePointIQNo + 1]) / 2.0));
		sampleQ *= sampleIQScaleFactor;
		sampleQ = std::min(1.0, std::max(-1.0, sampleQ));

		// Decode the RGB colour for this pixel
		double redNormalized;
		double greenNormalized;
		double blueNormalized;
		if (forceMonoOutputForColorDecoding)
		{
			redNormalized = sampleY;
			greenNormalized = sampleY;
			blueNormalized = sampleY;
		}
		else if (decodeAsYUV)
		{
			ConvertYUVToRGB(sampleY, sampleQ, sampleI, redNormalized, greenNormalized, blueNormalized);
		}
		else
		{
			ConvertYIQToRGB(sampleY, sampleI, sampleQ, redNormalized, greenNormalized, blueNormalized);
		}

		// Clamp the output colour components to the range 0.0 - 1.0
		lineResolutionOutputR[samplePointYNo] = std::min(std::max(0.0, redNormalized), 1.0);
		lineResolutionOutputG[samplePointYNo] = std::min(std::max(0.0, greenNormalized), 1.0);
		lineResolutionOutputB[samplePointYNo] = std::min(std::max(0.0, blueNormalized), 1.0);
	}

	// Rescale the output pixel data to match the required horizontal output resolution, and clip the pixel region to
	// the precise active scan boundaries, since we will have calculated sample points slightly past the beginning and
	// end of active scan.
	double sampleDataToSamplePointsYScaleFactor = (double)samplePointsYCount / (samplePointsYEndPos - samplePointsYStartPos);
	double activeScanStartPosInBuffer = (preciseLineStartPos - samplePointsYStartPos) * sampleDataToSamplePointsYScaleFactor;
	double activeScanEndPosInBuffer = (double)samplePointsYCount - ((samplePointsYEndPos - preciseLineEndPos) * sampleDataToSamplePointsYScaleFactor);
	CubicInterpolateCatmullRom(lineResolutionOutputR.data(), activeScanStartPosInBuffer, activeScanEndPosInBuffer, outputDataR);
	CubicInterpolateCatmullRom(lineResolutionOutputG.data(), activeScanStartPosInBuffer, activeScanEndPosInBuffer, outputDataG);
	CubicInterpolateCatmullRom(lineResolutionOutputB.data(), activeScanStartPosInBuffer, activeScanEndPosInBuffer, outputDataB);

	// Since we've just performed interpolation, we need to clamp our output colour components again to ensure they're
	// still in the range 0.0 - 1.0.
	for (size_t i = 0; i < outputPixelCount; ++i)
	{
		outputDataR[i] = std::min(std::max(0.0, outputDataR[i]), 1.0);
		outputDataG[i] = std::min(std::max(0.0, outputDataG[i]), 1.0);
		outputDataB[i] = std::min(std::max(0.0, outputDataB[i]), 1.0);
	}
	return true;
}
