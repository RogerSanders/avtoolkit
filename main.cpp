#define NOMINMAX
#include <windows.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "SyncDetector.h"
#include "FrameBuilder.h"
#include "spline.h"

//----------------------------------------------------------------------------------------
template<class T>
void ResampleLinear(const std::vector<T>& inputData, std::vector<T>& outputData)
{
	//This method is much more sophisticated than the one above, and will give the best
	//result possible from linear resampling when either upsampling or downsampling,
	//regardless of the respective dimensions of the source and target images.
	size_t oldSize = inputData.size();
	size_t newSize = outputData.size();
	float imageWidthConversionRatio = (float)oldSize / (float)newSize;

	//Generate each pixel in this line
	for(unsigned int xpos = 0; xpos < newSize; ++xpos)
	{
		//Calculate the beginning and end of the sample region on the X axis in the
		//source image, which is being mapped onto this pixel in the target image.
		//Note that because we're adding 1 to the current pixel location, the second
		//sample point may be past the end of the image, but this is ok, because this
		//sample point is really a limit. We protect against reading past the
		//boundaries of the image further below, and if this is attempted, this dud
		//sample will be assigned a weight of zero.
		float firstSamplePointX = (float)xpos * imageWidthConversionRatio;
		float lastSamplePointX = (float)(xpos + 1) * imageWidthConversionRatio;

		//Calculate the total domain, or length, of this sample region on the X axis.
		float totalDomainX = lastSamplePointX - firstSamplePointX;

		//Calculate the first and last pixels of interest from the source region
		unsigned int firstSamplePosX = (unsigned int)firstSamplePointX;
		unsigned int lastSamplePosX = (unsigned int)lastSamplePointX;

		//Calculate the total domain, or area, of the sample region in the source
		//image. We use this to normalize the sampled data back to an area of one
		//pixel at the end.
		float totalDomain = totalDomainX;

		//Combine sample values from the source image, on both the X and Y axis,
		//with their respective weightings.
		float finalSample = 0.0f;
		for(unsigned int currentSampleX = firstSamplePosX; currentSampleX <= lastSamplePosX; ++currentSampleX)
		{
			float sampleStartPointX = 0.0f;
			if(currentSampleX == firstSamplePosX)
			{
				sampleStartPointX = firstSamplePointX - (float)firstSamplePosX;
			}
			float sampleEndPointX = 1.0f;
			if(currentSampleX == lastSamplePosX)
			{
				sampleEndPointX = lastSamplePointX - (float)lastSamplePosX;
			}
			float sampleWeightX = sampleEndPointX - sampleStartPointX;

			float sample = (float)inputData[currentSampleX] / (float)std::numeric_limits<T>::max();
			finalSample += sample * sampleWeightX;
		}

		//Normalize the sample value back to a single pixel value, by dividing it
		//by the total area of the sample region in the source image.
		finalSample /= totalDomain;

		//Write the generated pixel to the image
		outputData[xpos] = (T)((finalSample * (float)std::numeric_limits<T>::max()) + (std::numeric_limits<T>::is_integer ? 0.5f : 0.0f));
	}
}

//----------------------------------------------------------------------------------------
template<class T>
void ResampleCubic(const std::vector<T>& inputData, std::vector<T>& outputData)
{
	std::vector<double> x(inputData.size());
	std::vector<double> y(inputData.size());
	for (unsigned int i = 0; i < inputData.size(); ++i)
	{
		x[i] = (double)i;
		y[i] = (double)inputData[i];
	}

	tk::spline s;
	s.set_points(x, y);

	for (unsigned int i = 0; i < outputData.size(); ++i)
	{
		double samplePos = ((double)i * ((double)inputData.size() / (double)outputData.size()));
		outputData[i] = (T)s(samplePos);
	}
}

//----------------------------------------------------------------------------------------
//##TODO## Put this in a class
//##TODO## This function does a very good job, but not as good as the colour burst phase offset detection approach used
//in ld-decode. We should consider that approach the primary method, and use our approach here as a secondary method.
//##TODO## Then again, it's questionable locking to colour burst is actually the right approach. A real monitor locks to
//the sync signal, and with a monochrome signal there might be no colour burst at all to lock to. Perform more
//comparisons which focus on the actual image region, and ignore the synchronization of the colour burst between lines,
//to determine if the colour burst appears to be a more or less reliable sync point.
void FindPreciseLineStartEndPos(const SyncDetector::SyncInfo& leadingSyncInfo, const SyncDetector::SyncInfo& followingSyncInfo, const tk::spline& lineSpline, double& preciseLineStartPos, double& preciseLineEndPos)
{
	//##FIX## Make constants like these configurable
	//##DEBUG## Pretty good for sonic 2, not as good when large min/max range
	//double syncSearchPosIncrement = 0.0001;
	//double syncAmplitudeMinTolerance = 0.2;
	//double slopeValueFlatToleranceAsPercentage = 0.3;

	double syncSearchPosIncrement = 0.0001;
	double syncAmplitudeMinTolerance = 0.125;
	double slopeValueFlatToleranceAsPercentage = 0.3;

	size_t lineSplineSampleRange = followingSyncInfo.endSampleNo - leadingSyncInfo.startSampleNo;

	// Find the point at which we cross the minimum threshold for ending the leading sync run
	double leadingSyncEndPosInSpline = (double)(leadingSyncInfo.endSampleNo - leadingSyncInfo.startSampleNo) / (double)lineSplineSampleRange;
	double leadingSyncEndMinimumThreshold = leadingSyncInfo.averageSyncLevel + (leadingSyncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance);
	double leadingSyncEndSearchPosInSpline = leadingSyncEndPosInSpline;
	while (lineSpline(leadingSyncEndSearchPosInSpline) < leadingSyncEndMinimumThreshold)
	{
		leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
	}

	// Find the point at which the leading sync run levels out to an acceptably flat slope, and mark it as the precise
	// line start pos.
	double leadingSyncLastSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
	double leadingSyncCurrentSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	double leadingSyncSlopeValueFlatTolerance = leadingSyncInfo.approxMinMaxSampleRange * slopeValueFlatToleranceAsPercentage;
	while (std::fabs(leadingSyncCurrentSlopeValue - leadingSyncLastSlopeValue) > leadingSyncSlopeValueFlatTolerance)
	{
		leadingSyncEndSearchPosInSpline += syncSearchPosIncrement;
		leadingSyncLastSlopeValue = leadingSyncCurrentSlopeValue;
		leadingSyncCurrentSlopeValue = lineSpline(leadingSyncEndSearchPosInSpline);
	}
	preciseLineStartPos = leadingSyncEndSearchPosInSpline;

	// Find the point at which we cross the minimum threshold for starting the following sync run
	double followingSyncStartPosInSpline = (double)(followingSyncInfo.startSampleNo - leadingSyncInfo.startSampleNo) / (double)lineSplineSampleRange;
	double followingSyncStartMinimumThreshold = followingSyncInfo.averageSyncLevel + (followingSyncInfo.approxMinMaxSampleRange * syncAmplitudeMinTolerance);
	double followingSyncStartSearchPosInSpline = followingSyncStartPosInSpline;
	while (lineSpline(followingSyncStartSearchPosInSpline) < followingSyncStartMinimumThreshold)
	{
		followingSyncStartSearchPosInSpline -= syncSearchPosIncrement;
	}

	// Find the point at which the following sync run levels out to an acceptably flat slope, and mark it as the precise
	// line end pos.
	double followingSyncLastSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	followingSyncStartSearchPosInSpline += syncSearchPosIncrement;
	double followingSyncCurrentSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	double followingSyncSlopeValueFlatTolerance = followingSyncInfo.approxMinMaxSampleRange * slopeValueFlatToleranceAsPercentage;
	while (std::fabs(followingSyncCurrentSlopeValue - followingSyncLastSlopeValue) > followingSyncSlopeValueFlatTolerance)
	{
		followingSyncStartSearchPosInSpline -= syncSearchPosIncrement;
		followingSyncLastSlopeValue = followingSyncCurrentSlopeValue;
		followingSyncCurrentSlopeValue = lineSpline(followingSyncStartSearchPosInSpline);
	}
	preciseLineEndPos = followingSyncStartSearchPosInSpline;
}

//----------------------------------------------------------------------------------------
template<class SampleType>
void WriteFrameToFile(const std::wstring& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo)
{
	const FrameBuilder::FieldInfo& fieldEntry = frameInfo.fieldInfo.front();
	unsigned int lineCount = fieldEntry.lineCount;
	//##FIX##
	//unsigned int pixelsPerLine = (unsigned int)((double)lineCount * (4.0 / 3.0) + 0.5);
	unsigned int pixelsPerLine = 910;
	double blankingPercentage = 0.1;
	double blankingLeadingPercentage = 0.95;

	unsigned int bitsPerPixel = 24;
	unsigned int bytesPerPixel = ((bitsPerPixel + 7) / 8);
	unsigned int lineByteCount = pixelsPerLine * 3;

	//Calculate the amount of padding on each line. Lines are padded out to DWORD
	//boundaries.
	unsigned int linePaddingByteCount = 0;
	if((lineByteCount % sizeof(DWORD)) != 0)
	{
		linePaddingByteCount = sizeof(DWORD) - (lineByteCount % sizeof(DWORD));
	}

	//Calculate the totals and offsets we need to write to the file
	unsigned int pixelDataOffset = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	unsigned int pixelDataSize = (unsigned int)((lineByteCount + linePaddingByteCount) * lineCount);
	unsigned int fileSize = pixelDataOffset + pixelDataSize;

	//Write the bitmap file header to the file
	BITMAPFILEHEADER fileHeader;
	char typeByte1 = 'B';
	char typeByte2 = 'M';
	fileHeader.bfType = ((WORD)typeByte2 << 8) | (WORD)typeByte1;
	fileHeader.bfSize = fileSize;
	fileHeader.bfReserved1 = 0;
	fileHeader.bfReserved2 = 0;
	fileHeader.bfOffBits = pixelDataOffset;
	//Write the bitmap info header to the file
	BITMAPINFOHEADER bitmapHeader;
	ZeroMemory(&bitmapHeader, sizeof(BITMAPINFOHEADER));
	bitmapHeader.biSize = sizeof(BITMAPINFOHEADER);
	bitmapHeader.biWidth = (LONG)pixelsPerLine;
	bitmapHeader.biHeight = -((LONG)lineCount);
	bitmapHeader.biPlanes = 1;
	bitmapHeader.biBitCount = (WORD)bitsPerPixel;
	bitmapHeader.biCompression = BI_RGB;
	bitmapHeader.biSizeImage = (DWORD)fileSize;
	bitmapHeader.biXPelsPerMeter = 0;
	bitmapHeader.biYPelsPerMeter = 0;
	bitmapHeader.biClrUsed = 0;
	bitmapHeader.biClrImportant = 0;

	//##DEBUG##
	//std::ofstream outfile(outputFilePath, std::ios_base::binary);
	//outfile.write((char*)&fileHeader, sizeof(fileHeader));
	//outfile.write((char*)&bitmapHeader, sizeof(bitmapHeader));

	//Write the pixel data to the file
	float sampleConversionFactor = ((float)std::numeric_limits<unsigned char>::max() / ((float)std::numeric_limits<SampleType>::max() + -(float)std::numeric_limits<SampleType>::min()));
	auto syncEventIterator = fieldEntry.syncEvents.cbegin();
	for (unsigned int lineNo = 0; lineNo < lineCount; ++lineNo)
	{
		//##TODO##
		while ((syncEventIterator != fieldEntry.syncEvents.cend()) && (syncEventIterator->type != SyncDetector::SyncType::Horizontal))
		{
			++syncEventIterator;
		}
		if (syncEventIterator == fieldEntry.syncEvents.cend())
		{
			std::cout << "Incorrect line count!\n";
			system("pause");
			break;
		}
		const SyncDetector::SyncInfo& syncEvent = *syncEventIterator;
		++syncEventIterator;
		const SyncDetector::SyncInfo& nextSyncEvent = (syncEventIterator != fieldEntry.syncEvents.cend()) ? *syncEventIterator : fieldEntry.followingSyncEvent;

		std::vector<SampleType> inputData(nextSyncEvent.endSampleNo - syncEvent.startSampleNo);
		for (unsigned int i = 0; i < inputData.size(); ++i)
		{
			inputData[i] = sampleData[syncEvent.startSampleNo + i];
		}

		//##DEBUG##
		//std::vector<SampleType> outputDataLinear(pixelsPerLine);
		//std::vector<SampleType> outputData(pixelsPerLine);
		//ResampleLinear(inputData, outputDataLinear);
		//ResampleCubic(inputData, outputData);

		//##DEBUG##
		std::vector<double> x(inputData.size());
		std::vector<double> y(inputData.size());
		for (unsigned int i = 0; i < inputData.size(); ++i)
		{
			x[i] = (double)i / (double)(inputData.size() - 1);
			y[i] = (double)inputData[i];
		}

		tk::spline s;
		s.set_points(x, y);

		double preciseLineStartPos;
		double preciseLineEndPos;
		FindPreciseLineStartEndPos(syncEvent, nextSyncEvent, s, preciseLineStartPos, preciseLineEndPos);
		//##DEBUG##
		double preciseLineWidth = (preciseLineEndPos - preciseLineStartPos);
		//std::wcout << "LineDetect: " << preciseLineStartPos << "\t" << preciseLineEndPos << "\t" << preciseLineWidth << '\n';

		std::vector<SampleType> outputData(pixelsPerLine);
		unsigned int outputSampleCountForLine = (unsigned int)outputData.size();
		unsigned int blankingLeadingSampleCount = (unsigned int)(blankingPercentage * blankingLeadingPercentage * (double)outputSampleCountForLine);
		unsigned int blankingFollowingSampleCount = (unsigned int)(blankingPercentage * (double)outputSampleCountForLine) - blankingLeadingSampleCount;
		unsigned int activeImageSampleCount = outputSampleCountForLine - (blankingLeadingSampleCount + blankingFollowingSampleCount);
		for (unsigned int i = 0; i < outputSampleCountForLine; ++i)
		{
			//##DEBUG##
			//double samplePos = (double)i * (double)(outputData.size() - 1);
			//double samplePos = (preciseLineStartPos - (preciseLineWidth * blankingLeadingPercentage * blankingPercentage)) + ((double)i * ((preciseLineWidth + (preciseLineWidth * blankingPercentage)) / outputData.size()));
			double samplePos;
			if (i < blankingLeadingSampleCount)
			{
				samplePos = (double)i * (preciseLineStartPos / blankingLeadingSampleCount);
			}
			else if (i < (blankingLeadingSampleCount + activeImageSampleCount))
			{
				samplePos = preciseLineStartPos + ((double)(i - blankingLeadingSampleCount) * ((preciseLineEndPos - preciseLineStartPos) / activeImageSampleCount));
			}
			else
			{
				samplePos = preciseLineEndPos + ((double)(i - (blankingLeadingSampleCount + activeImageSampleCount)) * ((preciseLineStartPos / blankingLeadingPercentage) / blankingFollowingSampleCount));
			}
			outputData[i] = (SampleType)s(samplePos);
		}

		//##DEBUG##
		//for (unsigned int pixelNo = 0; pixelNo < pixelsPerLine; ++pixelNo)
		//{
		//	for (unsigned int i = 0; i < 3; ++i)
		//	{
		//		unsigned char temp = (unsigned char)(((float)outputData[pixelNo] + -(float)std::numeric_limits<SampleType>::min()) * sampleConversionFactor);
		//		outfile.write((char*)&temp, 1);
		//	}
		//}
		//for (unsigned int i = 0; i < linePaddingByteCount; ++i)
		//{
		//	outfile.write("", 1);
		//}
	}
}

//----------------------------------------------------------------------------------------
template<class SampleType>
float SampleToIRE(SampleType sampleValue, float ireLevel0, float ireLevel100)
{
	return ((float)sampleValue - ireLevel0) * (100.0f / (ireLevel100 - ireLevel0));
}

//----------------------------------------------------------------------------------------
template<class SampleType>
SampleType IREToSample(float ire, float ireLevel0, float ireLevel100)
{
	return (SampleType)(((ire * ((ireLevel100 - ireLevel0) / 100.0f)) + ireLevel0) + (std::numeric_limits<T>::is_integer ? 0.5f : 0.0f));
}

//----------------------------------------------------------------------------------------
//##TODO##
//template<class SampleType>
//bool DetectSyncBurst()
//{
//}

//----------------------------------------------------------------------------------------
template<class SampleType>
void WriteFramesToFiles(const std::wstring& outputFolderPath, const std::wstring& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::list<FrameBuilder::FrameInfo>& frameInfo, unsigned int threadCount = 0)
{
	size_t frameCount = frameInfo.size();

	// Determine the number of threads to use for this operation
	if (threadCount <= 0)
	{
		unsigned int coreCount = std::thread::hardware_concurrency();
		threadCount = (coreCount > 0) ? coreCount : 4;
	}
	if (threadCount > frameCount)
	{
		threadCount = 1;
	}

	std::vector<FrameBuilder::FrameInfo> frameInfoAsVector(frameInfo.cbegin(), frameInfo.cend());

	size_t chunkCount = threadCount;
	size_t framesPerChunk = frameCount / chunkCount;
	std::vector<std::thread> workerThreads;
	for (unsigned int i = 0; i < chunkCount; ++i)
	{
		workerThreads.emplace_back(std::thread([&, i]
		{
			size_t frameNo = i * framesPerChunk;
			size_t lastFrameNoForChunk = frameNo + framesPerChunk;
			while (frameNo < lastFrameNoForChunk)
			{
				const FrameBuilder::FrameInfo& frameEntry = frameInfoAsVector[frameNo];
				////##DEBUG##
				//if (frameNo == 4)
				//{
				//	for (auto entry : frameEntry.fieldInfo.front().syncEvents)
				//	{
				//		//std::cout << "Sync: " << entry.startSampleNo << "\t" << entry.endSampleNo << "\t" << entry.averageSyncLevel << '\n';
				//		std::wcout << "Sync: " << entry.startSampleNo << "\t" << entry.endSampleNo << "\t" << entry.endSampleNo - entry.startSampleNo << "\t" << entry.averageSyncLevel << '\n';
				//	}
				//}
				//else
				//{
				//	++frameNo;
				//	continue;
				//}

				std::wstring outputFilePath = outputFolderPath + L"\\" + outputFileNameBase + std::to_wstring(frameNo) + L".bmp";
				WriteFrameToFile(outputFilePath, sampleData, frameEntry);
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
void ProcessVideo(const std::wstring& inputFilePath, const std::wstring& outputFolderPath, const std::wstring& outputFileNameBase, bool useSlidingWindow)
{
	std::ifstream infile(inputFilePath, std::ios_base::binary);
	infile.seekg(0, infile.end);
	size_t length = infile.tellg();
	infile.seekg(0, infile.beg);

	//##DEBUG##
	//length = 700000;

	std::vector<SampleType> tempBuffer;
	tempBuffer.resize(length / sizeof(SampleType));
	infile.read((char*)&tempBuffer[0], tempBuffer.size() * sizeof(tempBuffer[0]));

	//##DEBUG##
	std::cout << "Detecting sync pulses\n";

	SyncDetector syncDetector;
	syncDetector.enableMinMaxSlidingWindow = useSlidingWindow;
	std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = syncDetector.DetectSyncPulses(tempBuffer);

	////##DEBUG##
	//std::map<size_t, std::pair<size_t, std::list<size_t>>> runCountsByLength;
	//for (auto entry : syncPulseInfo)
	//{
	//	auto& countEntry = runCountsByLength[entry.endSampleNo - entry.startSampleNo];
	//	++countEntry.first;
	//	countEntry.second.push_back(entry.startSampleNo);
	//}
	//for (auto entry : runCountsByLength)
	//{
	//	std::wcout << "Length: " << entry.first << "\tCount: " << entry.second.first;
	//	if (entry.second.first <= 5)
	//	{
	//		std::wcout << '\t';
	//		for (size_t startSampleNo : entry.second.second)
	//		{
	//			std::wcout << startSampleNo << ' ';
	//		}
	//	}
	//	std::wcout << '\n';
	//}

	//##DEBUG##
	//system("pause");

	//##DEBUG##
	std::cout << "Extracting sync events\n";

	std::list<SyncDetector::SyncInfo> syncInfo = syncDetector.DetectSyncEvents(tempBuffer, syncPulseInfo);

	////##DEBUG##
	//SyncDetector::SyncType lastSyncType;
	//size_t syncCount = 0;
	//for (auto entry : syncInfo)
	//{
	//	if ((syncCount > 0) && (entry.type == lastSyncType))
	//	{
	//		++syncCount;
	//	}
	//	else
	//	{
	//		if (syncCount > 0)
	//		{
	//			std::cout << "Sync: " << (int)lastSyncType << ", " << syncCount << '\n';
	//		}
	//		lastSyncType = entry.type;
	//		syncCount = 1;
	//	}
	//}
	//std::cout << "Sync: " << (int)lastSyncType << ", " << syncCount << '\n';

	//##DEBUG##
	std::cout << "Detecting frames\n";

	// Collect the sync events into frames
	FrameBuilder frameBuilder;
	std::list<FrameBuilder::FrameInfo> frameInfo = frameBuilder.DetectFrames(tempBuffer, syncInfo);

	//##DEBUG##
	//std::cout << "Frames:\n";
	//for (auto entry : frameInfo)
	//{
	//	std::cout << "lineCount: " << entry.fieldInfo.front().syncEvents.size() << "\tLength: " << entry.fieldInfo.front().syncEvents.back().endSampleNo - entry.fieldInfo.front().syncEvents.front().startSampleNo << "\tStart: " << entry.fieldInfo.front().syncEvents.front().startSampleNo << "\tEnd: " << entry.fieldInfo.back().syncEvents.back().startSampleNo << "\n";
	//}

	//##DEBUG##
	std::cout << "Writing frames\n";
	auto start = std::chrono::high_resolution_clock::now();

	// Write frames to files
	WriteFramesToFiles(outputFolderPath, outputFileNameBase, tempBuffer, frameInfo);

	//##DEBUG##
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Processing time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << '\n';
}

//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
	std::wstring inputFilePath;
	std::wstring outputFolderPath;
	std::wstring outputFileNameBase;
	bool useSlidingWindow = false;

	//##DEBUG##
	//inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\CompExternalRunFailure.bin";
	//inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\FantasiaCompositeSignedShort.bin";
	//inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\CompExternalVSyncMerge.bin";
	//inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\FinalTrimTest.bin";
	
	inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\WindowsDecode\\FantasiaCompositeSigned.bin";
	outputFolderPath = L"D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputFantasia";
	outputFileNameBase = L"Fantasia";
	ProcessVideo<short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow);

	//inputFilePath = L"D:\\Emulation\\Roms\\MegaLD\\CompExternal.raw";
	//outputFolderPath = L"D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputSonic";
	//outputFileNameBase = L"SonicTest";
	//useSlidingWindow = true;
	//ProcessVideo<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow);

	//##DEBUG##
	system("pause");
	return 0;
}
