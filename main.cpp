#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "FileSystemInterop.h"
#include "SyncDetector.h"
#include "FrameBuilder.h"
#include "FrameConverter.h"

//----------------------------------------------------------------------------------------
template<class SampleType>
void WriteFramesToFiles(const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::list<FrameBuilder::FrameInfo>& frameInfo, unsigned int threadCount = 0)
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

	FrameConverter frameConverter;

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

				PathString outputFilePath = outputFolderPath + PathSeparatorChar + outputFileNameBase + ToPathString(std::to_string(frameNo) + ".bmp");
				frameConverter.WriteFrameToFile(outputFilePath, sampleData, frameEntry, 910);
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
void ProcessVideo(const PathString& inputFilePath, const PathString& outputFolderPath, const PathString& outputFileNameBase, bool useSlidingWindow)
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
#if defined(_WIN32) && defined(_UNICODE)
int wmain(int argc, PathChar* argv[])
#else
int main(int argc, PathChar* argv[])
#endif
{
	PathString inputFilePath;
	PathString outputFolderPath;
	PathString outputFileNameBase;
	bool useSlidingWindow = false;

	//##DEBUG##
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternalRunFailure.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\FantasiaCompositeSignedShort.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternalVSyncMerge.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\FinalTrimTest.bin");

	inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\WindowsDecode\\FantasiaCompositeSigned.bin");
	outputFolderPath = ToPathString(L"D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputFantasia");
	outputFileNameBase = ToPathString("Fantasia");
	ProcessVideo<short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow);

	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternal.raw");
	//outputFolderPath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputSonic");
	//outputFileNameBase = ToPathString("SonicTest");
	//useSlidingWindow = true;
	//ProcessVideo<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow);

	return 0;
}
