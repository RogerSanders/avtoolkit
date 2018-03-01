#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "FileSystemInterop.h"
#include "SyncDetector.h"
#include "FrameBuilder.h"
#include "FrameConverter.h"
#include "SplineHelpers.h"

//----------------------------------------------------------------------------------------
template<class SampleType>
void WriteFramesToFiles(const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::vector<FrameBuilder::FrameInfo>& frames, unsigned int lineWidthInPixels, unsigned int threadCount = 0)
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
				const FrameBuilder::FrameInfo& frameEntry = frames[frameNo];
				PathString outputFilePath = outputFolderPath + PathSeparatorChar + outputFileNameBase + ToPathString(std::to_string(frameNo) + ".bmp");
				frameConverter.WriteFrameToBMP(outputFilePath, sampleData, frameEntry, lineWidthInPixels);
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
void ProcessVideo(const PathString& inputFilePath, const PathString& outputFolderPath, const PathString& outputFileNameBase, bool useSlidingWindow, unsigned int lineWidthInPixels)
{
	//##TODO## Read the file in chunks (default 1GB) and combine chunks together
	std::ifstream infile(inputFilePath, std::ios_base::binary);
	infile.seekg(0, infile.end);
	size_t fileLength = infile.tellg();
	infile.seekg(0, infile.beg);

	//##DEBUG##
	//fileLength = 70000000;
	//fileLength = 500000000;

	std::vector<SampleType> fileData;
	fileData.resize(fileLength / sizeof(SampleType));
	infile.read((char*)&fileData[0], fileData.size() * sizeof(fileData[0]));

	//##DEBUG##
	std::cout << "Detecting sync pulses\n";

	SyncDetector syncDetector;
	syncDetector.enableMinMaxSlidingWindow = useSlidingWindow;
	std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = syncDetector.DetectSyncPulses(fileData);

	//##DEBUG##
	std::cout << "Extracting sync events\n";
	std::list<SyncDetector::SyncInfo> syncInfo = syncDetector.DetectSyncEvents(fileData, syncPulseInfo);

	//##DEBUG##
	std::cout << "Detecting frames\n";
	FrameBuilder frameBuilder;
	//##DEBUG##
	//frameBuilder.combineInterlacedFields = false;

	// Collect the sync events into frames
	std::list<FrameBuilder::FieldInfo> fields = frameBuilder.DetectFields(fileData, syncInfo);

	//##DEBUG##
	std::cout << "Detecting lines\n";
	std::vector<FrameBuilder::FrameInfo> frames = frameBuilder.DetectFrames(fields);
	frameBuilder.DetectLines(fileData, frames);

	//##DEBUG##
	std::cout << "Writing frames\n";

	// Write frames to files
	WriteFramesToFiles(outputFolderPath, outputFileNameBase, fileData, frames, lineWidthInPixels);
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
	unsigned int lineWidthInPixels = 930; //910 //3640 //7280

	//##DEBUG##
	auto start = std::chrono::high_resolution_clock::now();

	//##DEBUG##
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternalRunFailure.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\FantasiaCompositeSignedShort.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternalVSyncMerge.bin");
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\FinalTrimTest.bin");

	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\WindowsDecode\\FantasiaCompositeSigned.bin");
	//outputFolderPath = ToPathString(L"D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputFantasiaInterlaced3");
	//outputFileNameBase = ToPathString("Fantasia");
	//ProcessVideo<short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);

	inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternal.raw");
	outputFolderPath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputSonic3");
	outputFileNameBase = ToPathString("SonicTest");
	useSlidingWindow = true;
	ProcessVideo<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);

	//##DEBUG##
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Total time: " << (end - start).count() << "\n";

	//##DEBUG##
	system("pause");
	return 0;
}
