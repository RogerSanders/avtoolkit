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
void WriteFramesToFiles(const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::vector<FrameBuilder::FrameInfo>& frames, unsigned int lineWidthInPixels, size_t initialFrameNo, unsigned int threadCount = 0)
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
	//##DEBUG##
	auto start = std::chrono::high_resolution_clock::now();

	//##FIX## Make this configurable
	//##DEBUG##
	const size_t maxChunkSizeInBytes = 1024*1024*1024;
	//const size_t maxChunkSizeInBytes = (size_t)(1024*1024*1024) * 4;

	//##TODO## Read the file in chunks (default 1GB) and combine chunks together
	std::ifstream infile(inputFilePath, std::ios_base::binary);
	infile.seekg(0, infile.end);
	size_t fileLengthInBytes = infile.tellg();
	infile.seekg(0, infile.beg);

	//##DEBUG##
	//fileLengthInBytes = 70000000;
	//fileLengthInBytes = 500000000;

	std::vector<SampleType> fileData;
	size_t currentFilePos = 0;
	bool firstPass = true;
	size_t sampleCountInChunk = 0;
	size_t continuePos;
	size_t startFrameNo = 0;
	while (currentFilePos < fileLengthInBytes)
	{
		size_t chunkSizeInBytes = (fileLengthInBytes - currentFilePos);
		chunkSizeInBytes = (chunkSizeInBytes > maxChunkSizeInBytes) ? maxChunkSizeInBytes : chunkSizeInBytes;
		sampleCountInChunk = (chunkSizeInBytes / sizeof(SampleType));

		size_t fileReadStartPos;
		size_t sampleScanningStartPos;
		if (firstPass)
		{
			fileData.resize(sampleCountInChunk);
			sampleScanningStartPos = 0;
			fileReadStartPos = 0;
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

		//##DEBUG##
		std::cout << "Reading " << chunkSizeInBytes << " bytes from file pos " << currentFilePos << "\n";

		// Read in new data from the input file
		infile.read((char*)&fileData[fileReadStartPos], chunkSizeInBytes);
		currentFilePos += chunkSizeInBytes;

		//##DEBUG##
		std::cout << "Detecting sync pulses\n";

		SyncDetector syncDetector;
		syncDetector.enableMinMaxSlidingWindow = useSlidingWindow;
		std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = syncDetector.DetectSyncPulses(fileData, sampleScanningStartPos);

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
		//frameBuilder.DetectLines(fileData, frames, 1);

		//##DEBUG##
		std::cout << "Writing frames\n";

		// Write frames to files
		WriteFramesToFiles(outputFolderPath, outputFileNameBase, fileData, frames, lineWidthInPixels, startFrameNo);
		startFrameNo += frames.size();

		const FrameBuilder::FieldInfo& lastFieldInfo = frames.back().fields.back();
		const SyncDetector::SyncInfo& targetSyncInfo = lastFieldInfo.syncEvents[lastFieldInfo.syncEvents.size() / 2];
		continuePos = targetSyncInfo.startSampleNo;
		firstPass = false;
	}

	//##DEBUG##
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Total time: " << (end - start).count() << "\n";
}

//----------------------------------------------------------------------------------------
bool StringStartsWith(const PathString& targetString, const PathString& searchString, bool caseInsensitive = false)
{
	if (targetString.size() < searchString.size())
	{
		return false;
	}

	if (caseInsensitive)
	{
		for (size_t i = 0; i < searchString.size(); ++i)
		{
			if (toupper(targetString[i]) != toupper(searchString[i]))
			{
				return false;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < searchString.size(); ++i)
		{
			if (targetString[i] != searchString[i])
			{
				return false;
			}
		}
	}
	return true;
}

//----------------------------------------------------------------------------------------
bool StringEquals(const PathString& value1, const PathString& value2, bool caseInsensitive = false)
{
	if (value1.size() != value2.size())
	{
		return false;
	}

	if (caseInsensitive)
	{
		for (size_t i = 0; i < value1.size(); ++i)
		{
			if (toupper(value1[i]) != toupper(value2[i]))
			{
				return false;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < value1.size(); ++i)
		{
			if (value1[i] != value2[i])
			{
				return false;
			}
		}
	}
	return true;
}

//----------------------------------------------------------------------------------------
bool ExtractOptionNameFromArg(const std::vector<PathString>& optionPrefixes, const PathString& arg, PathString& optionName)
{
	for (const PathString& prefix : optionPrefixes)
	{
		if (StringStartsWith(arg, prefix))
		{
			optionName = arg.substr(prefix.size());
			return true;
		}
	}
	return false;
}

//----------------------------------------------------------------------------------------
enum class SampleType
{
	Int8,
	UInt8,
	Int16,
	UInt16,
	Int32,
	UInt32,
	Float32,
	Float64,
};

//----------------------------------------------------------------------------------------
#if defined(_WIN32) && defined(_UNICODE)
int wmain(int argc, PathChar* argv[])
#else
int main(int argc, PathChar* argv[])
#endif
{
	// Set our various options to their initial state
	PathString inputFilePath;
	PathString outputFolderPath;
	PathString outputFileNameBase;
	bool useSlidingWindow = false;
	bool outputHelp = false;
	unsigned int lineWidthInPixels = 930; //910 //3640 //7280
	SampleType sampleType;

	// Process our command line options
	std::vector<PathString> argPrefixes = { ToPathString("/"), ToPathString("--"), ToPathString("-") };
	int currentArg = 1;
	PathString optionName;
	while ((currentArg < argc) && ExtractOptionNameFromArg(argPrefixes, argv[currentArg], optionName))
	{
		++currentArg;
		if (StringEquals(optionName, ToPathString("?")) || StringEquals(optionName, ToPathString("h"), true) || StringEquals(optionName, ToPathString("help"), true))
		{
			outputHelp = true;
		}
		else if (StringEquals(optionName, ToPathString("i"), true) && (currentArg < argc))
		{
			inputFilePath = argv[currentArg++];
		}
		else if (StringEquals(optionName, ToPathString("o"), true) && (currentArg < argc))
		{
			outputFolderPath = argv[currentArg++];
		}
		else if (StringEquals(optionName, ToPathString("b"), true) && (currentArg < argc))
		{
			outputFileNameBase = argv[currentArg++];
		}
		else if (StringEquals(optionName, ToPathString("s"), true))
		{
			useSlidingWindow = true;
		}
		else if (StringEquals(optionName, ToPathString("t"), true) && (currentArg < argc))
		{
			PathString sampleTypeAsString = argv[currentArg++];
			if (StringEquals(sampleTypeAsString, ToPathString("Int8"), true))
			{
				sampleType = SampleType::Int8;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("UInt8"), true))
			{
				sampleType = SampleType::UInt8;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("Int16"), true))
			{
				sampleType = SampleType::Int16;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("UInt16"), true))
			{
				sampleType = SampleType::UInt16;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("Int32"), true))
			{
				sampleType = SampleType::Int32;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("UInt32"), true))
			{
				sampleType = SampleType::UInt32;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("Float32"), true))
			{
				sampleType = SampleType::Float32;
			}
			else if (StringEquals(sampleTypeAsString, ToPathString("Float64"), true))
			{
				sampleType = SampleType::Float64;
			}
		}
	}

	//##DEBUG##
	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\WindowsDecode\\FantasiaCompositeSigned.bin");
	//outputFolderPath = ToPathString(L"D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputFantasiaInterlaced3");
	//outputFileNameBase = ToPathString("Fantasia");
	//sampleType = SampleType::Int16;

	//inputFilePath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\CompExternal.raw");
	//outputFolderPath = ToPathString("D:\\Emulation\\Roms\\MegaLD\\TestVideoOutputSonic4");
	//outputFileNameBase = ToPathString("SonicTest");
	//sampleType = SampleType::UInt8;
	//useSlidingWindow = true;

	// Validate our command line options
	if (outputHelp || inputFilePath.empty() || outputFolderPath.empty())
	{
		std::cout << "Processes raw composite video data and converts it into other forms\n"
		             "\n"
		             "CompositeVideoDecode /i inputFilePath /o outputPath [/t sampleType] [/b baseOutputFileName] [/s]\n"
		             "\n"
		             "/i     Specifies the input file containing raw composite sample data\n"
		             "/o     Specifies the output folder or file (depending other on options) to write the data to\n"
		             "/b     When outputting multiple files to a folder, specifies the base name of the output files\n"
		             "/t     Specifies the sample data type in the input file\n"
		             "          Int8    - 8-bit signed integer\n"
		             "          UInt8   - 8-bit unsigned integer (default)\n"
		             "          Int16   - 16-bit signed integer\n"
		             "          UInt16  - 16-bit unsigned integer\n"
		             "          Int32   - 32-bit signed integer\n"
		             "          UInt32  - 32-bit unsigned integer\n"
		             "          Float32 - 32-bit floating point\n"
		             "          Float64 - 64-bit floating point\n"
		             "/s     Perform auto-ranging on the input sample data using an adaptive sliding window\n";
		return 0;
	}

	// Process the target composite video sample data
	//##TODO## These types are right for the sizes on Windows/Linux, but typedef them into guaranteed sizes in a header.
	//##TODO## Byte ordering is platform dependent currently. Consider if we want to formally support big endian and
	//little endian ordering in a platform independent manner.
	switch (sampleType)
	{
	case SampleType::Int8:
		ProcessVideo<signed char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::UInt8:
		ProcessVideo<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::Int16:
		ProcessVideo<short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::UInt16:
		ProcessVideo<unsigned short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::Int32:
		ProcessVideo<int>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::UInt32:
		ProcessVideo<unsigned int>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::Float32:
		ProcessVideo<float>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	case SampleType::Float64:
		ProcessVideo<double>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels);
		break;
	}

	//##DEBUG##
	system("pause");
	return 0;
}
