#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "Logger.h"
#include "StringHelpers.h"
#include "FileSystemInterop.h"
#include "SyncDetector.h"
#include "FrameBuilder.h"
#include "FrameConverter.h"

//----------------------------------------------------------------------------------------
template<class SampleType>
void WriteFramesToFiles(const Logger& logger, const PathString& outputFolderPath, const PathString& outputFileNameBase, const std::vector<SampleType>& sampleData, const std::vector<FrameBuilder::FrameInfo>& frames, unsigned int lineWidthInPixels, size_t initialFrameNo, unsigned int threadCount = 0)
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

	FrameConverter frameConverter(logger);
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
//##TODO## Wrap this up in a helper class
template<class SampleType>
void ConvertCompositeVideoToImages(const PathString& inputFilePath, const PathString& outputFolderPath, const PathString& outputFileNameBase, bool useSlidingWindow, unsigned int lineWidthInPixels, size_t maxChunkSizeInBytes)
{
	// Allocate a logger for this operation
	Logger logger;

	// Mark the current time so we can calculate our total decode time
	auto startTime = std::chrono::high_resolution_clock::now();

	// Open the target file, and calculate the size of the file in bytes.
	std::ifstream infile(inputFilePath, std::ios_base::binary);
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
		logger.Info("Reading {0} bytes from file pos {1}", chunkSizeInBytes, currentFilePos);
		infile.read((char*)&fileData[fileReadStartPos], chunkSizeInBytes);
		currentFilePos += chunkSizeInBytes;

		// Detect sync pulses in the sample data
		logger.Info("Detecting sync pulses");
		SyncDetector syncDetector(logger);
		syncDetector.enableMinMaxSlidingWindow = useSlidingWindow;
		std::list<SyncDetector::SyncPulseInfo> syncPulseInfo = syncDetector.DetectSyncPulses(fileData, sampleScanningStartPos);

		// Build sync events from the raw sync pulses
		logger.Info("Extracting sync events");
		std::list<SyncDetector::SyncInfo> syncInfo = syncDetector.DetectSyncEvents(fileData, syncPulseInfo);

		// Collect the sync events into frames
		logger.Info("Detecting frames");
		FrameBuilder frameBuilder(logger);
		std::list<FrameBuilder::FieldInfo> fields = frameBuilder.DetectFields(fileData, syncInfo);
		std::vector<FrameBuilder::FrameInfo> frames = frameBuilder.DetectFrames(fields);

		// Decode line information within each detected frame
		logger.Info("Detecting lines");
		frameBuilder.DetectLines(fileData, frames);

		// Write each detected frame out to an image file
		logger.Info("Writing frames");
		WriteFramesToFiles(logger, outputFolderPath, outputFileNameBase, fileData, frames, lineWidthInPixels, startFrameNo);
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

	// Log how long it took to decode the target file
	auto endTime = std::chrono::high_resolution_clock::now();
	auto totalTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	logger.Trace("Total time: {0}", totalTimeInMilliseconds.count());
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
	unsigned int lineWidthInPixels = 930;
	SampleType sampleType;
	//##FIX## Make this configurable
	size_t maxChunkSizeInBytes = 1024*1024*1024;

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
		ConvertCompositeVideoToImages<signed char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::UInt8:
		ConvertCompositeVideoToImages<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::Int16:
		ConvertCompositeVideoToImages<short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::UInt16:
		ConvertCompositeVideoToImages<unsigned short>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::Int32:
		ConvertCompositeVideoToImages<int>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::UInt32:
		ConvertCompositeVideoToImages<unsigned int>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::Float32:
		ConvertCompositeVideoToImages<float>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	case SampleType::Float64:
		ConvertCompositeVideoToImages<double>(inputFilePath, outputFolderPath, outputFileNameBase, useSlidingWindow, lineWidthInPixels, maxChunkSizeInBytes);
		break;
	}
	return 0;
}
