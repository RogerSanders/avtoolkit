#include <iostream>
#include <string>
#include <vector>
#include "Logger.h"
#include "StringHelpers.h"
#include "FileSystemInterop.h"
#include "VideoDecoder.h"

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
#if defined(_WIN32) && defined(_UNICODE)
int wmain(int argc, PathChar* argv[])
#else
int main(int argc, PathChar* argv[])
#endif
{
	// Allocate a logger for this program
	Logger logger;

	// Set our various options to their initial state
	PathString inputFilePath;
	PathString outputFolderPath;
	PathString outputFileNameBase;
	bool useSlidingWindow = false;
	bool outputHelp = false;
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

	// Create and configure the video decoder
	SyncDetector syncDetector(logger);
	FrameBuilder frameBuilder(logger);
	VideoDecoder videoDecoder(syncDetector, frameBuilder, logger);
	syncDetector.enableMinMaxSlidingWindow = useSlidingWindow;

	// Mark the current time so we can calculate our total decode time
	auto startTime = std::chrono::high_resolution_clock::now();

	// Process the target composite video sample data
	//##TODO## These types are right for the sizes on Windows/Linux, but typedef them into guaranteed sizes in a header.
	//##TODO## Byte ordering is platform dependent currently. Consider if we want to formally support big endian and
	//little endian ordering in a platform independent manner.
	switch (sampleType)
	{
	case SampleType::Int8:
		videoDecoder.ConvertCompositeVideoToImages<signed char>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::UInt8:
		videoDecoder.ConvertCompositeVideoToImages<unsigned char>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::Int16:
		videoDecoder.ConvertCompositeVideoToImages<short>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::UInt16:
		videoDecoder.ConvertCompositeVideoToImages<unsigned short>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::Int32:
		videoDecoder.ConvertCompositeVideoToImages<int>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::UInt32:
		videoDecoder.ConvertCompositeVideoToImages<unsigned int>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::Float32:
		videoDecoder.ConvertCompositeVideoToImages<float>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	case SampleType::Float64:
		videoDecoder.ConvertCompositeVideoToImages<double>(inputFilePath, outputFolderPath, outputFileNameBase);
		break;
	}

	// Log how long it took to decode the target file
	auto endTime = std::chrono::high_resolution_clock::now();
	auto totalTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	logger.Trace("Total time: {0}", totalTimeInMilliseconds.count());
	return 0;
}
