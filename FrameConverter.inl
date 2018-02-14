//----------------------------------------------------------------------------------------
// Frame conversion methods
//----------------------------------------------------------------------------------------
template<class SampleType>
void FrameConverter::WriteFrameToFile(const PathString& outputFilePath, const std::vector<SampleType>& sampleData, const FrameBuilder::FrameInfo& frameInfo, unsigned int pixelsPerLine) const
{
	const FrameBuilder::FieldInfo& fieldEntry = frameInfo.fieldInfo.front();
	unsigned int lineCount = fieldEntry.lineCount;

	unsigned int bitsPerPixel = 24;
	unsigned int lineByteCount = pixelsPerLine * 3;

	//Calculate the amount of padding on each line. Lines are padded out to DWORD
	//boundaries.
	unsigned int linePaddingByteCount = 0;
	if((lineByteCount % sizeof(unsigned int)) != 0)
	{
		linePaddingByteCount = sizeof(unsigned int) - (lineByteCount % sizeof(unsigned int));
	}

	//Calculate the totals and offsets we need to write to the file
	unsigned int pixelDataOffset = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);
	unsigned int pixelDataSize = (unsigned int)((lineByteCount + linePaddingByteCount) * lineCount);
	unsigned int fileSize = pixelDataOffset + pixelDataSize;

	//Write the bitmap file header to the file
	BitmapFileHeader fileHeader = {};
	char typeByte1 = 'B';
	char typeByte2 = 'M';
	fileHeader.bfType = ((unsigned short)typeByte2 << 8) | (unsigned short)typeByte1;
	fileHeader.bfSize = fileSize;
	fileHeader.bfReserved1 = 0;
	fileHeader.bfReserved2 = 0;
	fileHeader.bfOffBits = pixelDataOffset;

	//Write the bitmap info header to the file
	BitmapInfoHeader bitmapHeader = {};
	bitmapHeader.biSize = sizeof(BitmapInfoHeader);
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

		tk::spline lineAsSpline;
		lineAsSpline.set_points(x, y);

		double preciseLineStartPos;
		double preciseLineEndPos;
		FindPreciseLineStartEndPos(syncEvent, nextSyncEvent, lineAsSpline, preciseLineStartPos, preciseLineEndPos);
		double preciseLineWidth = (preciseLineEndPos - preciseLineStartPos);

		//##DEBUG##
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
				samplePos = preciseLineStartPos + ((double)(i - blankingLeadingSampleCount) * (preciseLineWidth / activeImageSampleCount));
			}
			else
			{
				samplePos = preciseLineEndPos + ((double)(i - (blankingLeadingSampleCount + activeImageSampleCount)) * ((preciseLineStartPos / blankingLeadingPercentage) / blankingFollowingSampleCount));
			}
			outputData[i] = (SampleType)lineAsSpline(samplePos);
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
// IRE conversion methods
//----------------------------------------------------------------------------------------
template<class SampleType>
float FrameConverter::SampleToIRE(SampleType sampleValue, float ireLevel0, float ireLevel100) const
{
	return ((float)sampleValue - ireLevel0) * (100.0f / (ireLevel100 - ireLevel0));
}

//----------------------------------------------------------------------------------------
template<class SampleType>
SampleType FrameConverter::IREToSample(float ire, float ireLevel0, float ireLevel100) const
{
	return (SampleType)(((ire * ((ireLevel100 - ireLevel0) / 100.0f)) + ireLevel0) + (std::numeric_limits<T>::is_integer ? 0.5f : 0.0f));
}

//----------------------------------------------------------------------------------------
// Resampling methods
//----------------------------------------------------------------------------------------
template<class SampleType>
void FrameConverter::ResampleLinear(const std::vector<SampleType>& inputData, std::vector<SampleType>& outputData) const
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
		outputData[xpos] = (SampleType)((finalSample * (float)std::numeric_limits<SampleType>::max()) + (std::numeric_limits<SampleType>::is_integer ? 0.5f : 0.0f));
	}
}

//----------------------------------------------------------------------------------------
template<class SampleType>
void FrameConverter::ResampleCubic(const std::vector<SampleType>& inputData, std::vector<SampleType>& outputData) const
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
		outputData[i] = (SampleType)s(samplePos);
	}
}
