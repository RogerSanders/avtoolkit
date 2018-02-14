//----------------------------------------------------------------------------------------------------------------------
// Frame detection methods
//----------------------------------------------------------------------------------------------------------------------
template<class T>
std::list<FrameBuilder::FrameInfo> FrameBuilder::DetectFrames(const std::vector<T>& sampleData, const std::list<SyncDetector::SyncInfo>& syncInfo) const
{
	// Build a list of fields in the current frame
	std::list<FieldInfo> fields;
	FieldInfo currentField;
	bool currentFieldPopulated = false;
	bool foundHSyncInCurrentField = false;
	unsigned int lineCount = 0;
	for (auto entry : syncInfo)
	{
		if ((entry.type == SyncDetector::SyncType::Vertical) && foundHSyncInCurrentField)
		{
			if (currentFieldPopulated)
			{
				currentField.lineCount = lineCount;
				currentField.endSampleNo = entry.startSampleNo;
				currentField.followingSyncEvent = entry;
				fields.push_back(std::move(currentField));
			}
			currentField = FieldInfo();
			currentFieldPopulated = true;
			foundHSyncInCurrentField = false;
			lineCount = 0;
		}

		if (entry.type == SyncDetector::SyncType::Horizontal)
		{
			++lineCount;
		}

		foundHSyncInCurrentField |= (entry.type == SyncDetector::SyncType::Horizontal);

		if (!currentFieldPopulated)
		{
			continue;
		}

		currentField.syncEvents.push_back(entry);
	}

	//##TODO## Combine the fields into frames
	std::list<FrameInfo> frames;
	for (auto entry : fields)
	{
		FrameInfo frameInfo;
		frameInfo.fieldInfo.push_back(std::move(entry));
		frames.push_back(std::move(frameInfo));
	}

	// Return the list of frames to the caller
	return frames;
}
