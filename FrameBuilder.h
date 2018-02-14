#ifndef __FRAMEBUILDER_H__
#define __FRAMEBUILDER_H__
#include <list>
#include "SyncDetector.h"

class FrameBuilder
{
public:
	//Structures
	struct FieldInfo
	{
		std::list<SyncDetector::SyncInfo> syncEvents;
		bool oddFrame;
		int lineCount;
		size_t endSampleNo;
		SyncDetector::SyncInfo followingSyncEvent;
	};
	struct FrameInfo
	{
		std::list<FieldInfo> fieldInfo;
	};

public:
	// Frame detection methods
	template<class SampleType>
	std::list<FrameInfo> DetectFrames(const std::vector<SampleType>& sampleData, const std::list<SyncDetector::SyncInfo>& syncInfo) const;
};

#include "FrameBuilder.inl"
#endif
