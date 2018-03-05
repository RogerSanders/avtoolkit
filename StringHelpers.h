#ifndef __STRINGHELPERS_H__
#define __STRINGHELPERS_H__
#include <string>
#include <cctype>

//----------------------------------------------------------------------------------------
template<class ElementType, class traits, class Alloc>
bool StringStartsWith(const std::basic_string<ElementType, traits, Alloc>& targetString, const std::basic_string<ElementType, traits, Alloc>& searchString, bool caseInsensitive = false)
{
	if (targetString.size() < searchString.size())
	{
		return false;
	}

	if (caseInsensitive)
	{
		for (size_t i = 0; i < searchString.size(); ++i)
		{
			if (std::toupper(targetString[i]) != std::toupper(searchString[i]))
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
template<class ElementType, class traits, class Alloc>
bool StringEquals(const std::basic_string<ElementType, traits, Alloc>& value1, const std::basic_string<ElementType, traits, Alloc>& value2, bool caseInsensitive = false)
{
	if (value1.size() != value2.size())
	{
		return false;
	}

	if (caseInsensitive)
	{
		for (size_t i = 0; i < value1.size(); ++i)
		{
			if (std::toupper(value1[i]) != std::toupper(value2[i]))
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

#endif
