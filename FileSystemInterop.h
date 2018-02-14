#ifndef __FILESYSTEMINTEROP_H__
#define __FILESYSTEMINTEROP_H__
#include <string>
#include <codecvt>

// Until we have C++17's filesystem library available, we handle platform-specific filesystem issues here.
#if defined(_WIN32) && defined(_UNICODE)
typedef std::wstring PathString;
typedef wchar_t PathChar;
#else
typedef std::string PathString;
typedef char PathChar;
#endif

#if defined(_WIN32) && defined(_UNICODE)
static const wchar_t PathSeparatorChar = L'\\';
static const wchar_t* PathSeparatorString = L"\\";
#else
static const char PathSeparatorChar = '\\';
static const wchar_t* PathSeparatorString = "\\";
#endif

#if defined(_WIN32) && defined(_UNICODE)
inline PathString ToPathString(const std::string& pathString)
{
	std::wstring_convert<std::codecvt_utf8<wchar_t>> conversion;
	return conversion.from_bytes(pathString);
}

inline PathString ToPathString(const std::wstring& pathString)
{
	return pathString;
}
#else
inline std::string ToPathString(const std::string& pathString)
{
	return pathString;
}

inline std::string ToPathString(const std::wstring& pathString)
{
	std::wstring_convert<std::codecvt_utf8<wchar_t>> conversion;
	return conversion.to_bytes(pathString);
}
#endif

#endif
