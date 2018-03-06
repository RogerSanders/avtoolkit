#include <iostream>
#include <sstream>
#include <cassert>

//----------------------------------------------------------------------------------------------------------------------
// Logging methods
//----------------------------------------------------------------------------------------------------------------------
void Logger::Log(Severity severity, const std::string& message) const
{
	LogInternal(severity, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Log(Severity severity, const std::string& formatString, Args... args) const
{
	LogInternal(severity, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Critical(const std::string& message) const
{
	LogInternal(Severity::Critical, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Critical(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Critical, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Error(const std::string& message) const
{
	LogInternal(Severity::Error, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Error(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Error, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Warning(const std::string& message) const
{
	LogInternal(Severity::Warning, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Warning(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Warning, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Info(const std::string& message) const
{
	LogInternal(Severity::Info, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Info(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Info, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Trace(const std::string& message) const
{
	LogInternal(Severity::Trace, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Trace(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Trace, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::Debug(const std::string& message) const
{
	LogInternal(Severity::Debug, message);
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::Debug(const std::string& formatString, Args... args) const
{
	LogInternal(Severity::Debug, formatString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::LogInternal(Severity severity, const std::string& message) const
{
	// Since we're about to use stdout, take a lock to make this operation thread-safe from this same logger.
	std::lock_guard<std::mutex> lock(_accessMutex);

	// Write the resolved string to stdout
	std::cout << message << std::endl;
}

//----------------------------------------------------------------------------------------------------------------------
template<class... Args>
void Logger::LogInternal(Severity severity, const std::string& formatString, Args... args) const
{
	// Convert all supplied arguments to a string representation
	const size_t argCount = sizeof...(Args);
	std::string argsResolved[argCount];
	ResolveArgs(argsResolved, args...);

	// Resolve the format string and argument values down to a single string
	std::string messageResolved = ResolveFormatString(formatString, argCount, argsResolved);

	// Log the resolved string with the requested severity level
	LogInternal(severity, messageResolved);
}

//----------------------------------------------------------------------------------------------------------------------
// Format string methods
//----------------------------------------------------------------------------------------------------------------------
std::string Logger::ResolveFormatString(const std::string& formatString, size_t argCount, const std::string* formatStringArgs) const
{
	// Allocate a variable to hold our resolved string content
	std::string resolvedString;
	resolvedString.reserve(formatString.size() * 2);

	// Parse the format string, storing literal text and argument values in our resolved string.
	size_t formatStringCurrentPos = 0;
	while (formatStringCurrentPos < formatString.size())
	{
		// Attempt to find the next argument insert position
		size_t formatStringStartPos = formatString.find_first_of('{', formatStringCurrentPos);
		if (formatStringStartPos == std::string::npos)
		{
			break;
		}
		size_t formatStringEndPos = formatString.find_first_of('}', formatStringStartPos);
		if (formatStringEndPos == std::string::npos)
		{
			break;
		}

		// Retrieve the argument index number
		size_t argIndex;
		std::string argIndexAsString = formatString.substr(formatStringStartPos + 1, formatStringEndPos - (formatStringStartPos + 1));
		argIndex = (size_t)std::stoi(argIndexAsString);

		// Ensure the argument index number is valid
		assert(argIndex < argCount);

		// Append any text leading up to this argument in the format string, and the string representation of the target
		// argument.
		resolvedString.append(formatString.data() + formatStringCurrentPos, formatStringStartPos - formatStringCurrentPos);
		resolvedString.append(formatStringArgs[argIndex]);

		// Advance past the argument in the format string
		formatStringCurrentPos = formatStringEndPos + 1;
	}
	resolvedString.append(formatString.data() + formatStringCurrentPos);

	// Return the resolved string to the caller
	return resolvedString;
}

//----------------------------------------------------------------------------------------------------------------------
template<class T>
void Logger::ResolveArgs(std::string* argAsString, T arg) const
{
	ResolveArg(arg, *argAsString);
}

//----------------------------------------------------------------------------------------------------------------------
template<class T, class... Args>
void Logger::ResolveArgs(std::string* argAsString, T arg, Args... args) const
{
	ResolveArg(arg, *(argAsString++));
	ResolveArgs(argAsString, args...);
}

//----------------------------------------------------------------------------------------------------------------------
template<class T>
void Logger::ResolveArg(T arg, std::string& argResolved) const
{
	std::ostringstream stringStream;
	stringStream << std::fixed << arg;
	argResolved = stringStream.str();
}

//----------------------------------------------------------------------------------------------------------------------
void Logger::ResolveArg(const std::wstring& arg, std::string& argResolved) const
{
	argResolved = std::string(arg.begin(), arg.end());
}
