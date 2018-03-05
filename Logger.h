#ifndef __LOGGER_H__
#define __LOGGER_H__
#include <string>
#include <mutex>

class Logger
{
public:
	// Enumerations
	enum class Severity
	{
		Critical,
		Error,
		Warning,
		Info,
		Trace,
		Debug,
	};

public:
	// Logging methods
	inline void Log(Severity severity, const std::string& message) const;
	template<class... Args>
	void Log(Severity severity, const std::string& formatString, Args... args) const;
	inline void Critical(const std::string& message) const;
	template<class... Args>
	void Critical(const std::string& formatString, Args... args) const;
	inline void Error(const std::string& message) const;
	template<class... Args>
	void Error(const std::string& formatString, Args... args) const;
	inline void Warning(const std::string& message) const;
	template<class... Args>
	void Warning(const std::string& formatString, Args... args) const;
	inline void Info(const std::string& message) const;
	template<class... Args>
	void Info(const std::string& formatString, Args... args) const;
	inline void Trace(const std::string& message) const;
	template<class... Args>
	void Trace(const std::string& formatString, Args... args) const;
	inline void Debug(const std::string& message) const;
	template<class... Args>
	void Debug(const std::string& formatString, Args... args) const;

private:
	// Logging methods
	inline void LogInternal(Severity severity, const std::string& message) const;
	template<class... Args>
	void LogInternal(Severity severity, const std::string& formatString, Args... args) const;

	// Format string methods
	inline std::string ResolveFormatString(const std::string& formatString, size_t argCount, const std::string* formatStringArgs) const;
	template<class T>
	void ResolveArgs(std::string* argAsString, T arg) const;
	template<class T, class... Args>
	void ResolveArgs(std::string* argAsString, T arg, Args... args) const;
	template<class T>
	void ResolveArg(T arg, std::string& argResolved) const;

private:
	mutable std::mutex _accessMutex;
};

#include "Logger.inl"
#endif
