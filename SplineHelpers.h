#ifndef __SPLINEHELPERS_H__
#define __SPLINEHELPERS_H__
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>

//##TODO## Clean this up, and migrate into files/types where appropriate.
//----------------------------------------------------------------------------------------
// Based on code from https://stackoverflow.com/a/23980479/3418322
template<class T>
class CubicPolynomial
{
public:
	CubicPolynomial() = default;
	CubicPolynomial(T x0, T x1, T t0, T t1)
	:_x0(x0), _x1(x1)
	{
		//Compute coefficients for a cubic polynomial
		//  p(s) = c0 + c1*s + c2*s^2 + c3*s^3
		//such that
		//  p(0) = x0, p(1) = x1
		// and
		//  p'(0) = t0, p'(1) = t1.
		_c[0] = x0;
		_c[1] = t0;
		_c[2] = (-3*x0) + (3*x1) - (2*t0) - t1;
		_c[3] = (2*x0) - (2*x1) + t0 + t1;
	}

	T Evaluate(T t)
	{
		T t2 = t * t;
		T t3 = t2 * t;
		return _c[0] + (_c[1]*t) + (_c[2]*t2) + (_c[3]*t3);
	}
	T Reverse(T targetValue, T searchTolerance)
	{
		// Calculate our allowable deviance values from the specified target
		double maxSampleDeviance = (std::abs(_x1 - _x0)) * searchTolerance;
		double syncEndMinimumSampleValue = targetValue - maxSampleDeviance;
		double syncEndMaximumSampleValue = targetValue + maxSampleDeviance;

		// Perform an iterative search for a position between samples x1 and x2 where we're within tolerance of the target
		// value
		double searchPos;
		double searchPosHigh = 1.0;
		double searchPosLow = 0.0;
		if (_x1 < _x0)
		{
			do
			{
				searchPos = (searchPosLow + ((searchPosHigh - searchPosLow) / 2.0));
				double curveVal = Evaluate(searchPos);

				if (curveVal < syncEndMinimumSampleValue)
				{
					searchPosHigh = searchPos;
					continue;
				}
				else if (curveVal > syncEndMaximumSampleValue)
				{
					searchPosLow = searchPos;
					continue;
				}
				break;
			}
			while (true);
		}
		else
		{
			do
			{
				searchPos = (searchPosLow + ((searchPosHigh - searchPosLow) / 2.0));
				double curveVal = Evaluate(searchPos);

				if (curveVal < syncEndMinimumSampleValue)
				{
					searchPosLow = searchPos;
					continue;
				}
				else if (curveVal > syncEndMaximumSampleValue)
				{
					searchPosHigh = searchPos;
					continue;
				}
				break;
			}
			while (true);
		}

		// Return the position to the caller
		return searchPos;
	}

private:
	T _c[4];
	T _x0;
	T _x1;
};

//----------------------------------------------------------------------------------------
template<class T>
T DistanceSquared(T p, T q)
{
	T dx = q - p;
	return dx*dx;
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomUniform(T x0, T x1, T x2, T x3)
{
	// Catmull-Rom with tension 0.5
	T tension = (T)0.5;
	T t1 = tension * (x2 - x0);
	T t2 = tension * (x3 - x1);
	return CubicPolynomial<T>(x1, x2, t1, t2);
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomUniform(const T* data)
{
	return CreateSplineCatmullRomUniform(data[0], data[1], data[2], data[3]);
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomNonUniform(T x0, T x1, T x2, T x3, T dt0, T dt1, T dt2)
{
	// compute tangents when parameterized in [t1,t2]
	T t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
	T t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;

	// rescale tangents for parametrization in [0,1]
	t1 *= dt1;
	t2 *= dt1;

	return CubicPolynomial<T>(x1, x2, t1, t2);
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomNonUniform(const T* data, T parameterization)
{
	T factor = parameterization / (T)2.0;
	T dt0 = std::pow(DistanceSquared(data[0], data[1]), factor);
	T dt1 = std::pow(DistanceSquared(data[1], data[2]), factor);
	T dt2 = std::pow(DistanceSquared(data[2], data[3]), factor);

	//##FIX## Determine what to do with this
	//// safety check for repeated points
	//if (dt1 < 1e-4f)    dt1 = 1.0f;
	//if (dt0 < 1e-4f)    dt0 = dt1;
	//if (dt2 < 1e-4f)    dt2 = dt1;

	return CreateSplineCatmullRomNonUniform(data[0], data[1], data[2], data[3], dt0, dt1, dt2);
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomCentripetal(const T* data)
{
	return CreateSplineCatmullRomNonUniform(data, (T)0.5);
}

//----------------------------------------------------------------------------------------
template<class T>
CubicPolynomial<T> CreateSplineCatmullRomChordal(const T* data)
{
	return CreateSplineCatmullRomNonUniform(data, (T)1.0);
}

//----------------------------------------------------------------------------------------
template<class T>
void CubicInterpolateCatmullRom(const T* data, double startPos, double endPos, std::vector<T>& outputData)
{
	double inputSampleCount = endPos - startPos;
	size_t outputSampleCount = outputData.size();
	double outputSampleToInputSampleRatio = (inputSampleCount / (double)outputSampleCount);

	size_t currentOutputPos = 0;
	size_t currentInputPosInSamples = (size_t)startPos;
	CubicPolynomial<double> cubicPolynomial = CreateSplineCatmullRomUniform((double)data[currentInputPosInSamples - 1], (double)data[currentInputPosInSamples], (double)data[currentInputPosInSamples + 1], (double)data[currentInputPosInSamples + 2]);
	while (currentOutputPos < outputData.size())
	{
		double newInputPos = startPos + ((double)currentOutputPos * outputSampleToInputSampleRatio);
		size_t newInputPosInSamples = (size_t)newInputPos;
		if (currentInputPosInSamples != newInputPosInSamples)
		{
			currentInputPosInSamples = newInputPosInSamples;
			cubicPolynomial = CreateSplineCatmullRomUniform((double)data[currentInputPosInSamples - 1], (double)data[currentInputPosInSamples], (double)data[currentInputPosInSamples + 1], (double)data[currentInputPosInSamples + 2]);
		}
		outputData[currentOutputPos++] = (T)cubicPolynomial.Evaluate(newInputPos - (double)newInputPosInSamples);
	}
}

//----------------------------------------------------------------------------------------
template<class T>
T CubicInterpolateCatmullRomSearch(T x0, T x1, T x2, T x3, T targetValue, T searchTolerance)
{
	// Ensure we've been passed a target value which is within the range of x1 and x2
	assert(((x1 <= x2) && (targetValue >= x1) && (targetValue <= x2)) || ((x1 >= x2) && (targetValue <= x1) && (targetValue >= x2)));

	// Construct a cubic polynomial to perform repeated cubic interpolation operations between the sample data
	CubicPolynomial<T> cubicPolynomial = CreateSplineCatmullRomUniform(x0, x1, x2, x3);

	// Calculate our allowable deviance values from the specified target
	double maxSampleDeviance = (std::abs(x2 - x1)) * searchTolerance;
	double syncEndMinimumSampleValue = targetValue - maxSampleDeviance;
	double syncEndMaximumSampleValue = targetValue + maxSampleDeviance;

	// Perform an iterative search for a position between samples x1 and x2 where we're within tolerance of the target
	// value
	double searchPos;
	double searchPosHigh = 1.0;
	double searchPosLow = 0.0;
	if (x2 < x1)
	{
		while (true)
		{
			searchPos = (searchPosLow + ((searchPosHigh - searchPosLow) / 2.0));
			double curveVal = cubicPolynomial.Evaluate(searchPos);

			if (curveVal < syncEndMinimumSampleValue)
			{
				searchPosHigh = searchPos;
				continue;
			}
			else if (curveVal > syncEndMaximumSampleValue)
			{
				searchPosLow = searchPos;
				continue;
			}
			break;
		}
	}
	else
	{
		while (true)
		{
			searchPos = (searchPosLow + ((searchPosHigh - searchPosLow) / 2.0));
			double curveVal = cubicPolynomial.Evaluate(searchPos);

			if (curveVal < syncEndMinimumSampleValue)
			{
				searchPosLow = searchPos;
				continue;
			}
			else if (curveVal > syncEndMaximumSampleValue)
			{
				searchPosHigh = searchPos;
				continue;
			}
			break;
		}
	}

	// Return the position to the caller
	return searchPos;
}

#endif
