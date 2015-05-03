
//    Copyright 2015 Rainer Gemulla
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.


// formatting.h
//
// formatting of time-related values

#ifndef RG_EVALUATION_FORMATTING_H
#define RG_EVALUATION_FORMATTING_H

#include <util/io.h>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- formatting information  
// ------------------------------------------------------------------------------------------------	

/** Holds enumerations used to determine output format of time-related values */
struct TimeFormat {
	/** Output formats for durations (e.g., 1ms) */
	enum Duration { 
		nanos,
		millis,
		micros,
		seconds
	};
	
	/** Output formats for frequencies (e.g., 1k) */
	enum Frequency { 
		plain,
		kilo,
		mega,
		giga
	};
};

// ------------------------------------------------------------------------------------------------
// -- manuipulators  
// ------------------------------------------------------------------------------------------------	

namespace evaluation { // hide

/** Describes formatting of rates (e.g., 1M/s) */
struct RateFormatDescriptor {
	TimeFormat::Frequency countScale;
	TimeFormat::Duration timeScale;
	
	RateFormatDescriptor() 
		: countScale(TimeFormat::plain), timeScale(TimeFormat::seconds) {};
		
	RateFormatDescriptor(TimeFormat::Frequency countScale_, TimeFormat::Duration timeScale_)
		: countScale(countScale_), timeScale(timeScale_) {};
		
	RateFormatDescriptor& operator=(RateFormatDescriptor other) {
		countScale = other.countScale;
		timeScale = other.timeScale;
		return *this;
	}
};

/** Manipulators for time-related values */
enum TimeManipulator { 
	duration,
	rate
};
	
static const int NO_MANIPULATORS = rate + 1;

/** Determines the index in the internal extensible array that corresponds to the given manipulator. 
 *
 * @param m type of manipulator
 * @return index of manipulator in the internal extensible array  
 */ 
inline int getHandle(TimeManipulator m) {
	static const int handles[evaluation::NO_MANIPULATORS] = {
		allocManipulator(TimeFormat::seconds), // duration
		allocManipulator(RateFormatDescriptor()) // rate
	};
	
	return handles[m];
};

} // namespace evaluation


// ------------------------------------------------------------------------------------------------
// -- utility functions for setting manipulators
// ------------------------------------------------------------------------------------------------	

/** Sets the style of formatting of durations (such as 1ms) */
inline Manipulator<TimeFormat::Duration> set_duration(const TimeFormat::Duration& d) {
   return Manipulator<TimeFormat::Duration>(evaluation::getHandle(evaluation::duration), d);
};

/** Sets the style of formatting of rates (such as 1M/s = millions per second) */
inline Manipulator<evaluation::RateFormatDescriptor> set_rate(
		const TimeFormat::Frequency& count, const TimeFormat::Duration per_time) {
	using evaluation::RateFormatDescriptor;
	
	return Manipulator<RateFormatDescriptor>(
		   evaluation::getHandle(evaluation::rate), 
		   RateFormatDescriptor(count, per_time));
};


// ------------------------------------------------------------------------------------------------
// -- methods for stream output  
// ------------------------------------------------------------------------------------------------	

/** Outputs a duration. */ 
template<typename Duration, typename CharT, typename Traits>
void outputDuration(std::basic_ostream<CharT, Traits>& out, const Duration& d) {
	using evaluation::getHandle;
	using evaluation::duration;
	
	switch (getManipulator<TimeFormat::Duration>(out, getHandle(duration))) {
	case TimeFormat::nanos:
		out << d.nanos() << "ns";
		break;
	case TimeFormat::micros:
		out << d.micros() << "us";
		break;
	case TimeFormat::millis:
		out << d.millis() << "ms";
		break;
	case TimeFormat::seconds: 
		out << d.seconds() << "s";
		break;
	default:
		out << "ILLEGAL DURATION FORMAT";
	}
}
	
/** Outputs a rate. */
template<typename Duration, typename CharT, typename Traits>
void outputRate(std::basic_ostream<CharT, Traits>& out, const Duration& d) {
	using evaluation::getHandle;
	using evaluation::rate;
	using evaluation::RateFormatDescriptor;
	
	// determine unit
	double multiplier = 1.;	
	string unit;
	RateFormatDescriptor rd = getManipulator<RateFormatDescriptor>(out, getHandle(rate)); 

	// count unit
	switch (rd.countScale) {
	case TimeFormat::giga:
		multiplier /= 1000000000.;
		unit = "G";
		break;
	case TimeFormat::mega:
		multiplier /= 1000000.;
		unit = "M";
		break;
	case TimeFormat::kilo:
		multiplier /= 1000.;
		unit = "k";
		break;
	case TimeFormat::plain: 
		break;
	default:
		out << "ILLEGAL FREQUENCY FORMAT";
		return;
	}
	
	// time unit
	switch (rd.timeScale) {
	case TimeFormat::nanos:
		multiplier *= 1000000000.;
		unit = "/ns";
		break;
	case TimeFormat::micros:
		multiplier *= 1000000.;
		unit = "/us";
		break;
	case TimeFormat::millis:
		multiplier *= 1000.;
		unit = "/ms";
		break;
	case TimeFormat::seconds: 
		unit += "/s";
		break;
	default:
		out << "ILLEGAL RATE FORMAT";
		return;
	}
	
	out << d.perSecond()*multiplier << unit;
}

} // namespace rg

#endif
