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


#ifndef RG_DURATION_H
#define RG_DURATION_H

#include <util/config.h>
#include <util/evaluation/formatting.h>
#include <util/portability/time.h>

#include <iostream>
#include <iomanip>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- Duration base class
// ------------------------------------------------------------------------------------------------	

/** This class represents a time duration. */
class Duration {
public:
	virtual ~Duration() { };
	
	// retrieve duration
	virtual double nanos() const = 0;
	double micros() const { return nanos()/1000.; }
	double millis() const { return micros()/1000.; }
	double seconds() const { return millis()/1000.; }
	double minutes() const { return seconds()/60.; }
	double hours() const { return minutes()/60.; }
	double days() const { return hours()/24.; }	
	
	// retrieve throughput
	double perNano() const { return 1./nanos(); }
	double perMicro() const { return 1./micros(); }
	double perMilli() const { return 1./millis(); }
	double perSecond() const { return 1./seconds(); }
	double perMinute() const { return 1./minutes(); }
	double perHour() const { return 1./hours(); }
	double perDay() const { return 1./days(); }
	
	// comparisons
	bool operator<(const Duration& o) const { return nanos() < o.nanos(); }
	bool operator<=(const Duration& o) const { return nanos() <= o.nanos(); }
	bool operator==(const Duration& o) const { return nanos() == o.nanos(); }
	bool operator>=(const Duration& o) const { return nanos() >= o.nanos(); }
	bool operator>(const Duration& o) const { return nanos() > o.nanos(); }
	bool operator!=(const Duration& o) const { return nanos() != o.nanos(); }
};


// ------------------------------------------------------------------------------------------------
// -- NanosDuration class
// ------------------------------------------------------------------------------------------------	

/** A time duration measured in nanoseconds.
 * 
 * @tparam NanoType type used to store nanosecond counter
 */
template<typename NanoType>
class NanoDuration : public Duration {
	NanoType m_nanos;
	
public:
	NanoDuration() : m_nanos(0) { };	
	NanoDuration(NanoType nanos) : m_nanos(nanos) { };
	
	double nanos() const { return m_nanos; }
	
	NanoDuration<NanoType>& operator=(NanoType nanos) { 
		m_nanos = nanos; 
		return *this;
	}	
	NanoDuration<NanoType>& operator+=(const Duration& o) { 
		m_nanos += (NanoType)o.nanos(); 
		return *this; 
	}
	NanoDuration<NanoType>& operator-=(const Duration& o) { 
		m_nanos -= (NanoType)o.nanos(); 
		return *this; 
	}
};

template<typename NanoType>
inline NanoDuration<NanoType> operator+(const NanoDuration<NanoType>& o1, const Duration& o2) { 
	return NanoDuration<NanoType>(o1) += o2; 
}

template<typename NanoType>
inline NanoDuration<NanoType> operator-(const NanoDuration<NanoType>& o1, const Duration& o2) { 
	return NanoDuration<NanoType>(o1) -= o2; 
}

// ------------------------------------------------------------------------------------------------
// -- Duration output
// ------------------------------------------------------------------------------------------------	

/** Formatted output of a duration */
template<typename Ch, typename Tr> 
std::basic_ostream<Ch,Tr>& operator<<(std::basic_ostream<Ch,Tr>& out, const Duration& t) {
	outputDuration(out, t);
	return out;   
}

} // namespace rg


#endif
