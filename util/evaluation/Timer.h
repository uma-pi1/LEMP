
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

#ifndef RG_TIMER_H
#define RG_TIMER_H

#include <util/config.h>

#include <util/portability/inttypes.h>
#include <util/portability/time.h>
#include <util/evaluation/duration.h>
#include <util/evaluation/formatting.h>
#include <util/evaluation/timeval.h>

#include <iostream>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- Timer class        
// ------------------------------------------------------------------------------------------------	

/** Utility class to measure wall-clock time. Everything is inlined to  minimize the overhead of 
 * time measurement. */
class Timer {
public:
	/** Standard constructor. The timer is not running after construction. */ 
	Timer();
	
	/** Resets the timer and starts it */
	void start();

	/** Stops the timer. The call to this method is ignored when the timer is not running. */
	void stop();

	/** Resumes the timer. The call to this method is ignored when the timer
	 * is already running. */
	void resume();

	/** Stops and resets the timer to 0. */
	void reset();
	
	/** Returns the elapsed time between all (non-ignored) start() and stop() calls. If the timer 
	 * is running, the time elapsed since the last call to start() is not included into the result. 
	 * 
	 * @return the elapsed time between all (non-ignored) start() and stop() calls
	 */
	NanoDuration<uint64_t> elapsedTime() const;
	
	/** Subtracts two timers. This method assumes that t1>t2. The returned timer is in 
	 * stopped state. */
	Timer operator-(const Timer& other);
	
private:
	/** Time since last call to start() or zeroed out if not running */
	timeval m_startTime;
	
	/** Total elapsed time */
	NanoDuration<uint64_t> m_elapsedTime;
};


// ------------------------------------------------------------------------------------------------
// -- Timer inline members      
// ------------------------------------------------------------------------------------------------

inline Timer::Timer() { 
	reset(); 
}

inline void Timer::resume() {
	// update startTime only when timer is not already running
	if (m_startTime.tv_sec == 0 && m_startTime.tv_usec == 0) {
		gettimeofday(&m_startTime, 0);
	}		
} 

inline void Timer::stop() {
	// get endTime quickly, then do checks
	timeval endTime;
	gettimeofday(&endTime, 0);
	
	// ignore call when timer is not running
	if (m_startTime.tv_sec == 0 && m_startTime.tv_usec == 0) {
		return;
	}
	
	// add to elapsed time
	m_elapsedTime += endTime - m_startTime;
	
	// clear starttime
	m_startTime.tv_sec = 0;
	m_startTime.tv_usec = 0;
}		

inline void Timer::reset() {
	m_startTime.tv_sec = 0;
	m_startTime.tv_usec = 0;
	m_elapsedTime = 0;
}

inline void Timer::start() {
	reset();
	resume();
}

inline NanoDuration<uint64_t> Timer::elapsedTime() const { 
	return m_elapsedTime;
}

inline Timer Timer::operator-(const Timer& other) {
	Timer t;
	t.m_elapsedTime = m_elapsedTime - other.m_elapsedTime;
	return t;
}


// ------------------------------------------------------------------------------------------------
// -- Timer output         
// ------------------------------------------------------------------------------------------------	

/** Formatted output of a Timer */
template<typename CharT, typename Traits> 
std::basic_ostream<CharT, Traits>& operator<<(
		std::basic_ostream<CharT, Traits>& out, const Timer& t) 
{
	return out << t.elapsedTime();
};

		
} // namespace rg

#endif
