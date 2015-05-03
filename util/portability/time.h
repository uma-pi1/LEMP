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


#ifndef RG_TIME_H
#define RG_TIME_H

#include <util/config.h>

// in windows: define timeval and gettimeofday() 
// otherwise: use sys/time.h
#if defined(_MSC_VER) || defined(_WINDOWS_)
	#include <windows.h>
	#include <time.h>
	
	#if !defined(_WINSOCK2API_) && !defined(_WINSOCKAPI_)
		struct timeval {
			long tv_sec;
            long tv_usec;
		};
   	#endif 
	
	inline int gettimeofday(struct timeval* tv, void* tzp) {
		union {
	         long long ns100;
	         FILETIME ft;
		} now;
		     
		GetSystemTimeAsFileTime (&now.ft);
		tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
		tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);

		return 0;
	}
#else
   #include <sys/time.h>
#endif 

#endif 
