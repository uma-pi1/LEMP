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


#ifndef RG_IO_MANIPULATOR_H
#define RG_IO_MANIPULATOR_H

#include <util/config.h>

#include <util/exception.h>

#include <iostream>
#include <map>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- public function declarations         
// ------------------------------------------------------------------------------------------------	

/** Allocates a new manipulator and returns a handle. The handle is equal to the index in 
 * the internal I/O stream array reserverd for the new manipulator. 
 * 
 * @tparam T type of data
 * @param defaultValue default value of the manipulator
 * @return handle for manipulator
 */ 
template<class T>
int allocManipulator(const T& defaultValue = T());

/** Updates the data associated with the specified manipulator in the specified stream.
 *
 * @tparam T type of data
 * @param stream the stream
 * @param handle handle of manipulator
 * @param data new data for the manipulator (will be copied into free memory)
 */ 
template<class T, class CharT, class Traits>
void setManipulator(std::basic_ios<CharT, Traits>& stream, int handle, const T& data);

/** Returns the data associated with the specified manipulator type in the specified stream.
 * The manipulator data must have been previously set with setManipulator.
 *
 * @tparam T type of data
 * @param stream the stream
 * @param handle handle of manipulator 
 * @return the data associated with the specified manipulator type
 * 
 * @throw IllegalStateException when manipulator has not been set
 */
template<class T, class CharT, class Traits>
const T& getManipulator(std::basic_ios<CharT, Traits>& stream, int handle);


// ------------------------------------------------------------------------------------------------
// -- Manipulator class        
// ------------------------------------------------------------------------------------------------	

/** Helper class that stores the data for stream manipulators */
template<class T>
class Manipulator {
	/** handle of manipulator */
	const int m_handle;
	
	/** data associated with it */
	const T& m_data;

public:
	explicit Manipulator(int handle, const T& data) : m_handle(handle), m_data(data) { };
  
	/** Update the manipulator and return the I/O stream. */ 
	template<class CharT, class Traits>
	void set(std::basic_ios<CharT, Traits>& io) const {
		setManipulator(io, m_handle, m_data);
	}
};


// ------------------------------------------------------------------------------------------------
// -- I/O wrappers
// ------------------------------------------------------------------------------------------------	

/** Use Manipulator objects to set manipulators */
template<class T, class CharT, class Traits>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& o, const Manipulator<T>& m) {
	m.set(o);
	return o;
};

/** Use Manipulator objects to set manipulators */
template<class T, class CharT, class Traits>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& i, const Manipulator<T>& m) {
	m.set(i);
	return i;
};


// ------------------------------------------------------------------------------------------------
// -- internal classes & functions         
// ------------------------------------------------------------------------------------------------	

namespace io { // hide 

/** Class for holding default manipulator values */
template<class T>
class DefaultManipulatorValues {
private:
	/** maps manipulator handles to default values */ 
	static typename std::map<int, T> defaultValues;
	
public:
	static const T& get(int handle) {
		typename std::map<int, T>::iterator it = defaultValues.find(handle);
		if (it == defaultValues.end()) RG_THROW(IllegalStateException, "default value not found");
		return it->second;
	}
	
	static void set(int handle, const T& defaultValue) {
		defaultValues[handle] = defaultValue; // copy
	}

	static void erase(int handle) {
		defaultValues.erase(handle);
	}
};

template<class T>
typename std::map<int, T> DefaultManipulatorValues<T>::defaultValues = std::map<int,T>();

/** Callback-method used to delete/copy manipulator strings when the stream is closed/copied.
 * 
 * @tparam T type of data
 * @param event type of event
 * @param stream the stream
 * @param index the index associated with the event 
 */
template<class T>
void manageManipulator(std::ios_base::event event, std::ios_base& stream, int handle) {
	T* p = static_cast<T*>(stream.pword(handle));
	
	switch (event) {
	case std::ios_base::erase_event: // delete data
		if (p) { 
			delete p;
		} 
		break;
	case std::ios_base::copyfmt_event:
		RG_THROW(NotImplementedException, "copying not implemented"); 
		break;
	default: 
		break; // imbue_event does not affect storage
	}
};

} // namespace io

// ------------------------------------------------------------------------------------------------
// -- function definitions (inlines and templates)
// ------------------------------------------------------------------------------------------------

template<class T>
int allocManipulator(const T& defaultValue = T()) {
	int handle = std::ios::xalloc();
	io::DefaultManipulatorValues<T>::set(handle, defaultValue);
	return handle;
}

template<class T, class CharT, class Traits>
void setManipulator(std::basic_ios<CharT, Traits>& stream, int handle, const T& data) {
	// get the current pointer 
	T* p = static_cast<T*>(stream.pword(handle));
	if (p) {
		// non-empty --> remove it from heap  
		delete p;
	} else { 
		// register a call-back
		stream.register_callback(io::manageManipulator<T>, handle);
	}
	
	// set the new pointer (copy)
	stream.pword(handle) = new T(data);
};

template<class T, class CharT, class Traits>
const T& getManipulator(std::basic_ios<CharT, Traits>& stream, int handle) {
	if (T* p = static_cast<T*>(stream.pword(handle))) return *p; 
	return io::DefaultManipulatorValues<T>::get(handle);
};

} // namespace rg

#endif
