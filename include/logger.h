/**
 *
 * @file logger.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 *
 * @author Mustafa Abduljabbar [mustafa.abduljabbar@kaust.edu.sa] and Mohammed Al Farhan [mohammed.farhan@kaust.edu.sa].
 *
 **/

#ifndef logger_h
#define logger_h
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <pthread.h>
#include <queue>
#include <stdint.h>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>

namespace bemfmm {
  //! Structure for pthread based tracer
  struct Tracer {
    pthread_t thread;                                           //!< pthread id
    double    begin;                                            //!< Begin timer of tracer
    double    end;                                              //!< End timer of tracer
    Tracer() {}                                                 //!< Constructor
  };

  //! Timer and Tracer logger
  namespace alogger {
    typedef std::map<std::string,double> Timer;                 //!< Map of timer event name to timed value
    typedef Timer::iterator              T_iter;                //!< Iterator of timer event name map
    typedef std::queue<Tracer>           Tracers;               //!< Queue of tracers
    typedef std::map<pthread_t,int>      ThreadMap;             //!< Map of pthread id to thread id

    Timer           beginTimer;                                 //!< Timer base value
    Timer           timer;                                      //!< Timings of all events
    Tracers         tracers;                                    //!< Tracers for all events
    pthread_mutex_t mutex;                                      //!< Pthread communicator

    int stringLength = 40;                                      //!< Max length of event name
    int decimal = 7;                                            //!< Decimal precision
    int precis = 17;                                            //!< Decimal precision for numerical results
    bool verbose = false;                                       //!< Print to screen
    const char * path = "./";                                   //!< Path to save files

    //! Timer function
    double get_time() {
      struct timeval tv;                                        // Time value
      gettimeofday(&tv, NULL);                                  // Get time of day in seconds and microseconds
      return double(tv.tv_sec)+double(tv.tv_usec)*1e-6;         // Combine seconds and microseconds and return
    }

    //! Cycle counter
    inline uint64_t get_cycle() {
      uint32_t low = 0, high = 0;                               // Define low and high 32 bits of cycle counter
#if !(__FUJITSU | _SX)
      asm volatile ("rdtsc" : "=a" (low), "=d" (high));         // Call rdtsc
#endif
      return (uint64_t(high) << 32) | uint64_t(low);            // Return 64 bit cycle counter
    }

    //! Cycle counter with thread ID
    inline uint64_t get_cycle(uint32_t * id) {
      uint32_t low = 0, high = 0;                               // Define low and high 32 bits of cycle counter
      if (!id) return 0;                                        // Count only for valid thread ID
#if !(__FUJITSU | _SX)
      asm volatile ("rdtscp" : "=a" (low), "=d" (high), "=c" (*id));// Call rdtscp
#endif
      return (uint64_t(high) << 32) | uint64_t(low);            // Return 64 bit cycle counter
    }

    //! Print message to standard output
    inline void printTitle(std::string title) {
      if (verbose) {                                            // If verbose flag is true
	title += " ";                                           //  Append space to end of title
	std::cout << "--- " << std::setw(stringLength)          //  Align string length
		  << std::left                                  //  Left shift
		  << std::setfill('-')                          //  Set to fill with '-'
		  << title << std::setw(10) << "-"              //  Fill until end of line
		  << std::setfill(' ') << std::endl;            //  Set back to fill with ' '
      }                                                         // End if for verbose flag
    }

    //! Start timer for given event
    inline void startTimer(std::string event) {
      beginTimer[event] = get_time();                           // Get time of day and store in beginTimer
    }

    //! Print timings of a specific event
    inline void printTime(std::string event, int force=0) {
      if (verbose||force) {                               // If verbose flag is true
	std::cout << std::setw(stringLength) << std::left       //  Set format
		  << event << " : " << std::setprecision(decimal) << std::fixed
		  << timer[event] << " s" << std::endl;         //  Print event and timer
      }                                                         // End if for verbose flag
    }
    
    template<typename ValT>
    void logItem(int stringLength, std::string item, ValT value) {
      std::cout << std::setw(stringLength) << std::fixed << std::left
      << item << " : " << value << std::endl;
    }

    void checkMatchingFiles(std::string file1, std::string file2) { 
      std::ifstream outfile1 (file1.c_str(),std::ios_base::binary);
      std::ifstream outfile2 (file2.c_str(),std::ios_base::binary);
      double diff  = 0.0;
      double norm1 = 0.0;
      if(outfile1.is_open() && outfile2.is_open()) {
        double num1, num2;
        uint32_t numsRead = 0;
        outfile1>>num1;
        outfile2>>num2;
        while(!outfile1.eof() && !outfile2.eof()) {     
          diff += (num1 - num2) * (num1 - num2);
          norm1 += num2*num2;
          outfile1>>num1;
          outfile2>>num2;
          numsRead++;
        }
        numsRead--;
        if(!outfile1.eof() || !outfile2.eof()) {
          std::cerr<< file1 << " and " << file2 <<" are not identical. Read " << numsRead << " numbers. File ";
          if(!outfile1.eof()) {
            std::cerr<<file1 << " ";
          }
          if(!outfile2.eof()) {
            std::cerr<<file2 << " ";
          }
          std::cerr<<" is not finished"<< std::endl;
          }
      } else {
         std::cerr<<"either " << file1 << " or " << file2 <<" cannot be openned"<<std::endl;
      }
      std::cout << std::setw(stringLength) << std::fixed <<std::left
                << file1 + "(err)" <<" : " << std::scientific << std::sqrt(diff/norm1) << std::endl;
      outfile1.close();
      outfile2.close();
    }

    template<typename Arr1, typename Arr2>
    void logNumericalError(Arr1 const& arr1, Arr2 const& arr2, size_t size, std::string const experiemnt) {
      double diff  = 0.0;
      double norm1 = 0.0;
      for (size_t i = 0; i < size; ++i) {
        diff += (std::abs(arr1[i]) - std::abs(arr2[i])) * (std::abs(arr1[i]) - std::abs(arr2[i]));
        norm1 += std::abs(arr2[i])*std::abs(arr2[i]);
      }
      if(verbose)
      	std::cout << std::setw(stringLength) << std::fixed << std::left
                << experiemnt + "(err)" <<" : " << std::scientific << std::sqrt(diff/norm1) << std::endl;

    }

    template<typename Arr1, typename Arr2, typename Arr3>
    void logNumericalError(Arr1 const& arr1, Arr2 const& arr2, Arr3 const& add, size_t size, std::string const experiemnt) {
      double diff  = 0.0;
      double norm1 = 0.0;
      for (size_t i = 0; i < size; ++i) {
        diff += (std::abs(arr1[add[i]]) - std::abs(arr2[add[i]])) * (std::abs(arr1[add[i]]) - std::abs(arr2[add[i]]));
        norm1 += std::abs(arr2[add[i]])*std::abs(arr2[add[i]]);
      }
      if(verbose)
      	std::cout << std::setw(stringLength) << std::fixed << std::left
                << experiemnt + "(err)" <<" : " << std::scientific << std::sqrt(diff/norm1) << std::endl;

    }
    //! Stop timer for given event
    double stopTimer(std::string event, int print=1, int force_print=0) {
      double endTimer = get_time();                             // Get time of day and store in endTimer
      timer[event] += endTimer - beginTimer[event];             // Accumulate event time to timer
      if ((verbose && print)||force_print) printTime(event, force_print);      // Print event and timer to screen
      return endTimer - beginTimer[event];                      // Return the event time
    }

    //! Stop and reset timer for given event
    double stopResetTimer(std::string event, int print=1, int force_print=0) {
      double endTimer = get_time();                             // Get time of day and store in endTimer
      timer[event] += endTimer - beginTimer[event];             // Accumulate event time to timer
      if ((verbose && print)||force_print) printTime(event, force_print);    // Print event and timer to screen
      timer.erase(event);                                       // Reset the timer event
      return endTimer - beginTimer[event];                      // Return the event time
    }

    //! Write timings of all events
    inline void writeTime(int mpirank=0) {
      std::stringstream name;                                   // File name
      name << path << "time" << std::setfill('0') << std::setw(6) // Set format
	   << mpirank << ".dat";                                // Create file name for timer
      std::ofstream timerFile(name.str().c_str());              // Open timer log file
      for (T_iter E=timer.begin(); E!=timer.end(); E++) {       // Loop over all events
	timerFile << std::setw(stringLength) << std::left       //  Set format
		  << E->first << " " << E->second << std::endl; //  Print event and timer
      }                                                         // End loop over all events
      timerFile.close();                                        // Close timer log file
      timer.clear();                                            // Clear timer
    }

    //! Erase single event in timer
    inline void resetTimer(std::string event) {
      timer.erase(event);                                       // Erase event from timer
    }

    //! Erase all events in timer
    inline void resetTimer() {
      timer.clear();                                            // Clear timer
    }

  };
}

#endif
