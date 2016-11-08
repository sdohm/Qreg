//
// logger.cpp
//
//  Created on: Aug 11,2014
//      Author: Sebastian Dohm <sebastian.dohm@uni-ulm.de>
//



/*
 * Copyright 2014-2016 Sebastian Dohm
 * 
 * This file is part of Qreg.
 *
 * Qreg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Qreg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Qreg.  If not, see <http://www.gnu.org/licenses/>.
 */



#include "logger.hpp"



namespace qreg
{



  Logger::Logger(const std::string& logFilePath):

    logFile_(logFilePath  ,
             std::ios::app )

  {


    //
    // DbC precondition
    //
    assert(logFile_.is_open());


  }



  Logger::~Logger(void)
  {


    //
    // DbC invariant
    //
    assert(logFile_.good());


    std::lock_guard<std::mutex> lock(threadSafety_);

    if(logFile_.is_open())
      logFile_.close();


  }



  auto
  Logger::message(const std::string& textHeader,
                  const std::string& textBody   )
  ->void
  {


    //
    // DbC invariant
    //
    assert(logFile_.good());


    //
    // auto
    // Logger::message(const std::string& textHeader,
    //                 const std::string& textBody   )
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(logFile_.is_open())
    {

      time_t        unixTime;
      struct tm* currentTime;

      unixTime    = time     (0        );
      currentTime = localtime(&unixTime);

      logFile_ << (currentTime->tm_year + 1900)
               << '-'
               << (currentTime->tm_mon  +    1)
               << '-'
               <<  currentTime->tm_mday
               << " "
               <<  currentTime->tm_hour
               << ":"
               <<  currentTime->tm_min
               << ":"
               <<  currentTime->tm_sec
               << "  ["
               << textHeader
               << "]"
               << std::endl
               << std::endl
               << textBody
               << std::endl
               << std::endl
               << std::endl                    ;

    }


    //
    // DbC invariant
    //
    assert(logFile_.good());


  }



}
