//
// logger.hpp
//
//  Created on: Nov 17,2014
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



#ifndef QREG_INTERFACE_LOGGER_
#define QREG_INTERFACE_LOGGER_



#include <cassert>
#include <fstream>
#include <mutex>
#include <string>



namespace qreg
{



  //
  // Provides functionality of logging to a text file.
  //
  class Logger final
  {



    public :


      //
      // Requires the path to a log file as an argument.
      //
      Logger(const std::string& logFilePath);

      //
      // Copy construction disabled.
      //
      Logger(const Logger& logger) = delete;

      //
      // Move construction disabled.
      //
      Logger(Logger&& logger) = delete;

      //
      // Copy assignment disabled.
      //
      auto
      operator=(const Logger& logger)
      ->Logger& = delete;

      //
      // Move assignment disabled.
      //
      auto
      operator=(Logger&& logger)
      ->Logger& = delete;

      //
      // Destructor.
      //
      ~Logger(void);


      //
      // Logs a text message with a header and a body. A time stamp is added.
      //
      auto
      message(const std::string& textHeader,
              const std::string& textBody   )
      ->void;


    private :


      //
      // Locking and Unlock the Logger for thread safety.
      //
      std::mutex threadSafety_;


      //
      // Stores an interface to the log file on disk.
      //
      std::ofstream logFile_;



  };



}



#endif // QREG_INTERFACE_LOGGER_
