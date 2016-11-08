//
// config.hpp
//
//  Created on: Nov 20, 2014
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



#ifndef QREG_INTERFACE_CONFIG_HPP_
#define QREG_INTERFACE_CONFIG_HPP_



#include "box.hpp"

#include <sstream>
#include <fstream>
#include <iostream>



namespace qreg
{



  //
  // The Config class handles parsing of the input file
  // and stores setup variables.
  //
  class Config final
  {



    public :


      //
      // Constructor.
      //
      Config(void);

      //
      // Copy construction disabled.
      //
      Config(const Config& config) = delete;

      //
      // Move construction enabled.
      //
      Config(Config&& config);

      //
      // Copy assignment disabled.
      //
      auto
      operator=(const Config& config)
      ->Config& = delete;

      //
      // Move assignment enabled.
      //
      auto
      operator=(Config&& config)
      ->Config&;

      //
      // Destructor.
      //
      ~Config(void);


      //
      // Parses the configuration file and stores gathered information.
      //
      auto
      initialize(const std::string& newConfigFilePath)
      ->void;


      //
      // Returns the path to the log file.
      //
      auto
      getLogFilePath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the MD binary.
      //
      auto
      getMdBinaryFilePath(void)const noexcept
      ->std::string;

      //
      // Returns the path to MD output folder.
      //
      auto
      getMdOutputFolderPath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the FF binary.
      //
      auto
      getFfBinaryFilePath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the FF output folder.
      //
      auto
      getFfOutputFolderPath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the QM binary.
      //
      auto
      getQmBinaryFilePath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the QM output folder.
      //
      auto
      getQmOutputFolderPath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the EVB binary.
      //
      auto
      getEvbBinaryFilePath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the EVB output folder.
      //
      auto
      getEvbOutputFolderPath(void)const noexcept
      ->std::string;

      //
      // Returns the path to the Qreg temporary files folder.
      //
      auto
      getTmpFolderPath(void)const noexcept
      ->std::string;

      //
      // Returns the MD header text.
      //
      auto
      getMdHeader(void)const noexcept
      ->std::string;

      //
      // Returns the FF header text.
      //
      auto
      getFfHeader(void)const noexcept
      ->std::string;

      //
      // Returns the QM header text.
      //
      auto
      getQmHeader(void)const noexcept
      ->std::string;

      //
      // Returns the EVB header text.
      //
      auto
      getEvbHeader(void)const noexcept
      ->std::string;


      //
      // Returns the geometry.
      //
      auto
      getGeometry(void)const noexcept
      ->Box;


      //
      // Returns N lines of M user defined commands
      // determining Qreg's work flow.
      //
      auto
      getNCommands(void)const noexcept
      ->std::vector< std::vector<std::string> >;


    private :


      //
      // Checks sanity i.e. Design by Contract compliance.
      //
      auto
      isSane(void)const noexcept
      ->bool;


      //
      // Locking and unlocking Config for thread safety.
      //
      std::mutex threadSafety_;


      //
      // Stores the path to the config file.
      //
      std::string configFilePath_;

      //
      // Stores the path to the log file.
      //
      std::string logFilePath_;

      //
      // Stores path to the MD binary.
      //
      std::string mdBinaryFilePath_;

      //
      // Stores the path to the MD output folder.
      //
      std::string mdOutputFolderPath_;

      //
      // Stores the path to the FF binary.
      //
      std::string ffBinaryFilePath_;

      //
      // Stores the path to the FF output folder.
      //
      std::string ffOutputFolderPath_;

      //
      // Stores the path to the QM binary.
      //
      std::string qmBinaryFilePath_;

      //
      // Stores the path to the QM output folder.
      //
      std::string qmOutputFolderPath_;

      //
      // Stores the path to the EVB binary.
      //
      std::string evbBinaryFilePath_;

      //
      // Stores the path to the EVB output folder.
      //
      std::string evbOutputFolderPath_;

      //
      // Stores the path to the Qreg temporary files folder.
      //
      std::string tmpFolderPath_;

      //
      // Stores the MD header text.
      //
      std::string mdHeader_;

      //
      // Stores the FF header text.
      //
      std::string ffHeader_;

      //
      // Stores the QM header text.
      //
      std::string qmHeader_;

      //
      // Stores the EVB header text.
      //
      std::string evbHeader_;


      //
      // Stores N lines of M user defined commands.
      // determining Qreg's work flow.
      //
      std::vector< std::vector<std::string> > nCommands_;


      //
      // Stores the geometry (if) defined in the configuration file.
      //
      Box geometry_;



  };



  inline auto
  Config::isSane(void)const noexcept
  ->bool
  {


    assert(!logFilePath_.empty());

    assert(!tmpFolderPath_.empty());

    assert(!nCommands_.empty());

    assert(!(mdBinaryFilePath_ .empty() &&
             ffBinaryFilePath_ .empty() &&
             qmBinaryFilePath_ .empty() &&
             evbBinaryFilePath_.empty()   ));

    assert(!(mdOutputFolderPath_ .empty() &&
             ffOutputFolderPath_ .empty() &&
             qmOutputFolderPath_ .empty() &&
             evbOutputFolderPath_.empty()   ));


    return true;


  }



}



#endif /* QREG_INTERFACE_CONFIG_HPP_ */

