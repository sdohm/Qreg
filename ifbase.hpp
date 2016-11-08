//
// ifbase.hpp
//
//  Created on: Nov 26,2014
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



#ifndef QREG_INTERFACE_IFBASE_HPP_
#define QREG_INTERFACE_IFBASE_HPP_



#include "box.hpp"
#include "idconverter.hpp"

// POSIX / C stuff
#include <cstdlib>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>



namespace qreg
{



  //
  // Traits ABC for interfaces with QMMM devices.
  //
  class IfBase
  {



    public:


      //
      // Pure virtual ABC.
      //
      IfBase(std::string newBinaryFilePath   = "",
             std::string newHeaderText       = "",
             std::string newOutputFolderPath = "");

      //
      // Copy construction disabled.
      //
      IfBase(const IfBase& ifBase) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfBase& ifBase)
      ->IfBase& = delete;

      //
      // Move construction enabled.
      //
      IfBase(IfBase&& ifBase);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfBase&& ifBase)
      ->IfBase&;


      //
      // Synthesized virtual destructor.
      //
      virtual
      ~IfBase(void) = default;


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfBase& lhs,
           IfBase& rhs )
      ->void;


      //
      // Virtual prototype for harvesting a calculation.
      //
      virtual auto
      harvest(Box& box)
      ->bool = 0;

      //
      // Virtual prototype for seeding a calculation.
      //
      virtual auto
      seed(const Box&          box            ,
           const unsigned int& quantumRegionNo,
           const signed int&   charge         ,
           const unsigned int& multiplicity   ,
           const double&       range           )
      ->void = 0;


      //
      // Checks if a calculation has been harvested already.
      //
      auto
      hasBeenHarvested(void)const noexcept
      ->bool;

      //
      // Checks if a calculation has finished.
      //
      virtual auto
      hasFinished(void)const noexcept
      ->bool;

      //
      // Checks whether this works on a specified quantum region.
      //
      auto
      hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
      ->bool;

      //
      // Checks if a calculation has terminated normally.
      //
      auto
      hasTerminatedNormally(void)const
      ->bool;


    protected:


      //
      // Returns energy conversion factor.
      //
      virtual auto
      getEnergyConversionFactor(void)const noexcept
      ->double = 0;

      //
      // Returns gradient conversion factor.
      //
      virtual auto
      getGradientConversionFactor(void)const noexcept
      ->double = 0;

      //
      // Returns message for normal termination.
      //
      virtual auto
      getNormalTerminationMessage(void)const noexcept
      ->std::string = 0;


      //
      // Locking and Unlock the Interface for thread safety.
      //
      std::mutex threadSafety_;


      //
      // Process ID.
      //
      pid_t processId_;


      //
      // Stores processing status.
      //
      bool hasBeenHarvested_;


      //
      // Stores unique quantum region.
      //
      unsigned int quantumRegionNo_;


      //
      // Stores path to binary.
      //
      std::string binaryFilePath_;

      //
      // Stores header text.
      //
      std::string headerText_;

      //
      // Stores input text.
      //
      std::string inputText_;

      //
      // Stores path to output file.
      //
      std::string outputFilePath_;

      //
      // Stores output file path.
      //
      std::string outputFolderPath_;

      //
      // Stores a timestamp.
      //
      std::string timeStamp_;


      //
      // Stores an ID conversion tool.
      //
      IdConverter<size_t,
                  size_t > idConverter_;



  };



}



#endif // qreg_iNTERFACE_IFBASE_HPP_
