//
// ifbase.cpp
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



#include "ifbase.hpp"



namespace qreg
{



  IfBase::IfBase(std::string newBinaryFilePath  ,
                 std::string newHeaderText      ,
                 std::string newOutputFolderPath ):

    processId_       (0                  ),
    hasBeenHarvested_(false              ),
    quantumRegionNo_ (0                  ),
    binaryFilePath_  (newBinaryFilePath  ),
    headerText_      (newHeaderText      ),
    inputText_       (                   ),
    outputFilePath_  (                   ),
    outputFolderPath_(newOutputFolderPath),
    timeStamp_       (                   ),
    idConverter_     (                   )

  {


    //
    // DbC PRE
    //

    assert(!binaryFilePath_    .empty());
    assert(!newOutputFolderPath.empty());


  }



  IfBase::IfBase(IfBase&& ifBase)
  {


    // Call move assignment operator.
    *this = std::move(ifBase);


  }



  auto
  IfBase::operator=(IfBase&& ifBase)
  ->IfBase&
  {



    //
    // Perform move operation.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    std::swap(       processId_,
              ifBase.processId_ );

    std::swap(       hasBeenHarvested_,
              ifBase.hasBeenHarvested_ );

    std::swap(       quantumRegionNo_,
              ifBase.quantumRegionNo_ );

    std::swap(       binaryFilePath_,
              ifBase.binaryFilePath_ );

    std::swap(       headerText_,
              ifBase.headerText_ );

    std::swap(       inputText_,
              ifBase.inputText_ );

    std::swap(       outputFilePath_,
              ifBase.outputFilePath_ );

    std::swap(       outputFolderPath_,
              ifBase.outputFolderPath_ );

    std::swap(       timeStamp_,
              ifBase.timeStamp_ );

    idConverter_ = ifBase.idConverter_;


    ifBase.~IfBase();


    return *this;



  }


  auto
  IfBase::hasBeenHarvested(void)const noexcept
  ->bool
  {


    //
    // auto
    // IfBase::hasBeenHarvested(void)const noexcept
    // ->bool
    //

    return hasBeenHarvested_;


  }



  auto
  IfBase::hasFinished(void)const noexcept
  ->bool
  {


    //
    // auto
    // IfBase::hasFinished(void)const noexcept
    // ->bool
    //

    return (kill(processId_,
                 0          ) != 0);


  }



  auto
  IfBase::hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
  ->bool
  {


    //
    // DbC PRE
    //
    assert(quantumRegionNo != 0);


    //
    // auto
    // IfBase::hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
    // ->bool
    //

    return (quantumRegionNo_ == quantumRegionNo);


  }



  auto
  IfBase::hasTerminatedNormally(void)const
  ->bool
  {


    //
    // auto
    // IfBase::hasTerminatedNormally(void)const
    // ->bool
    //

    if(!hasFinished())
      return false;

    else
    {
      std::string   line      ;
      std::ifstream outputFile;
      outputFile.open(outputFilePath_.c_str());
      if(outputFile.is_open())
      {
        while(std::getline(outputFile,
                           line       ))
        {
          if(line == getNormalTerminationMessage())
          {
            outputFile.close();
            return true;
          }
        }
        outputFile.close();
      }
    }


    return false;


  }



}
