//
// ifopengaff.hpp
//
//  Created on: Nov 19,2015
//      Author: Sebastian Dohm <sebastian.dohm@uni-ulm.de>
//



/*
 * Copyright 2015-2016 Sebastian Dohm
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



#ifndef QREG_INTERFACE_IFOPENGAFF_HPP_
#define QREG_INTERFACE_IFOPENGAFF_HPP



#include "ifbase.hpp"



namespace qreg
{



  //
  // Interface to OpenGAFF.
  //
  class IfOpenGaff final : public IfBase
  {



    public:


      //
      // Sets up interface to one OpenGAFF calculation.
      //
      IfOpenGaff(std::string newBinaryFilePath  ,
                 std::string newHeaderString    ,
                 std::string newOutputFolderPath )noexcept;

      //
      // Copy construction disabled.
      //
      IfOpenGaff(const IfOpenGaff& ifOpenGaff) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfOpenGaff& ifOpenGaff)
      ->IfOpenGaff& = delete;

      //
      // Move construction enabled.
      //
      IfOpenGaff(IfOpenGaff&& ifOpenGaff);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfOpenGaff&& ifOpenGaff)
      ->IfOpenGaff&;

      //
      // Destructor.
      //
      ~IfOpenGaff(void);


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfOpenGaff& lhs,
           IfOpenGaff& rhs )
      ->void;


      //
      // Harvests OpenGAFF calculation.
      //
      auto
      harvest(Box& box)
      ->bool override;

      //
      // Seeds OpenGAFF calculation.
      //
      auto
      seed(const Box&          box            ,
           const unsigned int& quantumRegionNo,
           const signed int&   charge         ,
           const unsigned int& multiplicity   ,
           const double&       range           )
      ->void override;


    private:


      //
      // Returns energy conversion factor.
      //
      auto
      getEnergyConversionFactor(void)const noexcept
      ->double override;

      //
      // Returns gradient conversion factor.
      //
      auto
      getGradientConversionFactor(void)const noexcept
      ->double override;

      //
      // Returns message for normal termination.
      //
      auto
      getNormalTerminationMessage(void)const noexcept
      ->std::string override;



  };



    inline auto
    IfOpenGaff::getEnergyConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOpenGaff::getEnergyConversionFactor(void)const noexcept
      // ->double
      //

      // no conversion
      return 1;


    }



    inline auto
    IfOpenGaff::getGradientConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOpenGaff::getGradientConversionFactor(void)const noexcept
      // ->double
      //

      // no conversion
      return 1;


    }



}



#endif // QREG_INTERFACE_IFOPENGAFF_HPP_
