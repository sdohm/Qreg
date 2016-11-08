//
// ifopenuff.hpp
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



#ifndef QREG_INTERFACE_IFOPENUFF_HPP_
#define QREG_INTERFACE_IFOPENUFF_HPP



#include "ifbase.hpp"



namespace qreg
{



  //
  // Interface to OpenUFF.
  //
  class IfOpenUff final : public IfBase
  {



    public:


      //
      // Sets up interface to one OpenUFF calculation.
      //
      IfOpenUff(std::string newBinaryFilePath  ,
                 std::string newHeaderString    ,
                 std::string newOutputFolderPath )noexcept;

      //
      // Copy construction disabled.
      //
      IfOpenUff(const IfOpenUff& ifOpenUff) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfOpenUff& ifOpenUff)
      ->IfOpenUff& = delete;

      //
      // Move construction enabled.
      //
      IfOpenUff(IfOpenUff&& ifOpenUff);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfOpenUff&& ifOpenUff)
      ->IfOpenUff&;

      //
      // Destructor.
      //
      ~IfOpenUff(void);


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfOpenUff& lhs,
           IfOpenUff& rhs )
      ->void;


      //
      // Harvests OpenUFF calculation.
      //
      auto
      harvest(Box& box)
      ->bool override;

      //
      // Seeds OpenUFF calculation.
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
    IfOpenUff::getEnergyConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOpenUff::getEnergyConversionFactor(void)const noexcept
      // ->double
      //

      // no conversion
      return 1;


    }



    inline auto
    IfOpenUff::getGradientConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOpenUff::getGradientConversionFactor(void)const noexcept
      // ->double
      //

      // no conversion
      return 1;


    }



}



#endif // QREG_INTERFACE_IFOPENUFF_HPP_
