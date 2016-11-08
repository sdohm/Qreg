//
// ifevb.hpp
//
//  Created on: Jul 1,2015
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



#ifndef QREG_INTERFACE_IFEVB_HPP_
#define QREG_INTERFACE_IFEVB_HPP_



#include "ifbase.hpp"



namespace qreg
{



  //
  // Interface to the EVB program.
  //
  class IfEvb final : public IfBase
  {



    public:


      //
      // Sets up interface to one EVB calculation.
      //
      IfEvb(std::string newBinaryFilePath  ,
            std::string newHeaderString    ,
            std::string newOutputFolderPath )noexcept;

      //
      // Copy construction disabled.
      //
      IfEvb(const IfEvb& ifEvb) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfEvb& ifEvb)
      ->IfEvb& = delete;

      //
      // Move construction enabled.
      //
      IfEvb(IfEvb&& ifEvb);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfEvb&& ifEvb)
      ->IfEvb&;

      //
      // Destructor.
      //
      ~IfEvb(void);


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfEvb& lhs,
           IfEvb& rhs )
      ->void;


      //
      // Harvests EVB calculation.
      //
      auto
      harvest(Box& box)
      ->bool override;

      //
      // Seeds EVB calculation.
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
    IfEvb::getEnergyConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfEvb::getEnergyConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree to E -23 J
      return 4.359810e+5;


    }



    inline auto
    IfEvb::getGradientConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfEvb::getGradientConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree to E -23 J
      return 4.359810e+5;


    }



}



#endif // QREG_INTERFACE_IFEVB_HPP_
