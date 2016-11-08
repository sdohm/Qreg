//
// ifuff.hpp
//
//  Created on: Mar 10,2015
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



#ifndef QREG_INTERFACE_IFUFF_HPP_
#define QREG_INTERFACE_IFUFF_HPP



#include "ifbase.hpp"



namespace qreg
{



  //
  // Interface to the UFF program.
  //
  class IfUff final : public IfBase
  {



    public:


      //
      // Sets up interface to one UFF calculation.
      //
      IfUff(std::string newBinaryFilePath  ,
            std::string newHeaderString    ,
            std::string newOutputFolderPath )noexcept;

      //
      // Copy construction disabled.
      //
      IfUff(const IfUff& ifUff) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfUff& ifUff)
      ->IfUff& = delete;

      //
      // Move construction enabled.
      //
      IfUff(IfUff&& ifUff);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfUff&& ifUff)
      ->IfUff&;

      //
      // Destructor.
      //
      ~IfUff(void);


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfUff& lhs,
           IfUff& rhs )
      ->void;


      //
      // Harvests UFF calculation.
      //
      auto
      harvest(Box& box)
      ->bool override;

      //
      // Seeds UFF calculation.
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
    IfUff::getEnergyConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfUff::getEnergyConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree to E -23 J
      return 4.3597482e+5;


    }



    inline auto
    IfUff::getGradientConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfUff::getGradientConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree per Bohr to E -23 J per Angstrom
      return 8.2386866e+5;


    }



}



#endif // QREG_INTERFACE_IFUFF_HPP_
