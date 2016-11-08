//
// iforc.hpp
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



#ifndef QREG_INTERFACE_IFORC_HPP_
#define QREG_INTERFACE_IFORC_HPP_



#include "ifbase.hpp"



namespace qreg
{



  //
  // Interface to the ORC program.
  //
  class IfOrc final : public IfBase
  {



    public:


      //
      // Sets up interface to one ORC calculation.
      //
      IfOrc(std::string newBinaryFilePath  ,
            std::string newHeaderString    ,
            std::string newOutputFolderPath )noexcept;

      //
      // Copy construction disabled.
      //
      IfOrc(const IfOrc& ifOrc) = delete;

      //
      // Copy assignement disabled.
      //
      auto
      operator=(const IfOrc& ifOrc)
      ->IfOrc& = delete;

      //
      // Move construction enabled.
      //
      IfOrc(IfOrc&& ifOrc);

      //
      // Move assignement enabled.
      //
      auto
      operator=(IfOrc&& ifOrc)
      ->IfOrc&;

      //
      // Destructor.
      //
      ~IfOrc(void);


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(IfOrc& lhs,
           IfOrc& rhs )
      ->void;


      //
      // Harvests ORC calculation.
      //
      auto
      harvest(Box& box)
      ->bool override;

      //
      // Seeds ORC calculation.
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
      // Returns delimiter message for beginning of charge.
      //
      static auto
      getChargeBeginMessage(void)noexcept
      ->std::string;

      //
      // Returns delimiter message for end of charge.
      //
      static auto
      getChargeEndMessage(void)noexcept
      ->std::string;

      //
      // Returns delimiter message for energy.
      //
      static auto
      getEnergyMessage(void)noexcept
      ->std::string;

      //
      // Returns delimiter message for end of gradient.
      //
      static auto
      getGradientEndMessage(void)noexcept
      ->std::string;

      //
      // Returns delimiter message for beginning of QM gradient.
      //
      static auto
      getQmGradientBeginMessage(void)noexcept
      ->std::string;

      //
      // Returns delimiter message for beginning of SE gradient.
      //
      static auto
      getSeGradientBeginMessage(void)noexcept
      ->std::string;

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
    IfOrc::getChargeBeginMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getChargeBeginMessage(void)noexcept
      // ->std::string
      //

      return "LOEWDIN ATOMIC CHARGES";


    }



    inline auto
    IfOrc::getChargeEndMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getChargeEndMessage(void)noexcept
      // ->std::string
      //

      return "LOEWDIN REDUCED ORBITAL CHARGES";


    }



    inline auto
    IfOrc::getEnergyMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getEnergyMessage(void)noexcept
      // ->std::string
      //

      return "FINAL";


    }



    inline auto
    IfOrc::getGradientEndMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getGradientEndMessage(void)noexcept
      // ->std::string
      //

      return "Norm";


    }



    inline auto
    IfOrc::getQmGradientBeginMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getQmGradientBeginMessage(void)noexcept
      // ->std::string
      //

      return "CARTESIAN GRADIENT";


    }



    inline auto
    IfOrc::getSeGradientBeginMessage(void)noexcept
    ->std::string
    {


      //
      // inline auto
      // IfOrc::getSeGradientBeginMessage(void)noexcept
      // ->std::string
      //

      return "The cartesian gradient:";


    }



    inline auto
    IfOrc::getEnergyConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOrc::getEnergyConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree to E -23 J
      return 4.359810e+5;


    }



    inline auto
    IfOrc::getGradientConversionFactor(void)const noexcept
    ->double
    {


      //
      // inline auto
      // IfOrc::getGradientConversionFactor(void)const noexcept
      // ->double
      //

      // Hartree to E -23 J
      return 4.359810e+5;


    }



}



#endif // QREG_INTERFACE_IFORC_HPP_

