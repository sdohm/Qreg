//
// partitioner.hpp
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



#ifndef QREG_PARTITIONER_HPP_
#define QREG_PARTITIONER_HPP_



#include "box.hpp"



namespace qreg
{



  //
  // Provides algorithms to partition quantum regions.
  //
  class Partitioner final
  {



    public:


      //
      // Policy class has no constructor.
      //
      Partitioner() = delete;

      //
      // Copy construction disabled.
      //
      Partitioner(const Partitioner& partitioner) = delete;

      //
      // Move construction disabled.
      //
      Partitioner(Partitioner&& partitioner) = delete;

      //
      // Copy assignment disabled.
      //
      auto
      operator=(const Partitioner& partitioner)
      ->Partitioner& = delete;

      //
      // Move assignment disabled.
      //
      auto
      operator=(Partitioner&& partitioner)
      ->Partitioner& = delete;

      //
      // Destruction disabled.
      //
      ~Partitioner(void) = delete;


      //
      // Linear buffer.
      //
      static auto
      computeLinearBuffer(Box&                box        ,
                          const unsigned int& markerNo   ,
                          const double&       innerRadius,
                          const double&       outerRadius )noexcept
      ->void;

      //
      // Smoothstep buffer.
      //
      static auto
      computeSmoothStepBuffer(Box&                box        ,
                              const unsigned int& markerNo   ,
                              const double&       innerRadius,
                              const double&       outerRadius )noexcept
      ->void;

      //
      // Smootherstep buffer.
      //
      static auto
      computeSmootherStepBuffer(Box&                box        ,
                                const unsigned int& markerNo   ,
                                const double&       innerRadius,
                                const double&       outerRadius )noexcept
      ->void;



  };



}



#endif // QREG_PARTITIONER_HPP_
