//
// partitioner.cpp
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



#include "partitioner.hpp"



namespace qreg
{



  auto
  Partitioner::computeLinearBuffer(Box&                box        ,
                                   const unsigned int& markerNo   ,
                                   const double&       innerRadius,
                                   const double&       outerRadius )noexcept
  ->void
  {


    //
    // DbC PRE
    //
   
    assert(innerRadius < outerRadius);
    assert(box.hasMarkerNo(markerNo));


    //
    // auto
    // Partitioner::computeLinearBuffer(Box&                box        ,
    //                                  const unsigned int& markerNo   ,
    //                                  const double&       innerRadius,
    //                                  const double&       outerRadius )noexcept
    // ->void
    //

    const auto nTargetAtomNos = box.getNAtomNosByMarkerNo(markerNo);

    assert(nTargetAtomNos.size() == 1);
    if(nTargetAtomNos.size() != 1)
      return;

    const auto targetAtomNo = nTargetAtomNos.front();

    auto nBufferedAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    outerRadius  );

    const auto nQmAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    innerRadius  );

    assert(!nQmAtomNos.empty());
    if(nQmAtomNos.empty())
      return;
    
    for(const auto& qmAtomNo:nQmAtomNos)
    {

      nBufferedAtomNos.erase
       (std::remove(nBufferedAtomNos.begin(),
                    nBufferedAtomNos.end  (),
                    qmAtomNo                 ),
                   nBufferedAtomNos.end()      );

    }

    for(const auto& bufferedAtomNo:nBufferedAtomNos)
    {


      assert((box.getXGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getYGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getZGradientGuess(bufferedAtomNo) != 0.0)   );

      if((box.getXGradient(bufferedAtomNo) != 0.0) &&
         (box.getYGradient(bufferedAtomNo) != 0.0) &&
         (box.getZGradient(bufferedAtomNo) != 0.0)   )
        box.setGradient(bufferedAtomNo,
         (box.getDistance(targetAtomNo                ,
                          bufferedAtomNo ) /
         (innerRadius - outerRadius      )  +
         outerRadius                 /
         (outerRadius - innerRadius)         ) *
         box.getXGradient(bufferedAtomNo)        +
         (1 - (box.getDistance(targetAtomNo  ,
                               bufferedAtomNo ) /
         (innerRadius - outerRadius)             +
         outerRadius /
         (outerRadius - innerRadius)              )) *
         box.getXGradientGuess(bufferedAtomNo)        ,
  
         (box.getDistance(targetAtomNo  ,
                          bufferedAtomNo ) /
         (innerRadius - outerRadius      )  +
         outerRadius                 /
         (outerRadius - innerRadius)         ) *
         box.getYGradient(bufferedAtomNo)        +
         (1 - (box.getDistance(targetAtomNo  ,
                               bufferedAtomNo ) /
         (innerRadius - outerRadius)             +
         outerRadius /
         (outerRadius - innerRadius)              )) *
         box.getYGradientGuess(bufferedAtomNo)        ,
  
         (box.getDistance(targetAtomNo  ,
                          bufferedAtomNo ) /
         (innerRadius - outerRadius      )  +
         outerRadius                 /
         (outerRadius - innerRadius)         ) *
         box.getZGradient(bufferedAtomNo)        +
         (1 - (box.getDistance(targetAtomNo  ,
                               bufferedAtomNo ) /
         (innerRadius - outerRadius)             +
         outerRadius /
         (outerRadius - innerRadius)              )) *
         box.getZGradientGuess(bufferedAtomNo)         );


    }


  }



  auto
  Partitioner::computeSmoothStepBuffer(Box&                box        ,
                                       const unsigned int& markerNo   ,
                                       const double&       innerRadius,
                                       const double&       outerRadius )noexcept
  ->void
  {


    //
    // DbC PRE
    //

    assert(innerRadius < outerRadius);
    assert(box.hasMarkerNo(markerNo));


    //
    // auto
    // Partitioner::computeSmoothStepBuffer(Box&                box        ,
    //                                      const unsigned int& markerNo   ,
    //                                      const double&       innerRadius,
    //                                      const double&       outerRadius )noexcept
    // ->void
    //

    const auto nTargetAtomNos = box.getNAtomNosByMarkerNo(markerNo);

    assert(nTargetAtomNos.size() == 1);
    if(nTargetAtomNos.size() != 1)
      return;

    const auto targetAtomNo = nTargetAtomNos.front();

    auto nBufferedAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    outerRadius  );

    const auto nQmAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    innerRadius  );

    assert(!nQmAtomNos.empty());
    if(nQmAtomNos.empty())
      return;

    for(const auto& qmAtomNo:nQmAtomNos)
    {

      nBufferedAtomNos.erase
       (std::remove(nBufferedAtomNos.begin(),
                    nBufferedAtomNos.end  (),
                    qmAtomNo                 ),
                   nBufferedAtomNos.end()      );

    }

    for(const auto& bufferedAtomNo:nBufferedAtomNos)
    {


      assert((box.getXGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getYGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getZGradientGuess(bufferedAtomNo) != 0.0)   );

      if((box.getXGradient(bufferedAtomNo) != 0.0) &&
         (box.getYGradient(bufferedAtomNo) != 0.0) &&
         (box.getZGradient(bufferedAtomNo) != 0.0)   )
      {


        const auto radius =
         (box.getDistance(bufferedAtomNo,
                          targetAtomNo   ) -
          outerRadius                       ) /
         (innerRadius - outerRadius         )  ;

        assert((box.getXGradientGuess(bufferedAtomNo) != 0.0) &&
               (box.getYGradientGuess(bufferedAtomNo) != 0.0) &&
               (box.getZGradientGuess(bufferedAtomNo) != 0.0)   );

        if((box.getXGradient(bufferedAtomNo) != 0.0) &&
           (box.getYGradient(bufferedAtomNo) != 0.0) &&
           (box.getZGradient(bufferedAtomNo) != 0.0)   )
          box.setGradient(bufferedAtomNo                ,
           (radius * radius * (3 - 2 * radius)) *
           box.getXGradient(bufferedAtomNo)            +
           (1 - (radius * radius * (3 - 2 * radius))) *
           box.getXGradientGuess(bufferedAtomNo)        ,

           (radius * radius * (3 - 2 * radius)) *
           box.getYGradient(bufferedAtomNo)            +
           (1 - (radius * radius * (3 - 2 * radius))) *
           box.getYGradientGuess(bufferedAtomNo)        ,

           (radius * radius * (3 - 2 * radius)) *
           box.getZGradient(bufferedAtomNo)            +
           (1 - (radius * radius * (3 - 2 * radius))) *
           box.getZGradientGuess(bufferedAtomNo)         );


      }


    }


  }



  auto
  Partitioner::computeSmootherStepBuffer(Box&                box        ,
                                         const unsigned int& markerNo   ,
                                         const double&       innerRadius,
                                         const double&       outerRadius )noexcept
  ->void
  {


    //
    // DbC PRE
    //

    assert(innerRadius < outerRadius);
    assert(box.hasMarkerNo(markerNo));


    //
    // auto
    // Partitioner::computeSmootherStepBuffer(Box&                box        ,
    //                                        const unsigned int& markerNo   ,
    //                                        const double&       innerRadius,
    //                                        const double&       outerRadius )noexcept
    // ->void
    //

    const auto nTargetAtomNos = box.getNAtomNosByMarkerNo(markerNo);

    assert(nTargetAtomNos.size() == 1);
    if(nTargetAtomNos.size() != 1)
      return;

    const auto targetAtomNo = nTargetAtomNos.front();

    auto nBufferedAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    outerRadius  );

    const auto nQmAtomNos = box.getNNeighborAtomNos(targetAtomNo,
                                                    innerRadius  );

    assert(!nQmAtomNos.empty());
    if(nQmAtomNos.empty())
      return;

    for(const auto& qmAtomNo:nQmAtomNos)
    {

      nBufferedAtomNos.erase
       (std::remove(nBufferedAtomNos.begin(),
                    nBufferedAtomNos.end  (),
                    qmAtomNo                 ),
                   nBufferedAtomNos.end()      );

    }

    for(const auto& bufferedAtomNo:nBufferedAtomNos)
    {


      assert((box.getXGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getYGradientGuess(bufferedAtomNo) != 0.0) &&
             (box.getZGradientGuess(bufferedAtomNo) != 0.0)   );

      if((box.getXGradient(bufferedAtomNo) != 0.0) &&
         (box.getYGradient(bufferedAtomNo) != 0.0) &&
         (box.getZGradient(bufferedAtomNo) != 0.0)   )
      {


        const auto radius =
         (box.getDistance(bufferedAtomNo,
                          targetAtomNo   ) -
          outerRadius                       ) /
         (innerRadius - outerRadius         )  ;

        assert((box.getXGradientGuess(bufferedAtomNo) != 0.0) &&
               (box.getYGradientGuess(bufferedAtomNo) != 0.0) &&
               (box.getZGradientGuess(bufferedAtomNo) != 0.0)   );

        if((box.getXGradient(bufferedAtomNo) != 0.0) &&
           (box.getYGradient(bufferedAtomNo) != 0.0) &&
           (box.getZGradient(bufferedAtomNo) != 0.0)   )
          box.setGradient(bufferedAtomNo                                          ,
           (radius * radius * radius * (radius * (radius * 6 - 15) + 10)) *
            box.getXGradient(bufferedAtomNo)                                     +
            (1 - (radius * radius * radius * (radius * (radius * 6 - 15) +10))) *
            box.getXGradientGuess(bufferedAtomNo)                                 ,

           (radius * radius * radius * (radius * (radius * 6 - 15) + 10)) *
            box.getYGradient(bufferedAtomNo)                                     +
            (1 - (radius * radius * radius * (radius * (radius * 6 - 15) +10))) *
            box.getYGradientGuess(bufferedAtomNo)                                 ,

           (radius * radius * radius * (radius * (radius * 6 - 15) + 10)) *
            box.getZGradient(bufferedAtomNo)                                     +
            (1 - (radius * radius * radius * (radius * (radius * 6 - 15) +10))) *
            box.getZGradientGuess(bufferedAtomNo)                                  );


      }


    }


  }



}

