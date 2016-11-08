//
// box.cpp
//
//  Created on: Aug 12, 2014
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



#include "box.hpp"



namespace qreg
{



  Box::Box(double newBoxXDimension,
           double newBoxXOrigin   ,
           double newBoxYDimension,
           double newBoxYOrigin   ,
           double newBoxZDimension,
           double newBoxZOrigin   ,
           double newTimeFrame    ,
           double newTotalEPot     )noexcept:

    boxXDimension_   (newBoxXDimension),
    boxXOrigin_      (newBoxXOrigin   ),
    boxYDimension_   (newBoxYDimension),
    boxYOrigin_      (newBoxYOrigin   ),
    boxZDimension_   (newBoxZDimension),
    boxZOrigin_      (newBoxZOrigin   ),
    dampeningFactor_ (100             ),
    temperature_     (298.15          ),
    timeFrame_       (newTimeFrame    ),
    totalEPot_       (newTotalEPot    ),

    nAtomNos_(),

    nAtoms_(),

    nEPots_                (),
    nEvbTotalEPots_        (),
    mnEvbAtomNos_          (),
    mnMarkerAtomNos_       (),
    mnQuantumRegionAtomNos_(),

    mnEvbEPots_(),

    evbNoCount_          (0),
    markerNoCount_       (0),
    quantumRegionNoCount_(0),

    atomCount_(0)

  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


  }



  Box::Box(const Box& box)noexcept
  {


    //
    // Design by Contract Invariant
    //
    assert(box.isSane());


    // Calling copy assignment operator.
    *this = box;


    //
    // Design by Contract Invariant
    //

    assert(box.isSane());
    assert(    isSane());


  }



  Box::Box(Box&& box)noexcept
  {

    
    //
    // Design by Contract Invariant
    //
    assert(box.isSane());


    // Calling move assignment operator.
    *this = std::move(box);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


  }



  auto
  Box::operator=(const Box& box)noexcept
  ->Box&
  {


    //
    // Design by Contract Invariant
    //
    assert(box.isSane());


    //
    // Perform copy assignment.
    //

    boxXDimension_    = box.boxXDimension_   ;
    boxXOrigin_       = box.boxXOrigin_      ;
    boxYDimension_    = box.boxYDimension_   ;
    boxYOrigin_       = box.boxYOrigin_      ;
    boxZDimension_    = box.boxZDimension_   ;
    boxZOrigin_       = box.boxZOrigin_      ;
    dampeningFactor_  = box.dampeningFactor_ ;
    temperature_      = box.temperature_     ;
    timeFrame_        = box.timeFrame_       ;
    totalEPot_        = box.totalEPot_       ;

    nAtomNos_ = box.nAtomNos_;

    nAtoms_ = box.nAtoms_;

    nEPots_                 = box.nEPots_                ;
    nEvbTotalEPots_         = box.nEvbTotalEPots_        ;
    mnEvbAtomNos_           = box.mnEvbAtomNos_          ;
    mnMarkerAtomNos_        = box.mnMarkerAtomNos_       ;
    mnQuantumRegionAtomNos_ = box.mnQuantumRegionAtomNos_;

    mnEvbEPots_ = box.mnEvbEPots_;

    evbNoCount_           = box.evbNoCount_          ;
    markerNoCount_        = box.markerNoCount_       ;
    quantumRegionNoCount_ = box.quantumRegionNoCount_;

    atomCount_ = box.atomCount_;


    //
    // Design by Contract Invariant
    //

    assert(box.isSane());
    assert(    isSane());


    return *this;


  }



  auto
  Box::operator=(Box&& box)noexcept
  ->Box&
  {


    //
    // Design by Contract Invariant
    //
    assert(box.isSane());


    //
    // Perform move assignment.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    std::swap(    boxXDimension_,
              box.boxXDimension_ );

    std::swap(    boxXOrigin_,
              box.boxXOrigin_ );

    std::swap(    boxYDimension_,
              box.boxYDimension_ );

    std::swap(    boxYOrigin_,
              box.boxYOrigin_ );

    std::swap(    boxZDimension_,
              box.boxZDimension_ );

    std::swap(    boxZOrigin_,
              box.boxZOrigin_ );

    std::swap(    dampeningFactor_,
              box.dampeningFactor_ );

    std::swap(    temperature_,
              box.temperature_ );

    std::swap(    timeFrame_,
              box.timeFrame_ );

    std::swap(    totalEPot_,
              box.totalEPot_ );

    std::swap(    nAtomNos_,
              box.nAtomNos_ );

    std::swap(    nAtoms_,
              box.nAtoms_ );

    std::swap(    nEPots_,
              box.nEPots_ );

    std::swap(    nEvbTotalEPots_,
              box.nEvbTotalEPots_ );

    std::swap(    mnEvbAtomNos_,
              box.mnEvbAtomNos_ );

    std::swap(    mnMarkerAtomNos_,
              box.mnMarkerAtomNos_ );

    std::swap(    mnQuantumRegionAtomNos_,
              box.mnQuantumRegionAtomNos_ );

    std::swap(    mnEvbEPots_,
              box.mnEvbEPots_ );

    std::swap(    evbNoCount_,
              box.evbNoCount_ );

    std::swap(    markerNoCount_,
              box.markerNoCount_ );

    std::swap(    quantumRegionNoCount_,
              box.quantumRegionNoCount_ );

    std::swap(    atomCount_,
              box.atomCount_ );

    box.~Box();


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    return *this;


  }



  Box::~Box(void)noexcept
  {



  }



  auto
  Box::cloneNonPeriodic(const double& range)const noexcept
  ->Box
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::cloneNonPeriodic(const double& range)const noexcept
    // ->Box
    //

    Box clone = *this;

    // Shift coords by randomly chosen atom coords.
    srand(static_cast<unsigned>(time(NULL)));
    const size_t centerAtomNo = rand() % getNAtomNos().size() + 1;
    const double xOrigin = getXCoordinate(centerAtomNo);
    const double yOrigin = getYCoordinate(centerAtomNo);
    const double zOrigin = getZCoordinate(centerAtomNo);
    for(const auto& atomNo:getNAtomNos())
      clone.setCoordinate(atomNo                                ,
                          clone.getXCoordinate(atomNo) - xOrigin,
                          clone.getYCoordinate(atomNo) - yOrigin,
                          clone.getZCoordinate(atomNo) - zOrigin );

    // Get rid of periodicity.
    clone.setBoxDimension(0.0,
                          0.0,
                          0.0 );

    // Add atoms by range.
    auto i = getNAtomNos().size();
    std::map<size_t,
             size_t > nExternalInternalAtomNos;
    for(const auto& atomNo:getNAtomNos())
    {


      const auto nullXCoord = clone.getXCoordinate(atomNo);
      const auto nullYCoord = clone.getYCoordinate(atomNo);
      const auto nullZCoord = clone.getZCoordinate(atomNo);

      const auto minusXCoord = nullXCoord - getBoxXDimension();
      const auto minusYCoord = nullYCoord - getBoxYDimension();
      const auto minusZCoord = nullZCoord - getBoxZDimension();

      const auto plusXCoord = nullXCoord + getBoxXDimension();
      const auto plusYCoord = nullYCoord + getBoxYDimension();
      const auto plusZCoord = nullZCoord + getBoxZDimension();

      const auto xMin = getBoxXOrigin() - range;
      const auto yMin = getBoxYOrigin() - range;
      const auto zMin = getBoxZOrigin() - range;

      const auto xMax = getBoxXOrigin() + getBoxXDimension() + range;
      const auto yMax = getBoxYOrigin() + getBoxYDimension() + range;
      const auto zMax = getBoxYOrigin() + getBoxZDimension() + range;

      auto augmentBox = [&](const double& newXCoord,
                            const double& newYCoord,
                            const double& newZCoord )
      ->void
      {
        if(((newXCoord > xMin) && (newXCoord < xMax)) &&
           ((newYCoord > yMin) && (newYCoord < yMax)) &&
           ((newZCoord > zMin) && (newZCoord < zMax))   )
        {

          ++i;
          nExternalInternalAtomNos.emplace(i     ,
                                           atomNo );
          clone.insertAtom(i                        ,
                           getChemicalSymbol(atomNo),
                           getElectronCount(atomNo) ,
                           newXCoord                ,
                           newYCoord                ,
                           newZCoord                 ); 

        }
      };

      // 0 0-Z
      augmentBox(nullXCoord, nullYCoord, minusZCoord);
      // 0 0+Z
      augmentBox(nullXCoord, nullYCoord, plusZCoord );
      // 0-Y 0
      augmentBox(nullXCoord, minusYCoord,nullZCoord );
      // 0-Y-Z
      augmentBox(nullXCoord, minusYCoord,minusZCoord);
      // 0-Y+Z
      augmentBox(nullXCoord, minusYCoord,plusZCoord );
      // 0+Y 0
      augmentBox(nullXCoord, plusYCoord, nullZCoord );
      // 0+Y-Z
      augmentBox(nullXCoord, plusYCoord, minusZCoord);
      // 0+Y+Z
      augmentBox(nullXCoord, plusYCoord, plusZCoord );
      //-X 0 0 
      augmentBox(minusXCoord,nullYCoord, nullZCoord );
      //-X 0-Z
      augmentBox(minusXCoord,nullYCoord, minusZCoord);
      //-X 0+Z
      augmentBox(minusXCoord,nullYCoord, plusZCoord );
      //-X-Y 0
      augmentBox(minusXCoord,minusYCoord,nullZCoord );
      //-X-Y-Z
      augmentBox(minusXCoord,minusYCoord,minusZCoord);
      //-X-Y+Z
      augmentBox(minusXCoord,minusYCoord,plusZCoord );
      //-X+Y 0
      augmentBox(minusXCoord,plusYCoord, nullZCoord );
      //-X+Y-Z
      augmentBox(minusXCoord,plusYCoord, minusZCoord);
      //-X+Y+Z
      augmentBox(minusXCoord,plusYCoord, plusZCoord );
      //+X 0 0
      augmentBox(plusXCoord, nullYCoord, nullZCoord );
      //+X 0-Z
      augmentBox(plusXCoord, nullYCoord, minusZCoord);
      //+X 0+Z
      augmentBox(plusXCoord, nullYCoord,plusZCoord  );
      //+X-Y 0
      augmentBox(plusXCoord, minusYCoord,nullZCoord );
      //+X-Y-Z
      augmentBox(plusXCoord, minusYCoord,minusZCoord);
      //+X-Y+Z
      augmentBox(plusXCoord, minusYCoord,plusZCoord );
      //+X+Y 0
      augmentBox(plusXCoord, plusYCoord,nullZCoord  );
      //+X+Y-Z
      augmentBox(plusXCoord, plusYCoord,minusZCoord );
      //+X+Y+Z
      augmentBox(plusXCoord, plusYCoord,plusZCoord  );


    }


    //
    // Design by Contract Postcondition
    //
    assert(isSane());    


    return clone;


  }



  auto
  Box::cloneShrinkGraphite(const double& range)const noexcept
  -> Box
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::cloneShrinkGraphite(const double& range)const noexcept
    // ->Box
    //

    Box clone = *this;

    std::vector<size_t> nGraphiteAtoms;
// DUMMY
if(range == 0.12345)
  {};
    auto graphiteCaseDetected = [&]()
    ->bool
    {
// DETECT GRAPHITE CASE
      return false;
    };

    auto cut = [&]()
    ->bool
    {
// REMOVE ALL QM STATES ON GRAPHITE C ATOMS OUTSIDE OF RANGE TO MARKED ATOM
      return false;
    };

    auto saturate = [&]()
    ->void
    {
// REPLACE OUTER GRAPHITE C ATOMS WITH H ATOMS VIA INSERTING ADDITIONAL ATOMS
    };

    if((graphiteCaseDetected()) &&
       (cut                 ())   )
      saturate();


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


    return clone;


  }



  auto
  Box::computeAllDistances(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::computeAllDistances(void)noexcept
    // ->void
    //

    std::vector<std::thread> nThreads;

    nThreads.reserve(atomCount_);

    // Loop over all atoms fixes DbC violation caused by
    // Box::computeAllDistancesForOneAtom.
    for(auto&& i:nAtomNos_)
      nThreads.push_back(std::thread
                         (&qreg::Box::computeAllDistancesForOneAtom,
                          this                                     ,
                          i                                        ,
                          nAtomNos_                                 ));

    for(auto&& i:nThreads)
      i.join();


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Box::computeNearestNeighborDistances(const double& distance)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(distance));
    assert(distance == distance);
    assert(distance > 0.0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::computeNearestNeighborDistances
    // (const unsigned int& nearestNeighborCount)noexcept
    // ->void
    //

    std::vector<std::thread> nThreads;

    nThreads.reserve(atomCount_);

    // Loop over all atoms fixes DbC violation caused by
    // Box::computeAllDistancesForOneAtom.
    for(auto&& i:nAtomNos_)
      nThreads.push_back(std::thread
                      (&qreg::Box::computeNearestNeighborDistancesForOneAtom,
                       this                                                 ,
                       i                                                    ,
                       distance                                              ));

    for(auto&& i:nThreads)
      i.join();


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Box::deleteAllDistances(const size_t& atomNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllDistances(const size_t& atomNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(auto&& i:nAtoms_)
      if((i.getAtomNo        (      ) != atomNo) &&
	 (i.hasNeighborAtomNo(atomNo)          )   )
        i.deleteDistance(atomNo);

    nAtoms_[atomNo - 1].deleteAllDistances();


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getNNeighborAtomNos().empty());

    // Atom with atomNo has no neighbors.
    assert(([this  ,
             atomNo ]()->bool
    {
      for(const auto& i:nAtoms_)
	if(i.hasNeighborAtomNo(atomNo))
	  return false;
      return true;
    })());

    assert(isSane());


  }



  auto
  Box::deleteAllDistances(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllDistances(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(auto&& i:nAtoms_)
      i.deleteAllDistances();


    //
    // Design by Contract Postcondition
    //

    // No atom has neighbors.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
	if(!i.getNNeighborAtomNos().empty())
	  return false;
      return true;
    })());

    assert(isSane());


  }



  auto
  Box::deleteAllEPots(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllEPots(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nEPots_.clear();


    //
    // Design by Contract Postcondition
    //

    assert(nEPots_.empty());

    assert(isSane());


  }



  auto
  Box::deleteAllEvbEPots(const unsigned int& evbNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEvbEPots(const unsigned int& evbNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnEvbEPots_.erase(evbNo);


    //
    // Design by Contract Postcondition
    //

    // Box:mnEvbEPots_ not assigned key value evbNo.
    assert(mnEvbEPots_.count(evbNo) == 0);

    assert(isSane());


  }



  auto
  Box::deleteAllEvbEPots(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllEvbEPots(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnEvbEPots_.clear();


    //
    // Design by Contract Postcondition
    //

    assert(mnEvbEPots_.empty());

    assert(isSane());


  }



  auto
  Box::deleteAllEvbGradients(const size_t& atomNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllEvbGradients(const size_t& atomNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(auto&& i:mnEvbAtomNos_)
    {

      i.second.erase(std::remove(i.second.begin(),
                                 i.second.end  (),
                                 atomNo           ), i.second.end());

      if(i.second.empty())
	mnEvbAtomNos_.erase(i.first);

    }

    evbNoCount_ = mnEvbAtomNos_.size();

    nAtoms_[atomNo - 1].deleteAllEvbGradients();


    //
    // Design by Contract Postcondition
    //

    // Box::mnEvbAtomNos_ assigned mapped values do not contain atomNo.
    assert(([this  ,
             atomNo ]()->bool
    {
      for(const auto& i:mnEvbAtomNos_)
	if(std::find(i.second.cbegin(),
	             i.second.cend  (),
	             atomNo            ) != i.second.cend())
	  return false;
      return true;
    })());

    assert(isSane());


  }



  auto
  Box::deleteAllEvbGradients(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllEvbGradients(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnEvbAtomNos_.clear();

    evbNoCount_ = mnEvbAtomNos_.size();

    for(auto&& i:nAtoms_)
      i.deleteAllEvbGradients();


    //
    // Design by Contract Postcondition
    //

    assert(evbNoCount_ == 0);

    assert(isSane());


  }



  auto
  Box::deleteAllEvbTotalEpots(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllEvbTotalEpots(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nEvbTotalEPots_.clear();


    //
    // Design by Contract Postcondition
    //

    assert(nEvbTotalEPots_.empty());

    assert(isSane());


  }



  auto
  Box::deleteAllMarkerNos(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllMarkerNos(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnMarkerAtomNos_.clear();

    markerNoCount_ = mnMarkerAtomNos_.size();

    for(auto&& i:nAtoms_)
      i.deleteAllMarkerNos();


    //
    // Design by Contract Postcondition
    //

    assert(markerNoCount_ == 0);

    assert(isSane());


  }



  auto
  Box::deleteAllQuantumRegionNos(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteAllQuantumRegionNos(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnQuantumRegionAtomNos_.clear();

    quantumRegionNoCount_ = mnQuantumRegionAtomNos_.size();

    for(auto&& i:nAtoms_)
      i.deleteAllQuantumRegionNos();


    //
    // Design by Contract Postcondition
    //

    assert(quantumRegionNoCount_ == 0);

    assert(isSane());


  }



  auto
  Box::deleteDistance(const size_t& atomNo        ,
                      const size_t& neighborAtomNo )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo         != 0     );
    assert(neighborAtomNo != 0     );
    assert(neighborAtomNo != atomNo);

    // Box::nAtomNos assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::nAtomNos assigned neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.cend());

    assert(nAtoms_[atomNo         - 1]
                  .hasNeighborAtomNo(neighborAtomNo));
    assert(nAtoms_[neighborAtomNo - 1]
                  .hasNeighborAtomNo(atomNo        ));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteDistance(const size_t& atomNo        ,
    //                     const size_t& neighborAtomNo )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo         - 1].deleteDistance(neighborAtomNo);
    nAtoms_[neighborAtomNo - 1].deleteDistance(atomNo        );


    //
    // Design by Contract Postcondition
    //

    assert(!nAtoms_[atomNo         - 1].hasNeighborAtomNo(neighborAtomNo));
    assert(!nAtoms_[neighborAtomNo - 1].hasNeighborAtomNo(atomNo        ));

    assert(isSane());


  }



  auto
  Box::deleteEPot(const unsigned int& quantumRegionNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    // Box:mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);

    // Box::nEPots_ assigned key value quantumRegionNo.
    assert(nEPots_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEPot(const unsigned int& quantumRegionNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nEPots_.erase(quantumRegionNo);


    //
    // Design by Contract Postcondition
    //

    // Box::nEPots_ not assigned key value quantumRegionNo.
    assert(nEPots_.count(quantumRegionNo) == 0);

    assert(isSane());


  }



  auto
  Box::deleteEvbEPot(const unsigned int& evbNo          ,
                     const unsigned int& quantumRegionNo )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo           != 0);
    assert(quantumRegionNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);

    // Box::mnEvbEPots_ assigned key value evbNo.
    assert(mnEvbEPots_.count(evbNo) != 0);

    // Box::mnEvbEPots_ assigned mapped value
    // assigned key value quantumRegionNo.
    assert(mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEvbEPot(const unsigned int& evbNo          ,
    //                    const unsigned int& quantumRegionNo )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnEvbEPots_.find(evbNo)->second.erase(quantumRegionNo);

    if(mnEvbEPots_.find(evbNo)->second.empty())
      mnEvbEPots_.erase(evbNo);


    //
    // Design by Contract Postcondition
    //

    // Requested EPot has been deleted.
    assert(([this           ,
             evbNo          ,
             quantumRegionNo ]()->bool
    {
      if((mnEvbEPots_.count(evbNo)                               != 0) &&
	 (mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) != 0)   )
	return false;
      return true;
    })());

    assert(isSane());


  }



  auto
  Box::deleteEvbGradient(const size_t&       atomNo,
                         const unsigned int& evbNo  )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtoms_ assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    assert(nAtoms_[atomNo - 1].hasEvbNo(evbNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEvbGradient(const size_t&       atomNo,
    //                        const unsigned int& evbNo  )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(mnEvbAtomNos_.find(evbNo)->second.empty())
      mnEvbAtomNos_.erase(evbNo);

    evbNoCount_ = mnEvbAtomNos_.size();

    nAtoms_[atomNo - 1].deleteEvbGradient(evbNo);


    //
    // Design by Contract Postcondition
    //

    assert(!nAtoms_[atomNo - 1].hasEvbNo(evbNo));

    assert(isSane());


  }



  auto
  Box::deleteEvbGradient(const unsigned int& evbNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEvbGradient(const unsigned int& evbNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    mnEvbAtomNos_.erase(evbNo);

    evbNoCount_ = mnEvbAtomNos_.size();

    for(auto&& i:nAtoms_)
      if(i.hasEvbNo(evbNo))
	i.deleteEvbGradient(evbNo);


    //
    // Design by Contract Postcondition
    //

    // Box::mnEvbAtomNos_ not assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) == 0);

    assert(isSane());


  }



  auto
  Box::deleteEvbTotalEPot(const unsigned int& evbNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    // Box::nEvbTotalEPots_ assigned key value evbNo.
    assert(nEvbTotalEPots_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteEvbTotalEPot(const unsigned int& evbNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nEvbTotalEPots_.erase(evbNo);


    //
    // Design by Contract Postcondition
    //

    // Box::nEvbTotalEPots_ not assigned key value evbNo.
    assert(nEvbTotalEPots_.count(evbNo) == 0);

    assert(isSane());


  }



  auto
  Box::deleteMarkerNo(const unsigned int& markerNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(markerNo != 0);

    // Box::mnMarkerAtomNos_ assigned key value markerNo.
    assert(mnMarkerAtomNos_.count(markerNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteMarkerNo(const unsigned int& markerNo)noexcept
    // ->void
    //
    
    std::lock_guard<std::mutex> lock(threadSafety_);

    for(const auto& i:mnMarkerAtomNos_.find(markerNo)->second)
      nAtoms_[i - 1].deleteMarkerNo(markerNo);

    mnMarkerAtomNos_.erase(markerNo);

    markerNoCount_ = mnMarkerAtomNos_.size();

    //
    // Design by Contract Postcondition
    //

    // Box::mnMarkerAtomNos_ not assigned key value markerNo.
    assert(mnMarkerAtomNos_.count(markerNo) == 0);

    assert(isSane());


  }



  auto
  Box::deleteQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::deleteQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(const auto& i:mnQuantumRegionAtomNos_.find(quantumRegionNo)->second)
      nAtoms_[i - 1].deleteQuantumRegionNo(quantumRegionNo);

    mnQuantumRegionAtomNos_.erase(quantumRegionNo);

    quantumRegionNoCount_ = mnQuantumRegionAtomNos_.size();


    //
    // Design by Contract Postcondition
    //

    // Box::mnQuantumRegionAtomNos_ not assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) == 0);

    assert(isSane());


  }



  auto
  Box::equalizeGradients(void)noexcept
  ->void
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::equalizeGradients(void)noexcept
    // ->void
    //

    for(auto&& atom:nAtoms_)
      atom.setGradient(0.0,
                       0.0,
                       0.0);


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Box::incrementNegativeCharge(const size_t& atomNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::incrementNegativeCharge(const size_t& atomNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].incrementNegativeCharge();


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Box::incrementPositiveCharge(const size_t& atomNo)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Atom::electronCount_ non-zero.
    assert(nAtoms_[atomNo - 1].getElectronCount() != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::incrementPositiveCharge(const size_t& atomNo)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].incrementPositiveCharge();


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Box::insertAtom(const size_t&       newAtomNo       ,
                  const std::string&  newChemicalName ,
                  const unsigned int& newElectronCount,
                  const double&       newXCoordinate  ,
                  const double&       newYCoordinate  ,
                  const double&       newZCoordinate   )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXCoordinate));
    assert(!std::isinf(newYCoordinate));
    assert(!std::isinf(newZCoordinate));

    assert(newXCoordinate == newXCoordinate);
    assert(newYCoordinate == newYCoordinate);
    assert(newZCoordinate == newZCoordinate);

    assert(newAtomNo != 0);

    // Box::nAtomNos_ not assigned newAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     newAtomNo          ) == nAtomNos_.cend());

    // newAtomNo consecutive number to present values of Box::nAtomNos_.
    assert(([this     ,
             newAtomNo ]()->bool
    {
      if(((nAtomNos_.empty()) &&
          (newAtomNo == 1   )                   ) ||
         ((!nAtomNos_.empty()               ) &&
          (newAtomNo == nAtomNos_.back() + 1)   )   )
        return true;
      return false; 
    })());

    assert(!newChemicalName.empty());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertAtom(const size_t&       newAtomNo       ,
    //                 const std::string&  newChemicalName ,
    //                 const unsigned int& newElectronCount,
    //                 const double&       newXCoordinate  ,
    //                 const double&       newYCoordinate  ,
    //                 const double&       newZCoordinate   )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    auto xCoordinate = newXCoordinate;
    auto yCoordinate = newYCoordinate;
    auto zCoordinate = newZCoordinate;

    // In case of periodic box.
    if((boxXDimension_ != 0.0) &&
       (boxYDimension_ != 0.0) &&
       (boxZDimension_ != 0.0)   )
    {
      // First atom added.
      if(newAtomNo == 1)
      {
        boxXOrigin_ = xCoordinate;
        boxYOrigin_ = yCoordinate;
        boxZOrigin_ = zCoordinate;
      }
      // Not on first atom added.
      else
      {
        // Check for box dimensions.
        while(xCoordinate >= boxXOrigin_ + boxXDimension_)
          xCoordinate -= boxXDimension_;
        while(xCoordinate < boxXOrigin_)
          xCoordinate += boxXDimension_;
        while(yCoordinate >= boxYOrigin_ + boxYDimension_)
          yCoordinate -= boxYDimension_;
        while(yCoordinate < boxYOrigin_)
          yCoordinate += boxYDimension_;
        while(zCoordinate >= boxZOrigin_ + boxZDimension_)
          zCoordinate -= boxZDimension_;
        while(zCoordinate < boxZOrigin_)
          zCoordinate += boxZDimension_;
      }
    }
    

    nAtomNos_.push_back(newAtomNo);


    if(newChemicalName == "Li")
      nAtoms_.push_back(Atom(newElectronCount,
                             3               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                              6.94           ,
                              0.0            ,
                              0.0            ,
                              1.82           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "Li"             ));

    else if(newChemicalName == "B")
      nAtoms_.push_back(Atom(newElectronCount,
                             5               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             10.81           ,
                              0.0            ,
                              0.0            ,
                              1.92           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "B"              ));

    else if(newChemicalName == "C")
      nAtoms_.push_back(Atom(newElectronCount,
                             6               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             12.011          ,
                              0.0            ,
                              0.0            ,
                              1.70           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "C"              ));

    else if(newChemicalName == "N")
      nAtoms_.push_back(Atom(newElectronCount,
                             7               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             14.0067         ,
                              0.0            ,
                              0.0            ,
                              1.55           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "N"              ));

    else if(newChemicalName == "O")
      nAtoms_.push_back(Atom(newElectronCount,
                             8               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             15.999          ,
                              0.0            ,
                              0.0            ,
                              1.52           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "O"              ));

    else if(newChemicalName == "F")
      nAtoms_.push_back(Atom(newElectronCount,
                             9               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             18.998403163    ,
                              0.0            ,
                              0.0            ,
                              1.47           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "F"              ));

    else if(newChemicalName == "Na")
      nAtoms_.push_back(Atom(newElectronCount,
                             11              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             22.98976928     ,
                              0.0            ,
                              0.0            ,
                              2.27           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "Na"             ));

    else if(newChemicalName == "Mg")
      nAtoms_.push_back(Atom(newElectronCount,
                             12              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             24.305          ,
                              0.0            ,
                              0.0            ,
                              1.73           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "Mg"             ));

    else if(newChemicalName == "Al")
      nAtoms_.push_back(Atom(newElectronCount,
                             13              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             26.9815385      ,
                              0.0            ,
                              0.0            ,
                              1.84           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "Al"             ));

    else if(newChemicalName == "Si")
      nAtoms_.push_back(Atom(newElectronCount,
                             14              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             28.085          ,
                              0.0            ,
                              0.0            ,
                              2.10           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "Si"             ));

    else if(newChemicalName == "P")
      nAtoms_.push_back(Atom(newElectronCount,
                             15              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             30.973761998    ,
                              0.0            ,
                              0.0            ,
                              1.80           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "P"              ));

    else if(newChemicalName == "S")
      nAtoms_.push_back(Atom(newElectronCount,
                             16              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             32.06           ,
                              0.0            ,
                              0.0            ,
                              1.80           ,
                             xCoordinate     ,
                              0.0            ,
                              0.0            ,
                             yCoordinate     ,
                              0.0            ,
                              0.0            ,
                             zCoordinate     ,
                              0.0            ,
                              0.0            ,
                             newChemicalName ,
                             "S"              ));

    else if(newChemicalName == "Pt")
      nAtoms_.push_back(Atom(newElectronCount,
                             78              ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             195.084         ,
                               0.0           ,
                               0.0           ,
                               1.75          ,
                             xCoordinate     ,
                               0.0           ,
                               0.0           ,
                             yCoordinate     ,
                               0.0           ,
                               0.0           ,
                             zCoordinate     ,
                               0.0           ,
                               0.0           ,
                             newChemicalName ,
                             "Pt"             ));

    else
      nAtoms_.push_back(Atom(newElectronCount,
                             1               ,
                             newAtomNo       ,
                             xCoordinate     ,
                             yCoordinate     ,
                             zCoordinate     ,
                             1.008           ,
                             0.0             ,
                             0.0             ,
                             1.20            ,
                             xCoordinate     ,
                             0.0             ,
                             0.0             ,
                             yCoordinate     ,
                             0.0             ,
                             0.0             ,
                             zCoordinate     ,
                             0.0             ,
                             0.0             ,
                             "H"             ,
                             "H"              ));

    atomCount_ = nAtoms_.size();


    //
    // Design by Contract Postcondition
    //

    // Box::nAtomNos_ assigned newAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     newAtomNo          ) != nAtomNos_.cend());

    assert(nAtoms_[newAtomNo - 1].getAtomNo() == newAtomNo);

    assert(isSane());


  }



  auto
  Box::insertDistance(const size_t& atomNo        ,
                      const size_t& neighborAtomNo,
                      const double& newDistance    )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newDistance));
    assert(newDistance == newDistance);

    assert(atomNo         != 0);
    assert(neighborAtomNo != 0);

    assert(atomNo != neighborAtomNo);

    // Box::nAtomNos_ assigned atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::nAtomNos_ assigned neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.cend());

    assert(newDistance >= 0.0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertDistance(const size_t& atomNo       ,
    //                     const size_t& neigborAtomNo,
    //                     const double& newDistance   )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].insertDistance(neighborAtomNo,
                                       newDistance    );

    nAtoms_[neighborAtomNo - 1].insertDistance(atomNo     ,
                                               newDistance );


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].hasNeighborAtomNo(neighborAtomNo));

    assert(nAtoms_[neighborAtomNo - 1].hasNeighborAtomNo(atomNo));

    assert(isSane());


  }



  auto
  Box::insertEPot(const unsigned int& quantumRegionNo,
                  const double&       newEPot         )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newEPot));
    assert(newEPot == newEPot);

    assert(quantumRegionNo != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertEPot(const unsigned int& quantumRegionNo,
    //                 const double&       newEPot         )noexcept
    // ->void
    //

    if(nEPots_.count(quantumRegionNo) != 0)
      deleteEPot(quantumRegionNo);

    std::lock_guard<std::mutex> lock(threadSafety_);

    nEPots_.emplace(quantumRegionNo,
                    newEPot         );


    //
    // Design by Contract Postcondition
    //

    // Box::nEPots_ assigned key value quantumRegionNo.
    assert(nEPots_.count(quantumRegionNo) != 0);

    // Box::nEPots_ assigned mapped value newEPot for key value quantumRegionNo.
    assert(nEPots_.find(quantumRegionNo)->second == newEPot);

    assert(isSane());


  }



  auto
  Box::insertEvbEPot(const unsigned int& evbNo          ,
                     const unsigned int& quantumRegionNo,
                     const double&       newEvbEPot      )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newEvbEPot));
    assert(newEvbEPot == newEvbEPot);

    assert(evbNo           != 0);
    assert(quantumRegionNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertEvbEPot(const unsigned int& evbNo          ,
    //                    const unsigned int& quantumRegionNo,
    //                    const double&       newEvbEPot      )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if((mnEvbEPots_.count(evbNo)                               != 0) &&
       (mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) != 0)   )
    {
      mnEvbEPots_.find(evbNo)->second.erase(quantumRegionNo);
      mnEvbEPots_.find(evbNo)->second.emplace(quantumRegionNo,
                                              newEvbEPot      );
    }

    if((mnEvbEPots_.count(evbNo)                               != 0) &&
       (mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) == 0)   )
      mnEvbEPots_.find(evbNo)->second.emplace(quantumRegionNo,
                                              newEvbEPot      );

    if(mnEvbEPots_.count(evbNo) == 0)
      mnEvbEPots_.emplace(evbNo                  ,
                          std::map<unsigned int,
                                   double       >
                          {
	                    {
	                      quantumRegionNo,
	                      newEvbEPot
	                    }
                          }                       );


    //
    // Design by Contract Postcondition
    //

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbEPots_.count(evbNo) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value evbNo assigned key value
    // quantumRegionNo.
    assert(mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) != 0);
    assert(mnEvbEPots_.find(evbNo)->second.find(quantumRegionNo)
                                   ->second
                                    == newEvbEPot               );


    assert(isSane());


  }



  auto
  Box::insertEvbGradient(const size_t&       atomNo         ,
                         const unsigned int& newEvbNo       ,
                         const double&       newEvbXGradient,
                         const double&       newEvbYGradient,
                         const double&       newEvbZGradient )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newEvbXGradient));
    assert(!std::isinf(newEvbYGradient));
    assert(!std::isinf(newEvbZGradient));

    assert(newEvbXGradient == newEvbXGradient);
    assert(newEvbYGradient == newEvbYGradient);
    assert(newEvbZGradient == newEvbZGradient);

    assert(atomNo   != 0);
    assert(newEvbNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertEvbGradient(const size_t&       atomNo         ,
    //                        const unsigned int& newEvbNo       ,
    //                        const double&       newEvbXGradient,
    //                        const double&       newEvbYGradient,
    //                        const double&       newEvbZGradient )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].insertEvbGradient(newEvbNo       ,
                                          newEvbXGradient,
                                          newEvbYGradient,
                                          newEvbZGradient );

    if((mnEvbAtomNos_.count(newEvbNo) != 0                       ) &&
       (std::find(mnEvbAtomNos_.find(newEvbNo)->second.cbegin(),
                  mnEvbAtomNos_.find(newEvbNo)->second.cend  (),
                  atomNo                                        )
         == mnEvbAtomNos_.find(newEvbNo)->second.cend()          )   )
      mnEvbAtomNos_.find(newEvbNo)->second.push_back(atomNo);

    if(mnEvbAtomNos_.count(newEvbNo) == 0)
      mnEvbAtomNos_.emplace(newEvbNo           ,
                            std::vector<size_t>
                             {
	                       atomNo
                             }                  );

    evbNoCount_ = mnEvbAtomNos_.size();


    //
    // Design by Contract Postcondition
    //

    // Box::mnEvbAtomNos_ assigned key value newEvbNo.
    assert(mnEvbAtomNos_.count(newEvbNo) != 0);

    // Box::nEvbAtomNos_ assigned key value newEvbNo
    // with mapped value atomNo.
    assert(std::find(mnEvbAtomNos_.find(newEvbNo)->second.cbegin(),
                     mnEvbAtomNos_.find(newEvbNo)->second.cend  (),
                     atomNo                                        )
            != mnEvbAtomNos_.find(newEvbNo)->second.cend()          );

    // Box::nAtoms_[atomNo - 1].hasEvbNo(newEvbNo) returns true.
    assert(nAtoms_[atomNo - 1].hasEvbNo(newEvbNo));

    // Box::nAtoms_[atomNo - 1].getEvbX/Y/ZGradient(newEvbNo)
    // returns newEvbX/Y/ZGradient.
    assert(nAtoms_[atomNo - 1].getEvbXGradient(newEvbNo) == newEvbXGradient);
    assert(nAtoms_[atomNo - 1].getEvbYGradient(newEvbNo) == newEvbYGradient);
    assert(nAtoms_[atomNo - 1].getEvbZGradient(newEvbNo) == newEvbZGradient);

    assert(isSane());


  }



  auto
  Box::insertEvbTotalEPot(const unsigned int& evbNo          ,
                          const double&       newEvbTotalEPot )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newEvbTotalEPot));
    assert(newEvbTotalEPot == newEvbTotalEPot);

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertEvbTotalEPot(const unsigned int& evbNo          ,
    //                         const double&       newEvbTotalEPot )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(nEvbTotalEPots_.count(evbNo) != 0)
      nEvbTotalEPots_.erase(evbNo);

    nEvbTotalEPots_.emplace(evbNo          ,
                            newEvbTotalEPot );


    //
    // Design by Contract Postcondition
    //

    // Box::nEvbTotalEPots_ assigned key value evbNo.
    assert(nEvbTotalEPots_.count(evbNo) != 0);

    // Box::nEvbTotalEPots_ assigned key value evbNo
    // with mapped value newEvbTotalEPot.
    assert(nEvbTotalEPots_.find(evbNo)->second == newEvbTotalEPot);

    assert(isSane());


  }



  auto
  Box::insertMarkerNo(const size_t&       atomNo     ,
                      const unsigned int& newMarkerNo )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo      != 0);
    assert(newMarkerNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertMarkerNo(const size_t&       atomNo     ,
    //                     const unsigned int& newMarkerNo )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if((mnMarkerAtomNos_.count(newMarkerNo) != 0                       ) &&
       (std::find(mnMarkerAtomNos_.find(newMarkerNo)->second.cbegin(),
                  mnMarkerAtomNos_.find(newMarkerNo)->second.cend  (),
                  atomNo                                              )
         == mnMarkerAtomNos_.find(newMarkerNo)->second.cend()          )   )
      mnMarkerAtomNos_.find(newMarkerNo)->second.push_back(atomNo);

    if(mnMarkerAtomNos_.count(newMarkerNo) == 0)
      mnMarkerAtomNos_.emplace(newMarkerNo        ,
                               std::vector<size_t>
                               {
	                         atomNo
                               }                   );

    markerNoCount_ = mnMarkerAtomNos_.size();

    nAtoms_[atomNo - 1].insertMarkerNo(newMarkerNo);

    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].hasMarkerNo(newMarkerNo));

    // Box::mnMarkerAtomNos_ assigned key value newMarkerNo.
    assert(mnMarkerAtomNos_.count(newMarkerNo) != 0);

    // Box::mnMarkerAtomNos_ assigned mapped value atomNo
    // for key value newMarkerNo.
    assert(std::find(mnMarkerAtomNos_.find(newMarkerNo)->second.cbegin(),
                     mnMarkerAtomNos_.find(newMarkerNo)->second.cend  (),
                     atomNo                                              )
            != mnMarkerAtomNos_.find(newMarkerNo)->second.cend()          );

    assert(isSane());


  }



  auto
  Box::insertQuantumRegionNo(const size_t&       atomNo            ,
                             const unsigned int& newQuantumRegionNo )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo             != 0);
    assert(newQuantumRegionNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::insertQuantumRegionNo(const size_t&       atomNo            ,
    //                            const unsigned int& newQuantumRegionNo )
    //                           noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if((mnQuantumRegionAtomNos_.count(newQuantumRegionNo) != 0      ) &&
       (std::find(mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                          ->second.cbegin()       ,
                  mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                          ->second.cend  ()       ,
                 atomNo                                            )
              == mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                         ->second.cend   ()         )   )
      mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                              ->second.push_back(atomNo);

    if(mnQuantumRegionAtomNos_.count(newQuantumRegionNo) == 0)
      mnQuantumRegionAtomNos_.emplace(newQuantumRegionNo ,
                                      std::vector<size_t>
                                      {
	                                atomNo
                                      }                  );

    quantumRegionNoCount_ = mnQuantumRegionAtomNos_.size();

    nAtoms_[atomNo - 1].insertQuantumRegionNo(newQuantumRegionNo);


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].hasQuantumRegionNo(newQuantumRegionNo));

    // Box::mnQuantumRegionAtomNos_ assigned key value newQuantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(newQuantumRegionNo) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned mapped value atomNo
    // for key value newQuantumRegionNo.
    assert(std::find(mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                             ->second.cbegin(),
                     mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                             ->second.cend  (),
                     atomNo                                          )
                  != mnQuantumRegionAtomNos_.find(newQuantumRegionNo)
                                             ->second.cend());

    assert(isSane());


  }



  auto
  Box::setBoxDimension(const double& newBoxXDimension,
                       const double& newBoxYDimension,
                       const double& newBoxZDimension )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newBoxXDimension));
    assert(!std::isinf(newBoxYDimension));
    assert(!std::isinf(newBoxZDimension));

    assert(newBoxXDimension == newBoxXDimension);
    assert(newBoxYDimension == newBoxYDimension);
    assert(newBoxZDimension == newBoxZDimension);

    assert(newBoxXDimension >= 0.0);
    assert(newBoxYDimension >= 0.0);
    assert(newBoxZDimension >= 0.0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setBoxDimension(const double& newBoxXDimension,
    //                      const double& newBoxYDimension,
    //                      const double& newBoxZDimension )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    boxXDimension_ = newBoxXDimension;
    boxYDimension_ = newBoxYDimension;
    boxZDimension_ = newBoxZDimension;

    boxXOrigin_ = 0.0;
    boxYOrigin_ = 0.0;
    boxZOrigin_ = 0.0;

    //
    // Design by Contract Postcondition
    //

    assert(boxXOrigin_ == 0.0);
    assert(boxYOrigin_ == 0.0);
    assert(boxZOrigin_ == 0.0);

    assert(boxXDimension_ == newBoxXDimension);
    assert(boxYDimension_ == newBoxYDimension);
    assert(boxZDimension_ == newBoxZDimension);

    assert(isSane());


  }



  auto
  Box::setCoordinate(const size_t& atomNo        ,
                     const double& newXCoordinate,
                     const double& newYCoordinate,
                     const double& newZCoordinate )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXCoordinate));
    assert(!std::isinf(newYCoordinate));
    assert(!std::isinf(newZCoordinate));

    assert(newXCoordinate == newXCoordinate);
    assert(newYCoordinate == newYCoordinate);
    assert(newZCoordinate == newZCoordinate);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setCoordinate(const size_t& atomNo        ,
    //                    const double& newXCoordinate,
    //                    const double& newYCoordinate,
    //                    const double& newZCoordinate )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    auto xCoordinate = newXCoordinate;
    auto yCoordinate = newYCoordinate;
    auto zCoordinate = newZCoordinate;

    // In case of periodic box.
    if((boxXDimension_ != 0.0) &&
       (boxYDimension_ != 0.0) &&
       (boxZDimension_ != 0.0)   )
    {
      // Check for box dimensions.
      while(xCoordinate >= boxXOrigin_ + boxXDimension_)
        xCoordinate -= boxXDimension_;
      while(xCoordinate < boxXOrigin_)
        xCoordinate += boxXDimension_;
      while(yCoordinate >= boxYOrigin_ + boxYDimension_)
        yCoordinate -= boxYDimension_;
      while(yCoordinate < boxYOrigin_)
        yCoordinate += boxYDimension_;
      while(zCoordinate >= boxZOrigin_ + boxZDimension_)
        zCoordinate -= boxZDimension_;
      while(zCoordinate < boxZOrigin_)
        zCoordinate += boxZDimension_;
    }

    nAtoms_[atomNo - 1].setCoordinate(xCoordinate,
                                      yCoordinate,
                                      zCoordinate );


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getXCoordinate() == xCoordinate);
    assert(nAtoms_[atomNo - 1].getYCoordinate() == yCoordinate);
    assert(nAtoms_[atomNo - 1].getZCoordinate() == zCoordinate);

    assert(isSane());


  }



  auto
  Box::setGradient(const size_t& atomNo      ,
                   const double& newXGradient,
                   const double& newYGradient,
                   const double& newZGradient )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXGradient));
    assert(!std::isinf(newYGradient));
    assert(!std::isinf(newZGradient));

    assert(newXGradient == newXGradient);
    assert(newYGradient == newYGradient);
    assert(newZGradient == newZGradient);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setGradient(const size_t& atomNo      ,
    //                  const double& newXGradient,
    //                  const double& newYGradient,
    //                  const double& newZGradient )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setGradient(newXGradient,
                                    newYGradient,
                                    newZGradient );


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getXGradient() == newXGradient);
    assert(nAtoms_[atomNo - 1].getYGradient() == newYGradient);
    assert(nAtoms_[atomNo - 1].getZGradient() == newZGradient);

    assert(isSane());


  }



  auto
  Box::setGradientGuess(const size_t& atomNo           ,
                        const double& newXGradientGuess,
                        const double& newYGradientGuess,
                        const double& newZGradientGuess )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXGradientGuess));
    assert(!std::isinf(newYGradientGuess));
    assert(!std::isinf(newZGradientGuess));

    assert(newXGradientGuess == newXGradientGuess);
    assert(newYGradientGuess == newYGradientGuess);
    assert(newZGradientGuess == newZGradientGuess);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setGradientGuess(const size_t& atomNo           ,
    //                       const double& newXGradientGuess,
    //                       const double& newYGradientGuess,
    //                       const double& newZGradientGuess )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setGradientGuess(newXGradientGuess,
                                         newYGradientGuess,
                                         newZGradientGuess );


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getXGradientGuess() == newXGradientGuess);
    assert(nAtoms_[atomNo - 1].getYGradientGuess() == newYGradientGuess);
    assert(nAtoms_[atomNo - 1].getZGradientGuess() == newZGradientGuess);

    assert(isSane());


  }



  auto
  Box::setElectronCount(const size_t&       atomNo          ,
                        const unsigned int& newElectronCount )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setElectronCount(const size_t&       atomNo          ,
    //                       const unsigned int& newElectronCount )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setElectronCount(newElectronCount);


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getElectronCount() == newElectronCount);

    assert(isSane());


  }



  auto
  Box::setMdSymbol(const size_t&      atomNo     ,
                   const std::string& newMdSymbol )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    assert(!newMdSymbol.empty());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setMdSymbol(const size_t&      atomNo     ,
    //                  const std::string& newMdSymbol )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setMdSymbol(newMdSymbol);


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getMdSymbol() == newMdSymbol);

    assert(isSane());


  }



  auto
  Box::setPartialCharge(const size_t& atomNo          ,
                        const double& newPartialCharge )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newPartialCharge));
    assert(newPartialCharge == newPartialCharge);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setPartialCharge(const size_t& atomNo          ,
    //                       const double& newPartialCharge )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setPartialCharge(newPartialCharge);


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getPartialCharge() == newPartialCharge);

    assert(isSane());


  }



  auto
  Box::setPartialChargeGuess(const size_t& atomNo               ,
                             const double& newPartialChargeGuess )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newPartialChargeGuess));
    assert(newPartialChargeGuess == newPartialChargeGuess);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setPartialChargeGuess(const size_t& atomNo               ,
    //                            const double& newPartialChargeGuess )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nAtoms_[atomNo - 1].setPartialChargeGuess(newPartialChargeGuess);


    //
    // Design by Contract Postcondition
    //

    assert(nAtoms_[atomNo - 1].getPartialChargeGuess()
            == newPartialChargeGuess                  );

    assert(isSane());


  }



  auto
  Box::setThermostat(const double& newTemperature    ,
                     const double& newDampeningFactor )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newTemperature));
    assert(newTemperature == newTemperature);
    assert(newTemperature >  0.0           );

    assert(!std::isinf(newDampeningFactor));
    assert(newDampeningFactor == newDampeningFactor);
    assert(newDampeningFactor >= 1.0               );

    assert(newDampeningFactor >= timeFrame_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setThermostat(const double& newTemperature    ,
    //                    const double& newDampeningFactor )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    temperature_     = newTemperature    ;
    dampeningFactor_ = newDampeningFactor;
    

    //
    // Design by Contract Postcondition
    //

    assert(temperature_     == newTemperature    );
    assert(dampeningFactor_ == newDampeningFactor);

    assert(isSane());


  }



  auto
  Box::setTimeFrame(const double& newTimeFrame)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newTimeFrame));
    assert(newTimeFrame == newTimeFrame);

    assert(newTimeFrame <= dampeningFactor_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setTimeFrame(const auto& newTimeFrame)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    timeFrame_ = newTimeFrame;


    //
    // Design by Contract Postcondition
    //

    assert(timeFrame_ == newTimeFrame);

    assert(isSane());


  }



  auto
  Box::setTotalEPot(const double& newTotalEPot)noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newTotalEPot));
    assert(newTotalEPot == newTotalEPot);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::setTotalEPot(const double& newTotalEPot)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    totalEPot_ = newTotalEPot;


    //
    // Design by Contract Postcondition
    //

    assert(totalEPot_ == newTotalEPot);

    assert(isSane());


  }



  auto
  Box::hasAtom(const size_t& atomNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //
    assert(atomNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasAtom(const size_t& atomNo)noexcept
    // ->bool
    //

    return (std::find(nAtomNos_.cbegin(),
                      nAtomNos_.cend  (),
                      atomNo             ) != nAtomNos_.cend());


  }



  auto
  Box::hasDistance(const size_t& atomNo        ,
                   const size_t& neighborAtomNo )const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo         != 0);
    assert(neighborAtomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::nAtomNos_ assigned value neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasDistance(const size_t& atomNo        ,
    //                  const size_t& neighborAtomNo )noexcept
    // ->bool
    //

    return (nAtoms_[atomNo - 1].hasNeighborAtomNo(neighborAtomNo));


  }



  auto
  Box::hasEPot(const unsigned int& quantumRegionNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    // Box::mnQuantumRegionAtomNos_ assigned quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasEPot(const unsigned int& quantumRegionNo)const noexcept
    // ->bool
    //

    return (nEPots_.count(quantumRegionNo));


  }



  auto
  Box::hasEvbEPot(const unsigned int& evbNo          ,
                  const unsigned int& quantumRegionNo )const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo           != 0);
    assert(quantumRegionNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_          .count(evbNo          ) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumrRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasEvbEPot(const unsigned int& evbNo          ,
    //                 const unsigned int& quantumRegionNo )const noexcept
    // ->bool
    //

    return (mnEvbEPots_                    .count(evbNo          ) &&
	    mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo)   );


  }



  auto
  Box::hasEvbNo(const size_t&       atomNo,
                const unsigned int& evbNo  )const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasEvbNo(const size_t& atomNo, const unsigned int& evbNo)noexcept
    // ->bool
    //

    return (nAtoms_[atomNo - 1].hasEvbNo(evbNo));


  }



  auto
  Box::hasEvbNo(const unsigned int& evbNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //
    assert(evbNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasEvbNo(const unsigned int& evbNo)noexcept
    // ->bool
    //

    return (mnEvbAtomNos_.count(evbNo));


  }



  auto
  Box::hasEvbTotalEPot(const unsigned int& evbNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    //Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasEvbTotalEPot(const unsigned int& evbNo)const noexcept
    // ->bool
    //

    return (nEvbTotalEPots_.count(evbNo));


  }



  auto
  Box::hasMarkerNo(const size_t&       atomNo  ,
                   const unsigned int& markerNo )const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo   != 0);
    assert(markerNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnMarkerAtomNos_ assigne dkey value markerNo.
    assert(mnMarkerAtomNos_.count(markerNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasMarkerNo(const size_t& atomNo, const unsigned int& markerNo)
    //                 noexcept
    // ->bool
    //

    return (nAtoms_[atomNo - 1].hasMarkerNo(markerNo));


  }



  auto
  Box::hasMarkerNo(const unsigned int& markerNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //
    assert(markerNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasMarkerNo(const unsigned int& markerNo)noexcept
    // ->bool
    //

    return (mnMarkerAtomNos_.count(markerNo));


  }



  auto
  Box::hasQuantumRegionNo(const size_t&       atomNo         ,
                          const unsigned int& quantumRegionNo )const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo          != 0);
    assert(quantumRegionNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasQuantumRegionNo(const size_t&       atomNo         ,
    //                         const unsigned int& quantumRegionNo )noexcept
    // ->bool
    //

    return (nAtoms_[atomNo - 1].hasQuantumRegionNo(quantumRegionNo));


  }



  auto
  Box::hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
  ->bool
  {


    //
    // Design by Contract Precondition
    //
    assert(quantumRegionNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::hasQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
    // ->bool
    //

    return (mnQuantumRegionAtomNos_.count(quantumRegionNo));


  }



  auto
  Box::getElectronCount(const size_t& atomNo)const noexcept
  ->unsigned int
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getElectronCount(const size_t& atomNo)const noexcept
    // ->unsigned int
    //

    return (nAtoms_[atomNo - 1].getElectronCount());


  }



  auto
  Box::getProtonCount(const size_t& atomNo)const noexcept
  ->unsigned int
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getProtonCount(const size_t& atomNo)const noexcept
    // ->unsigned int
    //

    return (nAtoms_[atomNo - 1].getProtonCount());


  }



  auto
  Box::getAtomCount(void)const noexcept
  ->size_t
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getAtomCount(void)const noexcept
    // ->size_t
    //

   return atomCount_;


  }



  auto
  Box::computeDistance(const size_t& atomNo        ,
                       const size_t& neighborAtomNo )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo         != 0);
    assert(neighborAtomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.end());

    // Box::nAtomNos_ assigned value neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.end());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::computeDistance(const size_t& atomNo, const size_t& neighborAtomNo)
    //                     noexcept
    // ->double
    //

    const auto aX  = nAtoms_[atomNo         - 1].getXCoordinate();
    const auto bX  = nAtoms_[neighborAtomNo - 1].getXCoordinate();

    const auto dx = std::min(std::abs(bX - aX)                            ,
                             std::min(std::abs(bX - aX - boxXDimension_),
                                      std::abs(bX - aX + boxXDimension_) ) );

    const auto aY  = nAtoms_[atomNo         - 1].getYCoordinate();
    const auto bY  = nAtoms_[neighborAtomNo - 1].getYCoordinate();

    const auto dy = std::min(std::abs(bY - aY)                            ,
                             std::min(std::abs(bY - aY - boxYDimension_),
                                      std::abs(bY - aY + boxYDimension_) ) );

    const auto aZ  = nAtoms_[atomNo         - 1].getZCoordinate();
    const auto bZ  = nAtoms_[neighborAtomNo - 1].getZCoordinate();

    const auto dz = std::min(std::abs(bZ - aZ)                            ,
                             std::min(std::abs(bZ - aZ - boxZDimension_),
                                      std::abs(bZ - aZ + boxZDimension_) ) );

    //
    // Design by Contract Postcondition
    //

    assert(((boxXDimension_ == 0.0) &&
            (boxYDimension_ == 0.0) &&
            (boxZDimension_ == 0.0)   ) ||
           ((dx <= boxXDimension_) &&
            (dy <= boxYDimension_) &&
            (dz <= boxZDimension_)    )   );

    assert(dx * dx + dy * dy + dz * dz >= 0.0);

    assert(isSane());


    return sqrt(dx * dx + dy * dy + dz * dz);


  }



  auto
  Box::computeVelocity(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::computeVelocity(const size_t& atomNo)const noexcept
    // ->double
    //

    const auto v = computeVelocityVector(atomNo);

    assert(v.size() == 3);

    return sqrt((v[0] * v[0]) +
                (v[1] * v[1]) +
                (v[2] * v[2])  );


  }



  auto
  Box::getBoxXDimension(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxXDimension(void)const noexcept
    // ->double
    //

    return boxXDimension_;


  }



  auto
  Box::getBoxXOrigin(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxXOrigin(void)const noexcept
    // ->double
    //

    return boxXOrigin_;


  }



  auto
  Box::getBoxYDimension(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxYDimension(void)const noexcept
    // ->double
    //

    return boxYDimension_;


  }



  auto
  Box::getBoxYOrigin(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxYOrigin(void)const noexcept
    // ->double
    //

    return boxYOrigin_;


  }



  auto
  Box::getBoxZDimension(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxZDimension(void)const noexcept
    // ->double
    //

    return boxZDimension_;


  }



  auto
  Box::getBoxZOrigin(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getBoxZOrigin(void)const noexcept
    // ->double
    //

    return boxZOrigin_;


  }



  auto
  Box::getDampeningFactor(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getDampeningFactor(void)const noexcept
    // ->double
    //

    return dampeningFactor_;


  }



  auto
  Box::getDistance(const size_t& atomNo        ,
                   const size_t& neighborAtomNo )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo         != 0);
    assert(neighborAtomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::nAtomNos_ assigned value neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.cend());

    assert(nAtoms_[atomNo - 1].hasNeighborAtomNo(neighborAtomNo));

    assert(nAtoms_[neighborAtomNo - 1].hasNeighborAtomNo(atomNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getDistance(const size_t& atomNo        ,
    //                  const size_t& neighborAtomNo )const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getDistance(neighborAtomNo);


  }



  auto
  Box::getEPot(const unsigned int& quantumRegionNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);

    // Box::nEPots_ assigned key value quantumRegionNo.
    assert(nEPots_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEPot(const unsigned int& quantumRegionNo)const noexcept
    // ->double
    //

    return nEPots_.find(quantumRegionNo)->second;


  }



  auto
  Box::getEvbEPot(const unsigned int& evbNo          ,
                  const unsigned int& quantumRegionNo )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo           != 0);
    assert(quantumRegionNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);

    // Box::mnEvbEPots_ assigned key value evbNo.
    assert(mnEvbEPots_.count(evbNo) != 0);

    // Box::mnEvbEPots_ assigned key value evbNo assigned mapped value
    // with key value quantumRegionNo.
    assert(mnEvbEPots_.find(evbNo)->second.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbEPot(const unsigned int& evbNo          ,
    //                 const unsigned int& quantumRegionNo )const noexcept
    // ->double
    //

    return mnEvbEPots_.find(evbNo)->second.find(quantumRegionNo)->second;


  }



  auto
  Box::getEvbTotalEPot(const unsigned int& evbNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);

    // Box::nEvbTotalEPots_ assigned key value evbNo.
    assert(nEvbTotalEPots_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbTotalEPot(const unsigned int& evbNo)const noexcept
    // ->double
    //

    return nEvbTotalEPots_.find(evbNo)->second;


  }



  auto
  Box::getEvbXGradient(const size_t&       atomNo,
                       const unsigned int& evbNo  )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo)!= 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbXGradient(const size_t&       atomNo,
    //                      const unsigned int& evbNo  )const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getEvbXGradient(evbNo);


  }



  auto
  Box::getEvbYGradient(const size_t&       atomNo,
                       const unsigned int& evbNo  )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo)!= 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbYGradient(const size_t&       atomNo,
    //                      const unsigned int& evbNo  )const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getEvbYGradient(evbNo);


  }



  auto
  Box::getEvbZGradient(const size_t&       atomNo,
                       const unsigned int& evbNo  )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo)!= 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbZGradient(const size_t&       atomNo,
    //                      const unsigned int& evbNo  )const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getEvbZGradient(evbNo);


  }



  auto
  Box::getLastDistance(const size_t& atomNo        ,
                       const size_t& neighborAtomNo )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo         != 0);
    assert(neighborAtomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::nAtomNos_ assigned value neighborAtomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     neighborAtomNo     ) != nAtomNos_.cend());

    assert(nAtoms_[atomNo - 1].hasNeighborAtomNo(neighborAtomNo));

    assert(nAtoms_[neighborAtomNo - 1].hasNeighborAtomNo(atomNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastDistance(const size_t& atomNo        ,
    //                      const size_t& neighborAtomNo )const noexcept
    // ->double
    //

    const auto v = getLastDistanceVector(atomNo        ,
                                         neighborAtomNo );

    return sqrt((v[0] * v[0]) +
                (v[1] * v[1]) +
                (v[2] * v[2])  );


  }



  auto
  Box::getLastXCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastXCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getLastXCoordinate();


  }



  auto
  Box::getLastXDistance(const size_t& atomNo0,
                        const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastXDistance(const size_t& atomNo0,
    //                       const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getLastXCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getLastXCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxXDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxXDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxXDimension_) < std::abs(b - a + boxXDimension_)) ?
        (b - a - boxXDimension_)                                            :
        (b - a + boxXDimension_)                                             ;


  }



  auto
  Box::getLastYCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastYCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getLastYCoordinate();


  }



  auto
  Box::getLastYDistance(const size_t& atomNo0,
                        const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastYDistance(const size_t& atomNo0,
    //                       const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getLastYCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getLastYCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxYDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxYDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxYDimension_) < std::abs(b - a + boxYDimension_)) ?
        (b - a - boxYDimension_)                                            :
        (b - a + boxYDimension_)                                             ;


  }



  auto
  Box::getLastZCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastZCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getLastZCoordinate();


  }



  auto
  Box::getLastZDistance(const size_t& atomNo0,
                        const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastZDistance(const size_t& atomNo0,
    //                       const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getLastZCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getLastZCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxZDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxZDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxZDimension_) < std::abs(b - a + boxZDimension_)) ?
        (b - a - boxZDimension_)                                            :
        (b - a + boxZDimension_)                                             ;


  }



  auto
  Box::getMass(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getMass(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getMass();


  }



  auto
  Box::getPartialCharge(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getPartialCharge(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getPartialCharge();


  }



  auto
  Box::getPartialChargeGuess(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getPartialChargeGuess(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getPartialChargeGuess();


  }



  auto
  Box::getTemperature(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getTemperature(void)const noexcept
    // ->double
    //

    return temperature_;


  }



  auto
  Box::getTimeFrame(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getTimeFrame(void)const noexcept
    // ->double
    //

    return timeFrame_;


  }



  auto
  Box::getTotalEPot(void)const noexcept
  ->double
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getTotalEPot(void)const noexcept
    // ->double
    //

    return totalEPot_;


  }



  auto
  Box::getVdwRadius(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getVdwRadius(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getVdwRadius();


  }



  auto
  Box::getXCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getXCoordinate();


  }



  auto
  Box::getXDistance(const size_t& atomNo0,
                    const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXDistance(const size_t& atomNo0,
    //                   const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getXCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getXCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxXDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxXDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxXDimension_) < std::abs(b - a + boxXDimension_)) ?
        (b - a - boxXDimension_)                                            :
        (b - a + boxXDimension_)                                             ;


  }



  auto
  Box::getXGradient(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXGradient(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getXGradient();


  }



  auto
  Box::getXGradientGuess(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXGradientGuess(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getXGradientGuess();


  }



  auto
  Box::getYCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getYCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getYCoordinate();


  }



  auto
  Box::getYDistance(const size_t& atomNo0,
                    const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getYDistance(const size_t& atomNo0,
    //                   const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getYCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getYCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxYDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxYDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxYDimension_) < std::abs(b - a + boxYDimension_)) ?
        (b - a - boxYDimension_)                                            :
        (b - a + boxYDimension_)                                             ;


  }



  auto
  Box::getYGradient(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getYGradient(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getYGradient();


  }



  auto
  Box::getYGradientGuess(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getYGradientGuess(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getYGradientGuess();


  }



  auto
  Box::getZCoordinate(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getZCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getZCoordinate();


  }



  auto
  Box::getZDistance(const size_t& atomNo0,
                    const size_t& atomNo1 )const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getZDistance(const size_t& atomNo0,
    //                   const size_t& atomNo1 )const noexcept
    // ->double
    //


    const auto a  = nAtoms_[atomNo0 - 1].getZCoordinate();
    const auto b  = nAtoms_[atomNo1 - 1].getZCoordinate();


    return ((std::abs(b - a) < std::abs(b - a - boxZDimension_)) &&
            (std::abs(b - a) < std::abs(b - a + boxZDimension_))   ) ?
      (b - a)                                                        :
      (std::abs(b - a - boxZDimension_) < std::abs(b - a + boxZDimension_)) ?
        (b - a - boxZDimension_)                                            :
        (b - a + boxZDimension_)                                             ;


  }



  auto
  Box::getZGradient(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getZGradient(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getZGradient();


  }



  auto
  Box::getZGradientGuess(const size_t& atomNo)const noexcept
  ->double
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getZGradientGuess(const size_t& atomNo)const noexcept
    // ->double
    //

    return nAtoms_[atomNo - 1].getZGradientGuess();


  }



  auto
  Box::updateCoordinates(void)noexcept
  ->double
  {



    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::updateCoordinates(void)noexcept
    // ->void
    //

    //
    // Compute all velocities and add up current total kinetic energy.
    //

    std::vector< std::vector<double> >nVelocities;

    double currentTotalEKin = 0.0;

    for(const auto& atom:nAtoms_)
    {


      // acceleration [Angstrom per Femtosecond sqaured]
      const auto acceleration = computeAcceleration(atom.getAtomNo());

      // velocity [Angstrom per Femtosecond]      
      auto velocity = computeVelocityVector(atom.getAtomNo());

      // Update velocity by acceleration.
      for(auto i  = 0;
               i != 3;
             ++i      )
        velocity[i] += (acceleration[i] * timeFrame_);

      nVelocities.push_back(velocity);


      //
      // Compute current Ekin.
      //

      // velocity [Angstrom per Femtosecond] 
      const auto velocityScalar =
       sqrt((velocity[0] * velocity[0]) +
            (velocity[1] * velocity[1]) +
            (velocity[2] * velocity[2])  );

      const auto eKin =
      (
        1E-17 * 1E23 *                // Ekin [E -23 Joule]
        atom.getMass() * 1.66053892 * // 1E-27 *
        velocityScalar *              // (1E-10/1E-15) *
        velocityScalar                // * (1E-10/1E-15)
                       /
        2.0          
      );

      currentTotalEKin += eKin;
      
    }


    //
    // Compute current temperature [Kelvin] from current Ekin [E -23 Joule].
    //

    assert(getAtomCount() != 0);

    const auto currentTemperature =
    (
       1E-23                                   *
       currentTotalEKin                        *
       2.0                                     /
       (
         static_cast<double> (getAtomCount()) *
         1.3806488E-23                        *
         3.0            
       )
    );


    //
    // Berendsen Thermostat.
    //

    const auto scalingFactor =
     (currentTemperature > 0.0)   ?
     sqrt
     (
       1.0 + 
       ((timeFrame_       /
         dampeningFactor_   )   *
       ((temperature_      /
         currentTemperature ) -
        1.0                    ) )
      )                           :
      (0.01)                       ;

    // Scale velocities.
    for(auto&& velocity:nVelocities)
      for(auto&& component:velocity)
        component *= scalingFactor;

    // Write out results.
    for(const auto& atomNo:nAtomNos_)
    {

      const auto xCoordinate = nAtoms_[atomNo - 1].getXCoordinate();
      const auto yCoordinate = nAtoms_[atomNo - 1].getYCoordinate();
      const auto zCoordinate = nAtoms_[atomNo - 1].getZCoordinate(); 

      setCoordinate
      (
        atomNo                                                 ,
        xCoordinate + (nVelocities[atomNo - 1][0] * timeFrame_),
        yCoordinate + (nVelocities[atomNo - 1][1] * timeFrame_),
        zCoordinate + (nVelocities[atomNo - 1][2] * timeFrame_)
      );

    }


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


    return currentTemperature;



  }



  auto
  Box::getChemicalSymbol(const size_t& atomNo)const noexcept
  ->std::string
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getChemicalSymbol(const size_t& atomNo)const noexcept
    // ->std::string
    //

    return nAtoms_[atomNo - 1].getChemicalSymbol();


  }



  auto
  Box::getMdSymbol(const size_t& atomNo)const noexcept
  ->std::string
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getMdSymbol(const size_t& atomNo)const noexcept
    // ->std::string
    //

    return nAtoms_[atomNo - 1].getMdSymbol();


  }



  auto
  Box::computeAcceleration(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    // force [E-23 Joule per Angstrom]
    auto acceleration = computeForce(atomNo);

    assert(acceleration.size() == 3);

    for(auto i  = 0;
             i != 3;
           ++i      )
    {
      // force [Joule per meter aka Newton]
      acceleration[i] *= (1e-23 * 1e10);
      // force to acceleration [Meter per second squared]
      acceleration[i] /= (nAtoms_[atomNo - 1].getMass() * 1.66053892e-27);
      // accelertion [Angstrom per Femtosecond squared]
      acceleration[i] *= (1e10 * (1e-15 * 1e-15));
    }


    //
    // Design by Contract Postcondition
    //
    assert(acceleration.size() == 3);


    return acceleration;


  }



  auto
  Box::computeForce(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    // force [E-23 Joule per Angstrom]
    auto force = nAtoms_[atomNo - 1].getGradient();

    if(force[0] == 0.0 &&
       force[1] == 0.0 &&
       force[2] == 0.0   )
      force = nAtoms_[atomNo - 1].getGradientGuess();

    for(auto&& component:force)
      component *= (-1);


    //
    // Design by Contract Postcondition
    //
    assert(force.size() == 3);

    return force;


  }



  auto
  Box::computeVelocityVector(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    const auto bX  = nAtoms_[atomNo - 1].getXCoordinate    ();
    const auto aX  = nAtoms_[atomNo - 1].getLastXCoordinate();

    const auto dx = ((std::abs(bX - aX)) < (std::min(std::abs(bX - aX - boxXDimension_),
                                                     std::abs(bX - aX + boxXDimension_) ))       ) ?
                    (bX - aX                                                                     ) :
                    ((std::abs(bX - aX - boxXDimension_)) < (std::abs(bX - aX + boxXDimension_)) ) ?
                    (bX - aX - boxXDimension_                                                    ) :
                    (bX - aX + boxXDimension_                                                    ) ;

    const auto bY  = nAtoms_[atomNo - 1].getYCoordinate    ();
    const auto aY  = nAtoms_[atomNo - 1].getLastYCoordinate();

    const auto dy = ((std::abs(bY - aY)) < (std::min(std::abs(bY - aY - boxYDimension_),
                                                     std::abs(bY - aY + boxYDimension_) ))       ) ?
                    (bY - aY                                                                     ) :
                    ((std::abs(bY - aY - boxYDimension_)) < (std::abs(bY - aY + boxYDimension_)) ) ?
                    (bY - aY - boxYDimension_                                                    ) :
                    (bY - aY + boxYDimension_                                                    ) ;

    const auto bZ  = nAtoms_[atomNo - 1].getZCoordinate    ();
    const auto aZ  = nAtoms_[atomNo - 1].getLastZCoordinate();

    const auto dz = ((std::abs(bZ - aZ)) < (std::min(std::abs(bZ - aZ - boxZDimension_),
                                                     std::abs(bZ - aZ + boxZDimension_) ))       ) ?
                    (bZ - aZ                                                                     ) :
                    ((std::abs(bZ - aZ - boxZDimension_)) < (std::abs(bZ - aZ + boxZDimension_)) ) ?
                    (bZ - aZ - boxZDimension_                                                    ) :
                    (bZ - aZ + boxZDimension_                                                    ) ;


    return std::vector<double>
    {
      dx / timeFrame_,
      dy / timeFrame_,
      dz / timeFrame_
    };


  }



  auto
  Box::getLastCoordinate(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //


    return std::vector<double>
    {
      nAtoms_[atomNo - 1].getLastXCoordinate(),
      nAtoms_[atomNo - 1].getLastYCoordinate(),
      nAtoms_[atomNo - 1].getLastZCoordinate(),
    };


  }



  auto
  Box::getCoordinate(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getCoordinate(const size_t& atomNo)const noexcept
    // ->double
    //


    return std::vector<double>
    {
      nAtoms_[atomNo - 1].getXCoordinate(),
      nAtoms_[atomNo - 1].getYCoordinate(),
      nAtoms_[atomNo - 1].getZCoordinate(),
    };


  }
      


  auto
  Box::getDistanceVector(const size_t& atomNo0,
                         const size_t& atomNo1 )const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),  
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getDistanceVector(const size_t& atomNo0,
    //                        const size_t& atomNo1 )const noexcept
    // ->std::vector<double>
    //


    return std::vector<double>
    {
      getXDistance(atomNo0,
                   atomNo1 ),
      getYDistance(atomNo0,
                   atomNo1 ),
      getZDistance(atomNo0,
                   atomNo1 ),
    };


  }



  auto
  Box::getEvbGradient(const size_t&       atomNo,
                      const unsigned int& evbNo  )const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);
    assert(evbNo  != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo)!= 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getEvbXGradient(const size_t&       atomNo,
    //                      const unsigned int& evbNo  )const noexcept
    // ->double
    //


    return std::vector<double>
    {
      nAtoms_[atomNo - 1].getEvbXGradient(evbNo),
      nAtoms_[atomNo - 1].getEvbYGradient(evbNo),
      nAtoms_[atomNo - 1].getEvbZGradient(evbNo)
    };


  }



  auto
  Box::getGradient(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXGradient(const size_t& atomNo)const noexcept
    // ->double
    //


    return std::vector<double>
    {
      nAtoms_[atomNo - 1].getXGradient(),
      nAtoms_[atomNo - 1].getYGradient(),
      nAtoms_[atomNo - 1].getZGradient()
    };


  }



  auto
  Box::getGradientGuess(const size_t& atomNo)const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getXGradientGuess(const size_t& atomNo)const noexcept
    // ->double
    //


    return std::vector<double>
    {
      nAtoms_[atomNo - 1].getXGradientGuess(),
      nAtoms_[atomNo - 1].getYGradientGuess(),
      nAtoms_[atomNo - 1].getZGradientGuess()
    };


  }



  auto
  Box::getLastDistanceVector(const size_t& atomNo0,
                             const size_t& atomNo1 )const noexcept
  ->std::vector<double>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo0 != 0);
    assert(atomNo1 != 0);

    // Box::nAtomNos_ assigned values atomNo0 and atomNo1.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo0            ) != nAtomNos_.cend());
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo1            ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getLastDistanceVector(const size_t& atomNo0,
    //                            const size_t& atomNo1 )const noexcept
    // ->std::vector<double>
    //


    return std::vector<double>
    {
      getLastXDistance(atomNo0,
                       atomNo1 ),
      getLastYDistance(atomNo0,
                       atomNo1 ),
      getLastZDistance(atomNo0,
                       atomNo1 ),
    };


  }



  auto
  Box::getNEvbNos(const size_t& atomNo)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNEvbNos(const size_t& atomNo)const noexcept
    // ->std::vector<unsigned int>
    //

    return nAtoms_[atomNo - 1].getNEvbNos();


  }



  auto
  Box::getNEvbNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNEvbNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nEvbNos;

    nEvbNos.reserve(evbNoCount_);

    for(const auto& i:mnEvbAtomNos_)
      nEvbNos.push_back(i.first);


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


    return nEvbNos;


  }



  auto
  Box::getNAtomNos(void)const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNAtomNos(void)const noexcept
    // ->std::vector<size_t>
    //

    return nAtomNos_;


  }



  auto
  Box::getNAtomNosByEvbNo(const unsigned int& evbNo)const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Box::mnEvbAtomNos_ assigned key value evbNo.
    assert(mnEvbAtomNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNAtomNosByEvbNo(const unsigned int& evbNo)const noexcept
    // ->std::vector<size_t>
    //

    return mnEvbAtomNos_.find(evbNo)->second;


  }



  auto
  Box::getNAtomNosByMarkerNo(const unsigned int& markerNo)const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Precondition
    //

    assert(markerNo != 0);

    // Box::mnMarkerAtomNos_ assigned key value markerNo.
    assert(mnMarkerAtomNos_.count(markerNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNAtomNosByMarkerNo(const unsigned int& marker)const noexcept
    // ->std::vector<size_t>
    //

    return mnMarkerAtomNos_.find(markerNo)->second;


  }



  auto
  Box::getNAtomNosByQuantumRegionNo(const unsigned int& quantumRegionNo)
                                   const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    // Box::mnQuantumRegionAtomNos_ assigned key value quantumRegionNo.
    assert(mnQuantumRegionAtomNos_.count(quantumRegionNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNAtomNosByQuantumRegionNo(const unsigned int&
    //                                     quantumRegionNo)const noexcept
    // ->std::vector<size_t>
    //

    return mnQuantumRegionAtomNos_.find(quantumRegionNo)->second;


  }



  auto
  Box::getNMarkerNos(const size_t& atomNo)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNMarkerNos(const size_t& atomNo)const noexcept
    // ->std::vector<unsigned int>
    //

    return nAtoms_[atomNo - 1].getNMarkerNos();


  }



  auto
  Box::getNMarkerNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNMarkerNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nMarkerNos;

    nMarkerNos.reserve(markerNoCount_);

    for(const auto& i:mnMarkerAtomNos_)
      nMarkerNos.push_back(i.first);


    //
    // Design by Contract Postcondition
    //

    assert(isSane());

    return nMarkerNos;


  }



  auto
  Box::getNNeighborAtomNos(const size_t& atomNo,
                           const double& range  )const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(range));
    assert(range == range);

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());

    assert(range > 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNNeighborAtomNos(const size_t& atomNo,
    //                          const double& range  )const noexcept
    // ->std::vector<size_t>
    //

    return nAtoms_[atomNo - 1].getNNeighborAtomNos(range);


  }



  auto
  Box::getNNeighborAtomNos(const size_t& atomNo)const noexcept
  ->std::vector<size_t>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNNeighborAtomNos(const size_t& atomNo)const noexcept
    // ->std::vector<size_t>
    //

    return nAtoms_[atomNo - 1].getNNeighborAtomNos();


  }



  auto
  Box::getNQuantumRegionNos(const size_t& atomNo)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Precondition
    //

    assert(atomNo != 0);

    // Box::nAtomNos_ assigned value atomNo.
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     atomNo             ) != nAtomNos_.cend());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNQuantumRegionNos(const size_t& atomNo)const noexcept
    // ->std::vector<unsigned int>
    //

    return nAtoms_[atomNo - 1].getNQuantumRegionNos();


  }



  auto
  Box::getNQuantumRegionNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Box::getNQuantumRegionNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nQuantumRegionNos;

    nQuantumRegionNos.reserve(quantumRegionNoCount_);

    for(const auto& i:mnQuantumRegionAtomNos_)
      nQuantumRegionNos.push_back(i.first);


    //
    // Design by Contract Postcondition
    //

    assert(isSane());


    return nQuantumRegionNos;


  }



}

