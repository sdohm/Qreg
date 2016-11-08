//
// atom.cpp
//
//  Created on: Aug 11, 2014
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



#include "atom.hpp"



namespace qreg
{



  Atom::Atom(unsigned int newElectronCount,
             unsigned int newProtonCount  ,

             size_t newAtomNo,

             double newLastXCoordinate   ,
             double newLastYCoordinate   ,
             double newLastZCoordinate   ,
             double newMass              ,
             double newPartialCharge     ,
             double newPartialChargeGuess,
             double newVdwRadius         ,
             double newXCoordinate       ,
             double newXGradient         ,
             double newXGradientGuess    ,
             double newYCoordinate       ,
             double newYGradient         ,
             double newYGradientGuess    ,
             double newZCoordinate       ,
             double newZGradient         ,
             double newZGradientGuess    ,

             std::string newChemicalSymbol,
             std::string newMdSymbol       )noexcept:

    electronCount_(newElectronCount),
    protonCount_  (newProtonCount  ),

    atomNo_(newAtomNo),

    lastXCoordinate_   (newLastXCoordinate   ),
    lastYCoordinate_   (newLastYCoordinate   ),
    lastZCoordinate_   (newLastZCoordinate   ),
    mass_              (newMass              ),
    partialCharge_     (newPartialCharge     ),
    partialChargeGuess_(newPartialChargeGuess),
    vdwRadius_         (newVdwRadius         ),
    xCoordinate_       (newXCoordinate       ),
    xGradient_         (newXGradient         ),
    xGradientGuess_    (newXGradientGuess    ),
    yCoordinate_       (newYCoordinate       ),
    yGradient_         (newYGradient         ),
    yGradientGuess_    (newYGradientGuess    ),
    zCoordinate_       (newZCoordinate       ),
    zGradient_         (newZGradient         ),
    zGradientGuess_    (newZGradientGuess    ),

    chemicalSymbol_(newChemicalSymbol),
    mdSymbol_      (newMdSymbol      ),

    nEvbNos_           (),
    nMarkerNos_        (),
    nQuantumRegionNos_ (),
    nEvbXGradients_    (),
    nEvbYGradients_    (),
    nEvbZGradients_    (),
    nNeighborDistances_(),
    nNeighborAtomNos_  (),

    evbNoCount_          (0),
    markerNoCount_       (0),
    neighborAtomNoCount_ (0),
    quantumRegionNoCount_(0)

  {


    //
    // Design by Contract Invariant
    //
    assert(isSane());


  }



  Atom::Atom(const Atom& atom)noexcept
  {


    // Calling copy assignment operator.
    *this = atom;


  }



  Atom::Atom(Atom&& atom)noexcept
  {


    // Calling move assignment operator.
    *this = std::move(atom);


  }



  auto
  Atom::operator=(const Atom& atom)noexcept
  ->Atom&
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(atom.isSane());


    //
    // Perform copy assignment.
    //

    electronCount_ = atom.electronCount_;
    protonCount_   = atom.protonCount_  ;

    atomNo_ = atom.atomNo_;

    lastXCoordinate_    = atom.lastXCoordinate_   ;
    lastYCoordinate_    = atom.lastYCoordinate_   ;
    lastZCoordinate_    = atom.lastZCoordinate_   ;
    mass_               = atom.mass_              ;
    partialCharge_      = atom.partialCharge_     ;
    partialChargeGuess_ = atom.partialChargeGuess_;
    vdwRadius_          = atom.vdwRadius_         ;
    xCoordinate_        = atom.xCoordinate_       ;
    xGradient_          = atom.xGradient_         ;
    xGradientGuess_     = atom.xGradientGuess_    ;
    yCoordinate_        = atom.yCoordinate_       ;
    yGradient_          = atom.yGradient_         ;
    yGradientGuess_     = atom.yGradientGuess_    ;
    zCoordinate_        = atom.zCoordinate_       ;
    zGradient_          = atom.zGradient_         ;
    zGradientGuess_     = atom.zGradientGuess_    ;

    chemicalSymbol_ = atom.chemicalSymbol_;
    mdSymbol_       = atom.mdSymbol_      ;

    nEvbNos_            = atom.nEvbNos_           ;
    nMarkerNos_         = atom.nMarkerNos_        ;
    nQuantumRegionNos_  = atom.nQuantumRegionNos_ ;
    nEvbXGradients_     = atom.nEvbXGradients_    ;
    nEvbYGradients_     = atom.nEvbYGradients_    ;
    nEvbZGradients_     = atom.nEvbZGradients_    ;
    nNeighborDistances_ = atom.nNeighborDistances_;
    nNeighborAtomNos_   = atom.nNeighborAtomNos_  ;

    evbNoCount_           = atom.evbNoCount_          ;
    markerNoCount_        = atom.markerNoCount_       ;
    neighborAtomNoCount_  = atom.neighborAtomNoCount_ ;
    quantumRegionNoCount_ = atom.quantumRegionNoCount_;


    //
    // Design by Contract Invariant
    //

    assert(atom.isSane());
    assert(     isSane());


    return *this;


  }



  auto
  Atom::operator=(Atom&& atom)noexcept
  ->Atom&
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(atom.isSane());


    //
    // Perform move assignment.
    //

    std::swap(     electronCount_,
              atom.electronCount_ );

    std::swap(     protonCount_,
              atom.protonCount_ );

    std::swap(     atomNo_,
              atom.atomNo_ );

    std::swap(     lastXCoordinate_,
              atom.lastXCoordinate_ );

    std::swap(     lastYCoordinate_,
              atom.lastYCoordinate_ );

    std::swap(     lastZCoordinate_,
              atom.lastZCoordinate_ );

    std::swap(     mass_,
              atom.mass_ );

    std::swap(     partialCharge_,
              atom.partialCharge_ );

    std::swap(     partialChargeGuess_,
              atom.partialChargeGuess_ );

    std::swap(     vdwRadius_,
              atom.vdwRadius_ );

    std::swap(     xCoordinate_,
              atom.xCoordinate_ );

    std::swap(     xGradient_,
              atom.xGradient_ );

    std::swap(     xGradientGuess_,
              atom.xGradientGuess_ );

    std::swap(     yCoordinate_,
              atom.yCoordinate_ );

    std::swap(     yGradient_,
              atom.yGradient_ );

    std::swap(     yGradientGuess_,
              atom.yGradientGuess_ );

    std::swap(     zCoordinate_,
              atom.zCoordinate_ );

    std::swap(     zGradient_,
              atom.zGradient_ );

    std::swap(     zGradientGuess_,
              atom.zGradientGuess_ );

    std::swap(     chemicalSymbol_,
              atom.chemicalSymbol_ );

    std::swap(     mdSymbol_,
              atom.mdSymbol_ );

    std::swap(     nEvbNos_,
              atom.nEvbNos_ );

    std::swap(     nMarkerNos_,
              atom.nMarkerNos_ );

    std::swap(     nQuantumRegionNos_,
              atom.nQuantumRegionNos_ );

    std::swap(     nEvbXGradients_,
              atom.nEvbXGradients_ );

    std::swap(     nEvbYGradients_,
              atom.nEvbYGradients_ );

    std::swap(     nEvbZGradients_,
              atom.nEvbZGradients_ );

    std::swap(     nNeighborDistances_,
              atom.nNeighborDistances_ );

    std::swap(     nNeighborAtomNos_,
              atom.nNeighborAtomNos_ );

    std::swap(     evbNoCount_,
              atom.evbNoCount_ );

    std::swap(     markerNoCount_,
              atom.markerNoCount_ );

    std::swap(     neighborAtomNoCount_,
              atom.neighborAtomNoCount_ );

    std::swap(     quantumRegionNoCount_,
              atom.quantumRegionNoCount_ );      

    atom.~Atom();


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    return *this;


  }



  auto
  Atom::operator!=(const Atom& atom)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    return (atomNo_ != atom.atomNo_);


  }



  auto
  Atom::operator<(const Atom& atom)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    return (atomNo_ < atom.atomNo_);


  }



  auto
  Atom::operator<=(const Atom& atom)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    return (atomNo_ <= atom.atomNo_);


  }



  auto
  Atom::operator>(const Atom& atom)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    return (atomNo_ > atom.atomNo_);


  }



  auto
  Atom::operator>=(const Atom& atom)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    return (atomNo_ >= atom.atomNo_);


  }



  Atom::~Atom(void)noexcept
  {



  }



  auto
  Atom::deleteAllDistances(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteAllDistances(void)noexcept
    // ->void
    //

    nNeighborAtomNos_  .clear();
    nNeighborDistances_.clear();

    neighborAtomNoCount_ = nNeighborAtomNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(neighborAtomNoCount_ == 0);

    assert(nNeighborAtomNos_  .empty());
    assert(nNeighborDistances_.empty());

    assert(isSane());


  }



  auto
  Atom::deleteAllEvbGradients(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteAllEvbGradients(void)noexcept
    // ->void
    //

    nEvbNos_.clear();

    nEvbXGradients_.clear();
    nEvbYGradients_.clear();
    nEvbZGradients_.clear();

    evbNoCount_ = nEvbNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(evbNoCount_ == 0);

    assert(nEvbNos_.empty());

    assert(nEvbXGradients_.empty());
    assert(nEvbYGradients_.empty());
    assert(nEvbZGradients_.empty());

    assert(isSane());


  }



  auto
  Atom::deleteAllMarkerNos(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteAllMarkerNos(void)noexcept
    // ->void
    //

    nMarkerNos_.clear();

    markerNoCount_ = nMarkerNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(markerNoCount_ == 0);

    assert(nMarkerNos_.empty());

    assert(isSane());


  }



  auto
  Atom::deleteAllQuantumRegionNos(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteAllQuantumRegionNos(void)noexcept
    // ->void
    //

    nQuantumRegionNos_.clear();

    quantumRegionNoCount_ = nQuantumRegionNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(quantumRegionNoCount_ == 0);

    assert(nQuantumRegionNos_.empty());

    assert(isSane());


  }



  auto
  Atom::deleteDistance(const size_t& neighborAtomNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(neighborAtomNo != 0      );
    assert(neighborAtomNo != atomNo_);

    // Atom::nNeighborDistances_ assigned neighborAtomNo.
    assert(nNeighborDistances_.count(neighborAtomNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteDistance(const size_t& neighborAtomNo)noexcept
    // ->void
    //

    for(std::multimap<double,
                      size_t >::iterator i  = nNeighborAtomNos_.begin();
                                         i != nNeighborAtomNos_.end  ();)
    {

      std::multimap<double,
                    size_t >::iterator eraseI = i++;

      if(eraseI->second == neighborAtomNo)
        nNeighborAtomNos_.erase(eraseI);

    }

    // Delete key value atomNo from Atom::nNeighborDistances_.
    nNeighborDistances_.erase(neighborAtomNo);

    neighborAtomNoCount_ = nNeighborAtomNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(neighborAtomNoCount_ == nNeighborAtomNos_.size());

    // Atom::nNeighborAtomNos not assigned mapped value neighborAtomNo.
    assert(([this          ,
             neighborAtomNo ]()->bool
    {
      for(const auto& i:nNeighborAtomNos_)
        if(i.second == neighborAtomNo)
          return false;
      return true;
    })());

    // Atom::nNeighborDistances not assigned key value neighborAtomNo.
    assert(nNeighborDistances_.count(neighborAtomNo) == 0);

    assert(isSane());


  }



  auto
  Atom::deleteEvbGradient(const unsigned int& evbNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Atom::nEvbNos_ assigned key value evbNo.
    assert(nEvbNos_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteEvbGradient(const unsigned int& evbNo)noexcept
    // ->void
    //

    // Delete key value evbNo from all associated containers.
    nEvbNos_       .erase(evbNo);
    nEvbXGradients_.erase(evbNo);
    nEvbYGradients_.erase(evbNo);
    nEvbZGradients_.erase(evbNo);

    evbNoCount_ = nEvbNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(evbNoCount_ == nEvbNos_.size());

    // Atom::nEvbNos_ not assigned key value evbNo.
    assert(nEvbNos_.count(evbNo) == 0);

    // Atom::nEvbX/Y/ZGradients_ not assigned key value evbNo.
    assert(nEvbXGradients_.count(evbNo) == 0);
    assert(nEvbYGradients_.count(evbNo) == 0);
    assert(nEvbZGradients_.count(evbNo) == 0);

    assert(isSane());


  }



  auto
  Atom::deleteMarkerNo(const unsigned int& markerNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(markerNo != 0);

    // Atom::nMarkerNos_ assigned key value markerNo.
    assert(nMarkerNos_.count(markerNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteMarkerNo(const unsigned int& MarkerNo)noexcept
    // ->void
    //

    nMarkerNos_.erase(markerNo);

    markerNoCount_ = nMarkerNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(markerNoCount_ == nMarkerNos_.size());

    // Atom::nMarkerNos_ not assigned key value markerNo.
    assert(nMarkerNos_.count(markerNo) == 0);

    assert(isSane());


  }



  auto
  Atom::deleteQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNo != 0);

    assert(nQuantumRegionNos_.count(quantumRegionNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::deleteQuantumRegion(const unsigned int& quantumRegion)noexcept
    // ->void
    //

    nQuantumRegionNos_.erase(quantumRegionNo);

    quantumRegionNoCount_ = nQuantumRegionNos_.size();


    //
    // Design by Contract Precondition
    //

    assert(quantumRegionNoCount_ == nQuantumRegionNos_.size());

    // Atom::nQuantumRegionNos_ not assigned key value quantumRegionNo.
    assert(nQuantumRegionNos_.count(quantumRegionNo) == 0);

    assert(isSane());


  }



  auto
  Atom::incrementNegativeCharge(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::incrementNegativeCharge(void)noexcept
    // ->void
    //

    ++electronCount_;


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Atom::incrementPositiveCharge(void)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(electronCount_ != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::incrementPositiveCharge(void)noexcept
    // ->void
    //

    --electronCount_;


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


  }



  auto
  Atom::insertDistance(const size_t& newNeighborAtomNo,
                       const double& newDistance       )noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(newNeighborAtomNo != 0      );
    assert(newNeighborAtomNo != atomNo_);

    assert(!std::isinf(newDistance));
    assert(newDistance == newDistance);
    assert(newDistance >= 0.0        );


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::insertDistance(const size_t& newNeighborAtomNo,
    //                      const double& newDistance       )noexcept
    // ->void
    //

    if(nNeighborDistances_.count(newNeighborAtomNo) != 0)
    {

      for(std::multimap<double,
                        size_t >::iterator i  = nNeighborAtomNos_.begin();
                                           i != nNeighborAtomNos_.end  ();)
      {

        std::multimap<double,
                      size_t >::iterator eraseI = i++;

        if(eraseI->second == newNeighborAtomNo)
          nNeighborAtomNos_.erase(eraseI);

      }

      nNeighborDistances_.erase(newNeighborAtomNo);

      neighborAtomNoCount_ = nNeighborAtomNos_.size();


    }

    nNeighborAtomNos_.emplace(newDistance      ,
                              newNeighborAtomNo );

    nNeighborDistances_.emplace(newNeighborAtomNo,
                                newDistance       );

    neighborAtomNoCount_ = nNeighborAtomNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(neighborAtomNoCount_ == nNeighborAtomNos_.size());

    // Atom::nNeighborAtomNos assigned key value newDistance.
    assert(nNeighborAtomNos_.count(newDistance) != 0);

    // Atom::nNeighborAtomNos assigned mapped value newNeighborAtomNo.
    assert(([this             ,
             newNeighborAtomNo,
             newDistance       ]()->bool
    {
      bool mappedValueExists = false;
      for(const auto& i:nNeighborAtomNos_)
	if((i.first  == newDistance      ) &&
	   (i.second == newNeighborAtomNo)   )
	  mappedValueExists = true;
      return mappedValueExists;
    })());

    // Atom::nNeighborDistances assigned key value newNeighborAtomNo.
    assert(nNeighborDistances_.count(newNeighborAtomNo) != 0);

    // Atom::nNeighborDistances assigned mapped value newDistance.
    assert(nNeighborDistances_.find(newNeighborAtomNo)->second == newDistance);

    assert(isSane());


  }



  auto
  Atom::insertEvbGradient(const unsigned int& newEvbNo       ,
                          const double&       newEvbXGradient,
                          const double&       newEvbYGradient,
                          const double&       newEvbZGradient )noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(newEvbNo != 0);

    assert(!std::isinf(newEvbXGradient));
    assert(!std::isinf(newEvbYGradient));
    assert(!std::isinf(newEvbZGradient));

    assert(newEvbXGradient == newEvbXGradient);
    assert(newEvbYGradient == newEvbYGradient);
    assert(newEvbZGradient == newEvbZGradient);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::insertEvbGradient(const unsigned int& newEvbNo ,
    //                         const double&       newEvbXGradient,
    //                         const double&       newEvbYGradient,
    //                         const double&       newEvbZGradient )noexcept
    // ->void
    //

    if(nEvbNos_.count(newEvbNo) != 0)
    {
      nEvbXGradients_.erase(newEvbNo);
      nEvbYGradients_.erase(newEvbNo);
      nEvbZGradients_.erase(newEvbNo);
    }

    nEvbNos_.emplace(newEvbNo);

    nEvbXGradients_.emplace(newEvbNo       ,
                            newEvbXGradient );

    nEvbYGradients_.emplace(newEvbNo       ,
                            newEvbYGradient );

    nEvbZGradients_.emplace(newEvbNo       ,
                            newEvbZGradient );

    evbNoCount_ = nEvbNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(evbNoCount_ == nEvbNos_.size());

    // Atom::nEvbNos assigned key value newEvbNo.
    assert(nEvbNos_.count(newEvbNo) != 0);

    // Atom::nEvbX/Y/ZGradients assigned key value newEvbNo.
    assert(nEvbXGradients_.count(newEvbNo) != 0);
    assert(nEvbYGradients_.count(newEvbNo) != 0);
    assert(nEvbZGradients_.count(newEvbNo) != 0);

    // Atom::nEvbX/Y/ZGradients assigned mapped value newEvbZGradient.
    assert(nEvbXGradients_.find(newEvbNo)->second == newEvbXGradient);
    assert(nEvbYGradients_.find(newEvbNo)->second == newEvbYGradient);
    assert(nEvbZGradients_.find(newEvbNo)->second == newEvbZGradient);

    assert(isSane());


  }



  auto
  Atom::insertMarkerNo(const unsigned int& newMarkerNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(newMarkerNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::insertMarkerNo(const unsigned int& newMarkerNo)noexcept
    // ->void
    //

    nMarkerNos_.emplace(newMarkerNo);

    markerNoCount_ = nMarkerNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(markerNoCount_ == nMarkerNos_.size());

    // Atom::nMarkerNos assigned key value newMarkerNo.
    assert(nMarkerNos_.count(newMarkerNo) != 0);

    assert(isSane());


  }



  auto
  Atom::insertQuantumRegionNo(const unsigned int& newQuantumRegion)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(newQuantumRegion != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::insertQuantumRegion(const unsigned int& newQuantumRegion)noexcept
    // ->void
    //

    nQuantumRegionNos_.emplace(newQuantumRegion);

    quantumRegionNoCount_ = nQuantumRegionNos_.size();


    //
    // Design by Contract Postcondition
    //

    assert(quantumRegionNoCount_ == nQuantumRegionNos_.size());

    // Atom::nQuantumRegionNos assigned key value newQuantumRegion.
    assert(nQuantumRegionNos_.count(newQuantumRegion) != 0);

    assert(isSane());


  }



  auto
  Atom::setAtomNo(const size_t& newAtomNo)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(newAtomNo != 0);

    // Atom::nNeighborDistances_ not assigned key value newAtomNo.
    assert(nNeighborDistances_.count(newAtomNo) == 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setId(const size_t& newAtomNo)noexcept
    // ->void
    //

    atomNo_ = newAtomNo;


    //
    // Design by Contract Postcondition
    //

    assert(atomNo_ == newAtomNo);

    assert(isSane());


  }



  auto
  Atom::setCoordinate(const double& newXCoordinate,
                      const double& newYCoordinate,
                      const double& newZCoordinate )noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXCoordinate));
    assert(!std::isinf(newYCoordinate));
    assert(!std::isinf(newZCoordinate));

    assert(newXCoordinate == newXCoordinate);
    assert(newYCoordinate == newYCoordinate);
    assert(newZCoordinate == newZCoordinate);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setCoordinate(const double& newXCoordinate,
    //                     const double& newYCoordinate,
    //                     const double& newZCoordinate )noexcept
    // ->void
    //

    lastXCoordinate_ = xCoordinate_;
    lastYCoordinate_ = yCoordinate_;
    lastZCoordinate_ = zCoordinate_;

    xCoordinate_ = newXCoordinate;
    yCoordinate_ = newYCoordinate;
    zCoordinate_ = newZCoordinate;


    //
    // Design by Contract Postcondition
    //

    assert(xCoordinate_ == newXCoordinate);
    assert(yCoordinate_ == newYCoordinate);
    assert(zCoordinate_ == newZCoordinate);

    assert(isSane());


  }



  auto
  Atom::setElectronCount(const unsigned int& newElectronCount)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setElectronCount(const unsigned int& newElectronCount)noexcept
    // ->void
    //

    electronCount_ = newElectronCount;


    //
    // Design by Contract Postcondition
    //

    assert(electronCount_ == newElectronCount);

    assert(isSane());


  }



  auto
  Atom::setGradient(const double& newXGradient,
                    const double& newYGradient,
                    const double& newZGradient )noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXGradient));
    assert(!std::isinf(newYGradient));
    assert(!std::isinf(newZGradient));

    assert(newXGradient == newXGradient);
    assert(newYGradient == newYGradient);
    assert(newZGradient == newZGradient);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setGradient(const double& newXGradient,
    //                   const double& newYGradient,
    //                   const double& newZGradient )noexcept
    // ->void
    //

    xGradient_ = newXGradient;
    yGradient_ = newYGradient;
    zGradient_ = newZGradient;


    //
    // Design by Contract Postcondition
    //

    assert(xGradient_ == newXGradient);
    assert(yGradient_ == newYGradient);
    assert(zGradient_ == newZGradient);

    assert(isSane());


  }



  auto
  Atom::setGradientGuess(const double& newXGradientGuess,
                         const double& newYGradientGuess,
                         const double& newZGradientGuess )noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newXGradientGuess));
    assert(!std::isinf(newYGradientGuess));
    assert(!std::isinf(newZGradientGuess));

    assert(newXGradientGuess == newXGradientGuess);
    assert(newYGradientGuess == newYGradientGuess);
    assert(newZGradientGuess == newZGradientGuess);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setGradientGuess(const double& newXGradientGuess,
    //                        const double& newYGradientGuess,
    //                        const double& newZGradientGuess )noexcept
    // ->void
    //

    xGradientGuess_ = newXGradientGuess;
    yGradientGuess_ = newYGradientGuess;
    zGradientGuess_ = newZGradientGuess;


    //
    // Design by Contract Postcondition
    //

    assert(xGradientGuess_ == newXGradientGuess);
    assert(yGradientGuess_ == newYGradientGuess);
    assert(zGradientGuess_ == newZGradientGuess);

    assert(isSane());


  }



  auto
  Atom::setMdSymbol(const std::string& newMdSymbol)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(!newMdSymbol.empty());


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setMdSymbol(const std::string& newMdSymbol)noexcept
    // ->void
    //

    mdSymbol_ = newMdSymbol;


    //
    // Design by Contract Postcondition
    //

    assert(mdSymbol_ == newMdSymbol);

    assert(isSane());


  }



  auto
  Atom::setPartialCharge(const double& newPartialCharge)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newPartialCharge));
    assert(newPartialCharge == newPartialCharge);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setPartialCharge(const double& newPartialCharge)noexcept
    // ->void
    //

    partialCharge_ = newPartialCharge;


    //
    // Design by Contract Postcondition
    //

    assert(partialCharge_ == newPartialCharge);

    assert(isSane());


  }



  auto
  Atom::setPartialChargeGuess(const double& newPartialChargeGuess)noexcept
  ->void
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(newPartialChargeGuess));
    assert(newPartialChargeGuess == newPartialChargeGuess);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::setPartialChargeGuess(const double& newPartialChargeGuess)noexcept
    // ->void
    //

    partialChargeGuess_ = newPartialChargeGuess;


    //
    // Design by Contract Postcondition
    //

    assert(partialChargeGuess_ == newPartialChargeGuess);

    assert(isSane());


  }



  auto
  Atom::hasEvbNo(const unsigned int& evbNo)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


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
    // Atom::hasEvbNo(const unsigned int& evbNo)const noexcept
    // ->bool
    //
    
    return nEvbNos_.count(evbNo);


  }



  auto
  Atom::hasMarkerNo(const unsigned int& markerNo)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


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
    // Atom::hasMarkerNo(const unsigned int& MarkerNo)const noexcept
    // ->bool
    //
    
    return nMarkerNos_.count(markerNo);


  }



  auto
  Atom::hasNeighborAtomNo(const size_t& neighborAtomNo)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(neighborAtomNo != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::hasneighborAtomNo(const size_t& neighborAtomNo)const noexcept
    // ->bool
    //
    
    return nNeighborDistances_.count(neighborAtomNo);


  }



  auto
  Atom::hasQuantumRegionNo(const unsigned int& quantumRegion)const noexcept
  ->bool
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //
    assert(quantumRegion != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::hasQuantumRegion(const unsigned int& quantumRegion)const noexcept
    // ->bool
    //
    
    return nQuantumRegionNos_.count(quantumRegion);


  }



  auto
  Atom::getElectronCount(void)const noexcept
  ->unsigned int
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getElectronCount(void)const noexcept
    // ->unsigned int
    //
    
    return electronCount_;


  }



  auto
  Atom::getProtonCount(void)const noexcept
  ->unsigned int
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getProtonCount(void)const noexcept
    // ->unsigned int
    //
    
    return protonCount_;


  }



  auto
  Atom::getAtomNo(void)const noexcept
  ->size_t
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getId(void)const noexcept
    // ->unsigned int
    //
    
    return atomNo_;


  }



  auto
  Atom::getDistance(const size_t& neighborAtomNo)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(neighborAtomNo != atomNo_);
    assert(neighborAtomNo != 0      );

    // Atom::nNeighborDistances assigned key value neighborAtomNo.
    assert(nNeighborDistances_.count(neighborAtomNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getDistance(const size_t& neighborAtomNo)const noexcept
    // ->double
    //
    
    return nNeighborDistances_.find(neighborAtomNo)->second;


  }



  auto
  Atom::getEvbXGradient(const unsigned int& evbNo)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Atom::nEvbXGradients assigned key value evbNo.
    assert(nEvbXGradients_.count(evbNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getEvbXGradient(const unsigned int& evbNo)const noexcept
    // ->double
    //
    
    return nEvbXGradients_.find(evbNo)->second;


  }



  auto
  Atom::getEvbYGradient(const unsigned int& evbNo)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Atom::nEvbYGradients assigned key value evbNo.
    assert(nEvbYGradients_.count(evbNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getEvbYGradient(const unsigned int& evbNo)const noexcept
    // ->double
    //
    
    return nEvbYGradients_.find(evbNo)->second;


  }



  auto
  Atom::getEvbZGradient(const unsigned int& evbNo)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Atom::nEvbZGradients assigned key value evbNo.
    assert(nEvbZGradients_.count(evbNo));


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getEvbZGradient(const unsigned int& evbNo)const noexcept
    // ->double
    //
    
    return nEvbZGradients_.find(evbNo)->second;


  }



  auto
  Atom::getLastXCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getLastXCoordinate(void)const noexcept
    // ->double
    //

    return lastXCoordinate_;


  }



  auto
  Atom::getLastYCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getLastYCoordinate(void)const noexcept
    // ->double
    //

    return lastYCoordinate_;


  }



  auto 
  Atom::getLastZCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getLastZCoordinate(void)const noexcept
    // ->double
    //

    return lastZCoordinate_;


  }



  auto
  Atom::getMass(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getMass(void)const noexcept
    // ->double
    //

    return mass_;


  }



  auto
  Atom::getPartialCharge(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getPartialCharge(void)const noexcept
    // ->double
    //
    
    return partialCharge_;


  }



  auto
  Atom::getPartialChargeGuess(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getPartialChargeGuess(void)const noexcept
    // ->double
    //
    
    return partialChargeGuess_;


  }



  auto
  Atom::getVdwRadius(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getVdwRadius(void)const noexcept
    // ->double
    //

    return vdwRadius_;


  }



  auto
  Atom::getXCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getXCoordinate(void)const noexcept
    // ->double
    //
    
    return xCoordinate_;


  }



  auto
  Atom::getXGradient(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getXGradient(void)const noexcept
    // ->double
    //
    
    return xGradient_;


  }



  auto
  Atom::getXGradientGuess(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getXGradientGuess(void)const noexcept
    // ->double
    //
    
    return xGradientGuess_;


  }



  auto
  Atom::getYCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getYCoordinate(void)const noexcept
    // ->double
    //
    
    return yCoordinate_;


  }



  auto
  Atom::getYGradient(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getYGradient(void)const noexcept
    // ->double
    //
    
    return yGradient_;


  }



  auto
  Atom::getYGradientGuess(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getYGradientGuess(void)const noexcept
    // ->double
    //
    
    return yGradientGuess_;


  }



  auto
  Atom::getZCoordinate(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getZCoordinate(void)const noexcept
    // ->double
    //
    
    return zCoordinate_;


  }



  auto
  Atom::getZGradient(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getZGradient(void)const noexcept
    // ->double
    //
    
    return zGradient_;


  }



  auto
  Atom::getZGradientGuess(void)const noexcept
  ->double
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getZGradientGuess(void)const noexcept
    // ->double
    //
    
    return zGradientGuess_;


  }



  auto
  Atom::getChemicalSymbol(void)const noexcept
  ->std::string
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getChemicalSymbol(void)const noexcept
    // ->std::string
    //
    
    return chemicalSymbol_;


  }



  auto
  Atom::getMdSymbol(void)const noexcept
  ->std::string
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getMdSymbol(void)const noexcept
    // ->std::string
    //
    
    return mdSymbol_;


  }



  auto
  Atom::getCoordinate(void)const noexcept
  ->std::vector<double>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getCoordinate(void)const noexcept
    // ->std::vector<double>
    //
    std::vector<double> coordinate = {xCoordinate_,
                                      yCoordinate_,
                                      zCoordinate_ };


    //
    // Design by Contract Postcondition
    //

    assert(coordinate.size() == 3);

    assert(coordinate[0] == xCoordinate_);
    assert(coordinate[1] == yCoordinate_);
    assert(coordinate[2] == zCoordinate_);

    assert(isSane());


    return coordinate;


  }



  auto
  Atom::getEvbGradient(const unsigned int& evbNo)const noexcept
  ->std::vector<double>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(evbNo != 0);

    // Atom::nEvbX/Y/ZGradients assigned key value evbNo.
    assert(nEvbXGradients_.count(evbNo) != 0);
    assert(nEvbYGradients_.count(evbNo) != 0);
    assert(nEvbZGradients_.count(evbNo) != 0);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getEvbGradient(const unsigned int& evbNo)const noexcept
    // ->std::vector<double>
    //
    std::vector<double> evbGradient = {nEvbXGradients_.find(evbNo)->second,
	                               nEvbYGradients_.find(evbNo)->second,
	                               nEvbZGradients_.find(evbNo)->second };


    //
    // Design by Contract Postcondition
    //

    assert(evbGradient.size() == 3);

    // EvbGradient assigned appropriate nEvbX/Y/ZGradients_ values.
    assert(evbGradient[0] == nEvbXGradients_.find(evbNo)->second);
    assert(evbGradient[1] == nEvbYGradients_.find(evbNo)->second);
    assert(evbGradient[2] == nEvbZGradients_.find(evbNo)->second);

    assert(isSane());


    return evbGradient;


  }



  auto
  Atom::getGradient(void)const noexcept
  ->std::vector<double>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getGradient(void)const noexcept
    // ->std::vector<double>
    //
    std::vector<double> gradient = {xGradient_,
                                    yGradient_,
                                    zGradient_ };


    //
    // Design by Contract Postcondition
    //

    assert(gradient.size() == 3);

    // Gradient assigned appropriate x/y/zGradient_ values.
    assert(gradient[0] == xGradient_);
    assert(gradient[1] == yGradient_);
    assert(gradient[2] == zGradient_);

    assert(isSane());


    return gradient;


  }



  auto
  Atom::getGradientGuess(void)const noexcept
  ->std::vector<double>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getGradientGuess(void)const noexcept
    // ->std::vector<double>
    //
    std::vector<double> gradientGuess = {xGradientGuess_,
	                                 yGradientGuess_,
	                                 zGradientGuess_ };


    //
    // Design by Contract Postcondition
    //

    assert(gradientGuess.size() == 3);

    // GradientGuess assigned appropriate x/y/zGradientGuess_ values.
    assert(gradientGuess[0] == xGradientGuess_);
    assert(gradientGuess[1] == yGradientGuess_);
    assert(gradientGuess[2] == zGradientGuess_);

    assert(isSane());


    return gradientGuess;


  }



  auto
  Atom::getNEvbNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getNEvbNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nEvbNos;

    nEvbNos.reserve(evbNoCount_);

    for(const auto& i:nEvbNos_)
      nEvbNos.push_back(i);


    //
    // Design by Contract Postcondition
    //

    assert(nEvbNos.size() == evbNoCount_);

    // NEvbNos assigned Atom::nEvbNos_ key values.
    assert(([this   ,
             nEvbNos ]()->bool
    {
      for(const auto& i:nEvbNos)
	if(nEvbNos_.count(i) == 0)
	  return false;
      return true;
    })());

    assert(isSane());


    return nEvbNos;


  }



  auto
  Atom::getNMarkerNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getNMarkerNoNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nMarkerNos;

    nMarkerNos.reserve(markerNoCount_);

    for(const auto& i:nMarkerNos_)
      nMarkerNos.push_back(i);


    //
    // Design by Contract Postcondition
    //

    assert(nMarkerNos.size() == markerNoCount_);

    // NMarkerNos assigned Atom::nMarkerNos_ key values.
    assert(([this      ,
             nMarkerNos ]()->bool
    {
      for(const auto& i:nMarkerNos)
	if(nMarkerNos_.count(i) == 0)
	  return false;
      return true;
    })());

    assert(isSane());


    return nMarkerNos;


  }



  auto
  Atom::getNNeighborAtomNos(const double& range)const noexcept
  ->std::vector<size_t>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Precondition
    //

    assert(!std::isinf(range));
    assert(range == range);
    assert(range >  0.0  );


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getnNeighborAtomNos(void)const noexcept
    // ->std::vector<size_t>
    //

    std::vector<size_t> nNeighborAtomNos;

    nNeighborAtomNos.reserve(64);

    for(auto&&  i         = nNeighborAtomNos_.cbegin()    ;
               (i        != nNeighborAtomNos_.cend  ()) &&
               (i->first <= range                     )   ;
              ++i                                          )
      nNeighborAtomNos.push_back(i->second);


    //
    // Design by Contract Postcondition
    //
    assert(isSane());


    return nNeighborAtomNos;


  }


  auto
  Atom::getNNeighborAtomNos(void)const noexcept
  ->std::vector<size_t>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getnNeighborAtomNos(void)const noexcept
    // ->std::vector<size_t>
    //

    std::vector<size_t> nNeighborAtomNos;

    nNeighborAtomNos.reserve(neighborAtomNoCount_);

    for(const auto& i:nNeighborAtomNos_)
      nNeighborAtomNos.push_back(i.second);


    //
    // Design by Contract Postcondition
    //

    assert(nNeighborAtomNos.size() == neighborAtomNoCount_);

    // NNeighborAtomNos assigned nNeighborAtomNos mapped values.
    assert(([this            ,
             nNeighborAtomNos ]()->bool
    {
      for(const auto& i:nNeighborAtomNos)
	if(nNeighborDistances_.count(i) == 0)
	  return false;
      return true;
    })());

    assert(isSane());


    return nNeighborAtomNos;


  }



  auto
  Atom::getNQuantumRegionNos(void)const noexcept
  ->std::vector<unsigned int>
  {


    std::lock_guard<std::mutex> lock(threadSafety_);


    //
    // Design by Contract Invariant
    //
    assert(isSane());


    //
    // auto
    // Atom::getNQuantumRegionNos(void)const noexcept
    // ->std::vector<unsigned int>
    //

    std::vector<unsigned int> nQuantumRegionNos;

    nQuantumRegionNos.reserve(quantumRegionNoCount_);

    for(const auto& i:nQuantumRegionNos_)
      nQuantumRegionNos.push_back(i);


    //
    // Design by Contract Postcondition
    //

    assert(nQuantumRegionNos.size() == quantumRegionNoCount_);

    // NQuantumRegionNos assigned nQuantumRegionNos values.
    assert(([this             ,
             nQuantumRegionNos ]()->bool
    {
      for(const auto& i:nQuantumRegionNos)
	if(nQuantumRegionNos_.count(i) == 0)
	  return false;
      return true;
    })());

    assert(isSane());


    return nQuantumRegionNos;


  }



}

