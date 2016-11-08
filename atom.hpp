//
// atom.hpp
//
//  Created on: Aug 11,2014
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



#ifndef QREG_ATOMBOX_ATOM_HPP_
#define QREG_ATOMBOX_ATOM_HPP_



#define NDEBUG



#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>



namespace qreg
{



  //
  // Describes an atom. Defaults to hydrogen.
  //
  class Atom final
  {



    public :


      //
      // Defines a hydrogen atom. AtomNo is "1".
      // Group properties are zero or empty.
      //
      Atom(unsigned int newElectronCount = 1,
           unsigned int newProtonCount   = 1,

           size_t newAtomNo = 1,

           double newLastXCoordinate    = 0.0  ,
           double newLastYCoordinate    = 0.0  ,
           double newLastZCoordinate    = 0.0  ,
           double newMass               = 1.008,
           double newPartialCharge      = 0.0  ,
           double newPartialChargeGuess = 0.0  ,
           double newVdwRadius          = 1.20 ,
           double newXCoordinate        = 0.0  ,
           double newXGradient          = 0.0  ,
           double newXGradientGuess     = 0.0  ,
           double newYCoordinate        = 0.0  ,
           double newYGradient          = 0.0  ,
           double newYGradientGuess     = 0.0  ,
           double newZCoordinate        = 0.0  ,
           double newZGradient          = 0.0  ,
           double newZGradientGuess     = 0.0  ,

           std::string newChemicalSymbol = "H",
           std::string newMdSymbol       = "H"  )noexcept;

      //
      // Copy construction enabled.
      //
      Atom(const Atom& atom)noexcept;

      //
      // Move construction enabled.
      //
      Atom(Atom&& atom)noexcept;

      //
      // Copy assignment enabled.
      //
      auto
      operator=(const Atom& atom)noexcept
      ->Atom&;

      //
      // Move assignment enabled.
      //
      auto
      operator=(Atom&& atom)noexcept
      ->Atom&;

      //
      // != operator available.
      //
      auto
      operator!=(const Atom& atom)const noexcept
      ->bool;

      //
      // < operator available.
      //
      auto
      operator<(const Atom& atom)const noexcept
      ->bool;

      //
      // <= operator available.
      //
      auto
      operator<=(const Atom& atom)const noexcept
      ->bool;

      //
      // > operator available.
      //
      auto
      operator>(const Atom& atom)const noexcept
      ->bool;

      //
      // >= operator available.
      //
      auto
      operator>=(const Atom& atom)const noexcept
      ->bool;
      
      //
      // Destructor.
      //
      ~Atom(void)noexcept;


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(Atom& lhs,
           Atom& rhs )
      ->void;


      //
      // Deletes all distances to neighbor atoms.
      // Linear complexity (neighbor atom count).
      //
      auto
      deleteAllDistances(void)noexcept
      ->void;

      //
      // Deletes all EVB gradients.
      // Linear complexity (EVB count).
      //
      auto
      deleteAllEvbGradients(void)noexcept
      ->void;

      //
      // Deletes all markers.
      // Linear complexity (marker count).
      //
      auto
      deleteAllMarkerNos(void)noexcept
      ->void;

      //
      // Deletes all quantum regions.
      // Linear complexity (quantum region count).
      //
      auto
      deleteAllQuantumRegionNos(void)noexcept
      ->void;

      //
      // Deletes one distance to another atom.
      // Logarithmic complexity (neighbor atom count).
      //
      auto
      deleteDistance(const size_t& atomNo)noexcept
      ->void;

      //
      // Deletes one EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      deleteEvbGradient(const unsigned int& evbNo)noexcept
      ->void;

      //
      // Deletes one marker.
      // Logarithmic complexity (marker count).
      //
      auto
      deleteMarkerNo(const unsigned int& markerNo)noexcept
      ->void;

      //
      // Deletes one quantum region.
      // Logarithmic complexity (quantum region count).
      //
      auto
      deleteQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
      ->void;

      //
      // Increases electron count by one.
      // Constant complexity.
      //
      auto
      incrementNegativeCharge(void)noexcept
      ->void;

      //
      // Decreases electron count by one.
      // Constant complexity.
      //
      auto
      incrementPositiveCharge(void)noexcept
      ->void;

      //
      // Inserts new distance to another atom.
      // Logarithmic complexity (neighbor atom count).
      //
      auto
      insertDistance(const size_t& neigborAtomNo,
                     const double& newDistance   )noexcept
      ->void;

      //
      // Inserts new EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      insertEvbGradient (const unsigned int& newEvbNo       ,
                         const double&       newEvbXGradient,
                         const double&       newEvbYGradient,
                         const double&       newEvbZGradient )noexcept
      ->void;

      //
      // Inserts new marker.
      // Logarithmic complexity (marker count).
      //
      auto
      insertMarkerNo(const unsigned int& newMarkerNo)noexcept
      ->void;

      //
      // Inserts new quantum region.
      // Logarithmic complexity (quantum region count).
      //
      auto
      insertQuantumRegionNo(const unsigned int& newQuantumRegionNo)noexcept
      ->void;

      //
      // Sets atom number.
      // Constant complexity.
      //
      auto
      setAtomNo(const size_t& newAtomNo)noexcept
      ->void;

      //
      // Sets coordinate.
      // Constant complexity.
      //
      auto
      setCoordinate(const double& newXCoordinate,
                    const double& newYCoordinate,
                    const double& newZCoordinate )noexcept
      ->void;

      //
      // Sets electron count.
      // Constant complexity.
      //
      auto
      setElectronCount(const unsigned int& newElectronCount)noexcept
      ->void;

      //
      // Sets QM gradient.
      // Constant complexity.
      //
      auto
      setGradient(const double& newXGradient,
                  const double& newYGradient,
                  const double& newZGradient )noexcept
      ->void;

      //
      // Sets MD gradient guess.
      // Constant complexity.
      //
      auto
      setGradientGuess(const double& newXGradientGuess,
                       const double& newYGradientGuess,
                       const double& newZGradientGuess )noexcept
      ->void;

      //
      // Sets MD symbol of an atom.
      // Constant complexity.
      //
      auto
      setMdSymbol(const std::string& newMdSymbol)noexcept
      ->void;

      //
      // Sets QM partial charge.
      // Constant complexity.
      //
      auto
      setPartialCharge(const double& newPartialCharge)noexcept
      ->void;

      //
      // Sets MD partial charge guess.
      // Constant complexity.
      //
      auto
      setPartialChargeGuess(const double& newPartialChargeGuess)noexcept
      ->void;


      //
      // Checks whether inclusion in an EVB is given.
      // Logarithmic complexity (EVB count).
      //
      auto
      hasEvbNo(const unsigned int& evbNo)const noexcept
      ->bool;

      //
      // Checks whether a marker is set.
      // Logarithmic complexity (marker count).
      //
      auto
      hasMarkerNo(const unsigned int& markerNo)const noexcept
      ->bool;

      //
      // Checks whether a neighbor atom exists.
      // Logarithmic complexity (neighbor atom count).
      //
      auto
      hasNeighborAtomNo(const size_t& neighborAtomNo)const noexcept
      ->bool;

      //
      // Checks whether inclusion in a quantum region is given.
      // Logarithmic complexity (quantum region count).
      //
      auto
      hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
      ->bool;


      //
      // Returns electron count.
      // Constant complexity.
      //
      auto
      getElectronCount(void)const noexcept
      ->unsigned int;

      //
      // Returns proton count.
      // Constant complexity.
      //
      auto
      getProtonCount(void)const noexcept
      ->unsigned int;


      //
      // Returns atom number.
      // Constant complexity.
      //
      auto
      getAtomNo(void)const noexcept
      ->size_t;


      //
      // Returns distance to neighbor.
      // Logarithmic complexity (neighbor atom count).
      //
      auto
      getDistance(const size_t& neighborAtomNo)const noexcept
      ->double;

      //
      // Returns x component of the EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbXGradient(const unsigned int& evbNo)const noexcept
      ->double;

      //
      // Returns y component of the EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbYGradient(const unsigned int& evbNo)const noexcept
      ->double;

      //
      // Returns z component of the EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbZGradient(const unsigned int& evbNo)const noexcept
      ->double;

      //
      // Returns x component of the previous coordinate.
      // Constant complexity.
      //
      auto
      getLastXCoordinate(void)const noexcept
      ->double;

      //
      // Returns y component of the previous coordinate.
      // Constant complexity.
      //
      auto
      getLastYCoordinate(void)const noexcept
      ->double;

      //
      // Returns z component of the previous coordinate.
      // Constant complexity.
      //
      auto
      getLastZCoordinate(void)const noexcept
      ->double;

      //
      // Returns mass.
      // Constant complexity.
      //
      auto
      getMass(void)const noexcept
      ->double;

      //
      // Returns QM partial charge.
      // Constant complexity.
      //
      auto
      getPartialCharge(void)const noexcept
      ->double;

      //
      // Returns MD partial charge guess.
      // Constant complexity.
      //
      auto
      getPartialChargeGuess(void)const noexcept
      ->double;

      //
      // Returns van der Waals radius.
      // Constant complexity.
      //
      auto
      getVdwRadius(void)const noexcept
      ->double;

      //
      // Returns x component of the coordinate.
      // Constant complexity.
      //
      auto
      getXCoordinate(void)const noexcept
      ->double;

      //
      // Returns x component of the QM gradient.
      // Constant complexity.
      //
      auto
      getXGradient(void)const noexcept
      ->double;

      //
      // Returns x component a MD the gradient guess.
      // Constant complexity.
      //
      auto
      getXGradientGuess(void)const noexcept
      ->double;

      //
      // Returns y component of the coordinate.
      // Constant complexity.
      //
      auto
      getYCoordinate(void)const noexcept
      ->double;

      //
      // Returns y component of the QM gradient.
      // Constant complexity.
      //
      auto
      getYGradient(void)const noexcept
      ->double;

      //
      // Returns y component of a MD gradient.
      // Constant complexity.
      //
      auto
      getYGradientGuess(void)const noexcept
      ->double;

      //
      // Returns z component of the coordinate.
      // Constant complexity.
      //
      auto
      getZCoordinate(void)const noexcept
      ->double;

      //
      // Returns z component of the QM gradient.
      // Constant complexity.
      //
      auto
      getZGradient(void)const noexcept
      ->double;

      //
      // Returns z component of a MD gradient guess.
      // Constant complexity.
      //
      auto
      getZGradientGuess(void)const noexcept
      ->double;


      //
      // Returns chemical symbol.
      // Constant complexity.
      //
      auto
      getChemicalSymbol(void)const noexcept
      ->std::string;

      //
      // Returns MD symbol.
      // Constant complexity.
      //
      auto
      getMdSymbol(void)const noexcept
      ->std::string;


      //
      // Returns coordinate.
      // Constant complexity.
      //
      auto
      getCoordinate(void)const noexcept
      ->std::vector<double>;

      //
      // Returns EVB gradient.
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbGradient(const unsigned int& evbNo)const noexcept
      ->std::vector<double>;

      //
      // Returns QM gradient.
      // Constant complexity.
      //
      auto
      getGradient(void)const noexcept
      ->std::vector<double>;

      //
      // Returns MD gradient guess.
      // Constant complexity.
      //
      auto
      getGradientGuess(void)const noexcept
      ->std::vector<double>;

      //
      // Returns previous coordinate.
      // Constant complexity.
      //
      auto
      getLastCoordinate(void)const noexcept
      ->std::vector<double>;

      //
      // Returns list of N EVB numbers.
      // Linear complexity (EVB count).
      //
      auto
      getNEvbNos(void)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N markers.
      // Linear complexity (marker count).
      //
      auto
      getNMarkerNos(void)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N neighbor atom numbers sorted by ascending distances
      // within given range.
      // Constant complexity (neighbor atom count) for range chosen within
      // reason. Linear complexity (neighbor atom count) otherwise.
      //
      auto
      getNNeighborAtomNos(const double& range)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N neighbor atom numbers sorted by ascending distances.
      // Linear complexity (neighbor atom count).
      //
      auto
      getNNeighborAtomNos(void)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N quantum regions.
      // Linear complexity (quantum region count).
      //
      auto
      getNQuantumRegionNos(void)const noexcept
      ->std::vector<unsigned int>;



    private :


      //
      // Checks sanity i.e. Design by Contract compliance.
      //
      auto
      isSane(void)const noexcept
      ->bool;


      //
      // Locking and unlocking the Atom for exception and thread safety.
      //
      mutable std::mutex threadSafety_;


      //
      // Stores electron count. Default "1". Assign non-negative value.
      //
      unsigned int electronCount_;

      //
      // Stores proton count. Default "1".
      // Assign non-negative non-zero value.
      //
      unsigned int protonCount_;


      //
      // Stores atom number. Default "1". Assign non-zero value.
      //
      size_t atomNo_;


      //
      // Stores x component of the previous coordinate in angstrom. Default "0".
      //
      double lastXCoordinate_;

      //
      // Stores y component of the previous coordinate in angstrom. Default "0".
      //
      double lastYCoordinate_;

      //
      // Stores z component of the previous coordinate in angstrom. Default "0".
      //
      double lastZCoordinate_;

      //
      // Stores mass in u. Default "1.008".
      //
      double mass_;

      //
      // Stores partial charge. Default "0".
      //
      double partialCharge_;

      //
      // Stores MD guess at partial charge. Default "0".
      //
      double partialChargeGuess_;

      //
      // Stores van der Waals radius. Default "1.20".
      //
      double vdwRadius_;

      //
      // Stores x component of the coordinate in angstrom. Default "0".
      //
      double xCoordinate_;

      //
      // Stores x component of the QM gradient in E-23 J per angstrom.
      // Default "0".
      //
      double xGradient_;

      //
      // Stores x component of a MD gradient guess in E-23 J per angstrom.
      // Default "0".
      //
      double xGradientGuess_;

      //
      // Stores y component of the coordinate in angstrom. Default "0".
      //
      double yCoordinate_;

      //
      // Stores y component of the QM gradient in E-23 J per angstrom.
      // Default "0".
      //
      double yGradient_;

      //
      // Stores y component of a MD gradient guess in E-23 J per angstrom.
      // Default "0".
      //
      double yGradientGuess_;

      //
      // Stores z component of the coordinate in angstrom. Default "0".
      //
      double zCoordinate_;

      //
      // Stores z component of the QM gradient in E-23 J per angstrom.
      // Default "0".
      //
      double zGradient_;

      //
      // Stores z component of a MD gradient guess in E-23 J per angstrom.
      // Default "0".
      //
      double zGradientGuess_;


      //
      // Stores chemical symbol such as "Fe" or "S". Default "H".
      // Assign non-empty string value. Default "H".
      //
      std::string chemicalSymbol_;

      //
      // Stores MD symbol such as "HB" or "H1". default "H".
      // Assign non-empty string value. Default "H".
      //
      std::string mdSymbol_;


      //
      // Stores N EVB numbers. Empty on creation. Assign non-zero values.
      //
      std::set<unsigned int> nEvbNos_;

      //
      // Stores N markers used to partition quantum regions.
      // Empty on creation. Assign non-zero values.
      //
      std::set<unsigned int> nMarkerNos_;

      //
      // Stores N identifiers symbolizing inclusion into quantum regions.
      // Empty on creation. Assign non-zero values.
      //
      std::set<unsigned int> nQuantumRegionNos_;


      //
      // Stores N x components of the EVB gradients in E-23 J per angstrom.
      // Empty on creation.
      //
      std::map<unsigned int,
               double       > nEvbXGradients_;

      //
      // Stores N y components of the EVB gradients in E-23 J per angstrom.
      // Empty on creation.
      //
      std::map<unsigned int,
               double       > nEvbYGradients_;

      //
      // Stores N z components of the EVB gradients in E-23 J per angstrom.
      // Empty on creation.
      //
      std::map<unsigned int,
               double       > nEvbZGradients_;

      //
      // Stores N distances to other atoms sorted by their atom numbers.
      // Empty on creation.
      //
      std::map<size_t,
               double > nNeighborDistances_;


      //
      // Stores N numbers of neighbor atoms sorted by proximity.
      // Empty on creation.
      //
      std::multimap<double,
                    size_t > nNeighborAtomNos_;


      //
      // Stores EVB count. Default "0".
      //
      decltype(nEvbNos_.size()) evbNoCount_;

      //
      // Stores marker count. Default "0".
      //
      decltype(nMarkerNos_.size()) markerNoCount_;

      //
      // Stores neighbor atom count. Default "0".
      //
      decltype(nNeighborAtomNos_.size()) neighborAtomNoCount_;

      //
      //  Stores quantum region count. Default "0".
      //
      decltype(nQuantumRegionNos_.size()) quantumRegionNoCount_;



  };



  inline auto
  Atom::isSane(void)const noexcept
  ->bool
  {


    //
    // Design by Contract Invariant
    //

    assert(atomNo_      != 0);
    assert(protonCount_ != 0);

    assert(!std::isinf(mass_));
    assert(mass_ == mass_);
    assert(mass_ >  0.0  );

    assert(!std::isinf(vdwRadius_));
    assert(vdwRadius_ == vdwRadius_);
    assert(vdwRadius_ >  0.0       );

    assert(!std::isinf(partialCharge_     ));
    assert(!std::isinf(partialChargeGuess_));

    assert(partialCharge_      == partialCharge_     );
    assert(partialChargeGuess_ == partialChargeGuess_);

    assert(!std::isinf(xCoordinate_));
    assert(!std::isinf(yCoordinate_));
    assert(!std::isinf(zCoordinate_));

    assert(xCoordinate_ == xCoordinate_);
    assert(yCoordinate_ == yCoordinate_);
    assert(zCoordinate_ == zCoordinate_);

    assert(!std::isinf(lastXCoordinate_));
    assert(!std::isinf(lastYCoordinate_));
    assert(!std::isinf(lastZCoordinate_));

    assert(lastXCoordinate_ == lastXCoordinate_);
    assert(lastYCoordinate_ == lastYCoordinate_);
    assert(lastZCoordinate_ == lastZCoordinate_);

    assert(!std::isinf(xGradient_));
    assert(!std::isinf(yGradient_));
    assert(!std::isinf(zGradient_));

    assert(xGradient_ == xGradient_);
    assert(yGradient_ == yGradient_);
    assert(zGradient_ == zGradient_);

    assert(!std::isinf(xGradientGuess_));
    assert(!std::isinf(yGradientGuess_));
    assert(!std::isinf(zGradientGuess_));

    assert(xGradientGuess_ == xGradientGuess_);
    assert(yGradientGuess_ == yGradientGuess_);
    assert(zGradientGuess_ == zGradientGuess_);

    assert(!chemicalSymbol_.empty());
    assert(!mdSymbol_      .empty());

    // Starting enumeration with one.
    assert(nEvbNos_           .count(0) == 0);
    assert(nMarkerNos_        .count(0) == 0);
    assert(nNeighborDistances_.count(0) == 0);
    assert(nQuantumRegionNos_ .count(0) == 0);

    assert(nEvbNos_       .size() == nEvbXGradients_.size());
    assert(nEvbXGradients_.size() == nEvbYGradients_.size());
    assert(nEvbYGradients_.size() == nEvbZGradients_.size());

    assert(nEvbNos_           .size() == evbNoCount_          );
    assert(nNeighborAtomNos_  .size() == neighborAtomNoCount_ );
    assert(nNeighborDistances_.size() == neighborAtomNoCount_ );
    assert(nMarkerNos_        .size() == markerNoCount_       );
    assert(nQuantumRegionNos_ .size() == quantumRegionNoCount_);

    // AtomNos are unique.
    assert(nNeighborDistances_.count(atomNo_) == 0);

    assert(nNeighborAtomNos_.size() == nNeighborDistances_.size());

    // All atomNos from nNeighborAtomNos_ also present in nNeighborDistances_.
    assert(([this]()->bool
    {
      for(const auto& i:nNeighborAtomNos_)
        if(nNeighborDistances_.count(i.second) == 0)
          return false;
      return true;
    })());

    // All atomNos from nNeighborDistances_ also present in nNeighborAtomNos_.
    assert(([this]()->bool
    {
      for(const auto& i:nNeighborDistances_)
        assert(([this,
                 i    ]()->bool
        {
	  for(const auto& j:nNeighborAtomNos_)
	    if(i.first == j.second)
	      return true;
	  return false;
        })());
      return true;
    })());

    // All distances from nNeighborAtomNos_ also present in nNeighborDistances_.
    assert(([this]()->bool
    {
      for(const auto& i:nNeighborDistances_)
        assert(([this,
                 i    ]()->bool
        {
	  for(const auto& j:nNeighborAtomNos_)
	    if(i.second == j.first)
	      return true;
	  return false;
        })());
      return true;
    })());

    // All distances from nNeighborDistances_ also present in nNeighborAtomNos_.
    assert(([this]()->bool
    {
      for(const auto& i:nNeighborDistances_)
        if(nNeighborAtomNos_.count(i.second) == 0)
          return false;
      return true;
    })());

    // EvbNos used consistently.
    assert(([this]()->bool
    {
      for(const auto& i:nEvbNos_)
        if((nEvbXGradients_.count(i) == 0) ||
           (nEvbYGradients_.count(i) == 0) ||
           (nEvbZGradients_.count(i) == 0)   )
          return false;
      return true;
    })());


    return true;


  }



}



#endif // QREG_ATOMBOX_ATOM_HPP_

