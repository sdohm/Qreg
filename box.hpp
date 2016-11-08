//
// box.hpp
//
//  Created on: Aug 12, 2014
//      Author: Sebastian Dohm<sebastian.dohm@uni-ulm.de>
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



#ifndef QREG_ATOMBOX_BOX_HPP_
#define QREG_ATOMBOX_BOX_HPP_



#include "atom.hpp"



#define NDEBUG



#include <algorithm>
#include <cmath>



namespace qreg
{



  //
  // Describes a collection of atoms in a MD environment..
  //
  class Box final
  {



     public :


      //
      // Creates an empty box without atoms and with zero box dimensions.
      //
      Box(double newBoxXDimension =      0.0,
          double newBoxXOrigin    =      0.0,
          double newBoxYDimension =      0.0,
          double newBoxYOrigin    =      0.0,
          double newBoxZDimension =      0.0,
          double newBoxZOrigin    =      0.0,
          double newTimeFrame     =      1.0,
          double newTotalEPot     =      0.0 )noexcept;

      //
      // Copy construction enabled.
      //
      Box(const Box& box)noexcept;

      //
      // Move construction enabled.
      //
      Box(Box&& box)noexcept;

      //
      // Copy assignment enabled.
      //
      auto
      operator=(const Box& box)noexcept
      ->Box&;

      //
      // Move assignment enabled.
      //
      auto
      operator=(Box&& box)noexcept
      ->Box&;

      //
      // Destructor.
      //
      ~Box(void)noexcept;


      //
      // ADL std::swap enabled.
      //
      auto friend
      swap(Box& lhs,
           Box& rhs )
      ->void;


      //
      // Creates nonperiodic clone of this box.
      //
      auto
      cloneNonPeriodic(const double& range)const noexcept
      ->Box;

      //
      // Creates clone of this box sans most of the graphite.
      //
      auto
      cloneShrinkGraphite(const double& range)const noexcept
      ->Box;


      //
      // Computes all distances.
      // Above-quadratic complexity (atom count).
      // Parallelization may improve performance.
      //
      auto
      computeAllDistances(void)noexcept
      ->void;

      //
      // Computes kinetics i.e. this is the thermostat.
      // Linear complexity (atom count).
      //
      auto
      computeKinetics(void)noexcept
      ->void;

      //
      // Recomputes some closest distances.
      // Linearithmetic complexity (atom count)
      // for distance chosen within reason.
      // Above-quadratic complexity (atom count) otherwise.
      // Parallelization may improve performance.
      //
      auto
      computeNearestNeighborDistances(const double& distance)noexcept
      ->void;

      //
      // Deletes all distances to one atom.
      // Linearithmic complexity (atom count).
      //
      auto
      deleteAllDistances(const size_t& atomNo)noexcept
      ->void;

      //
      // Deletes all distances.
      // Quadratic complexity (atom count).
      //
      auto
      deleteAllDistances(void)noexcept
      ->void;

      //
      // Deletes all EPots of all quantum regions.
      // Constant complexity (atom count).
      // Linear complexity (quantum region count).
      //
      auto
      deleteAllEPots(void)noexcept
      ->void;

      //
      // Deletes all EPots of one EVB. Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      deleteAllEvbEPots(const unsigned int& evbNo)noexcept
      ->void;

      //
      // Deletes all EPots of all EVBs. Constant complexity (atom count).
      // Linear complexity (EVB count).
      //
      auto
      deleteAllEvbEPots(void)noexcept
      ->void;

      //
      // Deletes all EVB gradients of one atom.
      // Constant complexity (atom count) for EVB sizes chosen within reason.
      // Linear complexity (atom count) otherwise.
      // Linear complexity (EVB count).
      //
      auto
      deleteAllEvbGradients(const size_t& atomNo)noexcept
      ->void;

      //
      // Deletes all EVB gradients. Linear complexity (atom count).
      // Linear complexity (EVB count).
      //
      auto
      deleteAllEvbGradients(void)noexcept
      ->void;

      //
      // Deletes all total EPots of all EVBs. Constant complexity (atom count).
      // Linear complexity (EVB count).
      //
      auto
      deleteAllEvbTotalEpots(void)noexcept
      ->void;

      //
      // Deletes all markers. Linear complexity (atom count).
      // Linear complexity (marker count).
      //
      auto
      deleteAllMarkerNos(void)noexcept
      ->void;

      //
      // Deletes all quantum regions. Linear complexity (atom count).
      // Linear complexity (quantum region count).
      //
      auto
      deleteAllQuantumRegionNos(void)noexcept
      ->void;

      //
      // Deletes one distance between two atoms.
      // Logarithmic complexity (atom count).
      //
      auto
      deleteDistance(const size_t& atomNo        ,
                     const size_t& neighborAtomNo )noexcept
      ->void;

      //
      // Deletes one EPot of one quantum region.
      // Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      deleteEPot(const unsigned int& quantumRegionNo)noexcept
      ->void;

      //
      // Deletes one EPot of one quantum region within one EVB.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      deleteEvbEPot(const unsigned int& evbNo          ,
                    const unsigned int& quantumRegionNo )noexcept
      ->void;

      //
      // Deletes one EVB gradient of one atom. Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      deleteEvbGradient(const size_t&       atomNo,
                        const unsigned int& evbNo  )noexcept
      ->void;

      //
      // Deletes one EVB gradient. Linear complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      deleteEvbGradient(const unsigned int& evbNo)noexcept
      ->void;

      //
      // Deletes total EPot of one EVB. Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      deleteEvbTotalEPot(const unsigned int& evbNo)noexcept
      ->void;

      //
      // Deletes one marker. Constant complexity (atom count).
      // Logarithmic complexity (marker count).
      //
      auto
      deleteMarkerNo(const unsigned int& markerNo)noexcept
      ->void;

      //
      // Deletes one quantum region. Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      deleteQuantumRegionNo(const unsigned int& quantumRegionNo)noexcept
      ->void;

      //
      // Equalizes all gradients. Linear complexity (atom count).
      //
      auto
      equalizeGradients(void)noexcept
      ->void;

      //
      // Increases the electron count of an atom by one.
      // Constant complexity (atom count).
      //
      auto
      incrementNegativeCharge(const size_t& atomNo)noexcept
      ->void;

      //
      // Decreases the electron count of an atom by one.
      // Constant complexity (atom count).
      //
      auto
      incrementPositiveCharge(const size_t& atomNo)noexcept
      ->void;

      //
      // Inserts atom.
      // Constant complexity (atom count).
      //
      auto
      insertAtom(const size_t&       newAtomNo       ,
                 const std::string&  newChemicalName ,
                 const unsigned int& newElectronCount,
                 const double&       newXCoordinate  ,
                 const double&       newYCoordinate  ,
                 const double&       newZCoordinate   )noexcept
      ->void;

      //
      // Inserts distance between two atoms.
      // Logarithmic complexity (atom count).
      //
      auto
      insertDistance(const size_t& atomNo       ,
                     const size_t& neigborAtomNo,
                     const double& newDistance   )noexcept
      ->void;

      //
      // Inserts EPot of a quantum region. Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      insertEPot(const unsigned int& quantumRegionNo,
                 const double&       newEPot         )noexcept
      ->void;

      //
      // Inserts EPot of a quantum region inside an EVB.
      // Logarithmic complexity (EVB count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      insertEvbEPot(const unsigned int& evbNo          ,
                    const unsigned int& quantumRegionNo,
                    const double&       newEvbEPot      )noexcept
      ->void;

      //
      // Inserts EVB gradient of an atom. Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      insertEvbGradient(const size_t&       atomNo         ,
                        const unsigned int& newEvbNo       ,
                        const double&       newEvbXGradient,
                        const double&       newEvbYGradient,
                        const double&       newEvbZGradient )noexcept
      ->void;

      //
      // Inserts total EPot of an EVB. Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      insertEvbTotalEPot(const unsigned int& evbNo          ,
                         const double&       newEvbTotalEPot )noexcept
      ->void;

      //
      // Inserts marker of an atom. Constant complexity (atom count).
      // Linear complexity (marker count).
      //
      auto
      insertMarkerNo(const size_t&       atomNo     ,
                     const unsigned int& newMarkerNo )noexcept
      ->void;

      //
      // Inserts quantum region of an atom. Constant complexity (atom count).
      // Linear complexity (quantum region count).
      //
      auto
      insertQuantumRegionNo(const size_t&       atomNo            ,
                            const unsigned int& newQuantumRegionNo )noexcept
      ->void;

      //
      // Sets box dimensions.
      // Constant complexity (atom count).
      //
      auto
      setBoxDimension(const double& newBoxXDimension,
                      const double& newBoxYDimension,
                      const double& newBoxZDimension )noexcept
      ->void;

      //
      // Sets coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      setCoordinate(const size_t& atomNo        ,
                    const double& newXCoordinate,
                    const double& newYCoordinate,
                    const double& newZCoordinate )noexcept
      ->void;

      //
      // Sets QM gradient of an atom.
      // Constant complexity (atom count).
      //
      auto
      setGradient(const size_t& atomNo      ,
                  const double& newXGradient,
                  const double& newYGradient,
                  const double& newZGradient )noexcept
      ->void;

      //
      // Sets MD gradient guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      setGradientGuess(const size_t& atomNo           ,
                       const double& newXGradientGuess,
                       const double& newYGradientGuess,
                       const double& newZGradientGuess )noexcept
      ->void;

      //
      // Sets electron count of an atom.
      // Constant complexity (atom count).
      //
      auto
      setElectronCount(const size_t&       atomNo          ,
                       const unsigned int& newElectronCount )noexcept
      ->void;

      //
      // Sets MD symbol of an atom.
      // Constant complexity (atom count).
      //
      auto
      setMdSymbol(const size_t&      atomNo     ,
                  const std::string& newMdSymbol )noexcept
      ->void;

      //
      // Sets QM partial charge of an atom.
      // Constant complexity (atom count).
      //
      auto
      setPartialCharge(const size_t& atomNo          ,
                       const double& newPartialCharge )noexcept
      ->void;

      //
      // Sets MD partial charge guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      setPartialChargeGuess(const size_t& atomNo               ,
                            const double& newPartialChargeGuess )noexcept
      ->void;

      //
      // Sets temperature and thermostat dampening factor.
      // Constant complexity (atom count).
      //
      auto
      setThermostat(const double& newTemperature    ,
                    const double& newDampeningFactor )noexcept
      ->void;

      //
      // Sets time steps for MD simulation.
      // Constant complexity (atom count).
      //
      auto
      setTimeFrame(const double& newTimeFrame)noexcept
      ->void;

      //
      // Sets total EPot.
      // Constant complexity (atom count).
      //
      auto
      setTotalEPot(const double& newTotalEPot)noexcept
      ->void;


      //
      // Checks whether specified atom exists.
      // Linear complexity (atom count).
      //
      auto
      hasAtom(const size_t& atomNo)const noexcept
      ->bool;

      //
      // Checks whether distance between two atom exists.
      // Logarithmic complexity (atom count).
      //
      auto
      hasDistance(const size_t& atomNo        ,
                  const size_t& neighborAtomNo )const noexcept
      ->bool;

      //
      // Checks whether a quantum region has EPot.
      // Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      hasEPot(const unsigned int& quantumRegionNo)const noexcept
      ->bool;

      //
      // Checks whether a quantum region within an EVB has EPot.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      hasEvbEPot(const unsigned int& evbNo          ,
                 const unsigned int& quantumRegionNo )const noexcept
      ->bool;

      //
      // Checks whether an atom is associated with an EVB.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      hasEvbNo(const size_t&       atomNo,
               const unsigned int& evbNo  )const noexcept
      ->bool;

      //
      // Checks whether an EVB exists.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      hasEvbNo(const unsigned int& evbNo)const noexcept
      ->bool;

      //
      // Checks whether an EVB has total EPot.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      hasEvbTotalEPot(const unsigned int& evbNo)const noexcept
      ->bool;

      //
      // Checks whether an atom has a marker.
      // Constant complexity (atom count).
      // Logarithmic complexity (marker count).
      //
      auto
      hasMarkerNo(const size_t&       atomNo  ,
                  const unsigned int& markerNo )const noexcept
      ->bool;

      //
      // Checks whether a marker exists.
      // Constant complexity (atom count).
      // Logarithmic complexity (marker count).
      //
      auto
      hasMarkerNo(const unsigned int& markerNo)const noexcept
      ->bool;

      //
      // Checks whether an atom is associated with a quantum region.
      // Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      hasQuantumRegionNo(const size_t&       atomNo         ,
                         const unsigned int& quantumRegionNo )const noexcept
      ->bool;

      //
      // Checks whether a quantum region exists.
      // Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      hasQuantumRegionNo(const unsigned int& quantumRegionNo)const noexcept
      ->bool;


      //
      // Returns electron count of an atom.
      // Constant complexity (atom count).
      //
      auto
      getElectronCount(const size_t& atomNo)const noexcept
      ->unsigned int;

      //
      // Returns proton count of an atom.
      // Constant complexity (atom count).
      //
      auto
      getProtonCount(const size_t& atomNo)const noexcept
      ->unsigned int;

      //
      // Returns atom count.
      // Constant complexity (atom count).
      //
      auto
      getAtomCount(void)const noexcept
      ->size_t;


      //
      // Computes distance [Angstrom] between two atoms.
      // Constant complexity (atom count).
      //
      auto
      computeDistance(const size_t& atomNo        ,
                      const size_t& neighborAtomNo )const noexcept
      ->double;

      //
      // Computes velocity [Angstrom per picosecond] of an atom.
      // Constant complexity (atom count).
      //
      auto
      computeVelocity(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns x component of the box dimension.
      // Constant complexity (atom count).
      //
      auto
      getBoxXDimension(void)const noexcept
      ->double;

      //
      // Returns x component of the box origin.
      // Constant complexity (atom count).
      //
      auto
      getBoxXOrigin(void)const noexcept
      ->double;

      //
      // Returns y component of the box dimension.
      // Constant complexity (atom count).
      //
      auto
      getBoxYDimension(void)const noexcept
      ->double;

      //
      // Returns y component of the box origin.
      // Constant complexity (atom count).
      //
      auto
      getBoxYOrigin(void)const noexcept
      ->double;

      //
      // Returns z component of the box dimension.
      // Constant complexity (atom count).
      //
      auto
      getBoxZDimension(void)const noexcept
      ->double;

      //
      // Returns z component of the box origin.
      // Constant complexity (atom count).
      //
      auto
      getBoxZOrigin(void)const noexcept
      ->double;

      //
      // Returns MD thermostat dampening factor.
      // Constant complexity (atom count).
      //
      auto
      getDampeningFactor(void)const noexcept
      ->double;

      //
      // Returns distance between two atoms.
      // Logarithmic complexity (atom count).
      //
      auto
      getDistance(const size_t& atomNo        ,
                  const size_t& neighborAtomNo )const noexcept
      ->double;

      //
      // Returns EPot of a quantum region.
      // Constant complexity (atom count).
      // Logarithmic complexity (atom count).
      //
      auto
      getEPot(const unsigned int& quantumRegionNo)const noexcept
      ->double;

      //
      // Returns EPot of a quantum region within an EVB.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      getEvbEPot(const unsigned int& evbNo          ,
                 const unsigned int& quantumRegionNo )const noexcept
      ->double;

      //
      // Returns total EPot of an EVB.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbTotalEPot(const unsigned int& evbNo)const noexcept
      ->double;

      //
      // Returns x component of an EVB gradient of an atom.
      // Constant complexity (atom count).
      // Constant complexity (EVB count).
      //
      auto
      getEvbXGradient(const size_t&       atomNo,
                      const unsigned int& evbNo  )const noexcept
      ->double;

      //
      // Returns y component of an EVB gradient of an atom.
      // Constant complexity (atom count).
      // Constant complexity (EVB count).
      //
      auto
      getEvbYGradient(const size_t&       atomNo,
                      const unsigned int& evbNo  )const noexcept
      ->double;

      //
      // Returns z component of an EVB gradient of an atom.
      // Constant complexity (atom count).
      // Constant complexity (EVB count).
      //
      auto
      getEvbZGradient(const size_t&       atomNo,
                      const unsigned int& evbNo  )const noexcept
      ->double;

      //
      // Returns previous distance between two atoms.
      // Logarithmic complexity (atom count).
      //
      auto
      getLastDistance(const size_t& atomNo        ,
                      const size_t& neighborAtomNo )const noexcept
      ->double;

      //
      // Returns x component of the previous coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getLastXCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns x component of the previous distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getLastXDistance(const size_t& atomNo0,
                       const size_t& atomNo1 )const noexcept
      ->double;

      //
      // Returns y component of the previous coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getLastYCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns y component of the previous distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getLastYDistance(const size_t& atomNo0,
                       const size_t& atomNo1 )const noexcept
      ->double;

      //
      // Returns z component of the previous coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getLastZCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns Z component of the previous distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getLastZDistance(const size_t& atomNo0,
                       const size_t& atomNo1 )const noexcept
      ->double;


      //
      // Returns mass of an atom.
      // Constant complexity (atom count).
      //
      auto
      getMass(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns QM partial charge of an atom.
      // Constant complexity (atom count).
      //
      auto
      getPartialCharge(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns MD partial charge guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      getPartialChargeGuess(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns temperature.
      // Constant complexity (atom count).
      //
      auto
      getTemperature(void)const noexcept
      ->double;

      //
      // Returns time frame for MD simulation.
      // Constant complexity (atom count).
      //
      auto
      getTimeFrame(void)const noexcept
      ->double;

      //
      // Returns total EPot.
      // Constant complexity (atom count).
      //
      auto
      getTotalEPot(void)const noexcept
      ->double;

      //
      // Returns van der Waals radius.
      // Constant complexity (atom count).
      //
      auto
      getVdwRadius(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns x component of the coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getXCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns x component of the distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getXDistance(const size_t& atomNo0,
                   const size_t& atomNo1 )const noexcept
      ->double;

      //
      // Returns x component of the QM gradient of an atom.
      // Constant complexity (atom count).
      //
      auto
      getXGradient(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns x component of the MD gradient guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      getXGradientGuess(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns y component of the coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getYCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns y component of the distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getYDistance(const size_t& atomNo0,
                   const size_t& atomNo1 )const noexcept
      ->double;

      //
      // Returns y component of the QM gradient of an atom.
      // Constant complexity (atom count).
      //
      auto
      getYGradient(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns y component of the MD gradient guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      getYGradientGuess(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns z component of the coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getZCoordinate(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns z component of the distance betwen two atoms.
      // Constant complexity (atom count).
      //
      auto
      getZDistance(const size_t& atomNo0,
                   const size_t& atomNo1 )const noexcept
      ->double;

      //
      // Returns z component of the QM gradient of an atom.
      // Constant complexity (atom count).
      //
      auto
      getZGradient(const size_t& atomNo)const noexcept
      ->double;

      //
      // Returns z component of the MD gradient guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      getZGradientGuess(const size_t& atomNo)const noexcept
      ->double;

      //
      // Advances one MD step ahead.
      // Returns temperature.
      // Linear complexity (atom count).
      //
      auto
      updateCoordinates(void)noexcept
      ->double;


      //
      // Returns chemical symbol of an atom.
      // Constant complexity (atom count).
      //
      auto
      getChemicalSymbol(const size_t& atomNo)const noexcept
      ->std::string;

      //
      // Returns MD symbol of an atom.
      // Constant complexity (atom count).
      //
      auto
      getMdSymbol(const size_t& atomNo)const noexcept
      ->std::string;


      //
      // Computes acceleration [Angstrom per Picosecond squared] on an atom.
      // onstant complexity (atom count).
      //
      auto
      computeAcceleration(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Computes force [Piconewton] on an atom.
      // Constant complexity (atom count).
      //
      auto
      computeForce(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Computes velocity vector [Angstrom per picosecond] of an atom.
      // Constant complexity (atom count).
      //
      auto
      computeVelocityVector(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Returns previous coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getLastCoordinate(const size_t& atomNo)const noexcept
      ->std::vector<double>;
 
      //
      // Returns coordinate of an atom.
      // Constant complexity (atom count).
      //
      auto
      getCoordinate(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Returns distance vector betwen two atom.
      // Constant complexity (atom count).
      //
      auto
      getDistanceVector(const size_t& atomNo0,
                        const size_t& atomNo1 )const noexcept
      ->std::vector<double>;

      //
      // Returns EVB gradient of an atom.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      getEvbGradient(const size_t&       atomNo,
                     const unsigned int& evbNo  )const noexcept
      ->std::vector<double>;

      //
      // Returns gradient of an atom.
      // Constant complexity (atom count).
      //
      auto
      getGradient(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Returns gradient guess of an atom.
      // Constant complexity (atom count).
      //
      auto
      getGradientGuess(const size_t& atomNo)const noexcept
      ->std::vector<double>;

      //
      // Returns previous distance vector betwen two atom.
      // Constant complexity (atom count).
      //
      auto
      getLastDistanceVector(const size_t& atomNo0,
                            const size_t& atomNo1 )const noexcept
      ->std::vector<double>;

      //
      // Returns list of N EVB numbers of an atom.
      // Constant complexity (atom count).
      // Linear complexity (EVB count).
      //
      auto
      getNEvbNos(const size_t& atomNo)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N EVB numbers.
      // Constant complexity (atom count).
      // Linear complexity (EVB count).
      //
      auto
      getNEvbNos(void)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N atom numbers.
      // Constant complexity (atom count).
      //
      auto
      getNAtomNos(void)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N atom numbers associated with an EVB.
      // Constant complexity (atom count).
      // Logarithmic complexity (EVB count).
      //
      auto
      getNAtomNosByEvbNo(const unsigned int& evbNo)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N atom numbers associated with a marker.
      // Constant complexity (atom count).
      // Logarithmic complexity (marker count).
      //
      auto
      getNAtomNosByMarkerNo(const unsigned int& markerNo)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N atom numbers associated with a quantum region.
      // Constant complexity (atom count).
      // Logarithmic complexity (quantum region count).
      //
      auto
      getNAtomNosByQuantumRegionNo(const unsigned int& quantumRegionNo)
                                  const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N markers of an atom.
      // Constant complexity (atom count).
      // Linear complexity (marker count).
      //
      auto
      getNMarkerNos(const size_t& atomNo)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N markers.
      // Constant complexity (atom count).
      // Linear complexity (marker count).
      //
      auto
      getNMarkerNos(void)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N neighbor atom numbers sorted by ascending distances
      // within given range of specified atom.
      // Constant complexity (atom count) for range chosen within
      // reason. Linear complexity (atom count) otherwise.
      //
      auto
      getNNeighborAtomNos(const size_t& atomNo,
                          const double& range  )const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N neighbor atom numbers associated with specified atom
      // sorted by ascending distances.
      // Linear complexity (atom count).
      //
      auto
      getNNeighborAtomNos(const size_t& atomNo)const noexcept
      ->std::vector<size_t>;

      //
      // Returns list of N quantum regions associated with specified atom.
      // Constant complexity (atom count).
      // Linear complexity (quantum region count).
      //
      auto
      getNQuantumRegionNos(const size_t& atomNo)const noexcept
      ->std::vector<unsigned int>;

      //
      // Returns list of N quantum regions.
      // Constant complexity (atom count).
      // Linear complexity (quantum region count).

      //
      auto
      getNQuantumRegionNos(void)const noexcept
      ->std::vector<unsigned int>;
  


    private :



      //
      // Computes all distances for one atom instance.
      // Causes temporary DbC sanity violation:
      // i.e. neighbor lists going out of sync!
      // Linearithmetic complexity (atom count).
      //
      auto
      computeAllDistancesForOneAtom(const size_t              atomNo ,
                                    const std::vector<size_t> nAtomNos )
                                   noexcept
      ->void;

      //
      // Recomputes closest distances for one atom instance.
      // Causes temporary DbC sanity violation:
      // i.e. neighbor lists going out of sync!
      // Logarithmic complexity (atom count)
      // for reasonable distance.
      // Up to above-cubic complexity (atom count) otherwise.
      //
      auto
      computeNearestNeighborDistancesForOneAtom
      (const size_t atomNo  ,
       const double distance )noexcept
      ->void;


      //
      // Checks sanity i.e. Design by Contract compliance.
      //
      auto
      isSane(void)const noexcept
      ->bool;


      //
      // Locking and unlocking the Box for multithreading.
      //
      mutable std::mutex threadSafety_;


      //
      // Stores x component of the box dimensions in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxXDimension_;

      //
      // Stores x component of the box origin in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxXOrigin_;

      //
      // Stores y component of the box dimensions in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxYDimension_;

      //
      // Stores y component of the box origin in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxYOrigin_;

      //
      // Stores z component of the box dimensions in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxZDimension_;

      //
      // Stores z component of the box origin in angstrom.
      // Assign non-negative values. Default "0".
      //
      double boxZOrigin_;

      //
      // Stores a dampening factor for MD thermostat.
      // Assign non-negative non-zero values. Default "1".
      //
      double dampeningFactor_;

      //
      // Stores temperature in Kelvin.
      // Assign non-negative non-zero values. Default "298.15".
      //
      double temperature_;

      //
      // Stores length of the time frame in femtoseconds.
      // Assign non-negative non-zero values. Default "1.0".
      //
      double timeFrame_;

      //
      // Stores total Epot in E-23 J. Default "0".
      //
      double totalEPot_;


      //
      // Stores N atom numbers. Assign non-zero continuous key values.
      // Empty on creation.
      //
      std::vector<size_t> nAtomNos_;


      //
      // Stores N atoms sorted by their numbers minus one.
      // Empty on creation.
      //
      std::vector<Atom>nAtoms_;


      //
      // Stores N Epots in E-23 J of quantum regions
      // sorted by numbers of quantum regions.
      // Assign non-zero key values. Empty on creation.
      //
      std::map<unsigned int,
               double       > nEPots_;

      //
      // Stores N total Epots in E-23 J for EVBs sorted by numbers of EVBs.
      // Assign non-zero key values. Empty on creation.
      //
      std::map<unsigned int,
               double       > nEvbTotalEPots_;

      //
      // Stores N atom numbers for M EVB numbers.
      //
      std::map<unsigned int       ,
               std::vector<size_t> > mnEvbAtomNos_;

      //
      // Stores N atom numbers for M markers.
      //
      std::map<unsigned int       ,
               std::vector<size_t> > mnMarkerAtomNos_;

      //
      // Stores N atom numbers for M quantum regions.
      //
      std::map<unsigned int       ,
               std::vector<size_t> > mnQuantumRegionAtomNos_;

      //
      // Stores N Epots in E-23 J of N quantum regions
      // sorted by N numbers of quantum regions for M EVBs
      // sorted by M numbers of M EVBs.
      // Assign non-zero key values. Empty on creation.
      //
      std::map<unsigned int           ,
               std::map<unsigned int,
                        double       > > mnEvbEPots_;


      //
      // Stores EVB number count. Default "0".
      //
      decltype(mnEvbAtomNos_.size()) evbNoCount_;

      //
      // Stores marker number count. Default "0".
      //
      decltype(mnMarkerAtomNos_.size()) markerNoCount_;

      //
      // Stores quantum region number count. Default "0".
      //
      decltype(mnQuantumRegionAtomNos_.size()) quantumRegionNoCount_;


      //
      // Stores atom count. Default "0".
      //
      decltype(nAtoms_.size()) atomCount_;



  };



  inline auto
  Box::computeAllDistancesForOneAtom(const size_t              atomNo  ,
                                     const std::vector<size_t> nAtomNos )
                                    noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //
    assert(atomNo != 0);


    const auto aX = nAtoms_[atomNo - 1].getXCoordinate();
    const auto aY = nAtoms_[atomNo - 1].getYCoordinate();
    const auto aZ = nAtoms_[atomNo - 1].getZCoordinate();


    //
    // DbC sanity disregarded for this method.
    // Neighbor lists will go out of sync temporarily.
    // Caller has responsibility to fix this.
    //
    for(const auto& i:nAtomNos)
      if(i != atomNo)
      {


        const auto bX = nAtoms_[i - 1].getXCoordinate();
        const auto bY = nAtoms_[i - 1].getYCoordinate();
        const auto bZ = nAtoms_[i - 1].getZCoordinate();

        const auto dX = std::min(std::abs(bX - aX)                            ,
                                 std::min(std::abs(bX - aX - boxXDimension_),
                                          std::abs(bX - aX + boxXDimension_) ) );

        const auto dY = std::min(std::abs(bY - aY)                            ,
                                 std::min(std::abs(bY - aY - boxYDimension_),
                                          std::abs(bY - aY + boxYDimension_) ) );

        const auto dZ = std::min(std::abs(bZ - aZ)                            ,
                                 std::min(std::abs(bZ - aZ - boxZDimension_),
                                          std::abs(bZ - aZ + boxZDimension_) ) );


        //
        // Design by Contract Postcondition
        //

        assert(((boxXDimension_ == 0.0) &&
                (boxYDimension_ == 0.0) &&
                (boxZDimension_ == 0.0)   ) ||
               ((dX <= boxXDimension_) &&
                (dY <= boxYDimension_) &&
                (dZ <= boxZDimension_)   )    );

        assert(dX * dX + dY * dY + dZ * dZ >= 0.0);


        std::lock_guard<std::mutex> lock(threadSafety_);

        nAtoms_[atomNo - 1].insertDistance(i                                ,
                                           sqrt(dX * dX + dY * dY + dZ * dZ) );


    }


  }



  inline auto
  Box::computeNearestNeighborDistancesForOneAtom
  (const size_t atomNo  ,
   const double distance )noexcept
  ->void
  {


    //
    // Design by Contract Precondition
    //
    assert(atomNo != 0);


    const auto aX = nAtoms_[atomNo - 1].getXCoordinate();
    const auto aY = nAtoms_[atomNo - 1].getYCoordinate();
    const auto aZ = nAtoms_[atomNo - 1].getZCoordinate();

    std::vector<size_t> nNeighborAtomNos = nAtoms_[atomNo - 1]
                                           .getNNeighborAtomNos(distance);


    //
    // DbC sanity disregarded for this method.
    // Neighbor lists will go out of sync temporarily.
    // Caller has responsibility to fix this.
    //
    for(const auto& i:nNeighborAtomNos)
      if(i != atomNo)
      {


        const auto bX = nAtoms_[i - 1].getXCoordinate();
        const auto bY = nAtoms_[i - 1].getYCoordinate();
        const auto bZ = nAtoms_[i - 1].getZCoordinate();

        const auto dX = std::min(std::abs(bX - aX)                            ,
                                 std::min(std::abs(bX - aX - boxXDimension_),
                                          std::abs(bX - aX + boxXDimension_) ) );

        const auto dY = std::min(std::abs(bY - aY)                            ,
                                 std::min(std::abs(bY - aY - boxYDimension_),
                                          std::abs(bY - aY + boxYDimension_) ) );

        const auto dZ = std::min(std::abs(bZ - aZ)                            ,
                                 std::min(std::abs(bZ - aZ - boxZDimension_),
                                          std::abs(bZ - aZ + boxZDimension_) ) );


        //
        // Design by Contract Postcondition
        //

        assert(((boxXDimension_ == 0.0) &&
                (boxYDimension_ == 0.0) &&
                (boxZDimension_ == 0.0)   ) ||
               ((dX <= boxXDimension_) &&
                (dY <= boxYDimension_) &&
                (dZ <= boxZDimension_)   )    );

        assert(dX * dX + dY * dY + dZ * dZ >= 0.0);


        std::lock_guard<std::mutex> lock(threadSafety_);

        nAtoms_[atomNo - 1].insertDistance(i                                ,
                                           sqrt(dX * dX + dY * dY + dZ * dZ) );


    }


  }



  inline auto
  Box::isSane(void)const noexcept
  ->bool
  {


    //
    // Design by Contract Invariant
    //

    assert((boxXDimension_ == 0.0 &&
            boxYDimension_ == 0.0 &&
            boxZDimension_ == 0.0   ) ||
           (boxXDimension_ > 0.0 &&
            boxYDimension_ > 0.0 &&
            boxZDimension_ > 0.0    )   );

    assert(!std::isinf(boxXDimension_));
    assert(!std::isinf(boxYDimension_));
    assert(!std::isinf(boxZDimension_));

    assert(boxXDimension_ == boxXDimension_);
    assert(boxYDimension_ == boxYDimension_);
    assert(boxZDimension_ == boxZDimension_);

    assert(!std::isinf(boxXOrigin_));
    assert(!std::isinf(boxYOrigin_));
    assert(!std::isinf(boxZOrigin_));

    assert(boxXOrigin_ == boxXOrigin_);
    assert(boxYOrigin_ == boxYOrigin_);
    assert(boxZOrigin_ == boxZOrigin_);

    assert(!std::isinf(temperature_));
    assert(temperature_ == temperature_);
    assert(temperature_ >= 0.0         );

    assert(!std::isinf(timeFrame_));
    assert(timeFrame_ == timeFrame_);
    assert(timeFrame_ >  0.0       );

    assert(!std::isinf(dampeningFactor_));
    assert(dampeningFactor_ == dampeningFactor_);
    assert(dampeningFactor_ >= timeFrame_      );

    assert(timeFrame_ <= dampeningFactor_);

    assert(!std::isinf(totalEPot_));
    assert(totalEPot_ == totalEPot_);


    // Starting enumeration with one.
    assert(mnEvbEPots_            .count(0) == 0);
    assert(mnEvbAtomNos_          .count(0) == 0);
    assert(mnMarkerAtomNos_       .count(0) == 0);
    assert(mnQuantumRegionAtomNos_.count(0) == 0);
    assert(nEPots_                .count(0) == 0);
    assert(nEvbTotalEPots_        .count(0) == 0);
    assert(std::find(nAtomNos_.cbegin(),
                     nAtomNos_.cend  (),
                     0                  ) == nAtomNos_.cend());

    assert(mnEvbAtomNos_          .size() == evbNoCount_          );
    assert(mnMarkerAtomNos_       .size() == markerNoCount_       );
    assert(mnQuantumRegionAtomNos_.size() == quantumRegionNoCount_);
    assert(nAtomNos_              .size() == atomCount_           );
    assert(nAtoms_                .size() == atomCount_           );

    assert(mnEvbEPots_    .size() <= evbNoCount_          );
    assert(nEPots_        .size() <= quantumRegionNoCount_);
    assert(nEvbTotalEPots_.size() <= quantumRegionNoCount_);

    // Box delimiters are being observed.
    assert(([this]()->bool
    {
      if((boxXDimension_ != 0.0) &&
         (boxYDimension_ != 0.0) &&
         (boxZDimension_ != 0.0)   )
      {
        for(const auto& i:nAtoms_)
        {
          if(((i.getXCoordinate() <  boxXOrigin_)                  ||
              (i.getXCoordinate() >= boxXOrigin_ + boxXDimension_)   ) ||
             ((i.getYCoordinate() <  boxYOrigin_)                  ||
              (i.getYCoordinate() >= boxYOrigin_ + boxYDimension_)   ) ||
             ((i.getZCoordinate() <  boxZOrigin_)                  ||
              (i.getZCoordinate() >= boxZOrigin_ + boxZDimension_)   )   )
            return false;
        }
      }
      
      return true;
    })());

    // NAtomNos_ assigned continuous numbers.
    assert(([this]()->bool
    {
      for(size_t i  = 0               ;
                 i != nAtomNos_.size();
               ++i                     )
        if(nAtomNos_[i] != i + 1)
          return false;
      return true;
    })());

    // Quantum region numbers used consistently.
    assert(([this]()->bool
    {
      for(const auto& i:nEPots_)
        if(mnQuantumRegionAtomNos_.count(i.first) == 0)
          return false;
      return true;
    })());

    // EVB numbers used consistently.
    assert(([this]()->bool
    {
      for(const auto& i:nEvbTotalEPots_)
        if(mnEvbAtomNos_.count(i.first) == 0)
          return false;
      return true;
    })());

    // Quantum region numbers used consistently with EVB numbers.
    assert(([this]()->bool
    {
      for(const auto& i:mnEvbEPots_)
      {
	if(mnEvbAtomNos_.count(i.first) == 0)
	  return false;
        for(const auto& j:i.second)
          if(mnQuantumRegionAtomNos_.count(j.first) == 0)
            return false;
      }
      return true;
    })());

    // All atom numbers from atoms also present in the box.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
      {
	if(std::find(nAtomNos_.cbegin(),
	             nAtomNos_.cend  (),
	             i.getAtomNo()      ) == nAtomNos_.cend())
	  return false;
      }
      return true;
    })());

    // All neighbor atom numbers from atoms also present in the box.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
      {
	for(const auto& j:i.getNNeighborAtomNos())
	  if(std::find(nAtomNos_.cbegin(),
		       nAtomNos_.cend  (),
		       j                  ) == nAtomNos_.cend())
	    return false;
      }
      return true;
    })());

    // All EVB numbers from atoms also present in the box.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
      {
        for(const auto& j:i.getNEvbNos())
          if((mnEvbAtomNos_.count(j) == 0)                        ||
             (std::find(mnEvbAtomNos_.find(j)->second.cbegin(),
                        mnEvbAtomNos_.find(j)->second.cend  (),
                        i.getAtomNo()                          )
              == mnEvbAtomNos_.find(j)->second.cend()           )   )
            return false;
      }
      return true;
    })());

    // All marker numbers from atoms also present in the box.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
      {
        for(const auto& j:i.getNMarkerNos())
          if((mnMarkerAtomNos_.count(j) == 0)                        ||
             (std::find(mnMarkerAtomNos_.find(j)->second.cbegin(),
                        mnMarkerAtomNos_.find(j)->second.cend  (),
                        i.getAtomNo()                             )
              == mnMarkerAtomNos_.find(j)->second.cend()           )   )
            return false;
      }
      return true;
    })());

    // All quantum region numbers from atoms also present in the box.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
      {
        for(const auto& j:i.getNQuantumRegionNos())
          if((mnQuantumRegionAtomNos_.count(j) == 0)                        ||
             (std::find(mnQuantumRegionAtomNos_.find(j)->second.cbegin(),
                        mnQuantumRegionAtomNos_.find(j)->second.cend  (),
                        i.getAtomNo()                                    )
              == mnQuantumRegionAtomNos_.find(j)->second.cend()           )   )
            return false;
      }
      return true;
    })());

    // All atom numbers from the box also present in atoms.
    assert(([this]()->bool
    {
      for(const auto& i:nAtomNos_)
	if(nAtoms_[i - 1].getAtomNo() != i)
	  return false;
      return true;
    })());


    // All EVB numbers from the box also present in atoms.
    #ifndef NDEBUG
      for(const auto& i:mnEvbAtomNos_)
        assert(([this,
                 i    ]()->bool
        {
          for(const auto& j:nAtoms_)
            if(j.hasEvbNo(i.first))
	      return true;
          return false;
        })());
    #endif
    assert(([this]()->bool
    {
      for(const auto& i:mnEvbAtomNos_)
	for(const auto& j:i.second)
	  if(!nAtoms_[j - 1].hasEvbNo(i.first))
	    return false;
      return true;
    })());

    // All marker numbers from the box also present in atoms.
    #ifndef NDEBUG
      for(const auto& i:mnMarkerAtomNos_)
        assert(([this,
                 i    ]()->bool
        {
          for(const auto& j:nAtoms_)
            if(j.hasMarkerNo(i.first))
              return true;
          return false;
        })());
    #endif
    assert(([this]()->bool
    {
      for(const auto& i:mnMarkerAtomNos_)
	for(const auto& j:i.second)
	  if(!nAtoms_[j - 1].hasMarkerNo(i.first))
	    return false;
      return true;
    })());

    // All quantum region numbers from the box also present in atoms.
    #ifndef NDEBUG
      for(const auto& i:mnQuantumRegionAtomNos_)
        assert(([this,
                 i    ]()->bool
        {
          for(const auto& j:nAtoms_)
            if(j.hasQuantumRegionNo(i.first))
              return true;
          return false;
        })());
    #endif
    assert(([this]()->bool
    {
      for(const auto& i:mnQuantumRegionAtomNos_)
        for(const auto& j:i.second)
	  if(!nAtoms_[j - 1].hasQuantumRegionNo(i.first))
	    return false;
      return true;
    })());

    // Neighbor atom lists consistent.
    assert(([this]()->bool
    {
      for(const auto& i:nAtoms_)
	for(const auto& j:i.getNNeighborAtomNos())
	  if(!nAtoms_[j - 1].hasNeighborAtomNo(i.getAtomNo()))
	    return false;
      return true;
    })());


    return true;


  }



}



#endif // QREG_ATOMBOX_BOX_HPP_//
