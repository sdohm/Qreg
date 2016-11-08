//
// idconverter.hpp
//
//  Created on: Aug 1,2015
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



#ifndef QREG_INTERFACE_IDCONVERTER_HPP_
#define QREG_INTERFACE_IDCONVERTER_HPP_



#include <cassert>
#include <map>
#include <mutex>



namespace qreg
{



  //
  // Enables ID conversion between different numbering schemes.
  //
  template<typename InternalId_T,
           typename ExternalId_T >
  class IdConverter final
  {



    public:


      //
      // Constructor.
      //
      IdConverter<InternalId_T,
                  ExternalId_T >(void)noexcept;

      //
      // Copy construction enabled.
      //
      IdConverter<InternalId_T,
                  ExternalId_T >(const IdConverter& idConverter)noexcept;

      //
      // Move construction disabled.
      //
      IdConverter<InternalId_T,
                  ExternalId_T >(IdConverter&& idConverter) = delete;

      //
      // Copy assignment enabled.
      //
      auto
      operator=(const IdConverter& idConverter)noexcept
      ->IdConverter<InternalId_T,
                    ExternalId_T >&;

      //
      // Move assignment disabled.
      //
      auto
      operator=(IdConverter&& idConverter)
      ->IdConverter<InternalId_T,
                    ExternalId_T >& = delete;

      //
      // Destructor.
      //
      ~IdConverter<InternalId_T,
                   ExternalId_T >(void)noexcept;



      //
      // Deletes all IDs.
      //
      auto
      deleteAllIds(void)noexcept
      ->void;

      //
      // Deletes an ID pair by external ID.
      //
      auto
      deleteExternalId(const ExternalId_T& externalId)noexcept
      ->void;

      //
      // Deletes an ID pair by Internal ID.
      //
      auto
      deleteInternalId(const InternalId_T& internalId)noexcept
      ->void;

      //
      // Inserts a pair of Internal ID and external ID.
      //
      auto
      insertIds(const InternalId_T& newInternalId,
                const ExternalId_T& newExternalId )noexcept
      ->void;


      //
      // Checks whether an external ID exists.
      //
      auto
      hasExternalId(const ExternalId_T& newExternalId)const noexcept
      ->bool;

      //
      // Checks whether an ID pair exists.
      //
      auto
      hasIdPair(const InternalId_T& newInternalId,
                const ExternalId_T& newExternalId )const noexcept
      ->bool;

      //
      // Checks whether an internal ID exists.
      //
      auto
      hasInternalId(const InternalId_T& newInternalId)const noexcept
      ->bool;


      //
      // Returns an internal ID by external ID.
      //
      auto
      getInternalId(const ExternalId_T& externalId)const noexcept
      ->InternalId_T;


      //
      // Returns an external ID by internal ID.
      //
      auto
      getExternalId(const InternalId_T& internalId)const noexcept
      ->ExternalId_T;


    private:


      //
      // Checks DbC sanity.
      //
      auto
      isSane(void)const noexcept
      ->bool;


      //
      // Locks and unlocks IdConverter for exception and thread safety.
      //
      std::mutex threadSafety_;

      
      //
      // Stores Internal IDs accessible via external IDs.
      //
      std::map<ExternalId_T,
               InternalId_T > nInternalIds_;


      //
      // Stores external IDs accessible via Internal IDs.
      //
      std::map<InternalId_T,
               ExternalId_T > nExternalIds_;



  };



  template<typename InternalId_T,
           typename ExternalId_T >
  IdConverter<InternalId_T,
              ExternalId_T >::IdConverter(void)noexcept:

    nInternalIds_(),
    nExternalIds_()

  {


    //
    // DbC INV
    //
    assert(isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  IdConverter<InternalId_T,
              ExternalId_T >::IdConverter(const IdConverter& idConverter)noexcept
  {


    //
    // DbC INV
    //
    assert(idConverter.isSane());


    // Calling copy assignment operator.
    *this = idConverter;


    //
    // DbC INV
    //

    assert(idConverter.isSane());
    assert(            isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  auto
  IdConverter<InternalId_T,
              ExternalId_T >::operator=(const IdConverter& idConverter)noexcept
  ->IdConverter&
  {


    //
    // DbC INV
    //
    assert(idConverter.isSane());


    //
    // Perform copy assignment.
    //

    nInternalIds_ = idConverter.nInternalIds_;
    nExternalIds_ = idConverter.nExternalIds_;


    //
    // DbC INV
    //

    assert(idConverter.isSane());
    assert(            isSane());


    return *this;


  }


  template<typename InternalId_T,
           typename ExternalId_T >
  IdConverter<InternalId_T,
              ExternalId_T >::~IdConverter(void)noexcept
  {



  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::deleteAllIds(void)noexcept
  ->void
  {


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::deleteAllIds(void)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nInternalIds_.clear();
    nExternalIds_.clear();


    //
    // DbC POST
    //

    assert(nInternalIds_.empty());
    assert(nExternalIds_.empty());

    assert(isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::deleteExternalId(const ExternalId_T& externalId)noexcept
  ->void
  {


    //
    // DbC PRE
    //
    assert(nInternalIds_.count(externalId) != 0);


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::deleteExternalId(const ExternalId_T& externalId)noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(typename std::map<InternalId_T,
                          ExternalId_T >::iterator i  = nExternalIds_.begin();
                                                   i != nExternalIds_.end  ();
                                                                              ) 
    {

      if(i->second == externalId)
        nExternalIds_.erase(i++);
      else
        ++i;

    }

    nInternalIds_.erase(externalId);


    //
    // DbC POST
    //

    assert(nInternalIds_.count(externalId) == 0);

    assert(isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::deleteInternalId(const InternalId_T& internalId)noexcept
  ->void
  {


    //
    // DbC PRE
    //
    assert(nExternalIds_.count(internalId) != 0);


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::deleteInternalId(const InternalId_T& internalId)noexcept
    // ->void

    std::lock_guard<std::mutex> lock(threadSafety_);

    for(typename std::map<ExternalId_T,
                          InternalId_T >::iterator i  = nInternalIds_.begin();
                                                   i != nInternalIds_.end  ();
                                                                              )
    {

      if(i->second == internalId)
        nInternalIds_.erase(i++);
      else
        ++i;

    }

    nExternalIds_.erase(internalId);


    //
    // DbC POST
    //

    assert(nExternalIds_.count(internalId) == 0);

    assert(isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::insertIds(const InternalId_T& newInternalId,
                                        const ExternalId_T& newExternalId )noexcept
  ->void
  {


    //
    //  DbC PRE
    //

    assert(nInternalIds_.count(newExternalId) == 0);
    assert(nExternalIds_.count(newInternalId) == 0);

    assert(newInternalId != 0);


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::insertIds(const InternalId_T& newInternalId,
    //                                       const ExternalId_T& newExternalId )noexcept
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    nInternalIds_.emplace(newExternalId,
                          newInternalId );
    nExternalIds_.emplace(newInternalId,
                          newExternalId );


    //
    // DbC POST
    //

    assert(nInternalIds_.count(newExternalId) != 0);
    assert(nExternalIds_.count(newInternalId) != 0);

    assert(isSane());


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::hasExternalId(const ExternalId_T& externalId)const noexcept
  ->bool
  {


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::hasExternalId(const ExternalId_T& externalId)const noexcept
    // ->bool
    //

    return (nInternalIds_.count(externalId) != 0);


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::hasIdPair(const InternalId_T& internalId,
                                        const ExternalId_T& externalId )const noexcept
  ->bool
  {


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::hasIdPair(const InternalId_T& internalId,
    //                                       const ExternalId_T& externalId )const noexcept
    // ->bool
    //

    return (getInternalId(externalId) == internalId);


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::hasInternalId(const InternalId_T& internalId)const noexcept
  ->bool
  {


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::hasInternalId(const InternalId_T& internalId)const noexcept
    // ->bool
    //

    return (nExternalIds_.count(internalId) != 0);


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::getInternalId(const ExternalId_T& externalId)const noexcept
  ->InternalId_T
  {


    //
    // DbC PRE
    //
    assert(nInternalIds_.count(externalId) != 0);


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::getInternalId(const ExternalId_T& externalId)const noexcept
    // ->TInternal
    //

    return nInternalIds_.find(externalId)->second;


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::getExternalId(const InternalId_T& internalId)const noexcept
  ->ExternalId_T
  {


    //
    // DbC PRE
    //
    assert(nExternalIds_.count(internalId) != 0);


    //
    // DbC INV
    //
    assert(isSane());


    //
    // template<typename InternalId_T,
    //          typename ExternalId_T >
    // inline auto
    // IdConverter<InternalId_T,
    //             ExternalId_T >::getExternalId(const InternalId_T& internalId)const noexcept
    // ->ExternalId_T
    //

    return nExternalIds_.find(internalId)->second;


  }



  template<typename InternalId_T,
           typename ExternalId_T >
  inline auto
  IdConverter<InternalId_T,
              ExternalId_T >::isSane(void)const noexcept
  ->bool
  {


    //
    // DbC INV
    //

    assert(nInternalIds_.size() == nExternalIds_.size());

    assert(([this]()->bool
    {

      bool isConsistent = true;

      for(const auto& idPair:nExternalIds_)
        if((nInternalIds_.count(idPair.second)        == 0           ) ||
           (nInternalIds_.find(idPair.second)->second != idPair.first)   )
          isConsistent = false;

      for(const auto& idPair:nInternalIds_)
        if((nExternalIds_.count(idPair.second)        == 0           ) ||
           (nExternalIds_.find(idPair.second)->second != idPair.first)   )
          isConsistent = false;

      return isConsistent ? true : false;

    })());


    return true;


  }



}



#endif // QREG_INTERFACE_IDCONVERTER_HPP_
