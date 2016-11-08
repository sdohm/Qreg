//
// qreg.hpp
//
//  Created on: Nov 17,2014
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



#ifndef QREG_QREG_HPP_
#define QREG_QREG_HPP_



#include "interface.hpp"
#include "partitioner.hpp"



#define NDEBUG
#undef  DEBUG_CONFIG
#undef  DEBUG_DISTANCES
#undef  DEBUG_MD_DETAILS
#undef  DEBUG_MD_ACCELERATION
#undef  DEBUG_MD_VELOCITIES



#define QREG_ORC          (1001)
#define QREG_HIT          (1002)
#define QREG_TFF          (1003)
#define QREG_EVB          (1004)
#define QREG_GFF          (1005)
#define QREG_UFF          (1006)

#define QREG_LINEARFIT    (2001)
#define QREG_SMOOTHSTEP   (2002)
#define QREG_SMOOTHERSTEP (2003)



namespace qreg
{



  //
  // Main program.
  //
  class Qreg final
  {



    public:


      //
      // Constructor requires path to config file.
      //
      Qreg(const std::string&  configFilePath              ,
           const unsigned int& newAutoToolChoice = QREG_ORC,
           const double&       newAutoCharge     = 1        );

      //
      // Copy construction disabled.
      //
      Qreg(const Qreg& qreg) = delete;

      //
      // Move construction enabled.
      //
      Qreg(Qreg&& qreg);

      //
      // Copy assignment disabled.
      //
      auto
      operator=(const Qreg& qreg)
      ->Qreg& = delete;

      //
      // Move assignment enabled.
      //
      auto
      operator=(Qreg&& qreg)
      ->Qreg&;

      //
      // Destructor.
      //
      ~Qreg(void);


      //
      // Run Qreg.
      //
      auto
      advance(Box& box)
      ->void;


    private:


      //
      // Adds attraction towards atomNo1 to atomNo2 in E -23 J.
      // Example:
      // Adding 1000 E -23 J to the gradient of atomNo2 towards atomNo1.
      // add ATTRACTION atomNo1 atomNo2 1000.0
      //
      auto
      addAttraction(Box&                            box    ,
                    const std::vector<std::string>& command,
                    Logger&                         log     )const
      ->void;


      //
      // Advances to the next MD step.
      // Usage:
      // advance MD
      //
      auto
      advanceMd(Box&    box,
                Logger& log )const
      ->void;


      //
      // Computes the gradient buffer.
      // Usage:
      // compute BUFFER
      //
      auto
      computeBuffer(Box&                box     ,
                    const unsigned int& markerNo,
                    Logger&             log      )const
      ->void;


      //
      // Computes all distances.
      // This requires considerable computational effort.
      // Usage:
      // compute DISTANCES
      //
      auto
      computeDistances(Box&    box,
                       Logger& log )const
      ->void;


      //
      // Retrieves gradient from automatically determined calculation of the current QM region. 
      // Usage:
      // harvest AUTO
      //
      auto
      harvestAuto(Box&                                     box         ,
                  std::map< unsigned int          ,
                            std::unique_ptr<IfEvb> >&      nIfEvbs     ,
                  std::map< unsigned int          ,
                            std::unique_ptr<IfUff> >&      nIfUffs     ,
                  std::map< unsigned int               ,
                            std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs,
                  std::map< unsigned int              ,
                            std::unique_ptr<IfOpenUff> >&  nIfOpenUffs ,
                  std::map< unsigned int          ,
                            std::unique_ptr<IfOrc> >&      nIfOrcs     ,
                  Logger&                                  log          )const
      ->void;


      //
      // Retrieves gradient from all EVB calculations.
      // Usage:
      // harvest EVB
      //
      auto
      harvestEvb(Box&                                box    ,
                 std::map< unsigned int          ,
                           std::unique_ptr<IfEvb> >& nIfEvbs,
                 Logger&                             log     )const
      ->void;


      //
      // Retrieves gradient from all force field calculations.
      // Usage:
      // harvest FF
      //
      auto
      harvestFf(Box&                                     box         ,
                std::map< unsigned int          ,
                          std::unique_ptr<IfUff> >&      nIfUffs     ,
                std::map< unsigned int               ,
                          std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs,
                std::map< unsigned int              ,
                          std::unique_ptr<IfOpenUff> >&  nIfOpenUffs ,
                Logger&                                  log          )const
      ->void;


      //
      // Retrieves gradient from all quantumchemical calculations.
      // Usage:
      // harvest QM
      //
      auto
      harvestQm(Box&                                box    ,
                std::map< unsigned int          ,
                          std::unique_ptr<IfOrc> >& nIfOrcs,
                Logger&                             log     )const
      ->void;


      //
      // Includes an atom into the current QM region.
      // Example:
      // Includes atom 42 into current QM region.
      // include ATOM 42
      //
      auto
      includeAtom(Box&                            box            ,
                  const unsigned int&             quantumRegionNo,
                  const std::vector<std::string>& command        ,
                  Logger&                         log             )const
      ->void;


      //
      // Includes all atoms into the current QM region.
      // Usage:
      // include ATOMS all
      //
      auto
      includeAtomsAll(Box&                box            ,
                      const unsigned int& quantumRegionNo,
                      Logger&             log             )const
      ->void;


      //
      // Includes all atoms in a given range around previously set markers
      // into the current QM region.
      // Example:
      // Includes atoms in range 12.3 Angstrom around all currently marked atoms
      // into the current QM region.
      // include ATOMS by_range 12.3
      //
      auto
      includeAtomsByRange(Box&                            box            ,
                          const unsigned int&             markerNo       ,
                          const unsigned int&             quantumRegionNo,
                          const std::vector<std::string>& command        ,
                          Logger&                         log             )const
      ->void;


      //
      // Includes a molecule into the current QM region.
      // Example:
      // Includes molecule based on atom 42 into current QM region.
      // include MOLECULE 42
      //
      auto
      includeMolecule(Box&                            box            ,
                      const unsigned int&             quantumRegionNo,
                      const std::vector<std::string>& command        ,
                      Logger&                         log             )const
      ->void;


      //
      // Includes all molecules in a given range around previously set markers
      // into the current QM region.
      // Example:
      // Includes moleculess in range 12.3 Angstrom around all currently marked atoms
      // into the current QM region.
      // include MOLECULES by_range 12.3
      //
      auto
      includeMoleculesByRange(Box&                            box            ,
                              const unsigned int&             markerNo       ,
                              const unsigned int&             quantumRegionNo,
                              const std::vector<std::string>& command        ,
                              Logger&                         log             )const
      ->void;


      //
      // Sets a marker on an atom to be processed by algorithms.
      // Example:
      // Marks atom 42.
      // mark ATOM 42
      //
      auto
      markAtom(Box&                            box       ,
               std::vector<unsigned int>&      nMarkerNos,
               unsigned int&                   markerNo  ,
               const std::vector<std::string>& command   ,
               Logger&                         log        )const   
      ->void;


      //
      // Marks up a Zundel ion by HIT types.
      // Usage:
      // mark ZUNDELION by_hit_types
      //
      auto
      markZundelIonByHitTypes(Box&                       box       ,
                              std::vector<unsigned int>& nMarkerNos,
                              unsigned int&              markerNo  ,
                              Logger&                    log        )const
      ->void;


      //
      // Modifies gradient of an atom by given amount.
      // Example:
      // Adds 1000.0 E -23 J to each component of the gradient of atom 42.
      // modify GRADIENT 42 1000.0 1000.0 1000.0
      //
      auto
      modifyGradient(Box&                            box    ,
                     const std::vector<std::string>& command,
                     Logger&                         log     )const
      ->void;


      //
      // Moves on to the next QM region.
      // Usage:
      // new QUANTUMREGION
      //
      auto
      newQuantumRegion(std::vector<unsigned int>& nQuantumRegions,
                       unsigned int&              quantumRegionNo,
                       Logger&                    log             )const
      ->void;


      //
      // Resets all markers.
      // Usage:
      // reset MARKERS
      //
      auto
      resetMarkers(Box&                       box       ,
                   std::vector<unsigned int>& nMarkerNos,
                   unsigned int&              markerNo  ,
                   Logger&                    log        )const
      ->void;


      //
      // Seeds automatically determined calculation of the current QM region.
      // Usage:
      // seed AUTO range distanceFactor
      //
      auto
      seedAuto(const Box&                               box            ,
               std::map< unsigned int          ,
                         std::unique_ptr<IfEvb> >&      nIfEvbs        ,
               std::map< unsigned int          ,
                         std::unique_ptr<IfUff> >&      nIfUffs        ,
               std::map< unsigned int               ,
                         std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs   ,
               std::map< unsigned int              ,
                         std::unique_ptr<IfOpenUff> >&  nIfOpenUffs    ,
               std::map< unsigned int          ,
                         std::unique_ptr<IfOrc> >&      nIfOrcs        ,
               const unsigned int                       quantumRegionNo,
               const std::vector<std::string>&          command        ,
               Logger&                                  log             )const
      ->void;


      //
      // Seeds an EVB calculation of the current QM region.
      // Usage:
      // seed EVB range distanceFactor
      //
      auto
      seedEvb(const Box&                          box            ,
              std::map< unsigned int,
                        std::unique_ptr<IfEvb> >& nIfEvbs        ,
              const unsigned int                  quantumRegionNo,
              const std::vector<std::string>&     command        ,
              Logger&                             log              )const
      ->void;


      //
      // Seeds a force field calculation of the current QM region.
      // Usage:
      // seed FF range distanceFactor
      //
      auto
      seedFf(const Box&                               box            ,
             std::map< unsigned int          ,
                       std::unique_ptr<IfUff> >&      nIfUffs        ,
             std::map< unsigned int               ,
                       std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs   ,
             std::map< unsigned int              ,
                       std::unique_ptr<IfOpenUff> >&  nIfOpenUffs    ,
             const unsigned int                       quantumRegionNo,
             const std::vector<std::string>&          command        ,
             Logger&                                  log             )const
      ->void;


      //
      // Seeds a quantumchemical calculation of the current QM region.
      // Usage:
      // seed QM range distanceFactor
      //
      auto
      seedQm(const Box&                          box            ,
             std::map< unsigned int,
                       std::unique_ptr<IfOrc> >& nIfOrcs        ,
             const unsigned int                  quantumRegionNo,
             const std::vector<std::string>&     command        ,
             Logger&                             log             )const
      ->void;


      //
      // Disables buffer region.
      // Usage:
      // set BUFFER inactive
      //
      auto
      setBufferInactive(Logger& log)
      ->void;
        

      //
      // Sets the inner radius of the spherical buffer region.
      // Example:
      // Sets inner radius to 5.0 Angstrom.
      // set BUFFER inner_radius 5.0
      //
      auto
      setBufferInnerRadius(const std::vector<std::string>& command,
                           Logger&                         log     )
      ->void;


      //
      // Sets the outer radius of the spherical buffer region.
      // Example:
      // Sets outer radius to 5.0 Angstrom.
      // set BUFFER outer_radius 5.0
      //
      auto
      setBufferOuterRadius(const std::vector<std::string>& command,
                           Logger&                         log     )
      ->void;


      //
      // Configures the buffer to use a linear fit function.
      // Usage:
      // set BUFFER using_linear_fit
      //
      auto
      setBufferUsingLinearFit(Logger& log)
      ->void;

      
      //
      // Configures the buffer to use a smootherstep function.
      // Usage:
      // set BUFFER using_smootherstep
      //
      auto
      setBufferUsingSmootherStep(Logger& log)
      ->void;


      //
      // Configures the buffer to use a smoothstep function.
      // Usage:
      // set BUFFER using_smoothstep
      //
      auto
      setBufferUsingSmoothStep(Logger& log)
      ->void;


      //
      // Disables EVB functionality.
      // Usage:
      // set EVB inactive
      //
      auto
      setEvbInactive(Logger& log)
      ->void;


      //
      // Sets the EVB program to be used as EVB implementation.
      // Usage:
      // set EVB using_evb
      //
      auto
      setEvbUsingEvb(Logger& log)
      ->void;


      //
      // Disables force field functionality.
      // Usage:
      // set FF inactive
      //
      auto
      setFfInactive(Logger& log)
      ->void;


      //
      // Sets the UFF program to be used as force field implementation.
      // Usage:
      // set FF using_uff
      //
      auto
      setFfUsingUff(Logger& log)
      ->void;


      //
      // Sets OpenGAFF to be used as force field implementation.
      // Usage:
      // set FF using_opengaff
      //
      auto
      setFfUsingOpenGaff(Logger& log)
      ->void;


      //
      // Sets OpenUFF to be used as force field implementation.
      // Usage:
      // set FF using_openuff
      //
      auto
      setFfUsingOpenUff(Logger& log)
      ->void;


      //
      // Enable standalone functionality.
      // Usage:
      // set GEOMETRY using_config_file
      //
      auto
      setGeometryUsingConfigFile(Box&    box,
                                 Logger& log )
      ->void;


      //
      // Switch control to MD frontend.
      // Usage:
      // set GEOMETRY using_md_frontend
      //
      auto
      setGeometryUsingMdFrontend(Logger& log)
      ->void;


      //
      // Sets the gradient of an atom manually.
      // Example:
      // Sets each component of the gradient of atom 42 to 1000.0.
      // set GRADIENT 42 1000.0 1000.0 1000.0
      //
      auto
      setGradient(Box&                            box    ,
                  const std::vector<std::string>& command,
                  Logger&                         log     )const
      ->void;


      //
      // Enables heat up period.
      // Example:
      // Heats up the system by 0.1 Kelvin per MD step.
      // set HEATUP 0.1
      //
      auto
      setHeatup(Box&                            box    ,
                const std::vector<std::string>& command,
                Logger&                         log     )const
      ->void;


      //
      // Disable external MD functionality.
      // Usage:
      // set MD inactive
      //
      auto
      setMdInactive(Logger& log)
      ->void;


      //
      // Sets the HIT program to be used as external MD implementation.
      // Usage:
      // set MD using_hit
      //
      auto
      setMdUsingHit(Logger& log)
      ->void;


      //
      // Disable neighbor list scheme for calculating distances.
      // Usage:
      // set NEIGHBORLISTS inactive
      //
      auto
      setNeighborListsInactive(Logger& log)
      ->void;


      //
      // Sets number of MD steps after which all distances are being reevaluated.
      // Example:
      // Enables recalculation of all distances every two steps.
      // set NEIGHBORLISTS step_count 2
      //
      auto
      setNeighborListsStepCount(const std::vector<std::string>& command,
                                Logger&                         log     )
      ->void;


      //
      // Sets number of nearest neighbor atoms to be included in
      // intermediary distance calculations.
      // Example:
      // Enables intermediary realculation of distances concerning
      // 16 nearest neighbors of each atom.
      // set NEIGHBORLISTS step_size 16
      //
      auto
      setNeighborListsStepSize(const std::vector<std::string>& command,
                               Logger&                         log     )
      ->void;


      //
      // Disable quantummechanical functionality.
      // Usage:
      // set QM inactive
      //
      auto
      setQmInactive(Logger& log)
      ->void;


      //
      // Sets the ORC program to be used as QM implementation.
      // Usage:
      // set QM using_orc
      //
      auto
      setQmUsingOrc(Logger& log)
      ->void;


      //
      // Sets the temperature and Berendsen dampening factor for MD simulation.
      // Setting the dampening factor to MD time frame results in simple velocity scaling.
      // Setting the dampening factor to higher values enables a Berendsen Thermostat;
      // its effect is fading away when the dampening factor reaches infinity.
      // Example:
      // Set temperature to 298.15 Kelvin and the Berendsen dampening factor to 400 Femtoseconds.
      // set THERMOSTAT 298.15 400
      //
      auto
      setThermostat(Box&                            box    ,
                    const std::vector<std::string>& command,
                    Logger&                         log     )const
      ->void;


      //
      // Sets the time frame for MD simulation.
      // Example:
      // Proceed MD simulation with 1.5 Femtosecond steps.
      // set TIMEFRAME 1.5
      //
      auto
      setTimeFrame(Box&                            box    ,
                   const std::vector<std::string>& command,
                   Logger&                         log     )const
      ->void;


      //
      // Updates the trajectory file.
      // Usage:
      // update TRAJECTORY
      //
      auto
      updateTrajectory(Box&    box,
                       Logger& log )const
      ->void;



      //
      // Locking and Unlock Qreg for thread safety.
      //
      std::mutex threadSafety_;


      //
      // Stores "true" or "false" whether an internal geometry from a config file is used or not.
      //
      bool internalGeometry_;


      //
      // Stores the choice of gradient engine.
      //
      unsigned int autoToolChoice_;

      //
      // Stores the choice of MD program.
      //
      unsigned int mdFrontend_;

      //
      // Stores the choice of FF program.
      //
      unsigned int ffBackend_;


      //
      // Stores the choice of QM program.
      //
      unsigned int qmBackend_;

      //
      // Stores the choice of EVB program.
      //
      unsigned int evbBackend_;

      //
      // Stores the run number.
      //
      unsigned int runNo_;

      //
      // Stores the number of steps between reevaluation of ALL distances.
      //
      unsigned int stepCount_;

      //
      // Stores the coice of buffer type.
      //
      unsigned int bufferType_;


      //
      // Stores total charge.
      //
      double autoCharge_;

      //
      // Stores the inner radius to be used for computing the buffer region.
      //
      double innerBufferRadius_;

      //
      // Stores the outer radius to be used for computing the buffer region.
      //
      double outerBufferRadius_;

      //
      // Stores the number of neighbor atoms considered when computing nearest neighbor distances.
      //
      double stepSize_;

      //
      // Stores the configuration options as an instance of the Config class.
      //
      Config configuration_;



  };



  inline auto
  Qreg::addAttraction(Box&                            box    ,
                      const std::vector<std::string>& command,
                      Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 5);


    //
    // inline auto
    // Qreg::addAttraction(Box&                            box    ,
    //                     const std::vector<std::string>& command,
    //                     Logger&                         log     )const
    // ->void
    //

    std::string message;

    if((box.getXGradient(std::stoi(command[3])) != 0.0) &&
       (box.getYGradient(std::stoi(command[3])) != 0.0) &&
       (box.getZGradient(std::stoi(command[3])) != 0.0)   )
      box.setGradient(std::stoi(command[3]),
                      std::stod(command[4]) *
       (box.getXCoordinate(std::stoi(command[3])) -
        box.getXCoordinate(std::stoi(command[2]))  )
       + box.getXGradient(std::stoi(command[3]))    ,
                      std::stod(command[4]) *
       (box.getYCoordinate(std::stoi(command[3])) -
        box.getYCoordinate(std::stoi(command[2]))  )
       + box.getYGradient(std::stoi(command[3]))    ,
                      std::stod(command[4]) *
       (box.getZCoordinate(std::stoi(command[3])) -
        box.getZCoordinate(std::stoi(command[2]))  )
       + box.getZGradient(std::stoi(command[3]))     );

    box.setGradientGuess(std::stoi(command[3]),
                         std::stod(command[4]) *
     (box.getXCoordinate(std::stoi(command[3])) -
      box.getXCoordinate(std::stoi(command[2]))  )
     + box.getXGradientGuess(std::stoi(command[3])),
                         std::stod(command[4]) *
     (box.getYCoordinate(std::stoi(command[3])) -
      box.getYCoordinate(std::stoi(command[2]))  )
     + box.getYGradientGuess(std::stoi(command[3])),
                         std::stod(command[4]) *
     (box.getZCoordinate(std::stoi(command[3])) -
      box.getZCoordinate(std::stoi(command[2]))  )
     + box.getZGradientGuess(std::stoi(command[3])) );

    message.append("Adding attraction of atom number ");
    message.append(command[3]                         );
    message.append(" to atom number "                 ); 
    message.append(command[2]                         );
    message.append(" by "                             );
    message.append(command[4]                         );
    message.append(" E -23 Joule per Angstrom."       );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::advanceMd(Box&    box,
                  Logger& log )const
  ->void
  {


    //
    // inline auto
    // Qreg::advanceMd(Box&    box,
    //                 Logger& log )const
    // ->void
    //

    std::string message;

    message.append("Advancing MD simulation ...");
    log.message("MD_DRIVER",
                message     );
    message.clear();

    message.append("Using MD step size of "          );
    message.append(std::to_string(box.getTimeFrame()));
    message.append(" Femtoseconds."                   );
    log.message("MD_DRIVER",
                message     );
    message.clear();

    #ifndef NDEBUG

    #ifdef DEBUG_MD_DETAILS
            
      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                          );
        message.append("Previous coordinate of atom ");
        message.append(std::to_string(atomNo)        );
        message.append(", "                          );
        message.append(box.getChemicalSymbol(atomNo) );
        message.append(": "                          );
        message.append
         (std::to_string(box.getXCoordinate(atomNo)) );
        message.append(" "                           );
        message.append
         (std::to_string(box.getYCoordinate(atomNo)) );
        message.append(" "                           );
        message.append
        (std::to_string(box.getZCoordinate(atomNo))  );
        message.append(" Angstrom"                   );

      }

      log.message("OLD_COORDINATES",
                  message  );
      message.clear();

      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                          );
        message.append("Previous velocity of atom "  );
        message.append(std::to_string(atomNo)        );
        message.append(", "                          );
        message.append(box.getChemicalSymbol(atomNo) );
        message.append(": "                          );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[0]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[1]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[2]));
        message.append(" Angstrom per Femtosecond"    );

      }

      log.message("OLD_VELOCITIES",
                  message          );
      message.clear();

      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                          );
        message.append("Force on atom "              );
        message.append(std::to_string(atomNo)        );
        message.append(", "                          );
        message.append(box.getChemicalSymbol(atomNo) );
        message.append(": "                          );
        message.append
         (std::to_string(box.computeForce(atomNo)[0]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeForce(atomNo)[1]));
        message.append(" "                           );
        message.append
        (std::to_string(box.computeForce(atomNo)[2]) );
        message.append(" E-23 J per Angstrom"        );

      }

      log.message("FORCES",
                  message  );
      message.clear();

    #endif

    #ifdef DEBUG_MD_ACCELERATION

      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                          );
        message.append("Acceleration on atom "       );
        message.append(std::to_string(atomNo)        );
        message.append(", "                          );
        message.append(box.getChemicalSymbol(atomNo) );
        message.append(": "                          );
        message.append
         (std::to_string(box.computeAcceleration(atomNo)[0]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeAcceleration(atomNo)[1]));
        message.append(" "                           );
        message.append
        (std::to_string(box.computeAcceleration(atomNo)[2]) );
        message.append(" Angstrom per Femtosecond squared");

      } 

      log.message("ACCELERATION",
                  message  );
      message.clear();

    #endif

    #endif

    const auto temperature = box.updateCoordinates();

    message.append("Berendsen Thermostat reports temperature of ");
    message.append(std::to_string(temperature)                   );
    message.append(" Kelvin."                                    );
    log.message("THERMOSTAT",
                message      );
      message.clear();

    #ifndef NDEBUG

    #ifdef DEBUG_MD_VELOCITIES

      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                          );
        message.append("Current velocity of atom  "  );
        message.append(std::to_string(atomNo)        );
        message.append(", "                          );
        message.append(box.getChemicalSymbol(atomNo) );
        message.append(": "                          );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[0]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[1]));
        message.append(" "                           );
        message.append
         (std::to_string(box.computeVelocityVector(atomNo)[2]));
        message.append(" Angstrom per Femtosecond"    );

      }

      log.message("NEW_VELOCITIES",
                  message          );
      message.clear();

    #endif

    #ifdef DEBUG_MD_DETAILS

      for(const auto& atomNo:box.getNAtomNos())
      {

        message.append("\n"                         );
        message.append("Current coordinate of atom ");
        message.append(std::to_string(atomNo)       );
        message.append(", "                         );
        message.append(box.getChemicalSymbol(atomNo));
        message.append(": "                         );
        message.append
         (std::to_string(box.getXCoordinate(atomNo)));
        message.append(" "                          );
        message.append
         (std::to_string(box.getYCoordinate(atomNo)));
        message.append(" "                          );
        message.append
        (std::to_string(box.getZCoordinate(atomNo)) );
        message.append(" Angstrom"                  );

      }

      log.message("NEW_COORDINATES",
                  message  );
      message.clear();

    #endif

    #endif

    message.append("Advancing MD simulation DONE");
    log.message("MD_DRIVER",
                message     );
    message.clear();


  }



  inline auto
  Qreg::computeBuffer(Box&                box     ,
                      const unsigned int& markerNo,
                      Logger&             log      )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(markerNo != 0);


    //
    // inline auto
    // Qreg::computeBuffer(Box&                box     ,
    //                     const unsigned int& markerNo,
    //                     Logger&             log      )const
    // ->void
    //

    if((bufferType_ == QREG_LINEARFIT   ) ||
       (bufferType_ == QREG_SMOOTHSTEP  ) ||
       (bufferType_ == QREG_SMOOTHERSTEP)   )
    {

      std::string message;

      message.append("Computing gradients in the buffer region ...");
      log.message("PARTITIONER",
                  message       );
      message.clear();

      if(bufferType_ == 0)
      {
        message.append("Buffering scheme unavailable.");
        log.message("PARTITIONER",
                    message       );
        message.clear();
      }

      else if(bufferType_ == QREG_LINEARFIT)
        Partitioner::computeLinearBuffer(box               ,
                                         markerNo          ,
                                         innerBufferRadius_,
                                         outerBufferRadius_ );
      else if(bufferType_ == QREG_SMOOTHSTEP)
        Partitioner::computeSmoothStepBuffer(box               ,
                                             markerNo          ,
                                             innerBufferRadius_,
                                             outerBufferRadius_ );
      else if(bufferType_ == QREG_SMOOTHERSTEP)
        Partitioner::computeSmootherStepBuffer(box               ,
                                               markerNo          ,
                                               innerBufferRadius_,
                                               outerBufferRadius_ );

      message.append("Computing gradients in the buffer region DONE");
      log.message("PARTITIONER",
                  message       );
      message.clear();

    }


  }



  inline auto
  Qreg::computeDistances(Box&    box,
                         Logger& log )const
  ->void
  {


    //
    // inline auto
    // Qreg::computeDistances(Box&    box,
    //                        Logger& log )const
    // ->void
    //

    std::string message;

    if((runNo_ % stepCount_ == 0) ||
       (runNo_              == 1)   )
    {

      message.append("Evaluating all distances ...");
      log.message("QREG"  ,
                   message );
      message.clear();
      box.computeAllDistances();

      #ifndef NDEBUG
      #ifdef DEBUG_DISTANCES

        const auto nAtomNos   = box.getNAtomNos();
        const auto nTargetNos = nAtomNos         ;

        message.append("TESTING DISTANCES");

        for(const auto& targetId:nTargetNos)
        {

          for(const auto& atomNo:nAtomNos)
          {

            if(targetId != atomNo)
            {
              message.append("\n"                         );
              message.append("DIST "                      );
              message.append(std::to_string(targetId)     );
              message.append(" to "                       );
              message.append(std::to_string(atomNo)       );
              message.append(": "                         );
              message.append
               (std::to_string(box.getDistance(targetId,
                                               atomNo   )));
            }

          }

        }

        log.message("DEBUG",
                    message );
        message.clear();

      #endif
      #endif

      message.append("Reevaluating all distances DONE");
      log.message("QREG" ,
                  message );
      message.clear();

    }

    else
    {

      message.append("Reevaluating distances between nearest neighbors up to ");
      message.append(std::to_string(stepSize_)                                );
      message.append(" Angstrom apart ..."                                    );
      log.message("QREG" ,
                  message );
      message.clear();

      box.computeNearestNeighborDistances(stepSize_);

      #ifndef NDEBUG
      #ifdef DEBUG_DISTANCES

        const auto nAtomNos   = box.getNAtomNos();
        const auto nTargetNos = nAtomNos         ;

        message.append("TESTING DISTANCES");

        for(const auto& targetId:nTargetNos)
        {

          for(const auto& atomNo:nAtomNos)
          {

            if(targetId != atomNo)
            {
              message.append("\n"                         );
              message.append("DIST "                      );
              message.append(std::to_string(targetId)     );
              message.append(" to "                       );
              message.append(std::to_string(atomNo)       );
              message.append(": "                         );
              message.append
               (std::to_string(box.getDistance(targetId,
                                               atomNo   )));
            }

          }

        }

        log.message("DEBUG",
                    message );
        message.clear();

      #endif
      #endif

      message.append("Reevaluating distances between nearest neighbors up to ");
      message.append(std::to_string(stepSize_)                                );
      message.append(" Angstrom apart DONE"                                   );
      log.message("QREG" ,
                  message );
      message.clear();

    }


  }



  inline auto
  Qreg::harvestAuto(Box&                                     box         ,
                    std::map< unsigned int          ,
                              std::unique_ptr<IfEvb> >&      nIfEvbs     ,
                    std::map< unsigned int          ,
                              std::unique_ptr<IfUff> >&      nIfUffs     ,
                    std::map< unsigned int               ,
                              std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs,
                    std::map< unsigned int              ,
                              std::unique_ptr<IfOpenUff> >&  nIfOpenUffs ,
                    std::map< unsigned int          ,
                              std::unique_ptr<IfOrc> >&      nIfOrcs     ,
                    Logger&                                  log          )const
  ->void
  {


    switch(autoToolChoice_)
    {


      case QREG_EVB:
        harvestEvb(box    ,
                   nIfEvbs,
                   log     );
        break;

      case QREG_TFF:
        harvestFf(box         ,
                  nIfUffs     ,
                  nIfOpenGaffs,
                  nIfOpenUffs ,
                  log          );
        break;

      case QREG_GFF:
        harvestFf(box         ,
                  nIfUffs     ,
                  nIfOpenGaffs,
                  nIfOpenUffs ,
                  log          );
        break;

      case QREG_UFF:
        harvestFf(box         ,
                  nIfUffs     ,
                  nIfOpenGaffs,
                  nIfOpenUffs ,
                  log          );
        break;

      case QREG_ORC:
        harvestQm(box    ,
                  nIfOrcs,
                  log     );
        break;

      default:
        harvestQm(box    ,
                  nIfOrcs,
                  log     );
        break;

    } 


  }



  inline auto
  Qreg::harvestEvb(Box&                                box    ,
                   std::map< unsigned int,
                             std::unique_ptr<IfEvb> >& nIfEvbs,
                   Logger&                             log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(!nIfEvbs.empty());


    //
    // inline auto
    // Qreg::harvestEvb(Box&                                box    ,
    //                  std::map< unsigned int,
    //                            std::unique_ptr<IfEvb> >& nIfEvbs,
    //                  Logger&                             log     )const
    // ->void
    //

    std::string message;

    message.append("Retrieving EVB gradient ...");
    log.message("QREG" ,
                message );
    message.clear();

    if(evbBackend_ == QREG_EVB)
    {

      for(auto&& ifEvb:nIfEvbs)
      {

        if(ifEvb.second->harvest(box) == true)
        {

          message.append("The gradient of QM region "            );
          message.append(std::to_string(ifEvb.first)             );
          message.append
           (" has been retrieved from an EVB calculation. Epot = ");
          message.append(std::to_string(box.getEPot(ifEvb.first)));
          message.append(" E -23 Joule."                         );
          log.message("IFEVB",
                      message );
          message.clear();

        }

        else
        {

          message.append("Abnormal termination of EVB reported. Dithering gradient.");
          log.message("IFEVB",
                      message );
          message.clear();

        }

        for(auto&& atomNo:box.getNAtomNosByQuantumRegionNo(ifEvb.first))
        {

          message.append("\n"                                    );
          message.append("Gradient of atom "                     );
          message.append(std::to_string(atomNo)                  );
          message.append(", "                                    );
          message.append(box.getChemicalSymbol(atomNo)           );
          message.append(": "                                    );
          message.append(std::to_string(box.getXGradient(atomNo)));
          message.append("   "                                   );
          message.append(std::to_string(box.getYGradient(atomNo)));
          message.append("   "                                   );
          message.append(std::to_string(box.getZGradient(atomNo)));
          message.append("   E -23 Joule per Angstrom"           );

        }

        log.message("GRADIENT",
                    message    );
        message.clear();

      }

    }

    message.append("Retrieving EVB gradient DONE");
    log.message("QREG" ,
                 message );
    message.clear();


  }



  inline auto
  Qreg::harvestFf(Box&                                     box         ,
                  std::map< unsigned int          ,
                            std::unique_ptr<IfUff> >&      nIfUffs     ,
                  std::map< unsigned int               ,
                            std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs,
                  std::map< unsigned int              ,
                            std::unique_ptr<IfOpenUff> >&  nIfOpenUffs ,
                  Logger&                                  log          )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(!((nIfUffs     .empty()) &&
             (nIfOpenGaffs.empty()) &&
             (nIfOpenUffs .empty())   ));

    //
    // inline auto
    // Qreg::harvestFf(Box&                                     box         ,
    //                 std::map< unsigned int          ,
    //                           std::unique_ptr<IfUff> >&      nIfUffs     ,
    //                 std::map< unsigned int               ,
    //                           std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs,
    //                 std::map< unsigned int              ,
    //                           std::unique_ptr<IfOpenUff> >&  nIfOpenUffs ,
    //                 Logger&                                  log          )const
    // ->void

    std::string message;

    message.append("Retrieving FF gradient ...");
    log.message("QREG" ,
                message );
    message.clear();

    if(ffBackend_ == QREG_TFF)
    {

      for(auto&& ifUff:nIfUffs)
      {

        if(ifUff.second->harvest(box) == true)
        {

          message.append("The gradient of QM region "            );
          message.append(std::to_string(ifUff.first)             );
          message.append
           (" has been retrieved from an UFF calculation. Epot = ");
          message.append(std::to_string(box.getEPot(ifUff.first)));
          message.append(" E -23 Joule."                         );
          log.message("IFUFF",
                      message );
          message.clear();

        }

        else
        {

          message.append("Abnormal termination of UFF reported. Dithering gradient.");
          log.message("IFUFF",
                      message );
          message.clear();

        }

        for(auto&& atomNo:box.getNAtomNosByQuantumRegionNo(ifUff.first))
        {

          message.append("\n"                                         );
          message.append("Gradient of atom "                          );
          message.append(std::to_string(atomNo)                       );
          message.append(", "                                         );
          message.append(box.getChemicalSymbol(atomNo)                );
          message.append(": "                                         );
          message.append(std::to_string(box.getXGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getYGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getZGradientGuess(atomNo)));
          message.append("   E -23 Joule per Angstrom"                );

        }

        log.message("GRADIENT",
                    message    );
        message.clear();

      }

    }

    if(ffBackend_ == QREG_GFF)
    {
      
      for(auto&& ifOpenGaff:nIfOpenGaffs)
      { 
        
        if(ifOpenGaff.second->harvest(box) == true)
        { 
          
          message.append("The gradient of QM region "          );
          message.append(std::to_string(ifOpenGaff.first)      );
          message.append
           (" has been retrieved from an OpenGAFF calculation.");
          log.message("IFOPENGAFF",
                      message      );
          message.clear();
        
        }
        
        else
        { 
          
          message.append("Abnormal termination of OpenGAFF reported. Dithering gradient.");
          log.message("IFOPENUFF",
                      message     );
          message.clear();
        
        }
        
        for(auto&& atomNo:box.getNAtomNosByQuantumRegionNo(ifOpenGaff.first))
        { 
          
          message.append("\n"                                         );
          message.append("Gradient of atom "                          );
          message.append(std::to_string(atomNo)                       );
          message.append(", "                                         );
          message.append(box.getChemicalSymbol(atomNo)                );
          message.append(": "                                         );
          message.append(std::to_string(box.getXGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getYGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getZGradientGuess(atomNo)));
          message.append("   E -23 Joule per Angstrom"                );
        
        }
        
        log.message("GRADIENT",
                    message    );
        message.clear();
        
      }
      
    }

    if(ffBackend_ == QREG_UFF)
    {

      for(auto&& ifOpenUff:nIfOpenUffs)
      {

        if(ifOpenUff.second->harvest(box) == true)
        {

          message.append("The gradient of QM region "         );
          message.append(std::to_string(ifOpenUff.first)      );
          message.append
           (" has been retrieved from an OpenUFF calculation.");
          log.message("IFOPENUFF",
                      message      );
          message.clear();

        }

        else
        {

          message.append("Abnormal termination of OpenUFF reported. Dithering gradient.");
          log.message("IFOPENUFF",
                      message     );
          message.clear();

        }

        for(auto&& atomNo:box.getNAtomNosByQuantumRegionNo(ifOpenUff.first))
        {

          message.append("\n"                                         );
          message.append("Gradient of atom "                          );
          message.append(std::to_string(atomNo)                       );
          message.append(", "                                         );
          message.append(box.getChemicalSymbol(atomNo)                );
          message.append(": "                                         );
          message.append(std::to_string(box.getXGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getYGradientGuess(atomNo)));
          message.append("   "                                        );
          message.append(std::to_string(box.getZGradientGuess(atomNo)));
          message.append("   E -23 Joule per Angstrom"                );

        }

        log.message("GRADIENT",
                    message    );
        message.clear();

      }

    }

    message.append("Retrieving FF gradient DONE");
    log.message("QREG" ,
                 message );
    message.clear();


  }



  inline auto
  Qreg::harvestQm(Box&                                box    ,
                  std::map< unsigned int,
                            std::unique_ptr<IfOrc> >& nIfOrcs,
                  Logger&                             log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(!nIfOrcs.empty());


    //
    // inline auto
    // Qreg::harvestQm(Box&                                box    ,
    //                 std::map< unsigned int,
    //                           std::unique_ptr<IfOrc> >& nIfOrcs,
    //                 Logger&                             log     )const
    // ->void
    //

    std::string message;

    message.append("Retrieving QM gradient ...");
    log.message("QREG" ,
                message );
    message.clear();

    if(qmBackend_ == QREG_ORC)
    {

      for(auto&& ifOrc:nIfOrcs)
      {

        if(ifOrc.second->harvest(box) == true)
        {

          message.append("The gradient of QM region "            );
          message.append(std::to_string(ifOrc.first)             );
          message.append
           (" has been retrieved from an ORC calculation. Epot = ");
          message.append(std::to_string(box.getEPot(ifOrc.first)));
          message.append(" E -23 Joule."                         );
          log.message("IFORC",
                      message );
          message.clear();

        }

        else
        {

          message.append("Abnormal termination of ORC reported. Dithering gradient.");
          log.message("IFORC",
                      message );
          message.clear();

        }

        for(auto&& atomNo:box.getNAtomNosByQuantumRegionNo(ifOrc.first))
        {

          message.append("\n"                                    );
          message.append("Gradient of atom "                     );
          message.append(std::to_string(atomNo)                  );
          message.append(", "                                    );
          message.append(box.getChemicalSymbol(atomNo)           );
          message.append(": "                                    );
          message.append(std::to_string(box.getXGradient(atomNo)));
          message.append("   "                                   );
          message.append(std::to_string(box.getYGradient(atomNo)));
          message.append("   "                                   );
          message.append(std::to_string(box.getZGradient(atomNo)));
          message.append("   E -23 Joule per Angstrom"           );

        }

        log.message("GRADIENT",
                    message    );
        message.clear();

        }

      }


      message.append("Retrieving QM gradient DONE");
      log.message("QREG" ,
                   message );
      message.clear();


  }



  inline auto
  Qreg::includeAtom(Box&                            box            ,
                    const unsigned int&             quantumRegionNo,
                    const std::vector<std::string>& command        ,
                    Logger&                         log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 3);


    //
    // inline auto
    // Qreg::includeAtom(Box&                            box            ,
    //                   const unsigned int&             quantumRegionNo,
    //                   const std::vector<std::string>& command        ,
    //                   Logger&                         log             )const
    // ->void
    //

    std::string message;

    message.append("Including atom ");
    message.append(command[2]);
    message.append(" ..."    );
    log.message("QREG" ,
                message );
    message.clear();


    box.insertQuantumRegionNo(std::stoi(command[2]),
                              quantumRegionNo       );

    message.append("Atom "                               );
    message.append(command[2]                            );
    message.append(" has been included in the QM region ");
    message.append(std::to_string(quantumRegionNo)       );
    message.append("."                                   );
    log.message("QREG" ,
                message );
    message.clear(); 



  }



  inline auto
  Qreg::includeAtomsAll(Box&                box            ,
                        const unsigned int& quantumRegionNo,
                        Logger&             log             )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(quantumRegionNo != 0);


    //
    // inline auto
    // Qreg::includeAtomsAll(Box&                box            ,
    //                       const unsigned int& quantumRegionNo,
    //                       Logger&             log             )const
    // ->void
    //

    std::string message;

    message.append("Including all atoms ...");
    log.message("QREG" ,
                message );
    message.clear();

    const auto nAtomNos = box.getNAtomNos();

    message.append("Atoms ");

    for(const auto& atomNo:nAtomNos)
    {

      box.insertQuantumRegionNo(atomNo         ,
                                quantumRegionNo );

      message.append(std::to_string(atomNo));
      message.append(" "                   );

    }

    message.append("have been included in the QM region ");
    message.append(std::to_string(quantumRegionNo)       );
    message.append("."                                   );
    log.message("QREG" ,
                message );
    message.clear();

    message.append("Including all atoms DONE");
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::includeAtomsByRange(Box&                            box            ,
                            const unsigned int&             markerNo       ,
                            const unsigned int&             quantumRegionNo,
                            const std::vector<std::string>& command        ,
                            Logger&                         log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 4);


    //
    // inline auto
    // Qreg::includeAtomsByRange(Box&                            box            ,
    //                          const unsigned int&             markerNo       ,
    //                          const unsigned int&             quantumRegionNo,
    //                          const std::vector<std::string>& command        ,
    //                          Logger&                         log             )const
    // ->void
    //

    std::string message;

    std::vector<size_t> nAtomNos;

    message.append("Including atoms in range ");
    message.append(command[3]                 );
    message.append(" Angstrom ..."            );
    log.message("QREG" ,
                message );
    message.clear();

    for(const auto& atomNo:box.getNAtomNosByMarkerNo(markerNo))
    {

      auto nNewAtomNos = box.getNNeighborAtomNos(atomNo               ,
                                                 std::stod(command[3]) );

      nNewAtomNos.push_back(atomNo);

      nAtomNos.insert(nAtomNos   .end  (),
                      nNewAtomNos.begin(),
                      nNewAtomNos.end  () );

      std::sort(nAtomNos.begin(),
                nAtomNos.end  () );

      nAtomNos.erase(std::unique(nAtomNos.begin(),
                                 nAtomNos.end  () ),
                     nAtomNos.end()                 );

    }

    message.append("Atoms ");
    for(const auto& atomNo:nAtomNos)
    {

      box.insertQuantumRegionNo(atomNo         ,
                                quantumRegionNo );

      message.append(std::to_string(atomNo));
      message.append(" "                   );

    }

    message.append("have been included in the QM region ");
    message.append(std::to_string(quantumRegionNo)       );
    message.append("."                                   );
    log.message("QREG" ,
                message );
    message.clear();

    message.append("Including atoms in range ");
    message.append(command[3]                 );
    message.append(" Angstrom DONE"           );        
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::includeMolecule(Box&                            box            ,
                        const unsigned int&             quantumRegionNo,
                        const std::vector<std::string>& command        ,
                        Logger&                         log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 4);


    //
    // inline auto
    // Qreg::includeMolecule(Box&                            box            ,
    //                       const unsigned int&             quantumRegionNo,
    //                       const std::vector<std::string>& command        ,
    //                       Logger&                         log             )const
    // ->void
    //

    std::string message;

    std::vector<size_t> nAtomNos =
    {
      static_cast<unsigned>(std::stol(command[2]))
    };
    auto nAtomNosSize = nAtomNos.size();

    message.append("Including molecule based on atom ");
    message.append(command[2]                         );
    message.append(" using distance factor "          );
    message.append(command[3]                         );
    message.append(" ..."                             );
    log.message("QREG" ,
                message );
    message.clear();

    do
    {

      std::vector<size_t> nNewAtomNos;

      nAtomNosSize = nAtomNos.size();

      for(const auto& atomNo:nAtomNos)
      {

        for(const auto& neighborAtomNo:box.getNNeighborAtomNos(atomNo,
                                                               8.0    ))
        {

          if(box.computeDistance(atomNo        ,
                                 neighborAtomNo ) < (std::stod(command[3]))             *
                                                    (box.getVdwRadius(atomNo        ) +
                                                     box.getVdwRadius(neighborAtomNo)  ) )
            nNewAtomNos.push_back(neighborAtomNo);

        }

      }

      nAtomNos.insert(nAtomNos   .end  (),
                      nNewAtomNos.begin(),
                      nNewAtomNos.end  () );

      std::sort(nAtomNos.begin(),
                nAtomNos.end  () );

      nAtomNos.erase(std::unique(nAtomNos.begin(),
                                 nAtomNos.end  () ),
                     nAtomNos.end()                 );   

    } while(nAtomNos.size() != nAtomNosSize);

    message.append("Atoms ");
    for(const auto& atomNo:nAtomNos)
    {

      box.insertQuantumRegionNo(atomNo         ,
                                quantumRegionNo );

      message.append(std::to_string(atomNo));
      message.append(" "                   );

    }
    message.append("have been included in the QM region ");
    message.append(std::to_string(quantumRegionNo)       );
    message.append("."                                   );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::includeMoleculesByRange(Box&                            box            ,
                                const unsigned int&             markerNo       ,
                                const unsigned int&             quantumRegionNo,
                                const std::vector<std::string>& command        ,
                                Logger&                         log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 5);


    //
    // inline auto
    // Qreg::includeMoleculesByRange(Box&                            box            ,
    //                               const unsigned int&             markerNo       ,
    //                               const unsigned int&             quantumRegionNo,
    //                               const std::vector<std::string>& command        ,
    //                               Logger&                         log             )const
    // ->void
    //

    std::string message;

    std::vector<size_t> nAtomNos;

    message.append("Including molecules in range "   );
    message.append(command[3]                        );
    message.append(" Angstrom using distance factor ");
    message.append(command[4]                        );
    message.append(" ..."                            );
    log.message("QREG" ,
                message );
    message.clear();

    for(const auto& atomNo:box.getNAtomNosByMarkerNo(markerNo))
    {

      auto nNewAtomNos = box.getNNeighborAtomNos(atomNo               ,
                                                 std::stod(command[3]) );

      nNewAtomNos.push_back(atomNo);

      nAtomNos.insert(nAtomNos   .end  (),
                      nNewAtomNos.begin(),
                      nNewAtomNos.end  () );

      std::sort(nAtomNos.begin(),
                nAtomNos.end  () );

      nAtomNos.erase(std::unique(nAtomNos.begin(),
                                 nAtomNos.end  () ),
                     nAtomNos.end()                 );

    }

    auto nAtomNosSize = nAtomNos.size();

    do
    {

      std::vector<size_t> nNewAtomNos;

      nAtomNosSize = nAtomNos.size();

      for(const auto& atomNo:nAtomNos)
      {

        for(const auto& neighborAtomNo:box.getNNeighborAtomNos(atomNo,
                                                               8.0    ))
        {

          if(box.computeDistance(atomNo        ,
                                 neighborAtomNo ) < (std::stod(command[4]))             *
                                                    (box.getVdwRadius(atomNo        ) +
                                                     box.getVdwRadius(neighborAtomNo)  ) )
            nNewAtomNos.push_back(neighborAtomNo);

        }

      }

      nAtomNos.insert(nAtomNos   .end  (),
                      nNewAtomNos.begin(),
                      nNewAtomNos.end  () );

      std::sort(nAtomNos.begin(),
                nAtomNos.end  () );

      nAtomNos.erase(std::unique(nAtomNos.begin(),
                                 nAtomNos.end  () ),
                     nAtomNos.end()                 );

    } while(nAtomNos.size() != nAtomNosSize);

    message.append("Atoms ");
    for(const auto& atomNo:nAtomNos)
    {

      box.insertQuantumRegionNo(atomNo         ,
                                quantumRegionNo );

      message.append(std::to_string(atomNo));
      message.append(" "                   );

    }

    message.append("have been included in the QM region ");
    message.append(std::to_string(quantumRegionNo)       );
    message.append("."                                   );
    log.message("QREG" ,
                message );
    message.clear();

    message.append("Including molecules in range ");
    message.append(command[3]                     );
    message.append(" Angstrom DONE"               );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::markAtom(Box&                            box       ,
                 std::vector<unsigned int>&      nMarkerNos,
                 unsigned int&                   markerNo  ,
                 const std::vector<std::string>& command   ,
                 Logger&                         log        )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 3);


    //
    // inline auto
    // Qreg::markAtom(Box&                            box       ,
    //                std::vector<unsigned int>&      nMarkerNos,
    //                unsigned int&                   markerNo  ,
    //                const std::vector<std::string>& command   ,
    //                Logger&                         log        )const
    // ->void
    //

    std::string message;

    ++markerNo;
    nMarkerNos.push_back(markerNo);

    box.insertMarkerNo(std::stoi(command[2]),
                       markerNo              );

    message.append("Marker "               );
    message.append(std::to_string(markerNo));
    message.append(" has been set on Atom ");
    message.append(command[2]              );
    message.append("."                     );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::markZundelIonByHitTypes(Box&                       box       ,
                                std::vector<unsigned int>& nMarkerNos,
                                unsigned int&              markerNo  ,
                                Logger&                    log        )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(markerNo != 0);


    //
    // inline auto
    // Qreg::markZundelIonByHitTypes(Box&                       box       ,
    //                               std::vector<unsigned int>& nMarkerNos,
    //                               unsigned int&              markerNo  ,
    //                               Logger&                    log        )const
    // ->void
    //

    std::string message;

    std::vector<unsigned long int> zundelIon;

    ++markerNo;
    nMarkerNos.push_back(markerNo);

    message.append("Marking Zundel ion by HIT types ...");
    log.message("QREG" ,
                message );
    message.clear(); 

    const auto nAtomNos = box.getNAtomNos();

    for(const auto& atomNo:nAtomNos)
      if((box.getMdSymbol(atomNo) == "H1") ||
         (box.getMdSymbol(atomNo) == "H2") ||
         (box.getMdSymbol(atomNo) == "HB") ||
         (box.getMdSymbol(atomNo) == "HM") ||
         (box.getMdSymbol(atomNo) == "O1") ||
         (box.getMdSymbol(atomNo) == "O2")   )
        zundelIon.push_back(atomNo);

    message.append("Atoms ");

    for(const auto& atomNo:zundelIon)
    {

      box.insertMarkerNo(atomNo  ,
                         markerNo );
      message.append(std::to_string(atomNo));
      message.append(" "                   );

    }

    message.append("identified as Zundel ion and marked as center ");
    message.append(std::to_string(markerNo)                        );
    message.append("."                                             );
    log.message("QREG" ,
                message );
    message.clear();

    message.append("Marking Zundel ion by HIT types DONE");
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::modifyGradient(Box&                            box    ,
                       const std::vector<std::string>& command,
                       Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 6);


    //
    // inline auto
    // Qreg::modifyGradient(Box&                            box    ,
    //                      const std::vector<std::string>& command,
    //                      Logger&                         log     )const
    // ->void
    //

    std::string message;

    if(box.getXGradient(std::stoi(command[2])) != 0.0 &&
       box.getYGradient(std::stoi(command[2])) != 0.0 &&
       box.getZGradient(std::stoi(command[2])) != 0.0   )
      box.setGradient(std::stoi(command[2]),
                      std::stod(command[3])
       + box.getXGradient(std::stoi(command[2])),
                      std::stod(command[4])
       + box.getYGradient(std::stoi(command[2])),
                      std::stod(command[5])
       + box.getZGradient(std::stoi(command[2])) );
    
    box.setGradientGuess(std::stoi(command[2]),
                         std::stod(command[3])
     + box.getXGradientGuess(std::stoi(command[2])),
                         std::stod(command[4])
     + box.getYGradientGuess(std::stoi(command[2])),
                         std::stod(command[5])
     + box.getZGradientGuess(std::stoi(command[2])) );
    
    message.append("Modifying gradient of atom number ");
    message.append(command[2]                          );
    message.append(" by "                              );
    message.append(command[3]                          );
    message.append(" "                                 );
    message.append(command[4]                          );
    message.append(" "                                 );
    message.append(command[5]                          );
    message.append(" E -23 Joule per Angstrom."        );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::newQuantumRegion(std::vector<unsigned int>& nQuantumRegions,
                         unsigned int&              quantumRegionNo,
                         Logger&                    log             )const
  ->void
  {


   //
   // inline auto
   // Qreg::newQuantumRegion(std::vector<unsigned int>& nQuantumRegions,
   //                        unsigned int&              quantumRegionNo,
   //                        Logger&                    log             )const
   // ->void
   //

   std::string message;

   ++quantumRegionNo;

   nQuantumRegions.push_back(quantumRegionNo);


   message.append("Creating new QM region "      );
   message.append(std::to_string(quantumRegionNo));
   message.append("."                            );
   log.message("QREG" ,
               message );
   message.clear();


  }



  inline auto
  Qreg::resetMarkers(Box&                       box       ,
                     std::vector<unsigned int>& nMarkerNos,
                     unsigned int&              markerNo  ,
                     Logger&                    log        )const
  ->void
  {


    //
    // inline auto
    // Qreg::resetMarkers(Box&                       box       ,
    //                   std::vector<unsigned int>& nMarkerNos,
    //                   unsigned int&              markerNo  ,
    //                   Logger&                    log        )const
    // ->void
    //

    std::string message;

    markerNo = 0;

    nMarkerNos.clear();

    box.deleteAllMarkerNos();

    message.append("Resetting all markers.");
    log.message("QREG" ,
                message );
    message.clear();


   //
   // DbC POST
   //
   assert(nMarkerNos.empty());  


  }



  inline auto
  Qreg::seedAuto(const Box&                               box            ,
                 std::map< unsigned int          ,
                           std::unique_ptr<IfEvb> >&      nIfEvbs        ,
                 std::map< unsigned int          ,
                           std::unique_ptr<IfUff> >&      nIfUffs        ,
                 std::map< unsigned int               ,
                           std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs   ,
                 std::map< unsigned int              ,
                           std::unique_ptr<IfOpenUff> >& nIfOpenUffs     ,
                 std::map< unsigned int          ,
                           std::unique_ptr<IfOrc> >&      nIfOrcs        ,
                 const unsigned int                       quantumRegionNo,
                 const std::vector<std::string>&          command        ,
                 Logger&                                  log             )const
  ->void
  {


    switch(autoToolChoice_)
    {


      case QREG_EVB:
        seedEvb(box            ,
                nIfEvbs        ,
                quantumRegionNo,
                command        ,
                log             );
        break;

      case QREG_TFF:
        seedFf(box            ,
               nIfUffs        ,
               nIfOpenGaffs   ,
               nIfOpenUffs    ,
               quantumRegionNo,
               command        ,
               log             );
        break;

      case QREG_GFF:
        seedFf(box            ,
               nIfUffs        ,
               nIfOpenGaffs   ,
               nIfOpenUffs    ,
               quantumRegionNo,
               command        ,
               log             );
        break;

      case QREG_UFF:
        seedFf(box            ,
               nIfUffs        ,
               nIfOpenGaffs   ,
               nIfOpenUffs    ,
               quantumRegionNo,
               command        ,
               log             );
        break;

      case QREG_ORC:
        seedQm(box            ,
               nIfOrcs        ,
               quantumRegionNo,
               command        ,
               log             );
        break;

      default:
        seedQm(box            ,
               nIfOrcs        ,
               quantumRegionNo,
               command        ,
               log             );
        break;

    }


  }



  inline auto
  Qreg::seedEvb(const Box&                          box            ,
                std::map< unsigned int,
                          std::unique_ptr<IfEvb> >& nIfEvbs        ,
                const unsigned int                  quantumRegionNo,
                const std::vector<std::string>&     command        ,
                Logger&                             log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 5);


    //
    // inline auto
    // Qreg::seedEvb(const Box&                          box            ,
    //               std::map< unsigned int,
    //                         std::unique_ptr<IfEvb> >& nIfEvbs        ,
    //               const unsigned int                  quantumRegionNo,
    //               const std::vector<std::string>&     command        ,
    //               Logger&                             log             )const
    // ->void
    //

    std::string message;

    if(evbBackend_ == QREG_EVB)
    {

      nIfEvbs.emplace(quantumRegionNo                                 ,
                      std::unique_ptr<IfEvb>
                 (new IfEvb(configuration_.getEvbBinaryFilePath  (),
                            configuration_.getEvbHeader          (),
                            configuration_.getEvbOutputFolderPath() )) );

      if(command[2] != "auto")
        nIfEvbs.find(quantumRegionNo)->second->seed(box                  ,
                                                    quantumRegionNo      ,
                                                    std::stoi(command[2]),
                                                    std::stoi(command[3]),
                                                    std::stod(command[4]) );

      else
        nIfEvbs.find(quantumRegionNo)->second->seed(box                                    ,
                                                    quantumRegionNo                        ,
                                                    static_cast<unsigned int> (autoCharge_),
                                                    std::stoi(command[3])                  ,
                                                    std::stod(command[4])                   );

      message.append("An EVB calculation of QM region ");
      message.append(std::to_string(quantumRegionNo)   );
      message.append(" has been submitted."            );
      log.message("IFEVB",
                  message );
      message.clear();

    }


  }



  inline auto
  Qreg::seedFf(const Box&                               box            ,
               std::map< unsigned int          ,
                         std::unique_ptr<IfUff> >&      nIfUffs        ,
               std::map< unsigned int               ,
                         std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs   ,
               std::map< unsigned int              ,
                         std::unique_ptr<IfOpenUff> >&  nIfOpenUffs    ,
               const unsigned int                       quantumRegionNo,
               const std::vector<std::string>&          command        ,
               Logger&                                  log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 5);


    //
    // inline auto
    // Qreg::seedFf(const Box&                               box            ,
    //              std::map< unsigned int          ,
    //                        std::unique_ptr<IfUff> >&      nIfUffs        ,
    //              std::map< unsigned int               ,
    //                        std::unique_ptr<IfOpenGaff> >& nIfOpenGaffs   ,
    //              std::map< unsigned int              ,
    //                        std::unique_ptr<IfOpenUff> >&  nIfOpenUffs    ,
    //              const unsigned int                       quantumRegionNo,
    //              const std::vector<std::string>&          command        ,
    //              Logger&                                  log             )const
    // ->void
    //

    std::string message;

    if(ffBackend_ == QREG_TFF)
    {

      nIfUffs.emplace(quantumRegionNo                                 ,
                      std::unique_ptr<IfUff>
                 (new IfUff(configuration_.getFfBinaryFilePath  (),
                            configuration_.getFfHeader          (),
                            configuration_.getFfOutputFolderPath() )) );

      if(command[2] != "auto")
        nIfUffs.find(quantumRegionNo)->second->seed(box                  ,
                                                    quantumRegionNo      ,
                                                    std::stoi(command[2]),
                                                    std::stoi(command[3]),
                                                    std::stod(command[4]) );

      else
        nIfUffs.find(quantumRegionNo)->second->seed(box                                    ,
                                                    quantumRegionNo                        ,
                                                    static_cast<unsigned int> (autoCharge_),
                                                    std::stoi(command[3])                  ,
                                                    std::stod(command[4])                   );

      message.append("An UFF calculation of QM region ");
      message.append(std::to_string(quantumRegionNo)   );
      message.append(" has been submitted."            );
      log.message("IFUFF",
                  message );
      message.clear();

    }

    if(ffBackend_ == QREG_GFF)
    {

      nIfOpenGaffs.emplace(quantumRegionNo                                 ,
                           std::unique_ptr<IfOpenGaff>
                 (new IfOpenGaff(configuration_.getFfBinaryFilePath  (),
                                 configuration_.getFfHeader          (),
                                 configuration_.getFfOutputFolderPath() )) );

      if(command[2] != "auto")
        nIfOpenGaffs.find(quantumRegionNo)->second->seed(box                  ,
                                                         quantumRegionNo      ,
                                                         std::stoi(command[2]),
                                                         std::stoi(command[3]),
                                                         std::stod(command[4]) );

      else
        nIfOpenGaffs.find(quantumRegionNo)->second->seed(box                                    ,
                                                         quantumRegionNo                        ,
                                                         static_cast<unsigned int> (autoCharge_),
                                                         std::stoi(command[3])                  ,
                                                         std::stod(command[4])                   );

      message.append("An OpenGAFF calculation of QM region ");
      message.append(std::to_string(quantumRegionNo)        );
      message.append(" has been submitted."                 );
      log.message("IFOPENGAFF",
                  message      );
      message.clear();

    }

    if(ffBackend_ == QREG_UFF)
    {

      nIfOpenUffs.emplace(quantumRegionNo                                ,
                          std::unique_ptr<IfOpenUff>
                 (new IfOpenUff(configuration_.getFfBinaryFilePath  (),
                                configuration_.getFfHeader          (),
                                configuration_.getFfOutputFolderPath() )) );

      if(command[2] != "auto")
        nIfOpenUffs.find(quantumRegionNo)->second->seed(box                  ,
                                                        quantumRegionNo      ,
                                                        std::stoi(command[2]),
                                                        std::stoi(command[3]),
                                                        std::stod(command[4]) );

      else
        nIfOpenUffs.find(quantumRegionNo)->second->seed(box                                    ,
                                                        quantumRegionNo                        ,
                                                        static_cast<unsigned int> (autoCharge_),
                                                        std::stoi(command[3])                  ,
                                                        std::stod(command[4])                   );

      message.append("An OpenUFF calculation of QM region ");
      message.append(std::to_string(quantumRegionNo)       );
      message.append(" has been submitted."                );
      log.message("IFOPENUFF",
                  message      );
      message.clear();

    }


  }



  inline auto
  Qreg::seedQm(const Box&                          box            ,
               std::map< unsigned int,
                         std::unique_ptr<IfOrc> >& nIfOrcs        ,
               const unsigned int                  quantumRegionNo,
               const std::vector<std::string>&     command        ,
               Logger&                             log             )const
  ->void
  {


    //
    // DbC PRE
    //

    assert(quantumRegionNo != 0);
    assert(command.size()  == 5);


    //
    // inline auto
    // Qreg::seedQm(const Box&                          box            ,
    //              std::map< unsigned int,
    //                        std::unique_ptr<IfOrc> >& nIfOrcs        ,
    //              const unsigned int                  quantumRegionNo,
    //              const std::vector<std::string>&     command        ,
    //              Logger&                             log             )const
    // ->void
    //

    std::string message;

    if(qmBackend_ == QREG_ORC)
    {

      nIfOrcs.emplace(quantumRegionNo                                ,
                      std::unique_ptr<IfOrc>
                 (new IfOrc(configuration_.getQmBinaryFilePath  (),
                            configuration_.getQmHeader          (),
                            configuration_.getQmOutputFolderPath() )) );

      if(command[2] != "auto")
        nIfOrcs.find(quantumRegionNo)->second->seed(box                  ,
                                                    quantumRegionNo      ,
                                                    std::stoi(command[2]),
                                                    std::stoi(command[3]),
                                                    std::stod(command[4]) );

      else
        nIfOrcs.find(quantumRegionNo)->second->seed(box                                    ,
                                                    quantumRegionNo                        ,
                                                    static_cast<unsigned int> (autoCharge_),
                                                    std::stoi(command[3])                  ,
                                                    std::stod(command[4])                   );

      message.append("An ORC calculation of QM region ");
      message.append(std::to_string(quantumRegionNo)   );
      message.append(" has been submitted."            );
      log.message("IFORC",
                  message );
      message.clear();

    }


  }



  inline auto
  Qreg::setBufferInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setBufferInactive(Logger& log)
    // ->void
    //

    std::string message;

    bufferType_ = 0;
    message.append("Buffering scheme disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setBufferInnerRadius(const std::vector<std::string>& command,
                             Logger&                         log     )
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 4);


    //
    // inline auto
    // Qreg::setBufferInnerRadius(const std::vector<std::string>& command,
    //                            Logger&                         log     )
    // ->void
    //

    std::string message;

    innerBufferRadius_ = std::stod(command[3]);
    message.append("Using a buffer region with inner radius ");
    message.append(std::to_string(innerBufferRadius_)        );
    message.append(" Angstrom."                              );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setBufferOuterRadius(const std::vector<std::string>& command,
                             Logger&                         log     )
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 4);


    //
    // inline auto
    // Qreg::setBufferOuterRadius(const std::vector<std::string>& command,
    //                            Logger&                         log     )
    // ->void
    //

    std::string message;

    outerBufferRadius_ = std::stod(command[3]);
    message.append("Using a buffer region with outer radius ");
    message.append(std::to_string(outerBufferRadius_)        );
    message.append(" Angstrom."                              );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setBufferUsingLinearFit(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setBufferUsingLinearFit(Logger& log)
    // ->void
    //

    std::string message;

    bufferType_ = QREG_LINEARFIT;
    message.append("Buffering scheme enabled using linear fit.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setBufferUsingSmootherStep(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setBufferUsingSmootherStep(Logger& log)
    // ->void
    //

    std::string message;

    bufferType_ = QREG_SMOOTHERSTEP;
    message.append("Buffering scheme enabled using Smootherstep.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setBufferUsingSmoothStep(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setBufferUsingSmoothStep(Logger& log)
    // ->void
    //

    std::string message;

    bufferType_ = QREG_SMOOTHSTEP;
    message.append("Buffering scheme enabled using Smoothstep.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setEvbInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setEvbInactive(Logger& log)
    // ->void
    //

    std::string message;

    evbBackend_ = 0;
    message.append("EVB backend disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setEvbUsingEvb(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setEvbUsingEvb(Logger& log)
    // ->void
    //

    std::string message;

    evbBackend_ = QREG_EVB;
    message.append("Using EVB backend EVB at "          );
    message.append(configuration_.getEvbBinaryFilePath());
    message.append("."                                  );       
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setFfInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setFfInactive(Logger& log)
    // ->void
    //

    std::string message;

    ffBackend_ = 0;
    message.append("Force field backend disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setFfUsingUff(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setFfUsingUff(Logger& log)
    // ->void
    //

    std::string message;

    ffBackend_ = QREG_TFF;
    message.append("Using force field backend UFF at ");
    message.append(configuration_.getFfBinaryFilePath());
    message.append("."                                 );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setFfUsingOpenGaff(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setFfUsingOpenGaff(Logger& log)
    // ->void
    //

    std::string message;

    ffBackend_ = QREG_GFF;
    message.append("Using force field backend OpenGAFF at ");
    message.append(configuration_.getFfBinaryFilePath()    );
    message.append("."                                     );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setFfUsingOpenUff(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setFfUsingOpenUff(Logger& log)
    // ->void
    //

    std::string message;

    ffBackend_ = QREG_UFF;
    message.append("Using force field backend OpenUFF at ");
    message.append(configuration_.getFfBinaryFilePath()   );
    message.append("."                                    );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setGeometryUsingConfigFile(Box&    box,
                                   Logger& log )
  ->void
  {


    //
    // inline auto
    // Qreg::setGeometryUsingConfigFile(Box&    box,
    //                                  Logger& log )
    // ->void
    //

    if(runNo_ == 1)
    {

      std::string message;

      internalGeometry_ = true;
      message.append("Loading geometry from config file ...");
      log.message("SETUP",
                  message );
      message.clear();

      box = configuration_.getGeometry();

      message.append("Loading geometry from config file DONE");
      log.message("SETUP",
                  message );
      message.clear();

    }


  }



  inline auto
  Qreg::setGeometryUsingMdFrontend(Logger& log)
  ->void
  {


   //
   // inline auto
   // Qreg::setGeometryUsingMdFrontend(Logger& log)
   // ->void
   //

   std::string message;

   internalGeometry_ = false;
   message.append("Using geometry provided by MD frontend.");
   log.message("SETUP",
               message );
   message.clear();


  }



  inline auto
  Qreg::setGradient(Box&                            box    ,
                    const std::vector<std::string>& command,
                    Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 6);


    //
    // inline auto
    // Qreg::setGradient(Box&                            box    ,
    //                   const std::vector<std::string>& command,
    //                   Logger&                         log     )const
    // ->void
    //

    std::string message;

    box.setGradient(std::stoi(command[2]),
                    std::stod(command[3]),
                    std::stod(command[4]),
                    std::stod(command[5]) );

    box.setGradientGuess(std::stoi(command[2]),
                         std::stod(command[3]),
                         std::stod(command[4]),
                         std::stod(command[5]) );

    message.append("Setting gradient of atom number ");
    message.append(command[2]                        );
    message.append(" to "                            );
    message.append(command[3]                        );
    message.append(" "                               );
    message.append(command[4]                        );
    message.append(" "                               );
    message.append(command[5]                        );
    message.append(" E -23 Joule per Angstrom."      );
    log.message("QREG" ,
                message );
    message.clear();


  }



  inline auto
  Qreg::setHeatup(Box&                            box    ,
                  const std::vector<std::string>& command,
                  Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 3);


    //
    // inline auto
    // setHeatup(Box&                            box    ,
    //           const std::vector<std::string>& command,
    //           Logger&                         log     )const
    // ->void
    //

    std::string message;

    if(box.getTemperature() > std::stod(command[2]) *
                              runNo_                 )
    {

      box.setThermostat(std::stod(command[2]) *
                        runNo_                 ,
                        1.0                     ); 
 
      message.append("Heating up system by ");
      message.append(command[2]);
      message.append(" Kelvin per MD step. Current target: ");
      message.append(std::to_string(std::stod(command[2]) *
                                    runNo_                 ));
      message.append(" Kelvin.");
      log.message("THERMOSTAT",
                  message );

    }

    message.clear();


  }



  inline auto
  Qreg::setMdInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setMdInactive(Logger& log)
    // ->void
    //

    std::string message;

    mdFrontend_ = 0;
    message.append("MD frontend disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setMdUsingHit(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setMdUsingHit(Logger& log)
    // ->void
    //

    std::string message;

    mdFrontend_ = QREG_HIT;
    message.append("Using MD frontend HIT at "         );
    message.append(configuration_.getMdBinaryFilePath());
    message.append("."                                 );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setNeighborListsInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setNeighborListsInactive(Logger& log)
    // ->void
    //

    std::string message;

    stepCount_ = 1;
    stepSize_  = 0;

    message.append("Neighborlist scheme for distance calculation disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setNeighborListsStepCount(const std::vector<std::string>& command,
                                  Logger&                         log     )
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 4);


    //
    // inline auto
    // Qreg::setNeighborListsStepCount(const std::vector<std::string>& command,
    //                                 Logger&                         log     )
    // ->void
    //

    std::string message;

    stepCount_ = std::stoi(command[3]);
    message.append("Reevaluating all distances every ");
    message.append(std::to_string(stepCount_)         );
    message.append(" steps."                          );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setNeighborListsStepSize(const std::vector<std::string>& command,
                                 Logger&                         log     )
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 4);


    //
    // inline auto
    // Qreg::setNeighborListsStepSize(const std::vector<std::string>& command,
    //                                Logger&                         log     )
    // ->void
    //

    std::string message;

    stepSize_ = std::stoi(command[3]);
    message.append("Reevaluating distances invloving neighbors up to ");
    message.append(std::to_string(stepSize_)                          );
    message.append(" Angstrom apart."                                 );
    log.message("SETUP",
                message );
    message.clear();


  }



  inline auto
  Qreg::setQmInactive(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setQmInactive(Logger& log)
    // ->void
    //

    std::string message;

    qmBackend_ = 0;
    message.append("QM backend disabled.");
    log.message("SETUP",
                message );
    message.clear();


  }




  inline auto
  Qreg::setQmUsingOrc(Logger& log)
  ->void
  {


    //
    // inline auto
    // Qreg::setQmUsingOrc(Logger& log)
    // ->void
    //

    std::string message;

    qmBackend_ = QREG_ORC;
    message.append("Using QM backend ORC at ");
    message.append(configuration_.getQmBinaryFilePath());
    message.append("."                                 );
    log.message("SETUP",
                message );
    message.clear();


  }


  inline auto
  Qreg::setThermostat(Box&                            box    ,
                      const std::vector<std::string>& command,
                      Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 4);


    //
    // inline auto
    // Qreg::setThermostat(Box&                            box    ,
    //                     const std::vector<std::string>& command,
    //                     Logger&                         log     )const
    // ->void
    //

    std::string message;

    box.setThermostat(std::stod(command[2]),
                      std::stod(command[3]) );

    message.append("Setting temperature to "                                  ); 
    message.append(command[2]                                                 );
    message.append(" Kelvin and the Berendsen Thermostat dampening factor to ");
    message.append(command[3]                                                 );
    message.append(" Femtoseconds."                                           );
    log.message("THERMOSTAT",
                message      );
    message.clear();


  }



  inline auto
  Qreg::setTimeFrame(Box&                            box    ,
                     const std::vector<std::string>& command,
                     Logger&                         log     )const
  ->void
  {


    //
    // DbC PRE
    //
    assert(command.size() == 3);


    //
    // inline auto
    // Qreg::setTimeframe(Box&                            box    ,
    //                    const std::vector<std::string>& command,
    //                    Logger&                         log     )const
    // ->void
    //

    std::string message;

    box.setTimeFrame(std::stod(command[2]));

    message.append("Setting MD time frame to ");
    message.append(command[2]                 );
    message.append(" Femtoseconds."           );
    log.message("MD_DRIVER",
                message     );
    message.clear();


  }



  inline auto
  Qreg::updateTrajectory(Box&    box,
                         Logger& log )const
  ->void
  {


    //
    // inline auto
    // Qreg::updateTrajectory(Box&    box,
    //                        Logger& log )const
    // ->void
    //

    std::string message;

    auto bulkTrajectoryFilePath = configuration_.getTmpFolderPath();
    bulkTrajectoryFilePath.append("/bulk.xyz");

    message.append("Updating trajectory files ...");
    log.message("QREG" ,
                message );
    message.clear();

    std::string trajectoryText;
    trajectoryText.append(std::to_string(box.getAtomCount()));
    trajectoryText.append("\n"                              );
    trajectoryText.append("Qreg Trajectory Step "           );
    trajectoryText.append(std::to_string(runNo_)            );
    trajectoryText.append(" "                               );
    trajectoryText.append(std::to_string(box.getTimeFrame()));
    trajectoryText.append(" Femtoseconds"                   );
    trajectoryText.append("\n"                              );

    for(const auto& atomNo:box.getNAtomNos())
    {

      trajectoryText.append(box.getChemicalSymbol(atomNo)             );
      trajectoryText.append("   "                                     );
      trajectoryText.append(std::to_string(box.getXCoordinate(atomNo)));
      trajectoryText.append("   "                                     );
      trajectoryText.append(std::to_string(box.getYCoordinate(atomNo)));
      trajectoryText.append("   "                                     );
      trajectoryText.append(std::to_string(box.getZCoordinate(atomNo)));
      trajectoryText.append("\n"                                      );

    }

    std::ofstream bulkTrajectoryFile;
    bulkTrajectoryFile.open(bulkTrajectoryFilePath.c_str(),
                            std::fstream::app              );

    if(bulkTrajectoryFile.good())
      bulkTrajectoryFile << trajectoryText;

    trajectoryText.clear();
    bulkTrajectoryFile.close();


    unsigned int quantumRegionNo = 1;
    while(box.hasQuantumRegionNo(quantumRegionNo))
    {

      auto quantumRegionTrajectoryFilePath = configuration_.getTmpFolderPath();
      quantumRegionTrajectoryFilePath.append("/qreg_"                       );
      quantumRegionTrajectoryFilePath.append(std::to_string(quantumRegionNo));
      quantumRegionTrajectoryFilePath.append(".xyz"                         );

      trajectoryText.append(std::to_string(
       box.getNAtomNosByQuantumRegionNo(quantumRegionNo).size()));
      trajectoryText.append("\n"                              );
      trajectoryText.append("Qreg Trajectory Step "           );
      trajectoryText.append(std::to_string(runNo_)            );
      trajectoryText.append(" "                               );
      trajectoryText.append(std::to_string(box.getTimeFrame()));
      trajectoryText.append(" Femtoseconds"                   );
      trajectoryText.append("\n"                              );

      for(const auto& atomNo:box.getNAtomNosByQuantumRegionNo(quantumRegionNo))
      {

        trajectoryText.append(box.getChemicalSymbol(atomNo)            );
        trajectoryText.append("   "                                    );
        trajectoryText.append(std::to_string(box.getXDistance(1     ,
                                                              atomNo )));
        trajectoryText.append("   "                                    );
        trajectoryText.append(std::to_string(box.getYDistance(1     ,
                                                              atomNo )));
        trajectoryText.append("   "                                    );
        trajectoryText.append(std::to_string(box.getZDistance(1     ,
                                                              atomNo )));
        trajectoryText.append("\n"                                     );

      }

      std::ofstream quantumRegionTrajectoryFile;
      quantumRegionTrajectoryFile.open(quantumRegionTrajectoryFilePath.c_str(),
                                       std::fstream::app                       );

      if(quantumRegionTrajectoryFile.good())
        quantumRegionTrajectoryFile << trajectoryText;

      trajectoryText.clear();
      quantumRegionTrajectoryFile.close();

      ++quantumRegionNo;

    }


    message.append("Updating trajectory files DONE");
    log.message("QREG" ,
                message );
    message.clear();


  }



}



#endif // QREG_QREG_HPP_
