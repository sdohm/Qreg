//
// qreg.cpp
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



#include "qreg.hpp"


namespace qreg
{



  Qreg::Qreg(const std::string&  configFilePath   ,
             const unsigned int& newAutoToolChoice,
             const double&       newAutoCharge     ):

    internalGeometry_ (true             ),
    autoToolChoice_   (newAutoToolChoice),
    mdFrontend_       (QREG_HIT         ),
    ffBackend_        (QREG_GFF         ),
    qmBackend_        (QREG_ORC         ),
    evbBackend_       (QREG_EVB         ),
    runNo_            ( 1               ),
    stepCount_        ( 3               ),
    bufferType_       (QREG_SMOOTHSTEP  ),
    autoCharge_       (newAutoCharge    ),
    innerBufferRadius_( 8.0             ),
    outerBufferRadius_(16.0             ),
    stepSize_         (16.0             ),
    configuration_    (                 )

  {



    configuration_.initialize(configFilePath);

    Logger log(configuration_.getLogFilePath());

    std::string message;

    message.append("Qreg 1.2.1 Copyright 2014-2016 Sebastian Dohm <sebastian.dohm@uni-ulm.de>\n\n");
    message.append("This program is free software: you can redistribute it and/or modify\n"       );
    message.append("it under the terms of the GNU General Public License as published by\n"       );
    message.append("the Free Software Foundation, either version 3 of the License, or\n"          );
    message.append("(at your option) any later version.\n\n"                                      );
    message.append("This program is distributed in the hope that it will be useful,\n"            );
    message.append("but WITHOUT ANY WARRANTY; without even the implied warranty of\n"             );
    message.append("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"              );
    message.append("GNU General Public License for more details.\n\n"                             );
    message.append("You should have received a copy of the GNU General Public License\n"          );
    message.append("along with this program.  If not, see <http://www.gnu.org/licenses/>."        );
    log.message("COPYRIGHT",
                message  );
    message.clear();

    #ifndef NDEBUG

      message.append("KNOWN BUGS\n\n"                                     );
      message.append("- forcefield lacks sufficient performance\n"        );
      message.append("\nRECENT CHANGES\n\n"                               );
      message.append("- several problems with periodicity solved\n"       );
      message.append("- support for POSIX interprocess management added\n");
      message.append("- Berendsen Thermostat added\n"                     );
      message.append("- restart files available\n"                        );
      message.append("- support for EVB added\n"                          );
      message.append("\nTHINGS TO COME\n\n"                               );
      message.append("- handling graphite\n"                              );
      message.append("- high performance forcefield"                      );
      log.message("README",
                  message );
      message.clear();
    
    #endif

    message.append("Configuration options have been loaded from ");
    message.append(configFilePath                                );
    message.append("."                                           );
    log.message("CONFIG",
                message  );
    message.clear();

    message.append("Logging to "                  );
    message.append(configuration_.getLogFilePath());
    message.append("."                            );
    log.message("LOGGER",
                message  );
    message.clear();

    #ifndef NDEBUG

    #ifdef DEBUG_CONFIG

      message.append("TESTING CONFIGURATION\n\n");

      for(const auto& line:configuration_.getNCommands())
      {
        message.append("CMD ");
        for(const auto& word:line)
        {
          message.append(word);
          message.append(" " );
        }
        message.append("\n");
      }

      message.append("configuration_.getLogFilePath returns: "         );
      message.append(configuration_.getLogFilePath()                   );
      message.append("\n"                                              );
      message.append("configuration_.getMdBinaryFilePath returns: "    );
      message.append(configuration_.getMdBinaryFilePath()              );
      message.append("\n"                                              );
      message.append("configuration_.getMdOutputFolderPath returns: "  );
      message.append(configuration_.getMdOutputFolderPath()            );
      message.append("\n"                                              );
      message.append("configuration_.getFfBinaryFilePath returns: "    );
      message.append(configuration_.getFfBinaryFilePath()              );
      message.append("\n"                                              );
      message.append("configuration_.getFfOutputFolderPath returns: "  );
      message.append(configuration_.getFfOutputFolderPath()            );
      message.append("\n"                                              );
      message.append("configuration_.getQmBinaryFilePath returns: "    );
      message.append(configuration_.getQmBinaryFilePath()              );
      message.append("\n"                                              );
      message.append("configuration_.getQmOutputFolderPath returns: "  );
      message.append(configuration_.getQmOutputFolderPath()            );
      message.append("\n"                                              );
      message.append("configuration_.getEvbBinaryFolderPath returns: " );
      message.append(configuration_.getEvbBinaryFilePath()             );
      message.append("\n"                                              );
      message.append("configuration_.getEvbOutputFolderPath returns: " );
      message.append(configuration_.getEvbOutputFolderPath()           );
      message.append("\n"                                              );
      message.append("configuration_.getTmpFolderPath returns: "       );
      message.append(configuration_.getTmpFolderPath()                 );
      message.append("\n"                                              );
      message.append("configuration_.getMdHeader returns:\n"           );
      message.append(configuration_.getMdHeader()                      );
      message.append("\n"                                              );
      message.append("configuration_.getFfHeader returns:\n"           );
      message.append(configuration_.getFfHeader()                      );
      message.append("\n"                                              );
      message.append("configuration_.getQmHeader returns:\n"           );
      message.append(configuration_.getQmHeader()                      );
      message.append("\n"                                              );
      message.append("configuration_getEvbHeader returns:\n"           );
      message.append(configuration_.getEvbHeader()                     );

      for(const auto& atomNo:configuration_.getGeometry().getNAtomNos())
      {

        message.append("\n"                                                  );
        message.append("atomNo: "                                            );
        message.append(std::to_string(atomNo)                                );
        message.append(" symbol: "                                           );
        message.append(configuration_.getGeometry().getChemicalSymbol(atomNo));
        message.append(" X: "                                                );
        message.append
         (std::to_string(configuration_.getGeometry().getXCoordinate(atomNo)));
        message.append(" Y: "                                                );
        message.append
         (std::to_string(configuration_.getGeometry().getYCoordinate(atomNo)));
        message.append(" Z: "                                                );
        message.append
         (std::to_string(configuration_.getGeometry().getZCoordinate(atomNo)));

      }

      log.message("DEBUG",
                  message );
      message.clear();

    #endif

    #endif
   


  }



  Qreg::Qreg(Qreg&& qreg)
  {


    // Calling move assignment operator.
    *this = std::move(qreg);


  }



  auto
  Qreg::operator=(Qreg&& qreg)
  ->Qreg&
  {


    //
    // Perform move assignment.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    std::swap(     internalGeometry_,
              qreg.internalGeometry_ );

    std::swap(     autoToolChoice_,
              qreg.autoToolChoice_ );

    std::swap(     mdFrontend_,
              qreg.mdFrontend_ );

    std::swap(     ffBackend_,
              qreg.ffBackend_ );

    std::swap(     qmBackend_,
              qreg.qmBackend_ );

    std::swap(     evbBackend_,
              qreg.evbBackend_ );

    std::swap(     runNo_,
              qreg.runNo_ );

    std::swap(     stepCount_,
              qreg.stepCount_ );

    std::swap(     bufferType_,
              qreg.bufferType_ );

    std::swap(     autoCharge_,
              qreg.autoCharge_ );

    std::swap(     innerBufferRadius_,
              qreg.innerBufferRadius_ );

    std::swap(     outerBufferRadius_,
              qreg.outerBufferRadius_ );

    std::swap(     stepSize_,
              qreg.stepSize_ );

    std::swap(     configuration_,
              qreg.configuration_ );

    qreg.~Qreg();


    return *this;


  }



  Qreg::~Qreg(void)
  {



  }



  auto
  Qreg::advance(Box& box)
  ->void
  {


    //
    // auto
    // Qreg::advance(Box& box)
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    unsigned int markerNo        = 0;
    unsigned int quantumRegionNo = 0;

    std::vector<unsigned int> nQuantumRegions;
    std::vector<unsigned int> nMarkerNos     ;

    std::string message;

    Logger log(configuration_.getLogFilePath());

    std::map <unsigned int          ,
              std::unique_ptr<IfEvb> > nIfEvbs;

    std::map <unsigned int          ,
              std::unique_ptr<IfOrc> > nIfOrcs;

    std::map <unsigned int          ,
              std::unique_ptr<IfUff> > nIfUffs;

    std::map <unsigned int               ,
              std::unique_ptr<IfOpenGaff> > nIfOpenGaffs;

    std::map <unsigned int              ,
              std::unique_ptr<IfOpenUff> > nIfOpenUffs;

    const auto nCommands = configuration_.getNCommands();

    message.append("Executing user commands from " );
    message.append(std::to_string(nCommands.size()));
    message.append(" lines of input. Run No. "     );
    message.append(std::to_string(runNo_)          );
    message.append(" ..."                          );
    log.message("QREG" ,
                message );
    message.clear();


    for(size_t commandNo  = 0               ;
               commandNo != nCommands.size();
             ++commandNo                     )
    {


      if(nCommands[commandNo].size() >= 1)
      {

        message.append("<< "                  );
        message.append(nCommands[commandNo][0]);
        message.append(" "                    );

        if(nCommands[commandNo].size() >= 2)
        {

          message.append(nCommands[commandNo][1]);
          message.append(" "                    );

        }

        if(nCommands[commandNo].size() >= 3)
        {

          message.append(nCommands[commandNo][2]);
          message.append(" "                    );

        }

        if(nCommands[commandNo].size() >= 4)
        {

          message.append(nCommands[commandNo][3]);
          message.append(" "                    );

        }

        if(nCommands[commandNo].size() >= 5)
        {

          message.append(nCommands[commandNo][4]);
          message.append(" "                    );

        }

        message.append(">>");
        log.message("COMMAND",
                    message   );
        message.clear();

      }


      if((nCommands[commandNo].size() == 5           ) &&
         (nCommands[commandNo][0]     == "add"       ) &&
         (nCommands[commandNo][1]     == "ATTRACTION")   )
        addAttraction(box                 ,
                      nCommands[commandNo],
                      log                  );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "advance") &&
              (nCommands[commandNo][1]     == "MD"     )   )
        advanceMd(box,
                  log );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "compute") &&
              (nCommands[commandNo][1]     == "BUFFER" )   )
        computeBuffer(box     ,
                      markerNo,
                      log      );

      else if((nCommands[commandNo].size() == 2          ) &&
              (nCommands[commandNo][0]     == "compute"  ) &&
              (nCommands[commandNo][1]     == "DISTANCES")   )
        computeDistances(box,
                         log );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "harvest") &&
              (nCommands[commandNo][1]     == "AUTO"   )   )
        harvestAuto(box         ,
                    nIfEvbs     ,
                    nIfUffs     ,
                    nIfOpenGaffs,
                    nIfOpenUffs ,
                    nIfOrcs     ,
                    log          );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "harvest") &&
              (nCommands[commandNo][1]     == "EVB"    )   )
        harvestEvb(box    ,
                   nIfEvbs,
                   log     );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "harvest") &&
              (nCommands[commandNo][1]     == "FF"     )   )
        harvestFf(box         ,
                  nIfUffs     ,
                  nIfOpenGaffs,
                  nIfOpenUffs ,
                  log          );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "harvest") &&
              (nCommands[commandNo][1]     == "QM"     )   )
        harvestQm(box    ,
                  nIfOrcs,
                  log     );

      else if((nCommands[commandNo].size() == 3        ) &&
              (nCommands[commandNo][0]     == "include") &&
              (nCommands[commandNo][1]     == "ATOM"   )   )
        includeAtom(box                 ,
                    quantumRegionNo     ,
                    nCommands[commandNo],
                    log                  );

      else if((nCommands[commandNo].size() == 3        ) &&
              (nCommands[commandNo][0]     == "include") &&
              (nCommands[commandNo][1]     == "ATOMS"  ) &&
              (nCommands[commandNo][2]     == "all"    )   )
        includeAtomsAll(box            ,
                        quantumRegionNo,
                        log             );

      else if((nCommands[commandNo].size() == 4         ) &&
              (nCommands[commandNo][0]     == "include" ) &&
              (nCommands[commandNo][1]     == "ATOMS"   ) &&
              (nCommands[commandNo][2]     == "by_range")   )
        includeAtomsByRange(box                 ,
                            markerNo            ,
                            quantumRegionNo     ,
                            nCommands[commandNo],
                            log                  );

      else if((nCommands[commandNo].size() == 4         ) &&
              (nCommands[commandNo][0]     == "include" ) &&
              (nCommands[commandNo][1]     == "MOLECULE")   )
        includeMolecule(box                 ,
                        quantumRegionNo     ,
                        nCommands[commandNo],
                        log                  );

      else if((nCommands[commandNo].size() == 5          ) &&
              (nCommands[commandNo][0]     == "include"  ) &&
              (nCommands[commandNo][1]     == "MOLECULES") &&
              (nCommands[commandNo][2]     == "by_range" )   )
        includeMoleculesByRange(box                 ,
                                markerNo            ,
                                quantumRegionNo     ,
                                nCommands[commandNo],
                                log                  );

      else if((nCommands[commandNo].size() == 3     ) &&
              (nCommands[commandNo][0]     == "mark") &&
              (nCommands[commandNo][1]     == "ATOM")   )
        markAtom(box                 ,
                 nMarkerNos          ,
                 markerNo            ,
                 nCommands[commandNo],
                 log                  );

      else if((nCommands[commandNo].size() == 2             ) &&
              (nCommands[commandNo][0]     == "mark"        ) &&
              (nCommands[commandNo][1]     == "ZUNDELION"   ) &&
              (nCommands[commandNo][2]     == "by_hit_types")   )
        markZundelIonByHitTypes(box       ,
                                nMarkerNos,
                                markerNo  ,
                                log        );

      else if((nCommands[commandNo].size() == 6         ) &&
              (nCommands[commandNo][0]     == "modify"  ) &&
              (nCommands[commandNo][1]     == "GRADIENT")   )
        modifyGradient(box                 ,
                       nCommands[commandNo],
                       log                  );

      else if((nCommands[commandNo].size() == 2              ) &&
              (nCommands[commandNo][0]     == "new"          ) &&
              (nCommands[commandNo][1]     == "QUANTUMREGION")   )
        newQuantumRegion(nQuantumRegions,
                         quantumRegionNo,
                         log             );

      else if((nCommands[commandNo].size() == 2        ) &&
              (nCommands[commandNo][0]     == "reset"  ) &&
              (nCommands[commandNo][1]     == "MARKERS")   )
        resetMarkers(box       ,
                     nMarkerNos,
                     markerNo  ,
                     log        );

      else if((nCommands[commandNo].size() == 6     ) &&
              (nCommands[commandNo][0]     == "seed") &&
              (nCommands[commandNo][1]     == "AUTO")   )
        seedAuto(box                 ,
                 nIfEvbs             ,
                 nIfUffs             ,
                 nIfOpenGaffs        ,
                 nIfOpenUffs         ,
                 nIfOrcs             ,
                 quantumRegionNo     ,
                 nCommands[commandNo], 
                 log                  );

      else if((nCommands[commandNo].size() == 6     ) &&
              (nCommands[commandNo][0]     == "seed") &&
              (nCommands[commandNo][1]     == "EVB" )   )
        seedEvb(box                 ,
                nIfEvbs             ,
                quantumRegionNo     ,
                nCommands[commandNo], 
                log                  );

      else if((nCommands[commandNo].size() == 6     ) &&
              (nCommands[commandNo][0]     == "seed") &&
              (nCommands[commandNo][1]     == "FF"  )   )
        seedFf(box                 ,
               nIfUffs             ,
               nIfOpenGaffs        ,
               nIfOpenUffs         ,
               quantumRegionNo     ,
               nCommands[commandNo], 
               log                  );

      else if((nCommands[commandNo].size() == 6     ) &&
              (nCommands[commandNo][0]     == "seed") &&
              (nCommands[commandNo][1]     == "QM"  )   )
        seedQm(box                 ,
               nIfOrcs             ,
               quantumRegionNo     ,
               nCommands[commandNo],
               log                  );

      else if((nCommands[commandNo].size() == 2         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "BUFFER"  ) &&
              (nCommands[commandNo][2]     == "inactive")   )
        setBufferInactive(log);

      else if((nCommands[commandNo].size() == 4             ) &&
              (nCommands[commandNo][0]     == "set"         ) &&
              (nCommands[commandNo][1]     == "BUFFER"      ) &&
              (nCommands[commandNo][2]     == "inner_radius")   )
        setBufferInnerRadius(nCommands[commandNo],
                             log                  );

      else if((nCommands[commandNo].size() == 4             ) &&
              (nCommands[commandNo][0]     == "set"         ) &&
              (nCommands[commandNo][1]     == "BUFFER"      ) &&
              (nCommands[commandNo][2]     == "outer_radius")   )
        setBufferOuterRadius(nCommands[commandNo],
                             log                  );

      else if((nCommands[commandNo].size() == 3                 ) &&
              (nCommands[commandNo][0]     == "set"             ) &&
              (nCommands[commandNo][1]     == "BUFFER"          ) &&
              (nCommands[commandNo][2]     == "using_linear_fit")   )
        setBufferUsingLinearFit(log);

      else if((nCommands[commandNo].size() == 3                   ) &&
              (nCommands[commandNo][0]     == "set"               ) &&
              (nCommands[commandNo][1]     == "BUFFER"            ) &&
              (nCommands[commandNo][2]     == "using_smootherstep")   )
        setBufferUsingSmootherStep(log);

      else if((nCommands[commandNo].size() == 3                 ) &&
              (nCommands[commandNo][0]     == "set"             ) &&
              (nCommands[commandNo][1]     == "BUFFER"          ) &&
              (nCommands[commandNo][2]     == "using_smoothstep")   )
        setBufferUsingSmoothStep(log);

      else if((nCommands[commandNo].size() == 3         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "EVB"     ) &&
              (nCommands[commandNo][2]     == "inactive")   )
        setEvbInactive(log);

      else if((nCommands[commandNo].size() == 3          ) &&
              (nCommands[commandNo][0]     == "set"      ) &&
              (nCommands[commandNo][1]     == "EVB"      ) &&
              (nCommands[commandNo][2]     == "using_evb")   )
        setEvbUsingEvb(log);

      else if((nCommands[commandNo].size() == 3         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "FF"      ) &&
              (nCommands[commandNo][2]     == "inactive")   )
        setFfInactive(log);

      else if((nCommands[commandNo].size() == 3          ) &&
              (nCommands[commandNo][0]     == "set"      ) &&
              (nCommands[commandNo][1]     == "FF"       ) &&
              (nCommands[commandNo][2]     == "using_uff")   )
        setFfUsingUff(log);

      else if((nCommands[commandNo].size() == 3               ) &&
              (nCommands[commandNo][0]     == "set"           ) &&
              (nCommands[commandNo][1]     == "FF"            ) &&
              (nCommands[commandNo][2]     == "using_opengaff")   )
        setFfUsingOpenGaff(log);

      else if((nCommands[commandNo].size() == 3               ) &&
              (nCommands[commandNo][0]     == "set"           ) &&
              (nCommands[commandNo][1]     == "FF"            ) &&
              (nCommands[commandNo][2]     == "using_openuff")   )
        setFfUsingOpenUff(log);

      else if((nCommands[commandNo].size() == 3                  ) &&
              (nCommands[commandNo][0]     == "set"              ) &&
              (nCommands[commandNo][1]     == "GEOMETRY"         ) &&
              (nCommands[commandNo][2]     == "using_config_file")   )
        setGeometryUsingConfigFile(box,
                                   log );

      else if((nCommands[commandNo].size() == 3                  ) &&
              (nCommands[commandNo][0]     == "set"              ) &&
              (nCommands[commandNo][1]     == "GEOMETRY"         ) &&
              (nCommands[commandNo][2]     == "using_md_frontend")   )
        setGeometryUsingMdFrontend(log);

      else if((nCommands[commandNo].size() == 6         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "GRADIENT")   )
        setGradient(box                 ,
                    nCommands[commandNo],
                    log                  );

      else if((nCommands[commandNo].size() == 3       ) &&
              (nCommands[commandNo][0]     == "set"   ) &&
              (nCommands[commandNo][1]     == "HEATUP")   )
        setHeatup(box                 ,
                  nCommands[commandNo],
                  log                  );

      else if((nCommands[commandNo].size() == 3         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "MD"      ) &&
              (nCommands[commandNo][2]     == "inactive")   )
        setMdInactive(log);

      else if((nCommands[commandNo].size() == 3          ) &&
              (nCommands[commandNo][0]     == "set"      ) &&
              (nCommands[commandNo][1]     == "MD"       ) &&
              (nCommands[commandNo][2]     == "using_hit")   )
        setMdUsingHit(log);

      else if((nCommands[commandNo].size() == 2              ) &&
              (nCommands[commandNo][0]     == "set"          ) &&
              (nCommands[commandNo][1]     == "NEIGHBORLISTS") && 
              (nCommands[commandNo][2]     == "inactive"     )   )
        setNeighborListsInactive(log);

      else if((nCommands[commandNo].size() == 4              ) &&
              (nCommands[commandNo][0]     == "set"          ) &&
              (nCommands[commandNo][1]     == "NEIGHBORLISTS") &&
              (nCommands[commandNo][2]     == "step_count"   )   )
        setNeighborListsStepCount(nCommands[commandNo],
                                  log                  );

      else if((nCommands[commandNo].size() == 4              ) &&
              (nCommands[commandNo][0]     == "set"          ) &&
              (nCommands[commandNo][1]     == "NEIGHBORLISTS") &&
              (nCommands[commandNo][2]     == "step_size"    )   )
        setNeighborListsStepSize(nCommands[commandNo],
                                 log                  );

      else if((nCommands[commandNo].size() == 3         ) &&
              (nCommands[commandNo][0]     == "set"     ) &&
              (nCommands[commandNo][1]     == "QM"      ) &&
              (nCommands[commandNo][2]     == "inactive")   )
        setQmInactive(log);

      else if((nCommands[commandNo].size() == 3          ) &&
              (nCommands[commandNo][0]     == "set"      ) &&
              (nCommands[commandNo][1]     == "QM"       ) &&
              (nCommands[commandNo][2]     == "using_orc")   )
        setQmUsingOrc(log);

      else if((nCommands[commandNo].size() == 4           ) &&
              (nCommands[commandNo][0]     == "set"       ) &&
              (nCommands[commandNo][1]     == "THERMOSTAT")   )
        setThermostat(box                 ,
                      nCommands[commandNo],
                      log                  );

      else if((nCommands[commandNo].size() == 3          ) &&
              (nCommands[commandNo][0]     == "set"      ) &&
              (nCommands[commandNo][1]     == "TIMEFRAME")   )
        setTimeFrame(box                 ,
                     nCommands[commandNo],
                     log                  );

      else if((nCommands[commandNo].size() == 2           ) &&
              (nCommands[commandNo][0]     == "update"    ) &&
              (nCommands[commandNo][1]     == "TRAJECTORY")   ) 
        updateTrajectory(box,
                         log );

      else if(nCommands[commandNo].size() != 0)
      {

        message.append("Unknown command.");
        log.message("COMMAND",
                    message   );
        message.clear();

      }



    }

    box.deleteAllEPots           ();
    box.deleteAllEvbEPots        ();
    box.deleteAllEvbGradients    ();
    box.deleteAllEvbTotalEpots   ();
    box.deleteAllMarkerNos       ();
    box.deleteAllQuantumRegionNos();
    box.equalizeGradients        ();

    message.append("Run No. "            );
    message.append(std::to_string(runNo_));
    message.append(" DONE"               );
    log.message("QREG" ,
                message );
    message.clear();


    ++runNo_;



  }



}
