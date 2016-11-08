//
// config.cpp
//
//  Created on: Nov 20, 2014
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



#include "config.hpp"



namespace qreg
{



  Config::Config(void):

    configFilePath_     ( ),
    logFilePath_        ( ),
    mdBinaryFilePath_   ( ),
    mdOutputFolderPath_ ( ),
    ffBinaryFilePath_   ( ),
    ffOutputFolderPath_ ( ),
    qmBinaryFilePath_   ( ),
    qmOutputFolderPath_ ( ),
    evbBinaryFilePath_  ( ),
    evbOutputFolderPath_( ),
    tmpFolderPath_      ( ),
    mdHeader_           ( ),
    ffHeader_           ( ),
    qmHeader_           ( ),
    evbHeader_          ( ),
    nCommands_          ( ),
    geometry_           ( )

  {



  }



  Config::Config(Config&& config)
  {


    // Calling move assignment operator.
    *this = std::move(config);


  }



  auto
  Config::operator=(Config&& config)
  ->Config&
  {


    //
    // Perform move assignment.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    std::swap(       configFilePath_,
              config.configFilePath_ );

    std::swap(       logFilePath_,
              config.logFilePath_ );

    std::swap(       mdBinaryFilePath_,
              config.mdBinaryFilePath_ );

    std::swap(       mdOutputFolderPath_,
              config.mdOutputFolderPath_ );

    std::swap(       ffBinaryFilePath_,
              config.ffBinaryFilePath_ );

    std::swap(       ffOutputFolderPath_,
              config.ffOutputFolderPath_ );

    std::swap(       qmBinaryFilePath_,
              config.qmBinaryFilePath_ );

    std::swap(       qmOutputFolderPath_,
              config.qmOutputFolderPath_ );

    std::swap(       evbBinaryFilePath_,
              config.evbBinaryFilePath_ );

    std::swap(       evbOutputFolderPath_,
              config.evbOutputFolderPath_ );

    std::swap(       tmpFolderPath_,
              config.tmpFolderPath_ );

    std::swap(       mdHeader_,
              config.mdHeader_ );

    std::swap(       ffHeader_,
              config.ffHeader_ );

    std::swap(       qmHeader_,
              config.qmHeader_ );

    std::swap(       evbHeader_,
              config.evbHeader_ );

    std::swap(       nCommands_,
              config.nCommands_ );

    std::swap(       geometry_,
              config.geometry_ );

    config.~Config();


    return *this;


  }



  Config::~Config(void)
  {



  }



  auto
  Config::initialize(const std::string& newConfigFilePath)
  ->void
  {


    //
    // DbC PRE
    //
    assert(!newConfigFilePath.empty());


    //
    // auto
    // Config::initialize(const std::string& newConfigFilePath)
    // ->void
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    std::ifstream configFile;

    std::string line           ;
    std::string restartFilePath;

    std::vector<std::string> nLines;

    std::vector<std::vector <std::string> > nCoordinates   ;
    std::vector<std::vector <std::string> > nNewCoordinates;

    off_t logLine            = -1;
    off_t mdHeaderBeginLine  = -1;
    off_t mdHeaderEndLine    = -1;
    off_t mdBinaryLine       = -1;
    off_t mdOutputLine       = -1;
    off_t ffHeaderBeginLine  = -1;
    off_t ffHeaderEndLine    = -1;
    off_t ffBinaryLine       = -1;
    off_t ffOutputLine       = -1;
    off_t qmHeaderBeginLine  = -1;
    off_t qmHeaderEndLine    = -1;
    off_t qmBinaryLine       = -1;
    off_t qmOutputLine       = -1;
    off_t evbHeaderBeginLine = -1;
    off_t evbHeaderEndLine   = -1;
    off_t evbBinaryLine      = -1;
    off_t evbOutputLine      = -1;
    off_t tmpLine            = -1;
    off_t geoBeginLine       = -1;
    off_t geoEndLine         = -1;
    off_t boxXLine           = -1;
    off_t boxYLine           = -1;
    off_t boxZLine           = -1;
    off_t cmdBeginLine       = -1;
    off_t cmdEndLine         = -1;

    configFilePath_ = newConfigFilePath;

    configFile.open(configFilePath_.c_str());

    if(configFile.bad())
    {

      std::cerr << "Qreg configuration file not available at "
                << configFilePath_
                << ", aborting."
                << std::endl                                  ;

      exit(-1);

    }

    if(configFile.is_open())
    {

      while(std::getline(configFile,
                         line       ))
        nLines.push_back(line);

      configFile.close();

    }

    assert(!nLines.empty());

    for(off_t i  = 0                                 ;
              i != static_cast<signed>(nLines.size());
            ++i                                       )
    {

      if(i != static_cast<signed>(nLines.size() - 1))
      {

        if(nLines[i] == "[LOG]")
          logLine = i;

        if(nLines[i] == "[MD_HEADER]")
          mdHeaderBeginLine = i;

        if(nLines[i] == "[MD_BIN]")
          mdBinaryLine = i;

        if(nLines[i] == "[MD_OUTPUT]")
          mdOutputLine = i;

        if(nLines[i] == "[FF_HEADER]")
          ffHeaderBeginLine = i;

        if(nLines[i] == "[FF_BIN]")
          ffBinaryLine = i;

        if(nLines[i] == "[FF_OUTPUT]")
          ffOutputLine = i;

        if(nLines[i] == "[QM_HEADER]")
          qmHeaderBeginLine = i;

        if(nLines[i] == "[QM_BIN]")
          qmBinaryLine = i;

        if(nLines[i] == "[QM_OUTPUT]")
          qmOutputLine = i;

        if(nLines[i] == "[EVB_HEADER]")
          evbHeaderBeginLine = i;

        if(nLines[i] == "[EVB_BIN]")
          evbBinaryLine = i;

        if(nLines[i] == "[EVB_OUTPUT]")
          evbOutputLine = i;

        if(nLines[i] == "[TMP]")
          tmpLine = i;

        if(nLines[i] == "[GEOMETRY]")
          geoBeginLine = i;

        if(nLines[i] == "[BOXXDIM]")
          boxXLine = i;

        if(nLines[i] == "[BOXYDIM]")
          boxYLine = i;

        if(nLines[i] == "[BOXZDIM]")
          boxZLine = i;

        if(nLines[i] == "[COMMANDS]")
          cmdBeginLine = i;

      }

      if(nLines[i] == "[/MD_HEADER]")
        mdHeaderEndLine = i;

      if(nLines[i] == "[/FF_HEADER]")
        ffHeaderEndLine = i;

      if(nLines[i] == "[/QM_HEADER]")
        qmHeaderEndLine = i;

      if(nLines[i] == "[/EVB_HEADER]")
        evbHeaderEndLine = i;

      if(nLines[i] == "[/GEOMETRY]")
        geoEndLine = i;

      if(nLines[i] == "[/COMMANDS]")
        cmdEndLine = i;


    }


    if(logLine != -1)
      logFilePath_ = nLines[logLine + 1];

    if(mdBinaryLine != -1)
      mdBinaryFilePath_ = nLines[mdBinaryLine + 1];

    if(mdOutputLine != -1)
      mdOutputFolderPath_ = nLines[mdOutputLine + 1];

    if(ffBinaryLine != -1)
      ffBinaryFilePath_ = nLines[ffBinaryLine + 1];

    if(ffOutputLine != -1)
      ffOutputFolderPath_ = nLines[ffOutputLine + 1];

    if(qmBinaryLine != -1)
      qmBinaryFilePath_ = nLines[qmBinaryLine + 1];

    if(qmOutputLine != -1)
      qmOutputFolderPath_ = nLines[qmOutputLine + 1];

    if(evbBinaryLine != -1)
      evbBinaryFilePath_ = nLines[evbBinaryLine + 1];

    if(evbOutputLine != -1)
      evbOutputFolderPath_ = nLines[evbOutputLine + 1];  

    if(tmpLine      != -1)
      tmpFolderPath_ = nLines[tmpLine + 1];

    if((boxXLine != -1) &&
       (boxYLine != -1) &&
       (boxZLine != -1)   )
      geometry_.setBoxDimension(std::stod(nLines[boxXLine + 1]),
                                std::stod(nLines[boxYLine + 1]),
                                std::stod(nLines[boxZLine + 1]) );

    if((geoBeginLine != -1) &&
       (geoEndLine   == -1)   )
    {  
      restartFilePath = nLines[geoBeginLine + 1];
      assert(!restartFilePath.empty());
    }
    

    if((mdHeaderBeginLine != -1               ) &&
       (mdHeaderEndLine   != -1               ) &&
       (mdHeaderEndLine   >  mdHeaderBeginLine)   )
    {

      for(off_t i  = mdHeaderBeginLine + 1;
                i != mdHeaderEndLine      ;
              ++i                          )
      {

        mdHeader_.append(nLines[i]);

        if(i != mdHeaderEndLine - 1)
          mdHeader_.append("\n");

      }

    }

    if((ffHeaderBeginLine != -1               ) &&
       (ffHeaderEndLine   != -1               ) &&
       (ffHeaderEndLine   >  ffHeaderBeginLine)   )
    {

      for(off_t i  = ffHeaderBeginLine + 1;
                i != ffHeaderEndLine      ;
              ++i                          )
      {

        ffHeader_.append(nLines[i]);

        if(i != ffHeaderEndLine - 1)
          ffHeader_.append("\n");

      }

    }

    if((qmHeaderBeginLine != -1               ) && 
       (qmHeaderEndLine   != -1               ) &&
       (qmHeaderEndLine   >  qmHeaderBeginLine)   )
    {

      for(off_t i  = qmHeaderBeginLine + 1;
                i != qmHeaderEndLine      ;
              ++i                          )
      {

        qmHeader_.append(nLines[i]);

        if(i != qmHeaderEndLine - 1)
          qmHeader_.append("\n");

      }

    }

    if((evbHeaderBeginLine != -1                ) &&
       (evbHeaderEndLine   != -1                ) &&
       (evbHeaderEndLine   >  evbHeaderBeginLine)   )
    {

      for(off_t i  = evbHeaderBeginLine + 1;
                i != evbHeaderEndLine      ;
              ++i                           )
      {

        evbHeader_.append(nLines[i]);

        if(i != evbHeaderEndLine - 1)
          evbHeader_.append("\n");

      }

    }

    if((cmdBeginLine != -1          ) &&
       (cmdEndLine   != -1          ) &&
       (cmdEndLine    > cmdBeginLine)   )
    {

      for(off_t i  = cmdBeginLine + 1;
                i != cmdEndLine      ;
              ++i                     )
      {

        std::stringstream cmdStream(nLines[i]);

        std::vector<std::string> cmdLine;
        std::string buffer;
        while(cmdStream >> buffer)
          cmdLine.push_back(buffer);

        nCommands_.push_back(cmdLine);

      }

    }

    if((geoBeginLine != -1          ) &&
       (geoEndLine   != -1          ) &&
       (geoEndLine   >  geoBeginLine)   )
    {

      for(off_t i  = geoBeginLine + 1;
                i != geoEndLine      ;
              ++i                     )
      {

        std::stringstream coordStream(nLines[i]);

        std::vector<std::string> coordLine;
        std::string buffer;
        while(coordStream >> buffer)
          coordLine.push_back(buffer);

        nCoordinates.push_back(coordLine);

      }

    }

    else if((geoBeginLine != -1) &&
            (geoEndLine   == -1)   )
    {


      std::ifstream restartFile;

      restartFile.open(restartFilePath.c_str());

      if(restartFile.bad())
      {

        std::cerr << "Qreg restart file not available at "
                  << restartFilePath
                  << ", aborting."
                  << std::endl                            ;

        exit(-1);

      }

      nLines.clear();
      if(restartFile.is_open())
      {

        while(std::getline(restartFile,
                           line        ))
          nLines.push_back(line);

        restartFile.close();

      }

      assert(!nLines.empty());

      std::vector<off_t> nStartLines;
      for(off_t i  = 0                                  ;
                i != static_cast<signed> (nLines.size());
              ++i                                        )
        if((nLines[i].size  (  ) >= 4     ) &&
           (nLines[i].substr(0,
                             4 ) == "Qreg")   )
          nStartLines.push_back(i);

      if(nStartLines.size() < 2)
      {

        std::cerr << "Qreg restart file malformated, aborting."
                  << std::endl                                 ;

        exit(-1);

      }

      off_t newCoordEndLine   = nLines.size     ()    ;
      off_t newCoordBeginLine = nStartLines.back() + 1;
      off_t oldCoordEndLine   = nStartLines.back() - 1;
      nStartLines.pop_back();
      off_t oldCoordBeginLine = nStartLines.back() + 1;

      for(auto i = oldCoordBeginLine;
               i < oldCoordEndLine  ;
             ++i                     )
      {
       
        std::stringstream coordStream(nLines[i]);

        std::vector<std::string> coordLine;
        std::string buffer;
        while(coordStream >> buffer)
          coordLine.push_back(buffer);

        nCoordinates.push_back(coordLine);

      }

      for(auto i = newCoordBeginLine;
               i < newCoordEndLine  ;
             ++i                     )
      {
        
        std::stringstream coordStream(nLines[i]);

        std::vector<std::string> coordLine;
        std::string buffer;
        while(coordStream >> buffer)
          coordLine.push_back(buffer);

        nNewCoordinates.push_back(coordLine);

      }

      if(nCoordinates.empty())
        nCoordinates = nNewCoordinates;

    }

    for(size_t i  = 0                  ;
               i != nCoordinates.size();
             ++i                        )
    {

      if(nCoordinates[i][0] == "Li")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             3                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "B")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             5                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "C")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             6                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "N")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             7                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "O")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             8                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "F")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             9                 , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "Na")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             11                , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "Mg")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             12                , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "Al")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             13                , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "Si")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             14                , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "P")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             15                , 
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "S")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             16                , 
                             std::stod(nCoordinates[i][1]), 
                             std::stod(nCoordinates[i][2]), 
                             std::stod(nCoordinates[i][3]) );

      else if(nCoordinates[i][0] == "Pt")
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             78                , 
                             std::stod(nCoordinates[i][1]), 
                             std::stod(nCoordinates[i][2]), 
                             std::stod(nCoordinates[i][3])  );

      else
        geometry_.insertAtom(i + 1             ,
                             nCoordinates[i][0],
                             1                 ,
                             std::stod(nCoordinates[i][1]),
                             std::stod(nCoordinates[i][2]),
                             std::stod(nCoordinates[i][3]) );


    }

    for(size_t i  = 0                     ;
               i != nNewCoordinates.size();
             ++i                           )
      geometry_.setCoordinate(i + 1,
                              std::stod(nNewCoordinates[i][1]),
                              std::stod(nNewCoordinates[i][2]),
                              std::stod(nNewCoordinates[i][3]) );



    //
    // Create folders.
    //

    if(!mdOutputFolderPath_.empty())
    {

      std::string cmd = "mkdir ";
      cmd.append(mdOutputFolderPath_);
      cmd.append(" > /dev/null 2>&1");
      if((std::system(cmd.c_str()) == -1    ) &&
         (errno                    != EEXIST)   )
      {

        std::cerr << "Unable to create folder for MD output."
                  << std::endl                               ;

        exit(-1);

      }

    }

    if(!ffOutputFolderPath_.empty())
    {
      std::string cmd = "mkdir ";
      cmd.append(ffOutputFolderPath_);
      cmd.append(" > /dev/null 2>&1");
      if((std::system(cmd.c_str()) == -1    ) &&
         (errno                    != EEXIST)   )
      {

        std::cerr << "Unable to create folder for force field output."
                  << std::endl                                        ;

        exit(-1);

      }
    }

    if(!qmOutputFolderPath_.empty())
    {

      std::string cmd = "mkdir ";
      cmd.append(qmOutputFolderPath_);
      cmd.append(" > /dev/null 2>&1");
      if((std::system(cmd.c_str()) == -1    ) &&
         (errno                    != EEXIST)   )
      {

        std::cerr << "Unable to create folder for QM output."
                  << std::endl                          ;

        exit(-1);

      }

    }

    if(!evbOutputFolderPath_.empty())
    {

      std::string cmd = "mkdir ";
      cmd.append(evbOutputFolderPath_);
      cmd.append(" > /dev/null 2>&1" );
      if((std::system(cmd.c_str()) == -1    ) &&
         (errno                    != EEXIST)   )
      {

        std::cerr << "Unable to create folder for EVB output."
                  << std::endl                                ;

        exit(-1);

      }

    }

    if(!tmpFolderPath_.empty())
    {

      std::string cmd = "mkdir ";
      cmd.append(tmpFolderPath_     );
      cmd.append(" > /dev/null 2>&1");
      if((std::system(cmd.c_str()) == -1    ) &&
         (errno                    != EEXIST)   )
      {

        std::cerr << "Unable to create folder for temporary files."
                  << std::endl                                     ;

        exit(-1);

      }

    }


    //
    // DbC POST
    //
    assert(isSane());


  }



  auto
  Config::getLogFilePath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getLogFilePath(void)const noexcept
    // ->std::string
    //

    return logFilePath_;


  }



  auto
  Config::getMdBinaryFilePath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getMdBinaryFilePath(void)const noexcept
    // ->std::string
    //
 
    return mdBinaryFilePath_;


  }



  auto
  Config::getMdOutputFolderPath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getMdOutputFolderPath(void)const noexcept
    // ->std::string
    //

    return mdOutputFolderPath_;


  }



  auto
  Config::getFfBinaryFilePath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getFfBinaryFilePath(void)const noexcept
    // ->std::string
    //

    return ffBinaryFilePath_;


  }



  auto
  Config::getFfOutputFolderPath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getFfOutputFolderPath(void)const noexcept
    // ->std::string
    //

    return ffOutputFolderPath_;


  }



  auto
  Config::getQmBinaryFilePath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getQmBinaryFilePath(void)const noexcept
    // ->std::string
    //

    return qmBinaryFilePath_;


  }



  auto
  Config::getQmOutputFolderPath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getQmOutputFolderPath(void)const noexcept
    // ->std::string
    //

    return qmOutputFolderPath_;


  }



  auto
  Config::getEvbBinaryFilePath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getEvbBinaryFilePath(void)const noexcept
    // ->std::string
    //

    return evbBinaryFilePath_;


  }



  auto
  Config::getEvbOutputFolderPath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getEvbOutputFolderPath(void)const noexcept
    // ->std::string
    //

    return evbOutputFolderPath_;


  }



  auto
  Config::getTmpFolderPath(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getTmpFolderPath(void)const noexcept
    // ->std::string
    //

    return tmpFolderPath_;


  }



  auto
  Config::getMdHeader(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getMdHeader(void)const noexcept
    // ->std::string
    //

    return mdHeader_;


  }



  auto
  Config::getFfHeader(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getFfHeader(void)const
    // ->std::string
    //

    return ffHeader_;


  }



  auto
  Config::getQmHeader(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getQmHeader(void)const noexcept
    // ->std::string
    //

    return qmHeader_;


  }



  auto
  Config::getEvbHeader(void)const noexcept
  ->std::string
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getEvbHeader(void)const noexcept
    // ->std::string
    //

    return evbHeader_;


  }



  auto
  Config::getGeometry(void)const noexcept
  ->Box
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getGeometry(void)
    // ->Box
    //

    return geometry_;


  }



  auto
  Config::getNCommands(void)const noexcept
  ->std::vector< std::vector<std::string> >
  {


    //
    // DbC PRE
    //
    assert(isSane());


    //
    // auto
    // Config::getNCommands(void)cons
    //  ->std::vector< std::vector<std::string> >
    //

    return nCommands_;


  }



}
