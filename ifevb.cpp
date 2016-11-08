//
// ifevb.cpp
//
//  Created on: Jul 1,2015
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



#include "ifevb.hpp"



namespace qreg
{



  IfEvb::IfEvb(std::string newBinaryFilePath  ,
               std::string newHeaderText      ,
               std::string newOutputFolderPath )noexcept:

    IfBase(newBinaryFilePath  ,
           newHeaderText      ,
           newOutputFolderPath )

  {


    //
    // DbC PRE
    //

    assert(!binaryFilePath_  .empty());
    assert(!headerText_      .empty());
    assert(!outputFolderPath_.empty());


  }



  IfEvb::IfEvb(IfEvb&& ifEvb)
  {


    // Call move assignment operator.
    *this = std::move(ifEvb);


  }



  auto
  IfEvb::operator=(IfEvb&& ifEvb)
  ->IfEvb&
  {


    //
    // Perform move operation.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    ifEvb.~IfEvb();


    return *this;


  }



  IfEvb::~IfEvb(void)
  {



  }



  auto
  IfEvb::getNormalTerminationMessage(void)const noexcept
  ->std::string
  {


    //
    // auto
    // IfEvb::getNormalTerminationMessage(void)const noexcept
    // ->std::string
    //

    return " EVB DONE";


  }



  auto
  IfEvb::harvest(Box& box)
  ->bool
  {


    //
    // DbC Invariant
    //
    assert(processId_ > 0);


    //
    // auto
    // IfEvb::harvest(Box& box)
    // ->bool
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(!hasFinished())
      waitpid(processId_,
              NULL      ,
              0          );

    if(hasTerminatedNormally())
    {


      //
      // Read output file.
      //

      std::ifstream outputFile;
      outputFile.open(outputFilePath_.c_str());

      std::string               line  ;
      std::vector <std::string> nLines;

      if(outputFile.is_open())
      {

        while(std::getline(outputFile,
                           line       ))
          nLines.push_back(line);

        outputFile.close();

      }


      //
      // Extract data.
      //

      const double ePot = std::stod(nLines[1]);

      box.insertEPot(quantumRegionNo_                  ,
                     ePot * getEnergyConversionFactor() );

      box.setTotalEPot(ePot * getEnergyConversionFactor());

      std::vector <std::vector<std::string> > nGradients;

      for(decltype(nLines.size()) i = 2                ;
                                  i < nLines.size() - 1;
                                ++i                     )
      {

        std::stringstream        gradStream(nLines[i]);
        std::string              buffer               ;
        std::vector<std::string> gradLine             ;

        while(gradStream >> buffer)
          gradLine.push_back(buffer);

        nGradients.push_back(gradLine);

      }

      for(decltype(nGradients.size()) i  = 0                ;
                                      i != nGradients.size(); 
                                    ++i                      )
      {

        assert(nGradients[i].size() == 3);
        if((std::stod(nGradients[i][0]) != 0.0) ||
           (std::stod(nGradients[i][1]) != 0.0) ||
           (std::stod(nGradients[i][2]) != 0.0)   )
          box.setGradient(idConverter_.getInternalId(i) ,
                          std::stod(nGradients[i][0]) *
                           getGradientConversionFactor(),
                          std::stod(nGradients[i][1]) *
                           getGradientConversionFactor(),
                          std::stod(nGradients[i][2]) *
                           getGradientConversionFactor() );

      }

      hasBeenHarvested_ = true;

      return hasBeenHarvested();


    }

    else
    {
      // Add random noise to the gradient.
      std::uniform_real_distribution<double> dist(-100.0,
                                                   100.0 );
      std::mt19937 rng;
      rng.seed(std::random_device{}());
      for(const auto& atomNo:box.getNAtomNosByQuantumRegionNo(quantumRegionNo_))
        box.setGradient(atomNo   ,
                        dist(rng),
                        dist(rng),
                        dist(rng) );
    }

    return false;


  }



  auto
  IfEvb::seed(const Box&          box            ,
              const unsigned int& quantumRegionNo,
              const signed int&   charge         ,
              const unsigned int& multiplicity   ,
              const double&       range           )
  ->void
  {


    //
    // DbC PRE
    //
    assert(quantumRegionNo != 0);


    //
    // auto
    // IfEvb::seed(const Box&          box            ,
    //             const unsigned int& quantumRegionNo,
    //             const signed int&   charge         ,
    //             const unsigned int& multiplicity   ,
    //             const double&       range          ,
    //             const double&       distanceFactor  )
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    quantumRegionNo_ = quantumRegionNo;


    //
    // Setup input text.
    //

    if(box.getNAtomNosByQuantumRegionNo(quantumRegionNo_).size() ==
       box.getNAtomNos()                                 .size()   )
    {

      Box clone = box.cloneNonPeriodic(range);

      const auto nAtomNos = clone.getNAtomNos();
      assert(!nAtomNos.empty());

      inputText_ = std::to_string(nAtomNos.size());

      inputText_.append("\n"                        );
      inputText_.append(headerText_                 );
      inputText_.append(" quasi periodic "          );
      inputText_.append(std::to_string(charge      ));
      inputText_.append(" "                         );
      inputText_.append(std::to_string(multiplicity));
      inputText_.append("\n"                        );

      for(decltype(nAtomNos.size()) i  = 0              ;
                                    i != nAtomNos.size();
                                  ++i                    )
      {

        idConverter_.insertIds(nAtomNos[i],
                                        i  );

        inputText_.append(clone.getChemicalSymbol(nAtomNos[i]));
        inputText_.append(" "                                 );

        inputText_.append(std::to_string(clone.getXDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append(" "                                             );
        inputText_.append(std::to_string(clone.getYDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append(" "                                             );
        inputText_.append(std::to_string(clone.getZDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append("\n"                                            );

      }

    }

    else
    {

      const auto nAtomNos = box.getNAtomNosByQuantumRegionNo(quantumRegionNo_);
      assert(!nAtomNos.empty());

      inputText_ = std::to_string(nAtomNos.size());

      inputText_.append("\n"                        );
      inputText_.append(headerText_                 );
      inputText_.append(" "                         );
      inputText_.append(std::to_string(charge      ));
      inputText_.append(" "                         );
      inputText_.append(std::to_string(multiplicity));
      inputText_.append("\n"                        );

      for(decltype(nAtomNos.size()) i  = 0              ;
                                    i != nAtomNos.size();
                                  ++i                    )
      {

        idConverter_.insertIds(nAtomNos[i],
                                        i  );

        inputText_.append(box.getChemicalSymbol(nAtomNos[i]));
        inputText_.append(" "                               );

        inputText_.append(std::to_string(box.getXDistance(nAtomNos[0],
                                                          nAtomNos[i] )));
        inputText_.append(" "                                           );
        inputText_.append(std::to_string(box.getYDistance(nAtomNos[0],
                                                          nAtomNos[i] )));
        inputText_.append(" "                                           );
        inputText_.append(std::to_string(box.getZDistance(nAtomNos[0],
                                                          nAtomNos[i] )));
        inputText_.append("\n"                                          );

      }

    }


    //
    // Input filename.
    //

    std::string inputFilePath = outputFolderPath_;

    inputFilePath.append("/qreg_"                        );
    inputFilePath.append(std::to_string(quantumRegionNo_));
    inputFilePath.append("_evb."                         );

    std::string purge = "rm -f ";
    purge.append("evb.tmp" );
    purge.append(" evb.out");
    if(std::system(purge.c_str()) != 0)
    {

      std::cerr << "Unable to call RM, aborting Qreg."
                << std::endl                          ;

      exit(-1);

    }

    inputFilePath.append("in");


    //
    // Output filename.
    //

    outputFilePath_ = "evb.out";


    //
    // Setup system call.
    //

    std::string systemCall = "exec ";

    systemCall.append(binaryFilePath_);
    systemCall.append(" -f "         );
    systemCall.append(inputFilePath  );


    //
    // Write to EVB input file.
    //

    std::ofstream inputFile;
    inputFile.open(inputFilePath.c_str());

    if(inputFile.good())
      inputFile << inputText_;

    inputFile.close();


    //
    // Execute system call.
    //

    pid_t pid;
    if((pid = fork()) == -1)
    {
      perror("fork");
      exit(-1);
    }

    if(pid == 0) // CHILD
    {
      execlp("/bin/sh"         ,
             "/bin/sh"         ,
             "-c"              ,
             systemCall.c_str(),
             (char*)NULL        );
      perror("exec");
      exit(-1);
    }

    else // PARENT
    {
      // Retrieve child process ID.
      processId_ = pid;
    }


  }




}
