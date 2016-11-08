//
// iforc.cpp
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



#include "iforc.hpp"



namespace qreg
{



  IfOrc::IfOrc(std::string newBinaryFilePath  ,
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



  IfOrc::IfOrc(IfOrc&& ifOrc)
  {


    // Call move assignment operator.
    *this = std::move(ifOrc);


  }



  auto
  IfOrc::operator=(IfOrc&& ifOrc)
  ->IfOrc&
  {


    //
    // Perform move operation.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    ifOrc.~IfOrc();


    return *this;


  }



  IfOrc::~IfOrc(void)
  {



  }



  auto
  IfOrc::getNormalTerminationMessage(void)const noexcept
  ->std::string
  {


    //
    // auto
    // IfOrc::getNormalTerminationMessage(void)const noexcept
    // ->std::string
    //

    return "                             ****ORCA TERMINATED NORMALLY****";


  }



  auto
  IfOrc::harvest(Box& box)
  ->bool
  {


    //
    // DbC Invariant
    //
    assert(processId_ > 0);


    //
    // auto
    // IfOrc::harvest(Box& box)
    // ->bool
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(!hasFinished())
      waitpid(processId_,
              NULL      ,
              0          );

    if(hasTerminatedNormally())
    {


      size_t energyLine      = 0;
      size_t gradBeginLine   = 0;
      size_t gradEndLine     = 0;
      size_t chargeBeginLine = 0;
      size_t chargeEndLine   = 0;

      double ePot = 0.0;


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

      for(decltype(nLines.size()) i  = 0            ;
                                  i != nLines.size();
                                ++i                  )
      {

        if(nLines[i].length() >= 5)
          if(nLines[i].substr(0,
                              5 ) == getEnergyMessage())
            energyLine = i;

        if(nLines[i] == getQmGradientBeginMessage())
          gradBeginLine = i + 3;

        if(nLines[i] == getSeGradientBeginMessage())
          gradBeginLine = i + 1;

        if(nLines[i].length() >= 3)
          if(nLines[i].substr(0,
                              4 ) == getGradientEndMessage())
            gradEndLine = i - 1;

        if(nLines[i] == getChargeBeginMessage())
          chargeBeginLine = i + 2;

        if(nLines[i] == getChargeEndMessage())
          chargeEndLine = i - 2;

      }

      if(energyLine != 0)
        ePot = std::stod(nLines[energyLine].substr(31,
                                                   10 ));

      box.insertEPot(quantumRegionNo_                  ,
                     ePot * getEnergyConversionFactor() );

      box.setTotalEPot(ePot * getEnergyConversionFactor());

      std::vector <std::vector<std::string> > nGradients;

      for(auto i = gradBeginLine;
               i < gradEndLine  ;
             ++i                 )
      {

        std::stringstream        gradStream(nLines[i]);
        std::string              buffer               ;
        std::vector<std::string> gradLine             ;

        while(gradStream >> buffer)
          gradLine.push_back(buffer);

        nGradients.push_back(gradLine);

      }

      std::vector <std::vector<std::string> > nCharges;

      for(auto i = chargeBeginLine;
               i < chargeEndLine  ;
             ++i                    )
      {

        std::string              buffer    ;
        std::vector<std::string> chargeLine;

        std::stringstream chargeStream(nLines[i]);

        while(chargeStream >> buffer)
          chargeLine.push_back(buffer);

        nCharges.push_back(chargeLine);

      }

      for(decltype(nGradients.size()) i  = 0                ;
                                      i != nGradients.size(); 
                                    ++i                      )
      {

        assert(nGradients[i].size() >= 6);
        if((!std::isinf(std::stod(nGradients[i][3]))                  ) &&
           (!std::isinf(std::stod(nGradients[i][4]))                  ) &&
           (!std::isinf(std::stod(nGradients[i][5]))                  ) &&
           (std::stod(nGradients[i][3]) == std::stod(nGradients[i][3])) &&
           (std::stod(nGradients[i][4]) == std::stod(nGradients[i][4])) &&
           (std::stod(nGradients[i][5]) == std::stod(nGradients[i][5])) &&
           (box.hasAtom(idConverter_.getInternalId(i))                )   )
          box.setGradient(idConverter_.getInternalId(i) ,
                          std::stod(nGradients[i][3]) *
                           getGradientConversionFactor(),
                          std::stod(nGradients[i][4]) *
                           getGradientConversionFactor(),
                          std::stod(nGradients[i][5]) *
                           getGradientConversionFactor() );

      }

      for(decltype(nCharges.size()) i  = 0              ; 
                                    i != nCharges.size();
                                  ++i                    )
      {

        assert(nCharges[i].size() >= 3);
        if(nCharges[i][2].compare(":") != 0)
          box.setPartialCharge(idConverter_.getInternalId(i),
                               std::stod(nCharges[i][2])     );
        else
        {
          assert(nCharges[i].size() >= 4);
          box.setPartialCharge(idConverter_.getInternalId(i),
                               std::stod(nCharges[i][3])     );
        }

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
  IfOrc::seed(const Box&          box            ,
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
    // IfOrc::seed(const Box&          box            ,
    //             const unsigned int& quantumRegionNo,
    //             const signed int&   charge         ,
    //             const unsigned int& multiplicity   ,
    //             const double&       range          ,
    //             const double&       distanceFactor   )
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

      inputText_ = headerText_;

      inputText_.append(" "                         );
      inputText_.append(std::to_string(charge      ));
      inputText_.append(" "                         );
      inputText_.append(std::to_string(multiplicity));
      inputText_.append("\n"                        );

      const auto nAtomNos = box.getNAtomNos();
      assert(!nAtomNos.empty());

      for(decltype(nAtomNos.size()) i  = 0              ;
                                    i != nAtomNos.size();
                                  ++i                    )
      {

        idConverter_.insertIds(nAtomNos[i],
                                        i  );

        inputText_.append(clone.getChemicalSymbol(nAtomNos[i]));
        inputText_.append(" "                                 );

        // non-atomic charges
        if(clone.getChemicalSymbol(nAtomNos[i]) == "Q")
        {
          inputText_.append(std::to_string(clone.
                                           getPartialChargeGuess(nAtomNos[i])));
          inputText_.append(" "                                               );
        }

        inputText_.append(std::to_string(clone.getXCoordinate(nAtomNos[i])));
        inputText_.append(" "                                              );
        inputText_.append(std::to_string(clone.getYCoordinate(nAtomNos[i])));
        inputText_.append(" "                                              );
        inputText_.append(std::to_string(clone.getZCoordinate(nAtomNos[i])));
        inputText_.append("\n"                                             );

      }

      inputText_.append("*");

    }

    else
    {

      Box clone = box.cloneShrinkGraphite(range);

      inputText_ = headerText_;

      inputText_.append(" "                         );
      inputText_.append(std::to_string(charge      ));
      inputText_.append(" "                         );
      inputText_.append(std::to_string(multiplicity));
      inputText_.append("\n"                        );

      const auto nAtomNos = box.getNAtomNosByQuantumRegionNo(quantumRegionNo_);

      assert(!nAtomNos.empty());

      for(decltype(nAtomNos.size()) i  = 0              ;
                                    i != nAtomNos.size();
                                  ++i                    )
      {

        idConverter_.insertIds(nAtomNos[i],
                                        i  );

        inputText_.append(clone.getChemicalSymbol(nAtomNos[i]));
        inputText_.append(" "                                 );

        // non-atomic charges
        if(clone.getChemicalSymbol(nAtomNos[i]) == "Q")
        {
          inputText_.append(std::to_string(clone.
                                           getPartialChargeGuess(nAtomNos[i])));
          inputText_.append(" "                                               );
        }

        inputText_.append(std::to_string(clone.getXDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append(" "                                              );
        inputText_.append(std::to_string(clone.getYDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append(" "                                             );
        inputText_.append(std::to_string(clone.getZDistance(nAtomNos[0],
                                                            nAtomNos[i] )));
        inputText_.append("\n"                                            );

      }

      inputText_.append("*");

    }


    //
    // Input filename.
    //

    std::string inputFilePath = outputFolderPath_;

    inputFilePath.append("/qreg_"                        );
    inputFilePath.append(std::to_string(quantumRegionNo_));
    inputFilePath.append("_orc."                         );

    std::string purge = "rm -f ";
    purge.append(inputFilePath);
    purge.append("*"          );
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

    outputFilePath_ = outputFolderPath_;

    outputFilePath_.append("/qreg_"                        );
    outputFilePath_.append(std::to_string(quantumRegionNo_));
    outputFilePath_.append("_orc.out"                      );


    //
    // Setup system call.
    //

    std::string systemCall = "exec ";

    systemCall.append(binaryFilePath_);
    systemCall.append(" "            );
    systemCall.append(inputFilePath  );
    systemCall.append(" > "          );
    systemCall.append(outputFilePath_);


    //
    // Write to ORC input file.
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

