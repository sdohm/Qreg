//
// ifuff.cpp
//
//  Created on: Mar 10,2015
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



#include "ifuff.hpp"



namespace qreg
{



  IfUff::IfUff(std::string newBinaryFilePath  ,
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



  IfUff::IfUff(IfUff&& ifUff)
  {


    // Call move assignment operator.
    *this = std::move(ifUff);


  }



  auto
  IfUff::operator=(IfUff&& ifUff)
  ->IfUff&
  {


    //
    // Perform move operation.
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    ifUff.~IfUff();


    return *this;


  }



  IfUff::~IfUff(void)
  {



  }



  auto
  IfUff::getNormalTerminationMessage(void)const noexcept
  ->std::string
  {


    //
    // auto
    // IfUff::getNormalTerminationMessage(void)const noexcept
    // ->std::string
    //

    return "UFF ended normally";


  }



  auto
  IfUff::harvest(Box& box)
  ->bool
  {


    //
    // DbC Invariant
    //
    assert(processId_ > 0);


    //
    // auto
    // IfUff::harvest(Box& box)
    // ->bool
    //

    std::lock_guard<std::mutex> lock(threadSafety_);

    if(!hasFinished())
      waitpid(processId_,
              NULL      ,
              0          );

    if(hasTerminatedNormally())
    {

      outputFilePath_ = outputFilePath_.substr(0                         ,
                                               outputFilePath_.size() - 4 );
      outputFilePath_.append("gradient");

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

      std::istringstream iss(nLines[1]);
      std::vector<std::string> tokens
      {
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{   }
      };

      box.insertEPot(quantumRegionNo_                                  ,
                     std::stod(tokens[6]) * getEnergyConversionFactor() );

      box.setTotalEPot(std::stod(tokens[6]) * getEnergyConversionFactor());

      for(size_t i  = (((nLines.size() - 3) / 2) + 2);
                 i != nLines.size() - 1              ;
               ++i                                    )
      {

        std::istringstream ss(nLines[i]);
        std::vector<std::string> gradientGuess
        {
          std::istream_iterator<std::string>{ss},
          std::istream_iterator<std::string>{  }
        };

        for(auto&& str:gradientGuess)
          std::replace(str.begin(),
                       str.end  (),
                       'D'        ,
                       'e'         );
        if((!std::isinf(std::stod(gradientGuess[0]))                  ) &&
           (!std::isinf(std::stod(gradientGuess[1]))                  ) &&
           (!std::isinf(std::stod(gradientGuess[2]))                  ) &&
           (std::stod(gradientGuess[0]) == std::stod(gradientGuess[0])) &&
           (std::stod(gradientGuess[1]) == std::stod(gradientGuess[1])) &&
           (std::stod(gradientGuess[2]) == std::stod(gradientGuess[2])) &&
           (box.hasAtom(idConverter_.getInternalId(i -
                             (((nLines.size() - 3) / 2) + 2)))        )   )
          box.setGradientGuess(idConverter_.getInternalId(i -
                             (((nLines.size() - 3) / 2) + 2)),
                               std::stod(gradientGuess[0]) *
                                getGradientConversionFactor(),
                               std::stod(gradientGuess[1]) *
                                getGradientConversionFactor(),
                               std::stod(gradientGuess[2]) *
                                getGradientConversionFactor() );

      }

      outputFilePath_ = outputFilePath_.substr(0                         ,
                                               outputFilePath_.size() - 8 );
      outputFilePath_.append(".out");

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
        box.setGradientGuess(atomNo   ,
                             dist(rng),
                             dist(rng),
                             dist(rng) );
    }


    return false;


  }



  auto
  IfUff::seed(const Box&          box            ,
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
    // IfUff::seed(const Box&          box            ,
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

        inputText_.append(std::to_string(clone.getXCoordinate(nAtomNos[i])));
        inputText_.append(" "                                              );
        inputText_.append(std::to_string(clone.getYCoordinate(nAtomNos[i])));
        inputText_.append(" "                                              );
        inputText_.append(std::to_string(clone.getZCoordinate(nAtomNos[i])));
        inputText_.append("\n"                                             );

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

        inputText_.append(std::to_string(box.getXCoordinate(nAtomNos[i])));
        inputText_.append(" "                                            );
        inputText_.append(std::to_string(box.getYCoordinate(nAtomNos[i])));
        inputText_.append(" "                                            );
        inputText_.append(std::to_string(box.getZCoordinate(nAtomNos[i])));
        inputText_.append("\n"                                           );

      }

    }


    //
    // Input filename.
    //

    outputFilePath_ = outputFolderPath_;
    outputFilePath_.append("/uff.out");

    std::string inputFilePath = outputFolderPath_;

    std::string purge = "rm -f ";
    purge.append(inputFilePath);
    purge.append("/control "      );
    purge.append(inputFilePath);
    purge.append("/coord "        );
    purge.append(inputFilePath);
    purge.append("/uffenergy "    );
    purge.append(inputFilePath);
    purge.append("/uffgradient "  );
    purge.append(inputFilePath);
    purge.append("/uffgradx "     );
    purge.append(inputFilePath);
    purge.append("/uffhessian0-0 ");
    purge.append(inputFilePath);
    purge.append("/uff.out "      );
    purge.append(inputFilePath);
    purge.append("/ufftopology "  );

    inputFilePath.append("/qreg_"                       );
    inputFilePath.append(std::to_string(quantumRegionNo));
    inputFilePath.append("_uff.in"                      );

    purge.append(inputFilePath);

    if(std::system(purge.c_str()) != 0)
    {

      std::cerr << "Unable to purge old UFF files, aborting Qreg."
                << std::endl                                      ;

      exit(-1);

    }


    //
    // Write to UFF input file.
    //

    std::ofstream inputFile;
    inputFile.open(inputFilePath.c_str());

    if(inputFile.good())
      inputFile << inputText_;

    inputFile.close();


    //
    // Setup system call.
    //

    std::string systemCall = "xyz2coord ";
    systemCall.append(inputFilePath    );
    systemCall.append(" >> "           );
    systemCall.append(outputFolderPath_);
    systemCall.append("/coord"         );

    if(std::system(systemCall.c_str()) != 0)
    {

      std::cerr << "XYZ reporting error state, aborting Qreg."
                << std::endl                                  ;

      exit(-1);

    }

    
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

      if(chdir(outputFolderPath_.c_str()) != 0)
      {

        std::cerr << "CHDIR reporting error state, aborting Qreg."
                  << std::endl                                    ;

        exit(-1);

      }

      systemCall.clear();
      systemCall = "exec ";
      systemCall.append(binaryFilePath_);
      systemCall.append(" > uff.out"   );

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
