//
// exe.cpp
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



#include "exe.hpp"



int main(int    argc,
         char** argv )
{



  static unsigned int cycleCount = 0;

  static std::string qregInputFilePath;


  switch(argc)
  {


    case 2:
    {

      cycleCount        = 1      ;
      qregInputFilePath = argv[1];

      break;

    }

    case 3:
    {

      cycleCount        = std::stoi(argv[2]);
      qregInputFilePath =           argv[1] ;
      break;

    }

    default:
    {

      std::cerr << "SYNTAX ERROR"
                << std::endl
                << "Usage: qreg inputfile [optional: cycle count]" 
                << std::endl                                      ;

      return -1;

      break ;

    }


  }


  static std::ifstream file(qregInputFilePath);

 
  //
  // Check Input file.
  //

  file.open(qregInputFilePath.c_str());

  if(file.bad())
  {

    std::cerr << "Qreg configuration file not available at "
              << qregInputFilePath
              << ", aborting Qreg."
              << std::endl                                  ; 

    return -1;

  }


  file.close();


  //
  // Initialize box.
  //

  static std::unique_ptr<qreg::Box> box(new qreg::Box);

  //
  // Initialize Qreg.
  //

  static std::unique_ptr<qreg::Qreg> qreg(new qreg::Qreg(qregInputFilePath,
                                                         QREG_ORC         ,
                                                         1                 ));


  //
  // Run main cycle.
  //

  for(unsigned int i  = 0         ;
                   i != cycleCount;
                 ++i               )
    qreg->advance(*box);


  return 0;



}
