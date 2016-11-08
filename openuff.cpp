//
// openuff.cpp
//
//  Created on: Nov 18,2015
//      Author: Sebastian Dohm <sebastian.dohm@uni-ulm.de>
//

/*
 * Copyright 2015 Sebastian Dohm
 * 
 * This file is part of OpenUFF.
 *
 * OpenUFF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenUFF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenUFF.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "openbabel/obconversion.h"
#include "openbabel/forcefield.h"
#include <cassert>
#include <iostream>
int main() {
  const auto conversionFactor = -1 * 694.77 / 4.184;
  std::string fNameIn = "uff.xyz";
  std::string fNameOut = "uff.grad";
  std::ifstream fIn(fNameIn);
  if(!fIn) {std::cerr << "Inputfile 'uff.xyz' unavailable." << std::endl; exit(-1);}
  if(fIn.is_open())
  {
    std::ofstream fOut(fNameOut       ,
                       std::ios::trunc );
    if(!fOut) {std::cerr << "Outputfile 'uff.grad' unavailable." << std::endl; exit(-1);}
    if(fOut.is_open())
    {
      auto mol  = OpenBabel::OBMol       ();
      auto conv = OpenBabel::OBConversion();
      const auto format = conv.FormatFromExt(fNameIn);
      if(!conv.SetInAndOutFormats(format,
                                  format )) {std::cerr << "Unknown file format." << std::endl; exit(-1);}
      conv.ReadFile(&mol   ,
                    fNameIn );
      auto ff = OpenBabel::OBForceField::FindForceField("UFF");
      if(!ff) {std::cerr << "Forcefield unavailable." << std::endl; exit(-1);}
      if(!ff->Setup(mol)) {std::cerr << "Forcefield setup failed." << std::endl; exit(-1);}
      ff->SteepestDescentInitialize();
      if(!ff->GetCoordinates(mol)) {std::cerr << "Forcefield geometry optimization failed." << std::endl; exit(-1);}
      const auto data = static_cast<OpenBabel::OBConformerData*>
       (mol.GetData(OpenBabel::OBGenericDataType::ConformerData));
      if(!data) {std::cerr << "No data received." << std::endl; exit(-1);}
      assert(data->GetForces().size() == 1);
      const auto nForces = data->GetForces()[0];
      fOut << "OpenUFF GRADIENT  ~  x  y  z  ~  [ E-23 J / Angstrom ]"
           << std::endl                                               ;
      for(const auto& force:nForces)
      {
        fOut << std::to_string(force.x() * conversionFactor)
             << "   "
             << std::to_string(force.y() * conversionFactor)
             << "   "
             << std::to_string(force.z() * conversionFactor)
             << std::endl                                   ;
      }
      fOut << "OpenUFF ended normally";
      fOut.close(); 
    }
    fIn.close();
  }
}
