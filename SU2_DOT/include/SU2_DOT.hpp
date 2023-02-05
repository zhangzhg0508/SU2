/*!
 * \file SU2_DOT.hpp
 * \brief Headers of the main subroutines of the code SU2_DOT.
 * \author F. Palacios, T. Economon
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../Common/include/parallelization/omp_structure.hpp"
#include "../../SU2_DEF/include/drivers/CDiscAdjDeformationDriver.hpp"

using namespace std;
