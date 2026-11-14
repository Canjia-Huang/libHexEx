/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "HexEx/Config/Export.hh"

namespace HexEx
{

/**
 * @brief Enables or disables the progress/debug output of the library.
 *
 * The extraction prints a fair amount of progress information to `std::cout`
 * (and some warnings to `std::cerr`). Applications that want a quiet run can
 * switch it off; it is enabled by default so that the command line tools and
 * the demos behave as before.
 */
HEXEX_EXPORT void setVerbose(bool enabled);

/// Whether the library's progress/debug output is enabled (the default).
HEXEX_EXPORT bool verbose();

} // namespace HexEx
