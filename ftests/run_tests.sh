#!/usr/bin/env -S bash
# -*- Mode:Shell-script; Coding:us-ascii-unix; fill-column:158 -*-
#########################################################################################################################################################.H.S.##
##
# @file      run_tests.sh
# @author    Mitch Richling http://www.mitchr.me/
# @date      2024-08-13
# @brief     Run all functional tests.@EOL
# @std       bash
# @see       
# @copyright 
#  @parblock
#  Copyright (c) 2024, Mitchell Jay Richling <http://www.mitchr.me/> All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
#  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  @endparblock
# @filedetails
#
#  NOTE: Run from the build directory!
#
#  WARNING: I developed this bash script for MSYS2 on WINDOWS!  It should work on UNIX'ish systems, but I have not tested it.
#
#########################################################################################################################################################.H.E.##

DO_BUILD='NO'
DO_RUN='YES'

FDIFF_PATH=~/world/my_prog/utils/float_diff.rb


if [ "$(basename $(pwd))" == "build" -a -e '../CMakeLists.txt' ]; then
  if [ "$DO_BUILD" == 'YES' ]; then
    echo " "
    echo " "
    echo "Building Tests Now"
    echo " "
    for f in ../ftests/*.vtu; do
      t=$(echo $(basename $f) | sed 's/\.vtu$//')
      make "$t"
    done
  fi
  if [ "$DO_RUN" == 'YES' ]; then
    echo " "
    echo " "
    echo "Running Tests Now"
    echo " "
    for f in ../ftests/*.vtu; do
      t=$(echo $(basename $f) | sed 's/\.vtu$//')
      b="$t"
      if [ -x "$b"'.exe' ]; then
        b="$t"'.exe'
      fi
      if [ -x "$b" ]; then
        ./"$b" > /dev/null
        if ruby "$FDIFF_PATH" "$t".vtu ../ftests/"$t".vtu >/dev/null; then
          echo "PASS: $t"
        else
          echo "FAIL: $t (invalid output)"
        fi
      else
        echo "FAIL: $t (no executable)"
      fi
    done
  fi
  echo " "
else
  echo "ERROR(run_tests.sh): Run this script from the build directory!"
fi
