#!/bin/bash
 ###############################################################################
 #
 # MatchMutalyzerOutputToBatchfile helps matching the Mutalyzer result files,
 # typically named in the format "Results_<timestamp>.txt" (for example,
 # "Results_140358894954.txt"), to the corresponding input file. Especially when
 # uploading multiple files, it's easy to lose track of which Mutalyzer output
 # file belonged to which input file.
 # This script will analyze the first 20 lines of the Mutalyzer results file and
 # the files used as input, and rename the results file to match the input file,
 # with a "_results.txt" suffix.
 #
 # Created     : 2013-08-21
 # Modified    : 2013-08-21
 # Version     : 0.1
 #
 # Copyright   : 2013-2015 Leiden University Medical Center; http://www.LUMC.nl/
 # Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 #
 #
 # This work is licensed under the Creative Commons
 # Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
 # copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
 # or send a letter to:
 # Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
 #
 ##############

TOTAL_ARGS=$#
COUNT=0

# Check: have more than one argument.
if [[ -z $2 ]]
then
  echo "Usage: $0 <mutalyzer_output_files> <batch_files>";
  exit 1
fi

for file in "$@"
do
  COUNT=$(($COUNT+1))

  # Determine if this is a mutalyzer output file or not.
  if [[ `head -n 1 "${file}" | grep Variant -c | cut -d ":" -f 2` == "1" ]]
  then
    # Gather "signature" from file, which is the first 20 lines.
    FILENAMES[$COUNT]="${file}"
    SIGNATURES[$COUNT]=`head -n 21 "${file}" | tail -n 20 | cut -f 1 | tr '\r' ' ' | tr '\n' ' ' | sed 's/  / /g'` # Normally doesn't contain \r, but well...
  else
    # Make same signature, and find it in the array.
    SIG=`head -n 20 "${file}" | tr '\r' ' ' | tr '\n' ' ' | sed 's/  / /g'`

    # Loop Mutalyzer files.
    for i in ${!FILENAMES[@]} # ! is needed to get $i to contain the key, not the value.
    do
      if [[ "${SIG}" == "${SIGNATURES[$i]}" ]]
      then
        echo "${FILENAMES[$i]} matches ${file}"
        FILE_NEW=`echo "${file}" | sed 's/.txt/_results.txt/'`
        mv -n "${FILENAMES[$i]}" "${FILE_NEW}"
        if [[ -f "${FILENAMES[$i]}" ]]
        then
          echo "${FILENAMES[$i]} could not be moved, because ${FILE_NEW} already exists."
        fi
        break
#      else
#        echo "${SIG} does not match ${SIGNATURES[$i]}"
      fi
    done
  fi

done


