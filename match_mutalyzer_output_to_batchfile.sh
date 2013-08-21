#!/bin/bash
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


