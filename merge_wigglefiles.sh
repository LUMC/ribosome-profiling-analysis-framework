#!/bin/bash
 ###############################################################################
 #
 # MergeWiggleFiles finds wiggle files in the current directory and merges any
 # wiggle files starting with the same part before the first dot (.) or
 # underscore (_) in the file name together, and creates a merged wiggle file.
 # It does so per strand (F and R), and separately for filtered and unfiltered
 # data. This information is taken from the file name.
 # NOTE: File names must contain "transcriptome_aligned" for the transcriptome
 # alignments, and "genome_aligned" for the genomic alignments.
 # The filtered transcriptome alignment files should end in
 # ".F.filtered.wig5" and ".R.filtered.wig5".
 # NOTE: The resulting Wiggle file, created by Wiggelen, uses spaces instead of
 # tabs as separator, regardless of what the original files are using. It is not
 # clear whether or not space or tab *should* be used; it seems that all
 # services accept both. No official standard seems to exist. Since we are using
 # tabs everywhere, because SAM tools are using them, we're replacing the spaces
 # with tabs.
 #
 # Created     : 2013-09-19
 # Modified    : 2015-03-26
 # Version     : 0.32
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

# Check: can not have argument.
if [[ ! -z $1 ]]
then
  echo "Usage: $0"
  echo "Do not use any arguments, MergeWiggleFiles will try to find"
  echo "  wiggle files in the current directory.";
  exit 1
fi

if [[ -z `ls | grep '.wig5$'` ]]
then
  echo "No wiggle files found in the current directory."
  exit 1;
fi

for sample in `ls *.wig5 | cut -d . -f 1 | cut -d _ -f 1 | sort | uniq`;
do
  for strand in F R;
  do
    if [[ -z `ls ${sample}*.wig5 | grep genome_aligned | grep ${strand}.wig5` ]]
    then
      echo "Could not find ${sample}*genome_aligned*.${strand}.wig5.";
      exit 1;
    elif [[ -z `ls ${sample}*.wig5 | grep transcriptome_aligned | grep ${strand}.wig5` ]]
    then
      echo "Could not find ${sample}*transcriptome_aligned*.${strand}.wig5.";
      exit 1;
    elif [[ -z `ls ${sample}*.wig5 | grep transcriptome_aligned | grep ${strand}.filtered.wig5` ]]
    then
      echo "Could not find ${sample}*transcriptome_aligned*.${strand}.filtered.wig5.";
      exit 1;
    fi

    rm -f "${sample}.merged_wiggle.${strand}.wig5";
    rm -f "${sample}.merged_wiggle.${strand}.filtered.wig5";
    wiggelen merge ${sample}*genome_aligned*.${strand}.wig5 ${sample}*transcriptome_aligned*.${strand}.wig5 | sed 's/^\([0-9]\+\) \([0-9]\+\)/\1\t\2/' > "${sample}.merged_wiggle.${strand}.wig5";
    echo -n ".";
    wiggelen merge ${sample}*genome_aligned*.${strand}.wig5 ${sample}*transcriptome_aligned*.${strand}.filtered.wig5 | sed 's/^\([0-9]\+\) \([0-9]\+\)/\1\t\2/' > "${sample}.merged_wiggle.${strand}.filtered.wig5";
    echo -n ".";
  done;
  # This process leaves .idx files laying around. Remove them by:
  rm ${sample}*.idx;
  echo -n " ";
done;
echo "";
