#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * This script analyzes the seq_gene.md file downloaded from:
 * ftp://ftp.ncbi.nih.gov/genomes/M_musculus/mapview/seq_gene.md.gz
 * and generates a list of transcripts that have a very short distance between
 * the ATG and the next intron. All transcripts where the distance is larger
 * than 10 codons (30 bases), are discarded.
 *
 * Created     : 2013-08-22
 * Modified    : 2013-08-22
 * Version     : 0.1
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.1',
        'maximum_distance' => 30, // In bases, the maximum distance to still be reported.
        'output_suffix' => '.distance_intron_to_ATG.txt',
    );

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) != 1) {
    die('Usage: ' . $sScriptName . ' path_to_seq_gene.md' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// Checking if we are allowed to create the output file.
$sFileOut = $aFiles[0] . $_SETT['output_suffix'];
if (file_exists($sFileOut)) {
    if (!is_writable($sFileOut)) {
        die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
    }
} elseif (!is_writable(dirname($sFileOut))) {
    die('Can not create ' . $sFileOut . ', aborting.' . "\n");
}





// Start by filtering using grep and cut, data directly into a variable (it's an 9.5MB string).
$aData = explode("\n", `grep CDS $aFiles[0] | grep GRCm38 | grep -v '|' | cut -f 3-5,12,14`);
$aTranscripts = array();
foreach ($aData as $sLine) {
    if (!trim($sLine)) {
        continue;
    }
    list($nPosStart, $nPosEnd, $sStrand, $sType, $sTranscript) = explode("\t", $sLine);
    if ($sType != 'CDS' || ($sStrand == '+' && isset($aTranscripts[$sTranscript]))) {
        continue;
    }
    $nLength = ($nPosEnd-$nPosStart)+1-3;
    $aTranscripts[$sTranscript] = $nLength;
}
unset($aData);



// We need to filter AFTER the array has been built, since getting the first exon in the array requires all values to be stored.
// It's just simpler to sort first, then stop outputting once we hit the threshold.
asort($aTranscripts, SORT_NUMERIC);
// FIXME: Also sort on transcript?



$fOut = @fopen($sFileOut, 'w');
if (!$fOut) {
    die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
}
$i = 0;
foreach ($aTranscripts as $sTranscript => $nLength) {
    if ($nLength > $_SETT['maximum_distance']) {
        break;
    }
    fputs($fOut, $sTranscript . "\t" . $nLength . "\n");
    $i ++;
}
print('Done, ' . $i . ' transcripts found within threshold of ' . $_SETT['maximum_distance'] . '.' . "\n");
?>
