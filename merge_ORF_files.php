#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * Merge ORF Files merges the ORF analysis files created by the "find_ORFs.php"
 * script. It creates a summary file, where only the genomic locations, the gene
 * and the coverage per sample is shown.
 *
 * Created     : 2013-10-08
 * Modified    : 2013-11-22
 * Version     : 0.2
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.2',
        'output_suffix' => '.merged_ORF_analyses.txt',
    );

echo 'Merge ORF Files v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 1) { // I guess you could run it with just one file...
    die('Usage: ' . $sScriptName . ' ORF_FILE1 [ORF_FILE2 [ORF_FILE3 [...]]]' . "\n\n");
}

// Check if all files can be read.
$nSamples = count($aFiles);
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// How should we call the new file?
if ($nSamples == 1) {
    $sFileOut = $aFiles[0] . $_SETT['output_suffix'];
} else {
    // Find prefix and suffix for file(s), to have an output file that matches the name.
    $nPrefixLength = $nSuffixLength = min(array_map('strlen', $aFiles)); // Length of the shortest file name in argument list.
    $sPrefix = substr($aFiles[0], 0, $nPrefixLength); // Limit prefix already to length of shortest file name.
    $sSuffix = substr($aFiles[0], -$nSuffixLength); // Limit suffix already to length of shortest file name.
    foreach ($aFiles as $sFile) {
        for ($i = 0; $i < $nPrefixLength; $i++) {
            if ($sPrefix{$i} != $sFile{$i}) {
                // No match!
                $sPrefix = substr($sPrefix, 0, $i);
                $nPrefixLength = strlen($sPrefix);
                break; // Go to next file.
            }
        }
    }
    foreach ($aFiles as $sFile) {
        for ($i = 1; $i < $nSuffixLength; $i++) {
            if (substr($sSuffix, -$i, 1) != substr($sFile, -$i, 1)) {
                // No match!
                $sSuffix = substr($sSuffix, -($i-1));
                $nSuffixLength = strlen($sSuffix);
                break; // Go to next file.
            }
        }
    }
    $sFileOut = $sPrefix . preg_replace('/^' . preg_quote($sPrefix, '/') . '(.+)' . preg_quote($sSuffix, '/') . '$/', '$1', $aFiles[0]) . '-' . preg_replace('/^' . preg_quote($sPrefix, '/') . '(.+)' . preg_quote($sSuffix, '/') . '$/', '$1', $aFiles[$nSamples-1]) . $sSuffix . $_SETT['output_suffix'];
}





// Checking if we are allowed to create the output file.
if (file_exists($sFileOut)) {
    if (!is_writable($sFileOut)) {
        die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
    }
} elseif (!is_writable(dirname($sFileOut))) {
    die('Can not create ' . $sFileOut . ', aborting.' . "\n");
}





// Now, loop the ORF analysis files, load them one by one in the memory.
$aData = array(); // chr => array(position => array(gene, sample1 => coverage, sample2 => coverage));
$aSamples = array();
foreach ($aFiles as $sFile) {
    $sSampleID = substr($sFile, $nPrefixLength, -$nSuffixLength);
    $aSamples[] = $sSampleID;
    $aORFFile = file($sFile, FILE_IGNORE_NEW_LINES);
    print('Parsing ' . $sFile . '... ');
    $sGene = '';
    $nPositions = 0;
    foreach ($aORFFile as $sLine) {
        if (!trim($sLine)) {
            continue;
        }
        // We're looking at the data, or just before.
        // If we don't have a gene yet, look for it.
        if (preg_match('/^(.+)\tPositions found:\t\d+\tPositions analyzed:\t\d+\tTSS found:\t\d+$/', $sLine, $aRegs)) {
            $sGene = $aRegs[1];
        }
        if ($sGene && preg_match('/^(chr(\d+|[XYM])):(\d+)\t(\d+)(\t[0-9*-]+)+$/', $sLine, $aRegs)) {
            // We have a gene, and we matched a data line.
            $nPositions ++;
            list(,$sChr, , $nPosition, $nCoverage) = $aRegs;
            if (!isset($aData[$sChr])) {
                $aData[$sChr] = array();
            }
            if (!isset($aData[$sChr][$nPosition])) {
                $aData[$sChr][$nPosition] = array($sGene);
            }
            $aData[$sChr][$nPosition][$sSampleID] = $nCoverage;
        }
    }
    print('done, loaded ' . $nPositions . ' positions.' . "\n");
}





$fOut = @fopen($sFileOut, 'w');
if (!$fOut) {
    die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
}
// Let user know we're working here...
print('Writing output to ' . $sFileOut . '... ');
fputs($fOut, '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n" .
    '# NOTE: When this script reports a coverage of 0, it simply means the position' . "\n" .
    '#   was not recognized as a translation start site in that sample.' . "\n" .
    '#   The actual measured coverage in that sample may not at all be 0.' . "\n" .
    '# Chromosome' . "\t" . 'Position');
foreach ($aSamples as $sSampleID) {
    fputs($fOut, "\t" . $sSampleID);
}
fputs($fOut, "\t" . 'Gene' . "\n");



// Now, sort the list of chromosomes properly, so we can loop through it to print the results.
ksort($aData);
foreach ($aData as $sChr => $aPositions) {
    // Sort the positions, too.
    ksort($aPositions);
    foreach ($aPositions as $nPosition => $aPosition) {
        fputs($fOut, $sChr . "\t" . $nPosition);
        foreach ($aSamples as $sSampleID) {
            if (isset($aPosition[$sSampleID])) {
                fputs($fOut, "\t" . $aPosition[$sSampleID]);
            } else {
                fputs($fOut, "\t" . '0');
            }
        }
        // Finally, put the gene symbol.
        fputs($fOut, "\t" . $aPosition[0] . "\n");
    }
}
print('done.' . "\n\n");
?>
