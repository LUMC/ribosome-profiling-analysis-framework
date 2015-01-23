#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * WIG2BATCHFILE converts wiggle files to the Mutalyzer batch format, such that
 * the position converter van be used.
 *
 * Created     : 2013-04-10
 * Modified    : 2014-07-10
 * Version     : 0.3
 *
 * Copyright   : 2013-2015 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 * Changelog   : 0.2
 *               While generating the batchfile, the wiggle file is already
 *               filtered for low coverage. This results in faster conversion,
 *               smaller files, and also less calculation time for Mutalyzer.
 *               0.3
 *               Fixed problem with skipping the chrom=NC_* header, adding its
 *               positions to the last used chromosome (usually chrY).
 *
 *
 * This work is licensed under the Creative Commons
 * Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
 * copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
 * or send a letter to:
 * Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
 *
 *************/

$_SETT =
    array(
        'version' => '0.3',
        'min_coverage' => 3,
    );

echo 'WIG2BATCHFILE v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 1) {
    die('Usage: ' . $sScriptName . ' WIG_FILE1 [WIG_FILE2 [WIG_FILE3 [...]]]' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// Find prefix for file(s), to have an output file that matches the name.
$nPrefixLength = min(array_map('strlen', $aFiles)); // Length of the shortest file name in argument list.
$sPrefix = substr($aFiles[0], 0, $nPrefixLength); // Limit prefix already to length of shortest file name.
foreach ($aFiles as $sFile) {
    for ($i = 0; $i < $nPrefixLength; $i ++) {
        if ($sPrefix{$i} != $sFile{$i}) {
            // No match!
            $sPrefix = substr($sPrefix, 0, $i);
            $nPrefixLength = strlen($sPrefix);
            break; // Go to next file.
        }
    }
}
$sPrefix = rtrim($sPrefix, '._-');
$sFileNameOut = $sPrefix . (!$sPrefix? '' : '_') . 'mutalyzer_batchfile.txt';

$fOut = @fopen($sFileNameOut, 'w');
if (!$fOut) {
    die('Unable to open file for writing.' . "\n\n");
}

$nLines = $nFiltered = 0;
foreach ($aFiles as $sFile) {
    $aFile = file($sFile);
    $sChrom = '';
    foreach ($aFile as $sLine) {
        if (preg_match('/^variableStep chrom=(.+)$/', $sLine, $aRegs)) {
            // Chromosome found.
//            if (!preg_match('/^chr([0-9]+|[XYM])$/', $aRegs[1])) { // FOR NOW, IGNORE chrM!!!
            if (!preg_match('/^chr([0-9]+|[XY])$/', $aRegs[1])) {
                echo 'Unrecognized chromosome: ' . $aRegs[1] . "\n";
                $sChrom = '';
            } else {
                $sChrom = $aRegs[1];
            }
            continue;
        }
        if ($sChrom) {
            list($nPos, $nCoverage) = explode("\t", $sLine);
            if ($nCoverage >= $_SETT['min_coverage']) {
                fputs($fOut, $sChrom . ':g.' . $nPos . 'del' . "\r\n");
                $nLines ++;
            } else {
                $nFiltered ++;
            }
        }
    }
}
die('Done, ' . $nLines . ' lines written (' . $nFiltered . ' filtered, coverage too low).' . "\n");
?>
