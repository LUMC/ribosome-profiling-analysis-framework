#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * PackSamFiles converts sam files into packed sam files, saving space and
 * grouping the reads together into transcript, position, coverage.
 *
 * Created     : 2013-08-13
 * Modified    : 2013-08-15
 * Version     : 0.1
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.1',
        'suffix' => '.packed',
        'terminal_width' => 150,
    );

echo 'PackSamFiles v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 1) {
    die('Usage: ' . $sScriptName . ' SAM_FILE1 [SAM_FILE2 [SAM_FILE3 [...]]]' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}





foreach ($aFiles as $sFile) {
    $nLine = 0;
    // SAM files are BIG. Very, very, very, BIG. Really. HUGE, actually. GBs. We need to read this line by line, to prevent a memory error.
    $fIn = @fopen($sFile, 'r');
    if (!$fIn) {
        die('Unable to open ' . $sFile . '.' . "\n\n");
    }
    $sFileSize = filesize($sFile);
    $sBytesRead = 0;
    $sFileOut = $sFile . $_SETT['suffix'];
    $fOut = @fopen($sFileOut, 'w');
    if (!$fOut) {
        die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
    }

    $aData = array(); // Will contain transcripts as keys, with an array (position => coverage) as value.
    while ($sLine = fgets($fIn)) {
        $nLine ++;
        $sBytesRead += strlen($sLine);
        $sLine = rtrim($sLine);
        list(, $sReference, $nPosition) = explode("\t", $sLine);
        list(, , , $sTranscript,) = explode('|', $sReference);
        if ($sTranscript && $nPosition) {
            if (!isset($aData[$sTranscript])) {
                $aData[$sTranscript] = array();
            }
            if (!isset($aData[$sTranscript][$nPosition])) {
                $aData[$sTranscript][$nPosition] = 0;
            }
            $aData[$sTranscript][$nPosition] ++;
        } else {
            die('Can\'t parse line ' . $nLine . ' in file ' . $sFile . '.' . "\n\n");
        }

        if (!($nLine % 50000)) {
            $nPercentageRead = round($sBytesRead/$sFileSize, 2);
            $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
            $lDone = round($nPercentageRead*$nAvailableWidth);
            print(str_repeat(chr(8), $_SETT['terminal_width']) .
                '[' . str_repeat('=', $lDone) . str_repeat(' ', $nAvailableWidth - $lDone) . '] ' . $nLine . ' ' . str_pad(round($nPercentageRead*100), 3, ' ', STR_PAD_LEFT) . '%');
        }
    }
    fclose($fIn);

    print("\n" .
          'Done reading ' . $nLine . ' lines, writing output... ');

    // Sorting the results would be nice, but perhaps really unneccesary.
    ksort($aData);
    foreach ($aData as $sTranscript => $aTranscript) {
        ksort($aTranscript);
        foreach ($aTranscript as $nPosition => $nCoverage) {
            fputs($fOut, $sTranscript . "\t" . $nPosition . "\t" . $nCoverage . "\n");
        }
    }
    fclose($fOut);
    print('Done, wrote ' . count($aData) . ' lines.' . "\n");
}
die('All files done.' . "\n");
?>
