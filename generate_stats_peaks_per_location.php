#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * Generates statistics; the number of peaks per location in a gene (5' UTR,
 * AUG, Downstream coding, 3' UTR, Multiple. It takes all analysis results files
 * from the find_ORFs.php script, and generates one result file per sample.
 *
 * Created     : 2014-01-08
 * Modified    : 2014-01-08
 * Version     : 0.1
 *
 * Copyright   : 2014 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.1',
        'output_suffix' => '.ORF_analysis_results.stats_peaks_per_location.txt',
        'categories' =>
        array(
            '5UTR',
            'AUG',
            'coding',
            '3UTR',
            'multiple',
        ),
    );

echo 'Stats: Peaks Per Location v.' . $_SETT['version'] . "\n";

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

// Go through input files, recognize types. Determine how we should call the new file.
$aSamples = array(); // Will store the info about the samples and which files we have on them.
if ($nSamples == 1) {
    $sFileOut = $aFiles[0] . $_SETT['output_suffix'];
} else {
    foreach ($aFiles as $nFile => $sFile) {
        // Toss possible statistics files or result files out.
        if (preg_match('/\.ORF_analysis_results[\._]stats(_peaks_per_location)?\./', $sFile)) {
            unset($aFiles[$nFile]);
            continue;
        }
        // Rest is matched.
        if (!preg_match('/^(.+)\.(F|R)(?:\..*)?\.ORF_analysis_results(_after_cutoff)?\.txt$/', $sFile, $aRegs)) {
            //              [1]   [2]                                 [3]
            // Unrecognized file, complain.
            die('Sorry, I do not understand the file name of the file \'' . $sFile . '\', aborting.' . "\n");
        }
        $sSampleID = $aRegs[1]; // Actually, basically the whole prefix until the strand info.
        $sStrand = $aRegs[2];
        $bCutOff = !empty($aRegs[3]);
        if (!isset($aSamples[$sSampleID])) {
            $aSamples[$sSampleID] =
                array(
                    'F' => array(), // We'll end up with two keys here hopefully: false = file with peaks before cutoff, true = file with peaks after cutoff.
                    'R' => array(),
                    'data' => array(
                        false => array_combine($_SETT['categories'], array_fill(0, count($_SETT['categories']), 0)),
                        true => array_combine($_SETT['categories'], array_fill(0, count($_SETT['categories']), 0)),
                    ),
                );
        }
        $aSamples[$sSampleID][$sStrand][$bCutOff] = $sFile;
    }
}





// Checking if we are allowed to create the output files.
foreach ($aSamples as $sSampleID => $aSample) {
    $sFileOut = $sSampleID . $_SETT['output_suffix'];
    $aSamples[$sSampleID]['file_out'] = $sFileOut;
    if (file_exists($sFileOut)) {
        if (!is_writable($sFileOut)) {
            die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
        }
    } elseif (!is_writable(dirname($sFileOut))) {
        die('Can not create ' . $sFileOut . ', aborting.' . "\n");
    }

    // Nicely sort the files, so we always parse them in the same order (before cutoff, after cutoff).
    ksort($aSamples[$sSampleID]['F']);
    ksort($aSamples[$sSampleID]['R']);
}





// Now, loop the ORF analysis files, load them one by one in the memory.
foreach ($aSamples as $sSampleID => $aSample) {
    foreach (array('F', 'R') as $sStrand) {
        $aFiles = $aSample[$sStrand];
        foreach ($aFiles as $bCutOff => $sFile) {
            $aORFFile = file($sFile, FILE_IGNORE_NEW_LINES);
            print('Parsing ' . $sFile . '... ');
            $nPositions = 0;
            foreach ($aORFFile as $sLine) {
                if (!trim($sLine)) {
                    continue;
                }
                if (preg_match('/^(chr(?:\d+|[XYM])):(\d+)\t(\d+)((?:\t[0-9*-]+)+)$/', $sLine, $aRegs)) {
                    //             [1]               [2]    [3]   [4]
                    // We have matched a data line.
                    $nPositions ++;
                    list(,$sChr, $nPosition, $nCoverage) = $aRegs;
                    $aPositions = array();
                    $aPositionsInGene = explode("\t", trim($aRegs[4]));
                    // Loop all positions in genes, and determine category. Store these categories.
                    foreach ($aPositionsInGene as $sPositionInGene) {
                        if ($sPositionInGene == '-') {
                            // There was no mapping on this transcript.
                            continue;
                        }
                        if ($sPositionInGene{0} == '-') {
                            if ($sPositionInGene < -12) {
                                $sCategory = '5UTR';
                            } elseif ($sPositionInGene >= -12 && $sPositionInGene <= -10) {
                                $sCategory = 'AUG';
                            } else {
                                $sCategory = 'coding';
                            }
                        } elseif ($sPositionInGene{0} == '*') {
                            $sCategory = '3UTR';
                        } else {
                            $sCategory = 'coding';
                        }
                        $aPositions[] = $sCategory;
                    }
                    // Now see how many different categories we have for this position.
                    $aPositions = array_unique($aPositions);
                    if (count($aPositions) == 1) {
                        // One transcript, or multiple but at least in the same category of position.
                        $sCategory = $aPositions[0];
                    } else {
                        // Multiple different positions. If one of them is AUG, we will assume AUG.
                        // Otherwise, we don't know what to do, and we call this 'multiple';
                        if (in_array('AUG', $aPositions)) {
                            $sCategory = 'AUG';
                        } else {
                            $sCategory = 'multiple';
                        }
                    }

                    // We have now determined the category. Store, and count.
                    // FIXME: We can also just inject directly into $aSample, so we don't need to reload $aSample later.
                    // Anyways we don't need this information outside of this loop.
                    $aSamples[$sSampleID]['data'][$bCutOff][$sCategory] ++;
                }
            }
            print('done, loaded ' . $nPositions . ' positions.' . "\n");
        }
    }
    $aSample = $aSamples[$sSampleID]; // Reload, since we're in a foreach and we're working on a copy of the array.





    $sFileOut = $aSample['file_out'];
    $fOut = @fopen($sFileOut, 'w');
    if (!$fOut) {
        die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
    }
    // Let user know we're working here...
    print('Writing output to ' . $sFileOut . '... ');
    fputs($fOut, '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n" .
        '# NOTE: Read start positions at the end of the coding region, less than 12 bp' . "\n" .
        '# away from the 3\'UTR, are counted as coding while in fact in reality they are' . "\n" .
        '# of course a peak in the 3\'UTR. This can not be detected however, because we' . "\n" .
        '# don\'t know the length of the coding region of the transcripts.' . "\n");

    foreach (array(false, true) as $bCutOff) {
        fputs($fOut, "\n\n" .
            '# Results for ORF start sites with ' . ($bCutOff? 'no' : 'a') . ' cutoff applied. Using files:' . "\n" .
            '# ' . $aSample['F'][$bCutOff] . "\n" .
            '# ' . $aSample['R'][$bCutOff] . "\n" .
            (!$bCutOff? '' :
                '# ' . $aSample['F'][!$bCutOff] . "\n" .
                '# ' . $aSample['R'][!$bCutOff] . "\n") .
            '# Category' . "\t" . 'Number of TSSs found' . "\n");

        foreach ($_SETT['categories'] as $sCategory) {
            fputs($fOut, $sCategory . "\t" . ($aSample['data'][$bCutOff][$sCategory] + (!$bCutOff? 0 : $aSample['data'][!$bCutOff][$sCategory])) . "\n");
        }
    }
    print('done.' . "\n\n");
}
?>
