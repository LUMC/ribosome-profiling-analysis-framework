#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * Merge ORF Files merges the ORF analysis files created by the "find_ORFs.php"
 * script. It creates a summary file, where only the genomic locations, the gene
 * and the coverage per sample is shown.
 *
 * Created     : 2013-10-08
 * Modified    : 2013-12-12
 * Version     : 0.3
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.3',
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
    $lPrefix = $lSuffix = min(array_map('strlen', $aFiles)); // Length of the shortest file name in argument list.
    $sPrefix = substr($aFiles[0], 0, $lPrefix); // Limit prefix already to length of shortest file name.
    $sSuffix = substr($aFiles[0], -$lSuffix); // Limit suffix already to length of shortest file name.
    foreach ($aFiles as $sFile) {
        for ($i = 0; $i < $lPrefix; $i++) {
            if ($sPrefix{$i} != $sFile{$i}) {
                // No match!
                $sPrefix = substr($sPrefix, 0, $i);
                $lPrefix = strlen($sPrefix);
                break; // Go to next file.
            }
        }
    }
    foreach ($aFiles as $sFile) {
        for ($i = 1; $i < $lSuffix; $i++) {
            if (substr($sSuffix, -$i, 1) != substr($sFile, -$i, 1)) {
                // No match!
                $sSuffix = substr($sSuffix, -($i-1));
                $lSuffix = strlen($sSuffix);
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
    $sSampleID = substr($sFile, $lPrefix, -$lSuffix);
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

// There are two modes in which we can merge the files. In case the actual replicates are given, we join normally.
// But in case the replicates were merged before, we need to use the peaks given by the input files, but read the
// individual Wiggle files of the replicates and show all sample's coverage's. The note about the peak calling
// (that 0 does not necessarily mean no coverage) is then no longer correct, and is hidden.
// We "detect" that the replicates have been merged, by looking at the length of the sample IDs. If those are 1,
// like "A" or "C", then we assume they are merged. Otherwise, with longer sample IDs (like A1, A2, C1, C2, etc),
// We assume they are not merged. We could also use a regular expression pattern, if we would be more flexible.
$bMergedReplicates = (strlen($aSamples[0]) == 1); // Sample IDs should be of the same length.

// If we've been merging, try and find the original Wiggle files.
if ($bMergedReplicates) {
    $aReplicates = array();
    foreach ($aFiles as $key => $sFile) {
        // Try and find the replicates.
        for ($i = 1; $i < 10; $i ++) {
            $sReplicateID = $aSamples[$key] . $i;
            $sReplicateFile = substr($sFile, 0, $lPrefix) . $sReplicateID . substr($sFile, -$lSuffix);
            // Remove .ORF_analysis_results.txt suffix.
            $sReplicateFile = preg_replace('/\.ORF_analysis_results(_after_cutoff)?\.txt$/', '', $sReplicateFile);
            if (preg_match('/wig5?$/', $sReplicateFile) && is_readable($sReplicateFile)) {
                $aReplicates[$sReplicateID] = $sReplicateFile;
            } else {
                break;
            }
        }
    }

    // Check if we found replicate Wiggle files, parse these, fill in coverages in $aData.
    if ($aReplicates) {
        foreach ($aReplicates as $sReplicateID => $sReplicateFile) {
            $aWiggleFile = file($sReplicateFile, FILE_IGNORE_NEW_LINES);
            $sChrom = '';
            foreach ($aWiggleFile as $sLine) {
                if (preg_match('/^variableStep chrom=(chr.+)$/', $sLine, $aRegs)) {
                    // Chromosome found.
                    if (!preg_match('/^chr([0-9]+|[XYM])$/', $aRegs[1])) {
                        $sChrom = '';
                    } else {
                        $sChrom = $aRegs[1];
                    }
                    continue;
                }
                if ($sChrom) {
                    // Check if this position is considered a peak.
                    list($nPos, $nCoverage) = explode("\t", $sLine);
                    $nPos = (int)$nPos;
                    if (isset($aData[$sChrom][$nPos])) {
                        $aData[$sChrom][$nPos][$sReplicateID] = $nCoverage;
                    }
                }
            }
        }

        // Overwrite $aSamples.
        $aSamples = array_keys($aReplicates);

    } else {
        // No replicates found, we just act as if we are not looking at merged files.
        $bMergedReplicates = false;
    }
}
fputs($fOut, '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n");

if ($bMergedReplicates) {
    fputs($fOut, '# Found and parsed Wiggle file: ' . implode("\n" . '# Found and parsed Wiggle file: ', $aReplicates) . "\n");
} else {
    fputs($fOut,
        '# NOTE: When this script reports a coverage of 0, it simply means the position' . "\n" .
        '#   was not recognized as a translation start site in that sample.' . "\n" .
        '#   The actual measured coverage in that sample may not at all be 0.' . "\n");
}

fputs($fOut, '# Chromosome' . "\t" . 'Position');
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
