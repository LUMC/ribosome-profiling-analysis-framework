#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * CodonPositionAnalyzer analyzes the Mutalyzer results, compared to the
 * original wiggle files, checks the strand of the transcript in a third file,
 * and presents the results.
 *
 * Created     : 2013-04-16
 * Modified    : 2014-07-14
 * Version     : 0.4
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 * Changelog   : 0.4
 *               Fixed problem with parsing the Wiggle files; we were skipping
 *               the chrom=NC_* header, adding its positions to the last used
 *               chromosome (usually chrY).
 *
 *************/

$_SETT =
    array(
        'version' => '0.4',
        'min_coverage' => 3,
        'max_upstream' => 500,
        'max_downstream' => 500,
    );

echo 'CodonPositionAnalyzer v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) != 4) {
    die('Usage: ' . $sScriptName . ' MUTALYZER_FILE WIGGLE_FILE GENE_LIST_FILE STRAND' . "\n\n");
}
$sArgStrand = array_pop($aFiles);
if ($sArgStrand == 'F') {
    $sArgStrand = '+';
} elseif ($sArgStrand == 'R') {
    $sArgStrand = '-';
} elseif (!in_array($sArgStrand, array('+','-'))) {
    die('Strand argument is invalid. Please choose from +, -, F or R.' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

function sortPositions ($key1, $key2) {
    // Sorts the positions.
    if ($key1{0} == '*' && $key2{0} == '*') {
        return (substr($key1, 1) < substr($key2, 1)? -1 : 1); // I know the values will never be equal.
    } elseif ($key1{0} == '*' && $key2{0} != '*') {
        return 1;
    } elseif ($key1{0} != '*' && $key2{0} == '*') {
        return -1;
    } else {
        return ($key1 < $key2? -1 : 1); // I know the values will never be equal.
    }
}



// Prepare wiggle file, read into memory.
$aWiggleFile = file($aFiles[1], FILE_IGNORE_NEW_LINES);
$aCoverages = array();
$nChroms = 0;
$sChrom = '';
print('Parsing Wiggle file... ');
foreach ($aWiggleFile as $sLine) {
    if (preg_match('/^variableStep chrom=(.+)$/', $sLine, $aRegs)) {
        // Chromosome found.
//        if (!preg_match('/^chr([0-9]+|[XYM])$/', $aRegs[1])) { // FOR NOW, IGNORE chrM!!!
        if (!preg_match('/^chr([0-9]+|[XY])$/', $aRegs[1])) {
            $sChrom = '';
        } else {
            $sChrom = $aRegs[1];
            $nChroms ++;
        }
        continue;
    }
    if ($sChrom) {
        list($nPos, $nCoverage) = explode("\t", $sLine);
        $aCoverages[$sChrom . ':g.' . $nPos . 'del'] = $nCoverage; // Faking the variant.
    }
}
print('done, loaded ' . $nChroms . ' chromosomes with ' . count($aCoverages) . ' positions in memory.' . "\n");

// Prepare gene file.
$aTranscriptFile = file($aFiles[2], FILE_IGNORE_NEW_LINES);
unset($aTranscriptFile[0]); // Header.
$aTranscripts = array();
print('Parsing gene file... ');
foreach ($aTranscriptFile as $sLine) {
    if (!trim($sLine) || $sLine{0} == '#') {
        continue;
    }
    list($sTranscript, $sStrand, $sGene) = explode("\t", $sLine);
    $aTranscripts[$sTranscript] = $sStrand; // We'll ignore the gene name for now.
}
print('done, loaded ' . count($aTranscripts) . ' transcripts in memory.' . "\n");





// Now, loop the Mutalyzer file, find coverages in array (if not present, complain), check strand in gene list, and summarize results.
$nFiltered = 0;
$nIntronicPositions = 0;
$nUnmappable = 0;
$nIntergenic = 0;
$aTranscriptsPerPosition = array(); // Stores how many positions have how many transcripts, pure statistics.
$aCodonPositions = array();
$aUnknownTranscripts = array();
$aMutalyzerResults = file($aFiles[0]);
unset($aMutalyzerResults[0]); // Header.
print('Parsing mutalyzer results file... ');
foreach ($aMutalyzerResults as $sLine) {
    $aLine = explode("\t", rtrim($sLine)); // Removing whitespace from the right.
    $sVariant = array_shift($aLine);

    // Get coverage.
    if (!isset($aCoverages[$sVariant])) {
        die("\n" .
            'Cannot find coverage for position ' . $sVariant . ', probably you selected the wrong Wiggle file for this Mutalyzer file?' . "\n");
    }
    $nCoverage = $aCoverages[$sVariant];
    // Filter for low coverage.
    if ($aCoverages[$sVariant] < $_SETT['min_coverage']) {
        $nFiltered ++;
        continue;
    }

    // We need at least 4 values; input, errors, chrom_var, var_on_transcript.
    $nTranscripts = count($aLine) - 2;
    if ($nTranscripts < 0) {
        $nTranscripts = 0;
    }
    // Prevent notices...
    if (!isset($aTranscriptsPerPosition[$nTranscripts])) {
        $aTranscriptsPerPosition[$nTranscripts] = 0;
    }
    $aTranscriptsPerPosition[$nTranscripts] ++;
    if ($nTranscripts >= 1) {
        $sError = array_shift($aLine);
        if ($sError) {
            // If we have an error, would we ever have an array of at least 4? I think not...
            die("\n" .
                'Position ' . $sVariant . ' somehow generated an error: ' . $sError . "\n");
        }
        array_shift($aLine); // We're ignoring the mapping on the chromosome, which was our input anyways.
        // What is left is an array with at least one mapping to a transcript.
        $b5UTR = false;
        $bCoding = false; // 'Coding region' also applies to -15 to -1.
        $bIntronic = false;
        $b3UTR = false;
        // Store all options first, then decide what to do depending on the resulting options.
        $aCodonOptions = array();
        foreach ($aLine as $sVOT) {
            if (preg_match('/^([NX]R_\d+)/', $sVOT)) {
                // Non-coding RNA... ignore it, even though it probably doesn't do much.
                continue;
            } elseif (!preg_match('/^([NX][RM]_\d+)\.\d+:(.+)/', $sVOT, $aRegs)) {
                die("\n" .
                    'Cannot parse variant ' . $sVOT . "\n");
            }

            $sTranscript = $aRegs[1];
            $sPosition = $aRegs[2];
            // Check strand and then check if position is in coding region. If so, store position!
            if (!isset($aTranscripts[$sTranscript])) {
                // This happens quite a lot...
                if (!in_array($sTranscript, $aUnknownTranscripts)) {
                    $aUnknownTranscripts[] = $sTranscript;
                }
                continue; // On to the next VOT for this position.
            }

            if ($sArgStrand == $aTranscripts[$sTranscript]) {
                // Correct strand!
                if (preg_match('/^[cn]\.(\*)?(\-?\d+)del$/', $sPosition, $aRegs)) {
                    // Got one!
                    $aCodonOptions[] = $aRegs[1] . $aRegs[2];
                    if ($aRegs[1]) {
                        // 3'UTR!
                        $b3UTR = true;
                    } elseif ($aRegs[2] < -15) {
                        $b5UTR = true; // Other values >= -15 && < 0 are very valuable and not regarded 5'UTR.
                    } else {
                        $bCoding = true;
                    }
                } elseif (preg_match('/^[cn]\.(\*)?(\-?\d+)([+-]\d+)del$/', $sPosition)) {
                    // Intronic, we must have got the wrong transcript here!
                    $bIntronic = true;
                } else {
                    // Euh....
                    die($sPosition);
                }
            }
        }

        // Now, decide, based on the options, what to do.
        // If we have no options, disregard position.
        if (!count($aCodonOptions)) {
            // Has not been mapped to any transcript on the correct strand!
            if ($bIntronic) {
                // We did find intronic mappings...
                $nIntronicPositions ++;
            } else {
                // Nothing on this strand...
                $nUnmappable ++;
            }
            continue; // Continue to next position.
        }

        // Now filter distances, more than 500bp up- or downstream is too much.
        foreach ($aCodonOptions as $sTranscript => $sPosition) {
            if ($sPosition{0} == '-' && $sPosition < -$_SETT['max_upstream']) {
                unset($aCodonOptions[$sTranscript]);
            } elseif ($sPosition{0} == '*' && substr($sPosition, 1) > $_SETT['max_downstream']) {
                unset($aCodonOptions[$sTranscript]);
            }
        }
        // Check if we filtered it out completely now.
        if (!count($aCodonOptions)) {
            // Too far from any known transcript!
            $nIntergenic ++;
            continue; // Continue to next position.
        }

        // Check if the 5' and 3' flags are still correct.
        $b5UTR = false;
        $b3UTR = false;
        foreach ($aCodonOptions as $sPosition) {
            if ($sPosition{0} == '-' && $sPosition < -15) {
                $b5UTR = true;
            } elseif ($sPosition{0} == '*') {
                $b3UTR = true;
            }
        }

        // If we have 5'UTR positions, but also positions in the coding region, the 5'UTR is taken out.
        if ($b5UTR && $bCoding) {
            foreach ($aCodonOptions as $nKey => $sPosition) {
                if ($sPosition{0} == '-' && $sPosition < -15) {
                    unset($aCodonOptions[$nKey]);
                }
            }
            $b5UTR = false;
        }

        // If we have 3'UTR positions, but also positions in the coding region, the 3'UTR is taken out.
        if ($bCoding && $b3UTR) {
            foreach ($aCodonOptions as $nKey => $sPosition) {
                if ($sPosition{0} == '*') {
                    unset($aCodonOptions[$nKey]);
                }
            }
            $b3UTR = false;
        }

        // If we have 5' and 3'UTR positions, but no positions in the coding region, compare the 5' and 3' positions.
        if ($b5UTR && !$bCoding && $b3UTR) {
            $nMin = -$_SETT['max_upstream'];
            $nMax = $_SETT['max_downstream'];
            // Determine the 5' and 3' positions closest to the translation start and end sites.
            foreach ($aCodonOptions as $sPosition) {
                if ($sPosition{0} == '*') {
                    // 3'
                    $nMax = min($nMax, substr($sPosition, 1));
                } else {
                    // 5'
                    $nMin = max($nMin, $sPosition);
                }
            }
            // Calculate the difference between each other.
            if (abs($nMin) < $nMax && ($nMax + $nMin) > 100) {
                // The 5' position is closer to the ATG than the 3' is to the TGA.
                foreach ($aCodonOptions as $sTranscript => $sPosition) {
                    if ($sPosition{0} == '*') {
                        unset($aCodonOptions[$sTranscript]);
                    }
                }
                $b3UTR = false;
            } elseif (abs($nMin) > $nMax && ($nMax + $nMin) < -100) {
                // The 3' position is closer to the TGA than the 5' is to the ATG.
                foreach ($aCodonOptions as $sTranscript => $sPosition) {
                    if ($sPosition{0} == '-') {
                        unset($aCodonOptions[$sTranscript]);
                    }
                }
                $b5UTR = false;
            }
        }



        // If everything is normal, report it.
        if ((!$b5UTR && $bCoding && !$b3UTR) || ($b5UTR && !$bCoding && !$b3UTR) || (!$b5UTR && !$bCoding && $b3UTR)) {
            foreach ($aCodonOptions as $sPosition) {
                if ($sPosition{0} != '*' && $sPosition >= -15) {
                    // Group reads together, store position in codon.
                    if ($sPosition < 0) {
                        // 'Coding region', but -15 to -1.
                        $nCodonPosition = ($sPosition%3)+4; // -6->4, -5->2, -4->3.
                        if ($nCodonPosition == 4) {
                            // Second position.
                            $nCodonPosition = 1;
                        }
                    } else {
                        // Real coding positions.
                        $nCodonPosition = $sPosition % 3; // 4->1, 5->2, 6->0.
                        if (!$nCodonPosition) {
                            // Third position.
                            $nCodonPosition = 3;
                        }
                    }
                    // Preventing notices...
                    if (!isset($aCodonPositions['%' . $nCodonPosition])) {
                        $aCodonPositions['%' . $nCodonPosition] = 0;
                    }
                    // Store value!
                    $aCodonPositions['%' . $nCodonPosition] += $aCoverages[$sVariant];

                    if ($sPosition > 6) {
                        // Not individually counted, goodbye...
                        continue; // To next codon position.
                    }
                }
                // Preventing notices...
                if (!isset($aCodonPositions[$sPosition])) {
                    $aCodonPositions[$sPosition] = 0;
                }
                $aCodonPositions[$sPosition] += $aCoverages[$sVariant];
            }
        } else {
            // Currently unhandled situation.
            print('Not implemented:' . "\n" .
                  '5UTR: ' . (int) $b5UTR . ', Coding: ' . (int) $bCoding . ', 3UTR: ' . (int) $b3UTR . ', Location: ' . $sVariant . "\n");
        }

    } elseif (count($aLine) == 1) {
        // Generated an error!
        $sError = $aLine[0];
        die('Position ' . $sVariant . ' somehow generated an error: ' . $sError . "\n");
    } else {
        // Unmappable...
        $nUnmappable ++;
    }
}
print('done!' . "\n");
sort($aUnknownTranscripts);
print('Total positions:' . "\t" . count($aMutalyzerResults) . "\n" .
      'Filtered for low coverage:' . "\t" . $nFiltered . "\n" .
      'Transcripts not found in gene file:' . "\t" . count($aUnknownTranscripts) . "\t" . implode(';', $aUnknownTranscripts) . "\n" .
      'Transcripts per position:' . "\n" .
      'Transcripts' . "\t" . 'Count' . "\n");
ksort($aTranscriptsPerPosition);
foreach($aTranscriptsPerPosition as $nTranscripts => $nPositions) {
    print($nTranscripts . "\t" . $nPositions . "\n");
}

print("\n" .
      'Positions left after filtering:' . "\t" . array_sum($aTranscriptsPerPosition) . "\n" .
      'Positions not mappable:' . "\t" . $nUnmappable . "\t" . 'Possible causes: no mapping by Mutalyzer, transcript is missing, or strand is wrong' . "\n" .
      'Positions intronic only:' . "\t" . $nIntronicPositions . "\t" . 'Possible causes: transcript is missing, or strand is wrong' . "\n" .
      'Positions too far from known genes:' . "\t" . $nIntergenic . "\t" . 'Possible causes: transcript is missing, newly discovered transcript, or strand is wrong' . "\n" .
      'Mapping to codon positions:' . "\n" .
      'Position' . "\t" . 'Total coverage' . "\n");
uksort($aCodonPositions, 'sortPositions');
foreach($aCodonPositions as $nPosition => $nCoverage) {
    print($nPosition . "\t" . $nCoverage . ($nPosition{0} != '%'? '' : "\t" . '(including -15 to -1)') . "\n");
}
print("\n");
exit(0);
?>
