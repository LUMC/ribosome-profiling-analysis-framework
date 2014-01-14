#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * ORF-Finder analyzes the Mutalyzer results, compared to the original wiggle
 * files, checks the strand of the transcript in a third file, and finds ORFs
 * that have not been annotated before.
 *
 * Created     : 2013-07-12
 * Modified    : 2014-01-14
 * Version     : 0.6 // 0.4 was version 2013-10-10, archived 2013-11-15.
 *
 * Copyright   : 2013-2014 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.6',
        'min_coverage' => 3,      // Positions with less coverage than this are ignored. NOTE: The Mutalyzer batch file has already been filtered for coverage lower than 3.
        'max_upstream' => 500,    // Maximum distance from known CDS look for ORFs.
        'max_downstream' => 500,  // Maximum distance from known CDS look for ORFs.
        'peak_finding' =>
        array(
            'min_coverage' => 10, // Candidate peaks need to have at least a coverage of 10.
            'codon_1st_pos_min_coverage_fraction' => 0.6, // The first position of a codon must have at least 60% of all coverage in that codon.
            'upstream_codons_to_check'   => 5, // How many codons upstream do we check to verify a candidate peak? (coverage should be higher than the max coverage in those codons)
            'downstream_codons_to_check' => 5, // How many codons downstream do we check to verify a candidate peak? (coverage should be higher than the max coverage in those codons)
            'downstream_codons_max_coverage_max_fraction' => 0.1,  // Downstream codons should not have a maximum coverage higher than 10% of the candidate peak.
            'candidate_coverage_min_fraction'             => 0.1,  // Each candidate peak must have at least 10% coverage relative to the highest candidate.
            'candidate_min_distance_reported_separately'  => 5000, // At what distance from the pORF start site should we report candidate peaks separately?
        ),
        'output_suffix' =>
        array(
            'stats'                => '.ORF_analysis_results_stats.txt',
            'results'              => '.ORF_analysis_results.txt',
            'results_after_cutoff' => '.ORF_analysis_results_after_cutoff.txt',
        ),
    );

echo 'ORF-Finder v.' . $_SETT['version'] . "\n";

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

// Checking if we are allowed to create the output files.
$aFilesOut = array();
foreach ($_SETT['output_suffix'] as $sType => $sSuffix) {
    $sFileOut = $aFiles[1] . $sSuffix;

    if (file_exists($sFileOut)) {
        if (!is_writable($sFileOut)) {
            die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
        }
    } elseif (!is_writable(dirname($sFileOut))) {
        die('Can not create ' . $sFileOut . ', aborting.' . "\n");
    }

    $fOut = @fopen($sFileOut, 'w');
    if (!$fOut) {
        die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
    }

    $aFilesOut[$sType] = array('name' => $sFileOut, 'handler' => $fOut);
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
$aCoverages = array(); // chr => array(position => coverage, ...)
$nChroms = 0;
$sChrom = '';
print('Parsing Wiggle file... ');
foreach ($aWiggleFile as $sLine) {
    if (preg_match('/^variableStep chrom=(chr.+)$/', $sLine, $aRegs)) {
        // Chromosome found.
//        if (!preg_match('/^chr([0-9]+|[XYM])$/', $aRegs[1])) { // FOR NOW, IGNORE chrM!!!
        if (!preg_match('/^chr([0-9]+|[XY])$/', $aRegs[1])) {
            $sChrom = '';
        } else {
            $sChrom = $aRegs[1];
            $nChroms ++;
            $aCoverages[$sChrom] = array();
        }
        continue;
    }
    if ($sChrom) {
        list($nPos, $nCoverage) = explode("\t", $sLine);
        $aCoverages[$sChrom][$nPos] = (int) $nCoverage;
    }
}
print('done.' . "\n");
fputs($aFilesOut['stats']['handler'], 'Loaded ' . $nChroms . ' chromosomes with ' . (count($aCoverages, true) - $nChroms) . ' positions in memory.' . "\n");

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
    $aTranscripts[$sTranscript] = array($sGene, $sStrand);
}
print('done.' . "\n");
fputs($aFilesOut['stats']['handler'], 'Loaded ' . count($aTranscripts) . ' transcripts in memory.' . "\n");





// Now, loop the Mutalyzer file, find coverages in array (if not present, complain), check strand in gene list, and summarize results.
$nFiltered = 0;
$nIntronicPositions = 0;
$nUnmappable = 0;
$nIntergenic = 0;
$aPositionsPerGene = array(); // Will contain genes as keys, with an array of positions, which is an array of positions per transcript and their mappings to the genome.
$aTranscriptsPerPosition = array(); // Stores how many positions have how many transcripts, pure statistics.
$aCodonPositions = array();
$aUnknownTranscripts = array();
$aMutalyzerResults = file($aFiles[0]);
unset($aMutalyzerResults[0]); // Header.
print('Parsing mutalyzer results file... ');
foreach ($aMutalyzerResults as $sLine) {
    $aLine = explode("\t", rtrim($sLine)); // Removing whitespace from the right.
    $sVariant = array_shift($aLine);
    list($sChr, $nPosition) = explode(';', preg_replace('/^(chr(?:[0-9]+|[XYM])):g\.(\d+)del$/', "$1;$2", $sVariant));

    // Get coverage.
    if (!isset($aCoverages[$sChr][$nPosition])) {
        die("\n" .
            'Cannot find coverage for position ' . $sVariant . ', probably you selected the wrong Wiggle file for this Mutalyzer file?' . "\n");
    }
    // Filter for low coverage.
    if ($aCoverages[$sChr][$nPosition] < $_SETT['min_coverage']) {
        // NOTE: The Mutalyzer batch file has already been filtered for coverage lower than 3.
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
            // Check if we have info on the transcript, then check strand and store positions!
            if (!isset($aTranscripts[$sTranscript])) {
                // This happens quite a lot...
                if (!in_array($sTranscript, $aUnknownTranscripts)) {
                    $aUnknownTranscripts[] = $sTranscript;
                }
                continue; // On to the next VOT for this position.
            }

            if ($sArgStrand == $aTranscripts[$sTranscript][1]) {
                // Correct strand!
                if (preg_match('/^[cn]\.(\*)?(\-?\d+)del$/', $sPosition, $aRegs)) {
                    // Got one!
                    $aCodonOptions[$sTranscript] = $aRegs[1] . $aRegs[2];
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

        // If we have 3'UTR positions, but also 5' UTR or coding region positions, remove the 3' positions.
        if ($b3UTR && ($b5UTR || $bCoding)) {
            foreach ($aCodonOptions as $sTranscript => $sPosition) {
                if ($sPosition{0} == '*') {
                    unset($aCodonOptions[$sTranscript]);
                }
            }
            $b3UTR = false;
        }



        // Now that we're done filtering, save all positions per gene, so we can loop through it and try and find patterns.
        foreach ($aCodonOptions as $sTranscript => $sPosition) {
            list($sGene, $sStrand) = $aTranscripts[$sTranscript];
            // Create gene array if it doesn't exist.
            if (!isset($aPositionsPerGene[$sGene])) {
                $aPositionsPerGene[$sGene] = array('chr' => $sChr, 'strand' => $sStrand, 'positions' => array(), 'unique_positions' => array(), 'unique_positions_analyzed' => array());
            }
            // If we don't know this transcript yet, add it to the list.
            if (!isset($aPositionsPerGene[$sGene]['positions'][$sTranscript])) {
                $aPositionsPerGene[$sGene]['positions'][$sTranscript] = array();
            }
            // Store position on transcript.
            $aPositionsPerGene[$sGene]['positions'][$sTranscript][$sPosition] = $nPosition;
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
$nGenesMapped = count($aPositionsPerGene);
print('done.' . "\n");
fputs($aFilesOut['stats']['handler'], 'Loaded ' . $nGenesMapped . ' genes with mappings.' . "\n");

sort($aUnknownTranscripts);
fputs($aFilesOut['stats']['handler'], 'Total positions:' . "\t" . count($aMutalyzerResults) . "\n" .
    'Filtered for low coverage:' . "\t" . $nFiltered . "\n" .
    'Transcripts not found in gene file:' . "\t" . count($aUnknownTranscripts) . "\t" . implode(';', $aUnknownTranscripts) . "\n" .
    'Transcripts per position:' . "\n" .
    'Transcripts' . "\t" . 'Count' . "\n");
ksort($aTranscriptsPerPosition);
foreach($aTranscriptsPerPosition as $nTranscripts => $nPositions) {
    fputs($aFilesOut['stats']['handler'], $nTranscripts . "\t" . $nPositions . "\n");
}

fputs($aFilesOut['stats']['handler'], "\n" .
    'Positions left after filtering:' . "\t" . array_sum($aTranscriptsPerPosition) . "\n" .
    'Positions not mappable:' . "\t" . $nUnmappable . "\t" . 'Possible causes: no mapping by Mutalyzer, transcript is missing, or strand is wrong' . "\n" .
    'Positions intronic only:' . "\t" . $nIntronicPositions . "\t" . 'Possible causes: transcript is missing, or strand is wrong' . "\n" .
    'Positions too far from known genes:' . "\t" . $nIntergenic . "\t" . 'Possible causes: transcript is missing, newly discovered transcript, or strand is wrong' . "\n" .
    'Positions left:' . "\t" . (array_sum($aTranscriptsPerPosition) - $nUnmappable - $nIntronicPositions - $nIntergenic) . "\n" .
    'Genes mapped to:' . "\t" . $nGenesMapped . "\n");
print('Now looking for peaks ');





// Now we'll look for candidate peaks. We'll start walking through the transcript, starting from 5' to 3', looking for peaks.
$nGenesDiscarded = 0;
$i = 0;
$aResults = array(); // Here we will store the list of found peaks per gene.
//$aPositionsPerGene = array('Csrp1' => $aPositionsPerGene['Csrp1']);
foreach ($aPositionsPerGene as $sGene => $aGene) {
    if (!(++$i%25)) {
        print('.');
    }

    foreach ($aGene['positions'] as $sTranscript => $aPositions) {
        if ($aGene['strand'] == '-' && count($aPositions) > 1) {
            uksort($aPositions, 'sortPositions');
        }
        // Store positions, for statistics.
        $aPositionsPerGene[$sGene]['unique_positions'] = array_merge($aPositionsPerGene[$sGene]['unique_positions'], array_values($aPositions));

        foreach ($aPositions as $sPosition => $nPosition) {
            if ($aCoverages[$aGene['chr']][$nPosition] < $_SETT['peak_finding']['min_coverage']) {
                // Coverage of this peak is not high enough for analysis.
                continue;
            }
            $b3UTR = (substr($sPosition, 0, 1) == '*'); // There should not be a lot of these.
            $b5UTR = (substr($sPosition, 0, 1) == '-');

            // Store position, for statistics.
            $aPositionsPerGene[$sGene]['unique_positions_analyzed'][] = $nPosition;

            // Step 1: Its coverage should be higher than that of the positions
            // located 3, 6, 9, 12 and 15 nucleotides upstream of it.
            for ($j = 1; $j <= $_SETT['peak_finding']['upstream_codons_to_check']; $j++) {
                if ($b3UTR) {
                    // FIXME: See note in README about this not working when crossing the stop codon position.
                    $sPositionToCheck = '*' . (substr($sPosition, 1) - (3*$j));
                } else {
                    $sPositionToCheck = $sPosition - (3*$j);
                    if (!$b5UTR && $sPositionToCheck <= 0) {
                        // We crossed the -1/1 border.
                        $sPositionToCheck--;
                    }
                }
                $nCoverage = (!isset($aPositions[$sPositionToCheck])? 0 : $aCoverages[$aGene['chr']][$aPositions[$sPositionToCheck]]);
                if ($aCoverages[$aGene['chr']][$nPosition] <= $nCoverage) {
                    // Coverage is NOT higher than coverage of previous position.
                    continue 2;
                }
            }

            // Step 2a: It should show a triplet periodicity...
            $nCoverageOfCodon = $aCoverages[$aGene['chr']][$nPosition];
            // The two other nucleotides in the possible codon, so the two nucleotides downstream of the position
            // we're analyzing, should not have an average coverage higher than 20% of the position's coverage.
            for ($j = 1; $j <= 2; $j++) {
                if ($b3UTR) {
                    $sPositionToCheck = '*' . ((int) substr($sPosition, 1) + $j);
                } else {
                    $sPositionToCheck = $sPosition + $j;
                    if ($b5UTR && $sPositionToCheck >= 0) {
                        // We crossed the -1/1 border.
                        $sPositionToCheck++;
                    }
                }
                $nCoverage = (!isset($aPositions[$sPositionToCheck])? 0 : $aCoverages[$aGene['chr']][$aPositions[$sPositionToCheck]]);
                $nCoverageOfCodon += $nCoverage;
            }
            if (($aCoverages[$aGene['chr']][$nPosition]/$nCoverageOfCodon) < $_SETT['peak_finding']['codon_1st_pos_min_coverage_fraction']) {
                // Coverage is NOT higher than coverage of previous position.
                continue;
            }

            // Step 2b: ... and a clear "harringtonine pattern".
            // The following 5 codons should also not have a higher maximum coverage than this position's coverage.
            // If any of the following 5 codons have a maximum coverage higher than 10% of the coverage of the position we're analyzing,
            // that codon must not show a conflicting triplet periodicity pattern (see 'codon_1st_pos_min_coverage_fraction' check above).
            for ($j = 1; $j <= $_SETT['peak_finding']['downstream_codons_to_check']; $j++) {
                $nCoverageOfCodon = 0;
                $nMaxCoverage = 0;
                $sPositionFirstOfCodon = '';
                for ($k = 0; $k <= 2; $k++) {
                    if ($b3UTR) {
                        $sPositionToCheck = '*' . ((int) substr($sPosition, 1) + (3*$j) + $k);
                    } else {
                        $sPositionToCheck = $sPosition + (3*$j) + $k;
                        if ($b5UTR && $sPositionToCheck >= 0) {
                            // We crossed the -1/1 border.
                            $sPositionToCheck++;
                        }
                    }
                    if (!$k) {
                        $sPositionFirstOfCodon = $sPositionToCheck;
                    }
                    $nCoverage = (!isset($aPositions[$sPositionToCheck])? 0 : $aCoverages[$aGene['chr']][$aPositions[$sPositionToCheck]]);
                    $nCoverageOfCodon += $nCoverage;
                    $nMaxCoverage = max($nMaxCoverage, $nCoverage);
                }
                if ($nMaxCoverage > $aCoverages[$aGene['chr']][$nPosition]) {
                    // Max coverage of downstream codon IS higher.
                    continue 2;
                } elseif ($nMaxCoverage > ($_SETT['peak_finding']['downstream_codons_max_coverage_max_fraction']*$aCoverages[$aGene['chr']][$nPosition])) {
                    // Max coverage of downstream codon is higher than 10% of the current codon. This codon should now NOT have a conflicting periodicity.
                    $nCoverage = (!isset($aPositions[$sPositionFirstOfCodon])? 0 : $aCoverages[$aGene['chr']][$aPositions[$sPositionFirstOfCodon]]);
                    if (($nCoverage/$nCoverageOfCodon) < $_SETT['peak_finding']['codon_1st_pos_min_coverage_fraction']) {
                        // Codon does NOT show correct triplet periodicity.
                        continue 2;
                    }
                }
            }

            // This position has passed all tests and will be reported.
            if (!isset($aResults[$sGene])) {
                $aResults[$sGene] = array('report_separately' => 0, 'positions' => array());
            }
            if (!isset($aResults[$sGene]['positions'][$nPosition])) {
                $aResults[$sGene]['positions'][$nPosition] = array();
            }
            $aResults[$sGene]['positions'][$nPosition][$sTranscript] = $sPosition;

            // Mark as potential false positive when too far downstream, but not in the 3'UTR,
            // because of the background in the coding region >= 5KB from the pORF start site.
            if (!$b5UTR && !$b3UTR && $sPosition >= $_SETT['peak_finding']['candidate_min_distance_reported_separately']) {
                $aResults[$sGene]['report_separately'] ++;
                $aResults[$sGene]['positions'][$nPosition]['report_separately'] = true;
            }

            // Continue 5 codons downstream.
            // FIXME: How to nicely skip the next 5 codons?
        }
    }
    $aPositionsPerGene[$sGene]['unique_positions'] = array_unique($aPositionsPerGene[$sGene]['unique_positions']);
    $aPositionsPerGene[$sGene]['unique_positions_analyzed'] = array_unique($aPositionsPerGene[$sGene]['unique_positions_analyzed']);

    // Step 3: Per gene, from all its found possible ORF starts, we will take the one with the highest coverage as a reference, and discard any
    // other candidate ORF starting points that do not have at least a coverage (on that position) of 10% of the reference (highest candidate).
    if (!isset($aResults[$sGene])) {
        continue;
    }

    $nMaxCoverage = 0;
    foreach ($aResults[$sGene]['positions'] as $nPosition => $aPosition) {
        $nMaxCoverage = max($nMaxCoverage, $aCoverages[$aGene['chr']][$nPosition]);
    }
    foreach ($aResults[$sGene]['positions'] as $nPosition => $aPosition) {
        if ($aCoverages[$aGene['chr']][$nPosition] < $_SETT['peak_finding']['candidate_coverage_min_fraction']*$nMaxCoverage) {
            // Coverage is not enough to be reported!
            if (isset($aResults[$sGene]['positions'][$nPosition]['report_separately'])) {
                $aResults[$sGene]['report_separately'] --;
            }
            unset($aResults[$sGene]['positions'][$nPosition]);
        }
    }

    // Step 4: Genes that have no candidate peaks left, are ignored.
    if (!count($aResults[$sGene]['positions'])) {
        unset($aResults[$sGene]);
    }
}
$nGenesLeft = count($aResults);
$nGenesDiscarded = count($aPositionsPerGene) - $nGenesLeft;
print(' done.' . "\n");
fputs($aFilesOut['stats']['handler'], 'Genes discarded, no peaks found:' . "\t" . $nGenesDiscarded . "\n" .
    'Genes left with found translation start sites:' . "\t" . $nGenesLeft . "\n");
////////////////////////////////////////////////////////////////////////////////
/*
+ For all genes remaining, the positions of the candidate ORF start sites are
  reported relative to all known transcripts of the gene in question. The
  original number of positions with high coverage (>10) is mentioned along with
  the number of remaining candidate ORF start sites.
+ Additionally, a separate file reports all candidate peaks found at a position
  higher than 5Kb from a known ORF start site, so they can easily be checked
  more thoroughly if they are false positives.

*/





fputs($aFilesOut['results_after_cutoff']['handler'], 'The following results are found in the coding region, at a distance to the TIS of at least ' . $_SETT['peak_finding']['candidate_min_distance_reported_separately'] . ' nucleotides.' . "\n\n");
foreach ($aResults as $sGene => $aGene) {
    $nOriCandidates = count($aPositionsPerGene[$sGene]['unique_positions']);
    $nCandidatesAnalyzed = count($aPositionsPerGene[$sGene]['unique_positions_analyzed']);
    $nResults = count($aResults[$sGene]['positions']);
    $aTranscripts = array_keys($aPositionsPerGene[$sGene]['positions']);
    if ($aGene['report_separately']) {
        // Gene has at least one position that needs to be reported separately.
        fputs($aFilesOut['results_after_cutoff']['handler'], $sGene . "\t" . 'Positions found:' . "\t" . $nOriCandidates . "\t" . 'Positions analyzed:' . "\t" . $nCandidatesAnalyzed . "\t" . 'TIS found:' . "\t" . $nResults . "\n" .
            'G_Position' . "\t" . 'Coverage' . "\t" . implode("\t", $aTranscripts) . "\n");
    }
    if ($aGene['report_separately'] < count($aGene['positions'])) {
        // Gene has at least one position that needs to be reported in the normal results file.
        fputs($aFilesOut['results']['handler'], $sGene . "\t" . 'Positions found:' . "\t" . $nOriCandidates . "\t" . 'Positions analyzed:' . "\t" . $nCandidatesAnalyzed . "\t" . 'TIS found:' . "\t" . $nResults . "\n" .
            'G_Position' . "\t" . 'Coverage' . "\t" . implode("\t", $aTranscripts) . "\n");
    }
    if ($sArgStrand == '+') {
        ksort($aGene['positions']);
    } else {
        krsort($aGene['positions']);
    }
    foreach ($aGene['positions'] as $nPosition => $aPosition) {
        $sFileHandler = (isset($aPosition['report_separately'])? 'results_after_cutoff' : 'results');
        fputs($aFilesOut[$sFileHandler]['handler'], $aPositionsPerGene[$sGene]['chr'] . ':' . $nPosition . "\t" . $aCoverages[$aPositionsPerGene[$sGene]['chr']][$nPosition]);
        foreach ($aTranscripts as $sTranscript) {
            if (isset($aPosition[$sTranscript])) {
                fputs($aFilesOut[$sFileHandler]['handler'], "\t" . $aPosition[$sTranscript]);
            } else {
                fputs($aFilesOut[$sFileHandler]['handler'], "\t-");
            }
        }
        fputs($aFilesOut[$sFileHandler]['handler'], "\n");
    }
    if ($aGene['report_separately']) {
        fputs($aFilesOut['results_after_cutoff']['handler'], "\n");
    }
    if ($aGene['report_separately'] < count($aGene['positions'])) {
        fputs($aFilesOut['results']['handler'], "\n");
    }
}

print("\n");
exit(0);
?>
