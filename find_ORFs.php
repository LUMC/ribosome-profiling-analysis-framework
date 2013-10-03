#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * ORF-Finder analyzes the Mutalyzer results, compared to the original wiggle
 * files, checks the strand of the transcript in a third file, and finds ORFs
 * that have not been annotated before.
 *
 * Created     : 2013-07-12
 * Modified    : 2013-10-03
 * Version     : 0.3
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.3',
        'min_coverage' => 3,
        'max_upstream' => 500,
        'max_downstream' => 500,
        'min_positions_per_gene' => 3, // A known gene with UORF has only 4 positions above minimum coverage, so we need to keep this low.
        'peak_finding' =>
            array(
                'neighbour_min_difference' => 10,  // Should be at least # higher than the genomic neighbour.
                'neighbour_min_factor' => 2,       // Should be at least #x higher than the genomic neighbour.
                'region_average_upstream' => 5,    // Region starts # bases upstream.
                'region_average_downstream' => 5,  // Region ends # bases downstream.
                'region_average_min_factor' => 10, // Should be at least #x higher than the average coverage of the region.
                'exonic_average_min_factor' => 2,  // Exonic reads should be at least #x higher than the average coverage of the exonic region.
                'exonic_report_max_peaks' => 2,    // Report the # highest peaks in the exonic region, discard the rest.
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
    if (preg_match('/^variableStep chrom=(chr.+)$/', $sLine, $aRegs)) {
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
    $aTranscripts[$sTranscript] = array($sGene, $sStrand);
}
print('done, loaded ' . count($aTranscripts) . ' transcripts in memory.' . "\n");





// Now, loop the Mutalyzer file, find coverages in array (if not present, complain), check strand in gene list, and summarize results.
$nFiltered = 0;
$nIntronicPositions = 0;
$nUnmappable = 0;
$nIntergenic = 0;
$aPositionsPerGene = array(); // Will contain genes as keys, with an array of transcripts and positions, which is an array of positions and their mappings on the transcripts.
$aTranscriptsPerPosition = array(); // Stores how many positions have how many transcripts, pure statistics.
$aCodonPositions = array();
$aUnknownTranscripts = array();
$aMutalyzerResults = file($aFiles[0]);
unset($aMutalyzerResults[0]); // Header.
print('Parsing mutalyzer results file... ');
foreach ($aMutalyzerResults as $sLine) {
    $aLine = explode("\t", rtrim($sLine)); // Removing whitespace from the right.
    $sVariant = array_shift($aLine);
    list($sChr, $nPosition) = explode(';', preg_replace('/^chr([0-9]+|[XYM]):g\.(\d+)del$/', "$1;$2", $sVariant));

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
        $bExonic = false;
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
            // Check strand and then check if position is exonic. If so, store position!
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
                        $bExonic = true;
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

        // If we have 5' and 3'UTR values, but no exonic values, compare the 5' and 3' positions.
        if ($b5UTR && !$bExonic && $b3UTR) {
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



        // Now that we're done filtering, save all positions per gene, so we can loop through it and try and find patters.
        foreach ($aCodonOptions as $sTranscript => $sPosition) {
            list($sGene, $sStrand) = $aTranscripts[$sTranscript];
            // Create gene array if it doesn't exist.
            if (!isset($aPositionsPerGene[$sGene])) {
                $aPositionsPerGene[$sGene] = array('chr' => $sChr, 'strand' => $sStrand, 'transcripts' => array(), 'positions' => array());
            }
            // If we don't know this transcript yet, add it to the list.
            if (!in_array($sTranscript, $aPositionsPerGene[$sGene]['transcripts'])) {
                $aPositionsPerGene[$sGene]['transcripts'][$sTranscript] = array('estimated_length' => 0, 'average_inframe_coverage' => 0);
            }
            // Create position's data array, if it didn't exist yet.
            if (!isset($aPositionsPerGene[$sGene]['positions'][$nPosition])) {
                $aPositionsPerGene[$sGene]['positions'][$nPosition] = array('coverage' => $nCoverage, 'mappings' => array());
            }
            $aPositionsPerGene[$sGene]['positions'][$nPosition]['mappings'][$sTranscript] = $sPosition;
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
print('done, ' . $nGenesMapped . ' genes with mappings.' . "\n");
// Filtering the gene list for low number of positions.
foreach ($aPositionsPerGene as $sGene => $aGene) {
    if (count($aGene['positions']) < $_SETT['min_positions_per_gene']) {
        unset($aPositionsPerGene[$sGene]);
    }
}
$nGenesMappedFiltered = count($aPositionsPerGene);

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
      'Positions left:' . "\t" . (array_sum($aTranscriptsPerPosition) - $nUnmappable - $nIntronicPositions - $nIntergenic) . "\n" .
      'Genes mapped to:' . "\t" . $nGenesMapped . "\n" .
      'Genes left after filtering:' . "\t" . $nGenesMappedFiltered . "\t" . 'Genes with too few positions (count < ' . $_SETT['min_positions_per_gene'] . ') were ignored' . "\n" .
      'Positions left:' . "\t");
$nPositions = 0;
foreach ($aPositionsPerGene as $aGene) {
    $nPositions += count($aGene['positions']);
}
print($nPositions . "\n" .
      'Now looking for peaks');



// Now we'll look for peaks. Several methods are used.
// We'll loop all positions (candidates) and check their absolute and relative distance to their direct neighbours and their surrounding region, of configurable size.
// If the thresholds are not reached, the candidate is discarded. -12 is never discarded, because we want to show it always.
$nGenesDiscarded = 0;
$i = 0;
foreach ($aPositionsPerGene as $sGene => $aGene) {
    if (!(++$i%10)) {
        print('.');
    }

    $aCandidates = array_keys($aGene['positions']);
    foreach ($aCandidates as $nKey => $nPosition) {
        // Step 1: We always accept -12 as a candidate.
        $aPosition = $aGene['positions'][$nPosition];
        foreach ($aPosition['mappings'] as $sPosition) {
            if (in_array($sPosition, array(-15, -12, -9, -6, -3))) {
                continue 2;
            }
        }

        $nCoverage = $aPosition['coverage'];
        // Step 2: Compare neighbour_min_difference and neighbour_min_facor.
        // FIXME? Note: Ignoring strand right now.
        $nUpstreamNeighbour = (isset($aGene['positions'][$nPosition-1])? $aGene['positions'][$nPosition-1]['coverage'] : 0);
        $nDownstreamNeighbour = (isset($aGene['positions'][$nPosition+1])? $aGene['positions'][$nPosition+1]['coverage'] : 0);
        if ($nCoverage < ($nUpstreamNeighbour + $_SETT['peak_finding']['neighbour_min_difference']) ||
            $nCoverage < ($nDownstreamNeighbour + $_SETT['peak_finding']['neighbour_min_difference']) ||
            $nCoverage < ($nUpstreamNeighbour * $_SETT['peak_finding']['neighbour_min_factor']) ||
            $nCoverage < ($nDownstreamNeighbour * $_SETT['peak_finding']['neighbour_min_factor'])) {
            // Not a candidate anymore.
            unset($aCandidates[$nKey]);
            continue;
        }

        // Step 3: Compare with surrounding region, using region_average_upstream, region_average_downstream, and region_average_min_factor.
        // Calculate coverage in surrounding sequence.
        // FIXME? This region is based on genomic sequence, not RNA sequence.
        // FIXME? Note: Ignoring strand right now.
        $nCoverageSumOfRegion = 0;
        for ($j = -$_SETT['peak_finding']['region_average_upstream']; $j <= $_SETT['peak_finding']['region_average_downstream']; $j ++) {
            if (!$j) {
                // We're at the peak itself now; we'll ignore it for the calculation of the average.
                continue;
            }
            $nCoverageSumOfRegion += (isset($aGene['positions'][$nPosition+$j])? $aGene['positions'][$nPosition+$j]['coverage'] : 0);
        }
        if ($nCoverage < (($nCoverageSumOfRegion/($_SETT['peak_finding']['region_average_upstream']+$_SETT['peak_finding']['region_average_downstream'])) * $_SETT['peak_finding']['region_average_min_factor'])) {
            // Not a candidate anymore.
            unset($aCandidates[$nKey]);
            continue;
        }
    }



    // Step 4: Filter out exonic background.
    // From all in-frame, exonic positions (1, 4, 7, 10, 13, etc) we calculate the average. All peaks from this region are then compared with this average.
    // If the factor is not high enough, the peak is discarded.
    // First, we need to guess the length of the transcript. We do so, by looping from the highest position backwards, and store the peak position with the highest exonic position.
    $aPositions = array_keys($aGene['positions']);
    if ($aGene['strand'] == '+') {
        // Gene on + strand. Start at highest genomic position.
        rsort($aPositions);
    }
    foreach ($aGene['transcripts'] as $sTranscript => $aTranscript) {
        $nLength = 0;
        $nTotalInframeCoverage = 0;
        foreach ($aPositions as $nPosition) {
            if (!isset($aGene['positions'][$nPosition]['mappings'][$sTranscript])) {
                // This position does not have a mapping on this transcript.
                continue;
            }
            $sPosition = $aGene['positions'][$nPosition]['mappings'][$sTranscript];
            if ($sPosition{0} == '*') {
                // 3' of the stop, not there yet.
                continue;
            } elseif ($sPosition < 0) {
                // We're upstream now. Went too far, break.
                break;
            } else {
                // Exonic region. Store length of we don't already have it, and take coverage if we're inframe.
                if (!$nLength) {
                    $nLength = $sPosition;
                }
                if (!(($sPosition-1)%3)) {
                    // In frame.
                    $nTotalInframeCoverage += $aGene['positions'][$nPosition]['coverage'];
                } else {
                    continue;
                }
                break;
            }
        }
        if ($nLength > 0) {
            // We actually found ourselves a length.
            $aPositionsPerGene[$sGene]['transcripts'][$sTranscript]['estimated_length'] = $nLength;
            // So maybe we also have some measured inframe coverage.
            $aPositionsPerGene[$sGene]['transcripts'][$sTranscript]['average_inframe_coverage'] = ($nTotalInframeCoverage/$nLength)*3; // Not 100.00% accurate, we don't know if $nLength is inframe.
        }
    }



    foreach ($aCandidates as $nKey => $nPosition) {
        /*
        // FIXME: filter now on these values.
                'exonic_average_min_factor' => 2,  // Exonic reads should be at least #x higher than the average coverage of the exonic region.
                'exonic_report_max_peaks' => 2,    // Report the # highest peaks in the exonic region, discard the rest.
        */


        // TEMPORARY!!! FILTER OUT EVERYTHING THAT IS WITHIN THE ORF ONLY.
        $aPosition = $aGene['positions'][$nPosition];
        $bExonic = true;
        foreach ($aPosition['mappings'] as $sPosition) {
            if (substr($sPosition, 0, 1) == '-') {
                // Don't delete it.
                $bExonic = false;
                break;
            }
        }
        if ($bExonic) {
            // Only mapped to 1 >= pos <= *.
            // Not a candidate anymore.
            unset($aCandidates[$nKey]);
        }
    }




    // FIXME: My theory is that long transcripts didn't have time to finish translation. When the ribosome is stopped at the ATG site, the others more downstream continue to translate the RNA.
    //   However, after 10 minutes the RNA isolation protocol is initiated. With long transcripts (e.g. Col3a1) we see that the background noise in the exonic region is increased greatly at the end of the transcript.
    //   We could perhaps estimate the size of the transcript by noting the peak furthest in the transcript, and then calculate the average coverage per base (exonic positions), taking the number of transcripts
    //   of this length into account, not blindly dividing by the total number of transcripts. This figure will probably show an image similar to Col3a1, and this can be used for adding an additional filtering
    //   step, where the background noise at the far downstream end of the transcript is compared to the average coverage in this region of the transcript.





    $aPositionsPerGene[$sGene]['candidates'] = $aCandidates;

    // No candidates left? Goodbye.
    if (!$aCandidates) {
        unset($aPositionsPerGene[$sGene]);
        $nGenesDiscarded ++;
    }
}
print(' done.' . "\n" .
      'Genes discarded, no peaks found:' . "\t" . $nGenesDiscarded . "\n" .
      'Genes left with candidate peaks:' . "\t" . count($aPositionsPerGene) . "\n");




//$aGenesToAnalyze = array('Atf4', 'Sulf1', 'Rab23', 'Pdcl3', 'Map4k4', 'Col3a1', 'Bzw1', 'Nop58', 'Bmpr2', 'Cyp20a1', 'Nrp2', 'Eef1b2', 'Igfbp2', 'Arpc2', 'Rqcd1', 'Dnajb2');
// Longest transcripts (> 4KB, 8500 transcripts, 5419 gene symbols):
//$aGenesToAnalyze = array('0610010B08Rik', '0610010F05Rik', '0610030E20Rik', '100043387', '1110003E01Rik', '1110003O08Rik', '1110007A13Rik', '1110012J17Rik', '1110018G07Rik', '1110028C15Rik', '1110037F02Rik', '1200009O22Rik', '1200011I18Rik', '1200016B10Rik', '1300001I01Rik', '1300002K09Rik', '1300010F03Rik', '1500010J02Rik', '1520402A15Rik', '1700013B16Rik', '1700017B05Rik', '1700021K19Rik', '1700025G04Rik', '1700028N14Rik', '1700029I01Rik', '1700034H15Rik', '1700049G17Rik', '1700052N19Rik', '1700081L11Rik', '1810012P15Rik', '1810013D10Rik', '1810013L24Rik', '1810041L15Rik', '1810043G02Rik', '1810048J11Rik', '1810074P20Rik', '2010106G01Rik', '2210009G21Rik', '2210012G02Rik', '2210018M11Rik', '2310003F16Rik', '2310008H09Rik', '2310021P13Rik', '2310035C23Rik', '2310044G17Rik', '2310057J16Rik', '2310067B10Rik', '2410089E03Rik', '2410131K14Rik', '2510009E07Rik', '2610002M06Rik', '2610008E11Rik', '2610030H06Rik', '2610034M16Rik', '2610101N10Rik', '2610109H07Rik', '2610507B11Rik', '2700049A03Rik', '2700050L05Rik', '2700078E11Rik', '2810030E01Rik', '2810046L04Rik', '2810408P10Rik', '2810474O19Rik', '2900026A02Rik', '2900064A13Rik', '2900092E17Rik', '3110082D06Rik', '3425401B19Rik', '4631416L12Rik', '4632411B12Rik', '4632427E13Rik', '4632428N05Rik', '4732418C07Rik', '4732471D19Rik', '4831426I19Rik', '4833424O15Rik', '4921513D23Rik', '4922501C03Rik', '4922501L14Rik', '4930402E16Rik', '4930402H24Rik', '4930407I10Rik', '4930417M19Rik', '4930422G04Rik', '4930432E11Rik', '4930453N24Rik', '4930473A06Rik', '4930506M07Rik', '4930534B04Rik', '4930539E08Rik', '4930555K19Rik', '4931406P16Rik', '4931408C20Rik', '4932431P20Rik', '4932438A13Rik', '4932441K18Rik', '4933403F05Rik', '4933407C03Rik', '4933407H18Rik', '4933411K20Rik', '4933413G19Rik', '4933426M11Rik', '4933427D14Rik', '4933432B09Rik', '4933439F18Rik', '5031439G07Rik', '5330417C22Rik', '5430411K18Rik', '5730419I09Rik', '5730455P16Rik', '5730494M16Rik', '5730590G19Rik', '5830418K08Rik', '6230409E13Rik', '6330403A02Rik', '6330403L08Rik', '6330408A02Rik', '6330409N04Rik', '6330439K17Rik', '6330503K22Rik', '6430527G18Rik', '6430548M08Rik', '6430704M03Rik', '6530418L21Rik', '6720456H20Rik', '7120451J01Rik', '8430419L09Rik', '8430427H17Rik', '9030409G11Rik', '9030418K01Rik', '9030420J04Rik', '9030425E11Rik', '9030624J02Rik', '9130019P16Rik', '9130404D08Rik', '9330101J02Rik', '9330182L06Rik', '9430020K01Rik', '9430031J16Rik', '9430038I01Rik', '9530053A07Rik', '9630014M24Rik', '9830001H06Rik', '9930013L23Rik', '9930021J03Rik', '9930111J21Rik1', '9930111J21Rik2', 'A130022J15Rik', 'A130042E20Rik', 'A230046K03Rik', 'A230051G13Rik', 'A2bp1', 'A2M', 'A330021E22Rik', 'A430105I19Rik', 'A430107O13Rik', 'A430110N23Rik', 'A530054K11Rik', 'A530088H08Rik', 'A630007B06Rik', 'A630033E08Rik', 'A730008H23Rik', 'A730011L01Rik', 'A830010M20Rik', 'A830018L16Rik', 'A930001N09Rik', 'A930038C07Rik', 'A930039A15Rik', 'AA408296', 'AA987161', 'aak1', 'AARS', 'AASDH', 'AATK', 'ABAT', 'ABCA1', 'Abca12', 'ABCA13', 'Abca14', 'Abca15', 'Abca16', 'Abca17', 'ABCA2', 'ABCA3', 'abca4', 'abca5', 'Abca6', 'ABCA7', 'Abca8a', 'Abca8b', 'Abca9', 'ABCB10', 'abcb11', 'Abcb1a', 'Abcb1b', 'ABCB4', 'ABCB7', 'ABCC1', 'ABCC10', 'ABCC12', 'Abcc2', 'Abcc3', 'abcc4', 'abcc5', 'ABCC6', 'ABCC8', 'abcc9', 'ABCD2', 'abhd13', 'ABHD2', 'ABI2', 'ABI3', 'ABI3BP', 'ABL1', 'ABLIM1', 'Ablim3', 'ABR', 'ABTB2', 'ACACA', 'Acacb', 'Acad10', 'ACAD11', 'ACAD9', 'ACAN', 'ACAP2', 'ACAP3', 'ACCS', 'ACE', 'Acer2', 'ACER3', 'ACIN1', 'ACLY', 'Acn9', 'ACOT11', 'ACP2', 'acpP', 'ACSL1', 'Acsl4', 'Acsl6', 'Acsm2', 'ACTR10', 'ACVR1B', 'ACVR1C', 'ACVR2A', 'ADAM10', 'ADAM11', 'adam12', 'Adam17', 'adam19', 'ADAM22', 'ADAM23', 'Adamts1', 'ADAMTS10', 'ADAMTS12', 'Adamts13', 'ADAMTS14', 'Adamts15', 'ADAMTS16', 'ADAMTS17', 'ADAMTS18', 'ADAMTS19', 'ADAMTS2', 'ADAMTS20', 'Adamts3', 'ADAMTS4', 'ADAMTS5', 'ADAMTS6', 'Adamts7', 'ADAMTS9', 'ADAMTSL1', 'ADAMTSL4', 'Adamtsl5', 'Adar', 'Adarb1', 'ADARB2', 'ADCY1', 'Adcy10', 'ADCY2', 'ADCY3', 'adcy5', 'Adcy6', 'ADCY7', 'adcy8', 'ADCY9', 'ADCYAP1R1', 'Add1', 'ADD3', 'ADNP', 'adnp2', 'Ado', 'Adora1', 'ADRA1A', 'Adra2b', 'Adrbk2', 'Adrm1', 'AEBP2', 'AEN', 'AF529169', 'AFAP1', 'AFF1', 'AFF2', 'AFF3', 'AFF4', 'AFG3L1', 'AGAP1', 'Agap2', 'Agbl3', 'AGBL4', 'agl', 'AGPAT3', 'AGPAT5', 'Agphd1', 'AGPS', 'AGRN', 'AGTPBP1', 'AHCTF1', 'AHCYL2', 'AHDC1', 'Ahi1', 'ahnak', 'AHNAK2', 'ahr', 'AHRR', 'AI118078', 'AI314180', 'AI464131', 'AI481877', 'AI593442', 'AI646023', 'AI661453', 'AI848100', 'AIM1', 'Aim1l', 'AK3L1', 'Ak3l2-ps', 'AKAP11', 'AKAP12', 'AKAP13', 'akap2', 'AKAP6', 'Akap9', 'Akna', 'AKNAD1', 'akt2', 'Akt3', 'Alcam', 'Aldh16a1', 'ALDH1L1', 'Aldh1l2', 'ALDH4A1', 'ALDH5A1', 'alg11', 'ALK', 'Alkbh5', 'ALMS1', 'ALPK1', 'alpk2', 'ALPK3', 'ALS2', 'ALS2CL', 'alx4', 'AMBRA1', 'AMIGO1', 'AMMECR1', 'AMMECR1L', 'amot', 'AMOTL1', 'Amotl2', 'AMPD3', 'AMZ1', 'ANAPC1', 'ANGPT1', 'ank1', 'ank2', 'ANK3', 'ankar', 'Ankfy1', 'ANKHD1', 'ANKIB1', 'Ankle2', 'ankrd11', 'ANKRD12', 'ANKRD17', 'ANKRD26', 'Ankrd27', 'ANKRD28', 'ANKRD33B', 'ANKRD34B', 'ANKRD36', 'ANKRD44', 'ANKRD50', 'ankrd52', 'ANKRD57', 'ANKRD6', 'Anks1', 'Anln', 'ANO1', 'Ano3', 'ANO4', 'ANO5', 'Ano6', 'ANTXR1', 'ANTXR2', 'aoc3', 'AOX1', 'Aox3', 'Aox3l1', 'AOX4', 'AP1B1', 'Ap1g1', 'Ap2a2', 'AP2B1', 'AP3B1', 'Ap3d1', 'ap3s2', 'Ap4e1', 'Apaf1', 'apbA1', 'APBB1IP', 'APC', 'APC2', 'Apeg3', 'apex2', 'APH1B', 'Apob', 'APPBP2', 'appl1', 'APTX', 'AQP4', 'aqr', 'Ar', 'ARAP1', 'Arap2', 'arap3', 'Arcn1', 'ARFGEF1', 'ARFGEF2', 'ARFIP2', 'ARHGAP11A', 'ARHGAP12', 'ARHGAP15', 'ARHGAP19', 'ARHGAP20', 'ARHGAP21', 'ARHGAP23', 'ARHGAP24', 'arhgap26', 'ARHGAP28', 'ARHGAP29', 'ARHGAP30', 'ARHGAP31', 'ARHGAP33', 'ARHGAP4', 'ARHGAP5', 'ARHGAP6', 'Arhgef1', 'arhgef10', 'ARHGEF10L', 'ARHGEF11', 'ARHGEF12', 'ARHGEF15', 'arhgef16', 'ARHGEF17', 'ARHGEF18', 'Arhgef2', 'arhgef5', 'ARHGEF6', 'Arhgef7', 'ARHGEF9', 'ARID1A', 'ARID1B', 'ARID2', 'arid3a', 'ARID4A', 'ARID4B', 'Arid5a', 'ARID5B', 'arl5a', 'ARL5B', 'ARMC2', 'armc8', 'ARMC9', 'ARNT', 'ARNT2', 'ARPC5L', 'ARPP19', 'ARRB1', 'arrdc3', 'Arvcf', 'ASAH2', 'asap1', 'Asap2', 'ASAP3', 'ASB1', 'ASB15', 'Asb6', 'ASB7', 'ASCC2', 'ASCC3', 'ASH1L', 'aspA', 'asph', 'ASPM', 'ASRGL1', 'Astn1', 'Astn2', 'Asxl1', 'ASXL2', 'ASXL3', 'ATAD2', 'atad2b', 'ATAD5', 'ate1', 'Atf2', 'Atf6', 'ATF7IP', 'Atg2a', 'ATG2B', 'ATG4B', 'ATL3', 'ATM', 'ATMIN', 'ATN1', 'ATP10B', 'ATP10D', 'Atp11a', 'Atp11b', 'ATP11C', 'Atp13a2', 'ATP13A3', 'ATP13A4', 'ATP13A5', 'Atp1a2', 'ATP1B4', 'ATP2A2', 'ATP2A3', 'ATP2B1', 'ATP2B2', 'Atp2b3', 'atp2b4', 'ATP2C1', 'ATP5G1', 'ATP6V0A1', 'Atp6v1a', 'ATP7A', 'Atp7b', 'ATP8A1', 'Atp8b1', 'ATP8B2', 'ATP8B3', 'ATP8B4', 'atp9b', 'ATR', 'ATRN', 'ATRNL1', 'Atrx', 'ATXN1', 'ATXN1L', 'ATXN2', 'ATXN3', 'ATXN7', 'Atxn7l1', 'AU022252', 'AU040320', 'AU040829', 'AVEN', 'AVL9', 'AVPR1B', 'AW146020', 'AW549877', 'AW551984', 'AW554918', 'AW555464', 'Axin2', 'AXL', 'AZIN1', 'B230120H23Rik', 'B230219D22Rik', 'B3galt1', 'B3galt2', 'B3galt5', 'B3galtl', 'b3gat2', 'b3gnt5', 'B4galt5', 'B4GALT6', 'B630005N14Rik', 'BACE1', 'Bach1', 'Bach2', 'Bag4', 'BAHCC1', 'BAHD1', 'BAI1', 'BAI2', 'BAI3', 'BAIAP2', 'BAMBI', 'BARD1', 'BAT2', 'Bat2l', 'Bat2l2', 'BAZ1A', 'baz1b', 'Baz2a', 'Baz2b', 'BBOX1', 'Bbs1', 'Bbx', 'BC005561', 'BC006779', 'BC016423', 'BC018507', 'BC019943', 'BC021767', 'BC021891', 'BC026590', 'BC027072', 'BC028471', 'BC030307', 'BC030336', 'BC031353', 'BC032203', 'BC046331', 'BC048546', 'BC049265', 'BC057079', 'BC067074', 'BC106179', 'BCAT1', 'bchE', 'BCL11A', 'BCL11B', 'BCL2', 'BCL2L11', 'BCL2L13', 'Bcl2l2', 'BCL7A', 'Bcl7c', 'BCL9', 'BCL9L', 'BCLAF1', 'bcor', 'BCORL1', 'BDNF', 'BDP1', 'BECN1', 'BEND3', 'BHLHE41', 'BICC1', 'BICD1', 'BICD2', 'birc6', 'Blm', 'BLNK', 'Bmf', 'bmp2k', 'BMPR1A', 'Bmpr1b', 'bms1', 'Bnc1', 'Bnc2', 'BOC', 'BOD1L', 'BPTF', 'BRAF', 'BRCA1', 'Brca2', 'BRCC3', 'BRD1', 'BRD2', 'BRD3', 'brd4', 'BRD8', 'brdt', 'BRI3BP', 'Brip1', 'BRPF1', 'BRPF3', 'BRSK2', 'BRWD1', 'Brwd3', 'BSN', 'Btaf1', 'BTBD11', 'BTBD12', 'Btbd3', 'BTBD7', 'BTBD9', 'btg1', 'BTNL9', 'btrc', 'Butr1', 'BZRAP1', 'C030039L03Rik', 'C030046E11Rik', 'C130022K22Rik', 'C130039O16Rik', 'C1QTNF5', 'C230081A13Rik', 'C230096C10Rik', 'C2CD2L', 'C2CD3', 'c2cd4c', 'C3', 'C330019G07Rik', 'C3AR1', 'C4a', 'C4b', 'C530008M17Rik', 'C77080', 'C77370', 'C79407', 'cabc1', 'CABIN1', 'cables2', 'CACHD1', 'Cacna1a', 'Cacna1b', 'CACNA1C', 'Cacna1d', 'CACNA1E', 'Cacna1f', 'CACNA1G', 'CACNA1H', 'CACNA1S', 'Cacna2d1', 'CACNA2D2', 'Cacna2d4', 'CACNB4', 'CACNG2', 'cad', 'CADM1', 'CADM2', 'CADM3', 'CADPS', 'cadps2', 'CALB1', 'CALCRL', 'CALD1', 'CALM1', 'CALM2', 'CALM3', 'Camk1d', 'Camk2a', 'Camk2b', 'CAMK2D', 'Camk2n1', 'camk4', 'camkk2', 'CAMSAP1', 'camsap1l1', 'CAMTA1', 'Camta2', 'CAND1', 'CAND2', 'Cant1', 'CAPN5', 'CAPN6', 'CAPN7', 'Caprin1', 'CARD11', 'CARD14', 'Card6', 'CASC1', 'CASC4', 'CASC5', 'CASK', 'CASKIN1', 'caskin2', 'Casp8ap2', 'Casr', 'CASZ1', 'CBFA2T2', 'Cbfa2t3', 'cbl', 'Cblb', 'CBLL1', 'Cbln3', 'cbx4', 'CBX5', 'cbx6', 'Cbx6-Nptxr', 'CC2D2A', 'Ccar1', 'ccbe1', 'ccdc108', 'Ccdc116', 'ccdc127', 'CCDC14', 'Ccdc141', 'CCDC144B', 'CCDC148', 'Ccdc157', 'Ccdc18', 'ccdc24', 'CCDC40', 'Ccdc47', 'CCDC50', 'CCDC52', 'CCDC59', 'CCDC6', 'CCDC72', 'CCDC85A', 'CCDC88A', 'Ccdc88b', 'CCDC88C', 'CCDC90A', 'Ccdc93', 'Ccnb3', 'CCND2', 'CCNG2', 'ccnl1', 'Ccnl2', 'CCNT2', 'CCNY', 'CCPG1', 'CCRN4L', 'CD109', 'CD163', 'CD22', 'CD28', 'CD2AP', 'CD300A', 'CD44', 'CD47', 'CD93', 'Cdan1', 'CDC14A', 'CDC14B', 'CDC23', 'CDC26', 'cdc27', 'CDC42BPA', 'CDC42BPB', 'CDC42BPG', 'cdc6', 'Cdc73', 'CDCP1', 'CDH1', 'CDH11', 'Cdh12', 'CDH13', 'CDH2', 'CDH20', 'Cdh23', 'CDH4', 'CDH7', 'CDH8', 'CDHR1', 'CDHR2', 'CDK12', 'CDK13', 'cdk14', 'CDK19', 'Cdk5r1', 'CDK5RAP2', 'CDK7', 'CDKL2', 'CDKN1B', 'cdon', 'Cds2', 'CDV3', 'CDYL2', 'CEBPG', 'Cebpz', 'CECR2', 'Celf1', 'CELF2', 'CELF5', 'celsr1', 'CELSR2', 'Celsr3', 'CENPE', 'CENPJ', 'CENPL', 'CENPO', 'Cep110', 'CEP120', 'CEP135', 'CEP152', 'CEP164', 'CEP170', 'CEP192', 'CEP250', 'CEP350', 'Cep68', 'Cep97', 'CERK', 'CFH', 'CFLAR', 'CFTR', 'CGGBP1', 'cgn', 'cgnl1', 'CHD1', 'chd3', 'CHD4', 'CHD5', 'CHD6', 'chd7', 'Chd8', 'CHD9', 'chdh', 'CHIC1', 'Chic2', 'CHID1', 'CHKA', 'CHL1', 'CHM', 'CHMP1B', 'Chn1', 'CHPF2', 'CHPT1', 'CHRD', 'CHRDL1', 'CHRM1', 'chrna1', 'Chrna4', 'CHRNB2', 'CHRNB3', 'CHST11', 'CHST15', 'chst2', 'CHST3', 'CHSY1', 'Cic', 'Ciita', 'CILP', 'Cilp2', 'CIT', 'ckap5', 'Clasp1', 'CLASP2', 'Clcn3', 'Clcn4-2', 'Clcn5', 'CLCN6', 'CLCN7', 'cldn18', 'CLDN19', 'Clec14a', 'CLEC16A', 'Clec1b', 'CLIC4', 'Clic5', 'clint1', 'CLIP1', 'CLIP2', 'CLIP4', 'CLMN', 'CLN8', 'CLOCK', 'clpB', 'CLPTM1', 'clspn', 'CLSTN1', 'Clstn2', 'Clstn3', 'CLTC', 'CMAH', 'CMYA5', 'CNGB1', 'CNGB3', 'CNKSR2', 'CNNM3', 'cnot1', 'CNOT4', 'CNOT6', 'cnr1', 'Cnst', 'CNTF', 'Cntln', 'Cntn1', 'CNTN2', 'CNTN3', 'Cntn4', 'CNTNAP1', 'CNTNAP2', 'CNTNAP3', 'CNTNAP4', 'Cntnap5a', 'Cntnap5b', 'cobL', 'Cobll1', 'Cog3', 'COG8', 'COL11A1', 'Col11a2', 'col12a1', 'col14a1', 'Col15a1', 'Col16a1', 'Col17a1', 'COL18A1', 'col19a1', 'COL1A1', 'COL1A2', 'COL20A1', 'COL23A1', 'COL24A1', 'COL25A1', 'Col27a1', 'COL28A1', 'COL2A1', 'COL3A1', 'Col4a1', 'col4a2', 'COL4A3', 'COL4A3BP', 'COL4A4', 'COL4A5', 'Col4a6', 'Col5a1', 'Col5a2', 'Col5a3', 'COL6A1', 'COL6A2', 'Col6a3', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COLEC10', 'copA', 'CopG', 'Corin', 'coro2a', 'CORO7', 'COX15', 'cox4nb', 'CP', 'CPA6', 'cpamd8', 'cpd', 'CPEB2', 'CPEB3', 'CPEB4', 'CPLX2', 'CPM', 'CPN2', 'CPNE1', 'CPNE3', 'CPNE5', 'CPNE9', 'Cps1', 'CPSF1', 'cpsf2', 'Cpsf6', 'Cpt1a', 'Cpt1b', 'CR2', 'Cramp1l', 'Crat', 'Crb1', 'CRB2', 'CREB1', 'CREB3L2', 'CREBBP', 'Crebzf', 'CREG2', 'Crim1', 'CRISPLD1', 'CRISPLD2', 'crk', 'Crocc', 'CRTC1', 'crtc3', 'CRYBG3', 'CSDE1', 'csf1', 'csf1r', 'CSF2RB', 'CSF2RB2', 'csgalnact2', 'CSMD1', 'CSMD2', 'CSMD3', 'CSNK1G1', 'Csnk1g3', 'Csnk2a1', 'cspg4', 'CSPP1', 'CSRNP2', 'CSRNP3', 'cstf2', 'CTDSP2', 'Ctdspl', 'Ctdspl2', 'CTNND1', 'CTNND2', 'CTR9', 'CTSB', 'Cttnbp2', 'CTTNBP2NL', 'cubn', 'CUL5', 'CUL7', 'Cul9', 'CUX1', 'CUX2', 'cxadr', 'CXCL12', 'CXXC4', 'CYB561D1', 'CYB5B', 'Cyb5rl', 'CYBB', 'cybrd1', 'CYFIP1', 'CYFIP2', 'CYLD', 'CYP1B1', 'CYP26B1', 'Cyp2j13', 'Cyp4x1', 'Cyp7a1', 'CYSLTR1', 'CYTIP', 'CYTSA', 'Cytsb', 'CYYR1', 'D030016E14Rik', 'D030074E01Rik', 'D10Bwg1379e', 'D130043K22Rik', 'D14Abb1e', 'D15Ertd621e', 'D15Wsu169e', 'D17H6S56E-3', 'D19Ertd386e', 'D2Ertd750e', 'D3Bwg0562e', 'D3Ertd751e', 'D430041D05Rik', 'D430042O09Rik', 'D5Ertd579e', 'D630003M21Rik', 'D630013G24Rik', 'D630037F22Rik', 'D630045J12Rik', 'D6Wsu116e', 'D730040F13Rik', 'D8Ertd82e', 'D930015E06Rik', 'D930020E02Rik', 'D930049A15Rik', 'DAAM1', 'DAAM2', 'dab1', 'DAB2', 'DAB2IP', 'Dach1', 'DACH2', 'DAG1', 'dagla', 'DAP3', 'DAPK1', 'dcaf11', 'DCAF5', 'dcaf7', 'Dcbld1', 'dcbld2', 'DCC', 'Dcdc2a', 'DCHS1', 'DCLK1', 'Dclk2', 'DCLRE1A', 'DCLRE1B', 'DCLRE1C', 'DCP1A', 'DCP1B', 'dcp2', 'DCTN1', 'DCUN1D3', 'dcun1d4', 'DCX', 'DDB1', 'DDHD1', 'DDHD2', 'Ddi2', 'DDX17', 'DDX19B', 'DDX21', 'ddx26b', 'DDX3X', 'DDX3Y', 'DDX42', 'ddx46', 'ddx51', 'Ddx60', 'dennd1a', 'DENND2A', 'DENND2D', 'DENND3', 'DENND4A', 'dennd4b', 'DENND4C', 'DENND5A', 'DENND5B', 'Depdc5', 'Depdc6', 'Derl2', 'Dfna5', 'DGAT2L6', 'DGCR2', 'Dgcr8', 'DGKB', 'DGKD', 'DGKE', 'DGKI', 'DGKK', 'dgkq', 'DGKZ', 'DHCR24', 'Dhfr', 'dhrs9', 'DHTKD1', 'DHX29', 'DHX30', 'DHX33', 'DHX34', 'Dhx38', 'Dhx57', 'DHX8', 'dhx9', 'Diap1', 'Diap2', 'Diap3', 'DICER1', 'DIDO1', 'dio2', 'dip2a', 'DIP2B', 'DIRAS2', 'DIRC2', 'DISP1', 'DISP2', 'DIXDC1', 'dlaT', 'DLC1', 'DLEC1', 'DLG1', 'Dlg2', 'DLG3', 'DLG5', 'dlgap1', 'DLGAP2', 'DLGAP3', 'DLGAP4', 'DLK1', 'dlx1', 'DMBT1', 'DMBX1', 'dmd', 'dmrta1', 'DMXL1', 'DMXL2', 'DNA2', 'Dnahc1', 'Dnahc10', 'DNAHC11', 'Dnahc17', 'Dnahc2', 'Dnahc3', 'Dnahc5', 'Dnahc6', 'Dnahc7a', 'Dnahc7b', 'Dnahc9', 'Dnaja1', 'Dnaja1-ps', 'Dnajb14', 'DNAJB5', 'DNAJC10', 'DNAJC11', 'DNAJC13', 'DNAJC14', 'Dnajc16', 'Dnajc18', 'DNAJC25', 'DNAJC27', 'DNAJC3', 'dnajc5', 'DNAJC6', 'Dnalc1', 'DNASE1L3', 'Dnase2a', 'DNHD1', 'DNM1', 'DNM1L', 'DNM3', 'Dnm3os', 'DNMBP', 'dnmt1', 'DNMT3A', 'DNMT3B', 'DOC2B', 'DOCK1', 'dock10', 'DOCK11', 'DOCK2', 'dock3', 'DOCK4', 'Dock5', 'Dock6', 'DOCK7', 'DOCK8', 'dock9', 'dopey1', 'Dopey2', 'Dot1l', 'DPM1', 'DPP10', 'DPP4', 'Dpp6', 'Dpp8', 'DPY19L1', 'Dpy19l3', 'DPY19L4', 'DPYD', 'Dpysl2', 'Dpysl3', 'DPYSL5', 'DRP2', 'DSC2', 'dsc3', 'dscam', 'Dscaml1', 'dse', 'Dsel', 'Dsg1a', 'Dsg1c', 'DSG2', 'DSG3', 'dsp', 'Dspp', 'DST', 'DSTYK', 'dtl', 'Dtna', 'dtx1', 'dtx2', 'dtx3l', 'Dtx4', 'DUOX1', 'Duox2', 'DUSP11', 'DUSP16', 'dusp18', 'DUSP27', 'DUSP3', 'DUSP8', 'DVL2', 'Dvwa', 'DYNC1H1', 'DYNC1I2', 'DYNC1LI2', 'DYNC2H1', 'DYRK1A', 'Dyrk1c', 'DYSF', 'DYTN', 'DZIP1', 'DZIP1L', 'DZIP3', 'E130112L23Rik', 'E130120F12Rik', 'E130308A19Rik', 'E130309D14Rik', 'E230008N13Rik', 'E2f2', 'E2F3', 'E2F7', 'E330009J07Rik', 'E330016A19Rik', 'E430025E21Rik', 'EAF1', 'Ears2', 'EBF1', 'EBF2', 'EBF3', 'ECE1', 'Echs1', 'ECT2', 'eda', 'EDA2R', 'EDC3', 'EDC4', 'EDEM1', 'edem3', 'edil3', 'EEA1', 'eef2k', 'EFCAB5', 'EFNA5', 'EFNB2', 'EFR3A', 'EFR3B', 'Egf', 'EGFLAM', 'egfr', 'Egr1', 'EHBP1', 'Ehbp1l1', 'Ehd2', 'Ehf', 'EHMT1', 'EHMT2', 'EIF1AY', 'Eif2ak1', 'EIF2AK2', 'eif2ak3', 'EIF2AK4', 'EIF2C1', 'EIF2C2', 'EIF2C3', 'EIF2C4', 'Eif2s2', 'eif3a', 'EIF3J', 'Eif4ebp3', 'Eif4g1', 'EIF4G2', 'Eif4g3', 'Eif5a2', 'Eif5b', 'ELAC1', 'ELAVL1', 'ELAVL2', 'ELAVL3', 'ELAVL4', 'ELF1', 'Elf4', 'ELFN2', 'ELK3', 'ELMO2', 'ELMOD2', 'ELOVL6', 'elovl7', 'Elp4', 'ELTD1', 'EME2', 'EML1', 'eml4', 'EML5', 'Eml6', 'En2', 'Enah', 'ENAM', 'Enc1', 'ENGASE', 'ENPEP', 'Enpp1', 'Enpp4', 'Enpp5', 'ENTPD1', 'Entpd5', 'ENTPD7', 'EP300', 'ep400', 'Epas1', 'Epb4.1', 'Epb4.1l1', 'Epb4.1l2', 'Epb4.1l3', 'Epb4.1l4b', 'Epb4.1l5', 'Epb4.2', 'Epb4.9', 'EPC2', 'EPHA10', 'EPHA3', 'EPHA4', 'Epha5', 'EPHA7', 'Epha8', 'ephb1', 'ephb2', 'EPHB3', 'EPHB4', 'EPM2AIP1', 'EPN2', 'eprs', 'EPS15', 'Eps8', 'Ept1', 'Erbb2', 'Erbb2ip', 'Erbb3', 'Erbb4', 'ERC1', 'ERC2', 'ERCC4', 'ERCC5', 'Ercc6', 'ERCC6L', 'ereg', 'ergic2', 'Eri1', 'ERI2', 'ERMAP', 'ermp1', 'ERN1', 'ERO1L', 'Ero1lb', 'ERP29', 'Esco1', 'ESPL1', 'ESPNL', 'Esr1', 'esrrb', 'esrrg', 'Esyt2', 'ESYT3', 'etaa1', 'Etl4', 'etnk1', 'ETS1', 'etv1', 'ETV3', 'ETV6', 'Evc', 'EVC2', 'EVI5', 'evpl', 'EVX2', 'EXD1', 'Exo1', 'exoc1', 'Exoc2', 'EXOC3', 'EXOC4', 'Exoc8', 'EXPH5', 'ext1', 'EXTL1', 'extl3', 'Eya1', 'eya3', 'EYA4', 'EZH1', 'F5', 'F730047E07Rik', 'F8', 'FADS6', 'FAF1', 'FAF2', 'FAIM2', 'Fam100a', 'FAM102A', 'Fam102b', 'FAM110B', 'FAM115A', 'FAM115C', 'Fam116a', 'Fam116b', 'fam117b', 'FAM120A', 'Fam120b', 'FAM120C', 'Fam122b', 'Fam123b', 'FAM123C', 'Fam125b', 'FAM126A', 'FAM126B', 'FAM129A', 'Fam131b', 'FAM135A', 'fam135b', 'Fam13a', 'FAM149A', 'FAM160A1', 'Fam160b1', 'FAM160B2', 'Fam164c', 'Fam167a', 'fam168a', 'fam168b', 'FAM169A', 'fam171a1', 'Fam171b', 'Fam178a', 'FAM179B', 'Fam184b', 'Fam186a', 'Fam188b', 'fam189a1', 'FAM190A', 'FAM193B', 'FAM196B', 'fam198b', 'fam199x', 'Fam19a2', 'FAM19A3', 'FAM20B', 'FAM38A', 'FAM40B', 'Fam46a', 'FAM46C', 'FAM49A', 'FAM53A', 'FAM53B', 'FAM53C', 'FAM54B', 'FAM57A', 'FAM59A', 'Fam5b', 'Fam63b', 'fam65a', 'FAM65B', 'fam73a', 'FAM78A', 'Fam78b', 'FAM83G', 'Fam83h', 'FAM84B', 'Fam98c', 'FANCA', 'Fancd2', 'fanci', 'FANCM', 'far1', 'FARP1', 'FASN', 'FASTKD2', 'Fat1', 'Fat2', 'FAT3', 'FAT4', 'fbf1', 'FBLN2', 'FBLN5', 'Fbn1', 'Fbn2', 'fbxl12', 'FBXL14', 'FBXL17', 'Fbxl18', 'Fbxl20', 'fbxl3', 'Fbxl5', 'Fbxo10', 'fbxo11', 'FBXO18', 'Fbxo28', 'FBXO3', 'FBXO30', 'FBXO32', 'FBXO38', 'Fbxo40', 'Fbxo42', 'FBXO43', 'Fbxo45', 'Fbxo7', 'Fbxw11', 'FBXW7', 'FBXW8', 'FCGBP', 'FCHO2', 'FCHSD1', 'fchsd2', 'Fem1a', 'fem1b', 'FEM1C', 'Fer1l4', 'FERMT1', 'fgd1', 'FGD3', 'FGD5', 'Fgd6', 'fgf10', 'Fgfr1', 'FGFR2', 'fgfr3', 'Fhad1', 'FHDC1', 'fhod1', 'FHOD3', 'FIBCD1', 'fign', 'FILIP1', 'Fip1l1', 'FKBP15', 'FKTN', 'FLG2', 'fliI', 'FLNA', 'FLNB', 'FLNC', 'FLRT1', 'FLRT2', 'FLT1', 'Flt4', 'Fmn1', 'Fmn2', 'FMNL2', 'FMNL3', 'FMO2', 'FMO5', 'FMR1', 'fn1', 'Fnbp1', 'FNBP4', 'FNDC1', 'FNDC3A', 'FNDC3B', 'Fndc3c1', 'Fndc7', 'FNIP1', 'FNIP2', 'Fosl2', 'FOXJ2', 'Foxj3', 'Foxk1', 'Foxk2', 'FOXM1', 'foxn2', 'foxo1', 'FOXO3', 'FOXP1', 'FOXP2', 'FOXRED2', 'FRAS1', 'FREM1', 'FREM2', 'Frem3', 'FRK', 'frmd3', 'Frmd4a', 'FRMD4B', 'Frmd5', 'FRMD6', 'FRMD8', 'FRMPD1', 'Frmpd4', 'frs2', 'FRY', 'FRYL', 'FSD1L', 'FSTL4', 'Fstl5', 'FTO', 'FTSJD2', 'FUBP1', 'FUCA2', 'Furin', 'fus', 'fut9', 'FXC1', 'Fyco1', 'FYTTD1', 'fzd1', 'fzd3', 'fzd5', 'FZD6', 'FZD7', 'G2E3', 'G3BP1', 'G3BP2', 'gab1', 'GAB2', 'GABBR1', 'Gabbr2', 'gabpa', 'gabpb1', 'GABPB2', 'Gabra1', 'GABRA4', 'GABRB3', 'Gabre', 'GABRG1', 'Gabrg2', 'Gabrg3', 'GABRQ', 'GAD2', 'GAK', 'GALNT10', 'GALNT13', 'GALNT5', 'Galnt7', 'GALNTL2', 'GANC', 'GAPVD1', 'GAS2L3', 'gas7', 'GATAD2A', 'Gatsl2', 'GBF1', 'GBP4', 'Gbp6', 'GBX2', 'Gcap14', 'gcc1', 'GCC2', 'GCLM', 'gcn1l1', 'GCNT1', 'GCNT2', 'gcnt3', 'GDA', 'GDF11', 'GEMIN5', 'GFM1', 'GFOD1', 'gfod2', 'GFPT1', 'gfra1', 'Gga2', 'Ggnbp1', 'GGT5', 'GHR', 'GHSR', 'gigyf1', 'GIGYF2', 'GIT2', 'GK5', 'glcci1', 'glcE', 'GLDN', 'GLE1', 'GLI2', 'Gli3', 'GLP2R', 'GLS', 'GLT25D2', 'GLTSCR1', 'glyctk', 'Gm10031', 'Gm10039', 'Gm10060', 'Gm10193', 'Gm1027', 'Gm11295', 'Gm11677', 'Gm11843', 'Gm11847', 'Gm12026', 'Gm12185', 'Gm12241', 'Gm12337', 'Gm12824', 'Gm12986', 'Gm13034', 'Gm13138', 'Gm13139', 'Gm13161', 'Gm13242', 'Gm13251', 'Gm13298', 'Gm13328', 'Gm1337', 'Gm13416', 'Gm13767', 'Gm13768', 'Gm13886', 'Gm13910', 'Gm14200', 'Gm14308', 'Gm14326', 'Gm14373', 'Gm14391', 'Gm14398', 'Gm14410', 'Gm14430', 'Gm14434', 'Gm14443', 'Gm14698', 'Gm14730', 'Gm15235', 'Gm15528', 'Gm15542', 'Gm15583', 'Gm1564', 'Gm1568', 'Gm15725', 'Gm15800', 'Gm15801', 'Gm16515', 'Gm1826', 'Gm1966', 'Gm2007', 'Gm2058', 'Gm2082', 'Gm2348', 'Gm2381', 'Gm2542', 'Gm2587', 'Gm2719', 'Gm2780', 'Gm2802', 'Gm2991', 'gm2a', 'Gm318', 'Gm3193', 'Gm3203', 'Gm336', 'Gm347', 'Gm3555', 'Gm3650', 'Gm3655', 'Gm3711', 'Gm382', 'Gm3864', 'Gm3888', 'Gm4002', 'Gm4070', 'Gm4349', 'Gm4475', 'Gm4479', 'Gm4500', 'Gm4521', 'Gm4631', 'Gm4672', 'Gm4674', 'Gm4724', 'Gm4776', 'Gm4799', 'Gm4972', 'Gm5043', 'Gm5064', 'Gm5109', 'Gm5117', 'Gm5182', 'Gm5435', 'Gm5469', 'Gm5490', 'Gm5550', 'Gm5595', 'Gm5607', 'Gm5778', 'Gm5780', 'Gm5815', 'Gm5820', 'Gm5827', 'Gm5828', 'Gm5831', 'Gm5856', 'Gm5878', 'Gm5896', 'Gm5898', 'Gm5951', 'Gm606', 'Gm608', 'Gm6086', 'Gm6153', 'Gm6158', 'Gm6159', 'Gm6162', 'Gm6190', 'Gm6206', 'Gm6254', 'Gm6335', 'Gm6413', 'Gm6493', 'Gm6506', 'Gm6528', 'Gm6641', 'Gm6710', 'Gm6758', 'Gm6766', 'Gm6776', 'Gm6793', 'Gm7004', 'Gm7054', 'Gm7236', 'Gm7281', 'Gm7282', 'Gm7298', 'Gm7308', 'Gm7451', 'Gm7498', 'Gm7526', 'Gm7551', 'Gm7712', 'Gm7743', 'Gm7749', 'Gm7816', 'Gm8000', 'Gm8007', 'Gm8086', 'Gm8177', 'Gm8188', 'Gm8258', 'Gm8261', 'Gm8291', 'Gm8545', 'Gm8616', 'Gm8750', 'Gm8752', 'Gm8783', 'Gm8818', 'Gm884', 'Gm8849', 'Gm8892', 'Gm8899', 'Gm8902', 'Gm8974', 'Gm8991', 'Gm8995', 'Gm9108', 'Gm9165', 'Gm9174', 'Gm9242', 'Gm9430', 'Gm9529', 'Gm9608', 'Gm9613', 'Gm9774', 'Gm9781', 'Gm98', 'Gm9845', 'Gm996', 'GMEB1', 'GMEB2', 'gmfb', 'GMPS', 'Gna13', 'gnal', 'GNAQ', 'GNB2L1', 'GNB4', 'gne', 'GNG12', 'Gnl3l', 'GNPTAB', 'Golga1', 'GOLGA2', 'GOLGA3', 'golga4', 'GOLGB1', 'Golim4', 'GON4L', 'GOPC', 'GOSR1', 'gpatch2', 'GPATCH8', 'GPC6', 'GPCPD1', 'GPD1L', 'Gpd2', 'GPM6B', 'GPN1', 'Gpr101', 'GPR107', 'Gpr112', 'Gpr116', 'GPR123', 'GPR124', 'GPR125', 'GPR126', 'GPR133', 'GPR146', 'Gpr149', 'GPR155', 'GPR156', 'GPR157', 'GPR158', 'Gpr17', 'gpr173', 'GPR174', 'GPR179', 'GPR22', 'Gpr26', 'GPR64', 'Gpr83', 'GPR89', 'GPR98', 'GPRASP1', 'Gprc5b', 'GPRIN1', 'Gramd1b', 'Gramd2', 'GRAMD4', 'GRB10', 'GREB1', 'GREB1L', 'GREM1', 'GRHL2', 'GRIA1', 'GRIA2', 'Gria3', 'GRIA4', 'GRID1', 'Grid2', 'grik2', 'Grik3', 'grik4', 'GRIN1', 'GRIN2A', 'GRIN2B', 'grin2c', 'GRIN2D', 'Grin3a', 'GRIP1', 'GRK1', 'GRLF1', 'GRM1', 'Grm4', 'Grm5', 'GRM6', 'GRM8', 'GRPEL2', 'Gse1', 'GSK3B', 'gtf2a1', 'gtf2i', 'GTF2IRD1', 'GTF3C1', 'GTF3C2', 'GTF3C4', 'GUCY1A2', 'GUCY1A3', 'Gucy2e', 'GUCY2F', 'Guf1', 'GVIN1', 'GXYLT1', 'GYK', 'Gyltl1b', 'gzf1', 'H13', 'H2AFY', 'H6pd', 'hadhb', 'hal', 'HAPLN1', 'HAS2', 'HAS3', 'haus6', 'Hc', 'Hcfc1', 'Hcn3', 'HDAC5', 'HDAC6', 'Hdac7', 'hdac9', 'HDC', 'HDGFRP3', 'HDLBP', 'HDX', 'HEATR1', 'heatr5a', 'HEATR5B', 'HEATR6', 'HEATR7A', 'HEATR7B1', 'Heatr7b2', 'Hectd1', 'Hectd2', 'hectd3', 'HECW1', 'Hecw2', 'HEG1', 'helB', 'Hells', 'helz', 'HEPH', 'HERC1', 'HERC2', 'Herc3', 'herc4', 'HERC5', 'heyl', 'HFM1', 'Hhip', 'Hic2', 'HIF1A', 'Hif1an', 'HIF3A', 'HIP1', 'Hip1r', 'hipk1', 'Hipk3', 'HIRA', 'Hivep1', 'hivep2', 'Hivep3', 'Hjurp', 'HK1', 'Hk2', 'HLCS', 'HLF', 'hltf', 'Hmcn1', 'Hmcn2', 'HMGCR', 'Hmgxb3', 'HMGXB4', 'Hnf4a', 'HNF4G', 'HNRNPA2B1', 'Hnrnpa3', 'Hnrnph1', 'HNRNPR', 'HNRNPUL2', 'HOMER1', 'homer2', 'homez', 'HOOK1', 'Hook3', 'HOXB4', 'HOXD11', 'HP1BP3', 'HPCAL4', 'Hps3', 'HPS5', 'HR', 'HRK', 'HRNR', 'Hs2st1', 'Hs3st3b1', 'HS6ST2', 'hsd17b7', 'HSF5', 'HSPA12A', 'hspa13', 'HSPA14', 'HSPA4', 'Hspa4l', 'Hspbap1', 'Hspg2', 'HSPH1', 'HTR1A', 'htr2c', 'HTR4', 'htr5a', 'HTT', 'Hunk', 'HUS1', 'HUWE1', 'HYDIN', 'Hyls1', 'hyou1', 'iars', 'iars2', 'Ibtk', 'ICK', 'icmt', 'IDE', 'ids', 'IFFO1', 'IFFO2', 'IFI44', 'IFIH1', 'IFT122', 'IFT140', 'IFT172', 'IFT80', 'IGBP1', 'IGDCC4', 'Igf1', 'IGF1R', 'IGF2', 'IGF2BP1', 'Igf2bp2', 'igf2bp3', 'Igf2r', 'IGFBP5', 'Igfn1', 'IGSF1', 'IGSF10', 'igsf3', 'IGSF9', 'IGSF9B', 'IKBKAP', 'IKBKB', 'ikbkg', 'ikzf1', 'ikzf2', 'IKZF4', 'IKZF5', 'IL16', 'IL17RD', 'IL18RAP', 'Il1r1', 'il1rap', 'IL1RAPL1', 'IL1RL1', 'IL28RA', 'IL2RA', 'Il4ra', 'IL6ST', 'IL9R', 'ILDR2', 'IMPAD1', 'IMPG2', 'INADL', 'INF2', 'Inf2q', 'ing5', 'Inhbb', 'Ino80', 'INO80D', 'INPP4A', 'INPP4B', 'Inpp5d', 'Inpp5e', 'INPP5F', 'INPPL1', 'INSR', 'INSRR', 'INTS1', 'INTS10', 'INTS2', 'INTS3', 'INTS6', 'INTS7', 'INTU', 'Invs', 'IP6K1', 'IPCEF1', 'IPMK', 'IPO11', 'IPO5', 'IPO8', 'Ipo9', 'IQCE', 'IQGAP2', 'Iqgap3', 'Iqsec1', 'IQSEC2', 'IQSEC3', 'IREB2', 'IRF2BP2', 'IRF4', 'irf6', 'Irgq', 'IRS1', 'irs2', 'IRS4', 'ITCH', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA3', 'ITGA4', 'ITGA5', 'itga6', 'ITGA7', 'ITGA8', 'ITGA9', 'ITGAD', 'ITGAL', 'ITGAM', 'Itgav', 'Itgax', 'ITGB3', 'ITGB4', 'ITGB6', 'ITGB8', 'ITGBL1', 'itih5', 'ITK', 'Itpkb', 'ITPR1', 'ITPR2', 'ITPR3', 'ITSN1', 'Itsn2', 'IWS1', 'JAG1', 'JAG2', 'JAK1', 'Jak2', 'JAK3', 'Jam2', 'jarid2', 'Jhdm1d', 'JMJD1C', 'JMJD4', 'JMJD6', 'JMJD7', 'JMY', 'JPH1', 'jph2', 'Jph4', 'Jrk', 'KALRN', 'KANK1', 'KANK2', 'Kank4', 'KAT2B', 'Kbtbd11', 'KBTBD7', 'KBTBD8', 'KCNA1', 'KCNA4', 'KCNA6', 'Kcnc1', 'Kcnc2', 'Kcnc3', 'KCND1', 'KCND2', 'Kcnd3', 'KCNG1', 'KCNH1', 'kcnh2', 'KCNH5', 'KCNH7', 'KCNH8', 'kcnj10', 'KCNJ15', 'KCNJ2', 'KCNJ3', 'Kcnj5', 'kcnk6', 'KCNMA1', 'Kcnmb1', 'KCNN1', 'KCNQ2', 'KCNQ5', 'KCNS2', 'KCNT1', 'KCNT2', 'KCNV1', 'KCNV2', 'Kcp', 'kctd12b', 'KCTD7', 'KDM1B', 'KDM2A', 'kdm2b', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'Kdm6a', 'KDM6B', 'KDR', 'kdsr', 'KIDINS220', 'KIF11', 'Kif13a', 'Kif13b', 'KIF14', 'KIF15', 'KIF16B', 'KIF1A', 'KIF1B', 'Kif1c', 'Kif20b', 'kif21a', 'Kif21b', 'Kif24', 'KIF26A', 'KIF27', 'KIF2A', 'kif3b', 'KIF3C', 'Kif4', 'Kif5a', 'Kif5b', 'KIF5C', 'Kif7', 'KIFAP3', 'KIRREL', 'KIT', 'Kitl', 'KL', 'KLF13', 'klf6', 'KLF7', 'klf8', 'KLHDC10', 'KLHDC5', 'Klhdc7a', 'Klhdc8a', 'Klhl1', 'KLHL14', 'KLHL15', 'KLHL18', 'KLHL20', 'KLHL23', 'KLHL24', 'KLHL29', 'KLHL31', 'KLHL32', 'KLHL9', 'KNDC1', 'KNTC1', 'KPNA1', 'kpna3', 'KPNA4', 'KPNA6', 'KPNB1', 'Kras', 'KRBA1', 'kremen1', 'KRIT1', 'Krt7', 'KSR1', 'Ktn1', 'KY', 'L1cam', 'L3mbtl3', 'LAMA1', 'Lama2', 'LAMA3', 'LAMA4', 'LAMA5', 'Lamb1-1', 'lamb2', 'Lamb3', 'LAMC1', 'LAMC2', 'LAMC3', 'lancl1', 'LARGE', 'LARP1', 'Larp4', 'LARP4B', 'LAS1L', 'lass5', 'LASS6', 'lats1', 'Lats2', 'LCA5', 'LCLAT1', 'Lcor', 'LCORL', 'Lct', 'Ldb3', 'Ldlr', 'LDOC1L', 'LEMD3', 'LENG8', 'LEPR', 'LETM1', 'lgi1', 'Lgi4', 'LGR4', 'LGR5', 'LHFPL2', 'lhx1', 'LHX4', 'Lhx9', 'LIF', 'Lifr', 'lig3', 'LIG4', 'LIMA1', 'LIMCH1', 'LIMD1', 'LIME1', 'LIMS1', 'Lin28b', 'LIN54', 'Lin7a', 'LIN7C', 'Lins2', 'Llgl1', 'Lman2', 'LMAN2L', 'LMBR1', 'LMBRD1', 'LMLN', 'Lmo7', 'LMTK2', 'LMTK3', 'lmx1b', 'LNPEP', 'LNX2', 'LOC100039037', 'LOC100039091', 'LOC100039744', 'LOC100043920', 'LOC100043982', 'LOC100043990', 'LOC100043992', 'LOC100043998', 'LOC100044008', 'LOC100044024', 'LOC100044039', 'LOC100044040', 'LOC100044041', 'LOC100044042', 'LOC100044059', 'LOC100044065', 'LOC100044106', 'LOC100044115', 'LOC100044118', 'LOC100044121', 'LOC100044124', 'LOC100044125', 'LOC100044133', 'LOC100044160', 'LOC100044161', 'LOC100044162', 'LOC100044163', 'LOC100044202', 'LOC100044208', 'LOC100044212', 'LOC100044222', 'LOC100044229', 'LOC100044261', 'LOC100044275', 'LOC100044315', 'LOC100044321', 'LOC100044322', 'LOC100044324', 'LOC100044325', 'LOC100044332', 'LOC100044338', 'LOC100044341', 'LOC100044356', 'LOC100044363', 'LOC100044374', 'LOC100044431', 'LOC100044468', 'LOC100044493', 'LOC100044500', 'LOC100044533', 'LOC100044547', 'LOC100044566', 'LOC100044642', 'LOC100044643', 'LOC100044686', 'LOC100044699', 'LOC100044742', 'LOC100044747', 'LOC100044760', 'LOC100044767', 'LOC100044776', 'LOC100044793', 'LOC100044842', 'LOC100044843', 'LOC100044882', 'LOC100044937', 'LOC100044968', 'LOC100044976', 'LOC100045018', 'LOC100045020', 'LOC100045032', 'LOC100045035', 'LOC100045048', 'LOC100045054', 'LOC100045092', 'LOC100045121', 'LOC100045132', 'LOC100045191', 'LOC100045202', 'LOC100045207', 'LOC100045216', 'LOC100045228', 'LOC100045283', 'LOC100045284', 'LOC100045296', 'LOC100045355', 'LOC100045356', 'LOC100045359', 'LOC100045371', 'LOC100045419', 'LOC100045420', 'LOC100045440', 'LOC100045448', 'LOC100045507', 'LOC100045522', 'LOC100045542', 'LOC100045585', 'LOC100045622', 'LOC100045628', 'LOC100045680', 'LOC100045749', 'LOC100045758', 'LOC100045772', 'LOC100045775', 'LOC100045780', 'LOC100045795', 'LOC100045884', 'LOC100045900', 'LOC100045958', 'LOC100045981', 'LOC100045983', 'LOC100045990', 'LOC100046012', 'LOC100046025', 'LOC100046035', 'LOC100046044', 'LOC100046048', 'LOC100046051', 'LOC100046056', 'LOC100046080', 'LOC100046137', 'LOC100046205', 'LOC100046208', 'LOC100046214', 'LOC100046241', 'LOC100046259', 'LOC100046304', 'LOC100046323', 'LOC100046333', 'LOC100046343', 'LOC100046525', 'LOC100046556', 'LOC100046566', 'LOC100046682', 'LOC100046690', 'LOC100046693', 'LOC100046698', 'LOC100046701', 'LOC100046704', 'LOC100046735', 'LOC100046744', 'LOC100046766', 'LOC100046781', 'LOC100046804', 'LOC100046825', 'LOC100046895', 'LOC100046926', 'LOC100046932', 'LOC100046953', 'LOC100047016', 'LOC100047029', 'LOC100047052', 'LOC100047069', 'LOC100047076', 'LOC100047082', 'LOC100047093', 'LOC100047109', 'LOC100047131', 'LOC100047134', 'LOC100047143', 'LOC100047147', 'LOC100047167', 'LOC100047183', 'LOC100047187', 'LOC100047194', 'LOC100047223', 'LOC100047237', 'LOC100047238', 'LOC100047284', 'LOC100047318', 'LOC100047353', 'LOC100047360', 'LOC100047385', 'LOC100047386', 'LOC100047391', 'LOC100047419', 'LOC100047441', 'LOC100047480', 'LOC100047481', 'LOC100047492', 'LOC100047516', 'LOC100047562', 'LOC100047575', 'LOC100047577', 'LOC100047589', 'LOC100047603', 'LOC100047616', 'LOC100047645', 'LOC100047670', 'LOC100047674', 'LOC100047733', 'LOC100047738', 'LOC100047776', 'LOC100047789', 'LOC100047790', 'LOC100047837', 'LOC100047843', 'LOC100047857', 'LOC100047863', 'LOC100047868', 'LOC100047880', 'LOC100047888', 'LOC100047896', 'LOC100047936', 'LOC100047937', 'LOC100047941', 'LOC100047960', 'LOC100047967', 'LOC100047973', 'LOC100047997', 'LOC100048010', 'LOC100048018', 'LOC100048049', 'LOC100048079', 'LOC100048082', 'LOC100048123', 'LOC100048139', 'LOC100048143', 'LOC100048165', 'LOC100048232', 'LOC100048250', 'LOC100048313', 'LOC100048362', 'LOC100048376', 'LOC100048415', 'LOC100048439', 'LOC100048476', 'LOC100048488', 'LOC100048526', 'LOC100048528', 'LOC100048537', 'LOC100048557', 'LOC100048559', 'LOC100048579', 'LOC100048600', 'LOC100048616', 'LOC100048709', 'LOC100048740', 'LOC100048751', 'LOC100048759', 'LOC100048814', 'LOC100048818', 'LOC100048820', 'LOC100048828', 'LOC100048832', 'LOC100048834', 'LOC100048845', 'LOC100048853', 'LOC100048858', 'LOC100048863', 'LOC100048879', 'LOC544757', 'LOC545974', 'LOC547345', 'LOC547380', 'LOC547385', 'LOC624231', 'LOC630401', 'LOC630507', 'LOC630567', 'LOC630671', 'LOC630776', 'LOC631094', 'LOC631406', 'LOC631806', 'LOC632263', 'LOC632530', 'LOC632664', 'LOC632805', 'LOC632847', 'LOC633311', 'LOC633677', 'LOC633778', 'LOC634102', 'LOC634379', 'LOC634417', 'LOC635773', 'LOC635918', 'LOC636537', 'LOC636856', 'LOC636969', 'LOC638024', 'LOC638050', 'LOC638642', 'LOC638935', 'LOC639171', 'LOC639433', 'LOC640441', 'LOC640793', 'LOC664986', 'LOC669888', 'LOC670614', 'LOC672215', 'LOC673550', 'LOC674168', 'LOC674193', 'LOC674236', 'LOC674368', 'LOC674654', 'LOC674677', 'LOC674773', 'LOC674888', 'LOC675056', 'LOC675366', 'LOC675521', 'LOC675560', 'LOC675709', 'LOC675933', 'LOC676013', 'LOC676256', 'LOC676463', 'LOC676464', 'LOC676881', 'LOC677069', 'LOC677213', 'LOC677282', 'LOC677366', 'LOC677447', 'LOC677502', 'LOC677560', 'LOC677582', 'LOC677611', 'Lonrf2', 'LOXHD1', 'loxl2', 'LOXL3', 'LPAR2', 'lpar4', 'lpgat1', 'LPHN1', 'LPHN2', 'LPHN3', 'lpin1', 'LPIN2', 'Lpl', 'lpp', 'LRAT', 'Lrba', 'Lrch1', 'LRCH2', 'LRCH3', 'Lrfn5', 'lrig1', 'LRIG2', 'LRIG3', 'LRIT1', 'lrp1', 'LRP10', 'LRP12', 'Lrp1b', 'LRP2', 'LRP3', 'LRP4', 'LRP5', 'LRP6', 'LRP8', 'lrpprc', 'LRRC14', 'LRRC15', 'Lrrc16a', 'LRRC16B', 'Lrrc3', 'LRRC40', 'LRRC49', 'LRRC56', 'lrrc58', 'lrrc7', 'LRRC8A', 'LRRC8C', 'LRRC8D', 'LRRK1', 'Lrrk2', 'Lrrtm2', 'LRRTM4', 'lrsam1', 'LRTM1', 'Lrtm2', 'lsamp', 'LSM11', 'LSS', 'Ltbp1', 'LTBP2', 'ltbp3', 'LTBP4', 'LUC7L', 'LUZP1', 'Ly75', 'Lypla1', 'LYST', 'LZTFL1', 'MACF1', 'macrod2', 'MADD', 'maf', 'MAFG', 'MAGEL2', 'MAGI1', 'MAGI2', 'MAGI3', 'MAGT1', 'MALT1', 'MAML1', 'MAML3', 'MAMLD1', 'MAMSTR', 'Man1a', 'MAN1A2', 'MAN1C1', 'MAN2A1', 'MAN2A2', 'Man2b1', 'MANEA', 'MAOA', 'MAP3K1', 'MAP3K12', 'MAP3K13', 'MAP3K14', 'MAP3K15', 'MAP3K2', 'MAP3K4', 'MAP3K5', 'MAP3K6', 'MAP3K7', 'MAP3K9', 'MAP4K2', 'map4k3', 'Map4k5', 'MAPK1', 'MAPK4', 'MAPK6', 'Mapk7', 'mapk8', 'Mapk8ip2', 'MAPK8IP3', 'Mapk9', 'MAPKAPK5', 'MAPKBP1', 'Mapre1', 'MAPRE2', 'mapt', 'MARCH1', 'March4', 'MARCH6', 'MARCH7', 'MARCH8', 'MARCKS', 'MARK1', 'MARK2', 'masp2', 'MAST1', 'MAST2', 'MAST4', 'MASTL', 'MBD5', 'MBD6', 'Mbnl1', 'mbnl2', 'Mbnl3', 'MBP', 'mbtd1', 'MBTPS1', 'MBTPS2', 'MCART1', 'MCC', 'MCCC1', 'MCF2', 'MCF2L', 'Mcm10', 'MCM9', 'MCPH1', 'MCTP1', 'MCTP2', 'MDC1', 'MDGA1', 'MDGA2', 'Mdm4-ps', 'MDN1', 'ME3', 'MECOM', 'mecp2', 'med1', 'MED12', 'MED12L', 'MED13', 'Med13l', 'Med14', 'MED23', 'MED24', 'MEF2A', 'MEF2C', 'MEF2D', 'MEGF10', 'Megf11', 'MEGF6', 'Megf8', 'MEGF9', 'MEI1', 'Meis2', 'MEOX2', 'MERTK', 'MESDC2', 'MET', 'metap2', 'Mett10d', 'METTL13', 'Mfap1a', 'Mfap1b', 'Mfap3', 'MFHAS1', 'mfi2', 'mfn1', 'MFN2', 'MFRP', 'MFSD11', 'Mfsd7b', 'MGA', 'MGAM', 'MGAT3', 'MGAT4A', 'Mgat5b', 'mgea5', 'MIA3', 'MIB1', 'MICAL2', 'MICAL3', 'MICALL1', 'Mid2', 'MIER1', 'Mier3', 'MINK1', 'MITF', 'MKI67', 'mkl1', 'MKL2', 'MKLN1', 'Mlec', 'MLH3', 'Mll1', 'MLL2', 'MLL3', 'MLL5', 'mllt10', 'MLLT3', 'mllt4', 'MLLT6', 'Mlph', 'MLXIP', 'MME', 'MMGT1', 'Mmp12', 'mmp15', 'mmp16', 'MMP17', 'mmp24', 'MMRN1', 'mn1', 'mnt', 'Mobkl1a', 'MOBKL1B', 'MOBKL2B', 'MON1B', 'MON2', 'Morc2a', 'MORC3', 'MORC4', 'Mov10l1', 'Mpa2l', 'MPDZ', 'Mpeg1', 'MPHOSPH9', 'MPP2', 'MPP3', 'mpp5', 'MPP7', 'Mprip', 'MPV17', 'mras', 'MRC1', 'Mrc2', 'mre11a', 'Mrgprb1', 'MRPL19', 'MRPL48', 'mrps10', 'MRPS25', 'MRPS26', 'MRPS35', 'Mrps5', 'MRVI1', 'msh6', 'MSI2', 'Msl1', 'MSL2', 'MST1R', 'Mtap1a', 'Mtap1b', 'Mtap2', 'Mtap4', 'Mtap6', 'Mtap9', 'MTF1', 'MTF2', 'Mthfr', 'mtif3', 'MTMR1', 'MTMR10', 'Mtmr12', 'mtmr3', 'mtmr4', 'MTOR', 'mtr', 'mtss1', 'MTSS1L', 'Mttp', 'MTUS1', 'Mtus2', 'MUC4', 'Muc5ac', 'MUC6', 'mug1', 'Mum1', 'Mum1l1', 'MXD1', 'MXI1', 'mybbp1a', 'mybl1', 'MYBPC3', 'MYH1', 'myh10', 'MYH11', 'MYH13', 'MYH14', 'myh2', 'Myh3', 'MYH4', 'Myh6', 'Myh7', 'Myh7b', 'MYH8', 'MYH9', 'MYLK', 'MYLK4', 'MYNN', 'Myo10', 'Myo15', 'Myo15b', 'Myo16', 'Myo18a', 'MYO18B', 'myo19', 'MYO1B', 'Myo1C', 'MYO1D', 'MYO1E', 'MYO3A', 'MYO3B', 'MYO5A', 'MYO5B', 'MYO5C', 'myo6', 'Myo7a', 'Myo7b', 'myo9b', 'Myocd', 'Myof', 'myom1', 'myom2', 'MYOM3', 'Mypn', 'MYRIP', 'Mysm1', 'MYST3', 'MYST4', 'Myt1', 'MYT1L', 'N28178', 'N4BP1', 'N4bp2', 'N4BP2L2', 'naa15', 'NAA25', 'NAA30', 'Naca', 'NACAD', 'Nacc1', 'NACC2', 'Naif1', 'Naip1', 'Naip1-rs1', 'Naip2', 'Naip3', 'Naip4', 'Naip5', 'Naip6', 'Naip7', 'NALCN', 'NAMPT', 'napB', 'NARF', 'NARG2', 'Nat10', 'NAT8L', 'NAV1', 'NAV2', 'NAV3', 'NBAS', 'NBEA', 'NBEAL1', 'NBEAL2', 'NBR1', 'ncam2', 'NCAN', 'NCAPD2', 'ncapd3', 'NCAPG2', 'NCEH1', 'NCK1', 'NCKAP1', 'NCKAP1L', 'NCKAP5', 'NCKAP5L', 'ncl', 'NCOA1', 'Ncoa2', 'NCOA3', 'NCOA6', 'NCOA7', 'NCOR1', 'ncor2', 'NCS1', 'NDFIP2', 'Ndnl2', 'Ndor1', 'NDRG3', 'NDST3', 'NDST4', 'NDUFB2', 'Ndufb7', 'NEB', 'Nebl', 'NECAB1', 'NEDD4', 'Nedd4l', 'Nedd9', 'NEGR1', 'nek1', 'nek9', 'Nell1', 'NEO1', 'NES', 'NETO2', 'Neu4', 'Neurl1a', 'NEURL1B', 'neurl4', 'NF1', 'Nf2', 'Nfasc', 'NFAT5', 'nfatc1', 'NFATC2', 'Nfatc3', 'NFE2L1', 'NFIA', 'NFIB', 'NFIC', 'nfix', 'NFKB1', 'NFKBIZ', 'Nfrkb', 'NFX1', 'NHLRC3', 'NHS', 'nhsl1', 'Nhsl2', 'nid1', 'NID2', 'nin', 'NINL', 'NIPAL1', 'NIPAL3', 'nipbl', 'nisch', 'Nkain3', 'NKAP', 'NKD1', 'NKD2', 'NKIRAS1', 'NKTR', 'nlgn1', 'Nlgn2', 'NLGN3', 'nlk', 'nlrc3', 'NLRC5', 'NLRP10', 'Nlrp12', 'Nlrp1b', 'NLRP3', 'Nlrp4a', 'Nmnat2', 'NMT1', 'NMT2', 'NOD1', 'NOD2', 'NOL6', 'NOL8', 'NOM1', 'Nomo1', 'NOS1', 'NOS2', 'NOS3', 'notch1', 'Notch2', 'NOTCH3', 'NOTCH4', 'NOVA1', 'nova2', 'NPAS2', 'NPAT', 'npc1', 'Npc1l1', 'npepps', 'NPHP3', 'NPHP4', 'NPHS1', 'Nploc4', 'NPNT', 'NPR1', 'NPR3', 'NPTX1', 'nptxr', 'NR1D2', 'NR2C2', 'NR2F2', 'NR3C1', 'NR3C2', 'nr4a3', 'Nr6a1', 'Nrap', 'NRAS', 'nrcam', 'NRD1', 'NRF1', 'NRIP1', 'Nrk', 'NRP1', 'NRP2', 'NRXN1', 'NRXN2', 'Nrxn3', 'NSD1', 'NSMAF', 'nsun2', 'NT5DC3', 'NTN1', 'NTN3', 'NTRK2', 'NUAK1', 'Nucks1', 'NUDT12', 'NUDT5', 'NUFIP2', 'NUMA1', 'Numbl', 'NUP133', 'NUP153', 'NUP160', 'NUP188', 'NUP205', 'Nup210', 'NUP210L', 'Nup214', 'NUP98', 'NUS1', 'nwd1', 'Nynrin', 'Nyx', 'Oas3', 'obscn', 'Obsl1', 'OCEL1', 'Ocln', 'OCRL', 'ODZ1', 'Odz2', 'Odz3', 'ODZ4', 'ofd1', 'Ogdh', 'OGFOD1', 'OGFRL1', 'Ogt', 'OLFM3', 'OLFML2A', 'Olfr1344', 'Olfr1443', 'Olfr613', 'ONECUT2', 'Onecut3', 'OPCML', 'OPHN1', 'OPRD1', 'ORC4L', 'OSBP', 'OSBP2', 'OSBPL11', 'OSBPL3', 'Osbpl6', 'OSBPL8', 'oscp1', 'OSMR', 'Otof', 'OTTMUSG00000016609', 'OTUD5', 'OTUD7B', 'oxr1', 'OXSR1', 'Oxtr', 'P2RX7', 'P2RY4', 'P4HA1', 'PACS1', 'PACS2', 'Pacsin1', 'PADI2', 'Pafah1b1', 'Pak2', 'Pak3', 'Pak6', 'pak7', 'palld', 'PALM2', 'pam', 'pan2', 'PAN3', 'Panc2', 'PANK2', 'PANK3', 'Papd5', 'Papln', 'Papola', 'papolb', 'PAPPA', 'PAQR3', 'PAQR8', 'Pard3', 'PARD3B', 'parg', 'PARM1', 'PARP14', 'PARVA', 'PASK', 'PATL1', 'Pax6', 'PAX7', 'PAX9', 'PAXIP1', 'pbrm1', 'Pbxip1', 'PCDH11X', 'PCDH12', 'PCDH15', 'PCDH17', 'PCDH18', 'Pcdh19', 'PCDH20', 'PCDH7', 'PCDH8', 'Pcdha@', 'PCDHA1', 'PCDHA10', 'PCDHA11', 'PCDHA12', 'PCDHA2', 'Pcdha4', 'PCDHA5', 'PCDHA6', 'PCDHA7', 'PCDHA9', 'PCDHAC1', 'PCDHAC2', 'PCDHB16', 'Pcdhb18', 'Pcdhb19', 'Pcdhb21', 'Pcdhg@', 'PCDHGA11', 'PCDHGA12', 'PCDHGA2', 'PCDHGA3', 'PCDHGA8', 'PCDHGA9', 'PCDHGB4', 'PCDHGC3', 'PCDHGC4', 'PCDHGC5', 'Pcf11', 'Pcgf3', 'PCIF1', 'PCLO', 'pcm1', 'PCMTD1', 'PCNT', 'pcnx', 'PCNXL2', 'PCNXL3', 'PCOLCE2', 'Pcsk2', 'PCSK5', 'Pcsk6', 'Pcx', 'PCYOX1', 'Pcyt1a', 'Pcyt1b', 'Pdcd11', 'PDCD6IP', 'PDE10A', 'PDE11A', 'PDE12', 'pde1a', 'Pde1c', 'Pde2a', 'PDE3A', 'PDE3B', 'Pde4a', 'Pde4b', 'PDE4D', 'PDE4DIP', 'PDE5A', 'Pde6a', 'Pde6b', 'PDE6H', 'PDE7A', 'Pde7b', 'pdgfra', 'Pdgfrb', 'PDIK1L', 'PDK1', 'PDLIM5', 'Pdp1', 'PDPK1', 'PDS5A', 'PDS5B', 'pdx1', 'PDXDC1', 'pdxK', 'Pdzd2', 'PDZD7', 'PDZD8', 'PDZRN3', 'Pear1', 'PEG3', 'PELI2', 'Per1', 'PER2', 'per3', 'Pet2', 'pex1', 'Pex13', 'PEX26', 'PFAS', 'Pfkfb2', 'pfkfb3', 'pfkp', 'PGAP1', 'PGM2L1', 'Pgm3', 'pgpep1', 'Pgr15l', 'PHACTR1', 'PHACTR2', 'Phactr3', 'PHACTR4', 'PHEX', 'PHF12', 'Phf15', 'PHF16', 'phf17', 'Phf2', 'phf20', 'PHF20L1', 'PHF21A', 'PHF3', 'PHF6', 'PHF8', 'PHIP', 'PHKA1', 'PHKA2', 'Phkb', 'PHLDB1', 'PHLDB2', 'PHLPP1', 'PHLPP2', 'PHRF1', 'PHTF1', 'phtf2', 'pi15', 'PI4KA', 'PIAS2', 'PICALM', 'PIGG', 'PIGM', 'Pign', 'pigo', 'PIGU', 'PIK3C2B', 'PIK3CA', 'Pik3cb', 'PIK3CD', 'PIK3CG', 'Pik3r1', 'PIK3R3', 'PIK3R4', 'Pik3r5', 'Pikfyve', 'PIP4K2B', 'PIP5K1C', 'PITPNC1', 'PITPNM1', 'pitpnm2', 'PITPNM3', 'piwil2', 'pja2', 'Pkd1', 'PKD1L2', 'Pkd1l3', 'Pkd2', 'Pkdrej', 'PKHD1', 'PKHD1L1', 'Pkib', 'pkn2', 'Pknox1', 'PKP1', 'PKP4', 'PLA2G4B', 'Pla2g4e', 'PLA2R1', 'PLAG1', 'PLAGL1', 'PLAGL2', 'PLB1', 'PLBD2', 'PLCB1', 'PLCB2', 'plcb3', 'Plcb4', 'PLCD4', 'PLCG1', 'PLCG2', 'PLCH1', 'PLCH2', 'PLCL2', 'PLCXD3', 'PLD1', 'PLD5', 'PLEC', 'plek', 'PLEKHA2', 'PLEKHA5', 'PLEKHA6', 'plekha7', 'PLEKHA8', 'PLEKHG1', 'Plekhg2', 'PLEKHG3', 'Plekhg4', 'PLEKHG5', 'PLEKHH1', 'PLEKHH2', 'Plekhm1', 'Plekhm2', 'PLEKHM3', 'PLIN4', 'PLOD1', 'PLP1', 'PLSCR4', 'PLXDC2', 'PLXNA1', 'Plxna2', 'Plxna3', 'PLXNA4', 'PLXNB1', 'plxnb2', 'PLXNB3', 'plxnc1', 'Plxnd1', 'PMEPA1', 'PML', 'PMPCB', 'PMS2', 'PNMA2', 'PNMAL2', 'PNPLA1', 'pnpla3', 'PNPLA6', 'PNPLA7', 'PODXL', 'Pofut1', 'Pogk', 'POGZ', 'polA1', 'Pole', 'POLG', 'POLK', 'POLQ', 'POLR1A', 'POLR2A', 'POLR2B', 'POLR3A', 'Polr3b', 'Polr3e', 'Polr3f', 'POM121', 'Pom121l2', 'POMT2', 'postn', 'pou2f1', 'POU6F1', 'POU6F2', 'PPARGC1A', 'PPFIA2', 'PPFIA3', 'PPFIBP1', 'PPHLN1', 'PPIG', 'Ppip5k1', 'ppip5k2', 'ppl', 'PPM1E', 'Ppm1f', 'ppm1h', 'Ppm1k', 'Ppm1l', 'PPP1CB', 'PPP1R10', 'PPP1R12B', 'PPP1R13B', 'PPP1R15B', 'PPP1R16B', 'Ppp1r2', 'Ppp1r3a', 'PPP1R3B', 'PPP1R3F', 'Ppp1r7', 'ppp1r9a', 'Ppp1r9b', 'PPP2CA', 'Ppp2r2a', 'PPP2R2C', 'Ppp2r3a', 'Ppp2r3d', 'PPP2R5C', 'PPP2R5E', 'ppp3ca', 'PPP6C', 'PPPDE1', 'pprc1', 'Pptc7', 'PRDM1', 'PRDM11', 'PRDM15', 'Prdm16', 'PRDM2', 'PREB', 'prex1', 'PREX2', 'PRICKLE1', 'prickle2', 'PRKAA1', 'prkaa2', 'PRKAB2', 'PRKACB', 'PRKAR2A', 'Prkca', 'Prkce', 'PRKCI', 'Prkcz', 'Prkd3', 'PRKDC', 'PRKX', 'PRLR', 'Prokr1', 'PROKR2', 'PROM2', 'Prosapip1', 'prpf19', 'prpf39', 'PRPF4', 'Prpf40a', 'prpf40b', 'prpf4b', 'Prpf8', 'PRR12', 'PRRC1', 'PRRG3', 'PRRX1', 'PRSS7', 'PRTG', 'PRUNE2', 'Prx', 'Psd2', 'PSD3', 'Psd4', 'Psg16', 'Psg29', 'Psme4', 'ptch1', 'PTCH2', 'Ptchd2', 'pten', 'Ptger1', 'PTGFR', 'Ptgfrn', 'PTGS2', 'PTK2', 'PTK2B', 'PTK7', 'ptp4a1', 'PTPDC1', 'PTPLAD2', 'ptpn1', 'PTPN11', 'PTPN13', 'PTPN14', 'Ptpn15', 'PTPN21', 'Ptpn23', 'Ptpn3', 'Ptpn9', 'PTPRB', 'Ptprc', 'Ptprd', 'PTPRE', 'ptprf', 'PTPRG', 'PTPRJ', 'PTPRM', 'Ptprn2', 'PTPRO', 'ptprq', 'Ptprs', 'PTPRT', 'ptpru', 'ptprz1', 'PUM1', 'PUM2', 'purA', 'pvrl1', 'PVRL3', 'PXDN', 'PYGO1', 'PZP', 'qk', 'QSER1', 'QSOX2', 'R3hdm1', 'R3HDM2', 'RAB11B', 'RAB11FIP1', 'Rab11fip2', 'RAB11FIP3', 'Rab11fip4', 'RAB11FIP5', 'Rab14', 'RAB22A', 'Rab23', 'RAB27B', 'Rab3gap1', 'RAB3GAP2', 'Rab43', 'RAB44', 'rab6b', 'RAB8A', 'RAB8B', 'RABEP1', 'RABGAP1', 'Rabgap1l', 'RABL3', 'RAD1', 'rad50', 'RAD54L', 'RAD54L2', 'RAG1', 'RAI1', 'RAI14', 'ralgapa1', 'RALGAPA2', 'RALGAPB', 'RALGDS', 'RALGPS1', 'ralgps2', 'RANBP10', 'RANBP17', 'Ranbp2', 'RANBP6', 'Rap1gap', 'RAP1GAP2', 'RAP2A', 'RAP2B', 'RAPGEF1', 'RAPGEF2', 'Rapgef3', 'Rapgef4', 'RAPGEF6', 'RAPGEFL1', 'Rasa1', 'RASA2', 'Rasa3', 'RASAL2', 'Rasef', 'RASGRF1', 'Rasgrf2', 'RASGRP1', 'rasgrp3', 'RASGRP4', 'RASSF2', 'RASSF4', 'rassf8', 'RAVER2', 'RB1', 'RB1CC1', 'RBBP4', 'RBBP6', 'Rbl1', 'RBL2', 'rbm16', 'RBM25', 'RBM27', 'RBM28', 'RBM33', 'RBM47', 'RBM9', 'RBMS1', 'rbms2', 'RBMS3', 'RBP3', 'Rbpj', 'RBPJL', 'RC3H2', 'Rcan3', 'RCOR1', 'rdx', 'RECK', 'Recql5', 'REEP3', 'RELN', 'REPS2', 'RERE', 'REST', 'RET', 'REV1', 'REV3L', 'REX2', 'rexo1', 'rfc1', 'RFFL', 'Rftn1', 'Rfwd2', 'RFWD3', 'RFX1', 'rfx3', 'RFX5', 'RFX7', 'RG9MTD2', 'RGAG4', 'RGL1', 'RGNEF', 'rgp1', 'rgs12', 'RGS17', 'rgs3', 'Rgs5', 'RGS7BP', 'Rgs8', 'RGS9BP', 'RHBDD2', 'Rhbdf1', 'Rhobtb1', 'Rhobtb2', 'RHOBTB3', 'Rhoh', 'RHOJ', 'rhoq', 'RHOT1', 'Rhox1', 'Rhox3b', 'RIC3', 'RICTOR', 'Rif1', 'RIMBP2', 'RIMKLA', 'RIMKLB', 'RIMS1', 'RIMS2', 'RIMS3', 'rims4', 'RIN1', 'RIN2', 'rlf', 'RLIM', 'Rmnd5a', 'RNASEL', 'rnasen', 'RNF11', 'Rnf111', 'RNF123', 'RNF139', 'Rnf144a', 'Rnf144b', 'RNF145', 'RNF146', 'Rnf150', 'RNF152', 'RNF157', 'RNF160', 'Rnf165', 'RNF169', 'RNF17', 'RNF19A', 'rnf20', 'RNF213', 'rnf216', 'RNF24', 'Rnf38', 'RNF40', 'Rnf41', 'RNF43', 'RNF44', 'Rnf8', 'RNFT1', 'Rnft2', 'Rngtt', 'rnmt', 'ROBO1', 'robo2', 'robo3', 'ROCK1', 'ROCK2', 'ROD1', 'ROR1', 'RORA', 'RORB', 'ROS1', 'Rp1', 'RP1L1', 'Rp2h', 'RPA1', 'rpap1', 'RPGRIP1', 'RPGRIP1L', 'RPH3A', 'Rprd1a', 'RPRD1B', 'Rprd2', 'RPS6KA2', 'Rps6ka3', 'RPS6KA5', 'RPS6KA6', 'Rptn', 'RPTOR', 'Rragd', 'RRBP1', 'Rreb1', 'Rrm2b', 'Rrp12', 'Rrp1b', 'RRP7A', 'RS1', 'RSAD1', 'RSAD2', 'RSBN1', 'RSBN1L', 'RSC1A1', 'rsf1', 'RSRC2', 'RTEL1', 'Rtf1', 'RTKN2', 'RTN3', 'rtn4', 'RTTN', 'RUFY2', 'Runx1', 'RUNX1T1', 'Runx2', 'RUSC2', 'RXFP1', 'rxfp3', 'RXRA', 'RYBP', 'RYR1', 'RYR2', 'RYR3', 'S100PBP', 'S1pr3', 'saal1', 'SACS', 'SAG', 'sall1', 'SALL2', 'sall3', 'sall4', 'Samd4', 'SAMD4B', 'SAMD8', 'SAMD9L', 'SAP130', 'SAPS3', 'SARDH', 'Sarm1', 'sart1', 'SASH1', 'Sass6', 'Satb1', 'SATB2', 'SBF1', 'Sbf2', 'Sbk1', 'Sbno1', 'SBNO2', 'SCAF1', 'Scai', 'SCAMP1', 'SCAP', 'scaper', 'Scd1', 'Scd2', 'SCML2', 'SCML4', 'Scn10a', 'SCN11A', 'Scn1a', 'Scn2a1', 'SCN3A', 'SCN3B', 'SCN4A', 'SCN4B', 'SCN5A', 'SCN7A', 'SCN8A', 'SCN9A', 'scrib', 'SCRN1', 'SCRT1', 'Scube1', 'SCUBE3', 'Scyl3', 'SDAD1', 'sdc3', 'SDF4', 'SDK1', 'sdk2', 'sec14l1', 'SEC16A', 'SEC16B', 'sec22c', 'SEC23A', 'SEC23IP', 'Sec24a', 'SEC24B', 'SEC24C', 'SEC31A', 'SEC63', 'secisbp2l', 'SEL1L', 'Sel1l3', 'SEMA3A', 'SEMA3C', 'Sema3d', 'Sema3e', 'SEMA3G', 'SEMA4D', 'SEMA4F', 'SEMA4G', 'SEMA5A', 'sema5b', 'SEMA6A', 'Sema6c', 'SEMA6D', 'SENP1', 'SENP2', 'SENP6', 'Senp7', 'sephs1', 'sepsecs', 'sept11', 'SEPT14', 'SEPT3', 'Sept6', 'SEPT8', 'SERAC1', 'serbp1', 'SERINC5', 'Serpini1', 'SERTAD2', 'Sestd1', 'SETBP1', 'SETD1A', 'SETD1B', 'Setd5', 'setd7', 'Setdb1', 'Setx', 'SEZ6', 'SEZ6L', 'Sf1', 'Sf3a1', 'sf3b1', 'sf3b3', 'SFI1', 'SFMBT1', 'SFMBT2', 'Sfrp1', 'SFRS1', 'SFRS14', 'Sfrs15', 'Sfrs17b', 'sfrs18', 'SFRS2IP', 'SFT2D2', 'Sgip1', 'sgms1', 'sgms2', 'SGOL2', 'sgpl1', 'SGSH', 'SGSM1', 'SGSM2', 'SH2B3', 'SH2D2A', 'SH3BP4', 'SH3D19', 'Sh3kbp1', 'SH3PXD2A', 'SH3PXD2B', 'sh3rf1', 'sh3rf2', 'SH3TC1', 'SH3TC2', 'SHANK1', 'SHANK2', 'SHANK3', 'Shc4', 'She', 'SHISA6', 'shisa7', 'shisa9', 'SHOC2', 'SHPRH', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'Shroom4', 'SIDT1', 'SIDT2', 'SIGLEC1', 'Sik1', 'SIK3', 'SIM1', 'SIN3A', 'Sin3b', 'SIPA1L1', 'SIPA1L2', 'sirt1', 'Sis', 'six1', 'Six4', 'SIX6', 'ski', 'Skil', 'Skp2', 'SLAIN2', 'SLC11A2', 'SLC12A1', 'SLC12A2', 'SLC12A5', 'SLC12A6', 'SLC12A7', 'SLC15A3', 'SLC16A1', 'SLC16A10', 'SLC16A12', 'SLC16A2', 'SLC16A6', 'SLC17A2', 'SLC17A3', 'SLC17A6', 'SLC17A8', 'Slc1a2', 'SLC1A3', 'SLC22A23', 'SLC22A6', 'SLC23A2', 'SLC24A1', 'SLC24A2', 'SLC24A3', 'SLC25A36', 'Slc25a37', 'SLC25A46', 'Slc26a2', 'SLC26A7', 'SLC26A9', 'SLC27A4', 'SLC28A3', 'Slc29a3', 'SLC2A12', 'SLC2A13', 'slc30a1', 'slc30a10', 'SLC30A4', 'SLC30A7', 'SLC31A1', 'slc34a2', 'SLC35A3', 'SLC35A5', 'SLC35B4', 'SLC35D1', 'SLC35E2', 'slc35f1', 'SLC36A1', 'SLC36A4', 'SLC37A2', 'SLC37A3', 'SLC38A1', 'SLC38A10', 'slc38a2', 'slc38a9', 'SLC39A10', 'slc39a14', 'SLC39A8', 'slc39a9', 'SLC41A1', 'SLC41A2', 'Slc43a2', 'slc44a1', 'SLC44A5', 'SLC45A4', 'SLC4A1', 'SLC4A10', 'SLC4A2', 'SLC4A3', 'SLC4A4', 'SLC4A5', 'Slc4a7', 'SLC4A8', 'SLC5A1', 'SLC5A3', 'SLC5A7', 'SLC5A8', 'Slc5a9', 'SLC6A1', 'SLC6A11', 'SLC6A17', 'SLC6A18', 'SLC6A19', 'SLC6A2', 'SLC6A5', 'slc6a6', 'SLC7A1', 'SLC7A11', 'SLC7A14', 'SLC7A2', 'SLC7A8', 'SLC8A3', 'slc9a1', 'SLC9A5', 'Slc9a6', 'SLC9A8', 'SLCO2A1', 'SLCO2B1', 'SLCO3A1', 'slco5a1', 'Slfn5', 'Slfn8', 'Slfn9', 'slit1', 'Slit2', 'SLIT3', 'SLITRK1', 'SLITRK4', 'SLITRK5', 'Slitrk6', 'SLK', 'SLMAP', 'SMAD3', 'Smad5', 'SMAD7', 'SMAD9', 'SMARCA1', 'smarca2', 'SMARCA4', 'SMARCA5', 'Smarcad1', 'SMARCAL1', 'SMARCB1', 'SMARCC1', 'smarcc2', 'Smc1a', 'SMC1B', 'SMC2', 'smc3', 'SMC4', 'smc5', 'SMC6', 'SMCHD1', 'SMCR7L', 'Smcr8', 'smek1', 'smek2', 'smg1', 'smg5', 'SMG6', 'SMG7', 'smpd2', 'Smpd3', 'SMPD4', 'SMURF1', 'SMURF2', 'Smyd3', 'SNAP25', 'SNAP91', 'SNAPC4', 'SNED1', 'Snhg11', 'SNPH', 'Snrk', 'SNRNP200', 'Snrpc', 'sntb2', 'SNTG1', 'Snx13', 'Snx18', 'SNX27', 'Snx30', 'SNX33', 'soat1', 'Sobp', 'socs2', 'Socs4', 'SOCS5', 'SOCS6', 'socs7', 'SON', 'SORBS2', 'SORCS1', 'Sorcs2', 'SORCS3', 'SORL1', 'Sort1', 'SOS1', 'SOS2', 'SOX12', 'Sox19', 'Sox4', 'SOX5', 'Sox6', 'SOX9', 'sp1', 'SP3', 'SP4', 'SP8', 'Sp9', 'Spag9', 'SPAST', 'SPATA13', 'SPATA2', 'SPEF2', 'speG', 'spen', 'SPG11', 'SPHK1', 'SPHKAP', 'spin1', 'Spin3', 'Spin4', 'SPINK5', 'SPIRE1', 'Spna1', 'Spna2', 'Spnb1', 'Spnb2', 'Spnb3', 'Spnb4', 'SPOCK2', 'SPON1', 'SPOPL', 'spred1', 'SPRED2', 'SPRED3', 'SPRY3', 'SPRY4', 'Sptlc2', 'SPTY2D1', 'Srcap', 'SRCIN1', 'Srebf1', 'srebf2', 'Srf', 'SRGAP1', 'srgap2', 'SRGAP3', 'srl', 'srpk1', 'SRPK2', 'SRRM1', 'SRRM2', 'srrm4', 'SS18', 'SS18L1', 'SSFA2', 'SSH1', 'SSH2', 'sspN', 'sspO', 'SSR1', 'ST18', 'ST3GAL1', 'ST3GAL2', 'ST3GAL6', 'ST5', 'ST6GAL1', 'ST6GAL2', 'ST6GALNAC5', 'ST7L', 'ST8SIA1', 'ST8SIA2', 'ST8SIA3', 'ST8SIA4', 'ST8SIA6', 'STAB1', 'STAB2', 'STAG1', 'STAG2', 'STAG3', 'STAM2', 'STAMBPL1', 'STARD13', 'STARD4', 'Stard8', 'STARD9', 'STAT1', 'STAT2', 'Stat3', 'Stat5b', 'Steap2', 'STEAP3', 'STIL', 'Stim2', 'STK10', 'Stk11ip', 'STK36', 'stk38l', 'STK4', 'STON2', 'STOX2', 'STRBP', 'STRC', 'strN', 'STRN3', 'STRN4', 'stt3b', 'STX17', 'STX1B', 'Stx8', 'Stxbp4', 'STXBP5', 'STXBP5L', 'Stxbp6', 'STYX', 'sufU', 'Sulf1', 'SULF2', 'SUN1', 'SUPT16H', 'SUPT6H', 'suv39h2', 'suv420h1', 'SUZ12', 'SV2A', 'SV2B', 'SVEP1', 'SVIL', 'SWAP70', 'sycp2', 'SYDE2', 'Sykb', 'SYMPK', 'SYN3', 'SYNCRIP', 'SYNE1', 'SYNE2', 'SYNGR1', 'SYNJ1', 'SYNJ2', 'SYNM', 'SYNPO', 'SYNPO2', 'SYNPO2L', 'SYNPR', 'SYNRG', 'Sypl', 'SYPL2', 'sys1', 'Syt1', 'Syt11', 'Syt2', 'SYT6', 'Syt7', 'SYTL2', 'SYVN1', 'Szt2', 'tab2', 'TAB3', 'TACC1', 'TACC2', 'Tacr1', 'TAF1', 'taf10', 'Taf2', 'taf3', 'Taf4a', 'tal1', 'TANC1', 'TANC2', 'Taok1', 'taok2', 'TAOK3', 'TARDBP', 'tas1r1', 'TBC1D1', 'TBC1D12', 'TBC1D14', 'TBC1D16', 'TBC1D2', 'Tbc1d24', 'TBC1D2B', 'tbc1d30', 'TBC1D5', 'TBC1D8', 'TBC1D8B', 'TBC1D9', 'tbc1d9b', 'tbcel', 'tbl1x', 'TBL1XR1', 'tbl3', 'tbx18', 'TBX20', 'tbx3', 'tceanc', 'Tceb3', 'tcerg1', 'TCF12', 'TCF20', 'tcf23', 'TCF4', 'Tcf7l2', 'Tcfap2b', 'Tcfcp2l1', 'Tcfe3', 'TCHH', 'TCP11L1', 'Tcte2', 'TDRD1', 'tdrd6', 'TDRD9', 'TEAD1', 'Tecpr1', 'TECTA', 'tef', 'TEK', 'TENC1', 'TEP1', 'Tet1', 'TET2', 'TET3', 'Tex14', 'TEX15', 'tex2', 'Tex9', 'TFAM', 'tfdp2', 'TFRC', 'Tg', 'TGFA', 'TGFB2', 'tgfbr1', 'Tgfbr2', 'Tgfbr3', 'TGFBRAP1', 'TGM3', 'Tgoln1', 'TGOLN2', 'tgs1', 'Thada', 'Thbs1', 'Thbs3', 'THEMIS', 'Thoc2', 'THRAP3', 'thrB', 'THSD1', 'THSD4', 'THSD7A', 'THSD7B', 'Tia1', 'Tial1', 'tiam1', 'TIAM2', 'Tie1', 'TIMELESS', 'TIMP3', 'TJP1', 'TJP2', 'TLE1', 'TLE2', 'TLE3', 'Tle4', 'Tlk1', 'TLK2', 'TLL1', 'TLN1', 'TLN2', 'tlr13', 'TLR3', 'TLR4', 'TM9SF3', 'TMBIM1', 'TMC1', 'TMC3', 'TMC7', 'Tmcc1', 'Tmcc2', 'tmcc3', 'tmco1', 'TMCO3', 'Tmco7', 'TMED8', 'TMEM104', 'tmem106b', 'Tmem127', 'TMEM131', 'TMEM132B', 'TMEM132C', 'Tmem132d', 'TMEM132E', 'TMEM151A', 'TMEM164', 'TMEM168', 'TMEM184A', 'Tmem194', 'Tmem2', 'TMEM200A', 'TMEM201', 'Tmem214', 'TMEM215', 'TMEM229A', 'tmem26', 'tmem33', 'TMEM44', 'TMEM47', 'Tmem48', 'TMEM5', 'TMEM56', 'tmem57', 'TMEM64', 'TMEM87A', 'tmem87b', 'TMEM8B', 'tmod2', 'TMPRSS11BNL', 'TMTC1', 'TMTC2', 'TMTC3', 'TMX3', 'TMX4', 'TNC', 'TNFAIP3', 'Tnfrsf11a', 'TNFRSF19', 'Tnfrsf22', 'TNFSF10', 'TNK2', 'tnks', 'TNKS1BP1', 'tnks2', 'TNN', 'tnpo1', 'Tnpo2', 'TNPO3', 'TNR', 'TNRC18', 'TNRC6A', 'TNRC6B', 'TNRC6C', 'TNS3', 'TNS4', 'TNXB', 'TOM1L1', 'TOM1L2', 'TOP2A', 'TOP2B', 'TOPBP1', 'TOR1AIP1', 'TOR1AIP2', 'Tor1b', 'Tox4', 'Tpcn1', 'tpp2', 'TPPP', 'tpr', 'traF2', 'TRAF3', 'TRAF6', 'TRAK1', 'TRAK2', 'Tram2', 'TRANK1', 'TRAPPC10', 'TRAPPC9', 'TREM2', 'TRERF1', 'TRHDE', 'TRIB2', 'trim2', 'trim23', 'TRIM24', 'TRIM25', 'TRIM26', 'trim33', 'TRIM36', 'TRIM37', 'TRIM44', 'TRIM45', 'TRIM56', 'Trim65', 'TRIM66', 'TRIM67', 'TRIM71', 'TRIM9', 'TRIO', 'Triobp', 'Trip11', 'Trip12', 'TRIP4', 'TRMT12', 'TRMT2B', 'Tro', 'Trove2', 'Trp53bp1', 'Trp53bp2', 'Trp53inp1', 'Trp53rk', 'Trp63', 'Trp73', 'trpa1', 'TRPC2', 'TRPC5', 'TRPM1', 'TRPM2', 'TRPM3', 'TRPM4', 'Trpm5', 'TRPM6', 'TRPM7', 'TRPM8', 'TRPS1', 'TRPV3', 'TRRAP', 'TRUB1', 'TSC1', 'TSC2', 'Tsc22d1', 'TSC22D2', 'Tsen2', 'TSHR', 'TSHZ1', 'TSHZ2', 'TSHZ3', 'TSPAN11', 'Tspan2', 'TSPAN7', 'TSPYL5', 'TTBK1', 'TTBK2', 'Ttc13', 'ttc14', 'ttc17', 'ttc21a', 'ttc21b', 'TTC25', 'TTC28', 'TTC3', 'TTC37', 'Ttc39b', 'Ttc7', 'Ttf1', 'TTF2', 'TTL', 'TTLL4', 'TTLL5', 'ttll7', 'TTN', 'ttpal', 'ttyh1', 'Ttyh3', 'TUB', 'TUBGCP4', 'tubgcp6', 'TUG1', 'TULP4', 'Twsg1', 'TXLNA', 'Txlnb', 'TXNDC11', 'Txndc16', 'Tyk2', 'tyro3', 'uaca', 'Uba1', 'UBA6', 'UBAP2', 'Ubap2l', 'UBASH3B', 'Ube2h', 'UBE2K', 'UBE2O', 'ube2v2', 'Ube3a', 'UBE3B', 'UBE3C', 'ube4a', 'Ube4b', 'ubfd1', 'Ubn1', 'ubn2', 'UBR1', 'UBR2', 'UBR3', 'UBR4', 'UBR5', 'UBXN2B', 'UBXN4', 'UBXN7', 'Uevld', 'UGCG', 'UGGT1', 'Uggt2', 'Ugt2b34', 'Uhmk1', 'UHRF1BP1', 'UHRF1BP1L', 'ULK1', 'Ulk2', 'Ulk4', 'Umodl1', 'UNC13A', 'Unc13b', 'Unc13c', 'UNC13D', 'UNC45B', 'Unc5b', 'UNC5C', 'Unc79', 'UNC80', 'Unkl', 'UPF1', 'Upf2', 'Upf3b', 'URB1', 'URB2', 'URGCP', 'USH2A', 'USP10', 'USP12', 'USP13', 'usp19', 'USP20', 'USP22', 'usp24', 'usp25', 'USP26', 'usp28', 'usp3', 'USP31', 'USP32', 'usp33', 'USP34', 'usp36', 'USP37', 'USP38', 'usp40', 'Usp42', 'Usp43', 'USP45', 'USP46', 'usp47', 'usp48', 'usp49', 'usp53', 'USP54', 'USP6NL', 'usp8', 'USP9X', 'USP9Y', 'USPL1', 'Ust', 'UTP15', 'Utp20', 'utrn', 'UTY', 'UXS1', 'Vangl2', 'vapB', 'vars', 'Vars2', 'Vash2', 'VAT1L', 'VAV2', 'VAV3', 'VCAN', 'vcl', 'Vcpip1', 'VDR', 'VEPH1', 'vezf1', 'VEZT', 'VIPR1', 'Vma21', 'vopp1', 'Vprbp', 'VPS13A', 'VPS13B', 'Vps13c', 'VPS13D', 'VPS18', 'VPS26B', 'vps33a', 'vps37a', 'VPS39', 'VPS4B', 'VPS54', 'vps8', 'VSIG10', 'vwa2', 'VWA3A', 'vwa5a', 'Vwa5b1', 'VWA5B2', 'VWC2L', 'vwf', 'WAPAL', 'WARS2', 'Wasf2', 'WASL', 'Wbp7', 'Wbscr17', 'Wdfy1', 'WDFY3', 'WDFY4', 'Wdhd1', 'WDR11', 'WDR13', 'WDR17', 'WDR19', 'Wdr25', 'Wdr3', 'WDR33', 'WDR35', 'WDR37', 'WDR44', 'WDR47', 'Wdr52', 'Wdr59', 'Wdr6', 'Wdr62', 'WDR65', 'WDR66', 'WDR7', 'Wdr72', 'WDR76', 'WDR77', 'WDR78', 'WDR81', 'wdr82', 'WDR86', 'WDR90', 'WDTC1', 'Whrn', 'Whsc1', 'WHSC1L1', 'WIPF1', 'WIPF2', 'wipi2', 'WISP1', 'WIZ', 'WNK1', 'WNK2', 'WNK3', 'WNK4', 'WNT4', 'WNT5A', 'WNT9B', 'wrn', 'WSCD2', 'WWC1', 'Wwc2', 'WWP1', 'WWTR1', 'X99384', 'XDH', 'XIAP', 'XIRP1', 'xirp2', 'XK', 'xkr8', 'XPC', 'xpo1', 'XPO6', 'XPO7', 'Xpot', 'XPR1', 'XRN1', 'YAP1', 'YEATS2', 'YES1', 'YIPF6', 'YLPM1', 'YME1L1', 'YPEL2', 'YTHDC1', 'Ythdc2', 'YTHDF3', 'YY2', 'zan', 'ZBED4', 'ZBTB1', 'zbtb11', 'ZBTB16', 'ZBTB20', 'Zbtb26', 'ZBTB33', 'ZBTB34', 'ZBTB38', 'ZBTB39', 'ZBTB4', 'Zbtb40', 'Zbtb41', 'ZBTB43', 'zbtb44', 'zbtb46', 'ZBTB5', 'Zbtb6', 'ZBTB7A', 'ZBTB7B', 'ZBTB7C', 'ZBTB8B', 'Zc3h11a', 'Zc3h12b', 'Zc3h12d', 'ZC3H13', 'Zc3h4', 'ZC3H6', 'ZC3H7B', 'ZC3HAV1', 'ZC3HAV1L', 'ZCCHC11', 'ZCCHC14', 'Zcchc18', 'ZCCHC2', 'zcchc24', 'ZCCHC4', 'ZCCHC6', 'ZCCHC7', 'zcchc8', 'ZDBF2', 'ZDHHC15', 'ZDHHC17', 'ZDHHC18', 'Zdhhc20', 'ZDHHC21', 'ZDHHC3', 'Zdhhc5', 'Zdhhc8', 'ZEB1', 'Zeb2', 'zer1', 'Zfand1', 'ZFAND5', 'Zfc3h1', 'Zfhx2', 'ZFHX3', 'ZFHX4', 'Zfml', 'ZFP106', 'Zfp113', 'ZFP12', 'Zfp120', 'Zfp142', 'Zfp148', 'Zfp157', 'Zfp174', 'Zfp191', 'Zfp192', 'Zfp217', 'Zfp229', 'Zfp238', 'Zfp26', 'Zfp260', 'Zfp275', 'Zfp276', 'ZFP28', 'Zfp280b', 'Zfp280c', 'Zfp280d', 'Zfp281', 'Zfp282', 'Zfp292', 'Zfp295', 'Zfp300', 'Zfp318', 'Zfp319', 'Zfp322a', 'Zfp324', 'Zfp326', 'Zfp329', 'Zfp334', 'Zfp335', 'Zfp354a', 'Zfp354c', 'Zfp365', 'Zfp395', 'Zfp397', 'Zfp40', 'Zfp407', 'Zfp423', 'Zfp451', 'Zfp456', 'Zfp46', 'Zfp462', 'Zfp473', 'Zfp503', 'Zfp507', 'Zfp516', 'Zfp518', 'Zfp518b', 'Zfp521', 'Zfp532', 'Zfp536', 'Zfp541', 'Zfp560', 'Zfp568', 'Zfp574', 'Zfp592', 'Zfp597', 'Zfp60', 'Zfp600', 'Zfp606', 'Zfp608', 'Zfp616', 'Zfp618', 'Zfp619', 'ZFP62', 'Zfp629', 'Zfp641', 'Zfp644', 'Zfp646', 'Zfp652', 'Zfp654', 'Zfp661', 'Zfp664', 'Zfp687', 'Zfp697', 'Zfp704', 'Zfp706', 'Zfp710', 'Zfp711', 'Zfp715', 'Zfp719', 'Zfp71-rs1', 'Zfp738', 'Zfp740', 'Zfp760', 'Zfp770', 'Zfp775', 'Zfp790', 'Zfp800', 'Zfp804a', 'Zfp81', 'Zfp827', 'Zfp828', 'Zfp831', 'Zfp84', 'ZFP91', 'Zfp91-Cntf', 'ZFP92', 'Zfp93', 'Zfpm2', 'zfr', 'ZFR2', 'ZFX', 'ZFYVE16', 'ZFYVE20', 'ZFYVE26', 'ZFYVE28', 'zfyve9', 'zgpat', 'ZHX1', 'ZHX2', 'zhx3', 'ZIC1', 'Zic5', 'Zik1', 'ZKSCAN1', 'Zkscan17', 'Zkscan3', 'zmat3', 'ZMAT4', 'Zmiz1', 'Zmiz2', 'zmym2', 'ZMYM3', 'Zmym4', 'ZMYM5', 'ZMYM6', 'zmynd11', 'Zmynd8', 'ZNF512B', 'Znfx1', 'znrf2', 'ZNRF3', 'ZPBP', 'ZRANB1', 'ZRANB3', 'ZRSR1', 'Zscan20', 'ZSCAN29', 'zswim4', 'ZSWIM5', 'Zswim6', 'ZXDB', 'ZXDC', 'ZYG11B', 'ZZEF1', 'ZZZ3');
// These genes have coverage >= 100 or show the distinctive pattern.
//$aGenesToAnalyze = array('Col3a1', 'Thbs1', 'Dclk1', 'Col1a2', 'Serbp1', 'Fbln2', 'Gnb2l1', 'Ppp2ca', 'Myo1c', 'Col1a1', 'Vcl', 'Slc5a3', 'Pdgfrb');
/*
$aGenesToAnalyze = array_map('strtolower', $aGenesToAnalyze);
$aPositionsToAnalyze = array();
foreach ($aPositionsPerGene as $sGene => $aGene) {
    if (in_array(strtolower($sGene), $aGenesToAnalyze)) {
        $aPositionsToAnalyze[$sGene] = $aGene;
    }
}
*/
//foreach ($aPositionsToAnalyze as $sGene => $aGene) {
foreach ($aPositionsPerGene as $sGene => $aGene) {
    $nOriCandidates = count($aGene['positions']);
    $nCandidates = count($aGene['candidates']);
    $aTranscripts = array_keys($aGene['transcripts']);
    print($sGene . "\t" . 'Positions analyzed:' . "\t" . $nOriCandidates . "\t" . 'Candidate peaks:' . "\t" . $nCandidates . "\n");
    print('G_Position' . "\t" . 'Coverage' . "\t" . implode("\t", $aTranscripts) . "\n");
    foreach ($aGene['candidates'] as $nPosition) {
        $aPosition = $aGene['positions'][$nPosition];
        print('chr' . $aGene['chr'] . ':' . $nPosition . "\t" . $aPosition['coverage']);
        foreach ($aTranscripts as $sTranscript) {
            if (isset($aPosition['mappings'][$sTranscript])) {
                print("\t" . $aPosition['mappings'][$sTranscript]);
            } else {
                print("\t-");
            }
        }
        print("\n");
    }
    print("\n");
}

print("\n");
exit(0);
?>
