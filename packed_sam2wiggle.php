#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * PackedSam2Wiggle converts packed sam files to wiggle files, mapping the reads
 * with transcriptome alignments to the genome using the transcript position
 * file produced by mm10_transcript_positions_create.php. This script creates 4
 * Wiggle files per SAM file; F unfiltered, F filtered (NR and XR removed), R
 * unfiltered, and R filtered.
 *
 * Created     : 2013-08-21
 * Modified    : 2013-09-19
 * Version     : 0.4
 *
 * Copyright   : 2013-2015 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
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
        'version' => '0.4',
        'suffix' => '.wig5',
    );

echo 'PackedSam2Wiggle v.' . $_SETT['version'] . "\n" .
     'PLEASE DO NOT USE THIS SCRIPT ON THE I DRIVE; IT\'S INCREDIBLY SLOW THERE.' . "\n\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 2) {
    die('Usage: ' . $sScriptName . ' TRANSCRIPT_POSITION_FILE PACKED_SAM_FILE1 [PACKED_SAM_FILE2 [...]]' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}





// Prepare transcript locations file, read into memory.
$sTranscriptPositionsFile = array_shift($aFiles);
$aTranscriptPositionsFile = file($sTranscriptPositionsFile, FILE_IGNORE_NEW_LINES);
$aTranscripts = array(); // transcript => array(chr, strand, array(array(exon_start, exon_end)))
$nTranscripts = 0;
print('Parsing Transcript locations file... ');
foreach ($aTranscriptPositionsFile as $nLine => $sLine) {
    $nLine ++;
    if (!trim($sLine) || $sLine{0} == '#') {
        continue;
    }
    if (preg_match('/^([NX][MR]_\d+\.\d+)\t(\d{1,2}|[XY])\t([+-])\t([\[\]0-9,]+)$/', $sLine, $aRegs)) {
        // Valid transcript position found.
        list(,$sTranscript, $sChr, $sStrand, $sExonPositions) = $aRegs;
        if (!($aExonPositions = json_decode($sExonPositions))) {
            die("\n" .
                'Can\'t parse line ' . $nLine . ' in file ' . $sTranscriptPositionsFile . '.' . "\n\n");
        }
        $aTranscripts[$sTranscript] = array($sChr, $sStrand, $aExonPositions);
        $nTranscripts ++;
    }
}
unset($aTranscriptPositionsFile);
// If the wrong file has been passed, we have no valid transcripts. Then it will
// make no sense at all to continue.
if (!count($aTranscripts)) {
    die("\n" .
        'Didn\'t find any valid transcript positions. Make sure you passed the correct transcript location file as the first argument.' . "\n\n");
}
print('done, loaded ' . $nTranscripts . ' transcript positions in memory.' . "\n");





// Mapping positions on several transcripts to the genome may result in multiple
// positions mapping to the same genomic location. So first we map everything,
// after that we'll sort and write the file.
foreach ($aFiles as $sFile) {
    $nLine = 0;
    // To save memory, we'll read the packed SAM files line by line (usually, they're about 6-14 MB).
    $fIn = @fopen($sFile, 'r');
    if (!$fIn) {
        die('Unable to open ' . $sFile . '.' . "\n\n");
    }
    $aFileNamesOut =
        array(
            '+' => array(
                $sFile . '.F' . $_SETT['suffix'],
                $sFile . '.F.filtered' . $_SETT['suffix']),
            '-' => array(
                $sFile . '.R' . $_SETT['suffix'],
                $sFile . '.R.filtered' . $_SETT['suffix']));
    $aFilesOut = array();
    // Data, per strand two arrays (non-filtered, filtered).
    $aData = array('+' => array(), '-' => array()); // strand => array(chromosome => array(position => array(coverage unfiltered, coverage filtered), ...), ...);
    foreach ($aFileNamesOut as $sStrand => $aStrandFiles) {
        $aFilesOut[$sStrand] = array();
        foreach ($aStrandFiles as $sStrandFile) {
            $f = @fopen($sStrandFile, 'w');
            if (!$f) {
                die('Unable to open file for writing: ' . $aFileNamesOut[$sStrand] . '.' . "\n\n");
            }
            $aFilesOut[$sStrand][] = $f;
        }
    }

    $nUnmapped = 0;
    print("\n" .
          'Reading ' . $sFile . '...' . "\n");
    while ($sLine = fgets($fIn)) {
        $nLine ++;
        $sLine = rtrim($sLine);
        list($sTranscript, $nPosition, $nCoverage) = explode("\t", $sLine);
        // Calculate the genomic position.
        // FIXME: Is this a smart idea?
        /*
        if (!isset($aTranscripts[$sTranscript])) {
            // Not found. This might be, because we have a newer, or older version. Try and find it?
            list($sTranscriptPrefix,$nVersion) = explode('.', $sTranscript, 2);
            for ($i = $nVersion; $i > 0; $i --) {
                if (isset($aTranscripts[$sTranscriptPrefix . '.' . $i])) {
                    $sTranscript = $sTranscriptPrefix . '.' . $i;
                    break;
                }
            }
        }
        // Has no effect.
        /*
        if (!isset($aTranscripts[$sTranscript])) {
            // Not found. This might be, because we have a newer, or older version. Try and find it?
            list($sTranscriptPrefix,$nVersion) = explode('.', $sTranscript, 2);
            for ($i = $nVersion; $i < $nVersion + 10; $i ++) {
                if (isset($aTranscripts[$sTranscriptPrefix . '.' . $i])) {
                    $sTranscript = $sTranscriptPrefix . '.' . $i;
                    break;
                }
            }
        }
        */
        if (isset($aTranscripts[$sTranscript])) {
            // Found!
            list($sChr, $sStrand, $aExonPositions) = $aTranscripts[$sTranscript];
            $sChr = 'chr' . $sChr;
            // Calculate real position on genome.
            if ($sStrand == '+') {
                $nOffset = $nPosition;
                foreach ($aExonPositions as $aExon) {
                    $nLength = $aExon[1] - $aExon[0] + 1;
                    if ($nOffset > $nLength) {
                        $nOffset -= $nLength;
                        continue;
                    } else {
                        $nPosition = $aExon[0] + $nOffset - 1;
                        break;
                    }
                }
            } else {
                $nOffset = $nPosition;
                for ($i = count($aExonPositions) - 1; $i >= 0; $i --) {
                    $aExon = $aExonPositions[$i];
                    $nLength = $aExon[1] - $aExon[0] + 1;
                    if ($nOffset > $nLength) {
                        $nOffset -= $nLength;
                        continue;
                    } else {
                        $nPosition = $aExon[1] - $nOffset + 1;
                        break;
                    }
                }
            }
            if (!isset($aData[$sStrand][$sChr])) {
                $aData[$sStrand][$sChr] = array();
            }
            if (!isset($aData[$sStrand][$sChr][$nPosition])) {
                $aData[$sStrand][$sChr][$nPosition] = array(0, 0);
            }
            $aData[$sStrand][$sChr][$nPosition][0] += $nCoverage; // Unfiltered, always count.
            if (preg_match('/^.M_/', $sTranscript)) {
                $aData[$sStrand][$sChr][$nPosition][1] += $nCoverage; // NM or XM are counted, rest (NR, XR) is not.
            }
        } else {
            // We chose to ignore the ones that we cannot find mappings of. Just count.
            $nUnmapped ++;
//            print('Missing position/strand information on transcript ' . $sTranscript . ' on line ' . $nLine . ' in file ' . $sFile . '.' . "\n");
//            die('Missing position/strand information on transcript ' . $sTranscript . ' on line ' . $nLine . ' in file ' . $sFile . '.' . "\n\n");
        }
    }
    fclose($fIn);

    print('Done reading ' . $nLine . ' lines (' . $nUnmapped . ' unmappable positions), writing output... ');

    // Write output files.
    $nLines = 0;
    foreach ($aFilesOut as $sStrand => $aStrandFiles) {
        // Sort the results on chromosome.
        ksort($aData[$sStrand]);
        // Write track headers in both files.
        fputs($aStrandFiles[0], 'track type=wiggle_0 name=' . $sFile . ' description=' . $sFile . ' visibility=full' . "\n");
        fputs($aStrandFiles[1], 'track type=wiggle_0 name=' . $sFile . ' description=' . $sFile . ' visibility=full' . "\n");
        $nLines += 2;
        foreach ($aData[$sStrand] as $sChr => $aPositions) {
            fputs($aStrandFiles[0], 'variableStep chrom=' . $sChr . "\n");
            fputs($aStrandFiles[1], 'variableStep chrom=' . $sChr . "\n");
            $nLines += 2;
            // Sort positions based on their numbers.
            ksort($aPositions, SORT_NUMERIC);
            foreach ($aPositions as $nPosition => $aCoverage) {
                fputs($aStrandFiles[0], $nPosition . "\t" . $aCoverage[0] . "\n");
                $nLines ++;
                // Store in filtered file only if we have coverage...
                if ($aCoverage[1]) {
                    fputs($aStrandFiles[1], $nPosition . "\t" . $aCoverage[1] . "\n");
                    $nLines ++;
                }
            }
        }
        fclose($aStrandFiles[0]);
        fclose($aStrandFiles[1]);
        print(($sStrand == '+'? 'F' : 'R') . $_SETT['suffix'] . ' done... ');
    }


    print('Done, wrote ' . $nLines . ' lines in total.' . "\n");
}
die('All files done.' . "\n");
?>
