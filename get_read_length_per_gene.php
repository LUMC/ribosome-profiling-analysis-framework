#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * GetReadLengthPerGene takes Gene IDs and determines the read length
 * distribution for its transcripts per sample, from both transcriptome and
 * genome alignment SAM files. It requires the mm10_gene_list.txt file to see
 * which transcripts belongs to which gene, the mm10_transcript_positions.txt
 * file to check the strand and the positions of the transcripts, and the SAM
 * files.
 * The script creates one file per transcript, the read lengths in the first
 * column, and the number of reads with this length per sample in the next
 * columns.
 * It also creates a summary file, with all read lengths summed up.
 *
 * Created     : 2015-01-13
 * Modified    : 2015-01-15
 * Version     : 0.1
 *
 * Copyright   : 2015 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.1',
        'output_file_format' => '{{GENE}}_{{TRANSCRIPT}}_read_length_distribution.txt',
        'terminal_width' => 120,
    );

echo 'GetReadLengthPerTranscript v.' . $_SETT['version'] . "\n";
//     'PLEASE DO NOT USE THIS SCRIPT ON THE I DRIVE; IT\'S INCREDIBLY SLOW THERE.' . "\n\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 4) {
    die('Usage: ' . $sScriptName . ' FILE_WITH_GENES_TO_ANALYSE GENE_LIST_FILE TRANSCRIPT_POSITION_FILE SAM_FILE1 [SAM_FILE2 [SAM_FILE3 [...]]]' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to read ' . $sFile . '.' . "\n");
    }
}





// First, load file with gene symbols to analyse.
$sGeneFile = array_shift($aFiles);
$aGeneFile = file($sGeneFile, FILE_IGNORE_NEW_LINES);
$aGenesToAnalyze = array(); // '<gene>' => array('<transcript>', '<transcript>', ...);
print('Parsing genes to analyze... ');
foreach ($aGeneFile as $sLine) {
    if (!trim($sLine) || $sLine{0} == '#') {
        // Comment, header.
        continue; // Skip this line, continue to next.
    }
    // We don't check the file much.
    $aGenesToAnalyze[$sLine] = array();
}
unset($aGeneFile);
print('done, loaded ' . count($aGenesToAnalyze) . ' genes in memory.' . "\n");





// Load transcript/gene info file.
$sGeneFile = array_shift($aFiles); // Shifts next argument off the array (the gene info file's name), we don't need it anymore after this.
$aGeneFile = file($sGeneFile, FILE_IGNORE_NEW_LINES);
$aTranscriptsToAnalyze = array(); // transcript => array(chr, strand, transcript_start, transcript_end)
print('Parsing gene/transcript information file... ');
foreach ($aGeneFile as $sLine) {
    if (!trim($sLine) || $sLine{0} == '#') {
        // Comment, header.
        continue; // Skip this line, continue to next.
    }
    list($sTranscript, $sStrand, $sGene) = explode("\t", $sLine);
    if (isset($aGenesToAnalyze[$sGene])) {
        $aGenesToAnalyze[$sGene][] = $sTranscript;
        $aTranscriptsToAnalyze[$sTranscript] = array();
    }
}
unset($aGeneFile);
// If the wrong file has been passed, we have no valid transcripts. Then it will
// make no sense at all to continue.
if (!count($aTranscriptsToAnalyze)) {
    die("\n" .
        'Didn\'t find any valid transcripts. Make sure you passed the correct gene list and gene information file as the first two arguments.' . "\n\n");
}
print('done, loaded ' . count($aTranscriptsToAnalyze) . ' transcripts in memory.' . "\n");





// Prepare transcript locations file, read into memory.
$sTranscriptPositionsFile = array_shift($aFiles);
$aTranscriptPositionsFile = file($sTranscriptPositionsFile, FILE_IGNORE_NEW_LINES);
$nTranscripts = 0;
print('Parsing transcript locations file... ');
foreach ($aTranscriptPositionsFile as $nLine => $sLine) {
    $nLine ++;
    if (!trim($sLine) || $sLine{0} == '#') {
        continue;
    }
    if (preg_match('/^([NX][MR]_\d+\.\d+)\t(\d{1,2}|[XY])\t([+-])\t([\[\]0-9,]+)$/', $sLine, $aRegs)) {
        // Valid transcript position found.
        list(,$sTranscriptWithVersion, $sChr, $sStrand, $sExonPositions) = $aRegs;
        $sTranscript = substr($sTranscriptWithVersion, 0, strpos($sTranscriptWithVersion . '.', '.'));
        if (isset($aTranscriptsToAnalyze[$sTranscript])) {
            if (!($aExonPositions = json_decode($sExonPositions))) {
                die("\n" .
                    'Can\'t parse line ' . $nLine . ' in file ' . $sTranscriptPositionsFile . '.' . "\n\n");
            }
            $nStart = min($aExonPositions[0]);
            $nEnd = max($aExonPositions[count($aExonPositions)-1]);
            $aTranscriptsToAnalyze[$sTranscript] = array('chr' . $sChr, $sStrand, $nStart, $nEnd);
            $nTranscripts ++;
        }
    }
}
unset($aTranscriptPositionsFile);
// If the wrong file has been passed, we have no valid transcripts. Then it will
// make no sense at all to continue.
if (!$nTranscripts) {
    die("\n" .
        'Didn\'t find any valid transcript positions. Make sure you passed the correct transcript location file as the third argument.' . "\n\n");
}
print('done, loaded ' . $nTranscripts . ' transcript positions in memory.' . "\n");
// If this number does not match the number of total transcripts, should we die. We're going to get errors, otherwise...
if (count($aTranscriptsToAnalyze) != $nTranscripts) {
    die('Not all transcripts we should analyze, have annotation. Please update the transcripts annotation file: ' . $sTranscriptPositionsFile . "\n\n");
}
print("\n");





$aData = array('' => array()); // 'NR_000001' => array('<length>' => array('<sample>' => '<coverage>'); (transcript '' is for a big summary)
$aSamples = array(); // To easily see which samples we have.
foreach ($aFiles as $sFile) {
    $nLine = 0;
    // SAM files are BIG. Very, very, very, BIG. Really. HUGE, actually. GBs. We need to read this line by line, to prevent a memory error.
    $fIn = @fopen($sFile, 'r');
    if (!$fIn) {
        die('Unable to open ' . $sFile . '.' . "\n\n");
    }
    $nFileSize = filesize($sFile);
    $nBytesRead = 0;

    // Try and see the sample name, and the type of file. This could be done more generally, but for now we just want to be fast.
    $sBasename = basename($sFile);
    if (preg_match('/^merged_(..).fastq.trunc(.+)_M25.sam$/', $sBasename, $aRegs)) {
        $sSample = $aRegs[1];
        if (strpos($aRegs[2], 'genome') !== false) {
            $sType = 'genome';
        } else {
            $sType = 'transcriptome';
        }
    } else {
        $sSample = $sBasename;
        if (strpos($sBasename, 'genome') !== false) {
            $sType = 'genome';
        } else {
            $sType = 'transcriptome';
        }
    }

    print('Parsing ' . $sSample . ' ' . $sType . ' (' . $sFile . ')' . "\n");
    $sSample = $sSample{0}; // Keep 'A' for 'A1' to 'A3'.
    $aSamples[] = $sSample;
    while ($sLine = fgets($fIn)) {
        $nLine ++;
        $nBytesRead += strlen($sLine);
        $sLine = rtrim($sLine);
        if ($sLine{0} == '@') {
            // The genomic alignment files come with a long list of headers, all starting with a @.
            continue;
        }
        if ($sType == 'transcriptome') {
            list(,, $sReference,,,,,,, $sRead) = explode("\t", $sLine);
            $lRead = strlen($sRead);
            list(,,, $sTranscriptWithVersion) = explode('|', $sReference);
            $sTranscript = substr($sTranscriptWithVersion, 0, strpos($sTranscriptWithVersion . '.', '.'));
            $aTranscriptsMatching = array($sTranscript); // Because we need a loop for the genomic SAM file.
        } else {
            list(, $nBitFlag, $sChr, $nStartPosition,,,,,, $sRead) = explode("\t", $sLine);
            $sStrand = ($nBitFlag & 16? '-' : '+');
            $lRead = strlen($sRead);
            $nEndPosition = $nStartPosition + $lRead;
            // Now check if any of the transcripts we're looking for, overlaps with this location.
            // One genomic location, may match multiple transcripts.
            $aTranscriptsMatching = array();
            foreach ($aTranscriptsToAnalyze as $sTranscript => $aTranscript) {
                list($sChrTranscript, $sStrandTranscript, $nTranscriptStart, $nTranscriptEnd) = $aTranscript;
                if ($sChr == $sChrTranscript && $sStrand == $sStrandTranscript && (($nStartPosition >= $nTranscriptStart && $nStartPosition <= $nTranscriptEnd) || ($nEndPosition >= $nTranscriptStart && $nEndPosition <= $nTranscriptEnd))) {
                    $aTranscriptsMatching[] = $sTranscript;
                }
            }
        }
        if (!$sRead) {
            die('Can\'t parse line ' . $nLine . ' in file ' . $sFile . '.' . "\n\n");
        }
        foreach ($aTranscriptsMatching as $sTranscript) {
            if (!isset($aTranscriptsToAnalyze[$sTranscript])) {
                // We only count transcripts that we are interested in!
                continue;
            }
            if (!isset($aData[$sTranscript])) {
                $aData[$sTranscript] = array();
            }
            if (!isset($aData[$sTranscript][$lRead])) {
                $aData[$sTranscript][$lRead] = array();
            }
            if (!isset($aData[$sTranscript][$lRead][$sSample])) {
                $aData[$sTranscript][$lRead][$sSample] = 0;
            }
            $aData[$sTranscript][$lRead][$sSample] ++;
        }

        if (!($nLine % 50000)) {
            $nPercentageRead = round($nBytesRead/$nFileSize, 2);
            $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
            $lDone = round($nPercentageRead*$nAvailableWidth);
            print(str_repeat(chr(8), $_SETT['terminal_width']) .
                '[' . str_repeat('=', $lDone) . str_repeat(' ', $nAvailableWidth - $lDone) . '] ' . $nLine . ' ' . str_pad(round($nPercentageRead*100), 3, ' ', STR_PAD_LEFT) . '%');
        }
    }
    $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
    print(str_repeat(chr(8), $_SETT['terminal_width']) .
        '[' . str_repeat('=', $nAvailableWidth) . '] ' . $nLine . ' 100%');
    fclose($fIn);

    print("\n" .
        'Done reading ' . $nLine . ' lines.' . "\n");
}



print("\n" .
      'All files done, writing output...' . "\n");
$aSamples = array_unique($aSamples); // Needed when we group replicates.

// Loop through genes, loop through transcripts, write output.
foreach ($aGenesToAnalyze as $sGene => $aTranscripts) {
    print($sGene . '...');
    if (!$aTranscripts) {
        print(' (no transcripts found)' . "\n");
        continue;
    }
    foreach ($aTranscripts as $sTranscript) {
        print(' ' . $sTranscript);
        if (!isset($aData[$sTranscript])) {
            print(' (no reads found)');
        } else {
            ksort($aData[$sTranscript]); // Sort on read length.
            $sFile = str_replace(array('{{GENE}}', '{{TRANSCRIPT}}'), array($sGene, $sTranscript), $_SETT['output_file_format']);
            $fOut = @fopen($sFile, 'w');
            fputs($fOut, 'read_length' . "\t" . implode("\t", $aSamples) . "\r\n");
            foreach ($aData[$sTranscript] as $lRead => $aCoverage) {
                // Prepare total summary, first.
                if (!isset($aData[''][$lRead])) {
                    $aData[''][$lRead] = 0;
                }
                // Then the outfile.
                fputs($fOut, $lRead);
                foreach ($aSamples as $sSample) {
                    // Total summary first.
                    $aData[''][$lRead] += (!isset($aCoverage[$sSample])? 0 : $aCoverage[$sSample]);
                    // Then the output.
                    fputs($fOut, "\t" . (!isset($aCoverage[$sSample])? 0 : $aCoverage[$sSample]));
                }
                fputs($fOut, "\r\n");
            }
            fclose($fOut);
            print(' (total coverage: ' . array_sum(array_map("array_sum", $aData[$sTranscript])) . ')');
        }
    }
    print("\n");
}

// Total summary.
print('Total summary...');
ksort($aData['']); // Sort on read length.
$sFile = str_replace(array('{{GENE}}', '{{TRANSCRIPT}}'), array('ALL', 'ALL'), $_SETT['output_file_format']);
$fOut = @fopen($sFile, 'w');
fputs($fOut, 'read_length' . "\t" . 'coverage' . "\r\n");
foreach ($aData[''] as $lRead => $nCoverage) {
    // Then the outfile.
    fputs($fOut, $lRead . "\t" . $nCoverage . "\r\n");
}
fclose($fOut);
print(' (' . array_sum($aData['']) . ' total coverage)' . "\n");
die('All files done.' . "\n");
?>
