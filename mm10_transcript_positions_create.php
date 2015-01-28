#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * This script generates the mm10_transcript_positions.txt file based on the
 * alignment of the transcriptome FASTA sequences to the GENOME. The results of
 * this alignment is stored in a SAM file, which should be passed to this script
 * as the first argument.
 *
 * Created     : 2013-08-22
 * Modified    : 2013-09-05
 * Version     : 0.1
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
        'version' => '0.1',
        'output' => 'mm10_transcript_positions.txt',
        'unsupported_transcripts_output' => 'transcriptome_alignment_unsupported_transcripts.txt',
        'terminal_width' => 150,
    );

echo 'CreateTranscriptPositions v.' . $_SETT['version'] . "\n" .
     'PLEASE DO NOT USE THIS SCRIPT ON A NETWORK DRIVE; IT CAN BE INCREDIBLY SLOW THERE.' . "\n\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) != 1) {
    die('Usage: ' . $sScriptName . ' SAM_FILE' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// Checking if we are allowed to create the output file.
$aFilesOut = array($_SETT['output'], $_SETT['unsupported_transcripts_output']);
foreach ($aFilesOut as $sFileOut) {
    if (file_exists($sFileOut)) {
        if (!is_writable($sFileOut)) {
            die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
        }
    } elseif (!is_writable(dirname($sFileOut))) {
        die('Can not create ' . $sFileOut . ', aborting.' . "\n");
    }
}
list($sFileOut, $sFileOutUnsupportedTranscripts) = $aFilesOut;





// Open the file, read out line by line.
$sFileIn = $aFiles[0];
$fIn = fopen($sFileIn, 'r');
if (!$fIn) {
    die('Unable to open file for reading: ' . $sFileIn . '.' . "\n\n");
}
$nFileSize = filesize($sFileIn);
$nBytesRead = 0;

$fOut = @fopen($sFileOut, 'w');
if (!$fOut) {
    die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
}
fputs($fOut, '# Generated ' . date('r') . ' by ' . $sScriptName . "\n" .
             '# Transcriptome alignment to the genome taken from ' . $sFileIn . '.' . "\n" .
             '# Transcript' . "\t" . 'Chr' . "\t" . 'Strand' . "\t" . 'Exon_positions' . "\n");
$fOutUnsupportedTranscripts = @fopen($sFileOutUnsupportedTranscripts, 'w');
if (!$fOutUnsupportedTranscripts) {
    die('Unable to open file for writing: ' . $sFileOutUnsupportedTranscripts . '.' . "\n\n");
}
fputs($fOutUnsupportedTranscripts, '# Generated ' . date('r') . ' by ' . $sScriptName . "\n" .
                                   '# Transcriptome alignment to the genome taken from ' . $sFileIn . '.' . "\n" .
                                   '# Transcript' . "\t" . 'Reason_for_failure' . "\n");



$aData = array(); // Will contain transcripts as keys, with an array (chromosome, strand, positions_encoded) as value.
$aUnsupportedTranscripts = array(); // Will contain transcripts as keys, and the reason for rejection as value (no alignment, no real chromosome, twice in file (=chimeric), bad mapping (Del, Ins)).
$nLine = 0;
$sTranscript = '';
while ($sLine = fgets($fIn)) {
    $nLine ++;
    $nBytesRead += strlen($sLine);
    $sLine = rtrim($sLine);
    if (!$sLine || $sLine{0} == '@') {
        continue;
    }
    $sPreviousTranscript = $sTranscript;
    if (substr_count($sLine, "\t") < 5) {
        // This doesn't look like the requested format...
        die('Unable to parse file: ' . $sFileIn . ', line ' . $nLine . ':' . "\n" . $sLine . "\n\n");
    }
    list($sReference, $nBitFlag, $sChromosome, $nPosition, $nQuality, $sCIGAR) = explode("\t", $sLine); // Ignoring all the other cols.
    list(, , , $sTranscript) = explode('|', $sReference);
    $sStrand = ($nBitFlag & 16? '-' : '+');

    // If we see this reference twice, we need to kill both instances.
    if ($sTranscript == $sPreviousTranscript) {
        $aUnsupportedTranscripts[$sTranscript] = 'chimeric';
        unset($aData[$sPreviousTranscript]);
        continue;
    }
    if ($nBitFlag & 4) {
        // Could not be mapped, skip.
        $aUnsupportedTranscripts[$sTranscript] = 'no_alignment';
        continue;
        /*
        flag:
        0 = forward
        1 = template having multiple segments in sequencing
        2 = each segment properly aligned according to the aligner
        4 = segment unmapped
        8 = next segment unmapped
        16 = reverse
         */
    }

    if ($sTranscript && preg_match('/^chr([XYM]|\d{1,2})$/', $sChromosome) && $nPosition && preg_match_all('/^(\d+[MIDNSHP])+$/', $sCIGAR)) {
        // All seem OK. Store basic info first.
        $sChromosome = substr($sChromosome, 3);
        $aData[$sTranscript] = array($sChromosome, $sStrand);

        // Convert the CIGAR string into small bits, for each section.
        preg_match_all('/(\d+)([MIDNSHP])/', $sCIGAR, $aMatches);

        // Now, loop the CIGAR sections to collect the exon's positions.
        $aExonPositions = array();
        $nCurrentPosition = (int) $nPosition;
        $nPrependBases = 0;
        $nAppendBases = 0;
        foreach ($aMatches[0] as $nSection => $sSection) {
            $n = $aMatches[1][$nSection];
            $s = $aMatches[2][$nSection];
            switch ($s) {
                case 'S':
                case 'H':
                    // Soft or Hard clipping; these bases don't align; the first or the last bases of the transcript.
                    // To compensate, we reduce the position with the size of the clipping.
                    // According to the format, S and H can only be used at the start or the end of the CIGAR string.
                    if ($nCurrentPosition == $nPosition) {
                        // At the start of the read.
                        $nPrependBases += $n;
                    } else {
                        $nAppendBases += $n;
                    }
                    break;
                case 'M':
                    // Match; these bases align. This is an exon, and should be counted that way.
                    $nStartPosition = $nCurrentPosition;
                    $nCurrentPosition += $n;
                    $nEndPosition = $nCurrentPosition - 1;
                    $aExonPositions[] = array($nStartPosition, $nEndPosition);
                    break;
                case 'N':
                    // No match; these bases don't align. This is an intron, and should be counted that way.
                    $nCurrentPosition += $n;
                    break;
                case 'D':
                case 'I':
                    // Not allowed; this indicates a change in sequence that we can not handle. This transcript is too different
                    // from the genomic sequence, and can not be used for transcriptome alignment. Report it!
                    if (!isset($aUnsupportedTranscripts[$sTranscript])) {
                        // Not reported before.
                        $aUnsupportedTranscripts[$sTranscript] = 'bad_alignment:' . $s;
                    } else {
                        // Reported before. Just append code.
                        $aUnsupportedTranscripts[$sTranscript] .= $s;
                    }
                    continue 2;
                default:
                    // Unsupported modifier!
                    die('Can\'t parse line ' . $nLine . ', CIGAR string contains unknown modifier ' . $s . '.' . "\n\n");
            }
        }
////////////////////////////////////////////////////////////////////////////////
/*
NR_029642.1        0    chr7     3219189     40    82M
NM_001104543.1    16    chr17    18243570    40    932M8736N124M2671N216M714N792M336N295M19025N156M14S

CIGAR: CIGAR string. The CIGAR operations are given in the following table (set ‘*’ if unavailable):
M alignment match (can be a sequence match or mismatch)
I insertion to the reference
D deletion from the reference
N skipped region from the reference
S soft clipping (clipped sequences present in SEQ)
H hard clipping (clipped sequences NOT present in SEQ)
P padding (silent deletion from padded reference)
= sequence match
X sequence mismatch

• H can only be present as the first and/or last operation.
• S may only have H operations between them and the ends of the CIGAR string.
• For mRNA-to-genome alignment, an N operation represents an intron. For other types of
alignments, the interpretation of N is not defined.
• Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
*/
        // Now handle the clipping.
        if ($nPrependBases) {
            // Clipping at the start. Extend first exon.
            $aExonPositions[0][0] -= $nPrependBases;
        }
        if ($nAppendBases) {
            // Clipping at the end. Extend last exon.
            $aExonPositions[count($aExonPositions)-1][1] += $nAppendBases;
        }
        $aData[$sTranscript][] = json_encode($aExonPositions);
    } elseif (preg_match('/^chr(([XYM]|\d{1,2})_.+_random|Un_.+)$/', $sChromosome)) {
        // Unrecognized ("fake") chromosome...
        $aUnsupportedTranscripts[$sTranscript] = 'weird_alignment:' . $sChromosome;
    } else {
        die("\n" .
            'Can\'t parse line ' . $nLine . ':' . "\n" . $sLine . "\n\n");
    }

    if (!($nLine % 1000)) {
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
    'Done reading ' . $nLine . ' lines, writing output... ');

// First, write unsupported transcripts.
ksort($aUnsupportedTranscripts, SORT_STRING);
foreach ($aUnsupportedTranscripts as $sTranscript => $sReason) {
    fputs($fOutUnsupportedTranscripts, $sTranscript . "\t" . $sReason . "\n");
}
fclose($fOutUnsupportedTranscripts);



// Now, the actual results.
ksort($aData, SORT_STRING);
foreach ($aData as $sTranscript => $aTranscript) {
    fputs($fOut, $sTranscript . "\t" . implode("\t", $aTranscript) . "\n");
}
fclose($fOut);

print('Done.' . "\n" .
      'Identified ' . count($aUnsupportedTranscripts) . ' unsupported transcripts, wrote data for ' . count($aData) . ' transcripts to ' . $sFileOut . '.' . "\n");
?>
