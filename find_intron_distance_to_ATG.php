#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * FindIntronDistanceToATG analyzes the seq_gene.md file downloaded from:
 * ftp://ftp.ncbi.nih.gov/genomes/M_musculus/mapview/seq_gene.md.gz
 * and generates a list of transcripts that have a very short distance between
 * the ATG and the next intron. All transcripts where the distance is larger
 * than 5 codons (15 bases), are discarded. This value is configurable.
 * Transcripts that have no UTR annotated, and have a safe distance between the
 * ATG and the first 3' intron, are reported elsewhere in a separate error file.
 *
 * Created     : 2013-08-22
 * Modified    : 2013-10-01
 * Version     : 0.2
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
        'version' => '0.2',
        'maximum_distance' => 15, // In bases, the maximum distance to still be reported.
        'output_suffix' => '.distance_intron_to_ATG.txt',
        'output_suffix_error' => '.distance_intron_to_ATG_errors.txt',
        'warnings' => array(
            'first_exon' => 'WARNING: 5\' distance is estimate; transcript may be longer',
            'no_UTR' => 'WARNING: No 5\' UTR annotated',
        ),
    );

echo 'FindIntronDistanceToATG v.' . $_SETT['version'] . "\n\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) != 1) {
    die('Usage: ' . $sScriptName . ' path_to_seq_gene.md' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// Checking if we are allowed to create the output file.
$sFileOut = $aFiles[0] . $_SETT['output_suffix'];
if (file_exists($sFileOut)) {
    if (!is_writable($sFileOut)) {
        die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
    }
} elseif (!is_writable(dirname($sFileOut))) {
    die('Can not create ' . $sFileOut . ', aborting.' . "\n");
}
$sFileOutErr = $aFiles[0] . $_SETT['output_suffix_error'];
if (file_exists($sFileOutErr)) {
    if (!is_writable($sFileOutErr)) {
        die('Can not overwrite ' . $sFileOutErr . ', aborting.' . "\n");
    }
} elseif (!is_writable(dirname($sFileOutErr))) {
    die('Can not create ' . $sFileOutErr . ', aborting.' . "\n");
}





// Start by filtering using grep and cut, data directly into a variable (it's an 9.5MB string).
print('Parsing gene/exon locations file... ');
$aData = explode("\n", `grep -E "UTR|CDS" $aFiles[0] | grep GRCm38 | grep -v '|' | cut -f 3-5,12,14`);
$aTranscripts = array(); // array(transcript => array('UTRs' => array(array(5', 3'), ...), 'CDS' => array(5', 3'), 'distances' => array(5' distance, 3' distance), 'warning' => true|false);
foreach ($aData as $sLine) {
    if (!trim($sLine)) {
        continue;
    }
    list($nPosStart, $nPosEnd, $sStrand, $sType, $sTranscript) = explode("\t", $sLine);
    if (!isset($aTranscripts[$sTranscript])) {
        // Ignore non-coding transcripts, they don't have an ATG.
        if (!preg_match('/^[NX]M_/', $sTranscript)) {
            continue;
        }
        $aTranscripts[$sTranscript] =
            array(
                'UTRs' => array(),
                'CDS' => array(),
                'distances' => array(),
                'warning' => 'first_exon', // Warning is on by default (because the - strand needs that).
            );
    }
    if ($sStrand == '+' && !$aTranscripts[$sTranscript]['distances']) {
        // Plus strand, and no distances calculated yet.
        // 5' UTRs are stored, and after finding the first CDS we're done.
        if ($sType == 'UTR') {
            $aTranscripts[$sTranscript]['UTRs'][] = array($nPosStart, $nPosEnd);
        } elseif ($sType == 'CDS') {
            // If there is no distance between this (first) CDS and the previous UTR, we're in the same exon.
            // The length of the previous UTR is then the upstream distance.
            $nUTRs = count($aTranscripts[$sTranscript]['UTRs']);
            if (!$nUTRs) {
                // WTF! No UTRs before the CDS???
                $aTranscripts[$sTranscript]['distances'][] = 0;
                $aTranscripts[$sTranscript]['warning'] = 'no_UTR';
            } elseif ($aTranscripts[$sTranscript]['UTRs'][$nUTRs-1][1] == ($nPosStart - 1)) {
                // UTR joins directly to this CDS.
                $aTranscripts[$sTranscript]['distances'][] = $aTranscripts[$sTranscript]['UTRs'][$nUTRs-1][1] - $aTranscripts[$sTranscript]['UTRs'][$nUTRs-1][0] + 1;
            } else {
                // UTR was in a different exon; exon starts with ATG.
                $aTranscripts[$sTranscript]['distances'][] = 0;
                $aTranscripts[$sTranscript]['warning'] = false;
            }
            $aTranscripts[$sTranscript]['distances'][] = $nPosEnd - $nPosStart + 1 - 3; // '0' means exon ends directly following the ATG.
            // The number of UTRs define whether we throw a warning about the ATG in the first exon.
            if ($nUTRs > 1) {
                $aTranscripts[$sTranscript]['warning'] = false;
            }
        }

    } elseif ($sStrand == '-') {
        // Minus strand.
        // We will ignore the first UTRs; every CDS gets stored, overwriting the previous one,
        // when first encountering an UTR after the CDS we can calculate the distances.
        if ($sType == 'UTR') {
            // Only if we already have seen a CDS will we continue.
            if (!$aTranscripts[$sTranscript]['CDS']) {
                continue;
            }
            // UTR following an CDS. If this is the first UTR, its length will
            // define the 5' distance. The length of the last CDS defines the 3'
            // distance. If more UTRs follow, the "first exon" warning is
            // dropped, unless this UTR does not connect perfectly to the CDS.
            if (!$aTranscripts[$sTranscript]['distances']) {
                // First UTR after the CDS.
                if ($aTranscripts[$sTranscript]['CDS'][1] == ($nPosStart - 1)) {
                    // CDS joins directly to this UTR.
                    $aTranscripts[$sTranscript]['distances'][] = $nPosEnd - $nPosStart + 1; // '0' means exon ends directly following the ATG.
                } else {
                    // CDS was in a different exon; exon starts with ATG.
                    $aTranscripts[$sTranscript]['distances'][] = 0;
                    $aTranscripts[$sTranscript]['warning'] = false;
                }
                $aTranscripts[$sTranscript]['distances'][1] = $aTranscripts[$sTranscript]['CDS'][1] - $aTranscripts[$sTranscript]['CDS'][0] + 1 - 3;
            } else {
                // We've seen the CDS, we've seen an UTR before. This is a
                // non-coding exon, so we can remove the warning.
                $aTranscripts[$sTranscript]['warning'] = false;
            }
        } elseif ($sType == 'CDS') {
            // We store CDS, overwriting the previous one if present.
            $aTranscripts[$sTranscript]['CDS'] = array($nPosStart, $nPosEnd);
        }
    }
}
print('Done, ' . count($aData) . ' lines parsed.' . "\n");
unset($aData);



// Now, we need to filter the data. Instead of removing transcripts from the big
// array based on their distance, I will just make a new list, easy to sort.
// Then I can remove $aTranscripts as a whole and then print the results.
$aData = array(); // array(transcript => array(5', 3', warning, transcript), ...);
foreach ($aTranscripts as $sTranscript => $aTranscript) {
    $aDistance = $aTranscript['distances'];
    if (!count($aDistance)) {
        // Some genes do not have UTRs annotated. But for the genes in reverse I
        // rely on the presence of the UTRs to be able to count the 5' distance.
        // So I need to repeat the calculation...

        // Only if we already have seen a CDS will we continue.
        if (!$aTranscript['CDS']) {
            // No CDSs either... :S
            die('Impossible transcript: ' . $sTranscript . "\n");
        }
        // No distance between start of transcript and CDS.
        $aDistance[] = 0;
        $aDistance[] = $aTranscript['CDS'][1] - $aTranscript['CDS'][0] + 1 - 3;
        $aTranscript['warning'] = 'no_UTR';
    }

    if (min($aDistance) <= $_SETT['maximum_distance']) {
        $aData[$sTranscript] = array_merge($aDistance, array($aTranscript['warning']));

        // 2013-09-25; OK, NOW I'm too lazy. I need to change the sorting algorithm
        // AGAIN because sometimes it's returning 0, which makes the order of the
        // elements undefined when the distances are exactly equal. However, that
        // makes comparing the files completely unreliable, if run again on a
        // different source file. So I always have to define the order. I need the
        // transcript ID for it, but I can't use it because uasort() does not
        // receive the key. So I just include the transcript ID as a value, which
        // means I'm using more memory than necessary, but whatever...
        $aData[$sTranscript][] = $sTranscript;
    }
}
unset($aTranscripts);



// Sorting is not too easy, and I need a custom function for it.
print('Sorting... ');
uasort($aData, function ($a, $b) {
    // Properly sort the results...

    // Determine if a AND b sort on 5' or 3'. if a on 5' and b on 3', sort using that.
    // otherwise, if sorting on the same thing, check the values.
    if ($a[0] <= $a[1] && $b[0] > $b[1]) {
        // 5' versus 3' preference.
        return -1;
    } elseif ($a[0] > $a[1] && $b[0] <= $b[1]) {
        // 3' versus 5' preference.
        return 1;
    } elseif ($a[0] <= $a[1]) {
        // 5' is shorter in A than 3', or equal.
        return ($a[0] > $b[0]? -1 : ($a[0] < $b[0]? 1 : ($a[1] > $b[1]? -1 : ($a[1] < $b[1]? 1 : strcmp($a[3], $b[3])))));

    } else {
        // 3' shorter.
        return ($a[1] < $b[1]? -1 : ($a[1] > $b[1]? 1 : ($a[0] < $b[0]? -1 : ($a[0] > $b[0]? 1 : strcmp($a[3], $b[3])))));
    }
});
print('Done.' . "\n");



print('Writing output... ');
$fOut = @fopen($sFileOut, 'w');
if (!$fOut) {
    die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
}
$fOutErr = @fopen($sFileOutErr, 'w');
if (!$fOutErr) {
    die('Unable to open file for writing: ' . $sFileOutErr . '.' . "\n\n");
}
$i = 0;
fputs($fOut, '#transcript' . "\t" . '5\' distance' . "\t" . '3\' distance' . "\n");
fputs($fOutErr, '#transcript' . "\t" . '5\' distance' . "\t" . '3\' distance' . "\n");
foreach ($aData as $sTranscript => $aTranscript) {
    // Unset unnecessary transcript ID.
    unset($aTranscript[3]);

    // Remove 'first_exon' warning when the 5' distance is big enough, anyways.
    if ($aTranscript[2] == 'first_exon' && $aTranscript[0] > $_SETT['maximum_distance']) {
        $aTranscript[2] = false;
    }
    // Transcripts without UTR and plenty of space 3' are moved to a different file, because we really can't say if they match our criteria or not.
    if ($aTranscript[2] == 'no_UTR' && $aTranscript[1] > $_SETT['maximum_distance']) {
        fputs($fOutErr, $sTranscript . "\t" . implode("\t", $aTranscript) . "\n");
        continue;
    }

    if (isset($_SETT['warnings'][$aTranscript[2]])) {
        // Replace code by human readable text.
        $aTranscript[2] = $_SETT['warnings'][$aTranscript[2]];
    }

    fputs($fOut, $sTranscript . "\t" . implode("\t", $aTranscript) . "\n");
    $i ++;
}
print('Done, ' . $i . ' transcripts found within threshold of ' . $_SETT['maximum_distance'] . '.' . "\n");
?>
