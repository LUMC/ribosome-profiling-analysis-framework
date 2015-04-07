#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * Takes labels and genomic positions from a file manually generated, imports
 * the transcriptomic coordinates from the ORF analysis results file(s), fetches
 * the sequence of the TSS until the end of the transcript, translates the
 * sequence, and reports the protein sequence until the first stop.
 *
 * Created     : 2014-04-25
 * Modified    : 2014-09-23
 * Version     : 0.31
 *
 * Copyright   : 2014-2015 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 * Changelog   : 0.3     2014-05-16
 *               Besides the full DNA sequence, the script now also shows the
 *               DNA sequence up to and including the first stop codon.
 *               0.31    2014-09-23
 *               Renamed "extended_5UTR" category to "unannotated_5UTR".
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
        'version' => '0.31',
        'output_suffix' => '.mRNA_sequence_report.txt',
        'ORF_results_suffix' => '.ORF_analysis_results.peaks_classification.txt',
        'NM_cache_dir' => '/data/NM_cache/',
        'terminal_width' => 100,
    );

function RPF_translateDNA ($sSequence)
{
    // Takes a DNA sequence and returns the protein sequence.
    static $aTranslationTable = array();

    if (!$sSequence || !is_string($sSequence)) {
        return false;
    }

    if (!$aTranslationTable) {
        $aAminoAcids =
            array(
                'A' => array('GCA','GCC','GCG','GCT'),
                'C' => array('TGC','TGT'),
                'D' => array('GAC','GAT'),
                'E' => array('GAA','GAG'),
                'F' => array('TTC','TTT'),
                'G' => array('GGA','GGC','GGG','GGT'),
                'H' => array('CAC','CAT'),
                'I' => array('ATA','ATC','ATT'),
                'K' => array('AAA','AAG'),
                'L' => array('CTA','CTC','CTG','CTT','TTA','TTG'),
                'M' => array('ATG'),
                'N' => array('AAC','AAT'),
                'P' => array('CCA','CCC','CCG','CCT'),
                'Q' => array('CAA','CAG'),
                'R' => array('AGA','AGG','CGA','CGC','CGG','CGT'),
                'S' => array('AGC','AGT','TCA','TCC','TCG','TCT'),
                'T' => array('ACA','ACC','ACG','ACT'),
                'V' => array('GTA','GTC','GTG','GTT'),
                'W' => array('TGG'),
                'Y' => array('TAC','TAT'),
                '*' => array('TAA','TAG','TGA'),
            );

        // Parse it into a easier format for us to use.
        $aTranslationTable = array();
        foreach ($aAminoAcids as $sAminoAcid => $aCodons) {
            foreach ($aCodons as $sCodon) {
                $aTranslationTable[$sCodon] = $sAminoAcid;
            }
        }
    }

    $sTranslatedSequence = '';

    // Loop through sequence in codons.
    $sSequence = strtoupper($sSequence);
    $lSequence = strlen($sSequence);
    for ($i = 0; $i < $lSequence; $i += 3) {
        $sCodon = substr($sSequence, $i, 3);
        if (isset($aTranslationTable[$sCodon])) {
            $sTranslatedSequence .= $aTranslationTable[$sCodon];
        } else {
            $sTranslatedSequence .= '?';
        }
    }
    return $sTranslatedSequence;
}





echo 'mRNA Sequence Report v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 1) { // We need at least one file...
    die('Usage: ' . $sScriptName . ' POSITIONS_FILE [PEAKS_CLASSIFICATION_FILE [PEAKS_CLASSIFICATION_FILE [...]]]' . "\n\n");
}

// I will always also check the current directory for ORF result files...
$h = opendir('.');
if ($h) {
    while (($sFile = readdir($h)) !== false) {
        if (is_file($sFile) && substr($sFile, -strlen($_SETT['ORF_results_suffix'])) == $_SETT['ORF_results_suffix']) {
            $aFiles[] = $sFile;
        }
    }
}
closedir($h);
$aFiles = array_unique($aFiles);

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}

// Check if the NM cache is available...
if (substr($_SETT['NM_cache_dir'], -1) != '/') {
    $_SETT['NM_cache_dir'] .= '/';
}
if (!is_readable($_SETT['NM_cache_dir'])) {
    die('Unable to open the NM cache. Please verify if the path is correct: ' . $_SETT['NM_cache_dir'] . '.' . "\n");
}
if (!is_writable($_SETT['NM_cache_dir'])) {
    die('Unable to write to the NM cache. Please verify if the path is correct: ' . $_SETT['NM_cache_dir'] . '.' . "\n");
}
// Open NM cache.
$aNMCache = array();
$h = opendir($_SETT['NM_cache_dir']);
if (!$h) {
    die('Unexpected error while reading the NM cache (' . $_SETT['NM_cache_dir'] . ').' . "\n");
}
closedir($h);
//while (($sFile = readdir($h)) !== false) {
//    if (is_file($_SETT['NM_cache_dir'] . $sFile) && preg_match('/^(NM_\d+\.\d+)\.gb$/', $sFile, $aRegs)) {
//        $aNMCache[$aRegs[1]] = $_SETT['NM_cache_dir'] . $sFile;
//    }
//}

$sInputFile = array_shift($aFiles);
$sFileOut = $sInputFile . $_SETT['output_suffix'];
if (file_exists($sFileOut)) {
    if (!is_writable($sFileOut)) {
        die('Can not overwrite ' . $sFileOut . ', aborting.' . "\n");
    }
} elseif (!is_writable(dirname($sFileOut))) {
    die('Can not create ' . $sFileOut . ', aborting.' . "\n");
}

$fFileOut = @fopen($sFileOut, 'w');
if (!$fFileOut) {
    die('Unable to open file for writing: ' . $sFileOut . '.' . "\n\n");
}
print('Writing output to ' . $sFileOut . '... ' . "\n");
fputs($fFileOut, '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n");





// Parse the mappings from the ORF result files...
$aMappings = array(); // chr => array(pos => array(nm => pos))
foreach ($aFiles as $sMappingFile) {
    $aMappingFile = file($sMappingFile, FILE_IGNORE_NEW_LINES);
    foreach ($aMappingFile as $sLine) {
        if (!trim($sLine) || $sLine{0} == '#') {
            continue;
        }
        $aLine = explode("\t", $sLine);
        if (count($aLine) == 10) {
            list($sChr, $nPosition) = explode(':', strtolower(trim($aLine[1])), 2); // strtolower() is needed because of typos.
            if (!isset($aMappings[$sChr])) {
                $aMappings[$sChr] = array();
            }
            if (!isset($aMappings[$sChr][$nPosition])) {
                $aMappings[$sChr][$nPosition] = array();
            }
            $aMappings[$sChr][$nPosition][$aLine[6]] = $aLine[7];
        }
    }
    // Mention which files we're using for getting the transcriptomic positions.
    fputs($fFileOut, '# Imported transcriptomic positions from ' . $sMappingFile . "\n");
}
fputs($fFileOut, "\n" .
                 '# PeakPosGenomic' . "\t" . 'RefSeqID' . "\t" . 'PosTrans+12' . "\t" . 'DNASeqToSTOP' . "\t" . 'DNASeq' . "\t" . 'ProtSeqToSTOP' . "\n");





// Read the input file, parse the headers and isolate the positions. Look for the transcriptomic mappings, fetch sequences, translate and report...
$aInputFile = file($sInputFile, FILE_IGNORE_NEW_LINES);
$sLabel = '';
$nLine = 0;
$nLines = count($aInputFile);
foreach ($aInputFile as $sLine) {
    $nLine ++;
    if (!trim($sLine)) {
        continue;
    }
    if ($sLine{0} == '>') {
        if ($sLabel) {
            // Wait... we already had the label... no position!
            die("\n" . '  Error: Label without a position: ' . $sLabel . "\n\n");
        }
        $sLabel = trim($sLine);
        continue;
    }

    list($sChr, $nPosition) = explode(':', strtolower(trim($sLine)), 2); // strtolower() is needed because of typos.
    if ($sChr && $nPosition) {
        // Now look for the mappings.
        if (!isset($aMappings[$sChr][$nPosition])) {
            die("\n" . '  Error: Couldn\'t find transcriptomic position for : ' . $sLabel . ":\n" . $sLine . "\n\n");
        }
    }

    fputs($fFileOut, $sLabel . "\n");
    // Multiple positions are possible, as long as there are different positions.
    foreach ($aMappings[$sChr][$nPosition] as $sTranscript => $nPositionOnTranscript) {
        // Do we already have the transcript parsed?
        if (!isset($aNMCache[$sTranscript])) {
            // File hasn't been parsed yet.
            $sNMFile = $_SETT['NM_cache_dir'] . $sTranscript . '.gb';
            if (!is_file($sNMFile)) {
                // In fact, it hasn't been downloaded yet!
                $fNM = fopen($sNMFile, 'w');
                $sNM = file_get_contents('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' . $sTranscript . '&rettype=gb');
                if (!$sNM) {
                    // Failed to download NM.
                    die("\n" .
                        'Failed to download NM sequence for ' . $sTranscript . "\n");
                }
                fputs($fNM, $sNM);
                fclose($fNM);
            } else {
                $sNM = file_get_contents($sNMFile);
            }

            // Parse NM, isolate sequence and isolate CDS position.
            $aNMCache[$sTranscript] = array();
            if (!preg_match('/^\s+CDS\s+(\d+)\.\.(\d+)$/m', $sNM, $aRegs)) {
//                                die("\n" .
//                                    'Failed to find CDS for ' . $sTranscript . "\n");
                $nCDSstart = $nCDSend = 0;
            } else {
                list(,$nCDSstart, $nCDSend) = $aRegs;
            }
            // Get sequence.
            list(,$sSequenceRaw) = preg_split('/^ORIGIN\s+$/m', $sNM, 2);
            $sSequence = rtrim(preg_replace('/[^a-z]+/', '', $sSequenceRaw), "\n/");
            $aNMCache[$sTranscript] = array($nCDSstart, $nCDSend, $sSequence);
        }



        // Isolate sequence starting from given position (note: still need to compensate + 12).
        list($nCDSstart, $nCDSend, $sSequence) = $aNMCache[$sTranscript];
        $bTranslatable = true;
        $sSequenceProtein = '';

        // Handle 3' UTR positions.
        if ($nPositionOnTranscript{0} == '*') {
            die("Yup, 3' UTR...");
        }

        // Compensate for read length.
        $nPositionOnTranscript += 12;
        if ($nPositionOnTranscript < 12) {
            $nPositionOnTranscript ++;
        }

        // Position was relative to TSS, make it relative to the sequence.
        if (!$nCDSstart) {
            // No annotated CDS found (could not parse)
            $bTranslatable = false;
            $sSequenceProtein = 'could_not_parse_CDS';
        } else {
            $nPositionInString = $nPositionOnTranscript + $nCDSstart - 2; // 0-based position.
            if ($nPositionInString < 0) {
                // We don't have the start of the sequence.
                $bTranslatable = false;
                if ($nCDSstart == 1) {
                    $sSequenceProtein = 'no_5UTR';
                } else {
                    // No upstream sequence, or not enough upstream sequence available.
                    $sSequenceProtein = 'unannotated_5UTR';
                }
            }
        }


        $sSequenceToTranslate = $sSequenceToTranslateToStop = '';
        if ($bTranslatable) {
            $sSequenceToTranslate = substr($sSequence, $nPositionInString);
            $sSequenceProtein = RPF_translateDNA($sSequenceToTranslate);

            // Shorten sequence, only show until the first stop.
            $sSequenceProtein = substr($sSequenceProtein, 0, strpos($sSequenceProtein . '*', '*')+1);

            // Then create a short mRNA sequence, up and until the first stop.
            $sSequenceToTranslateToStop = substr($sSequenceToTranslate, 0, strlen($sSequenceProtein)*3);
        }

        // Output...
        fputs($fFileOut, $sChr . ':' . $nPosition . "\t" . $sTranscript . "\t" . $nPositionOnTranscript . "\t" . $sSequenceToTranslateToStop . "\t" . $sSequenceToTranslate . "\t" . $sSequenceProtein . "\n");
    }

    // Clean up...
    $sLabel = '';

    // Print the progress.
    if (!($nLine % 5)) {
        $nPercentageRead = round($nLine/$nLines, 2);
        $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
        $lDone = round($nPercentageRead*$nAvailableWidth);
        print(str_repeat(chr(8), $_SETT['terminal_width']) .
            '[' . str_repeat('=', $lDone) . str_repeat(' ', $nAvailableWidth - $lDone) . '] ' . $nLine . ' ' . str_pad(round($nPercentageRead*100), 3, ' ', STR_PAD_LEFT) . '%');
    }
}

$nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
print(str_repeat(chr(8), $_SETT['terminal_width']) .
    '[' . str_repeat('=', $nAvailableWidth) . '] ' . $nLine . ' 100%' . "\n\n");
?>
