#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * This script generates a file with a list of genes with their coverage based
 * on the packed sam file, using the mm10_gene_list.txt file to retrieve the
 * gene symbol for the transcripts.
 *
 * Created     : 2013-08-29
 * Modified    : 2013-08-29
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
        'output_suffix' => '.coverage_per_gene.txt',
        'show_top_unknown_transcripts' => 5,
    );

print('PackedSam2CoveragePerGene v.' . $_SETT['version'] . "\n");

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles); // Shifts first argument off the array (this script's name), we don't need it anymore after this.
if (count($aFiles) < 2) {
    die('Usage: ' . $sScriptName . ' path_to_mm10_gene_list.txt PACKED_SAM_FILE1 [ PACKED_SAM_FILE2 [ ... ]]' . "\n\n");
}

// Check if all files can be read.
foreach ($aFiles as $sFile) {
    if (!is_readable($sFile)) {
        die('Unable to open ' . $sFile . '.' . "\n");
    }
}





// Load transcript/gene info file.
$sGeneFile = array_shift($aFiles); // Shifts next argument off the array (the gene info file's name), we don't need it anymore after this.
$aGeneFile = file($sGeneFile, FILE_IGNORE_NEW_LINES);
$aTranscripts = array();
print('Parsing gene file... ');
foreach ($aGeneFile as $sLine) {
    if (!trim($sLine) || $sLine{0} == '#') {
        // Comment, header.
        continue; // Skip this line, continue to next.
    }
    list($sTranscript, $sStrand, $sGene) = explode("\t", $sLine);
    $aTranscripts[$sTranscript] = $sGene;
}
print('done, loaded ' . count($aTranscripts) . ' transcripts in memory.' . "\n\n");



// Now, loop the SAM files, load them one by one in the memory.
$nSamples = count($aFiles);
$aCoveragePerGene = array();
$aTranscriptsNotFoundPerSample = array();
$aTranscriptsNotFoundWithCoveragePerSample = array();
foreach ($aFiles as $nFile => $sFile) {
    $aSAMFile = file($sFile, FILE_IGNORE_NEW_LINES);
    $aTranscriptsNotFoundWithCoverage = array();
    print('Parsing ' . $sFile . '... ');
    foreach ($aSAMFile as $sLine) {
        list($sTranscript, $nPosition, $nCoverage) = explode("\t", $sLine);
        // Take away version since it doesn't matter and the gene list doesn't have it either.
        $sTranscript = substr($sTranscript, 0, strpos($sTranscript . '.', '.'));
        if (isset($aTranscripts[$sTranscript])) {
            $sGene = $aTranscripts[$sTranscript];
            if (!isset($aCoveragePerGene[$sGene])) {
                $aCoveragePerGene[$sGene] = array_fill(0, $nSamples, 0);
            }
            $aCoveragePerGene[$sGene][$nFile] += $nCoverage;
        } else {
            // We don't know the transcript, it's not in our mm10_gene_list.txt file. Store it, so we can report it.
            if (!isset($aTranscriptsNotFoundWithCoverage[$sTranscript])) {
                $aTranscriptsNotFoundWithCoverage[$sTranscript] = 0;
            }
            $aTranscriptsNotFoundWithCoverage[$sTranscript] += $nCoverage;
        }
    }
    $aTranscriptsNotFoundPerSample[$nFile] = array_keys($aTranscriptsNotFoundWithCoverage); // Copy the list of transcripts not found, without coverage info.
    sort($aTranscriptsNotFoundPerSample[$nFile]); // Then sort this list.
    arsort($aTranscriptsNotFoundWithCoverage, SORT_NUMERIC); // Associated array Reverse Sort transcripts without genes on coverage.
    // Isolate the top # transcripts based on their coverage.
    $aTranscriptsNotFoundWithCoveragePerSample[$nFile] = array();
    for ($i = 0; $i < $_SETT['show_top_unknown_transcripts']; $i++) {
        list($sTranscript, $nCoverage) = each($aTranscriptsNotFoundWithCoverage);
        $aTranscriptsNotFoundWithCoveragePerSample[$nFile][$sTranscript] = $nCoverage;
    }
    print('done.' . "\n" .
        '  ' . count($aTranscriptsNotFoundWithCoverage) . ' transcripts were not found, details will be in the output file.' . "\n");
}



// How should we call the new file?
if ($nSamples == 1) {
    $sFileOut = $aFiles[0] . $_SETT['output_suffix'];
} else {
    // Find prefix and suffix for file(s), to have an output file that matches the name.
    $nPrefixLength = $nSuffixLength = min(array_map('strlen', $aFiles)); // Length of the shortest file name in argument list.
    $sPrefix = substr($aFiles[0], 0, $nPrefixLength); // Limit prefix already to length of shortest file name.
    $sSuffix = substr($aFiles[0], -$nSuffixLength); // Limit suffix already to length of shortest file name.
    foreach ($aFiles as $sFile) {
        for ($i = 0; $i < $nPrefixLength; $i++) {
            if ($sPrefix{$i} != $sFile{$i}) {
                // No match!
                $sPrefix = substr($sPrefix, 0, $i);
                $nPrefixLength = strlen($sPrefix);
                break; // Go to next file.
            }
        }
    }
    foreach ($aFiles as $sFile) {
        for ($i = 1; $i < $nSuffixLength; $i++) {
            if (substr($sSuffix, -$i, 1) != substr($sFile, -$i, 1)) {
                // No match!
                $sSuffix = substr($sSuffix, -($i-1));
                $nSuffixLength = strlen($sSuffix);
                break; // Go to next file.
            }
        }
    }
    $sFileOut = $sPrefix . preg_replace('/^' . preg_quote($sPrefix, '/') . '(.+)' . preg_quote($sSuffix, '/') . '$/', '$1', $aFiles[0]) . '-' . preg_replace('/^' . preg_quote($sPrefix, '/') . '(.+)' . preg_quote($sSuffix, '/') . '$/', '$1', $aFiles[$nSamples-1]) . $sSuffix . $_SETT['output_suffix'];
}

// Open the file and, looping through the samples, print the basic info.
$f = fopen($sFileOut, 'w');
if (!$f) {
    die('FAILED, can\'t create file.' . "\n");
}

fputs($f, '# PackedSam2CoveragePerGene v.' . $_SETT['version'] . "\n" .
          '# NOTE: Transcripts not found happens for transcripts that can not be mapped on the genome. UCSC does not know about these transcripts at all.' . "\n");

foreach ($aFiles as $nFile => $sFile) {
    $nSampleID = $nFile + 1;
    fputs($f, "\n" .
              '# SAMPLE' . str_pad($nSampleID, 2, '0', STR_PAD_LEFT) . "\t" . $sFile . "\n" .
              '# Transcripts not found (' . count($aTranscriptsNotFoundPerSample[$nFile]) . '):' . "\t" . implode(',', $aTranscriptsNotFoundPerSample[$nFile]) . "\n" .
              '# Top unknown transcripts ordered on coverage:' . "\n" .
              '# Transcript' . "\t"  . 'Total coverage' . "\n");

    foreach ($aTranscriptsNotFoundWithCoveragePerSample[$nFile] as $sTranscript => $nCoverage) {
        fputs($f, $sTranscript . "\t" . $nCoverage . "\n");
    }
}





// Let user know we're working here...
print('Writing output to ' . $sFileOut . '... ');

// Now, sort the list of genes properly, so we can loop through it to print the results.
$aGenesInOrder = array_map('array_sum', $aCoveragePerGene);
arsort($aGenesInOrder, SORT_NUMERIC);

// Print the header.
fputs($f, "\n" .
          '# Gene' . "\t" . 'Total coverage');
for ($i = 1; $i <= $nSamples; $i ++) {
    fputs($f, "\t" . 'SAMPLE ' . str_pad($i, 2, '0', STR_PAD_LEFT));
}
fputs($f, "\n");

// Now, loop through these genes, to print the results per sample.
foreach ($aGenesInOrder as $sGene => $nCoverage) {
    fputs($f, $sGene . "\t" . $nCoverage);
    for ($i = 0; $i < $nSamples; $i ++) {
        fputs($f, "\t" . $aCoveragePerGene[$sGene][$i]);
    }
    fputs($f, "\n");
}
print('done.' . "\n\n");
?>
