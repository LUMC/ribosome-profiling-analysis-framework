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
 * Version     : 0.1
 *
 * Copyright   : 2013 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 *************/

$_SETT =
    array(
        'version' => '0.1',
        'output_suffix' => '.coverage_per_gene.txt',
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
unset($aGeneFile[0]); // Remove header (first line).
$aTranscripts = array();
print('Parsing gene file... ');
foreach ($aGeneFile as $sLine) {
    list($sTranscript, $sStrand, $sGene) = explode("\t", $sLine);
    $aTranscripts[$sTranscript] = $sGene;
}
print('done, loaded ' . count($aTranscripts) . ' transcripts in memory.' . "\n\n");



// Now, loop the SAM files, load them one by one, then write the results in a new file.
foreach ($aFiles as $sFile) {
    $aSAMFile = file($sFile, FILE_IGNORE_NEW_LINES);
    $aCoveragePerGene = array();
    $aTranscriptsNotFound = array();
    print('Parsing ' . $sFile . '... ');
    foreach ($aSAMFile as $sLine) {
        list($sTranscript, $nPosition, $nCoverage) = explode("\t", $sLine);
        // Take away version since it doesn't matter and the gene list doesn't have it either.
        $sTranscript = substr($sTranscript, 0, strpos($sTranscript . '.', '.'));
        if (isset($aTranscripts[$sTranscript])) {
            $sGene = $aTranscripts[$sTranscript];
            if (!isset($aCoveragePerGene[$sGene])) {
                $aCoveragePerGene[$sGene] = 0;
            }
            $aCoveragePerGene[$sGene] += $nCoverage;
        } else {
            // We don't know the transcript, it's not in our mm10_gene_list.txt file. Store it, so we can report it and maybe fix it.
            $aTranscriptsNotFound[] = $sTranscript; // No check on uniqueness, this list will contain duplicates.
        }
    }
    $aTranscriptsNotFound = array_unique($aTranscriptsNotFound); // Mention every transcript just once!
    sort($aTranscriptsNotFound); // And sort it on the transcript names.
    print('done, writing output... ');
    arsort($aCoveragePerGene, SORT_NUMERIC); // Associated array Reverse Sort genes on coverage.
    // Now print the results to a new file.
    $f = fopen($sFile . $_SETT['output_suffix'], 'w');
    if (!$f) {
        die('FAILED, can\'t create file.' . "\n");
    }
    // Start by writing down which transcripts were not found.
    fputs($f, '# Transcripts not found (' . count($aTranscriptsNotFound) . '):' . "\t" . implode(',', $aTranscriptsNotFound) . "\n" .
              '# Gene' . "\t"  . 'Total coverage' . "\n");
    foreach ($aCoveragePerGene as $sGene => $nCoverage) {
        fputs($f, $sGene . "\t" . $nCoverage . "\n");
    }
    print('done.' . "\n" .
          '  ' . count($aTranscriptsNotFound) . ' transcripts were not found, details in output file.' . "\n");
}
print("\n");
?>
