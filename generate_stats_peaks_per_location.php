#!/usr/bin/php
<?php
/*******************************************************************************
 *
 * Generates statistics; the number of peaks per location in a gene (5' UTR,
 * Annotated TIS, Downstream coding, 3' UTR, Multiple. It takes all analysis
 * results files from the find_ORFs.php script, and generates one result file
 * per sample.
 *
 * Created     : 2014-01-08
 * Modified    : 2014-09-23
 * Version     : 0.51
 *
 * Copyright   : 2014 Leiden University Medical Center; http://www.LUMC.nl/
 * Programmer  : Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
 *
 * Changelog   : 0.5     2014-07-07
 *               Now ignoring notices when encountering NM reference sequence
 *               files that do not seem to have a sequence.
 *               0.51    2014-09-23
 *               Renamed "extended_5UTR" category to "unannotated_5UTR".
 *
 *************/

$_SETT =
    array(
        'version' => '0.51',
        'output_suffix' =>
        array(
            'stats' => '.ORF_analysis_results.stats_peaks_per_location.txt',
            'peak_classification' => '.ORF_analysis_results.peaks_classification.txt',
            'peak_classification_5UTR' => '.ORF_analysis_results.peaks_classification_5UTR.txt',
        ),
        'categories' =>
        array(
            '5UTR',
            'annotated_TIS',
            'coding',
            '3UTR',
            'multiple',
        ),
        'NM_cache_dir' => '/home/ifokkema/tmp/ele/new_run/NM_cache/',
        'terminal_width' => 150,
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





echo 'Stats: Peaks Per Location v.' . $_SETT['version'] . "\n";

$aFiles = $_SERVER['argv'];
$sScriptName = array_shift($aFiles);
if (count($aFiles) < 1) { // I guess you could run it with just one file...
    die('Usage: ' . $sScriptName . ' ORF_FILE1 [ORF_FILE2 [ORF_FILE3 [...]]]' . "\n\n");
}

// Check if all files can be read.
$nSamples = count($aFiles);
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

// Go through input files, recognize types. Determine how we should call the new file.
$aSamples = array(); // Will store the info about the samples and which files we have on them.
if ($nSamples == 1) {
    // FIXME: Is this making sense? Shouldn't we try and recognize the file anyways? Like this, we also don't check the file's readability and fill $aSamples.
    $sFileOut = $aFiles[0] . $_SETT['output_suffix']['stats'];
} else {
    foreach ($aFiles as $nFile => $sFile) {
        // Toss possible statistics files or result files out.
        foreach ($_SETT['output_suffix'] as $sSuffix) {
            if (substr($sFile, -(strlen($sSuffix))) == $sSuffix) {
                unset($aFiles[$nFile]);
                continue 2;
            }
        }
        // Also ignore stats.
        if (preg_match('/\.ORF_analysis_results_stats\.txt$/', $sFile)) {
            unset($aFiles[$nFile]);
            continue;
        }
        // Rest is matched.
        if (!preg_match('/^(.+)\.(F|R)(?:\..*)?\.ORF_analysis_results(_after_cutoff)?\.txt$/', $sFile, $aRegs)) {
            //              [1]   [2]                                 [3]
            // Unrecognized file, complain.
            die('Sorry, I do not understand the file name of the file \'' . $sFile . '\', aborting.' . "\n");
        }
        $sSampleID = $aRegs[1]; // Actually, basically the whole prefix until the strand info.
        $sStrand = $aRegs[2];
        $bCutOff = !empty($aRegs[3]);
        if (!isset($aSamples[$sSampleID])) {
            $aSamples[$sSampleID] =
                array(
                    'F' => array(), // File names. We'll end up with two keys here hopefully: false = file with peaks before cutoff, true = file with peaks after cutoff.
                    'R' => array(), // File names. We'll end up with two keys here hopefully: false = file with peaks before cutoff, true = file with peaks after cutoff.
                    'data' => array(
                        false => array_combine($_SETT['categories'], array_fill(0, count($_SETT['categories']), array(0, 0))), // Number of peaks, Total coverage.
                        true => array_combine($_SETT['categories'], array_fill(0, count($_SETT['categories']), array(0, 0))), // Number of peaks, Total coverage.
                    ),
                    'peak_data' => array_combine($_SETT['categories'], array_fill(0, count($_SETT['categories']), array())), // The raw peak data, to be displayed at the bottom of the file.
                    'peak_count' => 0,
                );
        }
        $aSamples[$sSampleID][$sStrand][$bCutOff] = $sFile;
    }
}





// Checking if we are allowed to create the output files.
$aFilesOut = array();
foreach ($aSamples as $sSampleID => $aSample) {
    $aSamples[$sSampleID]['file_out'] = array();
    foreach ($_SETT['output_suffix'] as $sType => $sSuffix) {
        $sFileOut = $sSampleID . $sSuffix;
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

        $aSamples[$sSampleID]['file_out'][$sType] = array('name' => $sFileOut, 'handle' => $fOut);

        // Nicely sort the files, so we always parse them in the same order (before cutoff, after cutoff).
        ksort($aSamples[$sSampleID]['F']);
        ksort($aSamples[$sSampleID]['R']);
    }
}





// Now, loop the ORF analysis files, load them one by one in the memory.
foreach ($aSamples as $sSampleID => $aSample) {
    foreach (array('F', 'R') as $sStrand) {
        $aFiles = $aSample[$sStrand];
        foreach ($aFiles as $bCutOff => $sFile) {
            $aORFFile = file($sFile, FILE_IGNORE_NEW_LINES);
            print('Parsing ' . $sFile . '... ');
            $sGene = '';
            $aTranscripts = array();
            $nPositions = 0;
            foreach ($aORFFile as $sLine) {
                if (!trim($sLine)) {
                    continue;
                }
                // We're looking at the data, or just before.
                // If we don't have a gene yet, look for it.
                if (preg_match('/^(.+)\tPositions found:\t\d+\tPositions analyzed:\t\d+\tT[IS]S found:\t\d+$/', $sLine, $aRegs)) {
                    $sGene = $aRegs[1];
                } elseif ($sGene && preg_match('/^G_Position\tCoverage((?:\tNM_[0-9]+(?:\.[0-9]+)?)+)$/', $sLine, $aRegs)) {
                    $aTranscripts = explode("\t", trim($aRegs[1]));
                } elseif ($sGene && preg_match('/^(chr(?:\d+|[XYM])):(\d+)\t(\d+)((?:\t[0-9*-]+)+)$/', $sLine, $aRegs)) {
                    //                             [1]                [2]    [3]  [4]
                    // We have matched a data line.
                    $nPositions ++;
                    list(,$sChr, $nPosition, $nCoverage) = $aRegs;
                    $aPositions = array();
                    $aPositionsInGene = explode("\t", trim($aRegs[4]));
                    // Loop all positions in genes, and determine category. Store these categories.
                    foreach ($aPositionsInGene as $key => $sPositionInGene) {
                        if ($sPositionInGene == '-') {
                            // There was no mapping on this transcript.
                            continue;
                        }
                        if ($sPositionInGene{0} == '-') {
                            if ($sPositionInGene < -12) {
                                $sCategory = '5UTR';
                            } elseif ($sPositionInGene >= -12 && $sPositionInGene <= -10) {
                                $sCategory = 'annotated_TIS';
                            } else {
                                $sCategory = 'coding';
                            }
                        } elseif ($sPositionInGene{0} == '*') {
                            $sCategory = '3UTR';
                        } else {
                            $sCategory = 'coding';
                        }
                        $aPositions[$key] = $sCategory;
                    }
                    // Now see how many different categories we have for this position.
                    // Apparently, without array_values(), sometimes we have no [0]. Weird.
                    $aPositionsUnique = array_values(array_unique($aPositions));
                     if (count($aPositionsUnique) == 1) {
                        // One transcript, or multiple but at least in the same category of position.
                        $sCategory = $aPositionsUnique[0];
                    } else {
                        // Multiple different positions. If one of them is annotated_TIS, we will assume annotated_TIS.
                        // Otherwise, we don't know what to do, and we call this 'multiple';
                        if (in_array('annotated_TIS', $aPositionsUnique)) {
                            $sCategory = 'annotated_TIS';
                        } else {
                            $sCategory = 'multiple';
                        }
                    }

                    // We have now determined the category. Store, and count.
                    // FIXME: We can also just inject directly into $aSample, so we don't need to reload $aSample later.
                    // Anyways we don't need this information outside of this loop.
                    $aSamples[$sSampleID]['data'][$bCutOff][$sCategory][0] ++;
                    $aSamples[$sSampleID]['data'][$bCutOff][$sCategory][1] += $nCoverage;

                    // v.0.3: Also store the raw peak data, so that we can show it. Only for when BEFORE the cut off.
                    if (!$bCutOff) {
                        // v.0.4: Changed the way the positions are stored; now in a big array.
                        // For 5'UTR peaks, we store *all* positions to be able to report *all* upstream sequences.
                        $aAllPositions = array();
                        if ($sCategory == '5UTR') {
                            foreach ($aPositionsInGene as $key => $sPositionInGene) {
                                if ($sPositionInGene == '-') {
                                    // There was no mapping on this transcript.
                                    continue;
                                }
                                if ($sPositionInGene{0} == '-' && $sPositionInGene < -12 && !isset($aAllPositions[$sPositionInGene])) {
                                    // 5'UTR position, that we haven't seen before (we're looking for unique positions; -15 twice is useless).
                                    $aAllPositions[$sPositionInGene] = array($sPositionInGene + 12, $aTranscripts[$key]);
                                }
                            }
                        }

                        // Pick transcript and report detailed information.
                        $key = (int) array_search($sCategory, $aPositions); // If the search returns false (category = 'multiple'), we'll get the first.
                        $sTranscript = $aTranscripts[$key];
                        $sPositionInGene = $aPositionsInGene[$key];
                        if (substr($sPositionInGene, 0, 1) == '*') {
                            $sPositionInGeneShifted = '*' . (substr($sPositionInGene, 1) + 12);
                        } elseif ($sPositionInGene >= -12 && $sPositionInGene < 0) {
                            $sPositionInGeneShifted = $sPositionInGene + 13; // We're skipping over the -1 -> 1 border here.
                        } else {
                            $sPositionInGeneShifted = $sPositionInGene + 12;
                        }

                        // v.0.4: Changed the way the positions are stored; now in a big array.
                        if (!$aAllPositions) {
                            // No positions stored yet, put in the selected one.
                            $aAllPositions[$sPositionInGene] = array($sPositionInGeneShifted, $sTranscript);
                        }

                        $aSamples[$sSampleID]['peak_data'][$sCategory][] = array($sChr . ':' . $nPosition, $sChr . ':' . ($sStrand == 'F'? $nPosition + 12 : $nPosition - 12), $sGene, $sStrand, $nCoverage, $aAllPositions);
                        $aSamples[$sSampleID]['peak_count'] ++;
                    }
                }
            }
            print('done, loaded ' . $nPositions . ' positions.' . "\n");
        }
    }
    $aSample = $aSamples[$sSampleID]; // Reload, since we're in a foreach and we're working on a copy of the array.





    // Let user know we're working here...
    print('Writing output to ' . $aSample['file_out']['stats']['name'] . '... ');
    fputs($aSample['file_out']['stats']['handle'], '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n" .
        '# NOTE: Read start positions at the end of the coding region, less than 12 bp' . "\n" .
        '# away from the 3\'UTR, are counted as coding while in fact in reality they are' . "\n" .
        '# of course a peak in the 3\'UTR. This can not be detected however, because we' . "\n" .
        '# don\'t know the length of the coding region of the transcripts.' . "\n");

    foreach (array(false, true) as $bCutOff) {
        fputs($aSample['file_out']['stats']['handle'], "\n\n" .
            '# Results for ORF start sites with ' . ($bCutOff? 'no' : 'a') . ' cutoff applied. Using files:' . "\n" .
            '# ' . $aSample['F'][$bCutOff] . "\n" .
            '# ' . $aSample['R'][$bCutOff] . "\n" .
            (!$bCutOff? '' :
                '# ' . $aSample['F'][!$bCutOff] . "\n" .
                '# ' . $aSample['R'][!$bCutOff] . "\n") .
            '# Category' . "\t" . 'Number of TISs found' . "\t" . 'Total coverage' . "\n");

        foreach ($_SETT['categories'] as $sCategory) {
            fputs($aSample['file_out']['stats']['handle'], $sCategory . "\t" . ($aSample['data'][$bCutOff][$sCategory][0] + (!$bCutOff? 0 : $aSample['data'][!$bCutOff][$sCategory][0])) . "\t" . ($aSample['data'][$bCutOff][$sCategory][1] + (!$bCutOff? 0 : $aSample['data'][!$bCutOff][$sCategory][1])) . "\n");
        }
    }
    print('done.' . "\n");



    // v.0.3: Print out the found TISs, sorted on category.
    $aTypes = array_keys($_SETT['output_suffix']);
    unset($aTypes[0]); // Stats removed.
    foreach ($aTypes as $sType) {
        print('Writing output to ' . $aSample['file_out'][$sType]['name'] . '...' . "\n");
        fputs($aSample['file_out'][$sType]['handle'], '# ' . $sScriptName . ' v.' . $_SETT['version'] . "\n" .
            '# NOTE: Read start positions at the end of the coding region, less than 12 bp' . "\n" .
            '# away from the 3\'UTR, are counted as coding while in fact in reality they are' . "\n" .
            '# of course a peak in the 3\'UTR. This can not be detected however, because we' . "\n" .
            '# don\'t know the length of the coding region of the transcripts.' . "\n\n\n" .
            ($sType != 'peak_classification'? '' :
                '# Showing all ORF start sites before the set cutoff (default: 5KB), sorted on category, strand and genomic position.' . "\n") .
            '# NOTE: The PosGenomic+12 field is the genomic position of the TIS, calculated by shifting the position of the start of the read (PeakPosGenomic)' . "\n" .
            '#       by 12 nucleotides downstream in the gene direction. It might be incorrect, since this kind of shifting does not compensate for introns.' . "\n" .
            '# The mentioned coverage is the summed coverage of the replicates.' . "\n" .
            ($sType != 'peak_classification_5UTR'? '' :
                '# The sequence from the found TIS until the annotated TIS is also given, and translated.' . "\n"));



        $aHeaders = array(
            'Category',
            'PeakPosGenomic',
            'PosGenomic+12',
            'Gene',
            'Strand',
            'Coverage',
            'RefSeqID',
            'PeakPosTrans',
            'PosTrans+12',
        );
        if ($sType == 'peak_classification') {
            $aHeaders[] = 'Motif';
        } else {
            $aHeaders[] = 'DNASeqToAUG';
            $aHeaders[] = 'ProtSeqToAUG';
        }
        $nHeaders = count($aHeaders);
        fputs($aSample['file_out'][$sType]['handle'],
            '# ' . implode("\t", $aHeaders) . "\n");

        // Print the peak data.
        $nLine = 0;
        foreach ($aSample['peak_data'] as $sCategory => $aData) {
            foreach ($aData as $aTIS) {
                $nLine ++;
                // For peak_classification_5UTR we only print 5UTR results...
                if ($sType == 'peak_classification_5UTR' && $sCategory != '5UTR') {
                    continue 2; // We will only have more non-5'UTR that follow.
                }

                // v.0.4: Changed the way the positions are stored; now in a big array.
                // For 5'UTR peaks, we report *all* positions and upstream sequences.
                $aAllPositions = array_pop($aTIS);

                // Map data with column names.
                array_unshift($aTIS, $sCategory);
                $aTIS = array_combine($aHeaders, array_pad($aTIS, $nHeaders, ''));

                foreach ($aAllPositions as $sPositionInGene => $aPosition) {
                    list($sPositionInGeneShifted, $sRefSeqID) = $aPosition;

                    if ($sCategory == 'multiple') {
                        // Remove the values we don't show for the Multiple group.
                        $aTIS['RefSeqID'] = $aTIS['PeakPosTrans'] = $aTIS['PosTrans+12'] = '';
                    } else {
                        $aTIS['RefSeqID'] = $sRefSeqID;
                        $aTIS['PeakPosTrans'] = $sPositionInGene;
                        $aTIS['PosTrans+12'] = $sPositionInGeneShifted;
                    }

                    // Now check if we already have the file in the cache; otherwise, download.
                    if ($aTIS['RefSeqID']) {
                        if (!isset($aNMCache[$aTIS['RefSeqID']])) {
                            // File hasn't been parsed yet.
                            $sNMFile = $_SETT['NM_cache_dir'] . $aTIS['RefSeqID'] . '.gb';
                            if (!is_file($sNMFile)) {
                                // In fact, it hasn't been downloaded yet!
                                $fNM = fopen($sNMFile, 'w');
                                $sNM = file_get_contents('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' . $aTIS['RefSeqID'] . '&rettype=gb');
                                if (!$sNM) {
                                    // Failed to download NM.
                                    die("\n" .
                                        'Failed to download NM sequence for ' . $aTIS['RefSeqID'] . "\n");
                                }
                                fputs($fNM, $sNM);
                                fclose($fNM);
                            } else {
                                $sNM = file_get_contents($sNMFile);
                            }

                            // Parse NM, isolate sequence and isolate CDS position.
                            $aNMCache[$aTIS['RefSeqID']] = array();
                            if (!preg_match('/^\s+CDS\s+(\d+)\.\.(\d+)$/m', $sNM, $aRegs)) {
//                                die("\n" .
//                                    'Failed to find CDS for ' . $aTIS['RefSeqID'] . "\n");
                                $nCDSstart = $nCDSend = 0;
                            } else {
                                list(,$nCDSstart, $nCDSend) = $aRegs;
                            }
                            // Get sequence.
                            @list(,$sSequenceRaw) = preg_split('/^ORIGIN\s+$/m', $sNM, 2); // Ignore notices unknown index 1.
                            $sSequence = rtrim(preg_replace('/[^a-z]+/', '', $sSequenceRaw), "\n/");
                            $aNMCache[$aTIS['RefSeqID']] = array($nCDSstart, $nCDSend, $sSequence);
                        }

                        list($nCDSstart, $nCDSend, $sSequence) = $aNMCache[$aTIS['RefSeqID']];
                        // Get Motif or upstream sequence.
                        if ($sType == 'peak_classification' && $nCDSstart) {
                            // Fetch motif.
                            if ($aTIS['PosTrans+12'] < 0) {
                                // Upstream.
                                if (!$nCDSstart) {
                                    // No annotated CDS found (could not parse)
                                    $sMotif = 'could_not_parse_CDS';
                                } elseif ($nCDSstart == 1) {
                                    $sMotif = 'no_5UTR';
                                } elseif (abs($aTIS['PosTrans+12']) > ($nCDSstart-1)) {
                                    // No upstream sequence, or not enough upstream sequence available.
                                    $sMotif = 'unannotated_5UTR';
                                } else {
                                    $sMotif = substr($sSequence, ($nCDSstart-1+$aTIS['PosTrans+12']), 3);
                                }
                            } else {
                                // Compensate 3'UTR reads.
                                if (substr($aTIS['PosTrans+12'], 0, 1) == '*') {
                                    $nPosMotif = substr($aTIS['PosTrans+12'], 1) + $nCDSend;
                                } else {
                                    $nPosMotif = $aTIS['PosTrans+12'];
                                }
                                $sMotif = substr($sSequence, ($nCDSstart+$nPosMotif-2), 3);
                            }
                            $aTIS['Motif'] = $sMotif;
                        } else {
                            // For 5'UTR (all we see here), get the whole upstream sequence.
                            if (!$nCDSstart) {
                                // No annotated CDS found (could not parse)
                                $sUpstreamSequence = 'could_not_parse_CDS';
                            } elseif ($nCDSstart == 1) {
                                $sUpstreamSequence = 'no_5UTR';
                            } elseif (abs($aTIS['PosTrans+12']) > ($nCDSstart-1)) {
                                // No upstream sequence, or not enough upstream sequence available.
                                $sUpstreamSequence = 'unannotated_5UTR';
                            } else {
                                $sUpstreamSequence = substr($sSequence, ($nCDSstart-1+$aTIS['PosTrans+12']), abs($aTIS['PosTrans+12']));
                            }
                            $aTIS['DNASeqToAUG'] = $sUpstreamSequence;

                            // Now, get it translated.
                            $sProteinSequence = RPF_translateDNA($sUpstreamSequence);
                            $aTIS['ProtSeqToAUG'] = $sProteinSequence;
                        }
                    }

                    fputs($aSample['file_out'][$sType]['handle'], implode("\t", $aTIS) . "\n");

                    // Only for 5'UTR classification, we show all. Otherwise, just the first will do.
                    if ($sType != 'peak_classification_5UTR') {
                        break;
                    }
                }

                if (!($nLine % 50)) {
                    $nPercentageRead = round($nLine/$aSample['peak_count'], 2);
                    $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
                    $lDone = round($nPercentageRead*$nAvailableWidth);
                    print(str_repeat(chr(8), $_SETT['terminal_width']) .
                        '[' . str_repeat('=', $lDone) . str_repeat(' ', $nAvailableWidth - $lDone) . '] ' . $nLine . ' ' . str_pad(round($nPercentageRead*100), 3, ' ', STR_PAD_LEFT) . '%');
                }
            }
        }

        $nAvailableWidth = $_SETT['terminal_width'] - 8 - strlen($nLine);
        print(str_repeat(chr(8), $_SETT['terminal_width']) .
            '[' . str_repeat('=', $nAvailableWidth) . '] ' . $nLine . ' 100%' . "\n");
    }
}
?>
