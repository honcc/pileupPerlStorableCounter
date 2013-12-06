#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to read the perl pileup storables generated from BAMToReadEndPerlStorable and count the reads on the features on GFF. It is multi-threaded for faster counting. It’ll output html as the results.
#
#	Input
#		--midPtPileupStorableIndexPath=		file path[compulsory]; path of the index file of the pileupStorable in midPt mode of to be counted;
#		--fullPileupStorableIndexPath=		file path[compulsory]; path of the index file of the pileupStorable in full mode of to be counted;
#		--gffPath=							file path[compulsory]; path of the reference GFF for gene annotation;
#		--fastaPath=						file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--maxThread=						positive integer [4]; max number of threads to be used;
#		--margin5End=						positive or negative interger [0]; the number of bases to extend (if positive) or trim (if negative) at the 5'end of the features; e.g. -30 meaning the 1st 30nt bases of the ranges (e.g. CDSRng) will be ignored; +30 will extend the 5'end by 30nt upstream;
#		--margin3End=						positive or negative interger [0]; the number of bases to extend (if positive) or trim (if negative) at the 5'end of the features; e.g. -30 meaning the last 30nt bases of the ranges (e.g. CDSRng) will be ignored; +30 will extend the 3'end by 30nt dnstream;
#		--minCov=							positive non-zero integer [1]; minimum coverage of a position defined as “covered”;
#		--outDir=							directory path ['./BAMToReadEndPerlStorable/']; output directory;
#
#	Usage
#		/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableCounter/v0.1/pileupPerlStorableCounter_v0.1.pl --midPtPileupStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/N26_footprint/BAMToReadEndPerlStorable/countMode.5.offset.25.baseComp.no/cntgCovPls/index.hsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927_TriTrypDB-4.2.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta
#
#	v0.1
#		[Wed 21 Aug 2013 20:43:41 CEST] debut;
#
#	v0.2
#		[Sat 24 Aug 2013 16:37:53 CEST] added --margin5End= --margin3End= option
#	
#	v0.3
#		[29/09/2013 17:48] added pileupStorableIndexPath changed to midPtPileupStorableIndexPath and fullPileupStorableIndexPath added, will be used to count coverage;
#
#	v0.4
#		[07/11/2013 11:13] will print an overall gene range count file with sense and antisense count with independent gene ID, for the purpose of inputing into DGE softwares;
#
#	v0.5
#		[12/11/2013 12:55] added metagene coverage plot and GC bias plot;
#		[14/11/2013 17:31] added minCov option;
#		[18/11/2013 13:37] added count contig based info subroutine, countContigReadAndGC
#		[19/11/2013 16:13] added HTML output
#		
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-11-19 17:06]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableCounter/v0.5/pileupPerlStorableCounter_v0.5.pl --midPtPileupStorableIndexPath=/Volumes/C_Analysis/NGS/results/E014_heatShock_polyA_basic_min35nt/E014_heatShock_HM1_T2H_polyA_1/BAMToReadEndPerlStorable/countMode.midPt.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --fullPileupStorableIndexPath=/Volumes/C_Analysis/NGS/results/E014_heatShock_polyA_basic_min35nt/E014_heatShock_HM1_T2H_polyA_1/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --maxThread=6 --margin5End=0 --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableCounter/v0.5/E014_heatShock_HM1_T2H_polyA_1/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableCounter/v0.5/pileupPerlStorableCounter_v0.5.pl
#	--midPtPileupStorableIndexPath=
#	--fullPileupStorableIndexPath=
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--maxThread=6
#	--margin5End=0
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableCounter/v0.5/E014_heatShock_HM1_T2H_polyA_1/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $global_scriptDir = dirname(rel2abs($0));
open DEBUGLOG, ">", "$global_scriptDir/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1607, readParameters|1939
#	secondaryDependOnSub: currentTime|928
#
#<section ID="startingTasks" num="0">
########################################################################## 
&printCMDLogOrFinishMessage("CMDLog");#->1607
my ($midPtPileupStorableIndexPath, $fullPileupStorableIndexPath, $gffPath, $fastaPath, $maxThread, $margin5End, $margin3End, $minCov, $outDir) = &readParameters();#->1939

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $metaGeneMinLen = 1000;
my $metaGeneMaxLen = 9999;
my $GCWinSize = 100;
my $metaGeneMinCov = 5;
my $paramTag = "M5End.$margin5End.M3End.$margin3End.minCov.$minCov.GCWinSize.$GCWinSize";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag/";
my $resultGenStatDir = "$resultDir/generalStats/"; push @mkDirAry, $resultGenStatDir;
my $resultContigCountDir = "$resultDir/contigCount/"; push @mkDirAry, $resultContigCountDir;
my $resultGeneCountDir = "$resultDir/geneCount/"; push @mkDirAry, $resultGeneCountDir;
my $resultGCPctDir = "$resultDir/GCPct/"; push @mkDirAry, $resultGCPctDir;
my $resultLenPctDir = "$resultDir/metaGenePlot/"; push @mkDirAry, $resultLenPctDir;
my $resultStorableDir = "$resultDir/storable/"; push @mkDirAry, $resultStorableDir;
my $resultHTMLDir = "$resultDir/html/"; push @mkDirAry, $resultHTMLDir;
my $ggplotDirHsh_ref = {};
my @ggplotFileTypeAry = qw /dat pdf R log/;
foreach my $fileType (@ggplotFileTypeAry) {$ggplotDirHsh_ref->{$fileType} = "$resultDir/ggplot/$fileType"; push @mkDirAry, $ggplotDirHsh_ref->{$fileType};}
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: checkGeneInfo|368, generateGeneByCntgHsh|946, getIndivCntgCovPlsPath|1300, readGFF_oneRNAPerGene|1786, readMultiFasta|1885, zipUnzipCntgCovInPlsPathHsh|2003
#	secondaryDependOnSub: currentTime|928, reportStatus|1982
#
#<section ID="processInputData" num="4">
########################################################################## 

my ($fastaHsh_ref) = &readMultiFasta($fastaPath);#->1885

my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->1786
&checkGeneInfo($geneInfoHsh_ref);#->368
my ($geneByCntgHsh_ref) = &generateGeneByCntgHsh($geneInfoHsh_ref);#->946

my ($midPtPileupStorablePathHsh_ref) = &getIndivCntgCovPlsPath($midPtPileupStorableIndexPath);#->1300
my ($fullPileupStorablePathHsh_ref) = &getIndivCntgCovPlsPath($fullPileupStorableIndexPath);#->1300
&zipUnzipCntgCovInPlsPathHsh('unzip', $midPtPileupStorablePathHsh_ref);#->2003
&zipUnzipCntgCovInPlsPathHsh('unzip', $fullPileupStorablePathHsh_ref);#->2003

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_countContigBasedInfo
#	primaryDependOnSub: countContigReadAndGC|396
#	secondaryDependOnSub: generateThreadHshWithRandomCntg|1273, reportStatus|1982
#
#<section ID="countContigBasedInfo" num="5">
my ($cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref) = &countContigReadAndGC($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $maxThread, $GCWinSize, $minCov, $fastaHsh_ref);#->396
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_countThePileupStorable
#	primaryDependOnSub: calculateNormalizedCount|290, countCtgryOnCntg|517, countGeneOnCntg|668, generateMetaGeneCoverageAndGCBias|1124
#	secondaryDependOnSub: countIndividualGeneRng|829, generateThreadHshWithRandomCntg|1273, getIndivCntgCovPlsPath|1300, reportStatus|1982
#
#<section ID="countThePileupStorable" num="6">
my ($geneIndivDataHsh_ref, $ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryGeneRngSumHsh_ref, $geneIndivCovIdxPathHsh_ref) = &countGeneOnCntg($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneInfoHsh_ref, $geneByCntgHsh_ref, $maxThread, $margin5End, $margin3End, $resultStorableDir, $minCov);#->668

my ($metaGeneLenPctCovSumHsh_ref, $GCPctCovHsh_ref, $metaGeneNumCountHsh_ref) = &generateMetaGeneCoverageAndGCBias($geneIndivCovIdxPathHsh_ref, $metaGeneMinLen, $metaGeneMaxLen, $metaGeneMinCov, $geneInfoHsh_ref, $GCWinSize, $fastaHsh_ref, $maxThread);#->1124

my ($ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $nonStructuralCount) = &countCtgryOnCntg($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $margin5End, $margin3End, $minCov);#->517
my ($mRNASCountUpperQuantile, $countCovDataHsh_ref, $geneRngLenHsh_ref) = &calculateNormalizedCount($geneIndivDataHsh_ref, $geneInfoHsh_ref, $nonStructuralCount, $resultStorableDir, $margin5End, $margin3End);#->290
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_metagene and GC plot
#	primaryDependOnSub: outputContigReadAndGC|1369, outputGCPctData|1455, outputMetaGeneCovData|1534
#	secondaryDependOnSub: ggplotXYLinesMultipleSamples|1334, reportStatus|1982
#
#<section ID="metagene and GC plot" num="7">
my $allResultFilePathHsh_ref = {};
&outputGCPctData($GCPctCovHsh_ref, $resultGCPctDir, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);#->1455
&outputMetaGeneCovData($metaGeneLenPctCovSumHsh_ref, $resultLenPctDir, $metaGeneMinLen, $metaGeneMaxLen, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);#->1534
&outputContigReadAndGC($resultContigCountDir, $cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);#->1369
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_printStatisticsLog
#	primaryDependOnSub: printFeatureCount|1640, printReadStatistics|1715
#	secondaryDependOnSub: reportStatus|1982
#
#<section ID="printStatisticsLog" num="8">
&printReadStatistics($ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $ctgryGeneRngSumHsh_ref, $resultGenStatDir, $minCov, $allResultFilePathHsh_ref);#->1715
&printFeatureCount($resultGeneCountDir, $geneInfoHsh_ref, $ctgryDirtnReadCountHsh_ref, $countCovDataHsh_ref, $geneRngLenHsh_ref, $allResultFilePathHsh_ref);#->1640
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_outputHTML
#	primaryDependOnSub: generateMasterHTML|967
#	secondaryDependOnSub: >none
#
#<section ID="outputHTML" num="9">
&generateMasterHTML($allResultFilePathHsh_ref, $resultHTMLDir, $resultDir);#->967
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 10_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1607
#	secondaryDependOnSub: currentTime|928
#
#<section ID="finishingTasks" num="10">
&printCMDLogOrFinishMessage("finishMessage");#->1607
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	HTML [n=1]:
#		generateMasterHTML
#
#	coverage [n=1]:
#		generateMetaGeneCoverageAndGCBias
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=7]:
#		checkGeneInfo, currentTime, printCMDLogOrFinishMessage
#		readGFF_oneRNAPerGene, readMultiFasta, readParameters
#		reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, generateGeneByCntgHsh, readGFF_oneRNAPerGene
#
#	ggplot [n=1]:
#		ggplotXYLinesMultipleSamples
#
#	multithread [n=1]:
#		generateThreadHshWithRandomCntg
#
#	plotInR [n=1]:
#		ggplotXYLinesMultipleSamples
#
#	reporting [n=1]:
#		currentTime
#
#	specific [n=4]:
#		countContigReadAndGC, outputContigReadAndGC, outputGCPctData
#		outputMetaGeneCovData
#
#	storable [n=2]:
#		getIndivCntgCovPlsPath, zipUnzipCntgCovInPlsPathHsh
#
#	unassigned [n=6]:
#		calculateNormalizedCount, countCtgryOnCntg, countGeneOnCntg
#		countIndividualGeneRng, printFeatureCount, printReadStatistics
#
#====================================================================================================================================================#

sub calculateNormalizedCount {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 6_countThePileupStorable|181
#	secondaryAppearInSection: >none
#	input: $geneIndivDataHsh_ref, $geneInfoHsh_ref, $margin3End, $margin5End, $nonStructuralCount, $resultStorableDir
#	output: $countCovDataHsh_ref, $geneRngLenHsh_ref, $mRNASCountUpperQuantile
#	toCall: my ($mRNASCountUpperQuantile, $countCovDataHsh_ref, $geneRngLenHsh_ref) = &calculateNormalizedCount($geneIndivDataHsh_ref, $geneInfoHsh_ref, $nonStructuralCount, $resultStorableDir, $margin5End, $margin3End);
#	calledInLine: 191
#....................................................................................................................................................#
	my ($geneIndivDataHsh_ref, $geneInfoHsh_ref, $nonStructuralCount, $resultStorableDir, $margin5End, $margin3End) = @_;
	
	&reportStatus("Normalizing read counts", 20, "\n");#->1982

	my $countCovDataHsh_ref = {};
	my $geneRngLenHsh_ref = {};
	
	#---get upperQuantile of mRNA sense count
	my @allmRNASenseCount = ();
	foreach my $geneID (keys %{$geneIndivDataHsh_ref->{'readCount'}{'geneRng'}}) {
		next if $geneInfoHsh_ref->{$geneID}{'ctgry'} !~ m/mRNA/;
		push @allmRNASenseCount, $geneIndivDataHsh_ref->{'readCount'}{'geneRng'}{$geneID}{'s'};
	}
	
	my $statObj = Statistics::Descriptive::Full->new();
	$statObj->add_data(@allmRNASenseCount);
	my $mRNASCountUpperQuantile = $statObj->percentile(75);
	
	foreach my $rngType (keys %{$geneIndivDataHsh_ref->{'readCount'}}) {
		foreach my $geneID (keys %{$geneIndivDataHsh_ref->{'readCount'}{$rngType}}) {
			
			#---set everything as -1
			$geneRngLenHsh_ref->{$rngType}{$geneID} = -1;
			foreach my $dirtn (qw/a s/) {
				$countCovDataHsh_ref->{$rngType}{$geneID}{'readCount'}{$dirtn} = $geneIndivDataHsh_ref->{'readCount'}{$rngType}{$geneID}{$dirtn};
				$countCovDataHsh_ref->{$rngType}{$geneID}{'covPos'}{$dirtn} = $geneIndivDataHsh_ref->{'covPos'}{$rngType}{$geneID}{$dirtn};
				$countCovDataHsh_ref->{$rngType}{$geneID}{'covSum'}{$dirtn} = $geneIndivDataHsh_ref->{'covSum'}{$rngType}{$geneID}{$dirtn};

				foreach my $dataType (qw/RPKM UQNor RPM covPct covPerNt/) {
					$countCovDataHsh_ref->{$rngType}{$geneID}{$dataType}{$dirtn} = -1;
				}
			}
			
			next if not $geneInfoHsh_ref->{$geneID}{$rngType};

			#---get the rng length
			my $rngLength = 0;
			for (my $i=0; $i < $#{$geneInfoHsh_ref->{$geneID}{$rngType}}; $i += 2) {
				$rngLength += ${$geneInfoHsh_ref->{$geneID}{$rngType}}[$i+1] - ${$geneInfoHsh_ref->{$geneID}{$rngType}}[$i];
			}

			$rngLength += $margin5End+$margin3End;
			
			$geneRngLenHsh_ref->{$rngType}{$geneID} = $rngLength;

			#---get the normalized count
			foreach my $dirtn (qw/a s/) {
				#print $rngType."\t".$geneID."\t".$rngLength."\t".$nonStructuralCount."\n";
				eval {$countCovDataHsh_ref->{$rngType}{$geneID}{'RPM'}{$dirtn} = sprintf "%.5f", ($geneIndivDataHsh_ref->{'readCount'}{$rngType}{$geneID}{$dirtn}/($nonStructuralCount/1000000));};
				eval {$countCovDataHsh_ref->{$rngType}{$geneID}{'UQNor'}{$dirtn} = sprintf "%.5f", $geneIndivDataHsh_ref->{'readCount'}{$rngType}{$geneID}{$dirtn}/$mRNASCountUpperQuantile;};
				eval {$countCovDataHsh_ref->{$rngType}{$geneID}{'RPKM'}{$dirtn} = sprintf "%.5f", ($geneIndivDataHsh_ref->{'readCount'}{$rngType}{$geneID}{$dirtn}/($nonStructuralCount/1000000))/($rngLength/1000);}; #---in case if $rngLength = 0
				eval {$countCovDataHsh_ref->{$rngType}{$geneID}{'covPct'}{$dirtn} = (100*$geneIndivDataHsh_ref->{'covPos'}{$rngType}{$geneID}{$dirtn}/$rngLength);}; #---in case if $rngLength = 0
				
				#---[29/09/2013 21:37] ah hoc DEBUG as sometime $rngLength < $geneIndivDataHsh_ref->{'covPos'}{$rngType}{$geneID}{$dirtn}
				$countCovDataHsh_ref->{$rngType}{$geneID}{'covPct'}{$dirtn} = 100 if $countCovDataHsh_ref->{$rngType}{$geneID}{'covPct'}{$dirtn}>100;
				$countCovDataHsh_ref->{$rngType}{$geneID}{'covPct'}{$dirtn} = sprintf "%.5f", $countCovDataHsh_ref->{$rngType}{$geneID}{'covPct'}{$dirtn};

				eval {$countCovDataHsh_ref->{$rngType}{$geneID}{'covPerNt'}{$dirtn} = sprintf "%.5f", ($geneIndivDataHsh_ref->{'covSum'}{$rngType}{$geneID}{$dirtn}/$rngLength);}; #---in case if $rngLength = 0
			}
		}
	}
	
	store($countCovDataHsh_ref, "$resultStorableDir/countCovDataHsh.pls");

	return ($mRNASCountUpperQuantile, $countCovDataHsh_ref, $geneRngLenHsh_ref);
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: 
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 159
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 0, "\n");#->1982
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 0, "\n");#->1982
	}
	
	return ();
}
sub countContigReadAndGC {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: generateThreadHshWithRandomCntg|1273, reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 5_countContigBasedInfo|171
#	secondaryAppearInSection: >none
#	input: $GCWinSize, $fastaHsh_ref, $fullPileupStorablePathHsh_ref, $maxThread, $midPtPileupStorablePathHsh_ref, $minCov
#	output: $cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref
#	toCall: my ($cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref) = &countContigReadAndGC($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $maxThread, $GCWinSize, $minCov, $fastaHsh_ref);
#	calledInLine: 176
#....................................................................................................................................................#
	my ($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $maxThread, $GCWinSize, $minCov, $fastaHsh_ref) = @_;
	
	&reportStatus("Counting contig based info", 10, "\n");#->1982
	my @cntgAry = (keys %{$midPtPileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->1273
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1982

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;

				my $cntgCountDataHsh_inThr_ref = {};
				my $wholeGenomeGCDataHsh_inThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $midPtCntgCovAry_ref = retrieve($midPtPileupStorablePathHsh_ref->{$cntg});
					my $fullCntgCovAry_ref = retrieve($fullPileupStorablePathHsh_ref->{$cntg});
					
					#---[18/11/2013 13:17]  generate GC bias data
					my %tmpBaseCountHsh = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
					my $covSum = 0;
					my %covPosHsh = ('+'=>0, '-'=>0);
					my %rdCountHsh = ('+'=>0, '-'=>0);
					
					foreach my $i (0..$#{$fullCntgCovAry_ref}) {

						my $base = substr $fastaHsh_ref->{$cntg}, $i, 1;
						$tmpBaseCountHsh{$base}++;

						my %tmpCovHsh = (
							'midPt'=>{'+'=>0, '-'=>0},
							'full'=>{'+'=>0, '-'=>0},
						);
						
						($tmpCovHsh{'midPt'}{'+'}, $tmpCovHsh{'midPt'}{'-'}) = split /,/, $midPtCntgCovAry_ref->[$i] if $midPtCntgCovAry_ref->[$i];
						($tmpCovHsh{'full'}{'+'}, $tmpCovHsh{'full'}{'-'}) = split /,/, $fullCntgCovAry_ref->[$i] if $fullCntgCovAry_ref->[$i];
						
						foreach my $strnd (qw/+ -/) {
							$covSum += $tmpCovHsh{'full'}{$strnd};
							$covPosHsh{$strnd}++ if $tmpCovHsh{'full'}{$strnd} >= $minCov;
							$rdCountHsh{$strnd} += $tmpCovHsh{'midPt'}{$strnd};
						}

						if ($i % $GCWinSize == 0 and $i > 0) {
							my $GCPct = sprintf "%.0f", 100*($tmpBaseCountHsh{'C'} + $tmpBaseCountHsh{'G'})/$GCWinSize;
							my $covPerNt = sprintf "%.2f", $covSum/$GCWinSize;
							push @{$wholeGenomeGCDataHsh_inThr_ref->{$GCPct}}, $covPerNt;
							%tmpBaseCountHsh = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
							$covSum = 0;
						}
					}
					
					foreach my $strnd (qw/+ -/) {
						my $covPct = sprintf "%.2f", 100*$covPosHsh{$strnd}/@{$fullCntgCovAry_ref};
						$cntgCountDataHsh_inThr_ref->{$cntg}{'cntgLength'} = @{$fullCntgCovAry_ref};
						$cntgCountDataHsh_inThr_ref->{$cntg}{'covPct'}{$strnd} = $covPct;
						$cntgCountDataHsh_inThr_ref->{$cntg}{'covPos'}{$strnd} = $covPosHsh{$strnd};
						$cntgCountDataHsh_inThr_ref->{$cntg}{'rdCount'}{$strnd} = $rdCountHsh{$strnd};
					}
					
					&reportStatus("Finished counting features on $cntgProc cntg", 20, "\r");#->1982
					
				}
				return ($cntgCountDataHsh_inThr_ref, $wholeGenomeGCDataHsh_inThr_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $cntgCountDataHsh_ref = {};
	my $wholeGenomeGCDataHsh_ref = {};

	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($cntgCountDataHsh_inThr_ref, $wholeGenomeGCDataHsh_inThr_ref) = $thr->join;

				foreach my $GCPct (keys %{$wholeGenomeGCDataHsh_inThr_ref}) {
					push @{$wholeGenomeGCDataHsh_ref->{$GCPct}}, @{$wholeGenomeGCDataHsh_inThr_ref->{$GCPct}};
				}

				foreach my $cntg (keys %{$cntgCountDataHsh_inThr_ref}) { 
					$cntgCountDataHsh_ref->{$cntg}{'cntgLength'} = $cntgCountDataHsh_inThr_ref->{$cntg}{'cntgLength'};
					foreach my $strnd (qw/+ -/) {
						$cntgCountDataHsh_ref->{$cntg}{'covPct'}{$strnd} = $cntgCountDataHsh_inThr_ref->{$cntg}{'covPct'}{$strnd};
						$cntgCountDataHsh_ref->{$cntg}{'covPos'}{$strnd} = $cntgCountDataHsh_inThr_ref->{$cntg}{'covPos'}{$strnd};
						$cntgCountDataHsh_ref->{$cntg}{'rdCount'}{$strnd} = $cntgCountDataHsh_inThr_ref->{$cntg}{'rdCount'}{$strnd};
					}
				}
				
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	return ($cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref);

}
sub countCtgryOnCntg {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateThreadHshWithRandomCntg|1273, reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 6_countThePileupStorable|181
#	secondaryAppearInSection: >none
#	input: $fullPileupStorablePathHsh_ref, $geneByCntgHsh_ref, $geneInfoHsh_ref, $margin3End, $margin5End, $maxThread, $midPtPileupStorablePathHsh_ref, $minCov, $resultStorableDir
#	output: $ctgryPosCountHsh_ref, $ctgryReadCountHsh_ref, $ctgryTotalPosHsh_ref, $nonStructuralCount, $totalCount
#	toCall: my ($ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $nonStructuralCount) = &countCtgryOnCntg($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $margin5End, $margin3End, $minCov);
#	calledInLine: 190
#....................................................................................................................................................#
	my ($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $margin5End, $margin3End, $minCov) = @_;
	
	#---define zero for all catgry
	my @cntgAry = (keys %{$midPtPileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->1273
	my $cntgProc :shared = 0;
	my %threadHsh = ();
	my $ambiguousCtgryCountHsh_ref = {};


	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 0, "\n");#->1982
		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781, the 
		
			sub {
				my ($cntgAry_ref) = @_;

				my $ctgryReadCountInThrHsh_ref = {};
				$ctgryReadCountInThrHsh_ref->{$geneInfoHsh_ref->{$_}{'ctgry'}} = 0 foreach (keys %{$geneInfoHsh_ref});
				$ctgryReadCountInThrHsh_ref->{'unannotated'} = 0;

				my $ctgryPosCountInThrHsh_ref = {};
				$ctgryPosCountInThrHsh_ref->{$geneInfoHsh_ref->{$_}{'ctgry'}} = 0 foreach (keys %{$geneInfoHsh_ref});
				$ctgryPosCountInThrHsh_ref->{'unannotated'} = 0;

				my $ctgryTotalPosInThrHsh_ref = {};
				$ctgryTotalPosInThrHsh_ref->{'unannotated'} = 0;
				
				foreach my $cntg (@{$cntgAry_ref}) {
					my $annotatedPos = 0;
					$cntgProc++;
					my $midPtCntgCovAry_ref = retrieve($midPtPileupStorablePathHsh_ref->{$cntg});
					my $fullCntgCovAry_ref = retrieve($fullPileupStorablePathHsh_ref->{$cntg});

					my %cntgAnnoHsh = ();
		
					if (exists $geneByCntgHsh_ref->{$cntg}) {
						foreach my $geneID (keys %{$geneByCntgHsh_ref->{$cntg}}) {
						
							my ($startMargin, $endMargin);
							if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
								($startMargin, $endMargin) = ($margin5End, $margin3End);
							} else {
								($startMargin, $endMargin) = ($margin3End, $margin5End);
							}
							
							foreach my $i (${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[0]-1-$startMargin..${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[1]-1+$endMargin) {
								if (not $cntgAnnoHsh{$i}) {
									$cntgAnnoHsh{$i} = $geneInfoHsh_ref->{$geneID}{'ctgry'};
								} else  {
									#---ambiguous if this location belongs two different categories
									if ($cntgAnnoHsh{$i} ne $geneInfoHsh_ref->{$geneID}{'ctgry'}) {
										$cntgAnnoHsh{$i} = 'ambiguous';
									}
								}
							}
						}
					}
					
					#---get the annoPos [29/09/2013 18:06]
					foreach my $i (keys %cntgAnnoHsh) {
						my $ctgry = $cntgAnnoHsh{$i};
						$annotatedPos++;
						$ctgryTotalPosInThrHsh_ref->{$ctgry}++;
					}

					my $unannotatedPos = @{$fullCntgCovAry_ref} - $annotatedPos;
					$ctgryTotalPosInThrHsh_ref->{'unannotated'} += $unannotatedPos;

					for my $i (0..$#{$midPtCntgCovAry_ref}) {
						if ($fullCntgCovAry_ref->[$i]) {
							my $ctgry = 'unannotated';
							$ctgry = $cntgAnnoHsh{$i} if $cntgAnnoHsh{$i};
							$ctgryPosCountInThrHsh_ref->{$ctgry}++ if sum((split /,/, $fullCntgCovAry_ref->[$i])) >= $minCov;
							if ($midPtCntgCovAry_ref->[$i]) {
								$ctgryReadCountInThrHsh_ref->{$ctgry} += sum((split /,/, $midPtCntgCovAry_ref->[$i]));
							}
						}
					}

					&reportStatus("$cntgProc cntg counted", 20, "\r");#->1982

				}

				#print $_."\n" foreach (keys %{$ctgryReadCountInThrHsh_ref});
				
				return ($ctgryReadCountInThrHsh_ref, $ctgryPosCountInThrHsh_ref, $ctgryTotalPosInThrHsh_ref);
			}
			,($cntgAry_ref)
		);
	}

	#---wait and collect ctgryReadCountInThrHsh_ref into ctgryReadCountHsh_ref
	my $ctgryReadCountHsh_ref = {};
	$ctgryReadCountHsh_ref->{$geneInfoHsh_ref->{$_}{'ctgry'}} = 0 foreach (keys %{$geneInfoHsh_ref});
	$ctgryReadCountHsh_ref->{'unannotated'} = 0;

	my $ctgryPosCountHsh_ref = {};
	$ctgryPosCountHsh_ref->{$geneInfoHsh_ref->{$_}{'ctgry'}} = 0 foreach (keys %{$geneInfoHsh_ref});
	$ctgryPosCountHsh_ref->{'unannotated'} = 0;
	
	my $ctgryTotalPosHsh_ref = {};
	$ctgryTotalPosHsh_ref->{$geneInfoHsh_ref->{$_}{'ctgry'}} = 0 foreach (keys %{$geneInfoHsh_ref});
	$ctgryTotalPosHsh_ref->{'unannotated'} = 0;

	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($ctgryReadCountInThrHsh_ref, $ctgryPosCountInThrHsh_ref, $ctgryTotalPosInThrHsh_ref) = $thr->join;
				$ctgryReadCountHsh_ref->{$_} += $ctgryReadCountInThrHsh_ref->{$_} foreach (keys %{$ctgryReadCountInThrHsh_ref});
				$ctgryPosCountHsh_ref->{$_} += $ctgryPosCountInThrHsh_ref->{$_} foreach (keys %{$ctgryPosCountInThrHsh_ref});
				$ctgryTotalPosHsh_ref->{$_} += $ctgryTotalPosInThrHsh_ref->{$_} foreach (keys %{$ctgryTotalPosInThrHsh_ref});
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	my $totalCount = 0;
	$totalCount += $ctgryReadCountHsh_ref->{$_} foreach (keys %{$ctgryReadCountHsh_ref});
	my $nonStructuralCount = $totalCount;
	eval {
		$nonStructuralCount -= $ctgryReadCountHsh_ref->{'rRNA'};
		$nonStructuralCount -= $ctgryReadCountHsh_ref->{'tRNA'};
	};

	&reportStatus("totalCount = $totalCount", 20, "\n");#->1982
	&reportStatus("nonStructuralCount = $nonStructuralCount", 20, "\n");#->1982

	store($ctgryReadCountHsh_ref, "$resultStorableDir/ctgryReadCountHsh.pls");
	store($ctgryPosCountHsh_ref, "$resultStorableDir/ctgryPosCountHsh.pls");
	store($ctgryTotalPosHsh_ref, "$resultStorableDir/ctgryTotalPosHsh.pls");

	return ($ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $nonStructuralCount);
}
sub countGeneOnCntg {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: countIndividualGeneRng|829, generateThreadHshWithRandomCntg|1273, reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 6_countThePileupStorable|181
#	secondaryAppearInSection: >none
#	input: $fullPileupStorablePathHsh_ref, $geneByCntgHsh_ref, $geneInfoHsh_ref, $margin3End, $margin5End, $maxThread, $midPtPileupStorablePathHsh_ref, $minCov, $resultStorableDir
#	output: $ctgryDirtnPosCountHsh_ref, $ctgryDirtnReadCountHsh_ref, $ctgryGeneRngSumHsh_ref, $geneIndivCovIdxPathHsh_ref, $geneIndivDataHsh_ref
#	toCall: my ($geneIndivDataHsh_ref, $ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryGeneRngSumHsh_ref, $geneIndivCovIdxPathHsh_ref) = &countGeneOnCntg($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneInfoHsh_ref, $geneByCntgHsh_ref, $maxThread, $margin5End, $margin3End, $resultStorableDir, $minCov);
#	calledInLine: 186
#....................................................................................................................................................#
	my ($midPtPileupStorablePathHsh_ref, $fullPileupStorablePathHsh_ref, $geneInfoHsh_ref, $geneByCntgHsh_ref, $maxThread, $margin5End, $margin3End, $resultStorableDir, $minCov) = @_;
	
	#---[12/11/2013 14:28] define the rngTypes to be counted
	my @rngTypeAry = (qw/exonRng CDSRng geneRng UTR5Rng UTR3Rng/);

	#---[12/11/2013 14:25] store the individual gene positions and coverage as storables for later use
	my @mkDirAry = ();
	my $indivGeneStroableDir = "$resultStorableDir/indivGene/";
	push @mkDirAry, $indivGeneStroableDir;
	foreach my $rngType (@rngTypeAry) {
		push @mkDirAry, "$indivGeneStroableDir/$rngType/";
	}
	foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

	my @cntgAry = (keys %{$midPtPileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->1273
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1982

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;

				my $geneIndivDataHsh_inThr_ref = {};
				my $geneIndivCovCntgHsh_inThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $midPtCntgCovAry_ref = retrieve($midPtPileupStorablePathHsh_ref->{$cntg});
					my $fullCntgCovAry_ref = retrieve($fullPileupStorablePathHsh_ref->{$cntg});
					my $indivGeneCovInCntgHsh_ref = {};
					
					next if not exists $geneByCntgHsh_ref->{$cntg};
					foreach my $geneID (keys %{$geneByCntgHsh_ref->{$cntg}}) {
						foreach my $rngType (@rngTypeAry) {
							my ($geneRngDataHsh_ref, $geneRngCovHsh_ref) = &countIndividualGeneRng($midPtCntgCovAry_ref, $fullCntgCovAry_ref, $rngType, $geneInfoHsh_ref, $geneID, $margin5End, $margin3End, $minCov);#->829

							if ($geneRngCovHsh_ref->{'idx'}) {
								$indivGeneCovInCntgHsh_ref->{$rngType}{$geneID} = $geneRngCovHsh_ref;
							}

							foreach my $dataType (keys %{$geneRngDataHsh_ref}) { 
								foreach my $dirtn (keys %{$geneRngDataHsh_ref->{$dataType}}) {
									$geneIndivDataHsh_inThr_ref->{$dataType}{$rngType}{$geneID}{$dirtn} = $geneRngDataHsh_ref->{$dataType}{$dirtn};
								}
							}
						}
					}
					
					&reportStatus("Finished counting features on $cntgProc cntg", 20, "\r");#->1982
					
					foreach my $rngType (keys %{$indivGeneCovInCntgHsh_ref}) {
						my $storableName = "$cntg.hsh.pls";
						my $absStorablePath = "$indivGeneStroableDir/$rngType/$storableName";
						store($indivGeneCovInCntgHsh_ref->{$rngType}, $absStorablePath);
						$geneIndivCovCntgHsh_inThr_ref->{$rngType}{$cntg} = $storableName;
					}
				}
				return ($geneIndivDataHsh_inThr_ref, $geneIndivCovCntgHsh_inThr_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $geneIndivDataHsh_ref = {};
	my $geneIndivCovCntgHsh_ref = {};
	my $geneIndivCovIdxPathHsh_ref = {};

	my $ctgryDirtnReadCountHsh_ref = {};
	my $ctgryDirtnPosCountHsh_ref = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($geneIndivDataHsh_inThr_ref, $geneIndivCovCntgHsh_inThr_ref) = $thr->join;

				foreach my $rngType (keys %{$geneIndivCovCntgHsh_inThr_ref}) {
					foreach my $cntg (keys %{$geneIndivCovCntgHsh_inThr_ref->{$rngType}}) { 
						$geneIndivCovCntgHsh_ref->{$rngType}{$cntg} = $geneIndivCovCntgHsh_inThr_ref->{$rngType}{$cntg};
					}
				}
				
				foreach my $dataType (keys %{$geneIndivDataHsh_inThr_ref}) { 
					foreach my $rngType (keys %{$geneIndivDataHsh_inThr_ref->{$dataType}}) {
						foreach my $geneID (keys %{$geneIndivDataHsh_inThr_ref->{$dataType}{$rngType}}) {
							foreach my $dirtn (keys %{$geneIndivDataHsh_inThr_ref->{$dataType}{$rngType}{$geneID}}) {
								$geneIndivDataHsh_ref->{$dataType}{$rngType}{$geneID}{$dirtn} = $geneIndivDataHsh_inThr_ref->{$dataType}{$rngType}{$geneID}{$dirtn};
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	#---[12/11/2013 18:23] store the geneIndivCovCntgHsh_ref in geneIndivCovIdxHshPath
	foreach my $rngType (@rngTypeAry) {
		my $dir = "$indivGeneStroableDir/$rngType/";

		if ($geneIndivCovCntgHsh_ref->{$rngType}) {
			my $geneIndivCovIdxHshPath = "$dir/index.hsh.pls";
			store($geneIndivCovCntgHsh_ref->{$rngType}, $geneIndivCovIdxHshPath);
			$geneIndivCovIdxPathHsh_ref->{$rngType} = $geneIndivCovIdxHshPath;
		} else {
			system "rm -fr $dir";
		}
	}
	
	#---[29/09/2013 21:05] get ctgry sense and antisense count 
	my %allGeneCtgryHsh = ();
	$allGeneCtgryHsh{$geneInfoHsh_ref->{$_}{'ctgry'}}++ foreach (keys %{$geneInfoHsh_ref});
	
	#---[29/09/2013 21:05] set the zero
	my $ctgryGeneRngSumHsh_ref = {};
	foreach my $ctgry (keys %allGeneCtgryHsh) {
		$ctgryGeneRngSumHsh_ref->{$ctgry} = 0;
		foreach my $dirtn (qw /a s/) {
			$ctgryDirtnReadCountHsh_ref->{$ctgry}{$dirtn} = 0;
			$ctgryDirtnPosCountHsh_ref->{$ctgry}{$dirtn} = 0;
		}
	}
	
	foreach my $geneID (keys %{$geneIndivDataHsh_ref->{'covPos'}{'geneRng'}}) {
		foreach my $dirtn (keys %{$geneIndivDataHsh_ref->{'covPos'}{'geneRng'}{$geneID}}) {
			$ctgryDirtnReadCountHsh_ref->{$geneInfoHsh_ref->{$geneID}{'ctgry'}}{$dirtn} += $geneIndivDataHsh_ref->{'readCount'}{'geneRng'}{$geneID}{$dirtn} if ($geneIndivDataHsh_ref->{'readCount'}{'geneRng'}{$geneID}{$dirtn} > 0);
			$ctgryDirtnPosCountHsh_ref->{$geneInfoHsh_ref->{$geneID}{'ctgry'}}{$dirtn} += $geneIndivDataHsh_ref->{'covPos'}{'geneRng'}{$geneID}{$dirtn} if ($geneIndivDataHsh_ref->{'covPos'}{'geneRng'}{$geneID}{$dirtn} > 0);
		}
	}
	
	#---[10/10/2013 14:13] get the gene rng some to calculate the dirtn percentage
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my @geneRngAry = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'geneRng'}};
		my $geneLength = $geneRngAry[-1]-$geneRngAry[0]+1;
		$ctgryGeneRngSumHsh_ref->{$ctgry} += $geneLength;
	}

	return ($geneIndivDataHsh_ref, $ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryGeneRngSumHsh_ref, $geneIndivCovIdxPathHsh_ref);
}
sub countIndividualGeneRng {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: countGeneOnCntg|668
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_countThePileupStorable|181
#	input: $fullCntgCovAry_ref, $geneID, $geneInfoHsh_ref, $margin3End, $margin5End, $midPtCntgCovAry_ref, $minCov, $rngType
#	output: $geneRngCovHsh_ref, $geneRngDataHsh_ref
#	toCall: my ($geneRngDataHsh_ref, $geneRngCovHsh_ref) = &countIndividualGeneRng($midPtCntgCovAry_ref, $fullCntgCovAry_ref, $rngType, $geneInfoHsh_ref, $geneID, $margin5End, $margin3End, $minCov);
#	calledInLine: 721
#....................................................................................................................................................#
	my ($midPtCntgCovAry_ref, $fullCntgCovAry_ref, $rngType, $geneInfoHsh_ref, $geneID, $margin5End, $margin3End, $minCov) = @_;
	
	my $geneRngDataHsh_ref = {};
	my $geneRngCovHsh_ref = {};
	my $cntgAryIndexAry_ref = [];
	
	my %strndToDirtnHsh = ();
	$strndToDirtnHsh{'+'}{'+'} = 's';
	$strndToDirtnHsh{'+'}{'-'} = 'a';
	$strndToDirtnHsh{'-'}{'+'} = 'a';
	$strndToDirtnHsh{'-'}{'-'} = 's';

	my $geneStrnd = $geneInfoHsh_ref->{$geneID}{'strnd'};
	
	foreach my $dataType (qw/readCount covPos covSum/) {
		foreach my $dirtn (qw/a s/) {
			$geneRngDataHsh_ref->{$dataType}{$dirtn} = -1;
		}
	}
	
	if (exists $geneInfoHsh_ref->{$geneID}{$rngType}) {

		foreach my $dataType (keys %{$geneRngDataHsh_ref}) {
			foreach my $dirtn (keys %{$geneRngDataHsh_ref->{$dataType}}) {
				$geneRngDataHsh_ref->{$dataType}{$dirtn} = 0;
			}
		}
		
		my @posAry = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{$rngType}};

		my ($startMargin, $endMargin);
		if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
			($startMargin, $endMargin) = ($margin5End, $margin3End);
		} else {
			($startMargin, $endMargin) = ($margin3End, $margin5End);
		}
		
		#---extend the margin if > 0;
		$posAry[0] -= $startMargin if $startMargin > 0;
		$posAry[-1] += $endMargin if $endMargin > 0;
		
		for (my $i=0; $i < $#posAry; $i += 2) {
			foreach my $cntgAryIndex ($posAry[$i]-1..$posAry[$i+1]-1) {
				push @{$cntgAryIndexAry_ref}, $cntgAryIndex;
			}
		}
		
		#---trim the margins if <0
		if ($startMargin < 0) {
			my $marginToTrim = -1*$startMargin;
			@{$cntgAryIndexAry_ref} = @{$cntgAryIndexAry_ref->[$marginToTrim..$#{$cntgAryIndexAry_ref}]};
		}

		if ($endMargin < 0) {
			my $marginToTrim = -1*$endMargin;
			@{$cntgAryIndexAry_ref} = @{$cntgAryIndexAry_ref->[0..$#{$cntgAryIndexAry_ref}-$marginToTrim]};
		}
		
		foreach my $cntgAryIndex (sort {$a <=> $b} @{$cntgAryIndexAry_ref}) {
			
			if ($fullCntgCovAry_ref->[$cntgAryIndex]) {
				my %tmpCovHsh = ();
				($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $fullCntgCovAry_ref->[$cntgAryIndex];
				foreach my $dataStrnd (qw/+ -/) {
					my $dirtn = $strndToDirtnHsh{$geneStrnd}{$dataStrnd};
					$tmpCovHsh{$dataStrnd} = 0 if not $tmpCovHsh{$dataStrnd};
					$geneRngDataHsh_ref->{'covSum'}{$dirtn} += $tmpCovHsh{$dataStrnd} if $tmpCovHsh{$dataStrnd} > 0;
					$geneRngDataHsh_ref->{'covPos'}{$dirtn} ++ if $tmpCovHsh{$dataStrnd} >= $minCov;
					push @{$geneRngCovHsh_ref->{'cov'}{$dirtn}}, $tmpCovHsh{$dataStrnd};
				}
				push @{$geneRngCovHsh_ref->{'idx'}}, $cntgAryIndex;
			}

			if ($midPtCntgCovAry_ref->[$cntgAryIndex]) {
				my %tmpCovHsh = ();
				($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $midPtCntgCovAry_ref->[$cntgAryIndex];
				foreach my $dataStrnd (qw/+ -/) {
					my $dirtn = $strndToDirtnHsh{$geneStrnd}{$dataStrnd};
					$tmpCovHsh{$dataStrnd} = 0 if not $tmpCovHsh{$dataStrnd};
					$geneRngDataHsh_ref->{'readCount'}{$dirtn} += $tmpCovHsh{$dataStrnd} if $tmpCovHsh{$dataStrnd} > 0;
				}
			}
		}
	}
	
	return ($geneRngDataHsh_ref, $geneRngCovHsh_ref);
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|1607, readGFF_oneRNAPerGene|1786, reportStatus|1982
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|91, 10_finishingTasks|230, 4_processInputData|149
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 1627, 1630, 1635, 1806, 1998
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateGeneByCntgHsh {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: $geneByCntgHsh_ref
#	toCall: my ($geneByCntgHsh_ref) = &generateGeneByCntgHsh($geneInfoHsh_ref);
#	calledInLine: 160
#....................................................................................................................................................#
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Generating geneByCntgHsh", 0, "\n");#->1982

	my $geneByCntgHsh_ref = {};
	$geneByCntgHsh_ref->{$geneInfoHsh_ref->{$_}{'cntg'}}{$_}++ foreach (keys %{$geneInfoHsh_ref});

	return ($geneByCntgHsh_ref);
}
sub generateMasterHTML {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 9_outputHTML|220
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $resultDir, $resultHTMLDir
#	output: none
#	toCall: &generateMasterHTML($allResultFilePathHsh_ref, $resultHTMLDir, $resultDir);
#	calledInLine: 225
#....................................................................................................................................................#

	my ($allResultFilePathHsh_ref, $resultHTMLDir, $resultDir) = @_;

	#$allResultFilePathHsh_ref->{'General_Statisitcs'}{'cntgCountData'} = $cntgCountDataPath;
	#$allResultFilePathHsh_ref->{'GC_vs_Coverage'}{'xls_Table'}{'wholeGenomeRng'}{'wholeGenomeData'} = $wholeGenomeGCDataPath;
	#$allResultFilePathHsh_ref->{'metaGene_Coverage'}{'xls_Table'}{$rngType}{$ctgry}{$dirtn} = $xlsPath;
	#$allResultFilePathHsh_ref->{'readCount_on_Genes'}{'Overall_Pooled_for_Differential_Expression'}{$dirtnType} = $overallXlsPathHsh{$dirtnType};

	my $masterIndexHTML = "$resultDir/master_index.html";

	my %htmlPathHsh = ();
	foreach my $resultType (keys %{$allResultFilePathHsh_ref}) {
		$htmlPathHsh{$resultType}{'abs'} = "$resultHTMLDir/$resultType.html";
		$htmlPathHsh{$resultType}{'rel'} = File::Spec->abs2rel($htmlPathHsh{$resultType}{'abs'}, $resultDir);
	}

	{#---[19/11/2013 16:13] General statistics
		my $resultType = 'General_Statisitcs';
		my $title = "General statisitcs";
		open (HTML, ">", $htmlPathHsh{$resultType}{'abs'});
		print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
		print HTML join "", ('<html>', "\n");
		print HTML join "", ("<head><title>$title</title></head>\n");
		print HTML join "", ('<body>', "\n");
		print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
		print HTML join "", ('</div><div>', "\n");
		print HTML join "", ("<h4>Click the following links to access the relevant information on $title</h4>\n");
		print HTML join "", ('<ul>', "\n");
		foreach my $dataType (sort keys %{$allResultFilePathHsh_ref->{$resultType}}) {
			my $absPath = $allResultFilePathHsh_ref->{$resultType}{$dataType};
			my $relPath = File::Spec->abs2rel($absPath, $resultHTMLDir);
			print HTML "<li><a href=\'$relPath\'>$dataType</a>\n";
		}
		print HTML join "", ('</body>', "\n");
		print HTML join "", ('</html>', "\n");
		close HTML;
	}

	{#---[19/11/2013 16:13] GC_vs_Coverage
		my $resultType = 'GC_vs_Coverage';
		my $title = "GC vs Coverage";
		open (HTML, ">", $htmlPathHsh{$resultType}{'abs'});
		print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
		print HTML join "", ('<html>', "\n");
		print HTML join "", ("<head><title>$title</title></head>\n");
		print HTML join "", ('<body>', "\n");
		print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
		print HTML join "", ('</div><div>', "\n");
		foreach my $xlsOrPdf (sort keys %{$allResultFilePathHsh_ref->{$resultType}}) {
			print HTML join "", ("<h2>$xlsOrPdf</h2>\n");
			print HTML join "", ('<ul>', "\n");
			foreach my $rng (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}}) {
				print HTML join "", ("<h4>range $rng</h4>\n");
				print HTML join "", ('<ul>', "\n");
				foreach my $ctgry (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}{$rng}}) {
					my $absPath = $allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}{$rng}{$ctgry};
					my $relPath = File::Spec->abs2rel($absPath, $resultHTMLDir);
					print HTML "<li><a href=\'$relPath\'>category $ctgry</a>\n";
				}
				print HTML join "", ('</ul>', "\n");
			}
			print HTML join "", ('</ul>', "\n");
		}
		print HTML join "", ('</body>', "\n");
		print HTML join "", ('</html>', "\n");
		close HTML;
	}

	{#---[19/11/2013 16:13] metaGene coverage plot
		my $resultType = 'metaGene_Coverage';
		my $title = "metaGene coverage plot";
		open (HTML, ">", $htmlPathHsh{$resultType}{'abs'});
		print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
		print HTML join "", ('<html>', "\n");
		print HTML join "", ("<head><title>$title</title></head>\n");
		print HTML join "", ('<body>', "\n");
		print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
		print HTML join "", ('</div><div>', "\n");
		foreach my $xlsOrPdf (sort keys %{$allResultFilePathHsh_ref->{$resultType}}) {
			print HTML join "", ("<h2>$xlsOrPdf</h2>\n");
			print HTML join "", ('<ul>', "\n");
			foreach my $rng (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}}) {
				print HTML join "", ("<h4>range $rng</h4>\n");
				print HTML join "", ('<ul>', "\n");
				foreach my $ctgry (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}{$rng}}) {
					foreach my $dirtn (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}{$rng}{$ctgry}}) {
						my $absPath = $allResultFilePathHsh_ref->{$resultType}{$xlsOrPdf}{$rng}{$ctgry}{$dirtn};
						my $relPath = File::Spec->abs2rel($absPath, $resultHTMLDir);
						print HTML "<li><a href=\'$relPath\'>category $ctgry direction $dirtn</a>\n";
					}
				}
				print HTML join "", ('</ul>', "\n");
			}
			print HTML join "", ('</ul>', "\n");
		}
		print HTML join "", ('</body>', "\n");
		print HTML join "", ('</html>', "\n");
		close HTML;
	}

	{#---[19/11/2013 16:13] metaGene coverage plot
		my $resultType = 'readCount_on_Genes';
		my $title = "read Count on Genes";
		open (HTML, ">", $htmlPathHsh{$resultType}{'abs'});
		print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
		print HTML join "", ('<html>', "\n");
		print HTML join "", ("<head><title>$title</title></head>\n");
		print HTML join "", ('<body>', "\n");
		print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
		print HTML join "", ('</div><div>', "\n");
		foreach my $rng (sort keys %{$allResultFilePathHsh_ref->{$resultType}}) {
			print HTML join "", ("<h4>range $rng</h4>\n");
			print HTML join "", ('<ul>', "\n");
			foreach my $ctgry (sort keys %{$allResultFilePathHsh_ref->{$resultType}{$rng}}) {
				my $absPath = $allResultFilePathHsh_ref->{$resultType}{$rng}{$ctgry};
				my $relPath = File::Spec->abs2rel($absPath, $resultHTMLDir);
				print HTML "<li><a href=\'$relPath\'>category $ctgry</a>\n";
			}
			print HTML join "", ('</ul>', "\n");
		}
		print HTML join "", ('</body>', "\n");
		print HTML join "", ('</html>', "\n");
		close HTML;
	}
	
	{
		my $title = "Data of perlStorablePileupCounter";
		open (HTML, ">", $masterIndexHTML);
		print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
		print HTML join "", ('<html>', "\n");
		print HTML join "", ("<head><title>Summary of differential expression analyses</title></head>\n");
		print HTML join "", ('<body>', "\n");
		print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
		print HTML join "", ('</div><div>', "\n");
		print HTML join "", ("<h4>Click the following links to access the relevant information</h4>\n");
		print HTML join "", ('<ul>', "\n");
		
		foreach my $resultType (keys %{$allResultFilePathHsh_ref}) {
			print HTML "<li><a href=\'$htmlPathHsh{$resultType}{'rel'}\'>$resultType</a>\n";
		}
		print HTML join "", ('</body>', "\n");
		print HTML join "", ('</html>', "\n");
		close HTML;
	}
}
sub generateMetaGeneCoverageAndGCBias {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: generateThreadHshWithRandomCntg|1273, getIndivCntgCovPlsPath|1300, reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 6_countThePileupStorable|181
#	secondaryAppearInSection: >none
#	input: $GCWinSize, $fastaHsh_ref, $geneIndivCovIdxPathHsh_ref, $geneInfoHsh_ref, $maxThread, $metaGeneMaxLen, $metaGeneMinCov, $metaGeneMinLen
#	output: $GCPctCovHsh_ref, $metaGeneLenPctCovSumHsh_ref, $metaGeneNumCountHsh_ref
#	toCall: my ($metaGeneLenPctCovSumHsh_ref, $GCPctCovHsh_ref, $metaGeneNumCountHsh_ref) = &generateMetaGeneCoverageAndGCBias($geneIndivCovIdxPathHsh_ref, $metaGeneMinLen, $metaGeneMaxLen, $metaGeneMinCov, $geneInfoHsh_ref, $GCWinSize, $fastaHsh_ref, $maxThread);
#	calledInLine: 188
#....................................................................................................................................................#
	my ($geneIndivCovIdxPathHsh_ref, $metaGeneMinLen, $metaGeneMaxLen, $metaGeneMinCov, $geneInfoHsh_ref, $GCWinSize, $fastaHsh_ref, $maxThread) = @_;
	
	my $metaGeneLenPctCovSumHsh_ref = {};
	my $metaGeneNumCountHsh_ref = {};
	my $GCPctCovHsh_ref = {};

	foreach my $rngType (keys %{$geneIndivCovIdxPathHsh_ref}) {

		my $geneIndivCovIdxHshPath = $geneIndivCovIdxPathHsh_ref->{$rngType};
		my ($geneIndivCovPlsPathHsh_ref) = &getIndivCntgCovPlsPath($geneIndivCovIdxHshPath);#->1300
		&reportStatus("Generating metagene plot for $rngType", 10, "\n");#->1982
		
		my @cntgAry = (keys %{$geneIndivCovPlsPathHsh_ref});
		my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->1273
		my $geneProc :shared = 0;
		my %threadHsh = ();

		foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
			my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
			my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
			&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1982

			#---spawn a new thread
			($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
				sub {
					my ($cntgAry_ref) = @_;
					my $metaGeneNumCountHsh_inThr_ref = {};
					my $metaGeneLenPctCovSumHsh_inThr_ref = {};
					my $GCPctCovByGeneHsh_inThr_ref = {};

					foreach my $cntg (@{$cntgAry_ref}) {
						my $indivGeneCovHsh_ref = retrieve($geneIndivCovPlsPathHsh_ref->{$cntg});
						system "gzip -f $geneIndivCovPlsPathHsh_ref->{$cntg}";
						&reportStatus("$geneProc genes processed", 10, "\r");#->1982

						foreach my $geneID (keys %{$indivGeneCovHsh_ref}) {
							$geneProc++;
							my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
							my $strnd = $geneInfoHsh_ref->{$geneID}{'strnd'};

							#---[12/11/2013 22:41] generate GC bias data
							my %tmpBaseCountHsh = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
							my $covSum = 0;
							foreach my $i (0..$#{$indivGeneCovHsh_ref->{$geneID}{'idx'}}) {
								my $idx = $indivGeneCovHsh_ref->{$geneID}{'idx'}->[$i];
								my $base = substr $fastaHsh_ref->{$cntg}, $idx, 1;
								$tmpBaseCountHsh{$base}++;
								$covSum += $indivGeneCovHsh_ref->{$geneID}{'cov'}{'s'}->[$i];

								if ($i % $GCWinSize == 0 and $i > 0) {
									my $GCPct = sprintf "%.0f", 100*($tmpBaseCountHsh{'C'} + $tmpBaseCountHsh{'G'})/$GCWinSize;
									my $covPerNt = sprintf "%.2f", $covSum/$GCWinSize;
									push @{$GCPctCovByGeneHsh_inThr_ref->{$geneID}{$GCPct}}, $covPerNt;
									%tmpBaseCountHsh = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
									$covSum = 0;
								}
							}
							
							#---[12/11/2013 22:41] generate metaGene lenPct vs covPct data
							my %tmpCovHsh = ();
							if ($strnd eq '+') {
								@{$tmpCovHsh{'s'}} = @{$indivGeneCovHsh_ref->{$geneID}{'cov'}{'s'}};
								@{$tmpCovHsh{'a'}} = @{$indivGeneCovHsh_ref->{$geneID}{'cov'}{'a'}};
							} elsif ($strnd eq '-') {
								@{$tmpCovHsh{'s'}} = reverse @{$indivGeneCovHsh_ref->{$geneID}{'cov'}{'s'}};
								@{$tmpCovHsh{'a'}} = reverse @{$indivGeneCovHsh_ref->{$geneID}{'cov'}{'a'}};
							} else {
								die "invalid gene strand\n";
							}
							
							my $rngLen = @{$tmpCovHsh{'s'}};
							next if ($rngLen < $metaGeneMinLen or $rngLen > $metaGeneMaxLen);
							foreach my $dirtn (keys %tmpCovHsh) {
								my $maxCov = max(@{$tmpCovHsh{$dirtn}});
								next if ($maxCov < $metaGeneMinCov);
								my %tmpLenCovPctHsh = ();
								foreach my $i (0..$#{$tmpCovHsh{$dirtn}}) {
									my $lenPct = sprintf "%.0f", 100*($i+1)/$rngLen;
									$lenPct++ if $lenPct == 0; #---[15/11/2013 10:51] force the min lenPct to be 1
									my $covPct = sprintf "%.2f", 100*$tmpCovHsh{$dirtn}->[$i]/$maxCov;
									push @{$tmpLenCovPctHsh{$lenPct}}, $covPct;
								}

								#---[12/11/2013 21:22] store the count
								$metaGeneNumCountHsh_inThr_ref->{$ctgry}{$dirtn}++;
				
								foreach my $lenPct (keys %tmpLenCovPctHsh) {
									my $avgCovPct = sprintf "%.2f", sum(@{$tmpLenCovPctHsh{$lenPct}})/@{$tmpLenCovPctHsh{$lenPct}};
									push @{$metaGeneLenPctCovSumHsh_inThr_ref->{$ctgry}{$dirtn}{$lenPct}}, $avgCovPct;
								}#---[12/11/2013 21:19] end of foreach my $lenPct (keys %{$tmpLenCovPctHsh})
							}#---[12/11/2013 21:20]  end of foreach my $dirtn (keys %{$tmpCovHsh})
						}#---[12/11/2013 21:20] end of foreach my $geneID (keys %{$indivGeneCovHsh_ref})
					}#---[12/11/2013 21:20] end of foreach my $cntg (keys %{$geneIndivCovPlsPathHsh_ref})
					return ($metaGeneLenPctCovSumHsh_inThr_ref, $metaGeneNumCountHsh_inThr_ref, $GCPctCovByGeneHsh_inThr_ref);
				}
				,($cntgAry_ref)
			);
		}#---[12/11/2013 22:27] end of foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref});
		
		while (keys %threadHsh) {
			foreach my $threadNum (keys %threadHsh) {
				my $thr = $threadHsh{$threadNum};
				if (not $thr->is_running()) {
					my ($metaGeneLenPctCovSumHsh_inThr_ref, $metaGeneNumCountHsh_inThr_ref, $GCPctCovByGeneHsh_inThr_ref) = $thr->join;

					foreach my $geneID (keys %{$GCPctCovByGeneHsh_inThr_ref}) {
						my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
						foreach my $GCPct (keys %{$GCPctCovByGeneHsh_inThr_ref->{$geneID}}) {
							push @{$GCPctCovHsh_ref->{'gene'}{$rngType}{$geneID}{$GCPct}}, @{$GCPctCovByGeneHsh_inThr_ref->{$geneID}{$GCPct}};
							push @{$GCPctCovHsh_ref->{'ctgry'}{$rngType}{$ctgry}{$GCPct}}, @{$GCPctCovByGeneHsh_inThr_ref->{$geneID}{$GCPct}};
						}
					}
					
					foreach my $ctgry (keys %{$metaGeneNumCountHsh_inThr_ref}) {
						foreach my $dirtn (keys %{$metaGeneNumCountHsh_inThr_ref->{$ctgry}}) {
							$metaGeneNumCountHsh_ref->{$rngType}{$ctgry}{$dirtn} = 0 if not $metaGeneNumCountHsh_ref->{$rngType}{$ctgry}{$dirtn};
							$metaGeneNumCountHsh_ref->{$rngType}{$ctgry}{$dirtn} += $metaGeneNumCountHsh_inThr_ref->{$ctgry}{$dirtn};
						}
					}

					foreach my $ctgry (keys %{$metaGeneLenPctCovSumHsh_inThr_ref}) {
						foreach my $dirtn (keys %{$metaGeneLenPctCovSumHsh_inThr_ref->{$ctgry}}) {
							foreach my $lenPct (keys %{$metaGeneLenPctCovSumHsh_inThr_ref->{$ctgry}{$dirtn}}) {
								push @{$metaGeneLenPctCovSumHsh_ref->{$rngType}{$ctgry}{$dirtn}{$lenPct}}, @{$metaGeneLenPctCovSumHsh_inThr_ref->{$ctgry}{$dirtn}{$lenPct}};
							}
						}
					}
					delete $threadHsh{$threadNum};
				}
			}
			sleep 1;
		}
	}#---[12/11/2013 21:20] end of foreach my $rngType (keys %{$geneIndivCovIdxPathHsh_ref})

	return ($metaGeneLenPctCovSumHsh_ref, $GCPctCovHsh_ref, $metaGeneNumCountHsh_ref);
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: countContigReadAndGC|396, countCtgryOnCntg|517, countGeneOnCntg|668, generateMetaGeneCoverageAndGCBias|1124
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_countContigBasedInfo|171, 6_countThePileupStorable|181
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 411, 532, 694, 1148
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return $randCntgInThreadHsh_ref;

}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|1982
#	appearInSub: generateMetaGeneCoverageAndGCBias|1124
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: 6_countThePileupStorable|181
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 162, 163, 1144
#....................................................................................................................................................#
	
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
		system ("gzip -fd $cntgCovPlsPath.gz") if -s "$cntgCovPlsPath.gz";
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	&reportStatus("pls path of $numCntg contig stored.", 0, "\n");#->1982
	
	return $cntgCovInPlsPathHsh_ref;
}
sub ggplotXYLinesMultipleSamples {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: outputContigReadAndGC|1369, outputGCPctData|1455, outputMetaGeneCovData|1534
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_metagene and GC plot|196
#	input: $RScriptPath, $XAXis, $YAxis, $YVariable, $dataPath, $extraArg, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: none
#	toCall: &ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);
#	calledInLine: 1445, 1526, 1600
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($YVariable, $YAxis, $XAXis)), "\n";
	foreach my $YCategory (sort keys %{$plotDataHsh_ref}) {
		foreach my $XVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$YCategory}}) {
			my $YVal = $plotDataHsh_ref->{$YCategory}{$XVal};
			print PLOTDATA join "", (join "\t", ($YCategory, $YVal, $XVal)), "\n";
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, colour=$YVariable)) + geom_line() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

}
sub outputContigReadAndGC {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: ggplotXYLinesMultipleSamples|1334
#	appearInSub: >none
#	primaryAppearInSection: 7_metagene and GC plot|196
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $cntgCountDataHsh_ref, $ggplotDirHsh_ref, $nonStructuralCount, $resultContigCountDir, $totalCount, $wholeGenomeGCDataHsh_ref
#	output: 
#	toCall: &outputContigReadAndGC($resultContigCountDir, $cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 204
#....................................................................................................................................................#
	my ($resultContigCountDir, $cntgCountDataHsh_ref, $wholeGenomeGCDataHsh_ref, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref) = @_;
	
	my @mkDirAry = ();
	my @ggplotFileTypeAry = qw /dat pdf R log/;
	my $sub_ggplotDirHsh_ref = {};
	foreach my $fileType (@ggplotFileTypeAry) {$sub_ggplotDirHsh_ref->{$fileType} = $ggplotDirHsh_ref->{$fileType}."/contigCount/"; push @mkDirAry, $sub_ggplotDirHsh_ref->{$fileType};}
	foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

	my $cntgCountDataPath = "$resultContigCountDir/cntgCountData.xls";
	$allResultFilePathHsh_ref->{'General_Statisitcs'}{'Count_Coverage_on_Contig'} = $cntgCountDataPath;
	
	open CNTGDATA, ">", $cntgCountDataPath;
	print CNTGDATA join "", (join "\t", (qw/cntg cntgLength covPct_+ covPos_+ rdCount_+ covPct_- covPos_- rdCount_-/)), "\n";
	foreach my $cntg (keys %{$cntgCountDataHsh_ref}) {
		my $cntgLength = $cntgCountDataHsh_ref->{$cntg}{'cntgLength'};
		my @outputAry = ($cntg, $cntgLength);
		foreach my $strnd (qw/+ -/) {
			my $covPct = $cntgCountDataHsh_ref->{$cntg}{'covPct'}{$strnd};
			my $covPos = $cntgCountDataHsh_ref->{$cntg}{'covPos'}{$strnd};
			my $rdCount = $cntgCountDataHsh_ref->{$cntg}{'rdCount'}{$strnd};
			push @outputAry, ($covPct, $covPos, $rdCount);
		}
		print CNTGDATA join "", (join "\t", (@outputAry)), "\n";
	}
	close CNTGDATA;

	my $nameTag = "wholeGenomeGCData";
	my $wholeGenomeGCDataPath = "$resultContigCountDir/$nameTag.xls";
	$allResultFilePathHsh_ref->{'GC_vs_Coverage'}{'xls_Table'}{'wholeGenomeRng'}{'wholeGenomeData'} = $wholeGenomeGCDataPath;

	open WHOLEGNGC, ">", $wholeGenomeGCDataPath;
	my $plotDataHsh_ref = {};
	print WHOLEGNGC join "", (join "\t", (qw/GCPct covMean covMeanPerTotalCount covMeanPerNSCount covSD covPct95 covPct75 covPct50 covPct25 covPct5 numSam/)), "\n";
	foreach my $GCPct (sort {$a <=> $b} keys %{$wholeGenomeGCDataHsh_ref}) {
		my $numSam = @{$wholeGenomeGCDataHsh_ref->{$GCPct}};
		if ($numSam >= 20) {
			my $statObj = Statistics::Descriptive::Full->new();
			$statObj->add_data(@{$wholeGenomeGCDataHsh_ref->{$GCPct}});
			my $covMean = $statObj->mean();
			my $covSD = $statObj->standard_deviation();
			my $covPct95 = $statObj->percentile(95);
			my $covPct75 = $statObj->percentile(75);
			my $covPct50 = $statObj->percentile(50);
			my $covPct25 = $statObj->percentile(25);
			my $covPct5 = $statObj->percentile(5);
			my $covMeanPerTotalCount = sprintf "%.5f", $covMean/($totalCount/1000000);
			my $covMeanPerNSCount = sprintf "%.5f", $covMean/($nonStructuralCount/1000000);
			print WHOLEGNGC join "", (join "\t", ($GCPct, $covMean, $covMeanPerTotalCount, $covMeanPerNSCount, $covSD, $covPct95, $covPct75, $covPct50, $covPct25, $covPct5, $numSam)), "\n";
			$plotDataHsh_ref->{'covMean'}{$GCPct} = $covMean;
			$plotDataHsh_ref->{'covPct75'}{$GCPct} = $covPct75;
			$plotDataHsh_ref->{'covPct25'}{$GCPct} = $covPct25;
		}
	}

	my $pdfPath = $sub_ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
	$allResultFilePathHsh_ref->{'GC_vs_Coverage'}{'pdf_Figure'}{'wholeGenomeRng'}{'wholeGenomeData'} = $pdfPath;
	my $dataPath = $sub_ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
	my $RScriptPath = $sub_ggplotDirHsh_ref->{'R'}."/$nameTag.R";
	my $logPath = $sub_ggplotDirHsh_ref->{'log'}."/$nameTag.log";
	my $xAxis = 'GCPercentage';
	my $YAxis = 'coverage';
	my $YVariable = 'covInterval';
	my $extraArg = '';
	my $height = 6;
	my $width = 14;
	&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1334

	print WHOLEGNGC join "", (join "\t", ("\ntotalCount = $totalCount")), "\n";
	print WHOLEGNGC join "", (join "\t", ("nonStructuralCount = $nonStructuralCount")), "\n";

	close WHOLEGNGC;

	return ();
}
sub outputGCPctData {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: ggplotXYLinesMultipleSamples|1334, reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 7_metagene and GC plot|196
#	secondaryAppearInSection: >none
#	input: $GCPctCovHsh_ref, $allResultFilePathHsh_ref, $ggplotDirHsh_ref, $nonStructuralCount, $resultGCPctDir, $totalCount
#	output: 
#	toCall: &outputGCPctData($GCPctCovHsh_ref, $resultGCPctDir, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 202
#....................................................................................................................................................#
	my ($GCPctCovHsh_ref, $resultGCPctDir, $totalCount, $nonStructuralCount, $ggplotDirHsh_ref, $allResultFilePathHsh_ref) = @_;

	&reportStatus("Outputing the GCPct data", 10, "\n");#->1982

	my @ggplotFileTypeAry = qw /dat pdf R log/;
	
	foreach my $rngType (keys %{$GCPctCovHsh_ref->{'ctgry'}}) {
		my @mkDirAry = ();
		my $xlsDir = "$resultGCPctDir/$rngType/"; push @mkDirAry, $xlsDir;
		my $sub_ggplotDirHsh_ref = {};
		foreach my $fileType (@ggplotFileTypeAry) {$sub_ggplotDirHsh_ref->{$fileType} = $ggplotDirHsh_ref->{$fileType}."/GCPct/$rngType/"; push @mkDirAry, $sub_ggplotDirHsh_ref->{$fileType};}
		foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

		foreach my $ctgry (keys %{$GCPctCovHsh_ref->{'ctgry'}{$rngType}}) {
			my $nameTag = "$rngType.GCPct.$ctgry";
			my $xlsPath = "$xlsDir/$nameTag.xls";

			$allResultFilePathHsh_ref->{'GC_vs_Coverage'}{'xls_Table'}{$rngType}{$ctgry} = $xlsPath;

			my $plotDataHsh_ref = {};
			open (GCOUT, ">", $xlsPath);
			print GCOUT join "", (join "\t", (qw/GCPct covMean covMeanPerTotalCount covMeanPerNSCount covSD covPct95 covPct75 covPct50 covPct25 covPct5 numSam/)), "\n";
			foreach my $GCPct (sort {$a <=> $b} keys %{$GCPctCovHsh_ref->{'ctgry'}{$rngType}{$ctgry}}) {
				my $numSam = @{$GCPctCovHsh_ref->{'ctgry'}{$rngType}{$ctgry}{$GCPct}};
				if ($numSam >= 20) {
					my $statObj = Statistics::Descriptive::Full->new();
					$statObj->add_data(@{$GCPctCovHsh_ref->{'ctgry'}{$rngType}{$ctgry}{$GCPct}});
					my $covMean = $statObj->mean();
					my $covSD = $statObj->standard_deviation();
					my $covPct95 = $statObj->percentile(95);
					my $covPct75 = $statObj->percentile(75);
					my $covPct50 = $statObj->percentile(50);
					my $covPct25 = $statObj->percentile(25);
					my $covPct5 = $statObj->percentile(5);
					my $covMeanPerTotalCount = sprintf "%.5f", $covMean/($totalCount/1000000);
					my $covMeanPerNSCount = sprintf "%.5f", $covMean/($nonStructuralCount/1000000);
					print GCOUT join "", (join "\t", ($GCPct, $covMean, $covMeanPerTotalCount, $covMeanPerNSCount, $covSD, $covPct95, $covPct75, $covPct50, $covPct25, $covPct5, $numSam)), "\n";
					$plotDataHsh_ref->{'covMean'}{$GCPct} = $covMean;
					$plotDataHsh_ref->{'covPct75'}{$GCPct} = $covPct75;
					$plotDataHsh_ref->{'covPct25'}{$GCPct} = $covPct25;
				}
			}

			print GCOUT join "", (join "\t", ("\ntotalCount = $totalCount")), "\n";
			print GCOUT join "", (join "\t", ("nonStructuralCount = $nonStructuralCount")), "\n";
			close GCOUT;
			
			my $pdfPath = $sub_ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			
			$allResultFilePathHsh_ref->{'GC_vs_Coverage'}{'pdf_Figure'}{$rngType}{$ctgry} = $pdfPath;
			
			my $dataPath = $sub_ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $RScriptPath = $sub_ggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $sub_ggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $xAxis = 'GCPercentage';
			my $YAxis = 'coverage';
			my $YVariable = 'covInterval';
			my $extraArg = '';
			my $height = 6;
			my $width = 14;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1334
			
		}
	}

	return ();
}
sub outputMetaGeneCovData {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: ggplotXYLinesMultipleSamples|1334
#	appearInSub: >none
#	primaryAppearInSection: 7_metagene and GC plot|196
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $ggplotDirHsh_ref, $metaGeneLenPctCovSumHsh_ref, $metaGeneMaxLen, $metaGeneMinLen, $resultLenPctDir
#	output: 
#	toCall: &outputMetaGeneCovData($metaGeneLenPctCovSumHsh_ref, $resultLenPctDir, $metaGeneMinLen, $metaGeneMaxLen, $ggplotDirHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 203
#....................................................................................................................................................#
	my ($metaGeneLenPctCovSumHsh_ref, $resultLenPctDir, $metaGeneMinLen, $metaGeneMaxLen, $ggplotDirHsh_ref, $allResultFilePathHsh_ref) = @_;
	
	my @ggplotFileTypeAry = qw /dat pdf R log/;

	foreach my $rngType (keys %{$metaGeneLenPctCovSumHsh_ref}) {
		foreach my $ctgry (keys %{$metaGeneLenPctCovSumHsh_ref->{$rngType}}) {
			my @mkDirAry = ();
			my $xlsDir = "$resultLenPctDir/$rngType/$ctgry/"; push @mkDirAry, $xlsDir;
			my $sub_ggplotDirHsh_ref = {};
			foreach my $fileType (@ggplotFileTypeAry) {$sub_ggplotDirHsh_ref->{$fileType} = $ggplotDirHsh_ref->{$fileType}."/metaGenePlot/$rngType/$ctgry/"; push @mkDirAry, $sub_ggplotDirHsh_ref->{$fileType};}
			foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
			my $nameTag = "$rngType.$ctgry.from.$metaGeneMinLen.to.$metaGeneMaxLen";
			my $plotDataHsh_ref = {};

			foreach my $dirtn (qw/a s/) {
				my $xlsPath = "$xlsDir/$nameTag.$dirtn.xls";
		
				$allResultFilePathHsh_ref->{'metaGene_Coverage'}{'xls_Table'}{$rngType}{$ctgry}{$dirtn} = $xlsPath;

				open (METAGENEOUT, ">", $xlsPath);
				print METAGENEOUT join "", (join "\t", (qw/lenPct covMean covSD covPct95 covPct75 covPct50 covPct25 covPct5 numSam/)), "\n";
				foreach my $lenPct (1..100) {
					my @tmpAry = (0);#---[15/11/2013 11:05] set a default value as 0 to prevent crash if the array is absent
					@tmpAry = @{$metaGeneLenPctCovSumHsh_ref->{$rngType}{$ctgry}{$dirtn}{$lenPct}} if exists $metaGeneLenPctCovSumHsh_ref->{$rngType}{$ctgry}{$dirtn}{$lenPct};
					my $numSam = @tmpAry;
					my $statObj = Statistics::Descriptive::Full->new();
					$statObj->add_data(@tmpAry);
					my $covMean = $statObj->mean();
					my $covSD = my $covPct95 = my $covPct75 = my $covPct50 = my $covPct25 = my $covPct5 = -999999;
					if ($numSam >= 20) {
						$covSD = $statObj->standard_deviation();
						$covPct95 = $statObj->percentile(95);
						$covPct75 = $statObj->percentile(75);
						$covPct50 = $statObj->percentile(50);
						$covPct25 = $statObj->percentile(25);
						$covPct5 = $statObj->percentile(5);
					}
					print METAGENEOUT join "", (join "\t", ($lenPct, $covMean, $covSD, $covPct95, $covPct75, $covPct50, $covPct25, $covPct5, $numSam)), "\n";
					$plotDataHsh_ref->{$dirtn}{$lenPct} = $covMean;
				}
			}
			
			my $pdfPath = $sub_ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";

			$allResultFilePathHsh_ref->{'metaGene_Coverage'}{'pdf_Figure'}{$rngType}{$ctgry}{'both_Directions'} = $pdfPath;

			my $dataPath = $sub_ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $RScriptPath = $sub_ggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $sub_ggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $xAxis = 'scaledGeneGength';
			my $YAxis = 'scaledReadCoverage';
			my $YVariable = 'direction';
			my $extraArg = '';
			my $height = 6;
			my $width = 14;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->1334
		}
	}

	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|928
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|91, 10_finishingTasks|230
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 97, 235
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->928
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->928
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->928
		print "=========================================================================\n\n";
	}
}
sub printFeatureCount {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 8_printStatisticsLog|209
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $countCovDataHsh_ref, $ctgryDirtnReadCountHsh_ref, $geneInfoHsh_ref, $geneRngLenHsh_ref, $resultGeneCountDir
#	output: 
#	toCall: &printFeatureCount($resultGeneCountDir, $geneInfoHsh_ref, $ctgryDirtnReadCountHsh_ref, $countCovDataHsh_ref, $geneRngLenHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 215
#....................................................................................................................................................#
	my ($resultGeneCountDir, $geneInfoHsh_ref, $ctgryDirtnReadCountHsh_ref, $countCovDataHsh_ref, $geneRngLenHsh_ref, $allResultFilePathHsh_ref) = @_;

	&reportStatus('Printing individual feature statistics', 10, "\n");#->1982
	
	system ("mkdir -pm 777 $resultGeneCountDir/overall/");

	my %overallXlsPathHsh = (
		'exonRng_both_dirtn'=>"$resultGeneCountDir/overall/overall.exonRng.count.both.dirtn.xls",
		'exonRng_sense_dirtn'=>"$resultGeneCountDir/overall/overall.exonRng.count.senseOnly.dirtn.xls",
		'exonRng_antisense_dirtn'=>"$resultGeneCountDir/overall/overall.exonRng.count.antisenseOnly.dirtn.xls",
	);

	open (OVERALLBOTHCOUNT, ">", $overallXlsPathHsh{'exonRng_both_dirtn'});
	open (OVERALLSENSECOUNT, ">", $overallXlsPathHsh{'exonRng_sense_dirtn'});
	open (OVERALLANTISENSECOUNT, ">", $overallXlsPathHsh{'exonRng_antisense_dirtn'});
	
	foreach my $dirtnType (keys %overallXlsPathHsh) {
		$allResultFilePathHsh_ref->{'readCount_on_Genes'}{'Overall_Pooled_for_Differential_Expression'}{$dirtnType} = $overallXlsPathHsh{$dirtnType};
	}

	my %FHHsh = ();
	foreach my $ctgry (keys %{$ctgryDirtnReadCountHsh_ref}) {
		my $subOutDir = "$resultGeneCountDir/$ctgry/";
		system ("mkdir -pm 777 $subOutDir");
		foreach my $rngType (keys %{$countCovDataHsh_ref}) {
			my $xlsPath = "$subOutDir/$rngType.$ctgry.xls";
			open $FHHsh{$ctgry}{$rngType}, ">", $xlsPath;
			$allResultFilePathHsh_ref->{'readCount_on_Genes'}{$rngType}{$ctgry} = $xlsPath;
			print {$FHHsh{$ctgry}{$rngType}} join "", ((join "\t", ('geneID', 'description', 'location', 'rngLength', 'covPct_s', 'covPct_a', 'covPerNt_s', 'covPerNt_a', 'covPos_s', 'covPos_a', 'covSum_s', 'covSum_a', 'readCount_s', 'readCount_a', 'RPKM_s', 'RPKM_a', 'RPM_s', 'RPM_a', 'UQNor_s', 'UQNor_a')), "\n");
		}
	}
	
	foreach my $rngType (keys %{$countCovDataHsh_ref}) {
		foreach my $geneID (keys %{$countCovDataHsh_ref->{$rngType}}) {
			my $description = $geneInfoHsh_ref->{$geneID}{'description'};
			my $location = $geneInfoHsh_ref->{$geneID}{'cntg'}.":".${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[0]."-".${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[1];
			my $rngLength = $geneRngLenHsh_ref->{$rngType}{$geneID};
			my @outputAry = ($geneID, $description, $location, $rngLength);
			foreach my $dataType (sort {lc $a cmp lc $b} keys %{$countCovDataHsh_ref->{$rngType}{$geneID}}) {
				foreach my $dirtn (sort {$b cmp $a} keys %{$countCovDataHsh_ref->{$rngType}{$geneID}{$dataType}}) {
					push @outputAry, $countCovDataHsh_ref->{$rngType}{$geneID}{$dataType}{$dirtn};
				}
			}
			print {$FHHsh{$geneInfoHsh_ref->{$geneID}{'ctgry'}}{$rngType}} join "", ((join "\t", @outputAry), "\n");
			if ($rngType eq 'exonRng') {
				foreach my $dirtn (qw/a s/) {
					my $RNAID = $dirtn."_".$geneID;
					my $count = $countCovDataHsh_ref->{$rngType}{$geneID}{'readCount'}{$dirtn};
					print OVERALLBOTHCOUNT join "", (join "\t", ($RNAID, $count, $geneID, $geneInfoHsh_ref->{$geneID}{'ctgry'})), "\n";

					if ($dirtn eq 's') {
						print OVERALLSENSECOUNT join "", (join "\t", ($geneID, $count, $geneID, $geneInfoHsh_ref->{$geneID}{'ctgry'})), "\n";
					}
					if ($dirtn eq 'a') {
						print OVERALLANTISENSECOUNT join "", (join "\t", ($geneID, $count, $geneID, $geneInfoHsh_ref->{$geneID}{'ctgry'})), "\n";
					}
				}
			}
		}
	}
	
	return ();
}
sub printReadStatistics {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 8_printStatisticsLog|209
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryDirtnReadCountHsh_ref, $ctgryGeneRngSumHsh_ref, $ctgryPosCountHsh_ref, $ctgryReadCountHsh_ref, $ctgryTotalPosHsh_ref, $minCov, $resultGenStatDir, $totalCount
#	output: 
#	toCall: &printReadStatistics($ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $ctgryGeneRngSumHsh_ref, $resultGenStatDir, $minCov, $allResultFilePathHsh_ref);
#	calledInLine: 214
#....................................................................................................................................................#
	my ($ctgryDirtnReadCountHsh_ref, $ctgryDirtnPosCountHsh_ref, $ctgryReadCountHsh_ref, $ctgryPosCountHsh_ref, $ctgryTotalPosHsh_ref, $totalCount, $ctgryGeneRngSumHsh_ref, $resultGenStatDir, $minCov, $allResultFilePathHsh_ref) = @_;
	
	&reportStatus('Printing read ctgry statistics', 10, "\n");#->1982
	
	my $ctgryCountSepStrandXlsPath = "$resultGenStatDir/category.SenseVsAntisense.count.minCov.$minCov.coverage.xls";
	$allResultFilePathHsh_ref->{'General_Statisitcs'}{'Category_ReadCount_on_Separate_Strands'} = $ctgryCountSepStrandXlsPath;
	
	open CTGRYSA, ">", $ctgryCountSepStrandXlsPath;
	print CTGRYSA join "", (join "\t", ('ctgry', 's_readCount', 'a_readCount', 'S_A_readCount_ratio', 's_covPct', 'a_covPct', 's_covPos', 'a_covPos', 'totalPos')), "\n";
	foreach my $ctgry (sort keys %{$ctgryDirtnReadCountHsh_ref}) {
		my $S_A_readCount_ratio = -1;
		my $totalPos = $ctgryGeneRngSumHsh_ref->{$ctgry};
		my $s_readCount = $ctgryDirtnReadCountHsh_ref->{$ctgry}{'s'};
		my $a_readCount = $ctgryDirtnReadCountHsh_ref->{$ctgry}{'a'};
		my $s_covPos = $ctgryDirtnPosCountHsh_ref->{$ctgry}{'s'};
		my $a_covPos = $ctgryDirtnPosCountHsh_ref->{$ctgry}{'a'};
		my $s_covPct = sprintf "%.5f", 100*$s_covPos/$totalPos;
		my $a_covPct = sprintf "%.5f", 100*$a_covPos/$totalPos;
		eval {$S_A_readCount_ratio = sprintf "%.5f", $s_readCount/$a_readCount;};
		print CTGRYSA join "", (join "\t", ($ctgry, $s_readCount, $a_readCount, $S_A_readCount_ratio, $s_covPct, $a_covPct, $s_covPos, $a_covPos, $totalPos)), "\n";
	}
	close CTGRYSA;

	my $ctgryCountPooledStrandXlsPath = "$resultGenStatDir/category.pooled.strand.absolute.percentage.xls";
	$allResultFilePathHsh_ref->{'General_Statisitcs'}{'Category_Total_ReadCount_on_Pooled_Strands'} = $ctgryCountPooledStrandXlsPath;

	open CTGRYPCT, ">", $ctgryCountPooledStrandXlsPath;
	print CTGRYPCT join "", (join "\t", ('ctgry', 'count', 'pct')), "\n";
	foreach my $ctgry (sort keys %{$ctgryReadCountHsh_ref}) {
		my $pct = sprintf "%.5f", 100*$ctgryReadCountHsh_ref->{$ctgry}/$totalCount;
		print CTGRYPCT join "", (join "\t", ($ctgry, $ctgryReadCountHsh_ref->{$ctgry}, $pct)), "\n";
	}
	print CTGRYPCT join "", (join "\t", ('total', $totalCount, '100')), "\n";
	close CTGRYPCT;
	
	my $ctgryCovPooledStrandXlsPath = "$resultGenStatDir/category.minCov.$minCov.coverage.xls";
	$allResultFilePathHsh_ref->{'General_Statisitcs'}{'Category_Coverage_on_Pooled_Strands'} = $ctgryCovPooledStrandXlsPath;

	open CTGRYCOV, ">", $ctgryCovPooledStrandXlsPath;
	print CTGRYCOV join "", (join "\t", ('ctgry', 'coveredPos', 'totalPos', 'coveredPct', 'genomePct')), "\n";
	#---get genomeSize [29/09/2013 18:54]
	my $genomeSize = 0;
	$genomeSize += $ctgryTotalPosHsh_ref->{$_} foreach (keys %{$ctgryTotalPosHsh_ref});

	foreach my $ctgry (sort keys %{$ctgryPosCountHsh_ref}) {
		my $coveredPos = $ctgryPosCountHsh_ref->{$ctgry};
		my $totalPos = $ctgryTotalPosHsh_ref->{$ctgry};
		my $genomePct = sprintf "%.5f", 100*$totalPos/$genomeSize;
		my $coveredPct = sprintf "%.5f", 100*$coveredPos/$totalPos;
		print CTGRYCOV join "", (join "\t", ($ctgry, $coveredPos, $totalPos, $coveredPct, $genomePct)), "\n";
	}
	my $genomeCovered = 0;
	$genomeCovered += $ctgryPosCountHsh_ref->{$_} foreach (keys %{$ctgryPosCountHsh_ref});
	my $coveredPct = sprintf "%.5f", 100*$genomeCovered/$genomeSize;
	print CTGRYCOV join "", (join "\t", ('total', $genomeCovered, $genomeSize, $coveredPct, '100')), "\n";
	close CTGRYCOV;
	
	return ();
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|928
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 158
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->928
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 156
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->1982
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/ *\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|91
#	secondaryAppearInSection: >none
#	input: none
#	output: $fastaPath, $fullPileupStorableIndexPath, $gffPath, $margin3End, $margin5End, $maxThread, $midPtPileupStorableIndexPath, $minCov, $outDir
#	toCall: my ($midPtPileupStorableIndexPath, $fullPileupStorableIndexPath, $gffPath, $fastaPath, $maxThread, $margin5End, $margin3End, $minCov, $outDir) = &readParameters();
#	calledInLine: 98
#....................................................................................................................................................#
	
	my ($midPtPileupStorableIndexPath, $fullPileupStorableIndexPath, $gffPath, $fastaPath, $maxThread, $margin5End, $margin3End, $minCov, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/pileupPerlStorableCounter/";
	$maxThread = 4;
	$margin5End = 0;
	$margin3End = 0;
	$minCov = 1;
	
	GetOptions 	("midPtPileupStorableIndexPath=s"  => \$midPtPileupStorableIndexPath,
				 "fullPileupStorableIndexPath=s" => \$fullPileupStorableIndexPath,
				 "gffPath=s"  => \$gffPath,
				 "fastaPath=s"  => \$fastaPath,
				 "maxThread:i"  => \$maxThread,
				 "margin5End:i"  => \$margin5End,
				 "margin3End:i"  => \$margin3End,
				 "minCov:i"  => \$minCov,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($gffPath, $fastaPath, $midPtPileupStorableIndexPath, $fullPileupStorableIndexPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($midPtPileupStorableIndexPath, $fullPileupStorableIndexPath, $gffPath, $fastaPath, $maxThread, $margin5End, $margin3End, $minCov, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|928
#	appearInSub: calculateNormalizedCount|290, checkGeneInfo|368, countContigReadAndGC|396, countCtgryOnCntg|517, countGeneOnCntg|668, generateGeneByCntgHsh|946, generateMetaGeneCoverageAndGCBias|1124, getIndivCntgCovPlsPath|1300, outputGCPctData|1455, printFeatureCount|1640, printReadStatistics|1715, readMultiFasta|1885, zipUnzipCntgCovInPlsPathHsh|2003
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|149, 5_countContigBasedInfo|171, 6_countThePileupStorable|181, 7_metagene and GC plot|196, 8_printStatisticsLog|209
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 303, 382, 390, 409, 418, 476, 541, 611, 658, 659, 701, 735, 959, 1145, 1155, 1169, 1329, 1468, 1653, 1728, 1903, 2018
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->928

	return ();
}
sub zipUnzipCntgCovInPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|1982
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|149
#	secondaryAppearInSection: >none
#	input: $cntgCovInPlsPathHsh_ref, $zipUnzip
#	output: none
#	toCall: &zipUnzipCntgCovInPlsPathHsh($zipUnzip, $cntgCovInPlsPathHsh_ref);
#	calledInLine: 164, 165
#....................................................................................................................................................#

	my ($zipUnzip, $cntgCovInPlsPathHsh_ref) = @_;
	
	foreach my $cntg (sort keys %{$cntgCovInPlsPathHsh_ref}) {
		&reportStatus("Trying to $zipUnzip cntg ary", 20, "\r");#->1982

		my $cntgCovPlsPath = "$cntgCovInPlsPathHsh_ref->{$cntg}";

		if ($zipUnzip eq 'unzip') {
			system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
		} else {
			system ("gzip -f $cntgCovPlsPath") if (-s "$cntgCovPlsPath");
		}
	}
	print "\n";
}

exit;
