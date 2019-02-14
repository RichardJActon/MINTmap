#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Parallel::Loops;
use Digest::MD5;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0) . '/lib';
use MINTmap::subs;
#############################################################################
### Options
#############################################################################

# Argument declarations
#######################################

my $help;
my $tRNAannoFile;
my $tRNAseqFile;
my $chrDir;
my $maxProcs = 2;#25
my $minLen = 16;
my $maxLen = 50;

my $subStr1 = 85;
my $subStr2 = -85;

my $tRNAseqSpaceFile = "./tRNAseqSpace.fa";
my $OtherAnnoFile = "./otherAnno.txt";
my $tRNAfragsFile = "./tRNAfrags.txt";
my $tRNAsubStrannoFile = "./tRNAsubStrannoFile.bed";

# my @chr = (1..22,"X","Y","M");
my @chr = ("T"); ##!!

# Argument getting
#######################################

GetOptions(
	"tRNAannoFile|i=s"		=> \$tRNAannoFile,
	"tRNAseqFile|s=s"		=> \$tRNAseqFile,
	"chrDir|c=s"			=> \$chrDir,
	"maxProcs|n=s"			=> \$maxProcs,
	"minLen=s"				=> \$minLen,
	"maxLen=s"				=> \$maxLen,
	"subStr1=s"				=> \$subStr1,
	"subStr2=s"				=> \$subStr2,
	"help|?" 				=> \$help
	) or pod2usage(1);

my $pl = Parallel::Loops -> new($maxProcs);

pod2usage(1) if $help;

#############################################################################
### body
#############################################################################

# considerations
# - splice order~
# - incl CCA?
# - size of upstream / downstream overlaps to allow

# Data Read-in
#######################################

# 100bp flank tRNA data hash ~ 
# hg19-tRNA-100bp53.bed
print STDERR "[".getTime()."] Reading tRNA annotation...\n";
my @tRNAannoLines;
open TRNAANNO, "$tRNAannoFile" or die "unable to open $tRNAannoFile";
while (my $line = <TRNAANNO>) {
	chomp $line;
	push @tRNAannoLines, $line;
}
close TRNAANNO or die "unable to close $tRNAannoFile";
print STDERR "[".getTime()."] Read tRNA annotation\n";
# 100bp flank fasta
# hg19-tRNA-100bp53.fa

my %tRNAseqs;
my $tRNAname;
print STDERR "[".getTime()."] Reading tRNA Sequences...\n";
open TRNASEQ, "$tRNAseqFile" or die "unable to open $tRNAseqFile";
while (my $line = <TRNASEQ>) {
	chomp $line;
	if ($line =~ /^>(.+)::/) {
		$tRNAname = $1;
		next;
	} else {
		$tRNAseqs{$tRNAname} .= uc $line;
	}
}
close TRNASEQ or die "unable to close $tRNAseqFile";
print STDERR "[".getTime()."] Read tRNA Sequences\n";

my %tRNAs = parsetRNAanno(\@tRNAannoLines);

foreach my $k (keys %tRNAs) { # reciprocal
	if (defined $tRNAseqs{$k}) {
		$tRNAs{$k}{seq} = $tRNAseqs{$k};
	} else {
		warn "$k is missing a sequence\n";
	}
	
}

# print STDOUT Dumper %tRNAs;
# print STDOUT Dumper %tRNAseqs;

# Substring generation 
#######################################

my $tRNAseqSpaceFasta;
print STDERR "[".getTime()."] Generating tRNA substrings...\n";
my @substrings;
foreach my $tRNAname (keys %tRNAs) {
	my $str;
	if (defined $tRNAs{$tRNAname}{seq}) {
		$str = substr $tRNAs{$tRNAname}{seq}, $subStr1, $subStr2; # introns...
		$tRNAseqSpaceFasta .= ">$tRNAname\n$str\n";
		## generate tRNAspaceSeqs.fa # MD5
	} else { 
		next 
	}
	# print "$k = $str\n";
	push @substrings, @{getAllSubstrings($minLen,$maxLen,$str)};
}
my %substr = map {$_ => ''} @substrings;
print STDERR "[".getTime()."] Generated tRNA substrings\n";

# Substring printing
#######################################

#tRNAsubStrTrans
# print STDERR Dumper %{tRNAsubStrAnno(\%tRNAs,$subStr1,$subStr2)};
# print STDERR Dumper tRNAsubStrTrans($tRNAs{'tRNA-Tyr-GTA-chr1-127'});

tRNAsubStrSave(\%tRNAs,$subStr1,$subStr2,$tRNAsubStrannoFile);

# Substring printing
#######################################

print STDERR "[".getTime()."] Saving tRNA substrings...\n";
my $tRNAseqsMD5 = Digest::MD5::md5_hex($tRNAseqSpaceFasta);
open(my $tsfh, ">$tRNAseqSpaceFile") or die $!;
#my $tRNAseqsMD5 = Digest::MD5->new->addfile($tsfh)->hexdigest;
print $tsfh $tRNAseqSpaceFasta;
close $tsfh or die $!;
print STDERR "[".getTime()."] Saved tRNA substrings\n";
# print Dumper @substrings;
# print Dumper %substr;

# Search Genome for unique substrings
#######################################

print STDERR "[".getTime()."] Matching ...\n";
#my @chr = ("H100k6"); ##!!
my %matches;
$pl -> share(\%matches);
$pl -> foreach(\@chr, sub {
	my $fwd = readChr($_, $chrDir);
	my $rev = transcribe($fwd);
	#my %matches;
	foreach my $k (keys %substr) {
		my @fwd = @{getCoords($k,$fwd)};
		#print "$k:\n$fwd\n";
		my @rev = @{getCoords($k,$rev)};
		if (scalar @fwd >= 1) {
			$matches{$k}{$_}{fwd} = \@fwd;
		} else {
			$matches{$k}{$_}{fwd} = "NA"
		}

		if (scalar @rev >= 1) {
			$matches{$k}{$_}{rev} = \@rev;
		} else {
			$matches{$k}{$_}{rev} = "NA"
		}
	}
});

print STDERR "[".getTime()."] Matched \n";
# print Dumper %matches;

# print matches
#######################################

# map {
# 	print STDOUT join "\t", @{$_};
# 	print STDOUT "\n";
# } @{matchesToBed(\%matches)};
print STDERR "[".getTime()."] Formatting Matched ...\n";
my $matchesString;
map {
	$matchesString .= join "\t", @{$_};
	$matchesString .= "\n";
} @{matchesToBed(\%matches)};
print STDERR "[".getTime()."] Formatted Matched \n";
#print Dumper @{matchesToBed(\%matches)};
# print STDOUT $matchesString;  

# intersect with tRNA gene regions here or with bedtools
print STDERR "[".getTime()."] Intersecting ...\n";
my @bedFrags = readpipe join("",
	"echo \"$matchesString\" | ",
	"bedtools intersect ",
	"-a stdin ",
	"-b $tRNAsubStrannoFile ", #$tRNAannoFile #..\/tmp-tRNAlookupTest.bed
	"-f 1 -c"
);
print STDERR "[".getTime()."] Intersected \n";

open TMP, ">./bedinterC.bed" or die $!;
map {print TMP $_} @bedFrags;
print TMP "\n";
close TMP or die $!;
# Annotate strings
#######################################
## generate other annotations
# 5' 3' etc - defer
print STDERR "[".getTime()."] Annotating ... \n";
# if >0 tRNA space exclusive?
# Y = 1, N = 0 or >1
# >1 in additional in other annotation?
#my @otherAnno;
my $otherAnnoStr;
map {
	my $line = $_;
	chomp $line;
	my @cols = split "\t", $line;
	if ($cols[6] == 1) {
	 	#push @otherAnno, join "\t", $cols[3], "Unique";
	 	$otherAnnoStr .= join("\t", $cols[3], "Unique"). "\n";
	} elsif ($cols[6] == 0) {
		#push @otherAnno, join "\t", $cols[3], "Non-unique";
		$otherAnnoStr .= join("\t", $cols[3], "Non-unique"). "\n";
	} else {
		#push @otherAnno, join "\t", $cols[3], "tRNAspace_multi-mapping";
		$otherAnnoStr .= join("\t", $cols[3], "tRNAspace_multi-mapping"). "\n";
	}
} @bedFrags;

print STDERR "[".getTime()."] Annotated ...\n";

# # print STDERR Dumper @otherAnno;

# print Annotated strings
#######################################

print STDERR "[".getTime()."] Saving Annotated \n";
my $OtherAnnoMD5 = Digest::MD5::md5_hex($otherAnnoStr);
open(my $oafh, ">$OtherAnnoFile") or die $!;
# my $OtherAnnoMD5 = Digest::MD5->new->addfile($oafh)->hexdigest;
# print $oafh join "\n", @otherAnno;
# print $oafh "\n";
print $oafh $otherAnnoStr;
close $oafh or die $!;
print STDERR "[".getTime()."] Saved Annotated \n";

# set exclusive values
#######################################
print STDERR "[".getTime()."] Marking Unique ...\n";
#my @tRNAfrags;
my $tRNAfrags;
# $pl -> share(\@tRNAfrags);
# $pl -> foreach(\@bedFrags, sub 
map {
	my $line = $_;
	chomp $line;
	my @cols = split "\t", $line;
	if ($cols[6] == 1) {
	 	#push @tRNAfrags, join "\t", $cols[3], "Y";
	 	$tRNAfrags .= join("\t", $cols[3], "Y"). "\n";
	} else {
		#push @tRNAfrags, join "\t", $cols[3], "N";
		$tRNAfrags .= join("\t", $cols[3], "N"). "\n";
	}
} @bedFrags;#);
print STDERR "[".getTime()."] Marked Unique \n";

# print fragments
#######################################
print STDERR "[".getTime()."] Saving Marked ...\n";
open(my $tffh, ">$tRNAfragsFile") or die $!;
print $tffh "#TRNASEQUENCES:$tRNAseqSpaceFile MD5SUM:$tRNAseqsMD5\n";
print $tffh "#OTHERANNOTATIONS:$OtherAnnoFile MD5SUM:$OtherAnnoMD5\n";
print $tffh $tRNAfrags;
# print $tffh join "\n", @tRNAfrags;
# print $tffh "\n";
close $tffh or die $!;
print STDERR "[".getTime()."] Saved Marked \n";


print STDERR "[".getTime()."] COMPLETE\n";

#$tRNAannoFile !! - need to account for substring

# my @uniques = 
# grep {
# 	chomp $_;
# 	my @cols = split("\t", $_);
# 	$cols[6] == 1; # only fragment (fully -f 1) overlapped by a target region exactly once
# } @bedFrags;

# map {
# 	print STDOUT join "\t", $_; 
# 	print STDOUT "\n"; 
# } @uniques;



# getAnnotations seq mod behaviour issues?


#############################################################################
### Usage
#############################################################################

__END__

=head1 SYNOPSIS

pretRNAfragGen.pl -i hg19-tRNAs.bed -s hg19-tRNAs.fa -c ~/genomes/hg19/

=head2 DESCRIPTION

...

=head1 ARGUMENTS

=over

=item C<-h | --help> - Displays this help.

=item C<-i> - tRNA data - .bed

=item C<-s> - tRNA seq .fa (same names and bed)

=item C<-c> - chromosome directory

=back

=cut
