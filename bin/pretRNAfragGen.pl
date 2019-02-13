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
print STDERR "Reading tRNA annotation...\n";
my @tRNAannoLines;
open TRNAANNO, "$tRNAannoFile" or die "unable to open $tRNAannoFile";
while (my $line = <TRNAANNO>) {
	chomp $line;
	push @tRNAannoLines, $line;
}
close TRNAANNO or die "unable to close $tRNAannoFile";
print STDERR "Read tRNA annotation\n";
# 100bp flank fasta
# hg19-tRNA-100bp53.fa

my %tRNAseqs;
my $tRNAname;
print STDERR "Reading tRNA Sequences...\n";
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
print STDERR "Read tRNA Sequences\n";

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
print STDERR "Generating tRNA substrings...\n";
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
#85,-85
my %substr = map {$_ => ''} @substrings;
print STDERR "Generated tRNA substrings\n";

print STDERR "Saving tRNA substrings...\n";
my $tRNAseqsMD5 = Digest::MD5::md5_hex($tRNAseqSpaceFasta);
open(my $tsfh, ">$tRNAseqSpaceFile") or die $!;
#my $tRNAseqsMD5 = Digest::MD5->new->addfile($tsfh)->hexdigest;
print $tsfh $tRNAseqSpaceFasta;
close $tsfh or die $!;
print STDERR "Saved tRNA substrings\n";
# print Dumper @substrings;
# print Dumper %substr;

# Search Genome for unique substrings
#######################################

# my @chr = (1..22,"X","Y","M");
my @chr = ("T"); ##!!

print STDERR "Matching \n";
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

print STDERR "Matched \n";
# print Dumper %matches;

# map {
# 	print STDOUT join "\t", @{$_};
# 	print STDOUT "\n";
# } @{matchesToBed(\%matches)};
print STDERR "Formatting Matched \n";
my $matchesString;
map {
	$matchesString .= join "\t", @{$_};
	$matchesString .= "\n";
} @{matchesToBed(\%matches)};
print STDERR "Formatted Matched \n";
#print Dumper @{matchesToBed(\%matches)};
# ~~?? intersect with tRNA gene regions here or with bedtools
# 
# print STDOUT $matchesString;  

#qx/echo "$matchesString" | bedtools intersect -a stdin -b ..\/tmp-tRNAlookupTest.bed -f 1 -c/;#$tRNAannoFile

print STDERR "Intersecting \n";
my @bedFrags = readpipe join("",
	"echo \"$matchesString\" | ",
	"bedtools intersect ",
	"-a stdin ",
	"-b $tRNAannoFile ", #$tRNAannoFile #..\/tmp-tRNAlookupTest.bed
	"-f 1 -c"
);
print STDERR "Intersected \n";
# my @uniques = 
# grep {
# 	chomp $_;
# 	my @cols = split("\t", $_);
# 	$cols[6] == 1; # only fragment (fully -f 1) overlapped by a target region exactly once
# } @bedFrags;

## generate other annotations # MD%
# 5' 3' etc - defer
print STDERR "Annotating \n";
my @otherAnno;
$pl -> share(\@otherAnno);
$pl -> foreach(\@bedFrags,sub {
	my $line = $_;
	chomp $line;
	my @cols = split "\t", $line;
	if ($cols[6] == 1) {
	 	push @otherAnno, join "\t", $cols[3],"Unique";
	} elsif ($cols[6] == 0) {
		push @otherAnno, join "\t", $cols[3],"Non-unique";
	} else {
		push @otherAnno, join "\t", $cols[3],"tRNAspace_multi-mapping";
	}
});
print STDERR "Annotated \n";

print STDERR "Saving Annotated \n";
open(my $oafh, ">$OtherAnnoFile") or die $!;
my $OtherAnnoMD5 = Digest::MD5->new->addfile($oafh)->hexdigest;
print $oafh join "\n", @otherAnno;
print $oafh "\n";
close $oafh or die $!;
print STDERR "Saved Annotated \n";


print STDERR "Marking Unique \n";
my @tRNAfrags;
$pl -> share(\@tRNAfrags);
$pl -> foreach(\@bedFrags,sub {
	my $line = $_;
	chomp $line;
	my @cols = split "\t", $line;
	if ($cols[6] == 1) {
	 	push @tRNAfrags, join "\t", $cols[3],"Y";
	} else {
		push @tRNAfrags, join "\t", $cols[3],"N";
	}
});
print STDERR "Marked Unique \n";

print STDERR "Saving Marked \n";
open(my $tffh, ">$tRNAfragsFile") or die $!;
print $tffh "#TRNASEQUENCES:$tRNAseqSpaceFasta MD5SUM:$tRNAseqsMD5\n";
print $tffh "#OTHERANNOTATIONS:$OtherAnnoFile MD5SUM:$OtherAnnoMD5\n";
print $tffh join "\n", @tRNAfrags;
print $tffh "\n";
close $tffh or die $!;
print STDERR "Saved Marked \n";
## generate list - head with 
# '#TRNASEQUENCES:filename.fa MD5SUM:nnn'
# '#OTHERANNOTATIONS:OtherAnnotations.txt MD5SUM:nnn'

#$tRNAannoFile !! - need to account for substring

# map {
# 	print STDOUT join "\t", $_; 
# 	print STDOUT "\n"; 
# } @uniques;

# if >0 tRNA space exclusive?
# Y = 1, N = 0 or >1
# >1 in additional in other annotation?

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
