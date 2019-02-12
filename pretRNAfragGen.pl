#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Parallel::Loops;

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

# Argument getting
#######################################

GetOptions(
	"tRNAannoFile|i=s"		=> \$tRNAannoFile,
	"tRNAseqFile|s=s"		=> \$tRNAseqFile,
	"chrDir|c=s"			=> \$chrDir,
	"maxProcs|n=s"			=> \$maxProcs,
	"help|?" 				=> \$help
	) or pod2usage(1);

my $pl = Parallel::Loops -> new($maxProcs);

pod2usage(1) if $help;

#############################################################################
### Fxn
#############################################################################

sub parsetRNAanno {
	my @lines = @{$_[0]};
	
	my %tRNAs;

	for (my $i = 0; $i < scalar @lines; $i++) {
		my @cols = split "\t", $lines[$i];
		my %tRNA;
		$tRNA{chr} = $cols[0];
		$tRNA{fStart} = $cols[1];
		$tRNA{fEnd} = $cols[2];
		
		$cols[3] =~ m/
			(?<nmt>\w*)-?tRNA-(?<aa>i?\w{3})
			(?:\w+)?-(?<anticodon>\w+)-\w+-\d+
			/x;
		$tRNA{nmt} = $+{nmt};
		$tRNA{aa} = $+{aa};
		$tRNA{anticodon} = $+{anticodon};

		$tRNA{strand} = $cols[5];
		$tRNA{tStart} = $cols[6];
		$tRNA{tEnd} = $cols[7];

		$tRNA{nExon} = $cols[9];

		#$tRNA{exonInfo1} = $cols[9];
		#$tRNA{exonInfo2} = $cols[10];
		$cols[10] =~ /(?<ep1>\d+),(?<ep2>\d*),?/;
		$tRNA{ep1} = $+{ep1};
		$tRNA{ep2} = $+{ep2};
		$cols[11] =~ /(?<ep3>\d+),(?<ep4>\d*),?/;
		$tRNA{ep3} = $+{ep3};
		$tRNA{ep4} = $+{ep4};

		# if not defined ep2/ep4
		if ($tRNA{ep2} ne "" || $tRNA{ep4} ne "") {
			$tRNA{is} = $tRNA{ep2} + 1;
			$tRNA{ie} = ($tRNA{ep4} - $tRNA{ep1}) + $tRNA{ep2};
		} else {
			$tRNA{is} = "NA";
			$tRNA{ie} = "NA";
		}
		
		#$tRNA{tRNAname} = $cols[2]; 
		$tRNAs{$cols[3]} = \%tRNA;
	}
	return %tRNAs;
}

sub getAllSubstrings {
	#@Ns = 16..50;
	# of length a though b
	my @Ns = $_[0]..$_[1];
	my $string = $_[2];
	my @substrings = map {
		$string =~ m/(?= (.{${_}}) )/gx;
		#$string => '';
	} @Ns;
	return \@substrings;
}

sub readChr {
	my $chr = $_[0];
	my $chrDir = $_[1];
	my $chrFile = "$chrDir\/chr$chr.fa.gz"; ####!!!! 
	my $str;
	open (my $fh, "zcat $chrFile |") or die "unable to open $chrFile";
	while (my $line = <$fh>) {
		chomp $line;
		if ($line =~ /^>/) {
			next;
		} else {
			$str .= uc $line;
		}
	}
	close $fh or die "unable to close $chrFile";
	return $str;
}

sub transcribe {
	my %lookup = (
		A => "T",
		C => "G",
		G => "C",
		T => "A",
		N => "N",
	);
	join "", map { $lookup{$_} } split "", $_[0];
}

sub getCoords {
	my ($search,$dict) = @_;
	# chr?
	my @matches;
	while ($dict =~ /(?= ${search} )/xg) {
		# $+[0] = $-[0] - lookahead
		# print "$dict\n$search\n";
		push @matches, [ $-[0], $-[0] + (length($search)-1) ]; 
	}
	return \@matches;
}

sub matchesToBed {
	my @matches;
	my %matches = %{$_[0]};
	foreach my $substr (sort keys %matches) {
		foreach my $chr (sort keys %{$matches{$substr}}) {
			foreach my $strand (sort keys %{$matches{$substr}{$chr}}) {
				my $res = $matches{$substr}{$chr}{$strand};
				if ($res ne 'NA') {
					my @coords = @{$res};
					for (my $i = 0; $i < scalar @coords; $i++) {
						my %strandLU = (fwd => "+", rev => "-");
						push @matches, [
							"chr$chr",
							$coords[$i][0],
							$coords[$i][1],
							$substr,
							"0",
							$strandLU{$strand}
						];
					}
				}
			}
		}
	}
	return \@matches;
}

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
my @tRNAannoLines;
open TRNAANNO, "$tRNAannoFile" or die "unable to open $tRNAannoFile";
while (my $line = <TRNAANNO>) {
	chomp $line;
	push @tRNAannoLines, $line;
}
close TRNAANNO or die "unable to close $tRNAannoFile";

# 100bp flank fasta
# hg19-tRNA-100bp53.fa

my %tRNAseqs;
my $tRNAname;
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

my @substrings;
foreach my $k (keys %tRNAs) {
	my $str;
	if (defined $tRNAs{$k}{seq}) {
		$str = substr $tRNAs{$k}{seq}, 85, -85; # introns...
	} else { 
		next 
	}
	# print "$k = $str\n";
	push @substrings, @{getAllSubstrings(16,50,$str)};
}
#85,-85
my %substr = map {$_ => ''} @substrings;

# print Dumper @substrings;
# print Dumper %substr;

# Search Genome for unique substrings
#######################################

# my @chr = (1..22,"X","Y","M");
my @chr = ("T"); ##!!

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

# print Dumper %matches;
# my @res = @{matchesToBed(\%matches)};
# for (my $i = 0; $i < scalar @res; $i++) {
# 	print STDOUT join "\t", @{$res[$i]};
# print STDOUT "\n";
# }
map {
	print STDOUT join "\t", @{$_};
	print STDOUT "\n";
} @{matchesToBed(\%matches)};

#print Dumper @{matchesToBed(\%matches)};
# ~~?? intersect with tRNA gene regions here or with bedtools
# 

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
