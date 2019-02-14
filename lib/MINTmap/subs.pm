package MINTmap::subs;
use warnings;
use strict;
use Digest::MD5;
use Sys::Hostname;
use POSIX qw(strftime);

use Exporter qw(import);

#our @EXPORT_OK = qw(transcribe);

	# loadLookupTable
	# loadtRNAfasta
	# loadOtherAnnotations
	# getAnnotations
	# generatesPlates
	# createOutput

	# tRNAsubStrTrans
	# tRNAsubStrAnno
#our @EXPORT_OK = qw(
our @EXPORT = qw(
	parsetRNAanno
	getAllSubstrings
	readChr
	transcribe
	getCoords
	matchesToBed
	writetRNAseqSpace
	getTime
	tRNAsubStrSave
);

#############################################################################
### xxx
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

		$tRNA{score} = $cols[4];
		$tRNA{strand} = $cols[5];
		$tRNA{tStart} = $cols[6];
		$tRNA{tEnd} = $cols[7];
		$tRNA{na} = $cols[8];
		$tRNA{nExon} = $cols[9];
		$tRNA{exonInfo1} = $cols[10];
		$tRNA{exonInfo2} = $cols[11];

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

sub tRNAsubStrTrans {
	my %tRNAs = %{$_[0]};
	my $substr1 = $_[1];
	my $substr2 = $_[2];
	if ($substr2 =~ /^-/) {
		$substr2 =~ s/^-//;
		$tRNAs{fStart} = $tRNAs{fStart} + ($substr1 + 0); 
		$tRNAs{fEnd} = $tRNAs{fEnd} - ($substr2 + 0);
	} else {
		$tRNAs{fStart} = $tRNAs{fStart} + ($substr1 + 0);
		$tRNAs{fEnd} = $tRNAs{fStart} + ($substr2 + 0);
	}
	return \%tRNAs;
}

sub tRNAsubStrAnno {
	my %tRNAs = %{$_[0]};
	my $substr1 = $_[1];
	my $substr2 = $_[2];
	foreach my $tRNAname (keys %tRNAs) {
		$tRNAs{$tRNAname} = tRNAsubStrTrans($tRNAs{$tRNAname},$substr1,$substr2);
	}
	return \%tRNAs;
}

sub tRNAsubStrSave {
	my %tRNAs = %{$_[0]};
	my $substr1 = $_[1];
	my $substr2 = $_[2];
	my $fileName = $_[3];
	my %proctRNAs = %{tRNAsubStrAnno(\%tRNAs, $substr1, $substr2)};
	open my $tRNAsubStrAnnoFH, ">$fileName" or die $!;
	foreach my $tRNAname (sort keys %proctRNAs) {
		print $tRNAsubStrAnnoFH join("\t",
			$proctRNAs{$tRNAname}{chr},
			$proctRNAs{$tRNAname}{fStart},
			$proctRNAs{$tRNAname}{fEnd},
			$tRNAname,
			$proctRNAs{$tRNAname}{score},
			$proctRNAs{$tRNAname}{strand},
			$proctRNAs{$tRNAname}{tStart},
			$proctRNAs{$tRNAname}{tEnd},
			$proctRNAs{$tRNAname}{na},
			$proctRNAs{$tRNAname}{nExon},
			$proctRNAs{$tRNAname}{exonInfo1},
			$proctRNAs{$tRNAname}{exonInfo2}
		);
		print $tRNAsubStrAnnoFH "\n";
	}
	close $tRNAsubStrAnnoFH or die $!;
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

sub getTime {
    return strftime "%H:%M:%S", localtime;
}

1;