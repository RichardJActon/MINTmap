package MINTmap::subs;
use warnings;
use strict;
use Digest::MD5;
use Sys::Hostname;

use Exporter qw(import);

#our @EXPORT_OK = qw(transcribe);

	# loadLookupTable
	# loadtRNAfasta
	# loadOtherAnnotations
	# getAnnotations
	# generatesPlates
	# createOutput
#our @EXPORT_OK = qw(
our @EXPORT = qw(
	parsetRNAanno
	getAllSubstrings
	readChr
	transcribe
	getCoords
	matchesToBed
	writetRNAseqSpace
);

# load pre-computed lookup table into memory
# sub loadLookupTable
# {
# 	my $hash_exclusive = $_[0];
# 	my $hash_notexclusive = $_[1];

# 	open my $ifh, "<$opt{l}" or die $!;

# 	# get md5sum of tRNA sequences from first line in lookup table.  If it doesn't match, report an error and exit
# 	my $md5check = <$ifh>; 
# 	chomp $md5check;
# 	$md5check =~ s/^.*MD5SUM://;
# 	if ($md5check ne $md5sum_trnasequences)
# 	{
# 		printf (STDERR "Error, exiting: md5sum defined in lookup table %s (%s) does not match that of the tRNA sequences %s (%s)\n", $opt{l}, $md5check, $opt{s}, $md5sum_trnasequences);
# 		exit (1);
# 	}

# 	$md5check = <$ifh>; 
# 	chomp $md5check;
# 	if (substr ($md5check, 0, 17) ne "#OTHERANNOTATIONS")
# 	{
# 		printf (STDERR "Error, the second line of the LOOKUP table must start with '#OTHERANNOTATIONS' even if the otherannotations file does not exist.  Exiting.\n");
# 		exit (1);
# 	}
# 	if ($using_otherannotations == 1) # only do the MD5CHECK if the file exists because it's optional
# 	{
# 		$md5check =~ s/^.*MD5SUM://;
# 		if ($md5check ne $md5sum_otherannotations)
# 		{
# 			printf (STDERR "Error, exiting: md5sum defined in lookup table %s (%s) does not match that of the tRF types %s (%s)\n", $opt{l}, $md5check, $opt{o}, $md5sum_otherannotations);
# 			exit (1);
# 		}
# 	}

# 	while (my $line = <$ifh>)
# 	{
# 		chomp $line; 
# 		my @splitline = split (/\t/, $line);
# 		my $fragseq = $splitline[0];

# 		# Note: a Y in col2 means the tRF sequence is exclusive to tRNA space and a N means it's not exclusive
# 		if ($splitline[1] eq "Y") 
# 		{
# 			$hash_exclusive->{$splitline[0]} = 0;
# 		}
# 		elsif ($splitline[1] eq "N")
# 		{
# 			$hash_notexclusive->{$splitline[0]} = 0;
# 		}
# 		else
# 		{
# 			printf (STDERR "Error, exiting: Lookup value of %s in %s not recognized\n", $splitline[1], $opt{l});
# 			exit (1);
# 		}
# 	}
# 	close ($ifh);
# }

# # load tRNA fasta into memory
# sub loadtRNAfasta
# {
# 	my $hash_tRNAseqs = $_[0];

# 	open my $ifh, "<$opt{s}" or die $!;
# 	$md5sum_trnasequences  = Digest::MD5->new->addfile($ifh)->hexdigest; # get md5sum of the tRNA sequence files and then open it
# 	seek $ifh, 0, 0;  # rewind the file after getting md5sum
# 	while (my $line = <$ifh>)
# 	{
# 		# read header from fasta file
# 		chomp $line; 
# 		if (!($line =~ m/^>/))
# 		{
# 			printf (STDERR "Error, exiting: Expected FASTA header but received %s\n", $line);
# 			exit (1);
# 		}
# 		my $trnaname = $line;
# 		$trnaname =~ s/^>//; # remove the > from the header name to get the tRNA name
# 		if (!($trnaname =~ m/^[-.|+_A-Za-z0-9]+$/))
# 		{
# 			printf (STDERR "Error, exiting: Only characters [-.|+_A-Za-z0-9] allowed in FASTA label but received %s.  Note, no whitespace is allowed\n", $trnaname);
# 			exit (1);
# 		}

# 		# read sequence from fasta file
# 		my $seq = <$ifh>;
# 		chomp $seq;
# 		if (!($seq =~ m/^[ATCGN]+$/))
# 		{
# 			printf (STDERR "Error, exiting: Only characters [ATCGN] allowed in FASTA sequence but received %s\n", $seq);
# 			exit (1);
# 		}
		
# 		# add to hash
# 		$hash_tRNAseqs->{$trnaname} = $seq;
# 	}
# 	close ($ifh);
# }

# # load other tRF annotations
# sub loadOtherAnnotations
# {
# 	my $hash_otherannotations = $_[0];

# 	open my $ifh, "<$opt{o}" or die $!;
# 	$md5sum_otherannotations = Digest::MD5->new->addfile($ifh)->hexdigest; # get md5sum of the tRNA sequence files and then open it
# 	seek $ifh, 0, 0;  # rewind the file after getting md5sum
# 	while (my $line = <$ifh>)
# 	{
# 		chomp $line;
# 		my @linesplit = split (/\t/, $line);
# 		# add to hash
# 		$hash_otherannotations->{$linesplit[0]} = $linesplit[1];
# 	}
# 	close ($ifh);
# }

# sub getAnnotations 
# {
# 	my $trfseq					= $_[0];
# 	my $hash_tRNAs				= $_[1];
# 	my $array_annotationoutput	= $_[2]; # store output here

# 	my $trfseq_len = length ($trfseq);

# 	# look for the tRFseq across all tRNA spliced sequences
# 	foreach my $trnaname (keys %{$hash_tRNAs})
# 	{
# 		my $trnaseq = $hash_tRNAs->{$trnaname};
# 		if (!defined $opt{z})
# 		{
# 			$trnaseq = $hash_tRNAs->{$trnaname} . "CCA"; # get tRNA spliced sequence and add CCA to it

# 			# check for -1 5' tRF since it's an exception
# 			my $trf5ptrunc = $trfseq;
# 			$trf5ptrunc =~ s/^.//; # remove the first NT from the tRF to check for 5' -1 post-modification possibilities
# 			if ($trf5ptrunc eq substr ($trnaseq, 0, length ($trf5ptrunc)))
# 			{
# 				push (@{$array_annotationoutput}, sprintf ("%s@%d%s.%d.%d", $trnaname, -1, substr ($trfseq, 0, 1), length ($trf5ptrunc), $trfseq_len));
# 			}
# 		}

# 		# now lets check the rest of them
# 		for (my $i = 0; $i < length ($trnaseq) - length ($trfseq) + 1; $i++)
# 		{
# 			if ($trfseq eq substr ($trnaseq, $i, $trfseq_len))
# 			{
# 				push (@{$array_annotationoutput}, sprintf ("%s@%d.%d.%d", $trnaname, $i + 1, $i + $trfseq_len, $trfseq_len));
# 			}
# 		}
# 	}

# 	if (scalar (@{$array_annotationoutput}) < 1)
# 	{
# 		printf (STDERR "Error, exiting: no annotations in sequence file %s found for tRF %s\n", $opt{s}, $trfseq);
# 		exit (1);
# 	}
# }

# # encode tRF sequences and store in output hash called MINTplates
# # e.g. java -cp ./MINTplates/ MINTcodes_tRF_license_plates -e MINTplates/example_sequences_to_encode.txt
# sub generatesPlates 
# {
# 	my $allexpressed_tRFs = $_[0];
# 	my $outputplates_hash = $_[1];

# 	my $host = hostname;
# 	my $hostprefix="$$\_$host"; # set prefix to unique host and pid incase process is ran multiple times on same machine
# 	open my $ofh, ">tmp.mintmap.seqstoencode.$hostprefix.txt" or die $!;
# 	foreach my $mykey (keys %{$allexpressed_tRFs}) 
# 	{
# 		chomp $mykey;
# 		printf ($ofh "%s\n", $mykey);
# 	}
# 	close ($ofh);

# 	# run java program
# 	printf ("Running Java-based MINTplates\n");
# 	`java -cp $opt{j} MINTcodes_tRF_license_plates -e tmp.mintmap.seqstoencode.$hostprefix.txt > tmp.mintmap.seqstoencode.encoded.$hostprefix.txt`;

# 	open my $ifh, "<tmp.mintmap.seqstoencode.encoded.$hostprefix.txt" or die $!;
# 	foreach my $line (<$ifh>)
# 	{
# 		if ($line =~ m/\t/) # only look at MINTplate output lines with tabs in them as the others are information
# 		{
# 			chomp $line;
# 			my @splitline = split (/\t/, $line);

# 			if ($opt{t} ne "tRF") # if fragment is not a tRF, rename the MINTplate prefix
# 			{
# 				$splitline[1] =~ s/^tRF/$opt{t}/;
# 			}
# 			$outputplates_hash->{$splitline[0]} = $splitline[1];
# 		}
# 	}
# 	close ($ifh);

# 	# cleanup
# 	unlink ("tmp.mintmap.seqstoencode.$hostprefix.txt");
# 	unlink ("tmp.mintmap.seqstoencode.encoded.$hostprefix.txt");
# }

# sub createOutput
# {
# 	my $read_hash = $_[0];
# 	my $annotation_hash = $_[1];
# 	my $otherannotations_hash = $_[2]; # currently used for tRF-type
# 	my $mintplates_hash = $_[3];
# 	my $total_frags_in_file = $_[4];
# 	my $tRFtypes = $_[5];
# 	my $filename = "";

# 	# output the tRFs in decreasing order of expression (double counting is not a problem because we are dealing with the raw reads)
# 	$filename = "$opt{p}-$scriptversion-$tRFtypes.expression.txt";
# 	printf ("Creating output file: %s\n", $filename);
# 	open my $ofh, ">$filename" or die $!;
# 	$filename = "$opt{p}-$scriptversion-$tRFtypes.expression.html";
# 	printf ("Creating output file: %s\n", $filename);
# 	open my $ofh_html, ">$filename" or die $!;
# 	printf ($ofh "MINTbase Unique ID\ttRF sequence\ttRF type(s)\tUnnormalized read counts\tRPM read counts (using all counts from this file[%d] in denominator)\tRPM read counts (using read count of input file from -f parameter[%d] in denominator)\tRPM read counts (using read count from optional -d parameter[%s] in denominator)\tSequence locations in tRNA space (comma deliminated)\n", $total_frags_in_file, $stat_totalstartingreads, defined ($opt{d}) ? $opt{d} : "na" );
# 	printf ($ofh_html "<html><head><title>%s expression</title></head>", $tRFtypes);
# 	printf ($ofh_html '<style> tr:nth-of-type(odd) { background-color:#ccc; } body { font-size: 18px; } table { font-size: 16px; } </style>');
	
# 	printf ($ofh_html '<body><p style="font-size:22px; display:inline">Table of %s for %s</p>', $tRFtypes, $opt{a});
# 	printf ($ofh_html '<br />Created by the <a target="_blank" href="http://cm.jefferson.edu">Computational Medicine Center</a> at <a target="_blank" href="http://www.jefferson.edu/">Thomas Jefferson University</a> using the MINTmap tool located <a target="_blank" href="http://cm.jefferson.edu/MINTcodes/">here</a>.<br />Please cite: Loher, P. <i>et al.</i> MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. <i>Sci. Rep.</i> 7, 41184; doi: 10.1038/srep41184 (2017).');
# 	printf ($ofh_html '<br /><br /><table style="width=100%%"><tr><td><b>MINTbase Unique ID</b><br />(sequence derived)</td><td><b>tRF sequence</b></td><td><b>tRF type(s)</b></td><td><center><b>Unnormalized<br />read counts</b></center></td><td><center><b>RPM read counts (using all counts from this file[%d] in denominator)</b></center></td><td><center><b>RPM read counts (using read count of input file from -f parameter[%d] in denominator)</b></center></td><td><center><b>RPM read counts (using read count from optional -d parameter[%s] in denominator)</b></center></td><td><center><b>MINTbase Summary Record</b></center></td><td><center><b>Sequence locations in tRNA space<br />(comma deliminated)</b></td></center></tr>', $total_frags_in_file, $stat_totalstartingreads, defined ($opt{d}) ? $opt{d} : "na");

# 	foreach my $mykey (sort { $read_hash->{$b} <=> $read_hash->{$a} or $a cmp $b } keys %{$read_hash}) 
# 	{
# 		if ($read_hash->{$mykey} != 0) # don't print things with 0 expression
# 		{
# 			my $unnorm_numreads = $read_hash->{$mykey};
# 			my $rpm_1 = (($unnorm_numreads/$total_frags_in_file)*1000000);
# 			my $rpm_2 = (($unnorm_numreads/$stat_totalstartingreads)*1000000);
# 			my $rpm_3 = "na";
# 			my $annotations = join (', ', @{$annotation_hash->{$mykey}});
# 			if (defined $opt{d})
# 			{
# 				$rpm_3 = sprintf ("%.2f", (($unnorm_numreads/$opt{d})*1000000));
# 			} 
		
# 			my $otherannotations = 'na';
# 			if ($using_otherannotations == 1)
# 			{
# 				$otherannotations = $otherannotations_hash->{$mykey};
# 			}

# 			printf ($ofh "%s\t%s\t%s\t%d\t%.2f\t%.2f\t%s\t%s\n", $mintplates_hash->{$mykey}, $mykey, $otherannotations, $unnorm_numreads, $rpm_1, $rpm_2, $rpm_3, $annotations);
# 			my $mintbase_html = "na";
# 			if ($opt{t} eq "tRF")
# 			{
# 				$mintbase_html = sprintf ('<a target="_blank" href="https://cm.jefferson.edu/MINTbase/InputController?g=%s&v=s&fs=%s">Summary</a>', $opt{a}, $mykey);
# 			}
# 			printf ($ofh_html '<tr><td>%s</td><td>%s</td><td>%s</td><td>%d</td><td>%.2f</td><td>%.2f</td><td>%s</td><td>%s</td><td>%s</td></tr>', $mintplates_hash->{$mykey}, $mykey, $otherannotations, $unnorm_numreads, $rpm_1, $rpm_2, $rpm_3, $mintbase_html, $annotations);
# 		}
# 	}
# 	printf ($ofh_html "</table></body></html>");
# 	close ($ofh);
# 	close ($ofh_html);

# 	$filename = "$opt{p}-$scriptversion-$tRFtypes.countsmeta.txt";
# 	printf ("Creating output file: %s\n", $filename);
# 	open $ofh, ">$filename" or die $!;
# 	printf ($ofh "Total reads in -f input file\tTotal unnormalized reads in %s\tPercent\n", $tRFtypes);
# 	printf ($ofh "%ld\t%ld\t%.2f%%\n", $stat_totalstartingreads, $total_frags_in_file, ($total_frags_in_file / $stat_totalstartingreads) * 100);
# 	close ($ofh);
# }

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

1;