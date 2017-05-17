#! /usr/bin/env perl
#Written by Dianne Velasco
use strict;
use warnings;

# User designates input and output files
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
# Try to open the file or die trying
open (IN, "$infile") or die "error reading $infile";
# Try to write fasta output file or die trying
open (OUT, ">$outfile") or die "error creating samples.fasta";
# While loop to read the file, transliterate the fastq identification format to fasta format
my $line = "";
my $seqid = "";
my $seq = "";
my $qualid = "";
my $qual = "";
my $i = 0;
while ($line = <IN>) {
	# Remove new line character from lines
	chomp ($line);
		if ($i == 0) {
			$seqid = $line;
			$i++;
		}
		elsif ($i == 1) {
			$seq = $line;
			$i++;
		}
                elsif ($i == 2) {
                        $qualid = $line; 
                        $i++;
                }
                elsif ($i == 3) {
                        $qual = $line; 
			print OUT "$seqid\t$seq\t$qualid\t$qual\n"; # print file with one line per read
                        $i = 0;
                }
		else {print "did not work!\n";
		}
}
close IN;
close OUT;
