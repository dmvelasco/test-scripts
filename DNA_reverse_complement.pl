#!/bin/perl

# Reverse compliment DNA sequence, including IUPAC ambiguity codes
# Usage: perl DNA_reverse_complement.pl < infile.fa > outfile.fa

open (FASTA, $file_name) or die "error $!";

# sequence, built by joining lines =~ /^[A-Z]+$/
my $sequence = '';

while (my $entry = <FASTA>)
{
    if ($entry =~ m/^[A-Z]+$/) {
        # Assemble the sequence from separate lines
        chomp($entry);
        $sequence .= $entry;
    }
    elsif ($entry =~ m/^\s*$/) { 
        # process and print the sequence and blank line, reset for next
        $sequence = reverse $sequence;
        $sequence =~ tr/ACGTRYSWKMBDHVNacgtryswkmbdhvn/TGCAYRSWMKVHDBNtgcayrswmkvhdbn/;
        print "$sequence\n";
        print "\n";
        $sequence = '';
    }
    else { # header
        print $entry;
    }
}

# Print the last sequence if the file didn't end with blank line    
if (length $sequence) {
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGUacgu/UGCAugca/;
    print "$sequence\n";
}
