#!/usr/bin/perl
 
# look for exact peptide matches in a protein fasta file
# Cache which proteins have which 5-aa sequences in them
# Then for a given peptide, look at the list of protein containing the first
# five aa, and do a simple substring match
 
use strict;
use warnings;
BEGIN {eval{ require Smart::Comments; import Smart::Comments ('##')};}
 
# hoh uses 1.2GB on uniprothum. 1.5 for hoa, but only 60s
# hostrings use 800M in 60s
my $Seed_Length = 4; 
# String in outfile when you don’t find a protein
my $Not_Found_String = "NOT_FOUND";
 
use Getopt::Long;
my ($peps_file, $fasta_file, $out_file);
my $Usage = "
perl find_pep_in_protein.pl -pep peptide_list.txt -fasta human.protein.faa -out results.txt
This will read a one-column list of peptides from peptide_list.txt
    and a FASTA protein file, and send results to results.txt
The result file is two tab-delimited columns: peptide, protein ID.
If the peptide wasn’t found, this column will say $Not_Found_String.
";
GetOptions(
"out_file=s" => \$out_file,
"peps_file=s" => \$peps_file,
"fasta_file=s" => \$fasta_file,
) or die $Usage;
die "Not enough arguments:\n$Usage\n"
    unless defined($out_file) && defined($peps_file) && defined($fasta_file);
 
my $prots_ref = read_fasta($fasta_file);
my %seed_proteins; # which proteins contain each small e.g., 5-aa seed seq?
my $seq_count = 0;
foreach my $seqid (keys %$prots_ref) {
    my $seq = $prots_ref->{$seqid}{sequence};
    my $max_l = length ($seq) - $Seed_Length;
    for (my $i = 0; $i <= $max_l; $i++) {
my $seed = substr($seq, $i, $Seed_Length);
# Assume repeats are rare enough to ignore duplicates of a given seed
$seed_proteins{$seed}.=":$seqid";
#push @{$seed_proteins{$seed}}, $seqid++;
    }
    $seq_count++;
    if ($seq_count % 1000 == 0) {
warn "$seq_count done. ",scalar(keys %seed_proteins), " seeds\n";
    }
}
 
open my $peps_fh, "<", $peps_file
    or die "Can't open peps file $peps_file: $!\n";
my %found;
my @not_found;
open my $out_fh, ">", $out_file
    or die "Can't open out file $out_file: $!\n";
 
#print "Starting search at "", scalar(localtime), "\n";
while (<$peps_fh>) {
    s/\r?\n//;
    my $peptide = $_;
    my $seed1 = substr($peptide, 0, $Seed_Length);
    my $seed2 = substr($peptide, -$Seed_Length);
    if ($seed1 eq $seed2) {
warn "$seed1 at beginning and end of $peptide\n";
my @prot_ids = split /:/, $seed_proteins{$seed1};
shift @prot_ids; # empty first id
my %done;
foreach my $prot_id (@prot_ids) {
   my $prot_seq = $prots_ref->{$prot_id}{sequence};
   if (index($prot_seq, $peptide) != -1) {
$found{$peptide} = $prot_id;
print $out_fh "$peptide\t$prot_id\n";
last;
   }
}
    }
    else {
my @prot_ids1 = split /:/, $seed_proteins{$seed1};
shift @prot_ids1;
my @prot_ids2_tmp = split /:/, $seed_proteins{$seed2};
shift @prot_ids2_tmp;
my %prot_ids2 = map {($_=>1)} @prot_ids2_tmp;
 
foreach my $prot_id1 (@prot_ids1) {
   if (exists $prot_ids2{$prot_id1}) { # both seeds in same protein
my $prot_seq = $prots_ref->{$prot_id1}{sequence};
if (index($prot_seq, $peptide) != -1) {
   $found{$peptide} = $prot_id1;
   print $out_fh "$peptide\t$prot_id1\n";
   last;
}
   }
}
    }
 
    if (! $found{$peptide}) {
print $out_fh "$peptide\t$Not_Found_String\n";
push @not_found, $peptide;
    }
}
close $peps_fh;
close $out_fh;
warn scalar(@not_found), " peptides NOT found, ", scalar(keys %found), " found\n";
#print “Ending search at “, scalar(localtime), “\n”;
 
 
exit;
 
################################################################################
sub read_fasta {
    my ($fasta_in) = @_;
    open my $fh, '<', $fasta_in or die "FASTA in $fasta_in: $!\n";
    my %seqs = ();
    my $current_id;
    while (<$fh>) {
s/\r?\n//;  # chomp even Windows \r\n
if (/>(\S+)(\s+(.*))?/) {
   $current_id = $1;
   my $desc = defined $2 ? $3: "";
   $seqs{$current_id} = {
sequence  => "",
description  => $desc,
id  => $current_id,
   }
} else {
   $seqs{$current_id}{sequence} .= $_;
}
    }
    close $fh;
    return \%seqs;
} 
