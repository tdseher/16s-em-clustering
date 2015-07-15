#!/usr/bin/perl -w
use strict;

# Removes hyphen characters from FASTA file input
# prints to STDOUT

while(<>)
{
    chomp;
    my $line = $_;
    unless ($line =~ m/^>/)
    {
        $line =~ s/[ -.]//g;
    }
    print "$line\n"
}
    