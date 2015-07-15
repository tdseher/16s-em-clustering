#!/usr/bin/perl
use strict;
use warnings;

my $count = 0;
while(<>)
{
    chomp;
    my $line = $_;
    $count += 1;
    
    if (($count % 4) == 1)
    {
        #if ($line =~ m/^\@/)
        print ">", substr($line, 1), "\n";
    }
    elsif (($count % 4) == 2)
    {
        print "$line\n";
    }
    
    #1 @header1
    #2 ACTGTGT
    #3 +header2
    #4 fhhacca
}
