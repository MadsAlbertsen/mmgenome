#!/usr/bin/env perl
###############################################################################
#
#    extract.fastq.for.reassembly.pl
#
#	Extracts a subset of sequences from a fastq file using a SAM file and 
#    the reference sequences.
#    
#    Copyright (C) 2012 Mads Albertsen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $inref;
my $insam;
my $infastq;
my $out;
my $fastqsplit;
my $samsplit;

$inref = &overrideDefault("ref.fa",'inref');
$insam = &overrideDefault("insam.sam",'insam');
$infastq = &overrideDefault("infastq.fastq",'infastq');
$out = &overrideDefault("out.fastq",'out');
$fastqsplit = &overrideDefault(" ",'fastqsplit');
$samsplit = &overrideDefault("_",'samsplit');
 
my %ref_id;
my %read_id;
my $count = 0;
my $print = 0;
my $tcount = 0;
my $ccount = 0;
my $mcount = 0;

######################################################################
# CODE HERE
######################################################################


## Extract the ID of each reference sequence

open(inref_fh, $inref) or die("Cannot read file: $inref\n");

while ( my $line = <inref_fh> ) {
	chomp $line;   	
	if ($line =~ m/>/) {
		$line =~ s/>//g;
		$ref_id{$line} = 1;
		$ccount++;
	}		
}

print "\nLoaded $ccount reference sequences from $inref\n\n";

close inref_fh;

## Extract the ID of the reads mapping to the subset of reference sequences

open(insam_fh, "$insam") or die("Cannot read file: $insam\n");

print "Looking for reads mapping to $inref in $insam\n\n";
$ccount = 0;

while ( my $line = <insam_fh> ) {
	chomp $line;   	
	if ($line !~ m/^\@/) {
		$tcount++;
		$ccount++;
		my @splitline = split("\t",$line);	
		if (exists($ref_id{$splitline[2]})){
              	my @splitline2 = split($samsplit,$splitline[0]);	
			$read_id{$splitline2[0]} = 1;
			$mcount++;
		}
		if ($tcount == 5000000){
			$tcount = 0;
			print "Looked though $ccount reads and $mcount mapped to the references\n";			
		}
	}
}

close insam_fh;

print "\n$mcount of $ccount reads mapped to $inref\n\n";

## Extract the reads mapping to the reference sequences

open(infastq_fh, "$infastq") or die("Cannot read file: $infastq\n");
open(out_fh, ">$out") or die("Cannot create file: $out\n");

print "Extracting mapped reads from the fastq file: $infastq\n\n";
$tcount = 0;
$ccount = 0;
$mcount = 0;

while ( my $line = <infastq_fh> ) {
	chomp $line; 
     $count++; 	
     if ($count == 1){
		$ccount++;
		$tcount++;
		my @splitline = split($fastqsplit, $line);	
		my $extract = $splitline[0];
		$extract =~ s/\@//;
		if (exists($read_id{$extract})){
			$print = 1;
			$mcount++;
		}
		if ($print == 1){ print out_fh "$line\n"; }	
	}
     if ($count == 2){
		if ($print == 1){ print out_fh "$line\n"; }	
	}
     if ($count == 3){
		if ($print == 1){ print out_fh "$line\n"; }
	}
     if ($count == 4){
		if ($print == 1){ print out_fh "$line\n"; }
		$print = 0;
		$count = 0;
	}
	if ($tcount == 5000000){
		$tcount = 0;
		print "Looked though $ccount reads and extracted $mcount\n";			
	}
}

close infastq_fh;
close out_fh;

print "\n$mcount of $ccount reads was extracted from $infastq\n\n";

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "insam|s:s", "inref|r:s","infastq|f:s","samsplit|x:s", "fastqsplit|y:s", "out|o:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

   extract.fastq.for.reassembly.pl

=head1 COPYRIGHT

   copyright (C) 2014 Mads Albertsen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

=head1 SYNOPSIS

extract.fastq.for.reassembly.pl  -s -r -f [-h -x -y -o]

 [-help -h]           Displays this basic usage information
 [-insam -s]          Mapping file in sam format.
 [-inref -r]          Reference sequences in fasta format.
 [-infastq -f]        Reads in fastq format.
 [-samsplit -x]       Delimiter used to distinguish read 1 and 2 in SAM files (default: "_").
 [-fastqsplit -y]     Delimiter used to distinguish read 1 and 2 in fastq files (default: " ").
 [-out -o]            Extracted reads in fastq format.
 
=cut
