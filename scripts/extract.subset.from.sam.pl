#!/usr/bin/env perl
###############################################################################
#
#    extract.subset.from.sam.using.list.pl
#
#	 Extracts a subset of sequences from a SAM file using a list of references
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
my $outputfile;

$inref = &overrideDefault("inref.txt",'inref');
$insam = &overrideDefault("insam.txt",'insam');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
 
my %extract;
my $ccount = 0;

######################################################################
# CODE HERE
######################################################################


## Extract the ID of each reference sequence

open(inref_fh, $inref) or die("Cannot read file: $inref\n");

while ( my $line = <inref_fh> ) {
	chomp $line;   	
	if ($line =~ m/>/) {
		$line =~ s/>//g;
		$extract{$line} = 1;
		$ccount++;
	}		
}

print "\nLoaded $ccount reference sequences from $inref\n\n";


close inref_fh;


open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");
open(INsam, "$insam") or die("Cannot read file: $insam\n");

my $count = 0;

while ( my $line = <INsam> ) {
	chomp $line;   	
	$count++;
	if ($count == 1){
		print OUT "$line\n";
	}
	else{
		if ($line =~ m/\@SQ/) { 	
			my @splitline = split("\t",$line);
			my @splitline1 = split(":",$splitline[1]);
			if (exists($extract{$splitline1[1]})){
				print OUT "$line\n";
			}
		}
		else{
			my @splitline = split("\t",$line);
			if (exists($extract{$splitline[2]})){
				print OUT "$line\n";
			}
		}
	}
}

close INsam;
close OUT;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "insam|s:s", "inref|r:s", "outputfile|o:s");
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

    vprobes.generateprobes.pl

=head1 COPYRIGHT

   copyright (C) 2012 Mads Albertsen

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

script.pl  -i [-h]

 [-help -h]           Displays this basic usage information
 [-insam -s]          Input sam file.
 [-inref -r]          Subset to be extracted in fasta format.
 [-outputfile -o]     Outputfile.
 
=cut
