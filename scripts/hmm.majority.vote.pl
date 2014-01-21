#!/usr/bin/env perl
###############################################################################
#
#    hmm.majority.vote.pl
#   
#    Version 2	 
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

my $inputfile;
my $outputfile;
my $noNA;

$inputfile = &overrideDefault("inputfile.txt",'inputfile');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
$noNA = &overrideDefault(1,'noNA');

my %phylum;
my %class;
my %phylum_sensus;
my %class_sensus;


######################################################################
# CODE HERE
######################################################################


open(IN, $inputfile) or die("Cannot read file: $inputfile\n");
open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

while ( my $line = <IN> ) {
	chomp $line;   	
	my @splitline = split(/\t/,$line);	
	my @tax = split(/;/,$splitline[1]);                        #1
	my @scaffold = split(/_/,$splitline[0]);                   #2
######################### Get phylum assignments ##################################
	if (!exists($phylum{$scaffold[0]})){
		if (exists($tax[4])){
			$phylum{$scaffold[0]} = $tax[3]; 
		}
		else{
			$phylum{$scaffold[0]} = "NA";
		}
	}
	else{
		if (exists($tax[4])){
			$phylum{$scaffold[0]} = $phylum{$scaffold[0]}.";".$tax[3];
		}
		else{
			$phylum{$scaffold[0]} = $phylum{$scaffold[0]}.";NA";
		}
	}

###################### Get class assignments #####################################
	if (!exists($class{$scaffold[0]})){
		if (exists($tax[5])){
			$class{$scaffold[0]} = "$tax[4]"; 
		}
		else{
			$class{$scaffold[0]} = "NA";
		}
	}
	else{
		if (exists($tax[5])){
			$class{$scaffold[0]} = "$class{$scaffold[0]};$tax[4]";
		}
		else{
			$class{$scaffold[0]} = "$class{$scaffold[0]};NA";
		}
	}	
	

}
print OUT "scaffold\tphylum\tcount\tall.phylum.assignments\tclass\tall.class.assignments\n";

###Phylum consensus assignments

foreach my $scaffold (keys %phylum){
	my @phylum_tax = split(/;/, $phylum{$scaffold});
	my %phylum_taxcons;
	foreach my $assignments (@phylum_tax){
		$phylum_taxcons{$assignments}++;
	}
	if ($noNA == 1){
		$phylum_taxcons{"NA"} = 0;                             
	}
	my $count = 0;
	my $temptax;
	my $tempcount;	
	my $printed = 0;
	
	foreach my $tax_name (sort { $phylum_taxcons{$b} <=> $phylum_taxcons{$a} } keys %phylum_taxcons){
		$count++;
		if ($count == 2){
			if ($phylum_taxcons{$tax_name} == $tempcount){	
				if ($noNA == 0){
					$phylum_sensus{$scaffold} = "NA";
				}
			}
			else{
				$phylum_sensus{$scaffold} = $temptax;
			}
		$printed = 1;
		}
		$temptax = "NA";
		if ($tax_name ne ""){
			$temptax = $tax_name;
		}
		$tempcount = $phylum_taxcons{$tax_name};
	}

	if ($printed == 0){
		if ($noNA == 0){	
			$phylum_sensus{$scaffold} = "NA";
		}
		else{
			if ($temptax ne "NA"){
				$phylum_sensus{$scaffold} = $temptax;
			}
		}
	}
}

### Class consensus assignements

foreach my $scaffold (keys %class){
	my @class_tax = split(/;/, $class{$scaffold});
	my %class_taxcons;
	foreach my $assignments (@class_tax){
		$class_taxcons{$assignments}++;
	}
	if ($noNA == 1){
		$class_taxcons{"NA"} = 0;                             
	}
	my $count = 0;
	my $temptax;
	my $tempcount;	
	my $printed = 0;
	
	foreach my $tax_name (sort { $class_taxcons{$b} <=> $class_taxcons{$a} } keys %class_taxcons){
		$count++;
		if ($count == 2){
			if ($class_taxcons{$tax_name} == $tempcount){	
				if ($noNA == 0){
					$class_sensus{$scaffold} = "NA";
				}
			}
			else{
				$class_sensus{$scaffold} = $temptax;
			}
		$printed = 1;
		}
		$temptax = "NA";
		if ($tax_name ne ""){
			$temptax = $tax_name;
		}
		$tempcount = $class_taxcons{$tax_name};
	}

	if ($printed == 0){
		if ($noNA == 0){	
			$class_sensus{$scaffold} = "NA";
		}
		else{
			if ($temptax ne "NA"){
				$class_sensus{$scaffold} = $temptax;
			}
		}
	}
}

foreach my $scaffold (keys %class){	
	if (exists($phylum_sensus{$scaffold})){
		my $class_out = "NA";
		my $class_all_out = "NA";
		if (exists($class_sensus{$scaffold})){
			$class_out = $class_sensus{$scaffold};
			$class_all_out = $class{$scaffold};
		}
	my $count = ($phylum{$scaffold} =~ tr/;//) + 1;
	print OUT "$scaffold\t$phylum_sensus{$scaffold}\t$count\t$phylum{$scaffold}\t$class_out\t$class_all_out\n";
	}	
}


close IN;
close OUT;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inputfile|i:s", "outputfile|o:s", "outlegend|l:s", "noNA|n:+", "level|l:s");
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

    hmm.majority.vote_v2.pl

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
 [-inputfile -i]      Input file 
 [-outputfile -o]     List 
 [-noNA -n]           Ignore ambigous assignments (flag, default yes).
 
=cut
