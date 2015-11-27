#!/usr/bin/env perl
###############################################################################
#
#    network.pl version 1.0
#    
#    Indicates connections between scaffolds and produces cytoscape 
#    files for an easy visualization 
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

my $infile;
my $outcon;
my $enddist;
my $minlength;
my $minconnections;
my $avgreadlength;
my $headersplit;

$infile = &overrideDefault("infile.sam",'infile');
$outcon = &overrideDefault("network.txt",'outcon');
$enddist = &overrideDefault("500",'enddist');
$minlength = &overrideDefault("3000",'minlength');
$minconnections = &overrideDefault("10",'minconnections');
$avgreadlength = &overrideDefault("126",'avgreadlength');
$headersplit = &overrideDefault("_",'headersplit');

my $line;
my $old;
my $new;
my $count = 0;
my $start = 0;
my $connectionorder;
my $readcount = 0;
my $printreadcount = 0;
my $out1;
my $out2;
my @splitline;
my @splitline1;
my @splitline2;
my @contigname;
my @contiglength;
my %contigs;
my %reads;
my %connections;
my %contigswhit;
my %condensedhit;

######################################################################
# CODE HERE
######################################################################

open(OUTconC, ">$outcon") or die("Cannot create connections file: $outcon\n");   
open(IN, $infile) or die("Cannot open the infile $infile\n");


if (my $is_even = $avgreadlength % 2 == 0){                                                                               #To be able to get a better estimate of where the read maps in the start and end / otherwise it would just have used the start position
	$avgreadlength = $avgreadlength/2;
}
else{
	$avgreadlength = ($avgreadlength+1)/2;
}


while ( $line = <IN> ) {
	chomp $line;                                                         												     
	$count++;
	@splitline = split(/\t/, $line); 	
	if ($line =~ m/\@SQ/) {                               												                  #if we are in the contig header area then retrive all contigs/scaffolds and store the name and length in the hash: contig		
			@contigname = split(/:/, $splitline[1]);                     												  #Retrive the contig name
			@contiglength = split(/:/, $splitline[2]);																	  #Retrive the contig length
			$contigs{$contigname[1]} = $contiglength[1];                   												  #Make a hash with key = "contig name" and value = "contig length"			
		}	
	else {
		if ($line !~ m/(\@PG|\@HD|\@SQ|\@RG)/) { 
			@splitline = split(/\t/, $line); 
			$start = 0;
			if ($splitline[1] != 19 and $splitline[1] != 35){															  #SAM Flags that indicate that these PE reads are maping as they are supposed to.. hence they are not interesting.. There are other PE flags that could be included to speed up the script.
				if (($splitline[3]+$avgreadlength) <= $enddist or ($splitline[3]+$avgreadlength) >= ($contigs{$splitline[2]}-$enddist)) {            #The read is required to hit within a certain distance from the contigs ends.. The middle postition of the read is used
					@splitline1 = split(/$headersplit/, $splitline[0]);													  #The read header - assumes "name_1xyz" or "name_2xyz"
					if (exists($reads{$splitline1[0]})){                      											  #If one of the PE reads has already been seen then add the hit to the match hash														
						@splitline2 = split(/\t/,$reads{$splitline1[0]});												  #get the contig name from the old read
						if ($splitline[3]<= $enddist){																      #If the new read hits in the start of the contig
							$new = $splitline[2]."s";                                                                     
							$start++;																					  #Counter used to indicate that both reads hits in the start - used in the "circular" feature prediction. Start = 2 equals not proper circular
						}
						else{
							$new = $splitline[2]."e";
						}
						if ($splitline2[1]<= $enddist){																      #If the old read hits in the start of the contig
							$old = $splitline2[0]."s";
							$start++;
						}
						else{
							$old = $splitline2[0]."e";
						}			
						if ($splitline[2] gt $splitline2[0]){															  #In order to create the right hashes we take the contig with the largest NAME to be the first as e.g. Astart_Aend+Bstart_Bend is equal Bend_Bstart+Aend_Astart and we do not want to create these duplicates..
							$connectionorder = $new.":".$old;
						}
						else{
							$connectionorder = $old.":".$new;
						}								
						if ($splitline[2] ne $splitline2[0]){                                                             #If two different contigs are matched by different reads
							if (exists($connections{$connectionorder})){												  #If there are connections already increase the count - else make the first connection
								$connections{$connectionorder}++;
							}
							else{							
								$connections{$connectionorder} = 1;                                                       
							}
						}
						else{			                                                                                  #If the same contig is hit by two reads
							if ($contigs{$splitline2[0]} >= $minlength and $start == 1){								  #Circular prediction: Make sure that its not just mapping 2x the start or end and that the contig is above a minimum length
								if (exists($connections{$connectionorder})){											  #if there is connections already increase it - else make the connection
									$connections{$connectionorder}++;
									}
								else{							
									$connections{$connectionorder} = 1;
								}
							}
						}
					}
					else{								
						$reads{$splitline1[0]} = "$splitline[2]\t$splitline[3]";										   #If the other PE read has not been seen then create the first instance of the pair
					}
				}
			}			
			else{																											
			}
			$readcount++;                                                                                                   #keep track of the number of reads look though
			$printreadcount++;
			if ($printreadcount == 1000000) {
				$printreadcount = 0;
				print "$readcount reads looked through\n";
			}		
		}	
	}
}

close IN;

###########################print the data ...                                                                                  

my @tempcon1 = keys %contigs;
@tempcon1 = sort @tempcon1;

my @tempcon = keys %connections;
@tempcon = sort @tempcon;
foreach my $con (@tempcon){
	if ($connections{$con} >= $minconnections){
		@splitline = split(/:/,$con);
		@splitline1 = split(//,$splitline[0]);
		my $end1 = $splitline1[-1];
		@splitline2 = split(//,$splitline[1]);
		my $end2 = $splitline2[-1];
		pop @splitline1;	
		pop @splitline2;
		$out1 = join("",@splitline1);
		$out2 = join("",@splitline2);
		$contigswhit{$out1} = 0; 
		$contigswhit{$out2} = 0; 
		my @corder;
		push (@corder, $out1);
		push (@corder, $out2);
		@corder = sort @corder;			
		my $ccon = "$corder[0]:$corder[1]";
		if (exists($condensedhit{$ccon})){
			$condensedhit{$ccon} = $condensedhit{$ccon} + $connections{$con};	
		}
		else{
			$ccon = "$corder[0]:$corder[1]";
			if (exists($condensedhit{$ccon})){
				$condensedhit{$ccon} = $condensedhit{$ccon} + $connections{$con};	
			}
			else{
			$condensedhit{$ccon} = $connections{$con};	
			}
		}			
	}
}

print OUTconC "scaffold1\tscaffold2\tconnections\n";

foreach my $con (keys %condensedhit){
	my @splitline = split(/:/,$con);
	$splitline[0] =~ tr/ //;
	print OUTconC "$splitline[0]\t$splitline[1]\t$condensedhit{$con}\n";
}
close OUTconC;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "infile|i:s", "outcon|o:s", "enddist|e:s", "minlength|m:s", "minconnections|f:s", "avgreadlength|a:s", "headersplit|s:s");
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

	network.pl

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

	Indicates connections between scaffolds and produces network 
	files for an easy visualization 

=head1 SYNOPSIS

cytoscapeviz.pl  -i [-o -e -m -f -a -s] : version 1.0

 [-help -h]           Displays this basic usage information
 [-infile -i]         SAM format infile
 [-outcon -o]         Connections file (default: network.txt)
 [-enddist -e]        Reads must hit within this distance of the end to be considered an end hit (default: 500)
 [-minlength -m]      The minimum lenght of a contig to infer closed circle (default: 3000)
 [-minconnections -f] The minimum number of connections to print the connection to the outfile (default: 10)
 [-avgreadlength -a]  The average readlength of the reads used to estimate read position (default: 125)
 [-headersplit -s]    The symbol used to split the header to make the read 1 and 2 headers identical (default: _)
	  
=cut
