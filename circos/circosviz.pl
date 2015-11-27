#!/usr/bin/env perl
###############################################################################
#
#    circosviz.pl version 2.0
#    
#    Indicates connections between contigs/scaffolds and produces circos 
#    files for visualization 
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
#use POSIX;            #to use the floor command

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

my $samfile;
my $fastafile;
my $enddist;
my $minlength;
my $countbinsize;
my $headersplit;
my $circoscbin;
my $minpesplit;
my $frc;

$samfile = &overrideDefault("data.sam",'samfile');
$fastafile = &overrideDefault("data.fasta",'fastafile');
$enddist = &overrideDefault("500",'enddist');
$minlength = &overrideDefault("3000",'minlength');
$countbinsize = &overrideDefault("500",'countbinsize');
$headersplit = &overrideDefault("_",'headersplit');
$circoscbin = &overrideDefault("10000",'circoscbin');
$minpesplit = &overrideDefault("1000",'minpesplit');
$frc = &overrideDefault("N",'frc');

my $old;
my $new;
my $count = 0;
my $connectionorder;
my $readcount = 0;
my $printreadcount = 0;
my $dummy = 0;
my $seq2;
my $contigcolor;
my $bandstart;
my $bandend;
my $firstn = 0;
my $contiglength;
my $bandcount;
my $countcends = 0;
my $countcdiff = 0;
my $countcsamewl = 0;
my $numcontigs = 0;
my $totallength = 0;
my $readl;
my @splitsam;
my @oldread;
my @splitreadname;
my @contigname;
my @contiglength;
my @seqarray;
my %contigs;
my %reads;
my %breads;
my %cov;
my %counts;
my %karyotypeh;
my %count_e;
my %count_d;
my %count_s;
my %link_e;
my %link_d;
my %link_s;

######################################################################
# CODE HERE
######################################################################

open(INsam, $samfile) or die("Cannot open $samfile\n");
open(INfasta, $fastafile) or die("Cannot open $fastafile\n");
open(OUTcends, ">circos.ends.txt") or die("Cannot create circos.ends.txt\n");   
open(OUTcdiff, ">circos.dcontigs.txt") or die("Cannot create circos.dcontigs.txt\n");   
open(OUTcsamewl, ">circos.scontigs.wl.txt") or die("Cannot create circos.scontigs.wl.txt\n");   
open(OUTckary, ">circos.karyotype.txt") or die("Cannot create circos.karyotype.txt\n");   
open(OUTccov, ">circos.coverage.txt") or die("Cannot create circos.coverage.txt\n");   
open(OUTcgc, ">circos.gc.txt") or die("Cannot create circos.gc.txt\n");   

print "Generating overview of pe connections in the SAM file.\n";

while ( my $line = <INsam> ) {
	chomp $line;       
	#$line =~ s/\r\n//g;	
	@splitsam = split(/\t/, $line); 	
	if ($line =~ m/\@SQ/) {                               												                                             # If we are in the contig header area then retrive all contigs/scaffolds and store the name and length in the hash: contig		
			@contigname = split(/:/, $splitsam[1]);                     												                             # Retrive the contig name
			@contiglength = split(/:/, $splitsam[2]);																	                             # Retrive the contig length
			$contigs{$contigname[1]} = $contiglength[1];                   												                             # Make a hash with key = "contig name" and value = "contig length"			
			$totallength = $totallength + $contiglength[1];                   											
			$numcontigs++;
			for (my $count = 1; $count <= $contiglength[1]; $count++)  {                                                                                           # Add coverage
				$cov{$contigname[1]}{$count} = 0;
			}
		}	
	else {
		if ($line !~ m/(\@PG|\@HD|\@SQ|\@RG)/) {                                                                                                     # If we are in the read section of the SAM file
			my $readname = $splitsam[0];
			my $flag = $splitsam[1];
			my $contig = $splitsam[2];
			my $startpos = $splitsam[3];
			my $readlength = length($splitsam[9]);
			my $endpos = $startpos + $readlength-1;
			my $mid_pos = sprintf("%.0f", $startpos+$readlength/2);						
			@splitreadname = split(/$headersplit/, $readname);	
			my $pename = $splitreadname[0];                                                                                                          # We extract the part of the header that is similar between paired-end reads.

			for (my $count = $startpos; $count <= $endpos; $count++)  {                                                                              # Add read coverage
				$cov{$contig}{$count}++;
			}
			
			if ($flag != 19 and $flag != 35){													                 		                             # SAM Flags that indicate that these PE reads are maping as they are supposed to.. hence they are not interesting..																				
				if (($endpos <= $enddist) or ($endpos >= ($contigs{$contig}-$enddist)) and !exists($breads{$pename})) {                          # The read is required to hit within a certain distance from the contigs ends.. The middle postition of the read is used					
					if (exists($reads{$pename})){	                         											                             # If one of the PE reads has already been seen; then add the hit to the match hash														
						@oldread = split(/\t/,$reads{$pename});	 
						my $old_contig = $oldread[0];
						my $old_startpos = $oldread[1];
						my $old_readlength = $oldread[2];
						my $old_endpos = $old_startpos + $old_readlength-1;
						my $old_mid_pos = sprintf("%.0f", $old_startpos+$old_readlength/2);
						delete $reads{$pename};
						if ($contig ne $old_contig){                                                                                                 # Good connection diff contigs																				
							$link_e{$readcount} = "$contig $startpos $endpos $old_contig $old_startpos $old_endpos id=$readcount";
							if (!exists($count_e{$contig}{$mid_pos})){$count_e{$contig}{$mid_pos} = 1;}else{$count_e{$contig}{$mid_pos}++;}
							if (!exists($count_e{$old_contig}{$old_mid_pos})){$count_e{$old_contig}{$old_mid_pos} = 1;}else{$count_e{$old_contig}{$old_mid_pos}++;}
							$countcends++;	
						}
						else{			                                                                                  
							if ($contigs{$old_contig} >= $minlength and abs($startpos-$old_startpos)>= $minpesplit){                                 # Good circular connection								
								$link_e{$readcount} = "$contig $startpos $endpos $old_contig $old_startpos $old_endpos id=$readcount";							
								if (!exists($count_e{$contig}{$mid_pos})){$count_e{$contig}{$mid_pos} = 1;}else{$count_e{$contig}{$mid_pos}++;}
								if (!exists($count_e{$old_contig}{$old_mid_pos})){$count_e{$old_contig}{$old_mid_pos} = 1;}else{$count_e{$old_contig}{$old_mid_pos}++;}
								$countcends++;	
							}
						}
					}
					else{								
						$reads{$pename} = "$contig\t$startpos\t$readlength";										                                 #If the other PE read has not been seen then create the first instance of the pair
					}
				}
				else{                                                                                                                                # One of the reads is a read with BAD values
					if (exists($reads{$pename}) or exists($breads{$pename})){                      			 
						if (exists($reads{$pename})){						
							@oldread = split(/\t/,$reads{$pename});												                                     # Retrive the other read
							delete $reads{$pename};
						}
						else{
							@oldread = split(/\t/,$breads{$pename});
							delete $breads{$pename};
						}
						my $old_contig = $oldread[0];
						my $old_startpos = $oldread[1];
						my $old_readlength = $oldread[2];
						my $old_endpos = $old_startpos + $old_readlength;
						my $old_mid_pos = sprintf("%.0f", $old_startpos+$old_readlength/2);
						if ($contig ne $old_contig){                                                                                                 # Bad connection different contigs
							$link_d{$readname} = "$contig $startpos $endpos $old_contig $old_startpos $old_endpos id=$readcount";							
							if (!exists($count_d{$contig}{$mid_pos})){$count_d{$contig}{$mid_pos} = 1;}else{$count_d{$contig}{$mid_pos}++;}
							if (!exists($count_d{$old_contig}{$old_mid_pos})){$count_d{$old_contig}{$old_mid_pos} = 1;}else{$count_d{$old_contig}{$old_mid_pos}++;}
							$countcdiff++;	
						}
						else{			                                                                                  
							if ($contigs{$old_contig} >= $minlength and abs($startpos-$old_startpos)>= $minpesplit){                                 # Bad connections same contigs
								$link_s{$readname} = "$contig $startpos $endpos $old_contig $old_startpos $old_endpos id=$readcount";							
								if (!exists($count_s{$contig}{$mid_pos})){$count_s{$contig}{$mid_pos} = 1;}else{$count_s{$contig}{$mid_pos}++;}
								if (!exists($count_s{$old_contig}{$old_mid_pos})){$count_s{$old_contig}{$old_mid_pos} = 1;}else{$count_s{$old_contig}{$old_mid_pos}++;}
								$countcsamewl++;	
							}
						}
					}
					else{								
						$breads{$pename} = "$contig\t$startpos\t$readlength";									    #If the other PE read has not been seen then create the first instance of the pair -
					}
				}
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
close INsam;


print "Generating coverage files.\n";

foreach my $contig_id (sort keys %cov){	                                                                                                                 # General read coverage for each contig
	my $binsize = 0;
	my $bincount = 0;
	for (my $pos_id = 1; $pos_id <= $contigs{$contig_id}; $pos_id++)  {                                                                                  # Run through the entire contig
		$binsize++;
		$bincount = $bincount + $cov{$contig_id}{$pos_id};
		if ($binsize == $circoscbin or $pos_id == $contigs{$contig_id}){
			my $start_pos = $pos_id - $binsize;
			my $avg_coverage = sprintf("%.3f", $bincount/$binsize);
			print OUTccov "$contig_id\t$start_pos\t$pos_id\t$avg_coverage\n";
			$binsize = 0;
			$bincount = 0;
		}
	}
}

close OUTccov;

foreach my $read_pair (keys %link_e){	                                                                                                             # Printing end links	
	my $bincount1 = 0;
	my $bincount2 = 0;
	my @readstats = split(/ /,$link_e{$read_pair});
	my $r1_start = sprintf("%.0f",($readstats[1]+$readstats[2])/2 - $countbinsize/2);	                                                             # Read 1
	my $r1_end = sprintf("%.0f",($readstats[1]+$readstats[2])/2 + $countbinsize/2);
	if ($r1_start < 1){
		$r1_start = 1;
	}
	if($r1_end > $contigs{$readstats[0]}){
		$r1_end = $contigs{$readstats[0]};
	}
	for (my $pos_id = $r1_start; $pos_id <= $r1_end; $pos_id++)  {                                                                                   # Extracting the local coverage of similar type links                                                           
		if(exists($count_e{$readstats[0]}{$pos_id})){	
			$bincount1 = $bincount1 + $count_e{$readstats[0]}{$pos_id};
		}
	}
	my $r2_start = sprintf("%.0f",($readstats[4]+$readstats[5])/2 - $countbinsize/2);	                                                             # Read 2
	my $r2_end = sprintf("%.0f",($readstats[4]+$readstats[5])/2 + $countbinsize/2);
	if ($r2_start < 1){
		$r2_start = 1;
	}
	if($r2_end > $contigs{$readstats[3]}){
		$r2_end = $contigs{$readstats[3]};
	}
	for (my $pos_id = $r2_start; $pos_id <= $r2_end; $pos_id++)  {                                                                                   
		if(exists($count_e{$readstats[3]}{$pos_id})){	
			$bincount2 = $bincount2 + $count_e{$readstats[3]}{$pos_id};
		}
	}	
	print OUTcends "$link_e{$read_pair},r1linkcount=$bincount1,r2linkcount=$bincount2\n";
}

close OUTcends;

foreach my $read_pair (keys %link_d){	                                                                                                             # Printing bad links on different contigs
	my $bincount1 = 0;
	my $bincount2 = 0;
	my @readstats = split(/ /,$link_d{$read_pair});
	my $r1_start = sprintf("%.0f",($readstats[1]+$readstats[2])/2 - $countbinsize/2);	                                                             # Read 1
	my $r1_end = sprintf("%.0f",($readstats[1]+$readstats[2])/2 + $countbinsize/2);
	if ($r1_start < 1){
		$r1_start = 1;
	}
	if($r1_end > $contigs{$readstats[0]}){
		$r1_end = $contigs{$readstats[0]};
	}
	for (my $pos_id = $r1_start; $pos_id <= $r1_end; $pos_id++)  {                                                                                   # Extracting the local coverage of similar type links                                                           
		if(exists($count_d{$readstats[0]}{$pos_id})){	
			$bincount1++;
		}
	}
	my $r2_start = sprintf("%.0f",($readstats[4]+$readstats[5])/2 - $countbinsize/2);	                                                             # Read 2
	my $r2_end = sprintf("%.0f",($readstats[4]+$readstats[5])/2 + $countbinsize/2);
	if ($r2_start < 1){
		$r2_start = 1;
	}
	if($r2_end > $contigs{$readstats[3]}){
		$r2_end = $contigs{$readstats[3]};
	}
	for (my $pos_id = $r2_start; $pos_id <= $r2_end; $pos_id++)  {                                                                                   
		if(exists($count_d{$readstats[3]}{$pos_id})){	
			$bincount2++;
		}
	}	
	print OUTcdiff "$link_d{$read_pair},r1linkcount=$bincount1,r2linkcount=$bincount2\n";
}

close OUTcdiff;

foreach my $read_pair (keys %link_s){	                                                                                                             # Printing bad links on same contig
	my $bincount1 = 0;
	my $bincount2 = 0;
	my @readstats = split(/ /,$link_s{$read_pair});
	my $r1_start = sprintf("%.0f",($readstats[1]+$readstats[2])/2 - $countbinsize/2);	                                                             # Read 1
	my $r1_end = sprintf("%.0f",($readstats[1]+$readstats[2])/2 + $countbinsize/2);
	if ($r1_start < 1){
		$r1_start = 1;
	}
	if($r1_end > $contigs{$readstats[0]}){
		$r1_end = $contigs{$readstats[0]};
	}
	for (my $pos_id = $r1_start; $pos_id <= $r1_end; $pos_id++)  {                                                                                   # Extracting the local coverage of similar type links                                                           
		if(exists($count_s{$readstats[0]}{$pos_id})){	
			$bincount1++;
		}
	}
	my $r2_start = sprintf("%.0f",($readstats[4]+$readstats[5])/2 - $countbinsize/2);	                                                             # Read 2
	my $r2_end = sprintf("%.0f",($readstats[4]+$readstats[5])/2 + $countbinsize/2);
	if ($r2_start < 1){
		$r2_start = 1;
	}
	if($r2_end > $contigs{$readstats[3]}){
		$r2_end = $contigs{$readstats[3]};
	}
	for (my $pos_id = $r2_start; $pos_id <= $r2_end; $pos_id++)  {                                                                                   
		if(exists($count_s{$readstats[3]}{$pos_id})){	
			$bincount2++;
		}
	}	
	print OUTcsamewl "$link_s{$read_pair},r1linkcount=$bincount1,r2linkcount=$bincount2\n";
}

close OUTcsamewl;

print "Generating karyotype file based on fasta file.\n";

while (my $line = <INfasta>)  {
	if ($line =~ m/>/) {
		chomp $line;
		#$line =~ s/\r\n//g;
		$line =~ s/\>//g; 	
		if ($dummy == 1){
			push (@seqarray, $seq2);
			my @splitline = split(/\t/,$seq2);
			$contiglength = length($splitline[1]);
			$karyotypeh{$splitline[0]}="chr - $splitline[0] $splitline[0] 0 $contiglength set1-7-qual-";
		}
		$seq2 = "$line\t";
		$dummy =1;
	}
	else {
		chomp $line;
		#$line =~ s/\r\n//g;
		$seq2 = $seq2.$line;
		}
}

close INfasta;

push (@seqarray, "$seq2");	          #to catch the last sequence
my @splitline = split(/\t/,$seq2);
$contiglength = length($splitline[1]);
$karyotypeh{$splitline[0]}="chr - $splitline[0] $splitline[0] 0 $contiglength set1-7-qual-";                                #Change the color scheme here if needed
foreach my $key (reverse sort {$contigs {$a} <=> $contigs {$b}} keys %contigs){
	$contigcolor++;
	if ($contigcolor == 8){																			    			        #Change the color scheme here if needed
		$contigcolor = 1;
	}
	print OUTckary "$karyotypeh{$key}$contigcolor\n";	
}
print OUTckary "\n";                                                                                                        #Just to seperate for the band part of the karyotype file



print "Generating gc file based on fasta file.\n";

foreach my $sequence (@seqarray){
	$counts{G} = 0;
	$counts{C} = 0;
	$counts{A} = 0;
	$counts{T} = 0;
	$firstn = 0;
	my @tseq = split("\t", $sequence);
	my @seq = split("", $tseq[1]);
	my $bincount = 0;
	my $lengthcount = 0;
	my $lastend = 0;
	foreach my $nucleotide (@seq) {
		$counts{$nucleotide}++;
		$bincount++;
		$lengthcount++;		
		if ($bincount == $circoscbin){
			my $gc = ($counts{G}+$counts{C})/($counts{G}+$counts{C}+$counts{A}+$counts{T})*100;
			print OUTcgc "$tseq[0] $lastend $lengthcount $gc\n";
			$lastend = $lengthcount;
			$counts{G} = 0;
			$counts{C} = 0;
			$counts{A} = 0;
			$counts{T} = 0;
			$bincount = 0;
		}
		if ($nucleotide eq "N"){
			if ($firstn == 0){
				$firstn = 1;
				$bandstart = $lengthcount;
			}				
		}
		else{
			if ($firstn == 1){
				$firstn = 0;
				$bandcount++;
				print OUTckary "band $tseq[0] band$bandcount band$bandcount $bandstart $lengthcount black\n";
			}
		}
	}
	my $gc = ($counts{G}+$counts{C})/($counts{G}+$counts{C}+$counts{A}+$counts{T})*100;
	print OUTcgc "$tseq[0] $lastend $lengthcount $gc\n";	
}



if ($frc ne "N"){
	print "Generating circos formated FRCbam file.";
	open(INfrc, $frc) or die("Cannot open $frc\n");
	open(OUTfrc, ">frc.tracks.txt") or die("Cannot create frc.tracks.txt\n");   
	while ( my $line = <INfrc> ) {
		chomp $line;   	
		#$line =~ s/\r\n//g;
		next if ($line =~ m/#/);
		my @splitline = split(/ /,$line);
		if ($line =~ m/STRECH_PE/) {$splitline[1] = 1};
		if ($line =~ m/COMPR_PE/) {$splitline[1] = 2};
		if ($line =~ m/LOW_NORM_COV_PE/) {$splitline[1] = 3};
		if ($line =~ m/LOW_COV_PE/) {$splitline[1] = 4};
		if ($line =~ m/HIGH_COV_PE/) {$splitline[1] = 5};
		if ($line =~ m/HIGH_NORM_COV_PE/) {$splitline[1] = 6};
		if ($line =~ m/HIGH_OUTIE_PE/) {$splitline[1] = 7};
		if ($line =~ m/HIGH_SINGLE_PE/) {$splitline[1] = 8};
		if ($line =~ m/HIGH_SPAN_PE/) {$splitline[1] = 9};
		if ($line =~ m/COMPR_MP/) {$splitline[1] = 10};
		if ($line =~ m/HIGH_OUTIE_MP/) {$splitline[1] = 11};
		if ($line =~ m/HIGH_SINGLE_MP/) {$splitline[1] = 12};
		if ($line =~ m/HIGH_SPAN_MP/) {$splitline[1] = 13};
		if ($line =~ m/STRECH_MP/) {$splitline[1] = 14};
		print OUTfrc "$splitline[0] $splitline[2] $splitline[3] $splitline[1]\n";
	}
	close INfrc;
	close OUTfrc;
}


########################### Print stats#################################

print "Total number of reads\n";
print "$readcount\n";
print "Good normal PE connections (%relative)\n";
print $readcount/2-$countcends-$countcsamewl-$countcdiff,"(100%)\n";
print "Good end connections:\n";
print "$countcends(",sprintf("%.3f",($countcends/(2*$numcontigs*$enddist))/(($readcount/2-$countcends-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";
print "Bad connection different contigs:\n";
print "$countcdiff(",sprintf("%.3f",($countcdiff/($totallength))/(($readcount-$countcends/2-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";
print "Bad connection same contig:\n";
print "$countcsamewl(",sprintf("%.3f",($countcsamewl/($totallength))/(($readcount-$countcends/2-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";

close OUTckary;
close OUTccov;	
close OUTcgc;	

exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "samfile|i:s","fastafile|f:s", "enddist|e:s", "minlength|m:s", "countbinsize|z:s", "headersplit|s:s", "circoscbin|b:s", "minpesplit|p:s", "frc|g:s" );
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

    circosviz.pl

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

	Indicates connections between contigs/scaffolds and produces circos 
	files for visualization

=head1 SYNOPSIS

circosviz.pl  -i -f [-e -m -z -s -b -p -g] : version 1.0

 [-help -h]           Displays this basic usage information
 [-samfile -i]        SAM formated mapping file
 [-fastafile -f]      Fasta file of sequences used in the mapping
 [-enddist -e]        Reads must hit within this distance of the end to be considered an end hit (default: 500)
 [-minlength -m]      The minimum lenght of a contig to infer closed circle (default: 3000)
 [-countbinsize -z]   The count bin size used to filter links in circos afterwards (default: 500)
 [-headersplit -s]    The symbol used to split the header to make the read 1 and 2 headers identical (default: _)
 [-circoscbin -b]     The windowlength used for the circos coverage file (default: 10000) 
 [-minpesplit -p]     Minimum distance between PE reads to be considered a split pe read (default: 1000)
 [-frc -g]            Add FRCbam Features using the Feature.gff output of FRCbam 	  
=cut
