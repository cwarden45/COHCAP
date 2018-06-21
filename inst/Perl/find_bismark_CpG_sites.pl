#NOTE: this script may not work with more recently processed samples

use warnings;
use strict;

$| =1;

my $os = $^O;
my $os_name;
if (($os eq "MacOS")||($os eq "darwin")||($os eq "linux"))
	{
		#Mac
		$os_name = "MAC";
	}#end if ($os eq "MacOS")
elsif ($os eq "MSWin32")
	{
		#PC
		$os_name = "PC";
	}#end if ($os eq "MacOS")
else
	{
		print "Need to specify folder structure for $os!\n";
		exit;
	}#end else
	
my $root = $ARGV[0];
my $combined_file = $ARGV[1];
my $GENCODE_file = $ARGV[2];
my $CpG_Islands = $ARGV[3];
my $COHCAP_annotation  = $ARGV[4];
my $shore_length = $ARGV[5];

my $CpG_Islands_annotated = "";

if($os_name eq "PC")
	{
		$CpG_Islands_annotated  = "$root\\UCSC_CpG_Islands_with_Genes.bed";
	}
elsif($os_name eq "MAC")
	{
		$CpG_Islands_annotated  = "$root/UCSC_CpG_Islands_with_Genes.bed";
	}

my @bed_files = ();
my @sample_names = ();

opendir DH, $root or die "Failed to open $root: $!";
my @files = readdir(DH);

foreach my $file (@files)
	{
		my $test_file = "";
		if($os_name eq "PC")
			{
				$test_file = "$root\\$file";
			}
		elsif($os_name eq "MAC")
			{
				$test_file = "$root/$file";
			}
		if(-f ($test_file) && ($file =~ /.bed$/))
			{
				my ($sampleID)= ($file =~ /(.*).bed/);
				push(@sample_names,$sampleID);
				push(@bed_files,$test_file);
			}#if(-f ($test_file) && ($file =~ /.bed$/))
	}#end foreach my $file (@files)

closedir(DH);

combine_files(\@sample_names, \@bed_files, $combined_file);

add_genes_to_CpG_sites($CpG_Islands, $CpG_Islands_annotated, $GENCODE_file);

create_annotation_file($combined_file, $COHCAP_annotation, $CpG_Islands_annotated, $shore_length);

unlink($CpG_Islands_annotated);

exit;


sub add_genes_to_CpG_sites
	{
		my ($inputfile, $outputfile, $gene_file)=@_;
		
		my %gene_hash = define_gene_hash($gene_file);
		
		open(OUTPUTFILE, ">$outputfile") || die("Could not open $outputfile!");
		print OUTPUTFILE "Chr\tStart\tStop\tGenes\n";
		
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				my $line = $_;
				chomp $line;
				my @line_info = split("\t",$line);
				my $test_chr = $line_info[0];
				my $test_start = $line_info[1] - 2000;
				my $test_stop = $line_info[2] + 2000;
				#print "$test_chr|$test_start|$test_stop\n";
				my $genes = "NA";
				if(defined($gene_hash{$test_chr}))
					{
						my $hash_ref = $gene_hash{$test_chr};
						my %chr_hash = %$hash_ref;
						#print scalar(keys %chr_hash),"\n";
						foreach my $id (keys %chr_hash)
							{
								my @hash_info = split(",",$chr_hash{$id});
								my $start = $hash_info[0];
								my $gene = $hash_info[1];
								#print "|$start|$gene|$test_start|$test_stop|\n";
								if(( $start >= $test_start) && ( $start <= $test_stop))
									{
										#print "$gene\n";
										if($genes eq "NA")
											{
												$genes = $gene
											}
										else
											{
												unless($genes =~ /$gene/)
													{
														$genes.= ",$gene";
													}
											}
									}#end if (($chr eq $test_chr) && ($test_pos >= $start) && ($test_pos <= $stop))
							}#end foreach my $window (keys %total_window_hash)
					}#end if(defined($window_hash{$test_chr}))
				else
					{
						next;
					}
				$line_info[3]=$genes;
				print OUTPUTFILE join("\t",@line_info),"\n";
			}#end while (<INPUTFILE>)
		close(INPUTFILE);

		close(OUTPUTFILE);
	}#end def find_illumina_window_overlap

sub define_gene_hash
	{
		my ($inputfile)=@_;
		
		my %hash;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 my $line = $_;
				 chomp $line;
				 my @line_info = split("\t",$line);
				 my $id=$line_info[1];
				 my $chr = $line_info[2];
				 my $strand = $line_info[3];
				 my $loc = -1;
				 if($strand eq "-")
					{
						$loc=$line_info[5];
					}
				else{
						$loc = $line_info[4];
					}
				my $gene = $line_info[12];				 
				 $hash{$chr}{$id}="$loc,$gene";
			}#end while (<REF>)
		close(INPUTFILE);
		
		return(%hash);
	}#end def
	
sub create_annotation_file
	{
		my ($inputfile, $outputfile, $window_file, $shore_length)=@_;
		
		my %window_hash = define_window_hash_v2($window_file, $shore_length);
		
		open(OUTPUTFILE, ">$outputfile") || die("Could not open $outputfile!");
		print OUTPUTFILE "SiteID\tChr\tLoc\tGene\tIsland\n";
		
		my $line_count=0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 $line_count++;
				 my $line = $_;
				 chomp $line;
				 if (($line_count > 1) && !($line =~ /NA/))
					{
						my @line_info = split("\t",$line);
						my $siteID = $line_info[0];
						my ($test_chr,$test_pos) = ($siteID =~ /(chr.+)\.(\d+)/);
						#print "$test_chr|$test_pos\n";
						if(defined($window_hash{$test_chr}))
							{
								my $hash_ref = $window_hash{$test_chr};
								my %chr_hash = %$hash_ref;
								foreach my $location (keys %chr_hash)
									{
										my ($start,$stop,$genes) = ($location =~ /(\d+)-(\d+):(.*)/);
										#print "$chr|$start|$stop\n";
										if(($test_pos >= $start) && ($test_pos <= $stop))
											{
												my $window = "$test_chr:$start-$stop";
												my $COHCAP_chr = $test_chr;
												$COHCAP_chr = s/chr//g;
												print OUTPUTFILE "$siteID\t$COHCAP_chr\t$test_pos\t$genes\t$window\n";
											}#end if (($chr eq $test_chr) && ($test_pos >= $start) && ($test_pos <= $stop))
									}#end foreach my $window (keys %total_window_hash)
							}#end if(defined($window_hash{$test_chr}))
						else
							{
								print "$test_chr doesn't have a chr hash!\n";
							}
					}#end  if ($line_count > 1)
				 my @line_info = split("\t",$line);
			}#end while (<INPUTFILE>)
		close(INPUTFILE);

		close(OUTPUTFILE);
	}#end def find_illumina_window_overlap

sub define_window_hash_v2
	{
		my ($inputfile, $shore_length)=@_;
		
		my %hash;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 my $line = $_;
				 chomp $line;
				 my @line_info = split("\t",$line);
				 if($line =~ /^chr/)
					{
						 my $chr = $line_info[0];
						 my $start = $line_info[1] - $shore_length;
						 my $stop = $line_info[2] + $shore_length;
						 my $genes = $line_info[3];
						 my $loc = " $start-$stop:$genes";				 
						 $hash{$chr}{$loc}=0;
					}
			}#end while (<REF>)
		close(INPUTFILE);
		
		return(%hash);
	}#end def
	
sub find_window_overlap
	{
		my ($inputfile, $outputfile, $window_file)=@_;
		
		my %total_window_hash = define_window_hash($window_file);
		my %window_count_hash = %total_window_hash;
		
		my $line_count=0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 $line_count++;
				 my $line = $_;
				 chomp $line;
				 if (($line_count > 1) && !($line =~ /NA/))
					{
						my @line_info = split("\t",$line);
						my ($test_chr,$test_pos) = ($line_info[0] =~ /(chr.+)\.(\d+)/);
						#print "$test_chr|$test_pos\n";
						my $beta = ($line_info[2] + $line_info[1])/2;
						foreach my $window (keys %total_window_hash)
							{
								my ($chr, $start,$stop) = ($window =~ /(.*):(\d+)-(\d+)/);
								#print "$chr|$start|$stop\n";
								if (($chr eq $test_chr) && ($test_pos >= $start) && ($test_pos <= $stop))
									{
										$total_window_hash{$window} = $total_window_hash{$window} + $beta;
										$window_count_hash{$window} = $window_count_hash{$window} + 1;
									}#end if (($chr eq $test_chr) && ($test_pos >= $start) && ($test_pos <= $stop))
							}#end foreach my $window (keys %total_window_hash)
					}#end  if ($line_count > 1)
				 my @line_info = split("\t",$line);
			}#end while (<INPUTFILE>)
		close(INPUTFILE);
		
		open(OUTPUTFILE, ">$outputfile") || die("Could not open $outputfile!");
		print OUTPUTFILE "Window\tCpG.count\tAverage.Beta\n";
		foreach my $window (keys %total_window_hash)
			{
				my $probe_count = $window_count_hash{$window};
				my $avg_beta = "NA";
				
				if($probe_count > 0)
					{
						$avg_beta = $total_window_hash{$window} / $window_count_hash{$window};
					}
				print OUTPUTFILE "$window\t$probe_count\t$avg_beta\n";
			}#end foreach my $window (keys %total_window_hash)
		close(OUTPUTFILE);
	}#end def find_illumina_window_overlap

sub define_window_hash
	{
		my ($inputfile)=@_;
		
		my %hash;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 my $line = $_;
				 chomp $line;
				 my @line_info = split("\t",$line);
				 my $loc = "$line_info[0]:$line_info[1]-$line_info[2]";				 
				 $hash{$loc}=0;
			}#end while (<REF>)
		close(INPUTFILE);
		
		return(%hash);
	}#end def
	
sub combine_files
	{
		my ($arr_ref1, $arr_ref2, $outputfile)=@_;
		my @samples = @$arr_ref1;
		my @files = @$arr_ref2;
		
		my %loc_hash;
		my $header =join("\t",@samples);
		
		#get locations
		for (my $i=0; $i<scalar(@files); $i++)
			{
				my $inputfile = $files[$i];
				open(INPUTFILE, $inputfile) || die("Could not open $inputfile !");
				while (<INPUTFILE>)
					{
						 my $line = $_;
						 chomp $line;
						 my @line_info = split("\t",$line);
						 my $loc = "$line_info[0].$line_info[1]";
						 my $beta = $line_info[3];
						 $loc_hash{$loc}=1;
					}#end while (<INPUTFILE>)
				close(INPUTFILE);
			}#end for (my $i=0; $i<scalar($arr_ref2); $i++)
		
		#combine data
		
		for (my $i=0; $i<scalar(@files); $i++)
			{
				my %temp_loc_hash = %loc_hash;
				my $inputfile = $files[$i];
				#print "$inputfile\t",scalar(keys%loc_hash),"\n";
				open(INPUTFILE, $inputfile) || die("Could not open $inputfile !");
				while (<INPUTFILE>)
					{
						 my $line = $_;
						 chomp $line;
						 $line =~ s/\n//g;
						 $line =~ s/\r//g;
						 my @line_info = split("\t",$line);
						 my $loc = "$line_info[0].$line_info[1]";
						 my $beta = $line_info[3];
						 delete($temp_loc_hash{$loc});
						 if($i==0)
							{
								$loc_hash{$loc}=$beta;
							}
						 else
							{
								$loc_hash{$loc}="$loc_hash{$loc}\t$beta";
							}
					}#end while (<INPUTFILE>)
				close(INPUTFILE);
				
				#print "$inputfile\t",scalar(keys%temp_loc_hash),"\n";
				foreach my $missing_loc (keys %temp_loc_hash)
					{
						 if($i==0)
							{
								$loc_hash{$missing_loc}="NA";
							}
						 else
							{
								$loc_hash{$missing_loc}="$loc_hash{$missing_loc}\tNA";
							}				
					}#end 
			}#end for (my $foreach my $missing_loc (keys %temp_loc_hash)i=0; $i<scalar($arr_ref2); $i++)
		
		open(OUTPUTFILE, ">$outputfile") || die("Could not open $outputfile!");
		print OUTPUTFILE "Location\t$header\n";
		foreach my $loc (keys %loc_hash)
			{
				print OUTPUTFILE "$loc\t$loc_hash{$loc}\n";
			}#end 
		close(OUTPUTFILE);
	}#end def combine_files
