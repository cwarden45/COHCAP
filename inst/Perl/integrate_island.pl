#Code written by Charles Warden (cwarden@coh.org, x60233)

use warnings;
use strict;
use diagnostics;
use Cwd 'abs_path'; 

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
	
my $expression_file = $ARGV[0];
my $methyl_file = $ARGV[1];
my $paired_file = $ARGV[2];

my %exp_hash = define_exp_hash($expression_file);
my @samples=();
		
open(OUTPUTFILE, ">$paired_file") || die("Could not open $paired_file!");
my $line_count=0;
open(INPUTFILE, $methyl_file) || die("Could not open $methyl_file!");
while (<INPUTFILE>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		$line_count++;
		if($line_count == 1)
			{
				my @line_info = split("\t",$line);
				my $island = shift(@line_info);
				my $gene = shift(@line_info);
				@samples=@line_info;
				print OUTPUTFILE "Island\tGene";
				for (my $i = 0; $i < scalar(@line_info); $i++)
					{
						my $sample = $line_info[$i];
						print OUTPUTFILE "\t$sample.Methyl\t$sample.Expr";
					}#end for (my $i = 2; $i < scalar(@line_info); $i++)
				print OUTPUTFILE "\n";
			}#end if($line_count == 1)
		else
			{
				my @line_info = split("\t",$line);
				my $island = shift(@line_info);
				my $gene = shift(@line_info);
				if(defined($exp_hash{$gene}))
					{
						print OUTPUTFILE "$island\t$gene";
						for (my $i = 0; $i < scalar(@line_info); $i++)
							{
								my $sample = $samples[$i];
								my $methyl = $line_info[$i];
								my $expr = "NA";
								if(defined($exp_hash{"$sample\t$gene"}))
									{
										$expr = $exp_hash{"$sample\t$gene"};
									}
								print OUTPUTFILE "\t$methyl\t$expr";
							}#end for (my $i = 2; $i < scalar(@line_info); $i++)
						print OUTPUTFILE "\n";
					}#end if(defined($island_hash{$island}))
			}#end else
	}#end while (<INPUTFILE>)
close(INPUTFILE);
close(OUTPUTFILE);

exit;
	
sub define_exp_hash
	{
		my ($inputfile) = @_;
		
		my %avg_hash;
		my %count_hash;
		my %total_hash;
		
		my @column_ids;
		
		my $line_count=0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 my $line = $_;
				 chomp $line;
				 $line =~ s/\r//g;
				 $line =~ s/\n//g;
				 $line_count++;
				if($line_count == 1)
					{
						$line =~ s/-/./g;
						my @line_info = split("\t",$line);
						my $label = shift(@line_info);
						@column_ids=@line_info;
					}
				elsif(($line_count > 1) && ($line =~ /\w+/))
					{
						my @line_info = split("\t",$line);
						my $label = shift(@line_info);
						
						#print "$label\n";				
						unless($label eq "")
							{
								#print "$label|$value|\n";
								#print join(",",@line_info),"\n";
								if(defined($total_hash{$label}))
									{
										my @temp_values = split("\t", $total_hash{$label});
										
										for (my $i=0; $i< scalar(@temp_values) ; $i++)
											{
												$temp_values[$i] = $temp_values[$i] + $line_info[$i];
											}#end for (my $i=0; $i< scalar(@temp_values) ; $i++)
										
										$total_hash{$label} = join("\t",@temp_values);
										$count_hash{$label} = $count_hash{$label} + 1;
									}#end if(defined($total_hash{$label}))
								else
									{
										$total_hash{$label} = join("\t",@line_info);
										$count_hash{$label} = 1;								
									}#end else
							}#end 
					}#end 		if(($line_count > 1) && ($line =~ /\w+/))
			}#end while (<INPUTFILE>)
		close(INPUTFILE);
		
		foreach my $index (keys %total_hash)
			{
				my @temp_arr = split("\t", $total_hash{$index});
				my $count = $count_hash{$index};
				
				for (my $i=0; $i< scalar(@temp_arr) ; $i++)
					{								
						my $sample = $column_ids[$i];
						if($sample =~ /^\d/)
							{
								$sample = "X$sample";
							}
						$avg_hash{$sample}=1;
						my $total_value = $temp_arr[$i];
						
						$avg_hash{"$sample\t$index"}=$total_value / $count;
					}#end for (my $i=0; $i< scalar(@temp_values) ; $i++)
				
				$avg_hash{$index}=join("\t",@temp_arr); 
			}#end foreach my $index (keys %total_hash)
			
		return(%avg_hash);
	}#end def define_average_hash