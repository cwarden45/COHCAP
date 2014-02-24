#Code written by Charles Warden (cwarden@coh.org, x60233)

use warnings;
use strict;
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
	
my $expr_file = $ARGV[0];
my $inputfile = $ARGV[1];
my $mDeU_file = $ARGV[2];
my $mUeD_file = $ARGV[3];
my $expr_pvalue_cutoff = $ARGV[4];
my $expr_fdr_cutoff = $ARGV[5];
my $expr_fc_cutoff = $ARGV[6];


my %gene_hash = define_gene_hash($expr_file);	
		
open(UP, ">$mDeU_file") || die("Could not open $mDeU_file!");
open(DOWN, ">$mUeD_file") || die("Could not open $mUeD_file!");

my $line_count=0;
open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
while (<INPUTFILE>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\n//g;
		$line =~ s/\r//g;
		$line_count++;
		if($line_count == 1)
			{
				print UP "$line\texpression.fold.change\texpression.pvalue\texpression.fdr\n";
				print DOWN "$line\texpression.fold.change\texpression.pvalue\texpression.fdr\n";
			}#end if($line_count == 1)
		else
			{
				my @line_info = split("\t",$line);
				my $gene = $line_info[1];
				my $methyl_status = $line_info[scalar(@line_info) - 1];
				if(defined($gene_hash{$gene}))
					{
						my $expr_text = $gene_hash{$gene};
						my @expr_info = split("\t",$expr_text);
						my $expr_fc = $expr_info[0];
						my $expr_pvalue = $expr_info[1];
						my $expr_fdr = $expr_info[2];
						if (($expr_pvalue <= $expr_pvalue_cutoff) && ($expr_fdr <= $expr_fdr_cutoff))
							{
								if (($methyl_status eq "Increased Methylation") && ($expr_fc <= -$expr_fc_cutoff))
									{
										print UP "$line\t$expr_text\n";
									}#end if ($methyl_status eq "Methylation Increase")
								elsif(($methyl_status eq "Decreased Methylation") && ($expr_fc >= $expr_fc_cutoff))
									{
										print DOWN "$line\t$expr_text\n"
									}#end elsif($methyl_status eq "Methylation Decrease")
							}#end if (($expr_pvalue < $expr_pvalue_cutoff) && ($expr_fdr < $expr_fdr_cutoff))
					}#end if(defined($gene_hash{$gene}) && ($num_sites > $num_sites_cutoff) && ($methyl_pvalue < $methyl_pvalue_cutoff) && ($methyl_fdr < $methyl_fdr_cutoff))
			}#end else
	}#end while (<INPUTFILE>)
close(INPUTFILE);
close(UP);
close(DOWN);

exit;

sub define_gene_hash
	{
		my ($inputfile) = @_;
		
		my %avg_hash;
		my %count_hash;
		my %total_hash;
		
		my $line_count=0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 my $line = $_;
				 chomp $line;
				 $line =~ s/\r//g;
				 $line =~ s/\n//g;
				 $line_count++;
				if(($line_count > 1) && ($line =~ /\w+/))
					{
						my @line_info = split("\t",$line);
						my $label = shift(@line_info);
						
						my $na_flag = 0;
						
						foreach my $value (@line_info)
							{
								if(($value eq "") || ($value eq "NA") || ($value eq "?"))
									{
										$na_flag=1;
									}
							}#endforeach my $value (@line_info)
						
						#print "$label\n";				
						unless(($label eq "") || $na_flag)
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
						$temp_arr[$i] = $temp_arr[$i] /$count;
					}#end for (my $i=0; $i< scalar(@temp_values) ; $i++)
				
				$avg_hash{$index}=join("\t",@temp_arr); 
			}#end foreach my $index (keys %total_hash)
			
		return(%avg_hash);
	}#end def define_average_hash