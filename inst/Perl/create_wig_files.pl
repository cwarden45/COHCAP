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
elsif (($os eq "MSWin32")||($os eq "msys"))
	{
		#PC
		$os_name = "PC";
	}#end if ($os eq "MacOS")
else
	{
		print "Need to specify folder structure for $os!\n";
		exit;
	}#end else
	
my $inputfile = $ARGV[0];
my $outputfolder = $ARGV[1];
		
my %output_indices;
		
#Define output indicies
my $line_count=0;
open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
while (<INPUTFILE>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_info = split("\t",$line);
		$line_count++;
		for (my $i=0; $i < scalar(@line_info); $i++)
			{
				if($line_info[$i] =~ /beta/)
					{
						my $outputfile = "$line_info[$i]".".wig";
						$output_indices{$outputfile}=$i;
					}#end if($line_info[$i] =~ /beta/)
			}#end for (my $i=0; $i < scalar(@line_info); $i++)
		last;
	}#end while (<INPUTFILE>)
close(INPUTFILE);
		
#move text to .wig file
foreach my $outputfile (keys %output_indices)
	{
		my $index = $output_indices{$outputfile};
		my $complete_file = "";
		if($os_name eq "PC")
			{
				$complete_file = "$outputfolder\\$outputfile";
			}
		elsif($os_name eq "MAC")
			{
				$complete_file = "$outputfolder/$outputfile";
			}
		open(OUTPUTFILE, ">$complete_file") || die("Could not open $complete_file!");
		my ($title) = ($outputfile =~ /(.*).wig$/);
		my $y_start = 0;
		my $y_stop = 1;
		if($outputfile =~ /delta/)
			{
				$y_start=-1;
			}
		print OUTPUTFILE "track type=wiggle_0 name=\"$title\" color=0,100,0 viewLimits=$y_start:$y_stop visibility=full autoScale=off\n";
				
		$line_count=0;
		my $start_index = 0;
		my $prev_chr = 0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				my $line = $_;
				chomp $line;
				$line =~ s/\r//g;
				$line =~ s/\n//g;
				my @line_info = split("\t",$line);
				$line_count++;	 
				if($line_count > 1)
					{
						my $chr = "chr$line_info[1]";
						my $pos = $line_info[2];
						my $value = $line_info[$index];

						if($chr ne $prev_chr)
							{
								print OUTPUTFILE "variableStep chrom=$chr\n";
							}#end if($chr ne $prev_chromosome)
												
						unless($value eq "NA")
							{
								print OUTPUTFILE "$pos\t$value\n";
							}#end unless($value eq "?")
											
						$prev_chr=$chr;	
					}#end else)
			}#end while (<INPUTFILE>)
		close(INPUTFILE);
		close(OUTPUTFILE);
	}#end foreach my $outputfile (keys %outputhash)

unlink($inputfile);

exit;