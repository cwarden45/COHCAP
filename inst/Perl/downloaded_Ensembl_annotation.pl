use warnings;
use strict;
use List::MoreUtils qw(uniq);

##NOTE: Function not implemented to run within R, since 1) you may need to modify code
##..... and 2) you have to download .gtf with gene symbols within annotations (this code tested for hg19 Ensembl .gtf)

##I was able to find .gtf for hg19 as follows (as of 6/20/2018)
#1) https://uswest.ensembl.org/index.html --> "Select a species" (Human, left side)
#2) Goes to https://uswest.ensembl.org/Homo_sapiens/Info/Index --> Select GRCh37 under "Other assemblies"
#3) Select "Downloads" (towards top of screen)
#4) Select "Databases" (towards right of screen)
#5) Select "GTF dumps" --> Gets to ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/

#For example, this code is designed to:
#a) define TSS range for "transcript" (although could also be done for gene, with one TSS)
#b) idetnify transcript/gene type is also in GTF (to filter pseudogene annotations)
#c) report gene in annotation and define region with TSS[length]_[gene] (for COHCAP, don't specify transcript, since slightly stagged TSS should start out as one region)

$| =1;

if (($ARGV[0] =~ /--help/) | ($ARGV[0] =~ /-h/))
	{
		die usage();
	}

my $input_table;
my $gtf;
my $tss_length = 2000;
my $output_table;

foreach my $arg (@ARGV)
	{
		if ($arg =~ /--beta=/)
			{
				#redefine outputfile
				($input_table) = ($arg =~ /--beta=(.*)/);
			}#end if ($arg =~ /--beta=/)

		if ($arg =~ /--gtf=/)
			{
				#redefine outputfile
				($gtf) = ($arg =~ /--gtf=(.*)/);
			}#end if ($arg =~ /--gtf=/)

		if ($arg =~ /--TSS=/)
			{
				#redefine outputfile
				($tss_length) = ($arg =~ /--TSS=(.*)/);
			}#end if ($arg =~ /--TSS=/)
			
		if ($arg =~ /--output=/)
			{
				#redefine outputfile
				($output_table) = ($arg =~ /--output=(.*)/);
			}#end if ($arg =~ /--output=/)
	}#end foreach my $arg (@ARGV)

unless(defined($input_table)){
	print "You didn't specify a methylation table with --beta=[project_methylation_table.txt] ~!\n";
	die usage();
}

unless(defined($gtf)){
	print "You didn't specify a GTF annotation with --gtf=[Homo_sapiens.GRCh37.87.gtf] ~!\n";
	die usage();
}

unless(defined($output_table)){
	print "You didn't specify an output annotation file with --output=[project_site_annotation.txt] ~!\n";
	die usage();
}
	
##create gene hash
print "Step #1: Creating gene hash from GTF\n";
my %transcript_TSS_hash;

open(IN, "$gtf") || die("Could not open $gtf!");

while(<IN>){
	my $line = $_;
	chomp $line;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	
	if (not($line =~ /^#/)){
		my @line_info = split("\t", $line);
		my $anno_type = $line_info[2];

		if ($anno_type eq "transcript"){
			my $anno_info = $line_info[8];
			my ($gene_symbol) = ($anno_info =~ /gene_name \"(\S+)\"/);
			my ($gene_type) = ($anno_info =~ /gene_biotype \"(\S+)\"/);
			
			if (not($gene_type =~ /pseudogene/)){
				#assume Bismark table is in UCSC format (so, need to add "chr" if parsing Ensembl file)
				my $chr = "chr".$line_info[0];
				my $start = int($line_info[3]);
				my $stop = int($line_info[4]);
				my $strand = $line_info[6];
				
				my $TSS;
				if($strand eq "+"){
					$TSS = $start;
				}elsif($strand eq "-"){
					$TSS = $stop;
				}
				
				if(defined($TSS)){					 
					 if(exists($transcript_TSS_hash{$chr}{$TSS})){
						if(not($transcript_TSS_hash{$chr}{$TSS} =~ /$gene_symbol/)){
							#print "Multiple genes have same TSS:\n";
							$transcript_TSS_hash{$chr}{$TSS} = $transcript_TSS_hash{$chr}{$TSS}.";".$gene_symbol;
							#print $transcript_TSS_hash{$chr}{$TSS},"\n";
						}#end if(not($transcript_TSS_hash{$chr}{$TSS} =~ /$gene_symbol/))
					 }else{
						$transcript_TSS_hash{$chr}{$TSS}=$gene_symbol;
					 }#end else
				}#end if(defined($TSS))
			}#end if (not($gene_type =~ /pseudogene/)
		}#end if ($anno_type eq "gene")
	}#end if (not($line =~ /^#/))
}#end while(<IN>)

close(IN);

##check for overlap with sites
print "Step #2: Checking for overlap with sites\n";
print "WARNING: Code may take a while if you don't have Targeted BS-Seq data, with a relatively small frequency of targeted sites\n";

open(OUT, "> $output_table") || die("Could not open $output_table!");
print OUT "SiteID\tChr\tLoc\tGene\tIsland\n";

open(IN, "$input_table") || die("Could not open $input_table!");

my $line_count = 0;

while(<IN>){
	my $line = $_;
	chomp $line;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	
	$line_count++;
	
	if ($line_count > 1){
		my @line_info = split("\t", $line);
		my @site_info = split(":",$line_info[0]);
		my $chr = $site_info[0];
		my $pos = int($site_info[1]);
		
		my $gene = "NA";
		my $island = "NA";
		
		my $chr_hash_ref;
		if(exists($transcript_TSS_hash{$chr})){
			$chr_hash_ref = $transcript_TSS_hash{$chr};
		}else{
			if($chr eq "chrM"){
				#discrepancy with Ensembl annotations
				#remember to add back in "chr", for earlier fix
				my $temp_chr = "chrMT";
				$chr_hash_ref = $transcript_TSS_hash{$temp_chr};
			}else{
				print "Issue mapping chromosome $chr in \n";
				exit;			
			}#end else
		}#end else
		my %chr_hash = %$chr_hash_ref;
		
		foreach my $chr_pos (keys %chr_hash){
			$chr_pos = int($chr_pos);
			my $test_start = $chr_pos - $tss_length;
			my $test_stop = $chr_pos + $tss_length;
			
			if(($pos >= $test_start) and ($pos <= $test_stop)){
				my $test_gene = $chr_hash{$chr_pos};
				if($gene eq "NA"){
					$gene = $test_gene;
				}elsif(not($gene =~ /$test_gene/)){
					$gene = $gene.";".$test_gene;				
				}
			}#end if(($pos >= $test_start) and ($pos <= $test_stop))
		}#end foreach my $chr_pos (keys %chr_hash)
		
		my @gene_arr = split(";",$gene);
		#use example from https://perlmaven.com/unique-values-in-an-array-in-perl
		@gene_arr = uniq @gene_arr;
		@gene_arr = sort @gene_arr;
		$gene = join(";",@gene_arr);
		
		if($gene ne "NA"){
			$island = "TSS$tss_length"."_$gene";
			$island =~ s/;/_/g;
		}#end if($gene ne "NA")
		
		print OUT "$line_info[0]\t$chr\t$pos\t$gene\t$island\n";
	}#end if ($line_count > 1)
}#end while(<IN>)

close(IN);

close(OUT);
	
exit;

sub usage{
  print <<EOF

  Usage: perl downloaded_Ensembl_annotation.pl --beta=[project]_methyl_table.txt --gtf=Homo_sapiens.GRCh37.87.gtf --TSS=2000 --output=[project]_Ensembl_hg19_site_annotation.txt

  --beta : Tab-delimited text file with Chr:Pos site ID in 1st column
				(can be created using COHCAP.BSSeq_V2.methyl.table() )
  
  --gtf : Ensembl Genomic Annotation (or similarly formatted .gtf)
				(includes gene symbol and some flags for filters)

  --TSS : Flanking length to define TSS
				
  --output : COHCAP Annotation table (for BS-Seq data, using observed sites)
							(to use with COHCAP.annotate() )  

EOF
    ;
  exit 1;
}