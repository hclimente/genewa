######################### TuRF-E ############################
# TuRF-E is a wrapper for creating an ensemble of TuRF   	#
# implemented in MDR package for SNP filtering      		#
# of case-control GWAS data.								#
# author:	Pengyi Yang										#
# email:	yangpy@it.usyd.edu.au							#
# date:		11 Dec. 2011								    #
# version:	1.1												#
# dependency:	mdr.jar										#
#############################################################

use strict;
use warnings;

unless (@ARGV) {
	usage();
}

my $cmdLen = @ARGV;
my $file = "null"; # input file
my $rankfile = "SNP_rank.txt";
my $output = "SNP_filtered.arff";
my $size = 50;	# aggregation size
my $top = 100;  # number of top SNPs to select
my %score;	# SNP score
my @selectedSNP;


# step 1: interrogate input parameters
for (my $i = 0 ; $i < $cmdLen; $i++) {
	if ($ARGV[$i] eq "-f") {
		$file = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-o") {
		$output = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-s") {
		$size = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-n") {
		$size = $ARGV[$i+1];
	}
  if ($ARGV[$i] eq "-t") {
    $top = $ARGV[$i+1];
  }
}

print "## parameters:\n";
print "input file: \"$file\"\n";
print "output file: \"$output\"\n";
print "ensemble size: $size\n";
print "number of SNPs to select: $top\n";


if ($file eq "null") {
	usage();
}

# step 2: execute the first filtering on the original dataset;
exeTuRF($file);


# step 3: generate permutated dataset and reapply TuRF
for (my $i = 0; $i < $size; $i++) {
	my $c = $i + 1;
	print "permutation ($c)\n";
	permutate($file);
	exeTuRF("replica.txt");
}

# step 4: clean up the temp files
system("rm replica.txt");
system("rm rank");


# step 5: rank SNPs according to the aggregation scores
my $count = keys(%score);
my @cList; # consensus list
foreach my $i (sort{$score{$a} <=> $score{$b}} keys(%score)) {
	my $value = $score{$i};
	$count--;
	#print "$i\t$value\n";
	my $rank = $count+ 1;
	$cList[$count] = "$i\t$value\t$rank\n";
}

# print consensus ranking
print "=========== consensus ranking ===========\n";
print @cList;
open (O, ">SNP_rank.txt");
print O "=========== consensus ranking ===========\n";
print O @cList;
close (O);

# step 6: output filtered SNP file in ARFF format for followup gene-gene interaction analysis
# header
open (O, ">".$output);
print O "\@relation SNP_top_$top\n\n";
open (I, $file);
my $header = <I>;
my @header = split(/\t/, $header);
my $length = @header;
my %headerIndex;
for (my $i = 0; $i < $length; $i++) {
   $headerIndex{$header[$i]} = $i;
}

my @selectedIndex;
for (my $i = 0; $i < $top; $i++) {
	my @snp = split(/\t/, $cList[$i]);
	print O "\@attribute $snp[0] numeric\n";
	$selectedIndex[$i] = $headerIndex{$snp[0]};
}

print O "\@attribute class {0, 1}\n";
print O "\n\@data\n";

while(<I>){
	chomp;
	my @tok = split(/\t/, $_);
	my $length = @tok;
	my $classLabel = $tok[$length - 1];

	for (my $i = 0; $i < $top; $i++) {
		my $idx = $selectedIndex[$i];
		print O "$tok[$idx],";
	}
	print O "$classLabel\n";
}

close(O);





##################################
##        subroutinues         ###
###############################################################

# usage information
sub usage {
	print "usage:\tperl TuRF-E.pl -f <file> [options]\n";
	print "\t  options\n";
	print "\t-o <string>: name of the output file (default = \"SNP_filtered.arff\")\n";
	print "\t-s <int>: aggregation size (default = 50)\n";
	print "\t-n <int>: select top n SNPs for further analysis (default = 100)\n";
	exit(1);
}

# permutation of dataset
sub permutate {
	my $file = shift; # file name

	my $sSize = -1;
	my $pSize = 0; # positive (case as 1)
	my $nSize = 0; # negative (control as 0)

	open(I, $file);
	while(<I>) {
		$sSize++;
		if ($_ =~ /0$/) {
			$nSize++;
		}
		if ($_ =~ /1$/) {
			$pSize++;
		}
	}
	close(I);

	open (I, $file);
	open (O, ">replica.txt");
	my $title = <I>;
	my %data;
	my $count = 0;
	print O "$title";
	while (<I>) {
		chomp;
		$data{$count} = $_;
		$count++;
	}
	close(I);

	# random sampling
	$count = 0;

	while ($count < ($sSize/2)) {
		my $r = int(rand($sSize));

		if ($data{$r}) {
			print O "$data{$r}\n";
			delete $data{$r};
			$count++;
		}
	}

	foreach my $i (values(%data)) {
		print O "$i\n";
	}
	close (O);
}

# execute TuRF and aggregate the scores
sub exeTuRF{
	my $file = shift; # file name
	system("java -Xmx1024m -jar ../../libs/turf/mdr.jar -filter=TURF $file > rank");

	open (I, "rank");
	while (<I>) {
		chomp;
		if ($_ =~ /.+\t.+\t\d+$/) {
			my @token = split(/\s+/, $_);
			if ($score{$token[0]}) {
				my $v = $score{$token[0]};
				$v += $token[1];
				$score{$token[0]} = $v;
			}
			else {
				$score{$token[0]} = $token[1];
			}
		}
	}
	close (I);
}
