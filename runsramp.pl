$fi=$ARGV[0];
$fo=$ARGV[1];
$mode=$ARGV[2];
if ($fi!~/\w+/)
{
print "Usage:\n1. Go to the install folder where this script is located\n 2. perl runsramp.pl \[INFILE\] \[OUTFILE\] \[MODE\], \[MODE\] should be specified as full or mature\n";
print "Example: perl runsramp.pl example.fa exampleout.txt full\n";
die;
}
if ($fo!~/\w+/)
{
print "Usage:\n1. Go to the install folder where this script is located\n 2. perl runsramp.pl \[INFILE\] \[OUTFILE\] \[MODE\], \[MODE\] should be specified as full or mature\n";
print "Example: perl runsramp.pl example.fa exampleout.txt full\n";
die;
}
if ($mode!~/\w+/)
{
print "Usage:\n1. Go to the install folder where this script is located\n 2. perl runsramp.pl \[INFILE\] \[OUTFILE\] \[MODE\], \[MODE\] should be specified as full or mature\n";
print "Example: perl runsramp.pl example.fa exampleout.txt full\n";
die;
}

if ($mode eq "full"){$wlb=30;$wlc=126;$kmax=3;$tresX=0.672;$tresH=0.600;$tresM=0.557;$tresL=0.528;}
elsif ($mode eq "mature"){$wlb=21;$wlc=55;$kmax=3;$tresX=0.676;$tresH=0.620;$tresM=0.584;$tresL=0.557;}
else {die "Error:\[MODE\] should be either \'full\' \(for the full-transcript\/pre-mRNA mode\) or \'mature\' \(for the mature mRNA\/cDNA mode\)\n";}
#Read Sequences;
@inseq=();@inhead=();$s=0;
open (IN,$fi) || die "Error: Cannot open input file\n";
while ($line=<IN>)
{
	chomp $line;
	if ($line=~/\>(.+)/){$s++;$inhead[$s]=substr($1,0,20);}
	else
	{
		$line=~s/T/U/g;
		$line=uc($line);
		if ($line=~/[^AUCGN]+/){die "Error: Non-standard charcter (other than A,C,G,T,U,N detected\n";}
		$seq[$s].=$line;
	}
}
close IN;
if ($s==0){die "Error: No FASTA input found\n";}
#Get fragments;

@frag=();@fragh=();@fragpis=();$f=0;
for ($q=1;$q<=$s;$q++)
{
	$flanking="N" x 1000;
	$seqsub=$flanking.$seq[$q].$flanking;
	$seqlen=length($seqsub);
	for ($i=1000;$i<=$seqlen-1002;$i++)
	{
		$motif=substr($seqsub,$i,5);
		if ($motif=~/[AGU][AG]AC[ACU]/)
		{
			$f++;
			$fragh[$f]=$inhead[$q];
			$frag[$f]=substr($seqsub,$i-998,2001);
			#$shortfrag=substr($seqsub,$i-20,45);
			#$shortfrag=~s/N/\-/g;
			$fragpis[$f]=$i-997;
			$i++;
		}
	}
}
#if ($f==0){die "Error: No m6A consensus motif found\n";}
print "There are $f DRACH consensus motifs found in the $s input sequences\n";
#Encoding;
$fb=$fi."\.binarytmp";
$fc=$fi."\.spectmp";
open (FB,">$fb");
$vv=($wlb*2+1)*4;
print FB "class";
for ($h=1;$h<=$vv;$h++){print FB " v$h";}
print FB "\n";
open (FC,">$fc");
$vv=4*4*($kmax+1);
print FC "class";
for ($h=1;$h<=$vv;$h++){print FC " v$h";}
print FC "\n";
for ($q=1;$q<=$f;$q++)
{
	$binary="0";
	for ($i=1000-$wlb;$i<=1000+$wlb;$i++)
	{
		$a=substr($frag[$q],$i,1);
		if ($a eq "A"){$binary.=" 1 0 0 0";}
		if ($a eq "U"){$binary.=" 0 1 0 0";}
		if ($a eq "G"){$binary.=" 0 0 1 0";}
		if ($a eq "C"){$binary.=" 0 0 0 1";}
		if ($a eq "N"){$binary.=" 0 0 0 0";}
	}
	print FB "$binary\n";
	%spec={};
	@word=qw(A U C G);
	$pp=substr($frag[$q],1000-$wlc,$wlc*2+1);$pp=~s/N//g;
	$pplen=length($pp);
	$skencode="0";
	for ($k=0;$k<=$kmax;$k++)
	{
		for ($i=0;$i<=$pplen-$k-2;$i++)
		{
			$j=$i+$k;
			$a=substr($pp,$i,1).$k.substr($pp,$i+$k+1,1);
			$spec{$a}=$spec{$a}+1/($pplen-$k-1);
		}
	}
	for ($k=0;$k<=$kmax;$k++)
	{
		for ($i=0;$i<=3;$i++)
		{
			for ($j=0;$j<=3;$j++)
			{
				$b=0;
				$a=$word[$i].$k.$word[$j];
				$b=int($spec{$a}*1000)/1000;
				$skencode.=" $b";
			}
		}
	}
	print FC "$skencode\n";
}
close FB,FC;
print "Encoding Finished\n";
if ($mode eq "full")
{
$Rcmdline="Rscript ./sramp-full.R ".$fi." >gc.tmp";
system($Rcmdline);
}
if ($mode eq "mature")
{
$Rcmdline="Rscript ./sramp-mature.R ".$fi." >gc.tmp";
system($Rcmdline);
}
$fpv=$fi."\.pvtmp";
open (FO,">$fo");
print FO "Seq_ID\tPosition\tScore(Binary)\tScore(Spectrum)\tScore(Combined)\tClassification\n";
open (PV,$fpv) || die "Error: The predictor reported nothing, please check the integrity of the downloaded files and run the program in the install folder again\n";
$pvline=<PV>;
for ($q=1;$q<=$f;$q++)
{
	$pvline=<PV>;chomp $pvline;@pvtemp=split(" ",$pvline);
	$deci="Non-m6A site";
	$pv1=substr($pvtemp[1],0,5);$pv2=substr($pvtemp[2],0,5);$pv3=substr($pvtemp[3],0,5);
	if ($pv3>=$tresL){$deci="m6A site (Low confidence)";}
	if ($pv3>=$tresM){$deci="m6A site (Moderate confidence)";}
	if ($pv3>=$tresH){$deci="m6A site (High confidence)";}
	if ($pv3>=$tresX){$deci="m6A site (Very high confidence)";}
	print FO "$fragh[$q]\t$fragpis[$q]\t$pv1\t$pv2\t$pv3\t$deci\n";
}
close FO,PV;
`rm $fb`;
`rm $fc`;
`rm $fpv`;
print "Done\!\n";
