#! /usr/bin/perl -w

use Bio::SeqIO;

my $usage = "[file] [format] [filtered file]\n";
my $file = shift or die $usage;
my $format = shift or die $usage;
my $result = shift or die $usage;


my $inseq = Bio::SeqIO->new( -file => "<$file", -format => $format,);
my $outseq = Bio::SeqIO->new( -file => ">$result", -format => $format,);
$i = 1;

while( my $seq = $inseq->next_seq) 
{
	print $seq->id()."\t$i\n";
	$seq->id("$i");
	$outseq->write_seq($seq);
	$i++;
}
		

