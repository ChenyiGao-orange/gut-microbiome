#!/usr/bin/perl -w

die "perl $0 <matrix_file> <map> <sortgroup> <whetherdiff(T/F)> <width> <height>\n" if(@ARGV!=6);
open IN, "$ARGV[0]";
my $head=<IN>;

my $i=0;
while(my $line=<IN>){
	my @dd=split/\t/,$line;
	$dd[0]=~s/\s+|^\[|\]|\(|\)/_/g;

	if($i%1==0){
		close OA;
		my $aa=int($i/1);
		open OA, ">$dd[0].tmp";
		print OA "$head$line";
	}
	else{
		print OA "$line";
	}
	$i++;
}
close IN;
close OA;
open IN2, "$ARGV[0]";
my $head2=<IN2>;
while(my $line2=<IN2>){
	        my @dd2=split/\t/,$line2;
		$dd2[0]=~s/\s+|^\[|\]|\(|\)/_/g;
`perl /mnt/sdb/lgq/bin/tools/row2column.pl  $dd2[0].tmp  > $dd2[0].tmp2`;
`rm  $dd2[0].tmp`;
`perl /mnt/sdb/lgq/bin/tools/boxplot/diff_box_merge.pl $ARGV[1] 1 $dd2[0].tmp2 0 > $dd2[0].xls`;
`rm  $dd2[0].tmp2 `;
`perl /mnt/sdb/lgq/bin/tools/boxplot/diff_box.pl $dd2[0].xls $ARGV[1] $dd2[0] $ARGV[2] $ARGV[3] $ARGV[4] $ARGV[5]` ;  
`rm $dd2[0].xls`; 
}
