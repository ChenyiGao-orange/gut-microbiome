#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use POSIX; 

my %opts;
my $v = "v2022-11-30";
#GuoqiLiu liuguoqii@163.com
my $help;
my $col=0;
my $row=0;




GetOptions (\%opts,"i1=s","o=s","i2=s","w=i","h=i","color=s","row!"=>\$row,"col!"=>\$col);#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i1*  	   <str>          cor file  
#                -i2*	   <str>          p file  
#                -o        <str>          default output.pdf
#		 -w        <int>          the width of the figure,default:6
#		 -h	   <int>          the height of the figure,default:6
#		 -color    <str>          two or more colors to ramp. eg: blue-white-red , black-red  or green-red-yellow  darkblue-darkgreen-yellow-darkred and so on. colors Seprated by "-". default : pheatmap color 
#               -row                      row cluster , default not cluster
#               -col                      column cluster, default not cluster
#
#                Example :                Usage:perl $0  -i1 cor.txt  -i2 p.txt 
#
#########################################################################################################################################################
USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2


die $usage if (!($opts{i1}));
die $usage if (!($opts{i2}));
$opts{o}=defined$opts{o}?$opts{o}:"output.pdf";
$opts{w}=defined$opts{w}?$opts{w}:6;
$opts{h}=defined$opts{h}?$opts{h}:6;
$opts{color}=defined$opts{color}?$opts{color}:"default";
$opts{row} ||= 0 ;
$opts{col} ||= 0 ;
#$opts{scale} ||= "none" ;




sub get_time {
	 my $hhh=shift;
	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
	 print  "\@\@$gettime====>>>>>>>$hhh\n";
}
&get_time("$0 program start running...") ;
if ($row == 0) {$row = "FALSE";} else {$row = "TRUE";}
if ($col == 0) {$col = "FALSE";} else {$col = "TRUE";}
#if ($opts{scale} ne "none" && $opts{scale} ne "row" && $opts{scale} ne "column" ) {
#         print "error: -scale only input row or column or not input none \n";
#         exit; 
#} 




my @col2;
my $ttmp;
my $tmp2;
if ($opts{color} ne "default") {@col2 =  split /-/, $opts{color} ; 


 $ttmp = join('","',@col2) ;
 $tmp2 = '"'.$ttmp.'"' ; 
}



open RCMD, ">cmd.r";

print RCMD "
library(pheatmap)
library(reshape2)
library(psych)


input1 <- read.table(file=\"$opts{i1}\",header=TRUE,sep=\"\\t\",check.names = FALSE,quote = \"\",comment.char = \"\",row.names=1)
input2 <- read.table(file=\"$opts{i2}\",header=TRUE,sep=\"\\t\",check.names = FALSE,quote = \"\",comment.char = \"\",row.names=1)

#input1 <- t(input1)
#input2 <- t(input2)

input2 <- input2[rownames(input1),]
input2 <- input2[,colnames(input1)]

#correlation <- corr.test(input1,input2)








#
#
cor <- input1

p1 <- pheatmap(cor,cluster_rows = $row ,cluster_cols = $col)
";
if ($row eq "TRUE") {

print RCMD "gn=rownames(cor)[p1\$tree_row[[\"order\"]]]\n";

	}
else {
print RCMD "gn=rownames(cor)\n";
    }
if ( $col eq "TRUE" ) { 
print RCMD 	"sn=colnames(cor)[p1\$tree_col[[\"order\"]]]\n";
		}
else {
print RCMD "sn=colnames(cor)\n";
}

print RCMD "

cor <- cor[gn,sn]


#est2 <- matrix(abs(rnorm(200)), 20, 10)
p <- input2
ok <- p
ok[ok<=0.001] <- \"***\"
ok[ok > 0.001 & ok <=0.01] <- \"**\"
ok[ok > 0.01 & ok <=0.05] <- \"*\"
ok[ok>0.05] <- \" \"

ok=ok[gn,sn]
";
if ($opts{color} ne "default")  {
       print  RCMD "pheatmap(cor, display_numbers = ok  ,number_color=\"black\", fontsize_number=15,color = colorRampPalette(colors = c($tmp2))(100),cluster_rows = $row ,cluster_cols = $col,width=$opts{w},height=$opts{h},filename=\"$opts{o}\",fontsize_row = 11)
\n";}
else {
	print RCMD "pheatmap(cor, display_numbers = ok  ,number_color=\"black\", fontsize_number=15,cluster_rows = $row ,cluster_cols = $col,width=$opts{w},height=$opts{h},filename=\"$opts{o}\",fontsize_row = 11)\n";
}
close RCMD; 
	`Rscript  cmd.r` ;
