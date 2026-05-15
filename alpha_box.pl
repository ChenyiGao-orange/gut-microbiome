#!/usr/bin/perl -w
###################################################################
#####Author : GuoqiLiu                                             #
#####Date   : 2024-04-17                                           #
#####Copyright (C) 2019~ Mingke Biotechnology (Hangzhou) Co., Ltd. #
#####Contact: liuguoqi@mingkebio.com                               #
#####Suppose: Auto diversity report  program                       #
#####step :                                                        #
#####1. diver_pipe_V3.pl generate basic report                     #
#####2. modify html and data                                       #
#####3. html2pdf report                                            #
#####4. Auto-send report ? whether or not                          #
#####Platform :                                                    #
##############Ubuntu 18.04 & Windows10 ,Perl v5.26 #################
####################################################################
#Log:
##V1 version 2024-05-28 

##
use strict;
use warnings;
use Getopt::Long;
use POSIX; 

my %opts;
my $v = "v2024-05-28";

my $help;

GetOptions (\%opts,"i=s","o=s","rep=s","w=f","h=f","m=s","g=s","color=s","rlc=f","point=s","clc=f","llc=f","lg=i","angle=i","test=s");#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i*  	   <str>             matrix file  
#                -m       <str>              map file
#                -rep*     <str>             ASV/otu_rep.fasta file or tree file  
#                -g        <str>             group sort list
#                -color    <str>             group color 
#                -point    <str>             T/F, T display point ,F not display point ,default:F   
#                -angle    <int>             xlab angle 0-360 default:0
#                -rlc     <float>            xlim size default:12
#                -clc     <float>            ylim size default:12
#                -llc     <float>            legend size default:12
#                -test     <str>             two group statistic test only can choice Wilcox-test,T-test,None,default:Wilcox-test
#
#		 -w        <float>           the width of the figure,default:4
#		 -h	   <float>           the height of the figure,default:5
#
#                Example :                Usage:perl $0  -i otu/ASV_table.xls   -rep ASV/otu_rep.fasta
#
#########################################################################################################################################################

USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2


die $usage if ( !(defined $opts{i}  && $opts{rep}) );

#die $usage if (!($opts{i}));
#die $usage if (!($opts{m}));
#die $usage if (!($opts{rep}));
#

#####                -llc     <float>            legend size
#die $usage if (!($opts{i2}));
#$opts{o}=defined$opts{o}?$opts{o}:"output";
$opts{w} ||= 4;
$opts{h} ||= 5;
$opts{g} ||= "F";
$opts{test} ||= "Wilcox-test";
#$opts{top} ||= 15 ;
$opts{color} ||= "F" ;
$opts{rlc} ||= 12 ;
$opts{clc} ||= 12 ;  
$opts{llc} ||= 12 ;
$opts{angle} ||= 0 ;
$opts{point} ||= "F" ;
$opts{m} ||= "F";
#$opts{color}=defined$opts{color}?$opts{color}:"default";
#$opts{row} ||= 0 ;
#$opts{col} ||= 0 ;
#$opts{scale} ||= "none" ;



my $qiime2_2="/mnt/sdb/lgq/bin/software/py38/envs/qiime2-2022.2/bin/qiime";
my @qiime2_arr=split/\//,$qiime2_2;
my $pathqiime2=join("/",@qiime2_arr[0..($#qiime2_arr-1)]);
# print "$pathqiime2\n"; 
sub get_time {
	 my $hhh=shift;
	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
	 print  "\@\@$gettime====>>>>>>>$hhh\n";
}
&get_time("$0 program start running...") ;


#print "$opts{angle} ------------------\n";

########################################################

use FindBin qw($Bin);
chomp (my $current_dir=`pwd`);


sub alpha {
	my ($qiime2,$asvtable,$rep) = @_;
	if (-e "feature_table.biom") {
		`rm feature_table.biom`;
	         `biom convert -i $asvtable   -o feature_table.biom --table-type "OTU table" --to-json`;
	 }
	 else {
		 `biom convert -i $asvtable   -o feature_table.biom --table-type "OTU table" --to-json`;
	 }

	 `rm tmp_feature_table.xls` if (-e "tmp_feature_table.xls") ;
	 `rm tree.nwk` if (-e "tree.nwk") ;
	 open I,$asvtable;
	 open O,">tmp_feature_table.xls";
         chomp (my $sam = <I>); 
	 my @samgroup = split/\t/,$sam; 
	 print O "ID\t",join("\t",@samgroup[1..$#samgroup]),"\n";
	 while (<I>){
		 print O $_;
	 }
	 close I;close O;
	 my @alpha_arr = ("chao1","observed_features","goods_coverage","pielou_e","simpson","shannon","ace");
	 if (-d "tmp_alpha") {`rm -r tmp_alpha`;}else {mkdir "tmp_alpha" ;}
	 if (-e "tree.nwk"){`rm tree.nwk`;}
	 # mkdir "tmp_alpha"  unless(-d "tmp_alpha");
	 open Oo,">alpha.sh";
	 #print Oo "export PATH=\"/mnt/sdb/lgq/bin/software/py38/envs/qiime2-2022.2/bin:\$PATH\"\n"; 
	 #$pathqiime2
	 print Oo "export PATH=\"$pathqiime2:\$PATH\"\n";
	 #close Oo;
	 #`sh config.sh`;
	 print Oo "$qiime2 tools  import --input-path feature_table.biom   --output-path tmp_alpha/feature_table.qza --type 'FeatureTable[Frequency]'  --input-format BIOMV100Format\n";
	 foreach my $i (@alpha_arr) {
		 #`$qiime2 tools  import --input-path feature_table.biom   --output-path tmp_alpha/feature_table.qza --type 'FeatureTable[Frequency]'  --input-format BIOMV100Format`;
		 print Oo "$qiime2    diversity  alpha --i-table tmp_alpha/feature_table.qza --p-metric $i  --o-alpha-diversity tmp_alpha/$i.qza \n";
		 print Oo "$qiime2 tools  extract --input-path tmp_alpha/$i.qza --output-path tmp_alpha/$i\n";
	 }
         
	 open TMP,"$rep";
	         my $tg=0;
		         while (<TMP>){if(/^>/) {$tg=1;} last;} 

	 if ($tg == 1) {
	 print Oo "$qiime2 tools import --input-path $rep  --output-path tmp_alpha/rep-seqs.qza  --type 'FeatureData[Sequence]' \n";
	 print Oo "$qiime2 phylogeny align-to-tree-mafft-fasttree  --i-sequences tmp_alpha/rep-seqs.qza --o-alignment tmp_alpha/aligned-rep-seqs.qza --o-masked-alignment tmp_alpha/masked-aligned-rep-seqs.qza --o-tree tmp_alpha/unrooted-tree.qza  --o-rooted-tree tmp_alpha/rooted-tree.qza\n";
           print Oo "$qiime2 tools  extract --input-path tmp_alpha/rooted-tree.qza --output-path tmp_alpha/rooted-tree\n";
	   print Oo "#ln -s tmp_alpha/rooted-tree/*/data/tree.nwk .\n";
	   print Oo "cp  -a  tmp_alpha/rooted-tree/*/data/tree.nwk .\n"; 
   }
   else {`cp $rep tree.nwk`;}
           print Oo "cp $Bin/pd.py .;\n";
	   print Oo "./pd.py \n";
           close Oo;
	   `sh alpha.sh`;
   }
#&alpha($qiime2_2,$opts{i},$opts{rep}); 

sub merge_table{
my ($table,$output) = @_;
open FILE,$table or die "$!\n";
chomp (my @files = <FILE>);
close FILE;

my (%abun,%tax,%sample);

foreach my $file(@files){
	open INF,$file or die "$!\n";
	my %num2sam;
	while(<INF>){
		chomp;
		my @temp = split /\t/;
		for(my $i = 1;$i < @temp;$i++){
			if($. == 1){
				$num2sam{$i} = $temp[$i];
				$sample{$temp[$i]}++;
			}else{
				$tax{$temp[0]}++;
				if(exists $abun{$temp[0]}{$num2sam{$i}}){
					$abun{$temp[0]}{$num2sam{$i}} += $temp[$i];
				}else{
					$abun{$temp[0]}{$num2sam{$i}} = $temp[$i];
				}
			}
		}
	}
	close INF;
}


open OUT,"> $output" or die "$!\n";
print OUT "Sample ID";
foreach my $sam(sort keys %sample){
	print OUT "\t$sam";
}
print OUT "\n";

foreach my $id(sort keys %tax){
	print OUT $id;
	foreach my $sam(sort keys %sample){
		if(exists $abun{$id}{$sam}){
			print OUT "\t$abun{$id}{$sam}";
		}else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;
}


sub samsum {
my ($table,$output) = @_;

my (%id2sam,%num);
open INT,$table or die "$!\n";
while(<INT>){
	chomp;
	my @temp = split /\t/,$_;
	my $num = @temp;
	#$num-- if($flag eq "T");
	for(my $i = 1;$i < $num;$i++){
		if($. == 1){
			$id2sam{$i} = $temp[$i];
		}else{
			$num{$id2sam{$i}} += $temp[$i];
		}
	}
}
close INT;


open OUT,"> $output" or die "$!\n";
print OUT "\tReads\n";
foreach my $sam(sort {$num{$a} <=> $num{$b}} keys %num){
		print OUT "$sam\t$num{$sam}\n";
	}
close OUT;
}




sub split_file {
	my ($input,$group) = @_; 
	open GROUP,$group;
	<GROUP>;
	my %grp=();
	while (<GROUP>) {
		chomp;
		my @group=split/\t/,$_;
		$grp{$group[0]} = $group[1];
	}
	close GROUP;
	open INPUT,$input;
	open TMP,">tmp";
	chomp(my $hh=<INPUT>);
	print TMP "$hh\tgroup\n";
	while (<INPUT>) {
		chomp;
		my @arr=split/\t/,$_;
		if (exists $grp{$arr[0]}) {
			print TMP "$_\t$grp{$arr[0]}\n";
		}
		else {
			print "error $arr[0] no group!!!!\n";
		}
	}
	close TMP;close INPUT;


	open Rscript,">cmd1.r";
	print Rscript "
	options(stringsAsFactors=FALSE)
	library(plyr)
	library(reshape2)
	library(ggplot2)
	data  <- read.table(\"tmp\",header = T,sep = \"\\t\",check.names=F,comment.char=\"\",quote=\"\")
	for (i in 2:(ncol(data)-1)) {
	             write.table(data[,c(1,i,ncol(data))],file=paste(colnames(data)[i],\".tsv\",sep=\"\"),sep=\"\\t\",quote=FALSE,row.names=FALSE)          }

";

`R --no-save < cmd1.r`;
}
#		  &split_file("qiime2_alpha.xls","group");




sub Rboxplot {
	my ($input,$group,$color,$sg,$point,$test,$angle,$rlc,$llc,$clc) = @_;
	#chao1.tsv,$opts{m},$opts{color},$opts{g},$opts{point},Wilcox-test/T-test,$opts{angle},$opts{rlc},$opts{llc},$opts{clc}
	#my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#0000cd","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#00ffff","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e","#ff7f50","#6495ed","#fff8dc","#dc143c","#00ffff");
	my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e");
	my $col;
	my @ttmp;
	#=head;
	##
	#
	my %group2=();
	open GROUP,$group ;
	<GROUP>;
	while(<GROUP>){
		chomp;
		my @a=split/\t/,$_;
		# push @{$group{$a[1]}},$a[0];
		$group2{$a[1]} =  0 ;  
	}
       close GROUP;
       my $mygroupfile = $color ;
       my $sort_group;
       my @group_m ; 
       if ($sg eq "F") {@group_m =  sort (keys %group2); }
       else { open TMM,$sg;chomp(@group_m =  <TMM>); }
       $sort_group = join("\",\"",@group_m);
       $sort_group = "\"".$sort_group."\"";
       if ($mygroupfile eq "F") {
	       open T,">color.txt";
	       foreach my $i (0..$#group_m) {
		       print T "$group_m[$i]\t$mycol[$i]\n";
		       push @ttmp,"\"$group_m[$i]\""."="."\"$mycol[$i]\"";
	       }
	       $col=join(",",@ttmp);
       }
       else {
	       open C,$mygroupfile;
	       while(<C>){
		       chomp;

		       my @ab=split/\t/,$_;
		       push @ttmp,"\"$ab[0]\""."="."\"$ab[1]\"";
		        }
			  $col=join(",",@ttmp);
		  }


 #my $mp = &delete_ns_compare($input) ;
my @label=split/\.tsv/,$input;
#my $mp = &delete_ns_compare($group,$input,$test) ;
my $mp;
#Wilcox-test/T-test
if ($test eq "Wilcox-test") {$test="wilcox.test"; $mp = &delete_ns_compare($group,$input,$test) ;}
if ($test eq "T-test") {$test="t.test";$mp = &delete_ns_compare($group,$input,$test) ;}
if ($test eq "None") {$mp="" ;}
open RSCRIPT,">cmd.r";
print RSCRIPT "
options(stringsAsFactors=FALSE)
library(vegan)
library(ggplot2)
#library(ggprism)
library(ggpubr)
library(reshape2)
library(plyr)
otu_raw <- read.table(file=\"$input\",sep=\"\\t\",header=T,check.names=FALSE)
colnames(otu_raw) <- c(\"sample\",\"value\",\"id\")
otu_raw\$id <- factor (otu_raw\$id,levels=c($sort_group))
#p <- ggplot(data=otu_raw,aes(x=id,y=value,fill=id))+geom_boxplot()+xlab(\"group\")+labs(fill=\"group\")
#p <-  ggboxplot(data=otu_raw, x=\"id\", y=\"value\", color = \"id\")+xlab(\"group\")+labs(color=\"group\")
p <-  ggboxplot(data=otu_raw, x=\"id\", y=\"value\", fill = \"id\")+xlab(\"\")+labs(color=\"group\")
if (\"$opts{point}\" != \"F\"){
p <- p+geom_jitter()}
p <- p+ylab(\"$label[0]\")
if (\"$test\" == \"None\") {
p2 <- p}else {
p2 <- p+stat_compare_means(method = \"$test\",method.args = list(\"two.sided\"),comparisons = list($mp),label= \"p.signif\")}


p2<-p2+scale_fill_manual(name=\"\",values=c($col),guide=FALSE)
if ($angle == 90) {      
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=0.5,size=$rlc,color=\"black\"))}else if ($angle == 0) {p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=0.5,size=$rlc,color=\"black\"))
}else {
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=1,size=$rlc,color=\"black\"))
}
p2 <- p2+theme( legend.title = element_text(size=$llc),legend.text = element_text(size=$llc),axis.text.y=element_text(size=$clc,color=\"black\"))

ggsave(p2,file = \"$label[0].pdf\",width = $opts{w},height=$opts{h})
png(file = \"$label[0].png\",width = $opts{w}*240,height=$opts{h}*200,res=300)
print(p2)
dev.off()
svg(file = \"$label[0].svg\",width = $opts{w},height=$opts{h})
print(p2)
dev.off()


";
`R --no-save < cmd.r`;




}




if ($opts{m} eq "F") {
	if (-e "alpha.txt") {
		print "alpha have done!!\n";
	}
        else {
	`rm -rf tmp_alpha`;
	&alpha($qiime2_2,$opts{i},$opts{rep});
	&samsum($opts{i},"reads.txt");
	`find $current_dir/tmp_alpha/*/*/data/alpha-diversity.tsv  > alpha.lst`;
	`find $current_dir/pd.txt >> alpha.lst`; 
	`find $current_dir/reads.txt >> alpha.lst`;
	&merge_table("alpha.lst","alpha.txt");
   }
}


else {

if (-e "alpha.txt") {
	&split_file("alpha.txt",$opts{m});
	&Rboxplot("chao1.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("goods_coverage.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("observed_features.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("pd.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("pielou_evenness.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("Reads.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("shannon_entropy.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("simpson.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
	&Rboxplot("ace.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
}
else {
	`rm -rf tmp_alpha`;
&alpha($qiime2_2,$opts{i},$opts{rep});
&samsum($opts{i},"reads.txt"); 
`find $current_dir/tmp_alpha/*/*/data/alpha-diversity.tsv  > alpha.lst`;
`find $current_dir/pd.txt >> alpha.lst`; 
`find $current_dir/reads.txt >> alpha.lst`;
&merge_table("alpha.lst","alpha.txt");
&split_file("alpha.txt",$opts{m});

&Rboxplot("chao1.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("goods_coverage.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("observed_features.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("pd.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("pielou_evenness.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("Reads.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("shannon_entropy.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("simpson.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});
&Rboxplot("ace.tsv",$opts{m},$opts{color},$opts{g},$opts{point},$opts{test},$opts{angle},$opts{rlc},$opts{llc},$opts{clc});


}

}
#die "perl $0 <otu.group.txt> <chao.group.txt>  <simpson.group.txt> <shannon.group.txt> <map.txt> <output.pdf> <sortfile>\n" if @ARGV!=7;
#otu|1,chao|2,PD|3,Shannon|4

=head;
open ORDER,$ARGV[6] ;
chomp (my @order = <ORDER>) ; 
close ORDER; 

my $ord = join ('","',@order) ;
$ord = '"'.$ord.'"';
=cut; 






########delete ns compare#######
sub delete_ns_compare {
         my ($groupfile,$file,$statistictest) = @_;



 #my %group=();
open I,$groupfile ;

<I> ;
my %group=();
while (<I>) {
	chomp;
	my @g = split/\t/,$_;
	$group{$g[1]} = 0;
}
close I;

my @group2 = keys (%group) ; 
#wilcox.test(re$value[re$id=="Gut"],re$value[re$id=="Skin"],alternative = "two.sided")
#
#
my @ok=();
my @ok2=();

foreach my $i (0..$#group2) {
	foreach my $j ($i+1..$#group2) {

                my $tmp=$group2[$i].",".$group2[$j] ; 
                push @ok2,$tmp;
		push @ok,"c(".'"'.$group2[$i].'"'.",".'"'.$group2[$j].'"'.")";
		#print $group2[$i],"==>>$group2[$j]\n";
		#open O,">$group2[$i]\_\_$group2[$j]\.sbgroup.tsv";
		#open II,$ARGV[0] ;
		#my $hh = <II>;
		#print O $hh;
		#while (<II>) {
		#	chomp;
		#	my @gg= split/\t/,$_;
		#	if ($gg[1] eq $group2[$i]) {
		#		print O join("\t",@gg),"\n" ;
		#	}
		#	if ($gg[1] eq $group2[$j]) {
		#		print O join("\t",@gg),"\n" ;
		#	}
		#}
		#close O;
		#close II;
		#print $group2[$i] ,"====>>>",$group2[$j],"\n";
	}
}

#print join(",",@ok),"\n"; 
########delete ns compare#######
#sub delete_ns_compare {
#           my $file = shift; 
           my $tmp = "";
           my @ok3 = (); 
           my $ok4 = ""; 
           open Rspt,">cmd0.r" ; 
           print Rspt "
           re  <- read.table(\"$file\",header=T,sep=\"\\t\",check.names = F)
	   colnames(re) <- c(\"sample\",\"value\",\"id\")
           ";
           my @group = ();
           foreach my $k (@ok2) {
                my @m = split/,/,$k;
                push @group,join("_",@m);
                 my $ttt = $m[0]."_".$m[1];
                  $ttt=~s/-/_/g;
                print Rspt "
           ##re  <- read.table(\"$file\",header=T,sep=\"\\t\",check.names = F)
           tmp  <- \"\"
	   if (\"$statistictest\" == \"t.test\"){
	       tmp <- t.test(re\$value[re\$id==\"$m[0]\"],re\$value[re\$id==\"$m[1]\"],alternative = \"two.sided\")
	       }else if (\"$statistictest\" == \"wilcox.test\") {
	        tmp <- wilcox.test(re\$value[re\$id==\"$m[0]\"],re\$value[re\$id==\"$m[1]\"],alternative = \"two.sided\")
	   }#else {stop(\"Please attention statistic test only can :  T-test or Wilcox-test !!!!!!!!!!!!!!!!!!!\")}
	   $ttt  <- c(\"$m[0]\",\"$m[1]\",tmp\$p.value)
           " ;
           }
           my $tmpp = join(",",@group) ;
                $tmpp=~s/-/_/g;
           print Rspt "
           tt <- as.data.frame(rbind($tmpp)) 
           write.table(tt,file = paste(\"$file\",\".tmp.txt\",sep=\"\"),sep=\"\\t\",quote=F)
           ";
           close Rspt;   
           `/usr/bin/R --no-save<cmd0.r` ;
	   # `rm cmd.r` ;
            open II,$file.".tmp.txt" ;
            <II> ;
            while (<II>) {
                   chomp;
                   my @a = split/\s+/,$_;
                   if($a[3] <= 0.05) {
                              push @ok3,"c(".'"'.$a[1].'"'.",".'"'.$a[2].'"'.")";
                                 }
           }
           close II;
           $ok4 = join(",",@ok3);
           return $ok4;

}
## Wilcox-test,T-test
#my $mp = &delete_ns_compare("group","sub_test.data","T-test") ;
#print "$mp\n";

if ($opts{m} ne  "F") {

my @file =glob ("*.tsv.tmp.txt") ;
open OUT,">alpha-test.txt";
print OUT "#Compare\tGroup1\tGroup2\tpvalue\talpha\n";
foreach my $i (@file) {
	my @split=split/\.tsv/,$i;
	open II,$i;
	<II>;
	while (<II>) {
		chomp;
		print OUT "$_\t$split[0]\n";
	}
	close II;
}
close OUT;
`rm *.tsv.tmp.txt`;
}
`rm *.r *.tsv`;
