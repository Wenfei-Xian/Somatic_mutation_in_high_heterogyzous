use strict;
use warnings;

open my$in1,"$ARGV[0]";
open my$out1,">$ARGV[0].Chr1.ID";
open my$out2,">$ARGV[0].Chr2.ID";
open my$out3,">$ARGV[0].Chr3.ID";
open my$out4,">$ARGV[0].Chr4.ID";
open my$out5,">$ARGV[0].Chr5.ID";
open my$out6,">$ARGV[0].Chr6.ID";
open my$out7,">$ARGV[0].Chr7.ID";
open my$out8,">$ARGV[0].Chr8.ID";
open my$out9,">$ARGV[0].Chr9.ID";
open my$out10,">$ARGV[0].Chr10.ID";
open my$out11,">$ARGV[0].Chr11.ID";
open my$out12,">$ARGV[0].Chr12.ID";
while(<$in1>){
	chomp;
	my$a=(split /\t/,$_)[1];
	if( $a eq "Chr1" ){
		print $out1 "$_\n";
	}
	elsif( $a eq "Chr2" ){
		print $out2 "$_\n";
	}
	elsif( $a eq "Chr3" ){
		print $out3 "$_\n";
	}
	elsif( $a eq "Chr4" ){
                print $out4 "$_\n";
        }
	elsif( $a eq "Chr5" ){
                print $out5 "$_\n";
        }
	elsif( $a eq "Chr6" ){
                print $out6 "$_\n";
        }
	elsif( $a eq "Chr7" ){
                print $out7 "$_\n";
        }
	elsif( $a eq "Chr8" ){
                print $out8 "$_\n";
        }
	elsif( $a eq "Chr9" ){
                print $out9 "$_\n";
        }
	elsif( $a eq "Chr10" ){
                print $out10 "$_\n";
        }
	elsif( $a eq "Chr11" ){
                print $out11 "$_\n";
        }
	elsif( $a eq "Chr12" ){
                print $out12 "$_\n";
        }
}
