use warnings;
use strict;

open IN1,"$ARGV[0]";
my%hash;
while(<IN1>){
	chomp;
	my@blast=(split /\t/,$_);
	next if( $blast[10] < 4050);
	if( exists $hash{$blast[0]} ){
	}
	else{
		$hash{$blast[0]}=$blast[0];
		my($read,$pair)=(split /_/,$blast[0])[0,1];
		print "$read\t$pair\t$blast[5]\t$blast[7]\t$blast[8]\n";
		#print "$_\n";
	}
}
