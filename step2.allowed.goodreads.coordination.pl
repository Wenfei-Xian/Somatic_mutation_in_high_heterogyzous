use strict;
use warnings;

sub binary_search {
    my ($array, $target) = @_;
    my ($low, $high) = (0, scalar(@$array) - 1);

    while ($low <= $high) {
        my $mid = int(($low + $high) / 2);
        my ($start, $end) = split /\t/, $array->[$mid];

        if ($target < $start) {
            $high = $mid - 1;
        } elsif ($target > $end) {
            $low = $mid + 1;
        } else {
            return $mid;  # 找到目标区间，返回索引
        }
    }

    return -1;  # 没有找到目标区间
}

my %hash;
my %hash_pos_hap;

open IN1, "$ARGV[0]" or die "Cannot open file $ARGV[0]: $!";  # output from step1
<IN1>;
while (<IN1>) {
    chomp;
    my @info = split /\t/;
    push @{$hash{$info[8]}}, "$info[9]\t$info[10]";
    $hash_pos_hap{"$info[8]\t$info[9]\t$info[10]"} = "$info[4]\t$info[5]\t$info[6]";
}
close IN1;

foreach my $chr (keys %hash) {
    @{$hash{$chr}} = sort {
        my ($a_start, $a_end) = split /\t/, $a;
        my ($b_start, $b_end) = split /\t/, $b;
        $a_start <=> $b_start;
    } @{$hash{$chr}};
}

open IN2, "samtools view $ARGV[1] $ARGV[2] |" or die "Cannot open pipe: $!";  # bam file
while (<IN2>) {
    chomp;
    my @sam = split /\t/;
    my ($read, $flag, $chr, $pos) = ($1, $2, $3, $4) if (/(\S+)\t(\d+)\t(\S+)\t(\d+)/);

    if (exists $hash{$chr}) {
        my $index = binary_search($hash{$chr}, $pos);
        if ($index != -1) {
            my ($start, $end) = split /\t/, $hash{$chr}[$index];
            if ($pos >= $start && $pos <= $end) {
                my $coor = "$chr\t$start\t$end";
		print "$read\t$hash_pos_hap{$coor}\n";
		#print "$read\t$hash_pos_hap{$coor}\t$_\n";
            }
        }
    }
}
close IN2;
