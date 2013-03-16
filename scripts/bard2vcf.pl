#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open :use_bgzip);
use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr( q$Rev: 0$, 4 ) / 10000;

=head1 NAME

bard2vcf.pl - script to convert bardCNV calls to VCFv4.1/TCGA1.2 format.

=head1 SYNOPSIS

  bard2vcf.pl --obs <bard observations file> --states <bard states file> --ref <reference fasta> --sample1 <sample_name> --sample2 <sample_name> [-notabix]

=head1 DESCRIPTION

This script reads in bardCNV output and writes a VCFv4.1/TCGA1.2 file containing the state calls when they are not diploid.  If output files end in '.gz' it will compress with bgzip and create tabix indexes. To prevent index creation pass -notabix option

=cut

#------------
# Begin MAIN
#------------

process_commandline();

my $bard_file = $Opt{'states'};
my $obs_file = $Opt{'obs'};
my $ref_file = $Opt{'ref'};
my $sample1_name = $Opt{'sample1'};
my $sample2_name = $Opt{'sample2'};
my $outfile = $Opt{'outfile'};

my $outfh = Open($outfile, 'w');

my $ref_length_string = `awk '{print \$1, \$2}' $ref_file.fai`;
my $rh_ref_length;
my @ref_entries = ();
foreach my $ref_line (split /\n/, $ref_length_string) {
    if ($ref_line =~ /^(\S+)\s(\d+)/) {
        $rh_ref_length->{$1} = $2;
        push @ref_entries, $1;
    }
    else {
        die "Misformatted string $ref_line in file $ref_file.fai!\n";
    }
}

my $mode = $Opt{'append'} ? 'a' : 'w';

print_header( $outfh, $sample1_name, $sample2_name ) unless ( $Opt{'noheader'} );
my ($rh_cnv_calls, $mean_rdrc) = read_cnv_states($bard_file);
read_observations($obs_file, $rh_cnv_calls, $mean_rdrc);
write_vcf_lines( $outfh, $rh_cnv_calls, $rh_ref_length, \@ref_entries);
undef $outfh;

if ( $outfile =~ /gz$/ && $Opt{tabix} ) {
    unlink "$outfile.tbi";  # tabix won't overwrite, so delete first
    my $cmd = " $Opt{tabix} -p vcf $outfile ";
    my $w = `$cmd`;
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( sample1 => 'Unknown', sample2 => 'Unknown', min_score => 60, max_pad => 50000 );
    GetOptions(
        \%Opt, qw(
            manual help+ version
            verbose states=s obs=s sample1=s sample2=s outfile=s ref=s
            append gzip noheader tabix! min_score=i min_somatic_score=i max_pad=i
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "mybard2vcf.pl, ", q$Revision: 4129 $, "\n";
	}

	if ( !exists( $Opt{states} ) ) {
		die "Specify bardCNV output filename with option --states.\n";
	}

	if ( !exists( $Opt{obs} ) ) {
		die "Specify bardCNV observation filename with option --obs.\n";
	}

	if ( !exists( $Opt{ref} ) ) {
		die "Specify reference fasta filename with option --ref.\n";
	}

	if ( !$Opt{outfile} ) {
		$Opt{outfile} = "$Opt{sample2}.mpv.vcf";
		if ($Opt{gzip}) {
			$Opt{outfile} .= ".gz";
		}
	}

	if ( !defined $Opt{tabix} ) {
		$Opt{tabix} = `which tabix`;
		chomp $Opt{tabix};

                if (!(defined $Opt{tabix})) {
                    print STDERR "Can\'t find tabix executable--will not be able to index VCF file!\n";
                }
	}
}

sub print_header {
        my $outfh = shift;
	my $sample1_name = shift;
	my $sample2_name = shift;

	print $outfh "##fileformat=VCFv4.1\n";
	print $outfh "##tcgaversion=1.2\n";
	my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime();
	$mon++;
	$year += 1900;
	printf $outfh "##fileDate=%d%02d%02d\n", $year, $mon, $mday;
	print $outfh "##center=\"NHGRI\"\n";
	print $outfh "##reference=GRCh37\n";
        print $outfh "##vcfProcessLog=<InputVCFSource=<bardCNV,bard2vcf.pl>,InputVCFVer=<0.1,0.1>>\n";

        print $outfh "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">\n";
        print $outfh "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";

	# included genotype id's:
	print $outfh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print $outfh "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        print $outfh "##FORMAT=<ID=MCN,Number=1,Type=Integer,Description=\"Minor allele copy number for imprecise events\">\n";
        print $outfh "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n";
        print $outfh "##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">\n";
	print $outfh "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position in the sample\">\n";
        print $outfh "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Depth of reads supporting alleles 0/1/2/3...\">\n"; 
        print $outfh "##FORMAT=<ID=BQ,Number=.,Type=Integer,Description=\"Average base quality for reads supporting alleles\">\n";
        print $outfh "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">\n";
        print $outfh "##FORMAT=<ID=SSC,Number=1,Type=Integer,Description=\"Somatic Score\">\n";

        print $outfh "##SAMPLE=<ID=$sample2_name>\n";
        print $outfh "##SAMPLE=<ID=$sample1_name>\n";

	# print required eight fields:
	print $outfh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample2_name\t$sample1_name\n";
}

sub read_cnv_states {
    my $state_file = shift;

    my $fh = Open($state_file);

    my $rh_state_call = {};
    my $rdr_sum = 0;
    my $rdr_total = 0;
    while (<$fh>) {
        if (/^(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($chr, $pos, $rdr, $alt_freq, $minor, $total, $score) =
                 ($1, $2, $3, $4, $5, $6, $7);
            my $somatic_score = ($minor==1 && $total==2) ? 0 : 
                      $score + 5*abs($total + $minor - 3.0);
            $rh_state_call->{$chr}->{$pos} = {'minor' => $minor,
                                              'total' => $total,
                                              'score' => $score,
                                              'somatic_score' => $somatic_score}; 
            if ($score == 90 && $total != 0) {
                $rdr_sum += $rdr/$total;
                $rdr_total++;
            }
        }
    }
    my $mean_rdrc = $rdr_sum/$rdr_total;
    return ($rh_state_call, $mean_rdrc);
}

sub read_observations {
    my $obs_file = shift;
    my $rh_states = shift;
    my $mean_rdrc = shift;

    my $fh = Open($obs_file);

    while (<$fh>) {
        if (/^(\S+)\s(\d+)\s(\S+)\s(\S+)$/) {
            my ($chr, $pos, $rdr, $alt_freq) = ($1, $2, $3, $4);
            if (!$rh_states->{$chr}->{$pos}) {
                $alt_freq = 1 - $alt_freq if ($alt_freq > 0.5);
                my $minor_est = int($alt_freq*$rdr/$mean_rdrc);
                $rh_states->{$chr}->{$pos} = {'minor' => $minor_est,
                                              'total' => int($rdr/$mean_rdrc),
                                              'score' => 99, 'somatic_score' => 90};
            }
        }
    }
}

sub write_vcf_lines {
    my $outfh = shift;
    my $rh_states = shift;
    my $rh_ref_length = shift;
    my $ra_ref_entries = shift;

    foreach my $chr (@ref_entries) {
        my @vcf_lines = ();
        my $soft_start = 1;
        my ($hard_start, $hard_end, $soft_end);
        my $score_sum = 0;
        my $somatic_score_sum = 0;
        my $state_positions = 0;
        my ($current_state_total, $current_state_minor);
        my $high_copy_flag;
        my $high_copy_sum = 0;
        my $high_copy_positions = 0;
        my $last_line; # keep track of last VCF line so we can fill non-high-quality space
        foreach my $pos (sort {$a <=> $b} keys %{$rh_states->{$chr}}) {
            my $state_total = $rh_states->{$chr}->{$pos}->{'total'}; 
            my $state_minor = defined ($rh_states->{$chr}->{$pos}->{'minor'}) ? $rh_states->{$chr}->{$pos}->{'minor'} : '.'; 
            my $score = $rh_states->{$chr}->{$pos}->{'score'}; 
            my $somatic_score = $rh_states->{$chr}->{$pos}->{'somatic_score'}; 
            next if ($somatic_score < $Opt{'min_score'});
            if (!defined($current_state_total) ) {
                $hard_start = $pos;
                $hard_end = $pos;
                $current_state_total = $state_total;
                $current_state_minor = $state_minor;
                $score_sum = $score;
                $somatic_score_sum = $somatic_score;
                $state_positions++;
                my $pad_start = ($hard_start > $Opt{'max_pad'}) ? $hard_start - $Opt{'max_pad'} : 1; 
                push @vcf_lines, vcf_line($chr, 0, $pad_start, $hard_start-1, 0, $current_state_total, $current_state_minor, 10, 10);
            }
            elsif (($current_state_total == $state_total && $current_state_minor eq $state_minor) || ($high_copy_flag && $score == 99)) { # extend hard boundaries
                $hard_end = $pos;
                $score_sum += $score;
                $somatic_score_sum += $somatic_score;
                $state_positions++;
                if ($high_copy_flag) {
                    $high_copy_sum += $state_total;
                    $high_copy_positions++;
                    $current_state_total = $state_total;
                    $current_state_minor = $state_minor;
                }
            }
            elsif (($current_state_total != $state_total) || ($current_state_minor != $state_minor)) { # end of last state, beginning of new!
                $soft_end = $pos;
                my $score_mean = ($state_positions) ? int($score_sum/$state_positions) : 0;
                my $somatic_score_mean = ($state_positions) ? int($somatic_score_sum/$state_positions) : 0;
                if ($high_copy_flag) {
                    $current_state_total = int($high_copy_sum/$high_copy_positions);
                    $current_state_minor = '.';
                    $score_mean = 10;
                    $somatic_score_mean = 90;
                }
                my $next_line = vcf_line($chr, $soft_start, $hard_start, $hard_end, $soft_end, $current_state_total, $current_state_minor, $score_mean, $somatic_score_mean);
                # do we need to fill space with low quality calls?
                if ($last_line) {
                    push @vcf_lines, low_qual_filler($last_line, $next_line);
                } 
                push @vcf_lines, $next_line; 
                $last_line = $next_line;
                # reset state: 
                $current_state_total = $state_total;
                $current_state_minor = $state_minor;
                $score_sum = $score;
                $somatic_score_sum = $somatic_score;
                $state_positions = 1;
                $soft_start = $hard_end;
                $hard_start = $pos;
                $hard_end = $pos;
                if ($score == 99) { # high copy number state
                    $high_copy_flag = 1;
                    $high_copy_sum += $state_total;
                    $high_copy_positions = 1;
                }
                else {
                    $high_copy_flag = 0;
                    $high_copy_sum = 0;
                    $high_copy_positions = 0;
                }
            }
        }
        # print last state for this chromosome:
        if ($state_positions) {
            $soft_end = $rh_ref_length->{$chr};
            my $score_mean = ($state_positions) ? int($score_sum/$state_positions) : 0;
            my $somatic_score_mean = ($state_positions) ? int($somatic_score_sum/$state_positions) : 0;
            my $next_line = vcf_line($chr, $soft_start, $hard_start, $hard_end, $soft_end, $current_state_total, $current_state_minor, $score_mean, $somatic_score_mean);
            if ($last_line) {
                push @vcf_lines, low_qual_filler($last_line, $next_line);
            }
            push @vcf_lines, $next_line;
            my $pad_end = ($soft_end > $hard_end + $Opt{'max_pad'}) ? $hard_end + $Opt{'max_pad'} : $soft_end;
            push @vcf_lines, vcf_line($chr, 0, $hard_end+1, $pad_end, 0, $current_state_total, $current_state_minor, 10, 10);
            $last_line = $next_line;
        }
        while (my $line = shift @vcf_lines) {
            my $cn_state = ($line =~ /GT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t(\S+)/) ? $1 : 'NA';
            my $old_end = ($line =~ /END=(\d+)/) ? $1 : 'NA';
            while ($vcf_lines[0] && 
                  $vcf_lines[0] =~ /END=(\d+).*GT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t$cn_state\s/ && 
                  $1<$old_end+$Opt{'max_pad'}) {
                my $add_line = shift @vcf_lines;
                my $new_end = ($add_line =~ /END=(\d+)/) ? $1 : 'NA';
                $line =~ s/END=\d+/END=$new_end/;
            }
            print $outfh "$line";
        }
    }
}

sub vcf_line {
    my $chr = shift;
    my $soft_start = shift;
    my $hard_start = shift;
    my $hard_end = shift;
    my $soft_end = shift;
    my $current_state_total = shift;
    my $current_state_minor = shift;
    my $score_mean = shift;
    my $somatic_score_mean = shift;

    return "$chr\t$hard_start\t.\t.\t<CNV>\t$somatic_score_mean\tPASS\tSOMATIC;END=$hard_end\tGT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t.:.:.:.:.:2:$current_state_minor:$current_state_total:$score_mean:$somatic_score_mean\t.:.:.:.:.:.:1:2:40:.\n";

}

sub low_qual_filler {
    my $last_line = shift;
    my $next_line = shift;

    my @lq_vcf_lines = ();
    my ($chr, $last_start, $last_end, $last_minor, $last_cn, $last_cnq, $next_start, $next_end, $next_minor, $next_cn, $next_cnq);
    if ($last_line =~ /^(\S+)\t(\d+).*END=(\d+)\tGT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t.:.:.:.:.:.:([^:]+):(\d+):(\d+):(\d+)/) {
        ($chr, $last_start, $last_end, $last_minor, $last_cn, $last_cnq) = ($1, $2, $3, $4, $5, $6);
    }
    else {
        die "Unrecognized format in VCF line $last_line\n";
    }
    if ($next_line =~ /^$chr\t(\d+).*END=(\d+)\tGT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t.:.:.:.:.:.:([^:]+):(\d+):(\d+):(\d+)/) {
        ($next_start, $next_end, $next_minor, $next_cn, $next_cnq) = ($1, $2, $3, $4, $5);
    }
    else {
        die "Unrecognized format in VCF line $next_line\n";
    }

    my $gap_length = $next_start - $last_end - 1;
    my $left_gap_portion = ($last_end - $last_start + 1)/($next_end - $next_start + $last_end - $last_start + 2);
    my $left_length = int($gap_length*$left_gap_portion);
    my $right_length = $gap_length - $left_length;
    my $left_start = $last_end + 1;
    my $left_end = $left_start + $left_length - 1;
    $left_end = $left_start + $Opt{'max_pad'} if ($left_end > $left_start + $Opt{'max_pad'});
    my $right_start = $left_end + 1;
    my $right_end = $right_start + $right_length - 1;
    $right_start = $right_end - $Opt{'max_pad'} if ($right_start < $right_end - $Opt{'max_pad'});
    if ($left_end - $left_start > 1) {
        push @lq_vcf_lines, "$chr\t$left_start\t.\t.\t<CNV>\t10\tPASS\tSOMATIC;END=$left_end\tGT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t.:.:.:.:.:2:$last_minor:$last_cn:10:10\t.:.:.:.:.:.:1:2:40:.\n";
    }
    if ($right_end - $right_start > 1) {
        push @lq_vcf_lines, "$chr\t$right_start\t.\t.\t<CNV>\t10\tPASS\tSOMATIC;END=$right_end\tGT:GQ:DP:AD:BQ:SS:MCN:CN:CNQ:SSC\t.:.:.:.:.:2:$next_minor:$next_cn:10:10\t.:.:.:.:.:.:1:2:40:.\n";
    }

    return @lq_vcf_lines;
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--sample>

Specify the sample name from which these genotypes come.  This will be used
inside the VCF file as the column header, and as the base name for the
default output file names.

=item B<--notabix>

By default when a output file ends in 'gz', tabix will create the corresponding index.
To override this behaviour specify --notabix.

=item B<--min_score> N

Minimum bardCNV score to record a "hard" position for a CNV (default 60).

=item B<--gzip>

Generate compressed output files (suffixed with .gz extension).  Has no
effect when C<--snv_outfile> or C<--div_outfile> are explicitly specified.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use.

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose.

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation.

=cut
