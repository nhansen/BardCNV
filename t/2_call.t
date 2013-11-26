# $Id$
# t/02_call.t - check that bardCNV successfully calls copy number alterations on a pair of bam files 
use strict;
use Test::More;
use Module::Build;

my $script = 'scripts/bardCNV.pl';

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => 4;

$ENV{PATH} = "./c:./scripts:$ENV{PATH}";

my $out;
my $testref = 't/ref.fasta';
my $testhet = 't/hetfile.bed';
my $normalbam = 't/normal_sample.bam';
my $tumorbam = 't/tumor_sample.bam';
system("perl -w -I lib $script --hetfile $testhet --ref $testref --outdir t/test1out --normalbam $normalbam --tumorbam $tumorbam > t/calltest1.out 2>&1");
$out = `wc -l t/test1out/states.out`;
like $out, qr/^447\s/, "$script states file";
$out = `grep 'mu=' t/test1out/best_model_file.txt`;
like $out, qr/\s0\.418411/, "$script mu value";
system("perl -w -I lib $script --hetfile $testhet --ref $testref --outdir t/test1out --normalbam $normalbam --tumorbam $tumorbam --diploid --skipcounts > t/calltest2.out 2>&1");
$out = `awk '\$6==1 && \$7==3 {print}' t/test1out/states.out | wc -l`;
like $out, qr/^446\s/, "$script diploid states file";
$out = `grep 'mu=' t/test1out/diploid_model_file.txt`;
like $out, qr/\s0\.455691/, "$script diploid mu value";
