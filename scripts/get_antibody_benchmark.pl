#!/usr/bin/perl

# get a list of all antibody complex structures
# search for antigen chain sequence(s), find unbound ones 
# search for antibody chain sequences, find unbound ones

use strict;

my $test_pdb = $ARGV[0];

my $PERCID_CUTOFF = 0.98;
my $COV_CUTOFF = 0.8;
my $RES_CUTOFF = 3.25;

my $data_dir = "/home/pierceb/databases/";
my $blast_exe = "/ibbr/ncbi-blast/bin/blastp";

my $sabdab_list = $data_dir . "all_antibody_structures_8_2020.txt";
my $seqres_file = $data_dir . "pdb_seqres_8_2020.fa";
my $pdb_res_file = $data_dir . "pdb_resolutions_8_2020.txt";

open(RES, $pdb_res_file) || die("unable to open file: $pdb_res_file\n");
my @res_lines = <RES>;
close(RES);

# parse those PDB resolutions
my %pdbres = ();
my $no_res_yet = 1;
for (my $i = 0; $i < @res_lines; $i++)
{
    if (substr($res_lines[$i], 0, 6) eq "IDCODE") { $i += 2; $no_res_yet = 0; }
    if ($no_res_yet == 1) { next; }
    (my $pdb, my $stuff, my $res) = split(" ", $res_lines[$i]);
    $pdbres{lc($pdb)} = $res;
}

open(SD, $sabdab_list) || die("unable to open file: $sabdab_list\n");
my @list_lines = <SD>;
close(SD);

open(SR, $seqres_file) || die("unable to open file: $seqres_file\n");
my @seqres_lines = <SR>;
close(SR);
chomp(@seqres_lines);

my %pdb_processed = ();
foreach my $line (@list_lines)
{
    my @fields = split("\t", $line);
    my $pdb = $fields[0];
    if ($pdb_processed{$pdb} ne "") { next; }
	$pdb_processed{$pdb} = 1;
    my $hchain = $fields[1];
    my $lchain = $fields[2];
    my $model = $fields[3];
    my $ag_chain = $fields[4];
    my $ag_type = $fields[5];
    if ($model == 1) { next; }
    my $res = $fields[16];
    my $method = $fields[17];
    if (($res > $RES_CUTOFF) || (substr($method, 0, 5) ne "X-RAY")) { next; }
    if (($pdbres{$pdb} > $RES_CUTOFF) || ($pdbres{$pdb} <= 0)) { next; } # this should not be necessary
    if (($hchain eq "NA") || ($lchain eq "NA") || ($ag_chain eq "NA")) { next; }
    if (($ag_type ne "protein") && (substr($ag_type, 0, 7) ne "protein")) { next; }
    if (($test_pdb ne "") && ($pdb ne lc($test_pdb))) { next; }
	# if ($hchain eq $lchain) { print "$pdb\n"; } # scFv?
    # now proceed to get the chain sequences
    my @ag_chains = split('\|', $ag_chain);
    for (my $i = 0; $i < @ag_chains; $i++) { $ag_chains[$i] =~ s/^\s+|\s+$//g; }
    # get the sequences
    my @seqs = ();
    my @chains = ($hchain, $lchain, @ag_chains);
    #for (my $i = 0; $i < @chains; $i++) { print $chains[$i]; }
    # get all the chain sequences
    my @chain_matches = ();
    for (my $i = 0; $i < @chains; $i++)
    {
	$chain_matches[$i] = "";
	my $chn = $chains[$i];
	my $seq = "";
	for (my $j = 0; $j < @seqres_lines; $j++)
	{
	    $line = $seqres_lines[$j];
	    if (substr($line, 0, 1) ne ">") { next; }
	    my $pdb2 = substr($line, 1, 4);
	    if ($pdb2 ne $pdb) { next; } 
	    my $chain = substr($line, 6, 1);
	    if ($chain ne $chn) { next; }
	    $j++;
	    while (($j < @seqres_lines) && (substr($seqres_lines[$j], 0, 1) ne ">")) { $seq .= $seqres_lines[$j]; $j++; }
	    last;
	}
	if ($seq eq "") { print "error finding sequence for $pdb $chn\n"; last; }
	open(TMP, "> temp.fa") || die("unable to open temp fa file\n");
	print TMP "> $pdb $chn\n";
	print TMP "$seq\n";

	# run BLAST
	my $blast_cmd = $blast_exe . " -query temp.fa -db $seqres_file";
	my $blout = `$blast_cmd`;
	
	my @blastout_lines = split("\n", $blout);
	my $tot_num = 0;
	my $match_num = 0;
	my $query_offset = -1;
	my $subj_offset = -1;
	my $curr_chn = "";
	my $curr_pdb = "";
	for (my $j = 0; $j < @blastout_lines; $j++)
	{
	    my $line = $blastout_lines[$j];
	    if ($line =~ /^> (\w\w\w\w)_(\w) mol:protein/) { $curr_pdb = $1; $curr_chn = $2; }
	    elsif ($line  =~ /^>(\w\w\w\w)_(\w) mol:protein/) { $curr_pdb = $1; $curr_chn = $2; } # other format
	    elsif ($line =~ /Identities = (\d+)\/(\d+) /)
	    {
		$match_num = $1;
		$tot_num = $2;
		my $perc_id = $match_num/$tot_num;
		if ($perc_id < $PERCID_CUTOFF) { next; }
		my $cov_perc = $tot_num/length($seq);
		if ($cov_perc < $COV_CUTOFF) { next; }
		if ($curr_pdb eq $pdb) { next; }
	
		# looks like we may have a winner
		if ($chain_matches[$i] ne "") { $chain_matches[$i] .= "\t"; }
		$chain_matches[$i] .= $curr_pdb . $curr_chn;
	    }
	}
    }
    for (my $i = 0; $i < @chains; $i++) { if ($chain_matches[$i] eq "") { next; }  } # at least one chain does not match

    # now let's see if we have any non-Ab or non-Ag chains in the PDBs
    # antibody
    my @hchain_matches = split("\t", $chain_matches[0]);
    my @lchain_matches = split("\t", $chain_matches[1]);
    my $ub_mab_matches = "";
    my %mab_checked = ();
    foreach my $hchns (@hchain_matches)
    {
	my $hpdb = substr($hchns, 0, 4);
	my $hchn = substr($hchns, 4, 1);
	if ($mab_checked{$hpdb} ne "") { next; }
	$mab_checked{$hpdb} = 1;
	if (($pdbres{$hpdb} > $RES_CUTOFF) || ($pdbres{$hpdb} <= 0)) { next; }
	my $hchains = "";
	my $lchains = "";

	foreach my $lchns (@lchain_matches)
	{
	    my $lpdb = substr($lchns, 0, 4);
	    my $lchn = substr($lchns, 4, 1);
	    if ($lpdb eq $hpdb) { $lchains .= $lchn; }
	}
	if ($lchains eq "") { next; }
	foreach my $hchns2 (@hchain_matches)
	{
	    my $hpdb2 = substr($hchns2, 0, 4);
	    my $hchn2 = substr($hchns2, 4, 1);
	    if ($hpdb2 eq $hpdb) { $hchains .= $hchn2; }
	}
	my $mabchains = $hchains . $lchains;
	my $num_matches = 0;
	my $mismatch_found = 0;

	# do all of the chains in the PDB correspond to antibody chains
	for (my $j = 0; $j < @seqres_lines; $j++)
	{
	    $line = $seqres_lines[$j];
	    if (substr($line, 0, 1) ne ">") { next; }
	    my $pdb2 = substr($line, 1, 4);
	    if ($pdb2 ne $hpdb) { next; } 
	    my $chain = substr($line, 6, 1);
	    $num_matches++;
	    if ($mabchains =~ /$chain/) { } # we have a match }
	    else { $mismatch_found = 1; last; }
	}
	if ($mismatch_found == 1) { next; }
	if ($num_matches != length($mabchains)) { print "error weird chain matches\n"; }
	if ($ub_mab_matches ne "") { $ub_mab_matches .= " "; }
	$ub_mab_matches .= $hpdb;
    }

    # antigen
    my $ub_ag_matches = "";
    my %ag_checked = ();
    my @tmp_matches = split("\t", $chain_matches[2]);
    foreach my $agchns (@tmp_matches)
     {
	my $agpdb = substr($agchns, 0, 4);
	my $agchn = substr($agchns, 4, 1);
	if ($ag_checked{$agpdb} ne "") { next; }
	$ag_checked{$agpdb} = 1;
	if (($pdbres{$agpdb} > $RES_CUTOFF) || ($pdbres{$agpdb} <= 0)) { next; }
	my $allagchns = "";
	for (my $j = 2; $j < @chain_matches; $j++)
	{
	    my $curr_agchns = "";
	    my @tmp_matches = split("\t", $chain_matches[$j]);
	    foreach my $agchns2 (@tmp_matches)
	    {
		my $agpdb2 = substr($agchns2, 0, 4);
		my $agchn2 = substr($agchns2, 4, 1);
		if ($agpdb2 eq $agpdb) { $curr_agchns .= $agchn2; }
	    }
	    if ($curr_agchns eq "") { $allagchns = ""; last; }
	    $allagchns .= $curr_agchns;
	}
	if ($allagchns eq "") { next; }
	my $num_matches = 0;
	my $mismatch_found = 0;
	for (my $j = 0; $j < @seqres_lines; $j++)
	{
	    $line = $seqres_lines[$j];
	    if (substr($line, 0, 1) ne ">") { next; }
	    my $pdb2 = substr($line, 1, 4);
	    if ($pdb2 ne $agpdb) { next; } 
	    my $chain = substr($line, 6, 1);
	    $num_matches++;
	    if ($allagchns =~ /$chain/) { } # we have a match }
	    else { $mismatch_found = 1; last; }
	}
	if ($mismatch_found == 1) { next; }
	if ($num_matches != length($allagchns)) { print "This is weird\n"; }
	if ($ub_ag_matches ne "") { $ub_ag_matches .= " "; }
	$ub_ag_matches .= $agpdb;
    }
    
    print "$pdb\t$ub_mab_matches\t$ub_ag_matches\n";
}


