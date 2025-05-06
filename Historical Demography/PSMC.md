# Using PSMC to infer historical demography 

    fq2psmcfa -q20 SRR25478317_eseal_sorted_subset.fq.gz > SRR25478317_eseal_sorted_subset.psmcfa

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o SRR25478317_eseal_sorted_subset.psmc SRR25478317_eseal_sorted_subset.psmcfa



    psmc_plot.pl -g 10 -u 1e-8 -X 10000000 -Y 16 SRR25478317_eseal_test SRR25478317_eseal_sorted_subset.psmc \ epstopdf SRR25478317_eseal_test.eps


    splitfa SRR25478317_eseal_sorted_subset.psmcfa > SRR25478317_eseal_sorted_subset_split_psmcfa_file.split.psmcfa


    seq 100 | xargs -P 4 -I {} psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o SRR25478317_eseal_sorted_subset_round-{}.psmc SRR25478317_eseal_sorted_subset_split_psmcfa_file.split.psmcfa

