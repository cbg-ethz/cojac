> **Notes:**
> - This directory contains Jupyter and Rstudio notebooks used in the [publication](#citation).
> - If you are interested in analysing wastewater samples, you should refer to repository [COWWID](https://github.com/cbg-ethz/cowwid).
>
>   It contains a detailed description of the latest up-to-date procedure used to generate the plots visible on
>   the [dashboard _Surveillance of SARS-CoV-2 genomic variants in wastewater_](https://bsse.ethz.ch/cbg/research/computational-virology/sarscov2-variants-wastewater-surveillance.html)
>   and displayed [on CoV-Spectrum](https://cov-spectrum.ethz.ch/story/wastewater-in-switzerland).


## Signature_mutations_in_patient_samples.ipynb

Notebook to assess the prevalence of signature mutation defining the lineages B.1.1.7 and 501.V2 in all 5758 non-B.1.1.7 and non-501.V2 consensus sequences from clinical samples collected in Switzerland before December 24.


## snv_count_wastewater3.ipynb

ShoRAH normally calls SNV by combining the results of three overlapping local windows.
Given the nature of our sample, we are lucky if a single windows can cover our mutations.
So this Notebook instead calls mutation based on the local haplotype on single ShoRAH windows.
We lose some confidence by only using a single windows per SVN, but because _variants of concern_ are detected based on multiple mutations simultaneously.


## mut-table.ipynb

Combines a list of SNV from the precedent Notebook with a reads counts, so we can also estimate mutation frequency.
