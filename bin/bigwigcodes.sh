samtools merge -@8 -o MAQCA_merged.bam MAQCA_1_dedup.bam MAQCA_2_dedup.bam

bamCoverage -b Merged_S7_S12.unique.bam -o Merged_S7_S12.Bin1.bw \
  --normalizeUsing CPM \
  --binSize 1 --ignoreDuplicates

computeMatrix scale-regions \
  -S Merged_S7_S12.Bin1.bw \
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_Merged_S7_S12_Bin1.gz \
  --outFileNameMatrix matrix_Merged_S7_S12_Bin1.tab \
  --smartLabels

plotProfile -m matrix_Merged_S7_S12_Bin1.gz -out geneBody_Merged_S7_S12_Bin1.png

bamCoverage -b MAQCA_merged.bam -o MAQCA_merged.bam.bw \
  --normalizeUsing CPM \
  --binSize 10 --ignoreDuplicates
# Strand-specific libraries? Add: --filterRNAstrand forward|reverse

bamCoverage -b MAQCA_2_dedup.bam -o MAQCA_2_dedup.bam.bw \
  --normalizeUsing CPM \
  --binSize 10 --ignoreDuplicates

computeMatrix reference-point \
  -S Merged_S7_S12.unique.bam.bw\
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  --metagene \
  --nanAfterEnd \
  --skipZeros \
  -o matrix_TSS_nanAfterEnd.gz \
  --outFileNameMatrix matrix_TSS_nanAfterEnd.tsv \
  --outFileSortedRegions regions_TSS_nanAfterEnd.bed

plotProfile -m matrix_MergedS7S12_nanAfterEnd.gz -out geneBody_MergedS7S12_nanAfterEnd.png

computeMatrix scale-regions \
  -S MAQCA_2_dedup.bam.bw \
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_MAQC2_QSP.gz \
  --outFileNameMatrix matrix_MAQC2_QSP.tab \
  --smartLabels



bamCoverage -b Microsample_1_50uL_20250409_Aligned.sortedByCoord.out.bam -o Microsample_1.bam.bw \
  --normalizeUsing CPM \
  --binSize 10 --ignoreDuplicates

computeMatrix scale-regions \
  -S Microsample_1.bam.bw \
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_QSP.gz \
  --outFileNameMatrix matrix_QSP.tab \
  --smartLabels

plotProfile -m matrix_QSP.gz -out geneBody_QSP.png


igvtools count ./testdata/Merged_S7_S12.unique.bam Merged_S7_S12.unique.bam.tdf $VSC_DATA_VO/PPOL/resources/repos/IGVTools/genomes/hg38.chrom.sizes
igvtools count MAQCA_1_dedup.bam MAQCA_1_dedup.bam.tdf $VSC_DATA_VO/PPOL/resources/repos/IGVTools/genomes/hg38.chrom.sizes
macs2 callpeak -t /user/gent/446/vsc44685/DataVO_dir/QSP_MAQCA/Samples/star/RNA021484_1Aligned.sortedByCoord.out_dedup.bam --outdir . -n CalledPeaks2_RNA021484_1 -g hs --bdg --keep-dup auto
macs2 callpeak -t MAQCA_merged.bam --outdir . -n CalledPeaks_MAQCA1 -g hs --bdg --keep-dup auto

macs2 callpeak \
  -t /user/gent/446/vsc44685/DataVO_dir/QSP_MAQCA/Samples/star/RNA021484_1Aligned.sortedByCoord.out_dedup.bam \
  -g hs -n RNA021484_1_3prime \
  --nomodel --extsize 75 --shift 0 \
  --call-summits -q 0.01 \
  --keep-dup auto \
  -B --SPMR


bedGraphToBigWig treat_pileup.sorted.bdg $VSC_DATA_VO/PPOL/resources/repos/IGVTools/genomes/hg38.chrom.sizes treat_pileup.bw


{"context":["the current pipeline with the option --makeplots generates coverage plots highlighting the selected peak","peaks are called using macs2 in the process MACS2_CALLPEAK"],
"new_feature":"add a new plot type that displays all called peaks within the gene span",
"details":["called peaks are written in the .narrowPeak file","peaks can be showed as horizontal bars in the bottom pannel of the plot","choose the best approach"]}


{
"version": "1.0",
"objective": "Analyze the provided R script that plots coverage and peaks and verify (or refute) the eight specific problem areas previously identified. Produce a concise, evidence-based report with exact code locations, why each issue matters, and concrete fixes.",
"inputs": {
"r\_script": "MakePlots_new.R",
"context": {
"expected\_behavior": "Top panel shows coverage; bottom panel shows exons, a primer/selected-peak lane, and all peaks from a MACS2 .narrowPeak file as horizontal bars spanning start–end."
}
},
"deliverables": {
"report": {
"format": "json",
"schema": {
"summary": "string",
"checks": \[
{
"id": "string",
"title": "string",
"status": "PASS|FAIL|WARN|NOT\_APPLICABLE",
"evidence": "string",
"code\_locations": \[
{
"excerpt": "string",
"line\_start": "integer",
"line\_end": "integer"
}
],
"why\_it\_matters": "string",
"fix": "string",
"severity": "critical|high|medium|low",
"confidence": "high|medium|low"
}
],
"optional\_improvements": \[
{
"title": "string",
"rationale": "string",
"suggested\_change": "string"
}
]
}
}
},
"checks": \[
{
"id": "C1\_fallback\_combination\_missing",
"title": "Fallback path saves only coverage panel (features/peaks never saved)",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"if \\(inherits\\(plt, "list"\\)\\)\s\*\\{\[\s\S]*ggsave\\(\[^,]+,\s*plt\\\$coverage\[\s\S]*\\}"
],
"success\_criteria": "When patchwork/cowplot are missing, the code still saves BOTH panels or errors clearly.",
"fix\_hint": "Use gridExtra::grid.arrange(plt\$coverage, plt\$features, ...) + ggsave(); or make patchwork a hard dependency; or save features panel too."
},
{
"id": "C2\_narrowpeak\_conversion\_order",
"title": "narrowPeak 0-based to 1-based conversion applied AFTER filtering",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"np <- np\\\[chr == gene\_chr .*\\]",
"np\\\$start <- pmax\\(np\\\$start \\+ 1"
],
"success\_criteria": "Coordinates converted to 1-based BEFORE any overlap filtering/clipping.",
"fix\_hint": "First: np\[, start := start + 1L]; then filter on gene bounds; then clip."
},
{
"id": "C3\_chr\_name\_mismatch",
"title": "Chromosome naming mismatch likely (e.g., 'chr1' vs '1')",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"np\\\[chr == gene\_chr",
"pb\\\$chr == gene\_chr"
],
"success\_criteria": "Seqnames normalized consistently across GTF, narrowPeak, and BED before comparisons.",
"fix\_hint": "Normalize via helper (sub('^chr','', x)) for both sides before filtering."
},
{
"id": "C4\_bed\_primer\_zero\_based",
"title": "BED primers not shifted from 0-based start",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"primer\_bed.*fread",
"setnames\\(pb, 1:3, c\\("chr","start","end"\\)\\)",
"pb\\\$start <- pmax\\("
],
"success\_criteria": "Primer start incremented by +1 BEFORE clipping.",
"fix\_hint": "Insert: pb\$start <- pb\$start + 1L  (before pmax/pmin)."
},
{
"id": "C5\_narrowpeak\_header\_track\_lines",
"title": "narrowPeak may include 'track' header or comments; fread without guard",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"fread\\(narrowpeak\_file,\s*header\s\*=\s*FALSE",
"comment\s*=\s\*"
],
"success\_criteria": "Either rtracklayer::import() is used or fread() skips comments and drops 'track' lines/NA rows.",
"fix\_hint": "Use comment='#' and drop grepl('^track', chr); or use rtracklayer::import."
},
{
"id": "C6\_selected\_peak\_not\_drawn\_lane",
"title": "Selected peak is only an overlay rectangle; not a dedicated bar lane",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"annotate\\("rect",\s*xmin = window\_start, xmax = window\_end",
"labels = c\\(valid\_transcripts, "Primer", "All Peaks"\\)"
],
"success\_criteria": "A distinct 'Selected Peak' lane exists with geom\_rect spanning window\_start–window\_end.",
"fix\_hint": "Add selected\_lane before primer; draw geom\_rect at that y; update scale\_y labels/breaks."
},
{
"id": "C7\_ggplot\_linewidth\_compat",
"title": "Use of 'linewidth' may drop layers on older ggplot2",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"linewidth\s*="
],
"success\_criteria": "Use size= for geom\_rect/geom\_segment or guard by ggplot2 version.",
"fix\_hint": "Replace linewidth= with size= for compatibility."
},
{
"id": "C8\_visibility\_zoom",
"title": "Peaks visually invisible across whole-gene x-range",
"detection\_strategy": "static",
"look\_for\_patterns": \[
"scale\_x\_continuous\\(limits = c\\(gene\_start, gene\_end\\)\\)",
"coord\_cartesian\\("
],
"success\_criteria": "Feature panel offers zoom to selected window or provides an inset/alternative zoom.",
"fix\_hint": "Add coord\_cartesian(xlim = c(window\_start, window\_end)) to 'feat' (optionally behind a flag)."
}
],
"additional\_validations": \[
{
"id": "V1\_qc\_optional",
"title": "QC subtitle robust across old/new formats",
"expectation": "No error if qc\_tsv missing or missing columns; subtitle informative when present."
},
{
"id": "V2\_bins\_and\_percent",
"title": "Coverage percent scales to TRUE max of BigWig over gene span",
"expectation": "score\_pct defined when cov\_df non-empty; 0–100 clamp applied."
},
{
"id": "V3\_legend\_behavior",
"title": "Legend visible only when all\_peaks\_df non-empty",
"expectation": "legend.position toggles correctly."
}
],
"analysis\_steps": \[
"Parse the R script and build a line-indexed map for precise references.",
"Run regex scans from checks.look\_for\_patterns. Capture first and secondary matches.",
"For each check, derive PASS/FAIL/WARN with brief justification and code excerpts.",
"Synthesize minimal code diffs (fix hints) tailored to matched locations.",
"Produce a final JSON report per deliverables.report.schema."
],
"constraints": {
"do\_not\_execute\_code": true,
"static\_analysis\_only": true,
"no\_external\_downloads": true
},
"severity\_guidance": {
"critical": "Causes peaks or feature panel not to render or to be omitted from output files.",
"high": "Likely to filter away valid peaks or misplace them (coordinate or chrom naming errors).",
"medium": "Visual correctness/compat issues that can make peaks appear missing.",
"low": "Quality-of-life or robustness improvements."
}
}


#batch
Completed at: 08-Sep-2025 16:04:13
Duration    : 5m 26s
CPU hours   : 0.1
Succeeded   : 6

#optimized
Completed at: 08-Sep-2025 15:59:54
Duration    : 5m 14s
CPU hours   : 0.1
Succeeded   : 6

#original
Completed at: 08-Sep-2025 16:06:32
Duration    : 5m 22s
CPU hours   : 0.1
Succeeded   : 6


#Salmon

for fq in *R2.fastq.gz
do
    sample=$(basename $fq _R2.fastq.gz)
    salmon quant -i /data/gent/vo/000/gvo00027/resources/Salmon_index/HomoSapiens/salmon_index \
      -l A \
      -r $fq \
      --fldMean 250 --fldSD 100 \
      -p 8 \
      -o salmon_out/${sample}
done


## gene body coverage only for dominant isoforms
#kallisto - getcdominat isoforms
python bin/dominat_isoform.py \
  -g /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  -q testdata/kallisto_res/kallisto_S7/abundance.tsv testdata/kallisto_res/kallisto_S8/abundance.tsv testdata/kallisto_res/kallisto_S9/abundance.tsv \
     testdata/kallisto_res/kallisto_S10/abundance.tsv testdata/kallisto_res/kallisto_S11/abundance.tsv testdata/kallisto_res/kallisto_S12/abundance.tsv \
  -o dominant_condition_kallisto.gtf


computeMatrix scale-regions \
  -S testdata/Merged_S7_S12.unique.bam.bw \
  -R testdata/dominant_condition_kallisto.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_dominant.gz \
  --outFileNameMatrix matrix_dominant.tab \
  --smartLabels

plotProfile -m matrix_dominant.gz -out geneBody_dominant.png

##splitting into long and short transcripts
awk '$3=="transcript" && $5-$4 < 1000 {print}' dominant_condition_kallisto.gtf > dominant_short.gtf
awk '$3=="transcript" && $5-$4 > 5000 {print}' dominant_condition_kallisto.gtf > dominant_long.gtf

#compute for short transcripts
computeMatrix scale-regions \
  -S Merged_S7_S12.unique.bam.bw \
  -R dominant_short.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_dominant_short.gz \
  --outFileNameMatrix matrix_dominant_short.tab \
  --smartLabels

plotProfile -m matrix_dominant_short.gz -out geneBody_dominant_short.png

#compute for long transcripts
computeMatrix scale-regions \
  -S Merged_S7_S12.unique.bam.bw \
  -R dominant_long.gtf \
  --regionBodyLength 100 \
  --beforeRegionStartLength 0 \
  --afterRegionStartLength 0 \
  --binSize 1 \
  --skipZeros \
  -o matrix_dominant_long.gz \
  --outFileNameMatrix matrix_dominant_long.tab \
  --smartLabels

plotProfile -m matrix_dominant_long.gz -out geneBody_dominant_long.png


#focus on 5' and 3' ends individually
computeMatrix reference-point \
  -S Merged_S7_S12.unique.bam.bw \
  -R dominant_condition_kallisto.gtf \
  --referencePoint TSS \
  -b 500 -a 500 \
  -o matrix_dominant_TSS.gz \
  --skipZeros

plotProfile -m matrix_dominant_TSS.gz -out geneBody_dominant_TSS.png --outFileNameData profile_TSS.tab

computeMatrix reference-point \
  -S Merged_S7_S12.unique.bam.bw \
  -R dominant_condition_kallisto.gtf \
  --referencePoint TES \
  -b 1000 -a 1000 \
  -o matrix_dominant_TES.gz \
  --skipZeros
  
plotProfile -m matrix_dominant_TES.gz -out geneBody_dominant_TES.png --outFileNameData profile_TES.tab


#with MAQCA_1 an TSS

computeMatrix reference-point \
  -S MAQCA_merged.bam.bw \
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --referencePoint TSS \
  -b 500 -a 500 \
  -o matrix_MAQCA_merged_TSS.gz \
  --skipZeros
plotProfile -m matrix_MAQCA_merged_TSS.gz -out geneBody_MAQCA_merged_TSS.png --outFileNameData profile_MAQCA_TSS.tab



computeMatrix reference-point \
  -S MAQCA_merged.bam.bw \
  -R /data/gent/vo/000/gvo00027/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf \
  --referencePoint TES \
  -b 500 -a 500 \
  -o matrix_MAQCA_merged_TES.gz \
  --skipZeros

plotProfile -m matrix_MAQCA_merged_TES.gz -out geneBody_MAQCA_merged_TES.png --outFileNameData profile_MAQCA_TES.tab

python bin/TSSratio.py \
  --tss testdata/profile_MAQCA_TSS.tab \
  --tes testdata/profile_MAQCA_TES.tab \
  -o ratios_MAQCA.txt
