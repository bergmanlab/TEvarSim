# TEvarSim

**TEvarSim** is a versatile genome simulation tool for generating polymorphic transposable element (TE) variants.   

Key features:
- Supports TE insertions, deletions, and LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs)
-	Simulates both real and random TE variants
-	Simulates both short- and long-read sequencing data
-	Supports population-scale genome simulation
-	Generates VCF files
-	Includes scripts to compare predicted vs. simulated variants, per genome (`Compare`) or across every simulated locus at once (`Evaluate`)
---

## Installation

The version of TEvarSim available through pip is outdated; install it from the github repo instead.

```bash
# gfatools and repeatmasker are used for *TEpan*.
# mason and pbsim3 are used for short-reads simulation and long-reads simulation, respectively.
# You may skip installing the software if you do not use the corresponding functionality.
conda create -n tevarsim -c bioconda gfatools repeatmasker mason pbsim3
conda activate tevarsim
#git clone TEvarSim
cd TEvarSim
pip install .
```
## Quick start
Example data can be found in the **testData** directory   

The number of each event class is specified explicitly with `--nINS` (insertions), `--nDEL`
(deletions), and `--nEXC` (LTR-LTR recombinations / excisions). Set at least one to a value > 0.

**1. Simulate 4 insertions and 2 deletions from known TEs**
```bash
tevarsim TEreal --knownINS MEI.fa --existingTEs rptmsk.out --CHR 21 --nINS 4 --nDEL 2
```
- `MEI.fa` is known pTE insertion, from paper [Logsdon, G.A. et al. Nature, 2025](https://www.nature.com/articles/s41586-025-09140-6)  
- `rptmsk.out` is the repeatmasker output for chr21_tiny.fa.

**1b. Add LTR-LTR recombinations (excisions)**
```bash
tevarsim TEreal --knownINS MEI.fa --existingTEs ltr_element.out --CHR 21 \
    --TEtype Gypsy --nDEL 1 --nEXC 1
```
- An excision converts a full-length LTR element (LTR-Internal-LTR) into a solo LTR.
- Full-length elements are identified from the fragment structure in a RepeatMasker `.out`
  file, so `--nEXC` requires `.out` input (not UCSC `.txt`). Include the relevant LTR family
  in `--TEtype`. `testData/ltr_element.out` is a small demonstration fixture.

**2. Simulate deletions plus random TE insertions**
```bash
tevarsim TErandom --consensus human_TE.fa --existingTEs rptmsk.out --ref chr21_tiny.fa --nINS 4 --nDEL 2 --CHR chr21,chr22
```
- `TEconsensus.fa` is human TE consensus sequences from Dfam
- Add `--nEXC N` (with a RepeatMasker `.out` `--existingTEs` file) to also simulate excisions.
- You can try adding `--CHR chr21,chr22` to restrict the analysis to just two 'chromosomes.'
- You can try adding `--regions regions.bed` or `--exclude regions.bed` to test out the functionality to restrict the simulation to certain regions or exclude them from the simulation.

**3. Simulate pTE from pangenome graph**
```bash
# Fetch pangenome graph from HPRC
curl https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-grch38.gfa.gz > hprc-v1.0-minigraph-grch38.gfa.gz
tevarsim TEpan --gfa hprc-v1.0-minigraph-grch38.gfa.gz --lib Homo_sapiens_DFAM.fa  --CHR chr21 --nINS 4 --nDEL 2
```
- `hprc-v1.0-minigraph-grch38.gfa.gz` is downloaded from [HPRC](https://data.humanpangenome.org/alignments)
- LTR-LTR recombination (`--nEXC`) is not yet supported for TEpan.

**4. Simulate 10 genomes with 6 pTE**  
```bash
tevarsim Simulate --ref chr21_tiny.fa --bed TEreal.bed --num 10 --pool MEI.fa
# if you want to generate sequence vairiations of the same TE between genomes, run below commonds
tevarsim Simulate --ref chr21_tiny.fa --bed TEreal.bed --num 10 --pool MEI.fa --diverse --diverse_config diverse.config
```
- `chr21_tiny.fa` is the reference sequence
- `real.bed` is the position of pTE positions that generated from `tevarsim TEreal`
- `diverse` : Introduce sequence diversity among individuals for the same TE event (which is suitable for evaluating methods that require a TE panel as input)
- `diverse_config` : A configuration file of parameters for introducing sequence diversity among individuals for the same TE event (optional)

**5. Generate sequencing reads from simulated genome** 
```bash
tevarsim Readsim --type short --genome Sim.fa --depth 10 
tevarsim Readsim --type long --genome Sim.fa --depth 10
```
- `type` : short reads or long reads

**6. Benchmark a caller against the simulation**
```bash
# one genome at a time: how well was this genome called?
tevarsim Compare --truth Sim.vcf --pred calls.vcf --truthID Sim_0 --predID Sim_0 -O compare_Sim_0
# every simulated locus at once, over all genomes: which were recovered, and where does it fail?
tevarsim Evaluate --truth Sim.vcf --pred calls.vcf -O evaluate
```
- `Compare` needs a genome named on both sides; `Evaluate` pairs the genomes by name itself and
  reports per-locus results aggregated across all of them

**7. The complete workflow** 
- We offered the complete workflow for the short/long-read based tool benchmarking. Please see the file `workflow.sh`

## Flowchart
![flowchart](https://github.com/JanMiao/TEvarSim/blob/main/chart.png)  
- The known TE deletion information can be obtained from [UCSC annotaion file (.txt)](https://genome.ucsc.edu/cgi-bin/hgTables) or [repeatmasker annotation (.out)  ](https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)
- The known TE insertion position can be obtained from our pre-built dataset (data/MEI_Callset_GRCh38.ALL.20241211.fasta). Any TE insertion sequence is acceptable , as long as the sequence ID follows the naming format **CHR-POS-ID**, e.g., **chr1-683234-AluSp#SINE/Alu**

## Usage
TEvarSim provides seven main command-line subcommands:
```bash
tevarsim <subcommand> [options]
```

### 1. TErandom
Generate pTE position from known deletion sites and random TE insertion.

**Required arguments:**
- `consensus` : Path to the TE consensus FASTA file. Headers may carry a `TSD=` tag after the
  sequence ID (see `Simulate --tsd-from-header`); it is propagated to the TE pool. The sequenceIDs in the FASTA header should be >TEname#class/superfamily, e.g., >AluY#SINE/Alu
- `ref` : Reference genome FASTA
  
**Optional arguments:**
- `existingTEs` : Known TE annotation file (RepeatMasker .out or UCSC .txt); required when `--nDEL` or `--nEXC` > 0
- `CHR` : Chromosome to simulate TE insertions on (e.g., chr21 or 21)
- `nINS` : Number of random TE insertions to simulate (default: 0)
- `nDEL` : Number of TE deletions to simulate from `--existingTEs` (default: 0)
- `nEXC` : Number of LTR-LTR recombinations (excisions of full-length LTR elements) to simulate; requires a RepeatMasker `.out` `--existingTEs` file (default: 0)
- `nNM` : Number of non-mobilizing SVs: deletions that cut into or swallow an `--existingTEs` element by no TE mechanism (default: 0). See below
- `nmKinds` : Which kinds of non-mobilizing SV, comma-separated, dealt out in equal shares (default: `internal,left,right,tight,wide`)
- `nmInternal / --nmStraddleFlank / --nmTightFlank / --nmWideFlank` : `MIN-MAX` lengths in bp: an internal deletion (100-1500, never more than half the element), the flank a left or right one takes (25-2000), and the flank a tight or wide one takes on each side (25-250, 250-5000)
- `nmMargin` : A breakpoint inside an element stays this many bp from its ends (default: 50)
- `outprefix` : Output prefix for TE pool FASTA (default: TErandom)
- `TEtype` : Which TE super families to be extracted from the TE deletion file (default: Alu, L1, and SVA). Specify the TE type by `--TEtype Alu --TEtype L1`
- `DELlen` : A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)
- `nMIN` : A minimum number of TE deletions for each TE super family to be simulated (default: 0)
- `TEdistance` : A minimum length of distance between two TE insertions (default: 500 bp)
- `nSV` : Number of background structural variants to simulate (default: 0)
- `regions` : Bed file with regions to include in the simulation.
- `exclude` : Bed file with regions to exclude from the simulation.
- `snp-rate` : SNP mutation rate per base (default: 0.02)
- `indel-rate` : INDEL mutation rate per base (default: 0.005)
- `indel-ins` : Proportion of indels that are insertions (default: 0.4)
- `indel-geom-p` : Geometric distribution parameter for indel lengths (default: 0.7)
- `truncated-ratio` : Proportion of sequences to truncate (default: 0.3)
- `truncated-max-length` : Maximum proportion of sequence to truncate (default: 0.5)
- `polyA-ratio` : Proportion of sequences to add polyA tail (default: 0.8)
- `polyA-min` : Minimum polyA length (default: 5)
- `polyA-max` : Maximum polyA length (default: 20)
- `seed` : Random seed (default: None)


**Non-mobilizing SVs.** A TE caller sees every structural variant that involves an element, not
only the ones a transposon made. `--nNM` adds such variants, deletions none of which any TE mechanism
makes, each on a reference element (an `--existingTEs` candidate that no deletion or excision took):

| kind | the deletion |
| :--- | :--- |
| `internal` | inside the element, neither end at the element's own: not a truncation, an excision or a presence/absence |
| `left` | from the left flank into the element |
| `right` | from inside the element into the right flank |
| `tight` | the whole element and 25-250 bp of flank either side, more than any TSD |
| `wide` | the whole element and 250-5000 bp of flank either side, taking whatever else lies there |

Each is written to the BED as `NM-<kind>-<host element>`. Simulate puts them into the genomes, draws
their allele frequencies from the insertion spec (they are new mutations, so rare), and writes them to
`<outprefix>.nonmobilizing.vcf` rather than to the TE truth. `Evaluate --nonmobilizing` then counts
the calls a caller made at them. They draw from a random stream of their own, so the deletions,
excisions and insertions are those a run without them would make.

### 2. TEreal
Automatically generate pTE positions from RepeatMasker or UCSC repeat annotations.

**Required arguments:**  
- `knownINS` : Known TE insertion file (The sequence ID should follow the naming format CHR-POS-ID, e.g., chr1-683234-AluSp#SINE/Alu)  
- `existingTEs` : Known TE annotation file (from RepeatMasker `.out` or UCSC `.txt`)  
- `CHR` : Chromosome used to simulate pTE  

**Optional arguments:**  
- `nINS` : Number of TE insertions to simulate from `--knownINS` (default: 0)
- `nDEL` : Number of TE deletions to simulate from `--existingTEs` (default: 0)
- `nEXC` : Number of LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs); requires a RepeatMasker `.out` `--existingTEs` file, and the LTR family must be included in `--TEtype` (default: 0)
- `DELlen` : A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)
- `nMIN` : A minimum number of TE deletions for each TE super family to be simulated (default: 0)
- `nSV` : Number of background structural variants to simulate (default: 0)
- `outprefix` : Output prefix for BED file (default: real)  
- `TEtype` : Which TE super families to be extracted from the TE deletion file (default: Alu, L1, and SVA). Specify the TE type by `--TEtype Alu --TEtype L1`
- `seed` : Random seed (default: None)  

### 3. TEpan
Generate pTE position from Pangenome graph.

**Required arguments:**  
- `gfa` : GFA file of the pangenome graph
- `lib` : RepeatMasker library file [pre-formatted libraries](https://www.repeatmasker.org/~cgoubert/GraffiTE_libraries/).
- `CHR` : Chromosome used to simulate pTE  

**Optional arguments:**  
- `outprefix` : Output prefix for BED file (default: TEpan)
- `nINS` : Number of TE insertions to simulate from the pangenome graph (default: 0)
- `nDEL` : Number of TE deletions to simulate from the pangenome graph (default: 0)
- `minLen` : Minimum length of structural variants to consider (default: 250)  
- `cov` : Minimum TE coverage to consider a structural variant as TE (0-1, default: 0.5)
- `tmpDir` : Temporary directory for intermediate files (default: tmp_TEpan)
- `TEtype` : TEs to be extracted from the TE deletion file, with the default set as LINE, SINE, LTR, and RC. Specify the TE type by `--TEtype LINE --TEtype SINE`
- `seed` : Random seed (default: None)  

### 4. Simulate
Simulate pTE insertions/deletions and generate VCF and modified genome FASTA.

**Required arguments:**
- `ref` : Reference genome FASTA  
- `te-pool` : TE pool FASTA  
- `bed` : BED file of TE positions (can be generated by `ppte TEreal`)  
- `num` : Number of simulated genomes

**Optional arguments:**  
- `diverse` : Introduce sequence diversity among individuals for the same TE event, which is suitable for evaluating methods that require a TE panel as input(default: False)
- `diverse_config` : Path to a configuration file for introducing sequence diversity among individuals for the same TE event (default: None; requires --diverse)
- `outprefix` : Output prefix (default: Sim)  
- `af-min / --af-max` : Bounds on the allele frequency of a simulated variant. Under the default
  `--af-dist uniform` they are the range it is drawn flat from; under `--af-dist exponential`
  they truncate the exponential (default: 0.1/0.9)
- `af-dist` : `uniform`, `exponential` or `reflected_exponential` (default: `uniform`). An exponential gives a population
  dominated by rare, largely private insertions with a thin tail of common ones, which is roughly
  what a real TE landscape looks like; a flat draw makes most insertions shared
- `af-mean` : Mean of that exponential. Required by `--af-dist exponential`, rejected without it
- `del-af-dist / --del-af-mean / --del-af-min / --del-af-max` : Allele frequency spec for
  DELETIONS only; each part defaults to its `--af-*` counterpart
- `allow-zero-carriers` : Keep events that no genome drew, instead of giving each one a carrier
- `tsd-min / --tsd-max` : Min/max TSD length, used for every element that does not state its
  own (default: 5/20)
- `tsd-from-header` : Take each insertion's TSD length from a `TSD=` tag on its pool FASTA header
- `sense-strand-ratio` : Proportion of sense-strand insertions (default: 0.5)  
- `bg-pi` : Background variation -- SNPs and short indels shared among the genomes by descent --
  at this pairwise diversity per bp (default: 0, none). See below
- `bg-sv-rate` : Background structural variants -- deletions and tandem duplications -- per genome per Mb
  against the reference, on the same genealogies (default: 0, none). See below
- `bg-sv-dup-frac / --bg-sv-min / --bg-sv-max` : Fraction that are duplications (0.5), and the log-uniform
  length range (50-10,000 bp)
- `bg-indel-frac / --bg-indel-max / --bg-margin / --bg-block` : Fraction of background variants
  that are indels (0.1), their longest length (50), how far they stay from any TE event (30 bp),
  and the block length given its own genealogy (50,000 bp)
- `seed` : Random seed (default: None)  

**Background variation.** Genomes that differ only by the events simulated never test whether a
caller can tell an element from everything else that varies at a locus. `--bg-pi` adds that
background: each `--bg-block` gets a random genealogy of the genomes (Kingman coalescent; one per
block stands in for recombination), and mutations fall on its branches at a rate that gives
pairwise diversity `--bg-pi` per bp, each carried by every genome below its branch, so genomes
share background by descent as real strains do. A `--bg-indel-frac` share are indels (half
insertions of random sequence, half deletions, geometric length up to `--bg-indel-max`). The
reference is the genealogy's root. Background stays `--bg-margin` bp clear of every TE event, so
TSDs and scored alleles are the reference's, and is written to `<outprefix>.background.vcf`, never
to the TE truth `<outprefix>.vcf`; TErandom `--nSV`'s synthetic insertions and deletions are
written there too. It draws from a random stream of its own: with or without it, the TE truth is
the same.

**Background SVs.** `--bg-sv-rate` adds deletions and tandem duplications that no TE made, at that
many per genome per Mb against the reference, on the same per-block genealogies as the background SNPs,
so genomes share them by descent. Lengths are log-uniform from `--bg-sv-min` to `--bg-sv-max`; they land
anywhere, over the reference's own elements as readily as elsewhere, but never within `--bg-margin` of a
TE event or of one another, nor over an N. A duplication is a copy of its stretch inserted right after it.
They are written to `<outprefix>.background_sv.vcf` (a duplication's record sits at the end of the copied
stretch, with `INFO/DUPSTART` at its start), never to the TE truth, and `Evaluate --background_sv` counts
the calls a caller made at them, by type and size.

**Per-element TSD.** TSD length is a property of the clade, not a constant: a cut-and-paste
transposon's TSD is set by the stagger between its transposase's two cuts, so hAT and P
duplicate 8 bp and Tc1/mariner 2 bp every time, while a non-LTR element inserting by
target-primed reverse transcription has no fixed stagger and produces a spread, and a Helitron
copies in by rolling-circle replication and duplicates nothing. One `--tsd-min`/`--tsd-max`
range for a whole library cannot say that, so an element can state its own in its library
header:

```
>copia#LTR/Copia TSD=5
>P-element#DNA/P TSD=8
>Bari1#DNA/Tc1-Mariner TSD=2
>Helitron#RC/Helitron TSD=0
>jockey#LINE/Jockey TSD=5-15
```

`TSD=N` is a fixed length and `TSD=N-M` an inclusive range drawn uniformly; `TSD=0` means no
duplication, which is a real value rather than a missing one. The tag lives in the description
-- everything after the sequence ID -- so it does not disturb the `TEname#class/superfamily`
grammar. `TErandom` copies it from `--consensus` onto every pool record it writes, including
the synthesised `NAME-FULL` records, which take the `-LTR` part's tag. A tag that is present
but malformed stops `TErandom`; that is the only command that reads the library, so it is the
only place a typo can be caught.

`Simulate --tsd-from-header` then uses it. Elements with no tag fall back to
`--tsd-min`/`--tsd-max` and are named in a single warning. It applies to insertions only:
deletions and excisions have no pool record and duplicate nothing, so their `INFO/TSD` remains
a draw from the global range and remains meaningless.

**Allele frequency.** `--af-dist exponential --af-mean M` draws each event's frequency from an
exponential of mean `M`, truncated to `[--af-min, --af-max]`. Truncating an exponential from
below also shifts it, so `--af-mean` is the mean *above* `--af-min` and the realised mean
frequency is about `--af-min + --af-mean`. Keep `--af-min` near `1/--num`, or at 0, unless you
mean to raise the floor. `--num` also bounds the resolution: with 20 genomes the smallest
non-zero frequency is 0.05, so a much smaller `--af-mean` than that just collapses onto it.

**Deletions measure the opposite thing from insertions.** For an insertion `GT=1` means the
genome *has* the element; for a deletion it means the genome *lacks* the reference's. So a
deletion's allele frequency is the frequency of an **absence**, and a reference element that
most genomes lack -- a reference-private insertion, seen from the other side -- is a deletion
at *high* frequency. Sharing one distribution between the two therefore gets one of them
backwards, which is what `--del-af-*` exists to fix.

An ordinary exponential cannot express it: asking for a high mean flattens the distribution
toward uniform rather than concentrating mass at the top. `reflected_exponential` is the same
exponential mirrored end for end within `[min, max]`, so the mass piles against the top and the
mean sits at `max` minus the mean. To give the reference its own private insertions:

```
tevarsim Simulate ... --nDEL 1800 \
    --af-dist exponential --af-mean 0.01 --af-min 0.01 \
    --del-af-dist reflected_exponential --del-af-max 1.0
```

Excisions stay on `--af-dist`: an excision is a derived event that arises and spreads from
rare, like an insertion, so only `DEL` carries the inverted sense. The split draw is taken only
when a `--del-af-*` option is actually given, so runs that do not use it are unaffected.

**Every event reaches a genome.** An event whose genotype draw gives no carrier used to be
dropped from the VCF and from every genome without comment, so at low frequencies most of the
pool evaporated -- exactly the regime that is supposed to produce private insertions. Such an
event is now given one carrier, a genome chosen at random, and the count is reported. This
applies to deletions and excisions as well as insertions: a requested event that reaches nobody
is the same silent loss whatever its type. The realised minimum allele frequency is therefore
`1/--num`. Pass `--allow-zero-carriers` for the old behaviour.

**Event types in the output VCF** (`INFO/TYPE`):
- `INS` : a TE insertion (point event; the inserted sequence comes from the TE pool).
- `DEL` : a TE deletion (a reference span is removed).
- `EXC` : an LTR-LTR recombination. `POS` is the 5' start of the full-length element, `REF` is
  the full `LTR-Internal-LTR` element, and `ALT` is the solo LTR left behind. `INFO/LTRLEN`
  gives the solo LTR length and `INFO/SVLEN` the (negative) net length change. Genotype 0
  keeps the full element; genotype 1 carries the solo LTR. Excision events are encoded in the
  BED with an `EXC-` ID prefix and a 7th column carrying the LTR length.


### 5. Readsim
Generate short or long reads from the simulated genome.

**Required arguments:**
- `type` : Type of reads to simulate (short or long)  
- `genome` : Reference genome file (FASTA) where reads will be simulated from  
- `depth` : Depth of simulated reads    

**Optional arguments(long reads only):**  
- `Lerror` : Sequencing error rate for long reads (default: 0.15)  
- `Lmean` : Average read length for long reads (default: 9000)   
- `Lstd` : Standard deviation of read length for long reads (default: 7000)    

**Optional arguments(short reads only):**  
- `length` : Read length (default: 150)  
- `Fmean` : Average fragment length (default: 300)   
- `Fstd` : Fragment length standard deviation (default: 30)    

**Optional arguments(general):**  
- `seed` : Random seed for reproducibility (default: None)

### 6. Compare
Compare predicted VCF to the simulated VCF.

**Required arguments:**  
- `truth` : Ground truth VCF (generated by `ppte simulate`)  
- `pred` : Predicted file
- `predType` : Type of the predicted file (VCF or BED, default: VCF)
- `outprefix` : the Matched variants in two VCF file  
- `truthID` : Sample ID in truth VCF  
- `predID` : Sample ID in predicted VCF  

**Optional arguments:**  
- `TEtype`: TE family in truth VCF to consider in the comparison. This is the family, not the
  superfamily — TY1, TY2, TY4 and TY5 are all LTR/Copia. The family is read out of the variant ID,
  ignoring the source locus and the polyA/strand annotation that a known-insertion pool decorates
  it with, so `--TEtype AluY` matches every insertion of `chr*-*-len*-AluY#SINE/Alu-polyA*-strand*`
- `INSonly`: Only compare insertions in truth VCF
- `nHap` : Ploidy level, e.g., 2 for Humans  (default: 2)  
- `max_dist` : Maximum allowed distance for variant matching (default: 100 bp)  

**Example**:
```bash
tevarsim compare --truth sim.vcf --pred variants.vcf --truthID Hap1_Hap2 --predID Sample
```
In simulation files, genomes are named Hap1, Hap2, etc.; for polyploids, combine haplotype IDs with `_` for one individual.  

### 7. Evaluate
Score a prediction against **every simulated locus at once**, aggregated over all genomes.

`Compare` answers "how well was *this* genome called?" — it takes one truth sample and one
prediction sample, so benchmarking a whole simulation means looping it over every genome and
stitching the per-genome outputs back together by hand. `Evaluate` turns the question around:
the unit of analysis is the simulated **locus**, not the genome. Every record in the truth VCF —
one mobile element and everything that happened to it — is scored once, using all the simulated
genomes at the same time, and the per-locus results are aggregated into overall metrics plus
breakdowns by event type, TE family, element size and allele frequency.

A locus is not an event: an element inserted in one generation and excised in a later one is one
locus carrying two events. Every count is of loci except the event-type breakdown, whose rows total
the events and so total more than the loci. No genome is named on the command line.

**Required arguments:**
- `truth` : Ground truth VCF (generated by `tevarsim Simulate`)
- `pred` : Predicted file

**Optional arguments:**
- `predType` : Type of the predicted file (VCF or BED, default: VCF). A BED prediction carries no
  allele sequences and no sample columns, so only detection and breakpoint accuracy are reported
- `outprefix` : Output prefix for the per-locus JSON and the HTML report (default: `evaluate`)
- `no_html` : Skip the HTML report; write only the text summary, the JSON and the per-locus files
- `sample_map` : Two-column file pairing truth and prediction samples (`truth_sample<TAB>pred_sample`).
  Default: pair samples by name
- `TEtype` : TE family in truth VCF to consider, e.g. `TY1-FULL` (case-insensitive). This is the
  family, not the superfamily — TY1, TY2, TY4 and TY5 are all LTR/Copia. Default: all families
- `INSonly` : Only evaluate insertions in truth VCF
- `nHap` : Number of haplotype columns per individual in the truth VCF; >1 merges them before
  pairing with the prediction (default: 1)
- `max_dist` : Maximum allowed distance (bp) to consider a prediction and a simulated locus the
  same locus (default: 100)
- `gt_len_tol` : Maximum allowed difference (bp) in allele length to consider two alleles the same
  (default: 50)
- `background_sv` : Simulate's `<outprefix>.background_sv.vcf`. Reports how many of the background SVs the
  prediction called, by type and size; scored as `--nonmobilizing` is
- `nonmobilizing` : Simulate's `<outprefix>.nonmobilizing.vcf`. Reports how many of the non-mobilizing
  SVs the prediction called, by kind: every unmatched call within `--max_dist` of one is charged to it.
  Those calls stay unmatched predictions, so they count against precision either way
- `size_bins` : Upper edges (bp) of the element-size strata; repeat the flag per edge
  (default: 100 500 1000 5000 10000)
- `af_bins` : Upper edges of the allele-frequency strata; repeat the flag per edge
  (default: 0.1 0.25 0.5 0.75)

**Example**:
```bash
tevarsim Evaluate --truth sim.vcf --pred variants.vcf --INSonly --TEtype TY1-FULL -O eval
```

**How genomes are paired.** Truth and prediction samples are paired by name, which is the pairing
a per-genome `Compare` loop is normally driven with (`-I $sample -J $sample`); as a convenience a
lone sample on each side is paired whatever it is called. Nothing is ever paired by column order —
a silently wrong pairing would report confident, meaningless genotype accuracy. Use `--sample_map`
when the two files name their genomes differently. If nothing can be paired, detection and
breakpoint accuracy are still reported and the carrier statistics are marked unavailable.

**Displaced excision calls.** A caller that anchors an SV where the two alleles begin to differ
puts an excision at the far end of the surviving solo LTR, because the LTR itself is identical on
both sides and collapses. The simulation writes `POS` at the start of the element, which is where
the solo LTR starts, so the two conventions differ by exactly the LTR length — 338 bp on Ty1, far
beyond any sane `--max_dist`. Left alone that scores one excision twice over: once as a miss and
once as a false positive.

Evaluate pairs such a call with its locus using `INFO/LTRLEN` (falling back to the length of the
excision allele), and records it as **displaced** rather than as a match. A displaced call earns
detection credit and nothing else: the locus counts as recovered, but the call is in the wrong
place, so it contributes no allele, breakpoint, carrier or genotype statistics, and its record
stays among the unmatched predictions — precision is not credited for it. The report breaks the
recovered count into `correctly placed` and `displaced`, and the per-locus file for the call
carries `displaced_match_for` so it is not mistaken for a spurious one.

**What is reported.** For each simulated locus:
- whether the prediction recovered the locus at all, and whether the allele it called is the same
  one that was simulated (by sequence, or by net length change within `--gt_len_tol`);
- the breakpoint offset and allele length error of the matched call, both **signed** — a call
  landing short of the simulated value reads negative and one past it positive, reported as a mean
  and SD so a systematic bias shows up as a mean away from zero instead of being folded in with
  random jitter;
- which genomes were called as carriers versus which genomes actually carry it (carrier TP/FP/FN);
- whether the genotype of each paired genome agrees, by allele content.

Events and prediction records are matched one-to-one, preferring allele agreement over breakpoint
proximity. That is what keeps a stacked pair of elements honest: two elements at one locus need two
prediction records for both to count as recovered, rather than both claiming the single record that
happens to sit there.

**Output.** `--outprefix` writes three things:

- `<outprefix>.html` — the report. See below; pass `--no_html` to skip it.

- `<outprefix>.json` — the run's parameters, the aggregate summary, and one object per simulated
  locus, so results can be reloaded for plotting or regression checks.
- `<outprefix>_loci/` — one self-contained JSON per locus, named `<chrom>_<pos>.json` after where
  the locus was simulated. Every simulated locus gets a file, including the ones the prediction
  never called, and so does every prediction that matched no simulated locus (`"kind":
  "unmatched_prediction"`, recording its nearest simulated locus so a call landing just outside
  `--max_dist` is distinguishable from one with nothing simulated near it). Two elements stacked at
  one position take a `-2`, `-3` suffix. Each file repeats enough of the run to be read on its own;
  strip its `run` and `kind` keys and what is left is exactly that locus's entry in the summary's
  `loci` list, which also carries a `locus_file` pointing back at it. Files left by an earlier run
  of the same prefix are cleared first, so the directory only ever describes the current run.

**The HTML report.** `<outprefix>.html` shows the shape the text summary can only average.
"Breakpoint offset: mean +0.07 bp, SD 0.60 bp" cannot separate a caller that is exact on every
locus but 30 bp out on two of them from one that scatters a base either side of the truth on all
of them; those are different problems, and only the distribution tells them apart. The page holds

- **breakpoint resolution** — the signed error of every correctly placed call as a histogram,
  centred so an exact call gets a bar of its own and symmetric about zero so a bias reads as an
  off-centre mass, beside the cumulative `exact / within N bp` readout;
- **allele length accuracy** — the same for the called allele's length, which is where a
  well-placed call that truncated the element shows up;
- **recovery by event size** — every simulated locus at the size it was simulated at, each column
  split into correctly placed / displaced / missed, so a sensitivity floor at a particular element
  size (the solo LTRs, typically) is visible instead of averaged away;
- **the breakdowns** — recovery rate, carrier F1 and genotype concordance as grouped bars beside
  the same sortable tables the text report prints;
- **truth and prediction side by side** — a second tab, reached by clicking either VCF's file name
  in the heading, or any locus named in a figure's table. Both files whole, in two panes. Clicking a record jumps
  the other pane to the record it was paired with and marks both; where a record has no partner,
  the other pane says so rather than sitting still, since "there is nothing there" and "the viewer
  did nothing" look identical in a pane that fails to move. Scrolling either pane carries the
  other to the same *relative* position, so two files of different lengths stay roughly aligned
  without pretending to a record-for-record correspondence they do not have. Records the filters
  dropped, and prediction records no sample calls, are kept and marked unscored — omitting them
  would show a smaller run than the one that happened. Allele sequences over 24 bp are shown as
  their first bases and their length; a 6 kb ALT per record would otherwise put megabytes of
  sequence in the page.

Each figure also carries the table of values behind it, naming the loci in every bar by
`chrom:pos` — the same key the per-locus JSON files are named after, so a bar leads to the loci
that made it, and each of those to its own file and to its record in the viewer.

The tabs are switched on by the page's own script: with JavaScript off it stays one long
document where every link is an ordinary anchor that still resolves, rather than stranding the
viewer behind a control that cannot be operated.

The figures are modelled on the Tangram paper's benchmarking figures (Wu et al. 2014,
*BMC Genomics* 15:795, Figures 1 and 2). Everything is generated in Python as inline SVG: there is
no plotting dependency to install and nothing is fetched when the page is opened, so the report
works over a mounted filesystem or copied off a cluster. It follows the browser's light/dark
setting.

**Note on the event-type breakdown.** A locus is counted under every event class in its history,
so the rows total the events rather than the loci: an element that was inserted and later excised is
one locus counted under both `INS` and `EXC`. This is read from `INFO/EVENTTYPE` (`Number=A`, one value
per ALT allele), not `INFO/TYPE` (`Number=1`), which names only the event that created the record —
by `TYPE` alone every excision of a simulated element would be filed under `INS` and the `EXC` row
would count only excisions of elements the reference already carried. A class is credited only where
the allele naming it has the matching shape: an `INS` allele has to add sequence relative to the
allele before it and an `EXC` allele has to remove it, by more than `--gt_len_tol`. Anything
declared without such an allele is left out of the count and listed in the report and in
`summary.unsupported_events`.

**Note on the TE family.** The family is parsed exactly as `Compare --TEtype` parses it, so the
labels in the "By TE family" table are the values that flag accepts: `AluY` for an insertion whose
pool sequence ID is `chr21-211282-len320-AluY#SINE/Alu-polyA23-strand+`, since neither the source
locus nor the polyA/strand annotation says which family the element belongs to. A "By TE
superfamily" table is added whenever the superfamily groups events that the family table does not
(the AluY / AluYa5 / AluYb8 of a human simulation all being `Alu`).

