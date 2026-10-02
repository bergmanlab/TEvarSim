# __main__.py
import argparse
import sys
import os
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Global interpreter lock.*")
from . import build_pool, TE_real, simulate, compare_vcf, evaluate, reads, TEpan
from . import __version__

class positive_int(int):
    def __new__(self, value):
        self = super().__new__(self, value)
        if self < 0:
            raise ValueError(f"Must be a positive integer: {value}")
        return self

class ratio(float):
    def __new__(self, value):
        self = super().__new__(self, value)
        if not 0 <= self <= 1:
            raise ValueError(f"Rate or ratio must be between 0 and 1: {value}")
        return self

class File_Path(str):
    def __new__(self, value: str):
        return super().__new__(self, os.path.abspath(value))
    def __bool__(self):
        try:
            return os.path.isfile(self) and os.stat(self).st_size != 0
        except (AttributeError, OSError, ValueError, TypeError):
            return False
class Existing_File_Path(File_Path):
    def __new__(self, value: str):
        path = super().__new__(self,value)
        if not path:
            raise ValueError(f"{value} is not an existing file path")
        return path
class Fasta_File_Path(Existing_File_Path):
    def __new__(self, value: str):
        path = super().__new__(self,value)
        check = False
        with open(path) as fin:
            # Iterate the file object rather than fin.readlines(). The two are identical
            # here -- the loop breaks at the first ">" either way -- but readlines()
            # materialises the WHOLE file as a list of strings before the loop starts, so
            # the break saves nothing. On a 3.1Gb human reference that is 5.2GB of resident
            # memory and 7 seconds, paid by every command that takes a --ref or a --pool,
            # before any work begins. It goes unnoticed on smaller genomes: the same code
            # costs ~230MB on a 137Mb Drosophila reference and reads as ordinary overhead.
            for line in fin:
                if line[0] == ">":
                    check = True
                    break
        if check is False:
            raise ValueError(f"{value} is not in fasta format")
        return path

def main():
    parser = argparse.ArgumentParser(prog="tevarsim", 
                                     description=f"TEvarSim: A genome simulation tool for transposable element (TE) variants\nVersion: {__version__}")
    parser.add_argument("--version", "-V", action="version", version=f"%(prog)s {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. TErandom
    TErandom_parser = subparsers.add_parser("TErandom", 
                               help="Generate pTE position from known deletion sites and random TE insertion")
    # Base
    TErandom_parser.add_argument("--ref", "-F", type=Fasta_File_Path, required=True, 
                    help="Reference genome FASTA file")
    TErandom_parser.add_argument("--consensus", "-C", type=Fasta_File_Path,  required=True,
                    help="Path to the TE consensus FASTA file. The sequenceIDs in the FASTA header should be >TEname#class/superfamily, e.g., >AluY#SINE/Alu")
    TErandom_parser.add_argument("--existingTEs", "-L", type=Existing_File_Path, required=False, 
                    help="Input known TE deletion file (RepeatMasker .out or UCSC .txt)")
    TErandom_parser.add_argument("--CHR", "-H", type=str,
                    help="Comma-separated list of chromosome(s) to simulate TE insertions on (e.g., chr21,chr22,chr23). Default: all")
    TErandom_parser.add_argument("--TEtype", "-e", type=str, action="append",
                    help="Which TE super families to be extracted from the TE deletion file. Default: all")
    TErandom_parser.add_argument("--nINS", "-N", type=int, default=0,
                    help="Number of random TE insertions to simulate (default: 0)")
    TErandom_parser.add_argument("--nDEL", type=int, default=0,
                    help="Number of TE deletions to simulate from --existingTEs (default: 0)")
    TErandom_parser.add_argument("--nEXC", type=int, default=0,
                    help="Number of LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs) to simulate from --existingTEs; requires a RepeatMasker .out file (default: 0)")
    TErandom_parser.add_argument("--outprefix", "-O", type=File_Path, default="TErandom", 
                    help="Output prefix for the generated TE pool FASTA file and the bed file (default: TErandom)")
    TErandom_parser.add_argument("--DELlen", type=int, default=100,
                    help="A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)")
    TErandom_parser.add_argument("--nMIN", type=int, default=0,
                    help="A minimum number of TE deletions for each TE super family to be simulated (default: 0)")
    TErandom_parser.add_argument("--TEdistance", type=positive_int, default=500,
                    help="A minimum length of distance between two TE insertions (default: 500)")
    TErandom_parser.add_argument("--nSV", type=int, default=0, help="Number of background structural variants to simulate (default: 0)")
    TErandom_parser.add_argument("--regions", type=Existing_File_Path, 
                    help="Bed file with regions to include in the simulation.")
    TErandom_parser.add_argument("--exclude", type=Existing_File_Path, 
                    help="Bed file with regions to exclude from the simulation.")
    TErandom_parser.add_argument("--sense-strand-ratio", type=ratio, default=None, 
                    help="Proportion of TE variants in the sense strand (default: no strand preference)")
    # SNP and INDEL
    TErandom_parser.add_argument("--snp-rate", "-S", type=ratio, default=0.02, 
                    help="SNP mutation rate per base (default: 0.02)")
    TErandom_parser.add_argument("--indel-rate", "-I", type=ratio, default=0.005, 
                    help="Indel mutation rate per base (default: 0.005)")
    TErandom_parser.add_argument("--indel-ins", "-r", type=ratio, default=0.4, 
                    help="Proportion of insertion events among INDELs (0-1, default: 0.4)")
    TErandom_parser.add_argument("--indel-geom-p", "-G", type=ratio, default=0.7, 
                    help="Parameter 'p' of geometric distribution for indel lengths (default: 0.7)")
    # Truncation
    TErandom_parser.add_argument("--truncated-ratio", "-T", type=ratio, default=0.3, 
                    help="Proportion of TE sequences to truncate (0-1, default: 0.3)")
    TErandom_parser.add_argument("--truncated-max-length", "-K", type=ratio, default=0.5, 
                    help="Maximum proportion of sequence length to truncate (0-1, default: 0.5)")
    # PolyA
    TErandom_parser.add_argument("--polyA-ratio", "-A", type=ratio, default=0.8, 
                    help="Proportion of TE sequences to add polyA tail (0-1, default: 0.8)")
    TErandom_parser.add_argument("--polyA-min", "-M", type=int, default=5, 
                    help="Minimum polyA tail length (default: 5)")
    TErandom_parser.add_argument("--polyA-max", "-X", type=int, default=20, 
                    help="Maximum polyA tail length (default: 20)")
    # Other
    TErandom_parser.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    TErandom_parser.set_defaults(func=build_pool.run)

    # 2. TEreal
    TEreal_parser = subparsers.add_parser("TEreal", 
                               help="Generate pTE position from Known TE insertion and deletion")
    # Input
    TEreal_parser.add_argument("--knownINS", "-K", type=Existing_File_Path, required=True, 
                    help="Input known TE insertion file (e.g., MEI_Callset)")
    TEreal_parser.add_argument("--existingTEs", "-L", type=Existing_File_Path, required=True, 
                    help="Input known TE deletion file (RepeatMasker .out or UCSC .txt)")
    TEreal_parser.add_argument("--TEtype", "-e", type=str, action="append",
                    help="TEs to be extracted from the TE deletion file, with the default set as LINE, SINE, LTR, and Helitron.")
    TEreal_parser.add_argument("--CHR", "-C", type=str, required=True,
                    help="Chromosome to simulate TE insertions on (e.g., chr21 or 21)")
    TEreal_parser.add_argument("--DELlen", type=int, default=100,
                    help="A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)")
    TEreal_parser.add_argument("--nMIN", type=int, default=0,
                    help="A minimum number of TE deletions for each TE super family to be simulated (default: 0)")
    TEreal_parser.add_argument("--nSV", type=int, default=0, help="Number of background structural variants to simulate (default: 0)")
    # Output
    TEreal_parser.add_argument("--outprefix", "-O", type=File_Path, default="TEreal", 
                    help="Output prefix for generated BED file (default: 'TEreal')")
    TEreal_parser.add_argument("--nINS", "-N", type=int, default=0,
                    help="Number of TE insertions to simulate from --knownINS (default: 0)")
    TEreal_parser.add_argument("--nDEL", type=int, default=0,
                    help="Number of TE deletions to simulate from --existingTEs (default: 0)")
    TEreal_parser.add_argument("--nEXC", type=int, default=0,
                    help="Number of LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs) to simulate from --existingTEs; requires a RepeatMasker .out file (default: 0)")
    # Other
    TEreal_parser.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    TEreal_parser.set_defaults(func=TE_real.run)

    # 3. TE from pangenome graph
    TEpan_parser = subparsers.add_parser("TEpan", 
                               help="Generate pTE position from Pangenome graph")
    # Input
    TEpan_parser.add_argument("--gfa", "-G", type=Existing_File_Path, required=True, 
                    help="GFA file of the pangenome graph")
    TEpan_parser.add_argument("--lib", "-L", type=Existing_File_Path, required=True, 
                    help="RepeatMasker library file")
    TEpan_parser.add_argument("--CHR", "-H", type=str, required=True, 
                    help="Chromosome to simulate TE insertions on (e.g., chr21 or 21)")
    TEpan_parser.add_argument("--minLen", "-I", type=int, default= 250,
                    help="Minimum length of structural variants to consider (default: 250)")
    TEpan_parser.add_argument("--TEtype", "-e", type=str, action="append",
                    help="TEs to be extracted from the RepeatMasker annotation file, with the default set as LINE, SINE, LTR, and Helitron.")
    TEpan_parser.add_argument("--cov", "-C", type=ratio, default= 0.5,
                    help="Minimum TE coverage to consider a structural variant as TE (0-1, default: 0.5)")
    TEpan_parser.add_argument("--tmpDir", "-T", type=str, default="tmp_TEpan", 
                    help="Temporary directory for intermediate files (default: tmp_TEpan)")
    TEpan_parser.add_argument("--nINS", "-N", type=int, default=0,
                    help="Number of TE insertions to simulate from the pangenome graph (default: 0)")
    TEpan_parser.add_argument("--nDEL", type=int, default=0,
                    help="Number of TE deletions to simulate from the pangenome graph (default: 0)")
    TEpan_parser.add_argument("--nEXC", type=int, default=0,
                    help="LTR-LTR recombination (--nEXC) is not yet supported for TEpan (must be 0)")
    # Output
    TEpan_parser.add_argument("--outprefix", "-O", type=File_Path, default="TEpan", 
                    help="Prefix for output files (VCF + modified genome FASTA)")
    TEpan_parser.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    TEpan_parser.set_defaults(func=TEpan.run)


    # 4. simulate
    Simulate_parser = subparsers.add_parser("Simulate", 
                               help="Simulate TE insertions/deletions and generate VCF and modified genome FASTA")
    # Input
    Simulate_parser.add_argument("--ref", "-F", type=Fasta_File_Path, required=True, 
                    help="Reference genome FASTA file")
    Simulate_parser.add_argument("--pool", "-P", type=Fasta_File_Path, required=True, 
                    help="FASTA file of TE sequences generated from 'TEpool'")
    Simulate_parser.add_argument("--bed", "-B", type=Existing_File_Path, required=True, 
                    help="BED file containing TE positions (can be generated by 'TEreal')")
    # Output
    Simulate_parser.add_argument("--outprefix", "-O", type=File_Path, default="Sim", 
                    help="Prefix for output files (VCF + modified genome FASTA)")
    # Options
    Simulate_parser.add_argument("--num", "-N", type=int, required=True, 
                    help="Number of genomes to simulate")
    Simulate_parser.add_argument("--diverse", "-I", action="store_true",
                    help="Introduce sequence diversity among individuals for the same TE event")
    Simulate_parser.add_argument("--diverse_config", "-c", type=Existing_File_Path,
                    help="Path to a configuration file for introducing sequence diversity among individuals for the same TE event")
    Simulate_parser.add_argument("--af-min", "-A", type=ratio, default=0.1, 
                    help="Minimum allele frequency for simulated TE variants (default: 0.1)")
    Simulate_parser.add_argument("--af-max", "-X", type=ratio, default=0.9, 
                    help="Maximum allele frequency for simulated TE variants (default: 0.9)")
    Simulate_parser.add_argument("--af-dist", type=str,
                    choices=["uniform", "exponential", "reflected_exponential"],
                    default="uniform",
                    help="Shape of the per-event allele frequency draw. 'uniform' draws flat "
                         "between --af-min and --af-max. 'exponential' draws from an exponential "
                         "of mean --af-mean truncated to [--af-min, --af-max], giving a "
                         "population dominated by rare, largely private insertions with a thin "
                         "tail of common ones (default: uniform)")
    Simulate_parser.add_argument("--af-mean", type=ratio, default=None,
                    help="Mean of the exponential drawn by --af-dist exponential; required by it "
                         "and rejected without it. Truncating below at --af-min shifts the "
                         "distribution, so the realised mean allele frequency is roughly "
                         "--af-min plus --af-mean -- keep --af-min near 1/--num, or 0, unless "
                         "you mean to raise the floor")
    Simulate_parser.add_argument("--del-af-dist", type=str,
                    choices=["uniform", "exponential", "reflected_exponential"], default=None,
                    help="Allele frequency distribution for DELETIONS only, defaulting to "
                         "--af-dist. Deletions need their own because the two event types "
                         "measure opposite things: for an insertion GT=1 means the genome HAS "
                         "the element, for a deletion it means the genome LACKS the "
                         "reference's, so a deletion's allele frequency is the frequency of an "
                         "ABSENCE. A reference element most genomes lack -- a reference-private "
                         "insertion seen from the other side -- is a deletion at HIGH "
                         "frequency, which is what 'reflected_exponential' gives: the same "
                         "exponential mirrored within [min, max] so the mass piles against the "
                         "top and the mean sits at max minus --del-af-mean. Excisions stay on "
                         "--af-dist, being derived events that spread from rare like insertions")
    Simulate_parser.add_argument("--del-af-mean", type=ratio, default=None,
                    help="Mean for --del-af-dist, defaulting to --af-mean. For "
                         "reflected_exponential it is the mean distance BELOW --del-af-max")
    Simulate_parser.add_argument("--del-af-min", type=ratio, default=None,
                    help="Lower truncation bound for deletions, defaulting to --af-min")
    Simulate_parser.add_argument("--del-af-max", type=ratio, default=None,
                    help="Upper truncation bound for deletions, defaulting to --af-max")
    Simulate_parser.add_argument("--allow-zero-carriers", action="store_true",
                    help="Keep events that no genome drew. By default an event whose genotype "
                         "roll gives no carrier is given one, a genome chosen uniformly at "
                         "random, so that every event in the BED reaches at least one genome; "
                         "without it a low --af-mean quietly discards most of the pool instead "
                         "of producing private insertions. The realised minimum allele "
                         "frequency is then 1/--num")
    Simulate_parser.add_argument("--tsd-min", "-M", type=int, default=5, 
                    help="Minimum TSD length (default: 5)")
    Simulate_parser.add_argument("--tsd-max", "-Y", type=int, default=20, 
                    help="Maximum TSD length (default: 20)")
    Simulate_parser.add_argument("--tsd-from-header", action="store_true",
                    help="Take each insertion's TSD length from a TSD= tag in its pool FASTA "
                         "header, e.g. '>copia#LTR/Copia_3SNP TSD=5' or '... TSD=5-15' (both "
                         "bounds inclusive; TSD=0 means no duplication). The tag comes from the "
                         "--consensus library and is carried onto the pool by TErandom. "
                         "Elements with no tag fall back to --tsd-min/--tsd-max. Applies to "
                         "insertions only: deletions and excisions duplicate nothing")
    Simulate_parser.add_argument("--sense-strand-ratio", "-S", type=ratio, default=0.5, 
                    help="Proportion of TE variants in the sense strand (default: 0.5)")
    # Other
    Simulate_parser.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    Simulate_parser.set_defaults(func=simulate.run)

    # 5. compare
    Compare_parser = subparsers.add_parser("Compare", help="Compare predicted VCF to ground truth VCF")
    # Input
    Compare_parser.add_argument("--truth", "-T", type=str, required=True, help="Ground truth VCF file")
    Compare_parser.add_argument("--pred", "-P", type=str, required=True, help="Predicted file to compare")
    Compare_parser.add_argument("--predType", "-p", type=str,  choices=["VCF", "BED"], default="VCF", 
                    help="Type of the predicted file (VCF or BED, default: VCF)")
    # Output
    Compare_parser.add_argument("--outprefix", "-O", type=str, required=True, help="Output matched TEs")
    # Options
    Compare_parser.add_argument("--truthID", "-I", type=str, required=True, help="Sample ID in the truth VCF")
    Compare_parser.add_argument("--predID", "-J", type=str, required=True, help="Sample ID in the predicted VCF")
    Compare_parser.add_argument("--TEtype", "-e", type=str, default=None,
                    help="TE family in truth VCF to consider in the comparison, e.g. TY1-FULL "
                         "(case-insensitive). This is the family, not the superfamily -- TY1, "
                         "TY2, TY4 and TY5 are all LTR/Copia. Default: all families")
    Compare_parser.add_argument("--INSonly", action="store_true",  help="Only compare insertions in truth VCF")
    Compare_parser.add_argument("--nHap", "-N", type=int, default=2,  help="Number of haplotypes in the genome (default: 2)")
    Compare_parser.add_argument("--max_dist", "-M", type=int, default=100,
                    help="Maximum allowed distance (bp) to consider two variants as matching")
    Compare_parser.add_argument("--gt_len_tol", type=int, default=50,
                    help="Maximum allowed difference (bp) in allele length to consider two "
                         "alleles as the same when comparing genotypes (default: 50)")
    Compare_parser.set_defaults(func=compare_vcf.run)

    # 6. Evaluate
    Evaluate_parser = subparsers.add_parser("Evaluate",
                               help="Score a prediction against every simulated event at once, "
                                    "aggregated over all genomes")
    # Input
    Evaluate_parser.add_argument("--truth", "-T", type=str, required=True, help="Ground truth VCF file")
    Evaluate_parser.add_argument("--pred", "-P", type=str, required=True, help="Predicted file to evaluate")
    Evaluate_parser.add_argument("--predType", "-p", type=str, choices=["VCF", "BED"], default="VCF",
                    help="Type of the predicted file (VCF or BED, default: VCF). A BED "
                         "prediction carries no genotypes, so only detection and breakpoint "
                         "accuracy are reported")
    # Output
    Evaluate_parser.add_argument("--outprefix", "-O", type=str, default="evaluate",
                    help="Output prefix for the per-event JSON and the HTML report "
                         "(default: evaluate)")
    Evaluate_parser.add_argument("--no_html", action="store_true",
                    help="Skip the HTML report; write only the text summary, the JSON and "
                         "the per-locus files")
    # Options
    Evaluate_parser.add_argument("--sample_map", "-S", type=Existing_File_Path, default=None,
                    help="Two-column file pairing truth and prediction samples "
                         "(truth_sample<TAB>pred_sample). Default: pair samples by name")
    Evaluate_parser.add_argument("--TEtype", "-e", type=str, default=None,
                    help="TE family in truth VCF to consider, e.g. TY1-FULL (case-insensitive). "
                         "This is the family, not the superfamily -- TY1, TY2, TY4 and TY5 are "
                         "all LTR/Copia. Default: all families")
    Evaluate_parser.add_argument("--INSonly", action="store_true", help="Only evaluate insertions in truth VCF")
    Evaluate_parser.add_argument("--nHap", "-N", type=int, default=1,
                    help="Number of haplotype columns per individual in the truth VCF; >1 "
                         "merges them before pairing with the prediction (default: 1)")
    Evaluate_parser.add_argument("--max_dist", "-M", type=int, default=100,
                    help="Maximum allowed distance (bp) to consider a prediction and a "
                         "simulated event the same locus (default: 100)")
    Evaluate_parser.add_argument("--gt_len_tol", type=int, default=50,
                    help="Maximum allowed difference (bp) in allele length to consider two "
                         "alleles as the same (default: 50)")
    Evaluate_parser.add_argument("--reference", "-R", type=Fasta_File_Path, default=None,
                    help="Reference FASTA the truth and prediction coordinates refer to (indexed, or in a "
                         "writable directory). Adds a sequence check: each correctly placed call's allele is "
                         "written into the reference and compared base-for-base (case-insensitive) with the "
                         "simulated haplotype, so an anchor shifted along a TSD still matches. The usual report "
                         "is followed by a second, sequence-exact one in which only identical calls count; it is "
                         "also written to <outprefix>.exact.json and .exact.html.")
    Evaluate_parser.add_argument("--carrier_info", type=str, default=None, metavar="KEY",
                    help="INFO field of the prediction that marks which ALTs are carrier alleles. With one "
                         "value per allele (REF first, e.g. miniME's ME_INFO, or ALTs only), an ALT whose value "
                         "is '.' is not a carrier allele and its genomes are noncarriers. On a single-ALT record "
                         "any other field marks the ALT by being present and not '.' (e.g. GraffiTE's "
                         "repeat_ids). Without it, a multi-allelic call's carriers are inferred from allele "
                         "lengths and the simulated neighbours each genome carries")
    Evaluate_parser.add_argument("--size_bins", type=int, action="append", default=None,
                    help="Upper edges (bp) of the event-size strata; repeat the flag per edge "
                         f"(default: {' '.join(map(str, evaluate.DEFAULT_SIZE_BINS))})")
    Evaluate_parser.add_argument("--af_bins", type=ratio, action="append", default=None,
                    help="Upper edges of the allele-frequency strata; repeat the flag per edge "
                         f"(default: {' '.join(map(str, evaluate.DEFAULT_AF_BINS))})")
    Evaluate_parser.set_defaults(func=evaluate.run)

    # 7. Read simulation
    Readsim_parser = subparsers.add_parser("Readsim", 
                               help="generate short or long reads from the simulated genome")
    # general
    Readsim_parser.add_argument("--type", "-T", choices=["short", "long"],
                    type=str, required=True, help="Simulate short reads or long reads")
    Readsim_parser.add_argument("--genome", "-G", type=str, required=True, 
                    help="The file contains genomes where reads simulated from")
    Readsim_parser.add_argument("--depth", "-P", type=int, required=True, 
                    help="Depth of simulated reads")
    #Readsim_parser.add_argument("--ngenome", "-N", type=int, default= 1, 
    #                help="Number of genomes for simultaneous reads simulation")
    #Readsim_parser.add_argument("--outprefix", "-O", type=str, required=True, 
    #                help="prefix of output files")    
    # long reads settings
    Readsim_parser.add_argument("--Lerror", "-E", type=ratio, default= 0.15, 
                    help="sequencing error rate for long reads (default: 0.15)")
    Readsim_parser.add_argument("--Lmean", "-M", type=int,  default= 9000,
                    help="average size of read length (only for long reads, default: 9000 bp)")
    Readsim_parser.add_argument("--Lstd", "-S", type=int, default= 7000,  
                    help="read length standard deviation (only for long reads, default: 7000 bp)")
    # short reads settings
    Readsim_parser.add_argument("--length", "-i", type=int, default= 150, 
                    help="read length (only for short reads, default: 150 bp)")
    Readsim_parser.add_argument("--Fmean", "-m", type=int, default= 500, 
                    help="average size of fragment length (only for short reads, default: 500 bp)")
    Readsim_parser.add_argument("--Fstd", "-s", type=int, default= 30,  
                    help="Fragment size standard deviation (only for short reads, default: 30 bp)")
    # random seed
    Readsim_parser.add_argument("--seed", "-D", type=int, default=None, 
                    help="Random seed for reproducibility (default: None)")
    Readsim_parser.set_defaults(func=reads.run)



    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    if args.command in ("TErandom", "TEreal", "TEpan"):
        if args.nINS < 0 or args.nDEL < 0 or args.nEXC < 0:
            parser.error("--nINS, --nDEL, and --nEXC must be non-negative.")
        if args.nINS + args.nDEL + args.nEXC == 0:
            parser.error("Nothing to simulate: set at least one of --nINS, --nDEL, --nEXC to a value > 0.")
    if args.command == "TErandom":
        if (args.nDEL > 0 or args.nEXC > 0) and not args.existingTEs:
            TErandom_parser.error("--existingTEs is required when --nDEL or --nEXC > 0.")
        if args.nEXC > 0 and not str(args.existingTEs).lower().endswith(".out"):
            TErandom_parser.error("--nEXC requires a RepeatMasker .out --existingTEs file (LTR fragment structure is needed to identify full-length elements).")
    if args.command == "TEreal" and args.nEXC > 0 and not str(args.existingTEs).lower().endswith(".out"):
            TEreal_parser.error("--nEXC requires a RepeatMasker .out --existingTEs file (LTR fragment structure is needed to identify full-length elements).")
    if args.command == "TEpan" and args.nEXC > 0:
        TEpan_parser.error("LTR-LTR recombination (--nEXC) is not yet supported for TEpan.")
    if args.command == "Simulate":
        if args.diverse_config and not args.diverse:
            parser.error("--diverse_config requires --diverse to be set")
        # An inverted TSD pair dies inside np.random.randint without naming either flag,
        # unlike the allele frequency pair, which Simulator checks for itself.
        if args.tsd_min > args.tsd_max:
            Simulate_parser.error("--tsd-min must be less than or equal to --tsd-max")
        if args.tsd_min < 0:
            Simulate_parser.error("--tsd-min must be non-negative")
        if args.af_dist in ("exponential", "reflected_exponential"):
            if args.af_mean is None or args.af_mean <= 0:
                Simulate_parser.error(f"--af-dist {args.af_dist} requires --af-mean > 0")
            # Truncating an exponential from below shifts it, so the realised mean cannot
            # fall below --af-min however small --af-mean is. Asking for one that does is
            # not an error, but it will not be honoured, so say so.
            if args.af_mean < args.af_min:
                print(f"[WARN] --af-mean {args.af_mean} is below --af-min {args.af_min}; the "
                      f"realised mean allele frequency cannot fall below --af-min. Lower "
                      f"--af-min (0 removes the floor entirely) to get the mean you asked for.",
                      file=sys.stderr)
        elif args.af_mean is not None:
            Simulate_parser.error("--af-mean only applies to an exponential --af-dist")
        # The deletion spec inherits part by part, so validate what it RESOLVES to rather
        # than what was typed: --del-af-dist exponential with no --del-af-mean is legal and
        # means "the insertion's mean", which only works if there is one.
        del_dist = args.del_af_dist if args.del_af_dist is not None else args.af_dist
        del_mean = args.del_af_mean if args.del_af_mean is not None else args.af_mean
        del_lo = args.del_af_min if args.del_af_min is not None else args.af_min
        del_hi = args.del_af_max if args.del_af_max is not None else args.af_max
        if del_dist in ("exponential", "reflected_exponential") and (del_mean is None or del_mean <= 0):
            Simulate_parser.error(f"--del-af-dist {del_dist} requires --del-af-mean > 0 "
                                  f"(or an --af-mean to inherit)")
        if del_lo > del_hi:
            Simulate_parser.error("--del-af-min must be less than or equal to --del-af-max")
        if args.del_af_mean is not None and del_dist == "uniform":
            Simulate_parser.error("--del-af-mean only applies to an exponential --del-af-dist")
    if args.command == "Evaluate":
        # Bin edges name the upper bound of each stratum, so out-of-order edges would
        # silently produce strata no event can ever fall into.
        for flag, edges in (("--size_bins", args.size_bins), ("--af_bins", args.af_bins)):
            if edges is not None and list(edges) != sorted(set(edges)):
                Evaluate_parser.error(f"{flag} must be given in strictly increasing order")
        if args.nHap < 1:
            Evaluate_parser.error("--nHap must be at least 1")
    args.func(args)


if __name__ == "__main__":
    main()
