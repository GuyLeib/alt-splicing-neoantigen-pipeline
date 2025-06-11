import gzip
import logging
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple, Any

import pandas as pd
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(funcName)s] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


DB = None
GENOME = None
EXPRESSION_DF = None


def initialize_global_resources(args: argparse.Namespace):
    """
    Loads and initializes global resources: genome, expression matrix, and GTF database.
    """
    global DB, GENOME, EXPRESSION_DF

    db_path = args.db_path
    if not db_path.exists():
        logging.info(f"GTF database not found at {db_path}. Creating it from {args.gtf_file}...")
        try:
            gffutils.create_db(
                str(args.gtf_file),
                dbfn=str(db_path),
                force=True,
                keep_order=True,
                merge_strategy="merge",
                sort_attribute_values=True,
                disable_infer_transcripts=True,
                disable_infer_genes=True
            )
            logging.info("Database created successfully.")
        except Exception as e:
            logging.error(f"Failed to create gffutils database: {e}", exc_info=True)
            raise
    DB = gffutils.FeatureDB(str(db_path), keep_order=True)

    logging.info("Loading reference genome into memory...")
    try:
        GENOME = SeqIO.to_dict(SeqIO.parse(args.reference_genome, "fasta"))
    except FileNotFoundError:
        logging.error(f"Reference genome file not found: {args.reference_genome}")
        raise
    logging.info("Reference genome loaded.")

    logging.info("Loading expression data...")
    try:
        EXPRESSION_DF = pd.read_csv(args.expression_file, sep='\t', compression='gzip', index_col='Name')
    except FileNotFoundError:
        logging.error(f"Expression file not found: {args.expression_file}")
        raise
    logging.info("Expression data loaded.")


# --- Core Functions ---

def get_expressed_genes_mask(df: pd.DataFrame, sample: str, min_expr: float, gene_col: str = 'gene_name') -> pd.Series:
    """
    Creates a boolean mask for a DataFrame based on gene expression levels.
    """
    if sample not in EXPRESSION_DF.columns:
        logging.warning(f"Sample '{sample}' not found in expression data. All events will be filtered out.")
        return pd.Series([False] * len(df), index=df.index)

    expression_map = EXPRESSION_DF[sample]
    gene_expression = df[gene_col].map(expression_map).fillna(0)

    return gene_expression >= min_expr


def extract_spladder_events(sample: str, event_type: str, data_dir: Path) -> pd.DataFrame:
    """
    Reads SplAdder output files, filtering for events with novel junctions using GFF.
    """
    confirmed_file = data_dir / f"{sample}_{event_type}_C2.confirmed.txt.gz"
    gff_file = data_dir / f"{sample}_{event_type}_C2.confirmed.gff3"

    if not confirmed_file.exists() or not gff_file.exists():
        return pd.DataFrame()

    try:
        gff_db = gffutils.create_db(str(gff_file), dbfn=":memory:", force=True, keep_order=True)
        novel_ids = {
            feature.attributes.get('ID', [None])[0]
            for feature in gff_db.all_features()
            if feature.attributes.get('HasNovelJunction', ['N'])[0] == 'Y'
        }
    except Exception as e:
        logging.warning(f"Could not parse GFF for {sample}, {event_type}: {e}")
        return pd.DataFrame()

    if not novel_ids:
        return pd.DataFrame()

    try:
        with gzip.open(confirmed_file, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
    except Exception as e:
        logging.warning(f"Could not read confirmed events file for {sample}, {event_type}: {e}")
        return pd.DataFrame()

    if 'event_id' not in df.columns:
        return pd.DataFrame()

    df_novel = df[df['event_id'].isin(novel_ids)].copy()
    if df_novel.empty:
        return df_novel

    df_novel['event_type'] = event_type
    df_novel['sample'] = sample
    return df_novel


def decide_event_usage(row: pd.Series, sample: str, psi_threshold: float) -> str:
    """
    Determines which isoform is predominantly used based on PSI value.
    """
    psi_col = f"{sample}:psi"
    psi_val = row.get(psi_col)

    if pd.isna(psi_val):
        return "unknown"

    usage_map = {
        'intron_retention': ('retained', 'spliced'),
        'exon_skip': ('inclusion', 'skipping'),
        'alt_3prime': ('alt_site', 'ref_site'),
        'alt_5prime': ('alt_site', 'ref_site'),
        'mutex_exons': ('exon2', 'exon3'),
    }

    event_type = row['event_type']
    alt_isoform, ref_isoform = usage_map.get(event_type, ('alt_isoform', 'ref_isoform'))

    return alt_isoform if psi_val >= psi_threshold else ref_isoform


def get_sequence(chrom: str, start: int, end: int, strand: str) -> str:
    """Extracts a genomic sequence from the pre-loaded FASTA dictionary."""
    if chrom not in GENOME:
        logging.warning(f"Chromosome '{chrom}' not found in reference genome. Skipping sequence.")
        return ""

    seq_record = GENOME[chrom]
    subseq = seq_record.seq[start - 1:end]

    return str(subseq.reverse_complement() if strand == '-' else subseq)


def translate_sequence(seq: str) -> str:
    """Translates a nucleotide sequence, stopping at the first stop codon."""
    try:
        protein = str(Seq(seq).translate(to_stop=True))
        return protein
    except Exception:
        return ""

def find_overlapping_cds_features(db: gffutils.FeatureDB, chrm: str, start: int, end: int, strand: str) -> List[Dict[str, Any]]:
    """
    Finds all CDS features within a given genomic region on a specific strand,
    and returns their overlapping coordinates and phase.
    """
    overlapping = []
    region_feats = db.region(seqid=chrm, start=start, end=end, strand=strand, completely_within=False)
    for feat in region_feats:
        if feat.featuretype == 'CDS':
            cds_start = max(start, feat.start)
            cds_end = min(end, feat.end)
            if cds_start <= cds_end:
                phase = int(feat.frame) if feat.frame is not None and feat.frame in ['0', '1', '2'] else 0
                overlapping.append({
                    'start': cds_start,
                    'end': cds_end,
                    'phase': phase,
                    'feature_id': feat.id
                })
    overlapping.sort(key=lambda x: x['start'])
    return overlapping

def build_codon_list(db: gffutils.FeatureDB, chrm: str, seg_start: int, seg_end: int, strand: str,
                     leftover_bases: str, current_frame: int) -> Tuple[List[str], str, int]:
    """
    Constructs a list of codons for a genomic segment, considering overlapping CDS features
    and maintaining the reading frame.
    """
    codons = []
    
    cds_feats = find_overlapping_cds_features(db, chrm, seg_start, seg_end, strand)

    segment_nucleotide_sequence = ""
    for subf in cds_feats:
        sub_seq_part = get_sequence(chrm, subf['start'], subf['end'], strand)
        segment_nucleotide_sequence += sub_seq_part
    
    combined_seq = leftover_bases + segment_nucleotide_sequence

    if current_frame >= len(combined_seq):
        return [], combined_seq, current_frame

    seq_to_translate = combined_seq[current_frame:]
    
    temp_leftover_bases = ""
    
    i = 0
    while i + 3 <= len(seq_to_translate):
        codon = seq_to_translate[i:i+3]
        codons.append(codon)
        i += 3
    
    temp_leftover_bases = seq_to_translate[i:]
    new_current_frame = len(temp_leftover_bases) % 3
    
    return codons, temp_leftover_bases, new_current_frame

def construct_spliced_protein(row: pd.Series) -> str:
    """
    Constructs the spliced CDS sequence and translates it, handling reading frames.
    """
    event_type = row['event_type']
    usage = row['usage']

    if usage == 'unknown':
        return ""

    chrom = row['chrm']
    strand = row['strand']

    segment_map = {
        'intron_retention': {
            'retained': [('e1_start', 'e1_end'), ('e2_start', 'e2_end'), ('e3_start', 'e3_end')],
            'spliced': [('e1_start', 'e1_end'), ('e3_start', 'e3_end')],
        },
        'exon_skip': {
            'inclusion': [('e1_start', 'e1_end'), ('e2_start', 'e2_end'), ('e3_start', 'e3_end')],
            'skipping': [('e1_start', 'e1_end'), ('e3_start', 'e3_end')],
        },
        'alt_3prime': {
            'alt_site': [('e1_start', 'e1_end'), ('e2_start', 'e2_end'), ('e3_start', 'e3_end')],
            'ref_site': [('e1_start', 'e1_end'), ('e3_start', 'e3_end')],
        },
        'alt_5prime': {
            'alt_site': [('e1_start', 'e1_end'), ('e2_start', 'e2_end'), ('e3_start', 'e3_end')],
            'ref_site': [('e1_start', 'e1_end'), ('e3_start', 'e3_end')],
        },
        'mutex_exons': {
            'exon2': [('e1_start', 'e1_end'), ('e2_start', 'e2_end'), ('e4_start', 'e4_end')],
            'exon3': [('e1_start', 'e1_end'), ('e3_start', 'e3_end'), ('e4_start', 'e4_end')],
        },
    }

    coord_pairs = segment_map.get(event_type, {}).get(usage)
    if not coord_pairs:
        logging.debug(f"No coordinate pairs found for event type {event_type} and usage {usage}")
        return ""

    segments_coords = []
    for start_col, end_col in coord_pairs:
        start = row.get(start_col)
        end = row.get(end_col)
        if pd.notna(start) and pd.notna(end):
            segments_coords.append((int(start), int(end)))

    segments_coords.sort(key=lambda p: p[0], reverse=(strand == '-'))
    
    all_codons = []
    leftover_bases = ""
    current_frame = 0 # Represents the number of bases needed to complete the current codon.

    # Establish initial reading frame from the gene's first CDS, if available.
    if segments_coords and DB is not None and 'gene_id' in row and pd.notna(row['gene_id']):
        try:
            gene_id = row['gene_id']
            first_gene_cds = next(DB.children(gene_id, featuretype='CDS', order_by='start', limit=1), None)
            if first_gene_cds and first_gene_cds.frame is not None and first_gene_cds.frame in ('0', '1', '2'):
                current_frame = int(first_gene_cds.frame)
        except gffutils.exceptions.FeatureNotFoundError:
            logging.debug(f"Gene ID {row['gene_id']} not found in annotation DB for initial frame.")
        
    for start, end in segments_coords:
        if DB is None:
            logging.error("GTF database (DB) is not initialized. Cannot construct protein sequence.")
            return ""
            
        segment_codons, new_leftover_bases, new_current_frame = \
            build_codon_list(DB, chrom, start, end, strand, leftover_bases, current_frame)
        
        all_codons.extend(segment_codons)
        leftover_bases = new_leftover_bases
        current_frame = new_current_frame

    final_nucleotide_sequence = "".join(all_codons) + leftover_bases
    
    return translate_sequence(final_nucleotide_sequence)


def process_sample(sample: str, args: argparse.Namespace) -> Path:
    """
    Runs the full processing pipeline for a single sample.
    """
    logging.info(f"Processing sample: {sample}")
    all_records = []
    output_dir = Path(args.output_dir)

    for event_type in args.event_types:
        df_evt = extract_spladder_events(sample, event_type, Path(args.data_dir))
        if df_evt.empty:
            continue
        logging.info(f"  - Found {len(df_evt)} novel '{event_type}' events.")

        expressed_mask = get_expressed_genes_mask(df_evt, sample, args.min_expression)
        df_expressed = df_evt[expressed_mask].copy()  
        if df_expressed.empty:
            logging.info(f"  - No events passed expression filter for '{event_type}'.")
            continue
        logging.info(f"  - {len(df_expressed)} events passed expression filter.")

        df_expressed['usage'] = df_expressed.apply(
            decide_event_usage, axis=1, sample=sample, psi_threshold=args.psi_threshold
        )

        df_expressed['protein_sequence'] = df_expressed.apply(construct_spliced_protein, axis=1)

        for _, row in df_expressed.iterrows():
            prot_seq = row['protein_sequence']
            if not prot_seq or len(prot_seq) < 8:
                continue

            record_id = f"{sample}|{row['event_id']}|{row['usage']}"
            psi_val = row.get(f'{sample}:psi', 'NA')
            psi_str = f"{psi_val:.2f}" if isinstance(psi_val, (int, float)) else "NA"
            description = (f"gene={row.get('gene_name', 'NA')} event_type={event_type} "
                           f"psi={psi_str}")

            all_records.append(SeqRecord(Seq(prot_seq), id=record_id, description=description))

    output_fasta_path = output_dir / f"{sample}_all_splicing_peptides.fasta"
    if all_records:
        with open(output_fasta_path, "w") as out_handle:
            SeqIO.write(all_records, out_handle, "fasta")
        logging.info(f"Wrote {len(all_records)} peptides for sample {sample} to {output_fasta_path}")
    else:
        logging.info(f"No valid peptides generated for sample: {sample}")

    return output_fasta_path


def main():
    """
    Main function to parse arguments and orchestrate the pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Process SplAdder output to find splicing-derived neoantigens.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # --- Input/Output Arguments ---
    parser.add_argument("-d", "--data-dir", type=Path, required=True, help="Directory containing SplAdder output files.")
    parser.add_argument("-o", "--output-dir", type=Path, required=True, help="Directory to save the output FASTA files.")
    parser.add_argument("-e", "--expression-file", type=Path, required=True, help="Path to the gzipped gene expression matrix (TSV format).")
    parser.add_argument("-g", "--gtf-file", type=Path, required=True, help="Path to the reference genome annotation GTF file.")
    parser.add_argument("-r", "--reference-genome", type=Path, required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--db-path", type=Path, default=Path("annotation.db"), help="Path to store/read the gffutils database.")

    # Pipeline Parameters 
    parser.add_argument("--min-expression", type=float, default=1.0, help="Minimum log10(TPM) expression value to consider a gene.")
    parser.add_argument("--psi-threshold", type=float, default=0.5, help="PSI threshold to determine dominant isoform usage.")
    parser.add_argument("--event-types", nargs='+', default=["exon_skip", "alt_3prime", "alt_5prime", "mutex_exons", "intron_retention"], help="List of SplAdder event types to process.")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel processes to use.")
    parser.add_argument("--sample", type=str, help="Process a single sample ID instead of all samples found in data-dir.")
    args = parser.parse_args()

    # --- Setup ---
    args.output_dir.mkdir(parents=True, exist_ok=True)
    initialize_global_resources(args)

    if args.sample:
        samples = [args.sample]
        logging.info(f"Running pipeline for a single sample: {args.sample}")
    else:
        samples = sorted(list(set(
            f.name.split('_')[0] for f in args.data_dir.glob('*_C2.confirmed.txt.gz')
            if f.name.split('_')[0] != 'merge'
        )))
        logging.info(f"Found {len(samples)} samples to process.")

    if not samples:
        logging.warning("No samples found to process. Exiting.")
        return

    # Execution 
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_sample, sample, args) for sample in samples]
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"A worker process failed: {e}", exc_info=True)

    logging.info("Pipeline finished.")


if __name__ == "__main__":
    main()