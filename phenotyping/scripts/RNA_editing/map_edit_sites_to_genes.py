
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

def load_tss_from_gtf(gtf: Path) -> pd.DataFrame:
    """Exactly your existing load_tss(), but returning a DataFrame with
    columns ['#chr','tss','gene_id'] so we can do closest‐TSS lookup."""
    anno = read_gtf(gtf)
    if type(anno).__module__ == 'polars.dataframe.frame':
        anno = anno.to_pandas()
    genes = anno.loc[anno['feature']=='gene', ['seqname','start','end','strand','gene_id']]
    # TSS coordinate:
    genes['tss'] = np.where(genes['strand']=='+', genes['start'], genes['end'])
    genes['#chr'] = genes['seqname']
    return genes[['#chr','tss','gene_id']]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("edit_sites_bed", help="BED of all editing sites (chr, start, end, site_id)")
    p.add_argument("ref_gtf",       help="GTF file for building TSS map")
    p.add_argument("output_tsv",    help="site<TAB>gene_id output (closest TSS)")
    return p.parse_args()

def main():
    args = parse_args()

    # 1) load sites
    sites = pd.read_csv(args.edit_sites_bed, sep="\t",
                        names=["#chr","start","end","site_id"], usecols=[0,1,2,3])

    # 2) load TSS
    tss = load_tss_from_gtf(Path(args.ref_gtf))

    # 3) merge & pick closest
    merged = sites.merge(tss, on="#chr", how="inner")
    merged["dist"] = (merged["start"] - merged["tss"]).abs()
    best = (
        merged
        .sort_values(["site_id","dist"])
        .drop_duplicates("site_id", keep="first")
        [["site_id","gene_id"]]
    )
    # 4) write
    best.to_csv(args.output_tsv, sep="\t", index=False, header=False)

if __name__=="__main__":
    main()
