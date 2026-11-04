#!/usr/bin/env python3
"""Convert each deduplicated disomy HiCUP BAM to a raw-count Juicer contact map."""
import argparse,json,subprocess
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("config");a=ap.parse_args()
    # Read the work directory and resource settings from the configuration.
    c=json.loads(Path(a.config).read_text());work=Path(c["work_dir"]).expanduser().resolve()
    def run(cmd):print(cmd,flush=True);subprocess.run(list(map(str,cmd)),check=True)
    # Convert each disomy replicate independently; do not pool the donors.
    for sample in ("16216","17759"):
        # Step 01 produces one BAM after duplicate removal across all its lanes.
        bams=list((work/"hicup"/sample/"output").glob("*.hicup.bam"))
        if len(bams)!=1:raise ValueError(f"Expected one library-deduplicated HiCUP BAM for {sample}; found {len(bams)}")
        bam=bams[0];dest=work/"maps"/sample
        if dest.exists():raise FileExistsError(f"Refusing to overwrite {dest}")
        dest.mkdir(parents=True);run(["samtools","quickcheck","-v",bam])
        # Name sorting preserves HiCUP's validated contact set; no extra MAPQ/duplicate filter.
        sorted_bam=dest/f"{sample}.name_sorted.hicup.bam"
        run(["samtools","sort","-n","-@",c.get("threads",16),"-o",sorted_bam,bam])
        # Convert the paired BAM alignments to Juicer pre-format contact records.
        run(["hicup2juicer",sorted_bam]);pre=Path(str(sorted_bam)+".prejuicer")
        if not pre.is_file() or not pre.stat().st_size:raise ValueError(f"Missing converter output: {pre}")
        # Build the raw-count hg38 map that HiC-DC+ will read in step 03.
        hic=dest/f"{sample}.hic"
        run(["juicer_tools",f"-Xmx{c.get('java_heap','48g')}","pre","--threads",c.get("threads",16),pre,hic,"hg38"])
        if not hic.is_file() or not hic.stat().st_size:raise ValueError(f"Missing map: {hic}")
if __name__=="__main__":main()
