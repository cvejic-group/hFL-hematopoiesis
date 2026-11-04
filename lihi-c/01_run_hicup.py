#!/usr/bin/env python3
"""Align/filter/deduplicate each biological library with HiCUP (Python standard library)."""
import argparse,csv,gzip,json,shutil,subprocess
from pathlib import Path

def main():
    # Read the configuration and the table listing paired FASTQs.
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config");ap.add_argument("samples");ap.add_argument("--check-only",action="store_true")
    a=ap.parse_args();c=json.loads(Path(a.config).read_text())
    work=Path(c["work_dir"]).expanduser().resolve(); fasta=Path(c["fasta"]).expanduser().resolve()
    if not fasta.is_file():raise ValueError(f"Missing reference FASTA: {fasta}")
    # Group lanes by donor; reject missing files and repeated inputs.
    rows=list(csv.DictReader(Path(a.samples).open(),delimiter="\t"));groups={}
    if not rows:raise ValueError("Sample sheet is empty")
    seen=set()
    for r in rows:
        sample=r["sample"]
        if sample not in {"16216","17759"}:raise ValueError(f"Disomy donors only: {sample}")
        pair=[Path(r[k]).expanduser().resolve() for k in ("r1","r2")]
        for f in pair:
            if not f.is_file() or not f.stat().st_size:raise ValueError(f"Missing/empty FASTQ: {f}")
            if f in seen:raise ValueError(f"Repeated FASTQ: {f}")
            seen.add(f)
        groups.setdefault(sample,[]).append(pair)
    if set(groups)!={"16216","17759"}:raise ValueError("Provide both disomy libraries")
    # --check-only validates paths without starting alignment or writing outputs.
    if a.check_only:
        print("Validated FASTA and manifest:", {k:len(v) for k,v in groups.items()});return
    # Print external commands and stop immediately if a tool fails.
    def run(cmd):
        print("Running:",cmd,flush=True);subprocess.run([str(x) for x in cmd],check=True)
    # Build the hg38 alignment index once, unless an existing index was supplied.
    refs=work/"references";refs.mkdir(parents=True,exist_ok=True)
    index=Path(c.get("bowtie2_index",str(refs/"hg38"))).expanduser().resolve()
    if not any(Path(str(index)+suffix).exists() for suffix in (".1.bt2",".1.bt2l")):
        index.parent.mkdir(parents=True,exist_ok=True);run(["bowtie2-build",fasta,index])
    # Use a supplied HindIII digest, or create one from this reference FASTA.
    if c.get("digest"):
        digest=Path(c["digest"]).expanduser().resolve()
    else:
        # Reference-specific empty directory prevents selecting an unrelated digest.
        digest_dir=refs/"HindIII_digest";digest_dir.mkdir(exist_ok=True)
        candidates=list(digest_dir.glob("Digest_*.txt"))
        if not candidates:
            subprocess.run(["hicup_digester","--genome","hg38","--re1","A^AGCTT,HindIII",str(fasta)],cwd=digest_dir,check=True)
            candidates=list(digest_dir.glob("Digest_*.txt"))
        if len(candidates)!=1:raise ValueError("Expected exactly one HindIII digest; set digest explicitly")
        digest=candidates[0]
    if not digest.is_file():raise ValueError(f"Missing digest: {digest}")
    # Process the two biological libraries separately. Existing runs are protected.
    for sample,pairs in groups.items():
        run_dir=work/"hicup"/sample
        if run_dir.exists():raise FileExistsError(f"Refusing to overwrite existing HiCUP run: {run_dir}")
        out=run_dir/"output";out.mkdir(parents=True)
        inputs=run_dir/"inputs";inputs.mkdir()
        # Concatenate lanes before HiCUP so duplicate removal is across the library.
        merged=[]
        for mate in range(2):
            dest=inputs/f"{sample}_R{mate+1}.fastq.gz"
            with dest.open("wb") as target:
                for pair in pairs:
                    src=pair[mate]
                    if src.name.endswith(".gz"):
                        with src.open("rb") as f:shutil.copyfileobj(f,target)
                    else:
                        with src.open("rb") as f, gzip.GzipFile(fileobj=target,mode="wb") as gz:shutil.copyfileobj(f,gz)
            merged.append(dest)
        # HiCUP aligns with Bowtie2, filters invalid ligation products, and removes duplicates.
        # Keep its BAM, original QC summaries, generated configuration, and command log.
        conf=run_dir/"hicup.conf"
        conf.write_text(f"Outdir: {out}\nThreads: {int(c.get('threads',16))}\nQuiet: 0\nKeep: 1\nZip: 1\nBowtie2: {shutil.which('bowtie2') or 'bowtie2'}\nIndex: {index}\nDigest: {digest}\nFormat: Sanger\n\n"+"\n".join(map(str,merged))+"\n")
        with (run_dir/"command.log").open("w") as log:
            subprocess.run(["hicup","--config",str(conf)],stdout=log,stderr=subprocess.STDOUT,check=True)
    print("HiCUP finished:",work/"hicup")
if __name__=="__main__":main()
