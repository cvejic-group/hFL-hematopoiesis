
## Files

* `get_metacell_wnn.py`: get metacells, similar idea as <https://github.com/Teichlab/skeletal_dev_atlas/blob/main/scenicplus/scanpy_helpers.py>

* `lmmPG.R`: use linear model to identify PG links


## PG detection

### HSC


```sh
# metacell
python HSC.00.get_mc_bySample.py

# PG
sbatch HSC.01.getPG.sh
# dev diff PG
sbatch HSC.02.getDevDiffPG.sh
# merge
cell=HSC
awk 'FNR>1 || NR==1' $cell/PG/lmmPG_{1..500}.tsv > $cell.lmmPG.tsv
awk 'FNR>1 || NR==1' $cell/PG_DevDiff/lmmPG_{1..500}.tsv > $cell.lmmPG_DevDiff.tsv
```


### Other cell types

```sh
# get metacell
python get_mc_byCelltype.py

# prep data for other cell types
nohup Rscript 01.prep_pg_data.R > 01.prep_pg_data.log &
# on lumi
for ct in GP Granulocyte MEMP-t MEMP MEP MEMP-Mast-Ery MEMP-Ery Early-Ery Late-Ery MEMP-MK MK MastP-t MastP Mast MDP Monocyte Kupffer cDC1 cDC2 pDC ASDC LMPP LP Cycling-LP PreProB ProB-1 ProB-2 Large-PreB Small-PreB IM-B NK ILCP T
do
  sbatch --job-name=$ct --output=logs/${ct}_%A_%a.out 02.submit_pg.sh $ct
done

# merge
for cell in GP Granulocyte MEMP-t MEMP MEP MEMP-Mast-Ery MEMP-Ery Early-Ery Late-Ery MEMP-MK MK MastP-t MastP Mast MDP Monocyte Kupffer cDC1 cDC2 pDC ASDC LMPP LP Cycling-LP PreProB ProB-1 ProB-2 Large-PreB Small-PreB IM-B NK ILCP T
do
awk 'FNR>1 || NR==1' $cell/PG/lmmPG_{1..100}.tsv > $cell/$cell.lmmPG.tsv
done
```


