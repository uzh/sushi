# ScMultiOmics SUSHI app deployment

## Prerequisite

Install `~/git/ezRun-scmultiomics` (branch `scmultiomics`) into the user library:

```bash
module load Dev/R/4.5.0
Rscript -e 'devtools::install("~/git/ezRun-scmultiomics", upgrade = FALSE)'
```

For a permanent install, merge the `scmultiomics` branch of `~/git/ezRun` to `master`
once Phase 1-4 changes are reviewed.

## Test instance (fgcz-h-083)

```bash
ssh fgcz-h-083
cp ~/git/sushi/master/lib/ScMultiOmicsApp.rb \
   /srv/sushi/<your_test_instance>/master/lib/
touch /srv/sushi/<your_test_instance>/master/tmp/restart.txt
```

Open the test SUSHI URL → verify "ScMultiOmics" appears under SingleCell.

## Production (fgcz-h-082)

```bash
scp ~/git/sushi/master/lib/ScMultiOmicsApp.rb fgcz-h-082:/tmp/
ssh trxcopy@fgcz-h-082 \
  'cp /tmp/ScMultiOmicsApp.rb /srv/sushi/production/master/lib/ && \
   touch /srv/sushi/production/master/tmp/restart.txt'
```

## Course + Demo (fgcz-h-081)

Per CLAUDE.md: course + demo SUSHI both live on fgcz-h-081 — apply to both:

```bash
scp ~/git/sushi/master/lib/ScMultiOmicsApp.rb fgcz-h-081:/tmp/
ssh trxcopy@fgcz-h-081 \
  'for inst in course_sushi demo_sushi; do
     cp /tmp/ScMultiOmicsApp.rb /srv/sushi/$inst/master/lib/
     touch /srv/sushi/$inst/master/tmp/restart.txt
   done'
```

## Submission flow

1. **CellRanger Multi or ARC input**: select the original CellRangerMulti / CellRangerARCCount
   dataset; required columns are `Name`, `Species`, `refBuild`, `CountMatrix`.
   Optional: `Report` / `SC Cluster Report` linking to a previous ScSeurat output
   (so we can reuse the annotated `scData.qs2`). Without that link, ezMethodScMultiOmics
   will refuse to run unless `SCDataOrigin = BDRhapsody`.
2. **BD Rhapsody input**: SCDataOrigin column must equal `BDRhapsody`. The app reads
   the pre-built `*_Seurat.rds` from `BDRhapsodyPath` (or the parent of `CountMatrix`).
3. **Standalone VDJ**: add `VDJTPath` and/or `VDJBPath` columns pointing at the
   directory containing `filtered_contig_annotations.csv`; auto-discovery is overridden.

## Post-deployment smoke checklist

- [ ] CellRanger Multi: re-run on p31662 P8_AT_PBMCs - expect Overview + ADT + VDJ + WNN tabs.
- [ ] CellRanger ARC: re-run on p39132 IC104_Plate_6217_1_ARC - expect Overview + ATAC + WNN.
- [ ] BD Rhapsody: re-run on p39179 394581_1-CD34run1 - expect Overview + ADT.
- [ ] exploreSC: load `scMultiData.qs2` for each - RNA UMAP must render unchanged.
