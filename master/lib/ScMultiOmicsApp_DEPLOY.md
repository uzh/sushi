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
- [ ] ScSeuratCombine: pick two ScMultiOmics output rows in SUSHI - the `SC Seurat` column
      must resolve to `<report>/scData.qs2` (a symlink to `scMultiData.qs2`); the combine
      job must load both via `qs2::qs_read` and `JoinLayers(..., assay = "RNA")` cleanly.

## Gotchas

These are things that bit us during the v7 → v9 iteration on `p31662 P8_AT_PBMCs`. Capture
here so the next person doesn't relearn them.

### exploreSC silent crash with empty error

If the `?data=` URL points to a file that does not exist at any `urlDataRoot`
(`/srv/gstore/projects` or `/srv/GT/analysis/course_sushi/public/gstore/projects`),
`execute_import()` triggers `show_error_modal("File Not Found")` + `stopApp()` and returns
`NULL`. `inputDataReactive` becomes `NULL`, then `server-umap.R` runs
`req(!is.null(NULL[["scData"]]))` → `req(FALSE)` → throws a `shiny.silent.error` whose
`$message` is `""`. The outer `tryCatch` in `app.R` then logged
`Error loading server components:` (with nothing after the colon), and the user only saw
a stalled spinner / "Connection reset".

**Fix landed in shinyproxy_apps `613b407`:**
- `parse_url_and_setup_dirs()` now logs `File not found. Tried: <both candidate paths>`.
- `execute_import()` has an explicit `is.na(dataDir)` guard.
- The outer `tryCatch` skips `shiny.silent.error` instead of dressing it up as a real
  error.

**Diagnose:** session stderr lives at `/srv/shinyproxy/logs/exploreSC_<proxy-id>_<date>_stderr.log`.
Find the proxy-id from the address bar (the cookie / URL), or grep `/srv/shinyproxy/logs/shinyproxy.log`.

### `fgcz-shiny.uzh.ch` is on fgcz-c-051

Distinct from `shiny.fgcz.uzh.ch` (private, port 8080) and `shiny-public.fgcz.uzh.ch`
(public, port 8081), which are documented in the `shinyproxy-deployment` skill. The
`fgcz-shiny.uzh.ch` reverse proxy resolves to fgcz-c-051 — same server, different port
and config (`/srv/shinyproxy/private.yml`).

### ADT DotPlot: never `scale_y_discrete(limits = ...)` with `coord_flip()`

`Seurat::DotPlot()` puts features on the data-x axis. After `coord_flip()` the visual
y-axis shows the features, but `scale_y_discrete()` always refers to the **data** y-axis
(the cluster ident). Setting `limits = rev(adt_features_ordered)` therefore tries to
restrict the cluster axis to feature names — every data point gets filtered out and the
plot renders empty (PNG generated, no circles).

**Right way to control the order:** pass `features =` in the desired pre-flip order
and stop there. Empirically (verified by extracting the rendered PNG from v9, then v10):
with `Seurat::DotPlot(...) + coord_flip()`, `features[1]` lands at the **bottom** of the
y-axis and `features[N]` at the top — opposite of what you'd guess. So to render
real markers alphabetical at the top with isotype controls grouped at the bottom:

```r
is_iso <- grepl("(?i)isotype", adt_features)
# Build the desired top-to-bottom order, then rev() it, because coord_flip() inverts.
adt_features_top_to_bottom <- c(gtools::mixedsort(adt_features[!is_iso]),
                                gtools::mixedsort(adt_features[is_iso]))
DotPlot(obj, features = rev(adt_features_top_to_bottom), assay = "ADT") + coord_flip()
```

Bug history: v9 passed `c(real, iso)` directly and got reverse-alphabetical real markers
on top with isotypes still above them — the rendered output looked random because the
plot inverts your input. v10 fixed it with the explicit `rev()` and renamed the
intermediate variable so the orientation is obvious from the code.

### `ezInteractiveTableRmd` in `results='asis'` chunks loses its dependencies

`print(htmltools::tagList(ezInteractiveTableRmd(mk, ...)))` inside an `asis` chunk emits
the widget HTML but does not register its DT JS/CSS htmlwidget deps. Symptom: the
rendered HTML has only ~2 occurrences of `DataTable` (just the title strings) instead of
the ~20+ a working interactive table produces, and the table is invisible / stripped.

**Fix:** put each modality's dotplot + table in a child Rmd
(`_scMultiOmics_wnn_deg_subtab.Rmd`) and call it via `knitr::knit_child()` from a small
asis loop in the parent. The child's chunks are non-asis, so the htmlwidget's deps load
normally.

### PBMC `scData.qs2` rarely has a celltype column

`ezRun::pickCellTypeColumn()` checks ~25 candidate names (`predicted.celltype.l2`,
`CyteType`, `azimuth_pan_human`, etc.). Standard ScSeurat output has none, so the
"by cell type" panels render the placeholder "No cell-type annotation column on the
object". For PBMC data, run `Azimuth::RunAzimuth(obj, reference = "pbmcref")` upstream
of render to populate `predicted.celltype.l2` — `pickCellTypeColumn()` will find it
automatically and the celltype panels populate.

`Azimuth` and `SeuratData::pbmcref` are both pre-installed in `Dev/R/4.5.0`.

### scSeuratCombine compatibility = symlink + dataset column

`ScSeuratCombineApp.rb` reads `input$getColumn("SC Seurat")` and expects each path to
end at a `scData.qs2`. `ezMethodScMultiOmics` saves the multi-modal object as
`scMultiData.qs2`, so the combine job can't find anything by default.

**Fix landed in this branch:**
- `app-ScMultiOmics.R` creates `scData.qs2 → scMultiData.qs2` (relative symlink)
  in the report dir after `makeRmdReport` returns.
- `ScMultiOmicsApp.rb#next_dataset` adds `'SC Seurat [Link]' => File.join(report_file,
  "scData.qs2")` so the output dataset is consumable by ScSeuratCombine /
  ScSeuratCombinedLabelClusters without any user intervention.

### `g-req copynow` is queued

The CLI prints `successfully copied` immediately, but the file lands in
`/srv/gstore/projects/...` ~15-30 s later (g-req drops the request into a queue
processed by `gstore-list`). Don't `ls` for it right after the copy command — wait or
poll. For overwriting an existing file, use `g-req copynow -f` (force); `g-req -w copy`
to a destination that already exists fails with "Destination path already exists".

### Tiny smoke objects (`P8_AT_PBMCs` = 2000 cells × 5000 genes)

The smoke fixtures in `/srv/GT/analysis/p31662/Analyses_Paul/scmultiomics_smoke/` are
deliberately downsampled. Tabs that gate on `nrow(obj[[assay]]) >= 2L` may silently
skip on these. When something looks "missing" in a smoke render, first check the
fixture dimensions before changing the template.

### exploreSC link in the Rmd needs the full gstore-relative path

The Rmd's setup chunk used to derive the link via:

```r
rp <- output$getColumn("Report")
if (!grepl("^p\\d+/", rp)) {
  proj <- regmatches(getwd(), regexec("/(p\\d+)/", getwd()))[[1]][2]
  rp   <- file.path(proj, basename(rp))     # BUG: drops intermediate dirs
}
```

This collapses any path with intermediate directories (e.g. `p31662/Analyses_Paul/<run>/<report_dir>`)
down to `p31662/<report_dir>/scMultiData.qs2`, which 404s in exploreSC. Symptom:
the rendered HTML's "Single Cell Explorer" link sends users to a path under
`/srv/gstore/projects/p31662/<report_dir>/scMultiData.qs2` that does not exist;
exploreSC then logs `File not found. Tried: ...` (see the silent-crash gotcha
above) and the app refuses to load.

**Fix:** the Rmd now declares two render params and prefers them over the cwd
heuristic:

```yaml
params:
  gstore_data_path: null   # gstore-relative dir holding scMultiData.qs2
  explore_url:      null   # full URL override
```

```r
explore_url <- if (!is.null(params$explore_url) && nzchar(params$explore_url)) {
  params$explore_url
} else if (!is.null(params$gstore_data_path) && nzchar(params$gstore_data_path)) {
  paste0("https://fgcz-shiny.uzh.ch/app/exploreSC/?data=",
         sub("/+$", "", params$gstore_data_path), "/scMultiData.qs2")
} else if (!is.null(output) && "Report" %in% output$colNames) {
  # ...legacy cwd heuristic, only fires when caller didn't supply a path...
}
```

The render driver (`v8_render.R`-style script for ad-hoc renders, or
`ezMethodScMultiOmics` for SUSHI-driven runs) passes the eventual gstore
destination:

```r
rmarkdown::render(
  input  = "ScMultiOmics.Rmd",
  ...,
  params = list(gstore_data_path =
    "p31662/Analyses_Paul/scmultiomics_smoke_2026-04-28/P8_AT_PBMCs_ScMultiOmicsReport")
)
```

When patching this for an existing rendered HTML (e.g. a delivery already on
gstore), copy the file out, `sed` the wrong URL to the right one, then
`g-req copynow -f` it back into place — same recipe as any other gstore
overwrite.
