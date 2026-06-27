#!/usr/bin/env ruby
# encoding: utf-8

require 'sushi_fabric'
require_relative 'global_variables'
include GlobalVariables


class DEAApp < SushiFabric::SushiApp
  def initialize
    super
    @name = 'DEA'
    @analysis_category = 'Proteomics'
    @description =<<-EOS
prolfqua differential expression analysis (DEA) on a DIANN quantification.<br/>
Reuses the <a href="https://gitlab.bfabric.org/proteomics/slurmworker">A414_DEA</a>
prolfqua logic; runs <code>prolfqua_dea.sh</code> in the prolfquapp container
(apptainer on a genomics node).<br/>
Input is the DIANN <b>per-sample (grandchild)</b> dataset, which carries each
sample's [Factor] columns plus the shared <code>DIANN Quant [File]</code>
directory. Pick a grouping factor and the two groups to compare; alternatively
supply CONTROL / Contrast columns on the dataset for factorial / custom designs.
    EOS
    @params['process_mode'] = 'DATASET'
    # Gate to the DIANN per-sample dataset (it also carries the [Factor] columns
    # used by the grouping/sampleGroup/refGroup selectors).
    @required_columns = ['Name', 'DIANN Quant [File]']
    @required_params  = ['grouping', 'sampleGroup', 'refGroup']

    # ---- Compute + identity ----
    @params['cores']   = '1'  ; @params['cores',   'context'] = 'slurm'
    @params['ram']     = '16' ; @params['ram',     'context'] = 'slurm'
    @params['scratch'] = '20' ; @params['scratch', 'context'] = 'slurm'
    @params['node']    = ['fgcz-c-050']
    @params['node', 'context'] = 'slurm'
    @params['name']    = 'DEA_Result'
    @params['mail']    = ''

    # ---- Contrast definition (Limma-style; single contrast) ----
    # grouping names the [Factor] column to model; sampleGroup/refGroup pick the
    # two levels. The ezpyz_dea glue synthesizes the CONTROL column from these.
    @params['grouping']    = ''
    @params['grouping', 'description'] = 'name of the [Factor] column to group by (blank = first [Factor])'
    @params['sampleGroup'] = ''
    @params['sampleGroup', 'description'] = 'treatment level; must differ from refGroup'
    @params['refGroup']    = ''
    @params['refGroup', 'description'] = 'reference/control level; must differ from sampleGroup'

    # ---- prolfqua DEA knobs (mirror the A414_DEA executable; readable names,
    #      the glue re-prefixes to the 03|..09| keys prolfqua reads) ----
    @params['model']                = ['lm_impute', 'lm', 'saint', 'Extra']
    @params['model', 'description'] = 'primary DEA model; Extra defers to model_extra'
    @params['model_extra']          = ['rfit', 'firth', 'firth_nested', 'ropeca_nested']
    @params['model_extra', 'description'] = 'used only when model = Extra'
    @params['Normalization']        = ['vsn', 'none', 'robscale']
    @params['Difference_threshold'] = ['1', '0.6', '2']
    @params['FDR_threshold']        = ['0.1', '0.05', '0.25']
    @params['REVpattern']           = ['NONE', '^REV_|^rev_', '^REV|^rev', 'REV_', 'rev_']
    @params['REVpattern', 'description'] = 'decoy regex; NONE disables'
    @params['CONpattern']           = ['NONE', '^CON|^zz', '^zz|^CON|^Cont_', '^CON', '^zz', '^Cont_']
    @params['CONpattern', 'description'] = 'contaminant regex; NONE disables'
    @params['PeptideLevel']         = ['false', 'true']
    @params['PeptideLevel', 'description'] = 'true selects the peptide-level prolfqua adapter'

    @modules = ['Dev/R', 'Dev/pixi']  # pixi runs the ezpyz_dea app; container is apptainer
    @inherit_columns = ['Order Id']
  end

  def next_dataset
    @comparison = "#{@params['sampleGroup']}--over--#{@params['refGroup']}"
    @params['comparison'] = @comparison
    @params['name'] = @comparison
    # ezpyz_dea writes the prolfqua result into <comparison>/ in the job cwd, so
    # SUSHI's g-req copies it to gstore by basename (matches these [File] paths).
    report_dir = File.join(@result_dir, @comparison)
    {
      'Name'          => @comparison,
      'Report [Link]' => File.join(report_dir, 'index.html'),
      'Report [File]' => report_dir
    }.merge(extract_columns(colnames: @inherit_columns))
  end

  def commands
    # Python entry point: from ezpyz_dea.app import EzAppDea
    # (run_PyApp resolves 'DEA' -> module ezpyz_dea + class EzAppDea)
    run_PyApp('DEA', pixi_enabled: true)
  end
end

if __FILE__ == $0

end
