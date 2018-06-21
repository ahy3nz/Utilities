from msibi_utils.plot_fit import plot_all_fits
from msibi_utils.animate_rdf import animate_all_pairs_states

logfile = 'charmm_bi/6-12-18.log'
plot_all_fits(logfile,use_agg=True,ylims=(0,1))  
animate_all_pairs_states(logfile,"charmm_bi/rdfs", 
        potentials_dir="charmm_bi/potentials",
        rdf_dir="charmm_bi/rdfs", 
        potentials2_dir="charmm_bulk_debi_flhe_dspc/potentials",
        rdf2_dir="charmm_bulk_debi_flhe_dspc/rdfs", use_agg=True)
