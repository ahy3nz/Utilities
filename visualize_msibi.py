from msibi_utils.plot_fit import plot_all_fits
from msibi_utils.animate_rdf import animate_all_pairs_states

logfile = 'from_accre/12-4.log'
plot_all_fits(logfile,use_agg=True,ylims=(-1,1))  
animate_all_pairs_states(logfile,"from_accre/rdfs", potentials_dir="from_accre/potentials",
        rdf_dir="from_accre/rdfs", use_agg=True)
