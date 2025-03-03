import argparse
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from plotnine import *
import pysam
from tqdm import tqdm
import random


class MartiStats:
    class NestedDict(defaultdict):
        def __call__(self):
            return MartiStats.NestedDict(self.default_factory)

    def __init__(self, bam_paths, experiment_names, max_reads_per_experiment=None, plot_sample_size=None):
    	self.bams = bam_paths
    	self.experiments = experiment_names
    	self.max_reads_per_experiment = max_reads_per_experiment
    	self.plot_sample_size = plot_sample_size
    	self.default_fig_size = (10, 10)
    	self.read_count = defaultdict(int)
    	self.class_profile_count = MartiStats.NestedDict(MartiStats.NestedDict(int))
    	self.class_count = MartiStats.NestedDict(MartiStats.NestedDict(int))
    	self.read_len = MartiStats.NestedDict(MartiStats.NestedDict(list))
    	self.structure_profile_count = MartiStats.NestedDict(MartiStats.NestedDict(MartiStats.NestedDict(int)))
    	self.adapter_lev = MartiStats.NestedDict(MartiStats.NestedDict(MartiStats.NestedDict(list)))
    	self.poly_lev = MartiStats.NestedDict(MartiStats.NestedDict(MartiStats.NestedDict(list)))
    	self.poly_len = MartiStats.NestedDict(MartiStats.NestedDict(list))
    	self.cdna_len = MartiStats.NestedDict(MartiStats.NestedDict(list))
    	self.start_adapter_offset = MartiStats.NestedDict(MartiStats.NestedDict(MartiStats.NestedDict(list)))
    	self.end_adapter_offset = MartiStats.NestedDict(MartiStats.NestedDict(MartiStats.NestedDict(list)))
    	self.collect_info()
        
    ########################
	#  BAM parsing / stats #
	########################
	#-------------------------------------------------------------------------------
    def collect_info(self):
        for exp_name, bam_name in zip(self.experiments, self.bams):
            input_bam = pysam.AlignmentFile(bam_name, "rb", check_sq=False, ignore_truncation=True)
            print("Loading ", bam_name)
            for read in tqdm(input_bam):
                if self.max_reads_per_experiment and self.read_count[exp_name] >= self.max_reads_per_experiment: break
                self.read_count[exp_name] += 1
                class_profile = read.get_tag("lb")
                self.class_profile_count[exp_name][class_profile] += 1
                for c in class_profile.split(","):
                	self.class_count[exp_name][c] += 1
                self.read_len[exp_name][class_profile].append(len(read.seq))
                if class_profile in ["TooShort", "LowRq"]: continue
                self.structure_profile_count[exp_name][class_profile][read.get_tag("st")] += 1
                structural_hits = read.get_tag("ch").split(',')
                for i, entry in enumerate(structural_hits):
                    name, start, end, lev = entry.split(":")
                    if any(s in name for s in ["adapterA", "adapterB"]):
                        self.adapter_lev[exp_name][class_profile][name].append(int(lev))
                        if i == 0:
                            self.start_adapter_offset[exp_name][class_profile][name].append(int(start))
                        elif i == len(structural_hits) - 1:
                            self.end_adapter_offset[exp_name][class_profile][name].append(len(read.seq) - int(end) - 1)
                    if name in ["polyA", "polyT"]:
                        self.poly_lev[exp_name][class_profile][name].append(int(lev))
                        self.poly_len[exp_name][class_profile].append(int(end) - int(start) + 1)
                    if name == "cDNA":
                        self.cdna_len[exp_name][class_profile].append(int(end) - int(start) + 1)
                        
    def summarize(self):
    	s = ""
    	for experiment in self.experiments:
    		s += "-----------------------------------------------------------------------\n"
    		s += "Number of reads in experiment %s: %d\n" % (experiment, self.read_count[experiment])
    		for c in self.class_profile_count[experiment]:
    			s += "-----------------------------%s------------------------------------\n" % c
    			s += "Number of reads: %d\n" % (self.class_profile_count[experiment][c])
    			s += "****Read structure****\n"
    			for st in self.structure_profile_count[experiment][c]:
    				s += "%s: %d\n" % (st, self.structure_profile_count[experiment][c][st])
    			s += "****Read length****\n"
    			s += str(pd.Series(self.read_len[experiment][c]).describe()) + "\n" 
    			if c in self.cdna_len[experiment]:
    				s += "****cDNA length****\n"
    				s += str(pd.Series(self.cdna_len[experiment][c]).describe()) + "\n"   
    			if c in self.poly_len[experiment]:
    				s += "****polyA length****\n"
    				s += str(pd.Series(self.poly_len[experiment][c]).describe()) + "\n"   	
    			s += "****Levenshtein distance****\n"
    			for adapter in self.adapter_lev[experiment][c]:
    				s += "-- " + adapter + ":\n"
    				s += str(pd.Series(self.adapter_lev[experiment][c][adapter]).describe()) + "\n"   
    			for poly in self.poly_lev[experiment][c]:
    				s += "-- " + poly + ":\n"
    				s += str(pd.Series(self.poly_lev[experiment][c][poly]).describe()) + "\n" 
    	return s   
    				
    							
    
    #-------------------------------------------------------------------------------

	###################
	#  Plotting utils #
	###################
	#-------------------------------------------------------------------------------
    def get_theme(self):
        return theme(
            axis_line=element_line(size=2, color="black"),
            axis_ticks_length=5, axis_ticks_pad=20,
            text=element_text(size=8, color="black"),
            title=element_text(size=8, color="black"),
            panel_grid_major_x=element_blank(),
            panel_grid_minor_x=element_blank(),
            panel_grid_major_y=element_line(size=1, color="grey", linetype="dotted"),
            panel_grid_minor_y=element_line(size=1, color="grey", linetype="dotted"),
            legend_title=element_blank(),
            panel_background=element_rect(size=2, colour="black", fill="white"),
            strip_background_x=element_rect(size=2, colour="white", fill='lightsteelblue'),
            strip_background_y=element_rect(size=2, colour="white", fill='lightsteelblue'),
            axis_text_x=element_text(rotation=45, hjust=1),
            figure_size=self.default_fig_size,
        )
        
    def class_label_func(self, s):
    	return '\n'.join(s.split(","))   
    #-------------------------------------------------------------------------------

	#################
	#  Plot library #
	#################
	#-------------------------------------------------------------------------------
    def plot_class_counts(self, show_combined): 	
        counts_dict = self.class_profile_count if show_combined else self.class_count
        data = []
        for experiment in self.experiments:
            for profile in counts_dict[experiment]:
            	data.append([experiment, profile, self.class_label_func(profile), counts_dict[experiment][profile], 
            			     100*counts_dict[experiment][profile]/self.read_count[experiment]])
        df = pd.DataFrame(data, columns=["experiment", "class_profile", "class_profile_wrap", "count", "percent"])
        return [((ggplot(df)
                + facet_wrap(facets="~experiment", ncol=4, scales="free_y")
                + aes(x='class_profile', y='percent', fill='class_profile')
                + geom_bar(stat='identity')#, width=0.7)
                + scale_y_log10()
                + self.get_theme()
                + xlab("Artifacts") + ylab("log10(% reads)")
                + theme(legend_position="none")
                + ggtitle("Artifact profile (by experiment)" + ("" if show_combined else " -- counting reads assigned to each individual artifact class")))),
                ((ggplot(df)
                + facet_wrap(facets="~class_profile_wrap", ncol=4, scales="free_y")#, labeller=self.class_label_func)
                + aes(x='experiment', y='percent')
                + geom_bar(stat='identity', fill='darkseagreen')
                + self.get_theme()
                + xlab("Experiments") + ylab("% reads")
                + theme(legend_position="none")
                + ggtitle("Artifact profile (by artifact)" + ("" if show_combined else " -- counting reads assigned to each individual artifact class"))))]


    def plot_lens(self):
        data_rlen = []
        data_poly_len = []
        data_cdna_len = []
        for experiment in self.experiments:
            for c in self.class_profile_count[experiment]:
                for l in random.sample(self.read_len[experiment][c],
                                       min(self.plot_sample_size, len(self.read_len[experiment][c]))):
                    data_rlen.append([experiment, c, l])
                for l in random.sample(self.poly_len[experiment][c],
                                       min(self.plot_sample_size, len(self.poly_len[experiment][c]))):
                    data_poly_len.append([experiment, c, l])
                for l in random.sample(self.cdna_len[experiment][c],
                                       min(self.plot_sample_size, len(self.cdna_len[experiment][c]))):
                    data_cdna_len.append([experiment, c, l])
        df_rlen = pd.DataFrame(data_rlen, columns=["experiment", "class_profile",  "len"])
        df_poly_len = pd.DataFrame(data_poly_len, columns=["experiment", "class_profile",  "len"])
        df_cdna_len = pd.DataFrame(data_cdna_len, columns=["experiment", "class_profile",  "len"])
        return [(ggplot(df_rlen)
                + facet_wrap(facets="~experiment", ncol=4)
                + aes(x='class_profile', y='len', fill='class_profile')
                + geom_boxplot()
                + self.get_theme()
                + xlab("") + ylab("log10(length (bp))")
                + theme(legend_position="none")
                + scale_y_log10()
                + ggtitle("Read length profile (by artifact)")
                ),
                (ggplot(df_cdna_len)
                + facet_wrap(facets="~experiment", ncol=4)
                + aes(x='class_profile', y='len', fill='class_profile')
                + geom_boxplot()
                + self.get_theme()
                + xlab("") + ylab("log10(length (bp))")
                + theme(legend_position="none")
                + scale_y_log10()
                + ggtitle("cDNA length profile (by artifact)")
                ),
                (ggplot(df_poly_len)
                + facet_wrap(facets="~experiment", ncol=4)
                + aes(x='class_profile', y='len', fill='class_profile')
                + geom_boxplot()
                + self.get_theme()
                + xlab("") + ylab("log10(length (bp))")
                + theme(legend_position="none")
                + scale_y_log10()
                + ggtitle("polyA length profile (by artifact)")
                )]

    
    def plot_lev(self):
        data_adapters = []
        data_poly = []
        for experiment in self.experiments:
        	for c in self.class_profile_count[experiment]:
        		for adapter in self.adapter_lev[experiment][c]:
        			for lev in random.sample(self.adapter_lev[experiment][c][adapter], min(self.plot_sample_size, len(self.adapter_lev[experiment][c][adapter]))):
        				data_adapters.append([experiment, c, adapter, lev, lev])
        		for poly in self.poly_lev[experiment][c]:
        			for lev in random.sample(self.poly_lev[experiment][c][poly], min(self.plot_sample_size, len(self.poly_lev[experiment][c][poly]))):
        				data_poly.append([experiment, c, poly, lev, lev])
        df_adapter = pd.DataFrame(data_adapters, columns=["experiment", "class_profile", "adapter",  "lev", "order"])
        df_poly = pd.DataFrame(data_poly, columns=["experiment", "class_profile", "adapter",  "lev", "order"])
        plots = []
        for experiment in self.experiments:
        	plots.append((ggplot(df_adapter[df_adapter['experiment']==experiment])
                + facet_grid(rows="class_profile", cols="adapter", scales="free")
                + geom_bar(aes(x='lev', y='..prop..'), fill='darkseagreen') 
                + self.get_theme()
                + ggtitle("Adapter Levenshtein distance" + ": experiment " + experiment)
                + xlab("Levenshtein distance") + ylab("% reads")
                + theme(strip_text_y=element_text(rotation=360, hjust=1))
                + theme(text=element_text(size=6, color="black"))
                ))
        for experiment in self.experiments:
        	plots.append((ggplot(df_poly[df_poly['experiment']==experiment])
                + facet_grid(rows="class_profile", cols="adapter", scales="free")
                + geom_bar(aes(x='lev', y='..prop..'), fill='darkseagreen') 
                + self.get_theme()
                + ggtitle("polyA/polyT Levenshtein distance" + ": experiment " + experiment)
                + xlab("Levenshtein distance") + ylab("% reads")
                + theme(strip_text_y=element_text(rotation=360, hjust=1))
                + theme(text=element_text(size=6, color="black"))
                ))
        return plots
                
                
    def plot_offset(self):
        data_start = []
        data_end = []
        for experiment in self.experiments:
        	for c in self.class_profile_count[experiment]:
        		for adapter in self.start_adapter_offset[experiment][c]:
        			for offset in random.sample(self.start_adapter_offset[experiment][c][adapter], min(self.plot_sample_size, len(self.start_adapter_offset[experiment][c][adapter]))):
        				data_start.append([experiment, c, adapter, offset])
        		for adapter in self.end_adapter_offset[experiment][c]:
        			for offset in random.sample(self.end_adapter_offset[experiment][c][adapter], min(self.plot_sample_size, len(self.end_adapter_offset[experiment][c][adapter]))):
        				data_end.append([experiment, c, adapter, offset])
        df_start = pd.DataFrame(data_start, columns=["experiment", "class_profile", "adapter",  "offset"])
        df_end = pd.DataFrame(data_end, columns=["experiment", "class_profile", "adapter",  "offset"])
        plots = []
        for experiment in self.experiments:
        	plots.append((ggplot(df_start[df_start['experiment']==experiment])
                + facet_grid(rows="class_profile", cols="adapter", scales="free")
                + geom_bar(aes(x='offset', y='..count..'), fill='cadetblue') 
                + self.get_theme()
                + scale_y_log10()
                + ggtitle("Left adapter offset" + ": experiment " + experiment)
                + xlab("Offset from terminal") + ylab("log10(Number of reads)")
                + theme(strip_text_y=element_text(rotation=360, hjust=1))
                + theme(text=element_text(size=6, color="black"))
                ))
        for experiment in self.experiments:
        	plots.append((ggplot(df_end[df_end['experiment']==experiment])
                + facet_grid(rows="class_profile", cols="adapter", scales="free")
                + geom_bar(aes(x='offset', y='..count..'), fill='cadetblue') 
                + self.get_theme()
                + scale_y_log10()
                + ggtitle("Right adapter offset" + ": experiment " + experiment)
                + xlab("Offset from terminal") + ylab("log10(Number of reads)")
                + theme(strip_text_y=element_text(rotation=360, hjust=1))
                + theme(text=element_text(size=6, color="black"))
                ))
        return plots

    #-------------------------------------------------------------------------------

def main():
    print("*********************************")
    print("*  marti: classification stats  *")
    print("*********************************")

    # ------ CLI ------
    parser = argparse.ArgumentParser(description='Plot marti classification results')
    parser.add_argument('--bams', required=True, nargs='+', help='List of BAM files classified by marti')
    parser.add_argument('--experiments', required=True, nargs='+',
                        help='List of unique experiment names corresponding to each BAM')
    parser.add_argument('--output_dir', required=True, help='Output directory for the report')
    parser.add_argument('--max_reads_per_experiment', type=int, default=None, help='Maximum reads to load from each BAM')
    parser.add_argument('--plot_sample_size', default=100000, help='Maximum number of data points to include in distribution plots')
    args = parser.parse_args()
    # -----------------
    print("Processing BAM file(s)...")
    marti_stats = MartiStats(args.bams, args.experiments, args.max_reads_per_experiment, args.plot_sample_size)
    print("Generating plots...")
    plots = []
    plots.extend(marti_stats.plot_class_counts(True))
    plots.extend(marti_stats.plot_class_counts(False))
    plots.extend(marti_stats.plot_lens())
    plots.extend(marti_stats.plot_lev())
    plots.extend(marti_stats.plot_offset())
    with PdfPages(args.output_dir + '/marti_stats.pdf') as report_pdf:
        print("Saving plots...")
        for plot in plots:
            report_pdf.savefig(plot.draw())
        print("Saving basic stats...")
        with open(args.output_dir + '/marti_stats.txt', "w") as report_txt:
        	report_txt.write(marti_stats.summarize())
    print("Marti PDF report has been written to: ", args.output_dir + "/marti_stats.pdf")
    print("Marti additional stats have been written to: ", args.output_dir + "/marti_stats.txt")
    
if __name__ == "__main__":
    main()
