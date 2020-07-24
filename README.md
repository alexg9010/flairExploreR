# flairExplorer

flairExplorer is a package for the visualization and interactive exploration of 
results generated from the [flair pipeline](https://github.com/BrooksLabUCSC/flair).

It provides a bash script to run a series of commands from the flair pipeline in a batch mode
and a shiny based dashboard to visualize results in a aggregated manner. 

While the flair pipeline itself includes some plotting scripts to visualize isoform usage, 
other data generated in the process can be best inspected in IGV, for example the productivity bed file. 
And there was no script to visualize the differential analysis results in a comprehensive manner, 
such as a volcano plot. This is why I came up with the interactve dashboard,
to visualize results in a gene centric manner, allowing to directly connect the outcome of all three
analysis, differential gene and transcript usage and differential splicing.

![Gene Centric Summary plot](misc/figure_panel.png)

In addition all the resulsts of the differential analysis are provided as interactive tables 
with an additional overall gene filter at shorthand. 

![Interactive Result tables](misc/table_panel.png)



## Installation


