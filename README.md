# shinydotplot
Simple code to generate an interactive dotplot to explore cell type-specific expression profiles

Instructions for testing R shiny interactive dotplot

1) Download and install required R packages: shiny, Matrix

2) Copy gene expression matrix and cluster identity vector to "starting_data" folder as rda objects. Currently, there are two dummy data files in there, with randomly generated gene expression values and cluster identities. The expression matrix should have genes in rows and cells in columns, and the cluster identity vector should have names that match the cell names in the expression matrix.

3) Run the generate_files.r script. This script will create an annotation matrix and a matrix with fractions of cells expressing each gene in each cluster. The order in which clusters should be displayed in the widget needs to be specified in this file, on line 10 (variable "cluster_order"). Make sure to create a directory called "generated_files" before running this script, or else it will not be able to save the output matrices correctly.

4) Once the generate_files.r script has been run, open an instance of Rstudio and set the directory to the shiny_dotplot_template folder. Include the shiny library (using the require(shiny)) command, and run the widget by typing runApp() in the console. This should pop up the widget and allow for plotting/DE gene identification/etc.

5) Extra modifications: built into the widget is an option to show "all" clusters, "neuron" clusters, or "glia" clusters. "All" will select all clusters, but the "neuron" and "glia" functions are hardwired to search for "GABA", "Glut", "Ex", or "Inh" in the cluster names in order to identify neuronal and non-neuronal clusters. To change these categories (or introduce new ones), change lines 37-41 & lines 162-181.
