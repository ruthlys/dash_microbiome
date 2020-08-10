# load necessary packages
library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(plotly)
library(dplyr)
library(base64)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(tidyverse)
library(picante)
library(data.table)
library(RColorBrewer)

#################################### PERFORM ANALYSIS AND CREATE DATA FRAMES #############################

# # load data and convert to matrix
read_csv <- function(x){
  df <- read.csv(x, sep="\t", header=T, row.names=1)
  as.matrix(df)
}

otubac <- read_csv("~/Dropbox/Work/INRS/Data/NI_experiment/Files/16S_ITS/16s_otu.txt")
taxbac <- read_csv("~/Dropbox/Work/INRS/Data/NI_experiment/Files/16S_ITS/16s_tax.txt")

# load phylogenetic information
read_tree <- function(x){
  read.tree(x)
}

treebac = read_tree("~/Dropbox/Work/INRS/Data/NI_experiment/Files/16S_ITS/16S.nwk")

# create sample data (and glimpse on structure (will have to be uploaded)
create_samdata <- function(x){
  sample_data(data.frame(
    SoilType = rep(c("Bulk soil","Rhizosphere"), each = 2),
    Date = rep(c("20190723", "20190905")),
    Treatment = rep(c("Control","Nitrapyrin"), each = 24),
    row.names=sample_names(x),
    stringsAsFactors=FALSE
  ))
}

# convert OTU, TAX, sample data and tree file into one global phyloseq objects
to_phyloseq <- function(otu,tax,read_tree){
  OTU <- otu_table(otu, taxa_are_rows = TRUE)
  TAX <- tax_table(tax)
  physeq <- phyloseq(OTU,TAX)
  sam_data <- create_samdata(physeq)
  phyloseq(OTU, TAX, sam_data, read_tree)
}

psbac <- to_phyloseq(otubac, taxbac, treebac)

# prepare data
psbac1 = prune_taxa(taxa_sums(psbac) > 0, psbac)
# subtract the number of OTUs in original (physeq) with number of OTUs in new phyloseq object
psbac1a = prune_taxa(taxa_sums(psbac1) > 1, psbac1)
psbac2 = subset_taxa(psbac1a, Phylum !="Unclassified") %>% subset_taxa(Genus !="Unclassified")
# transform abundance counts to fractional abundance
psbac2.rel = transform_sample_counts(psbac2, function(x) x /sum(x)) 
# filter mean abundance above 1% mean rel abundance (loss of low abundance AOB!)
psbac2.rel = filter_taxa(psbac2.rel, function(x) mean(x) > 0.001, TRUE)
# calculate percentage 
psbac2.rel = transform_sample_counts(psbac2.rel, function(x) 100 * x/sum(x))

# # calculate alpha diversity
# alphabac = estimate_richness(psbac1, measures=c("Shannon", "InvSimpson"))
# # get the metadata out as separate object
# alphabac.meta = meta(psbac1)
# # add the rownames as a new colum for integration later
# alphabac.meta$sam_name = rownames(alphabac.meta)
# # add the rownames to diversity table
# alphabac$sam_name = rownames(alphabac)
# # merge these two data frames into one
# alphabac.df = merge(alphabac.meta, alphabac, by = "sam_name")
# # calculate diversity (Faith's PD)
# psbac1.div = as.data.frame(psbac1@otu_table)
# psbac1.tree = psbac1@phy_tree
# # calculate PD and SR using picante
# psbac1.div = pd(t(psbac1.div), psbac1.tree, include.root=F) # t(otu_table) transposes the table for use in picante
# psbac1.div = setDT(psbac1.div, keep.rownames = "sam_name")
# # add results of PD to div.df file
# alphabac.df$Phylogenetic_Diversity = psbac1.div$PD
# colnames(alphabac.df) = c("sam_name", "SoilType", "Date", "Treatment", "Shannon", "InvSimpson", "Phylogenetic Diversity")
# alphabac.df = reshape2::melt(alphabac.df)
# 
# # continue processing data and bar charts of top 10 phyla and genera
# # remove singletons
# psbac1a = prune_taxa(taxa_sums(psbac1) > 1, psbac1)
# psbac2 = subset_taxa(psbac1a, Phylum !="Unclassified") %>% subset_taxa(Genus !="Unclassified")
# # transform abundance counts to fractional abundance
# psbac2.rel = transform_sample_counts(psbac2, function(x) x /sum(x))
# # filter mean abundance above 1% mean rel abundance (loss of low abundance AOB!)
# psbac2.rel = filter_taxa(psbac2.rel, function(x) mean(x) > 0.001, TRUE)
# # calculate percentage
# psbac2.rel = transform_sample_counts(psbac2.rel, function(x) 100 * x/sum(x))
# # convert to df
# psbac.df = psmelt(psbac2.rel)
# 
# # agglomerate taxa at phylum level and extract top 10 phyla for 16S
# bacglomp = tax_glom(psbac2.rel, taxrank = "Phylum")
# bactop10p = sort(tapply(taxa_sums(bacglomp), tax_table(bacglomp)[, "Phylum"], sum), decreasing = TRUE)[1:10]
# bactop10p = subset_taxa(bacglomp, Phylum %in% names(bactop10p))
# # agglomerate taxa at genus level and extract top 10 genera for 16S
# bacglomg = tax_glom(psbac2.rel, taxrank = "Genus")
# bactop10g = sort(tapply(taxa_sums(bacglomg), tax_table(bacglomg)[, "Genus"], sum), decreasing = TRUE)[1:10]
# bactop10g = subset_taxa(bacglomg, Genus %in% names(bactop10g))
# 
# # create dataframe from phyloseq object
# create_dataframe <- function(dataframe){
#   data.table(psmelt(dataframe))
# }
# 
# top10p.df <- create_dataframe(bactop10p)
# top10g.df <- create_dataframe(bactop10g)
# 
# # alternatively save and read in files
# fwrite(alphabac.df, file="/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/alphabac.csv", row.names = FALSE)
# fwrite(psbac.df, file="/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/psbac_melt.csv", row.names = FALSE)
# fwrite(top10p.df, file="/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/bactop10p.csv", row.names = FALSE)
# fwrite(top10g.df, file="/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/bactop10g.csv", row.names = FALSE)

####################################################################################################

#################################### APP START $ LOAD FILES ########################################

app <- Dash$new()
# Initiate application

# load files
read_csv <- function(x){
  df <- data.frame(fread(x, header=T))
}

# read in files that were previously analyzed
alphabac.df<- read_csv("/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/alphabac.csv")
psbac.df <- read_csv("/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/psbac_melt.csv")
top10p.df <- read_csv("/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/bactop10p.csv")
pca <- read.csv("/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/pca_16S.csv")
pca$Date <- as.factor(pca$Date)

# extract unique elements
soil <- unique(alphabac.df$SoilType)
date <- unique(alphabac.df$Date)
alphadiv_var <- unique(alphabac.df$variable)
yaxis_type <- c("log", "linear")

# Define the number of colors 
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

createBox <- function(alphadiv_file, yaxis_type){
  # create graph parameters 
  plot_ly(
    data = alphadiv_file, 
    x = ~variable, y = ~value,
    color = ~Treatment, colors = "Dark2",
    type = "box", alpha = 0.7) %>% 
    layout(
      boxmode = "group", 
      title = "Alpha diversity across treatments", 
      xaxis = list(title = 'Alpha Diversity Measure'),
      yaxis = list(title = 'Value', type = yaxis_type))
}

# create 16S stacked bar chart
createStackedbar <- function(stackedbar_file){
  stackedbar_file %>% mutate(Phylum=fct_reorder(Phylum, Abundance, .fun='mean', .desc=TRUE)) %>%
  plot_ly(
    x = ~Treatment, y = ~Abundance,
    color = ~Phylum, colors = "Dark2", alpha = 1,
    name = stackedbar_file$Phylum) %>% 
    layout(
    title = "Relative abundance of Phyla", 
    yaxis = list(title = "Relative abundance"), # ticksuffix = "%"
    barmode = "stack")
}

# # create plot of top10 phyla and genera 
# createBoxtop10p <- function(top10p_file, yaxis_type){
#   # sort in desc order of mean abundance 
#   top10p_file %>% mutate(Phylum=fct_reorder(Phylum, Abundance, .fun='mean', .desc=TRUE)) %>% 
#   # create graph parameters 
#   plot_ly(
#     x = ~Phylum, y = ~Abundance,
#     color = ~Treatment, colors = "Dark2",
#     type = "box", alpha = 0.7) %>% 
#     layout(
#       boxmode = "group", 
#       title = "Top 10 Phyla in decreasing order of mean abundance", 
#       xaxis = list(title = 'Phylum'),
#       yaxis = list(title = 'Relative abundance', type = yaxis_type))
# }

# # create 3d pca plot 
# create3dScatter <- function(pca_file){
#   plot_ly(
#     pca_file, 
#     x = ~x, 
#     y = ~y, 
#     z = ~z, 
#     color = ~Treatment, 
#     colors = "Dark2", 
#     mode = "markers", 
#     type = "scatter3d") %>% 
#     layout(
#       title = "Principal Coordinates Analysis (PCoA) of Volatile Organic Compounds (VOCs)",
#       scene = list(xaxis = list(title = 'PC1 (37.1%)'),
#                    yaxis = list(title = 'PC2 (13.3%)'),
#                    zaxis = list(title = 'PC3 (5.7%)')))
# }

# create PCoA plot 
m2 <- list(
  l = 0,
  r = 0,
  b = 20,
  t = 30
  # pad = 4
)

t <- list(
  family = "Open Sans",
  size = 12,
  color = "#444444",
  m = list(
    t = "25px"))


createNet <- function(phyloseq_obj) {
  ig <- make_network(phyloseq_obj, max.dist=0.3)
  p <- plot_network(ig, phyloseq_obj, 
                    color= "Date", 
                    shape = "Treatment",
                    alpha = 0.8,
                    line_weight=0.4, label = NULL) + scale_colour_brewer(palette="Dark2")
  p[["labels"]][["colour"]] <- NULL
  p[["labels"]][["shape"]] <- NULL
  p1 <- ggplotly(p) %>% style(showlegend = FALSE) %>% layout(title = "Microbiome network", margin = m2, font = t)
  p1
}

x <- list(
  title = "Axis 1"
)
y <- list(
  title = "Axis 2"
)

createPca <- function(pca_file) {
  plot_ly(pca_file, 
                 x= ~Axis.1, 
                 y= ~Axis.2, 
                 color = ~Treatment, 
                 colors = "Dark2", 
                 mode = "markers", 
                 type = "scatter",
                 alpha = 1,
                 marker = list(size = 13)) %>%
    layout(title = 'PCoA plot', xaxis = x, yaxis = y)
}

####################################################################################################

#################################### CREATE LAYOUT #################################################

app$layout(
  htmlDiv(
    list(
      # create and style title
      htmlDiv(
        list(
          htmlH2(
            id = "title",
                 'Interactive Microbiome Viewer', 
                 style = list(
                   textAlign = 'left'
                 )
               ),
          htmlH5(
            'Statistical analysis and visualization of microbiome data', 
            style = list(
              textAlign = 'left'
            )
          )
        ), 
      className = "banner"
      ),

      # create top container and graph child, set style 
      htmlDiv(
        list(
          htmlDiv(
            id = "top-graphs",
            style = list(
              display = "flex", 
              borderRadius = '5px',
              justifyContent = "space-evenly", 
              width = "100%",
              # height = "50px",
              margin = "0 auto"
            ),
            children = list(
              # create tabs
              htmlDiv(
                id = "control-tabs",
                className = "control-tabs",
                children = list(
                  dccTabs(
                    id="tabs", 
                    # style = list(colors = list(primary = "white")),
                    value = "what-is",
                    children = list(
                      dccTab(
                        label = "About",
                        value = "what-is",
                        children = htmlDiv(
                          className = "control-tab",
                          children = list(
                            htmlH4(
                              className = "what-is", 
                              children = "What is Alpha Diversity?"
                            ),
                            dccMarkdown(
                              'Alpha diversity describes the diversity within a particular sample or environment 
                              (for more info have a look at [Amy D. Willis paper](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)).
                              Different alpha diversity metrices can be selected to explore alpha diversity.'
                            ),
                            htmlH4(
                              className = "what-is", 
                              children = "What is Beta Diversity?"
                            ),
                            dccMarkdown(
                              'Beta diversity measures changes in diversity of species from one environment or sample to another. 
               This app performs principal coordinate analysis (PCoA) of a distance matrix, including two correction methods for negative eigenvalues.'
                            )
                          )
                        )
                      ),
                      
                      dccTab(
                        label = "Parameters",
                        value = "parameters",
                        children = htmlDiv(
                          className = "control-tab",
                          children = list(
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select alpha diversity"
                                ),
                                dccDropdown(
                                  id ="d_alphadiv",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px"),
                                  options = lapply(alphadiv_var, # create list of unique elements (alphadiv_var) and create labels using lapply 
                                                   function(alphadiv_var) {
                                                     list(label = alphadiv_var,
                                                          value = alphadiv_var)
                                                   }
                                  ),
                                  value = c('Shannon', 'InvSimpson'),
                                  multi = TRUE
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select Date"
                                ),
                                dccDropdown(
                                  id = "d_date-box",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px"),
                                  # style = list(width = "45%"),
                                  options = lapply(list("20190723", "20190905"),
                                                   function(x){
                                                     list(label = x, value = x)
                                                   }
                                  ),
                                  value = '20190723'
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select Soil compartment"
                                ),
                                dccDropdown(
                                  id = "d_soiltype",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px", width = "100%"),
                                  options = lapply(list("Rhizosphere", "Bulk soil"),
                                                   function(x){
                                                     list(label = x, value = x)
                                                   }
                                  ),
                                  value = 'Rhizosphere'
                                )
                              )
                            ), 
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Choose scale"
                                ),
                                dccRadioItems(
                                  id ="r_datatype",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px", width = "50%"),
                                  options=lapply(list("log", "linear"),
                                                 function(x){
                                                   list(label = x, value = x)
                                                 }
                                  ),
                                  value = "log",
                                  labelStyle = list(display = 'inline-block')
                                )
                              )
                            )
                          )
                        )
                      )
                    )))),
              
              htmlDiv(
                id = "left-top-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "alphadiv-16S-graph",
                        figure = createBox(alphabac.df, yaxis_type),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container top left
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "2rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "38%", 
                  float = "none", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              ),
              
              htmlDiv(
                id = "right-top-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "stacked-bar-16S",
                        figure = createStackedbar(psbac.df),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container top right
                style = list(
                  marginTop = "10px",
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "2rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "35%",
                  float = "none",
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              )
            )
          )
        )
      ),
      
      # create top container and graph child, set style 
      htmlDiv(
        list(
          htmlDiv(
            id = "bottom-graphs",
            style = list(
              display = "flex",
              borderRadius = '5px',
              justifyContent = "space-evenly", 
              width = "100%",
              margin = "0 auto"
            ),
            children = list(
              htmlDiv(
                id = "left-bottom-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "network",
                        figure = createNet(psbac2.rel),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container bottom left
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "3rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "48%", 
                  float = "none", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              ),

              htmlDiv(
                id = "right-bottom-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "3d-graph",
                        figure = createPca(pca)
                      )
                    )
                  )
                ),
                # style container bottom right
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "3rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "48%", 
                  float = "none", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              )
            )
          )
        )
      )
    )
  )
)

####################################################################################################

#################################### CALLBACKS #####################################################

app$callback(
  output=list(id="alphadiv-16S-graph", property="figure"),
  list(input(id='d_alphadiv', property='value'),
    input(id='d_date-box', property='value'),
    input(id='d_soiltype', property='value'), 
    input(id='r_datatype', property='value')),
  function(selected_variable, selected_date, selected_soiltype, yaxis_datatype){
    filtered_variable <- alphabac.df %>% dplyr::filter(variable %in% c(selected_variable), Date %in% c(selected_date), 
                                                       SoilType %in% c(selected_soiltype))
    createBox(filtered_variable, yaxis_datatype)
  }
)

app$callback(
  output(id='stacked-bar-16S', property='figure'),
  list(input(id='d_date-box', property='value'),
       input(id='d_soiltype', property='value')),
  function(selected_date, selected_soiltype){
    filtered_date <- psbac.df %>% dplyr::filter(Date %in% c(selected_date), SoilType %in% c(selected_soiltype)) 
    createStackedbar(filtered_date)
  }
)

app$callback(
  output=list(id="3d-graph", property="figure"),
  list(input(id='d_date-box', property='value')),
  function(selected_date){
    filtered_variable <- pca %>% dplyr::filter(Date %in% c(selected_date))
    createPca(filtered_variable)
  }
)

# run app
app$run_server(debug=F, threaded=TRUE, showcase = T)

