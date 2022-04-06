# function to plot GO term enrichment

go_plot <- function(genes, background) {
  require (gprofiler2)
  
  # use gprofiler2 function to calculate enrichment
  go <- gost(query = genes, organism = "athaliana", custom_bg = background, user_threshold = 0.05, 
             sources = "GO")$result
  
  go$log10_pvalue <- -log10(go$p_value)
  
  # set nice order
  go$term_name <- factor(go$term_name, 
                         levels = as.factor(go$term_name [order(go$source, go$log10_pvalue, decreasing = F)]))
  
  g <- ggplot(go, aes(log10_pvalue, term_name, fill = source)) +
          geom_bar(stat="identity", position=position_dodge()) +
          theme_classic() + 
          scale_fill_tableau()
  
  return (g)
}
