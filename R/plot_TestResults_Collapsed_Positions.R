#' Plot collapsed tests results for positions
#'
#' Plot collapsed tests results for Gene position and Relative to Island position
#'
#' @param x list with datataframes to plot
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_TestResults_Collapsed_Positions <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   x.df <- bind_rows(x, .id = "test")
   #..# x.df[1:3] <- lapply(x.df[1:3], function(xf) as.numeric(levels(xf))[xf])

   p <- ggplot(x.df, aes(x = colnames(x.df[2]), y = x.df[3], fill = test)) +
      geom_bar(stat="identity",   position=position_dodge(), width = 0.5) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      xlab("") +
      ylab("")


   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      ggplot2::ggsave(paste0(file.path( outputdir),"/",
                             unique(x.df$test)[which(length(unique(x.df$test))==min(length(unique(x.df$test))))],
                             "_",
                             filename,
                             ".pdf"),
                      p)
   }

   return(p)
}
