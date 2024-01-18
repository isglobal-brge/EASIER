#' Plot collapsed tests results
#'
#' Plot collapsed tests results
#'
#' @param x list with datataframes to plot
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_TestResults_Collapsed <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   x.df <- bind_rows(x, .id = "test")
   #.. 04/01/2022 ..# x.df[3:6] <- lapply(x.df[3:6], function(xf) as.numeric(levels(xf))[xf])
   x.df[,3:6] <- lapply(x.df[,3:6], function(xf) as.numeric(xf))
   colnames(x.df)[2] <- "Position"

   p <- ggplot(x.df, aes(x = Position, y = OR, fill = test)) +
      geom_bar(stat="identity",   position=position_dodge(), width = 0.8) +
      geom_errorbar(aes(ymin=OR.inf, ymax=OR.sup), width=0.2, position=position_dodge(.5)) +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed") # +
      # xlab("") +
      # ylab("")


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
