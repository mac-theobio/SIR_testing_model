library(McMasterPandemic)
library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)
library(latex2exp)
library(rootSolve)
library(cowplot) ## for ggdraw()/rectangle-hacking

library(shellpipes)

texfile <- "modeldefs.tex"
unlink(texfile)

rpcall("sir_plot.Rout sir_plot.R params.rda SIRfunctions.rda")

loadEnvironments()

## setwd("/home/ag/projects/SIRmodel/codes/")

# #################################
# Making Dataframe part
# #################################
# R0 contour plots
n_out <- 5 ## facets (rows/cols)
n_in <- 41 ## grid N within facets

#choose weight W_S 0.3 for TTI testing plan and 1 for the random testing.
W_S_random <- 1
W_S_targeted <- 0.3
## tweak to get non-overlapping axis labels:
##  adjust grid to be slightly bigger than ideal tick marks
inv_omega_min <- 0.5 ##inverse of omega in days
inv_omega_max <- 12
inv_omega_brks <- c(1,5,10)
rho_min <- 0
rho_max <- 0.013
rho_brks <- c(0,0.005,0.01)

inv_omega_fast_min <- 0.5
inv_omega_fast_max <- 5
inv_omega_fast_brks <- c(1,2,4)
rho_fast_min <- 0
rho_fast_max <- (1/inv_omega_fast_max)-0.000001
rho_fast_brks <- c(0,0.1,0.17)

catt("% parameter values\n")
latexout(params[["beta"]],"betaparam", file="modeldefs.tex")
latexout(params[["gamma"]],"gammaparam", file="modeldefs.tex")
latexout(1/params[["gamma"]],"invgammaparam", file="modeldefs.tex")
latexout(params[["beta"]]/params[["gamma"]],"Rnumparam", file="modeldefs.tex")
latexout(params[["beta"]]-params[["gamma"]],"rparam", file="modeldefs.tex")


background <- "grey90"
# make dataframe
df_random <- make_params_dat(params = params,
                eta_ws=0,eta_we=1, ## so theta_w
                eta_cs=0,eta_ce=1, ## so theta_c
                omega_s=1/inv_omega_max,omega_e=1/inv_omega_min,
                rho_s=rho_min,rho_e=rho_max)

df_targeted <- make_params_dat(params = update(params,W_S=W_S_targeted),
                eta_ws=0,eta_we=1, ## so theta_w
                eta_cs=0,eta_ce=1, ## so theta_c
                omega_s=1/inv_omega_max,omega_e=1/inv_omega_min,
                rho_s=rho_min,rho_e=rho_max)

## high rho and omega
df_random_h <- make_params_dat(params = params,
                             eta_ws=0,eta_we=1, ## so theta_w
                             eta_cs=0,eta_ce=1, ## so theta_c
                             omega_s=1/inv_omega_fast_max,omega_e=1/inv_omega_fast_min,
                             rho_s=rho_fast_min,rho_e=rho_fast_max)


df_targeted_h <- make_params_dat(params = update(params,W_S=W_S_targeted),
                               eta_ws=0,eta_we=1, ## so theta_w
                               eta_cs=0,eta_ce=1, ## so theta_c
                               omega_s=1/inv_omega_fast_max,omega_e=1/inv_omega_fast_min,
                               rho_s=rho_fast_min,rho_e=rho_fast_max)


##important contour, ie R0=1 thus threshold=1 when plotting R0 contours, or Delta(R0=1)
threshold <- 1-(params[["gamma"]]/params[["beta"]]) ## corresponding to R0=1
# #################################
# Plotting Part
# #################################
# Setting the plots
## hacked (not very robustly) to chain label_both() and label_parsed() ...
label_special <- function (labels, multi_line=FALSE, sep = "== ") {
  value <- ggplot2:::label_value(labels, multi_line = multi_line)
  variable <- ggplot2:::label_variable(labels, multi_line = multi_line)
  variable[[1]] <- gsub("_(.)$","[\\1]",variable[[1]])
  ## not using multiple faceting variables on each margin, so forget paste
  out <- Map(paste, variable, value, sep = sep)
  out[[1]] <- lapply(out[[1]], function(expr) c(parse(text=expr)))
  out
}

## truncate numeric values to (default) 2 digits
hack_breaks <- function(x,digits=2) {
  regex <- sprintf("0.([0-9]{%d})[0-9]+",digits)
  gsub(regex,"0.\\1",x)
}

n_contour <- 8 ## number of desired contours or bins
bins <- cut(sort(unique(c(df_random$Delta,df_targeted$Delta))), n_contour)
brks <- levels(bins)
brks_vec <- unique(as.numeric(unlist(lapply(strsplit(brks, ","), function(x) gsub("\\(|]", "", x)))))
brks_vec <- sort(ifelse(brks_vec[]< 0,0,brks_vec[]), decreasing = TRUE) ## replace neg with 0 and reorder

df_temp <- df_random
# df_temp <- df_targeted
p1 <- (ggplot(df_temp,aes(x=1/omega,y=rho,z=Delta))
    + theme_bw()
    + xlab(TeX('$\\1/omega$, mean test return time (day)'))
    + ylab(TeX('$\\rho$, testing intensity (1/day per 1000)'))
)
## There might be a solution for the tick label collision with expand_limits(), but that's hard to do across facets
## see https://stackoverflow.com/questions/41575045/avoiding-axis-tick-label-collision-in-faceted-ggplots
p1_temp <- (p1
  # + geom_contour_filled(breaks=brks)
  + geom_contour_filled(breaks=brks_vec)
  ##+ geom_contour(breaks=threshold,alpha=0.5,colour="black",lty="11",lwd=1)
  + facet_grid(theta_w~theta_c, labeller=label_special)
  + scale_x_continuous(expand=expansion(c(0,0)), breaks=inv_omega_brks,limits=c(inv_omega_min,inv_omega_max))
  ## scale by 1000, no digits after decimal
  + scale_y_continuous(expand=expansion(c(0,0)), breaks=rho_brks,
                       labels=scales::number_format(scale=1000,accuracy=1))
  + scale_fill_viridis_d(name=parse(text="Delta"),drop=FALSE,
                         labels=hack_breaks)
  + geom_rect(data=df_temp, fill=ifelse(df_temp$theta_c < df_temp$theta_w,background,"NA" ),
              color= NA,
              ymin=-1,
              ymax=10,
              xmin=-Inf,
              xmax=Inf)
  + theme(panel.spacing=grid::unit(0,"lines"))
)
print(p1_temp)

# #################################
# 1. Plot the Random Testing Scenario:
# #################################

ggsave(p1_temp + ggtitle(TeX('Random testing, $w_S=w_I=1$')) +
       theme(legend.position = "none"),
       filename = "R0contour_random.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 2. Plot the Targeted Testing Scenario:
# #################################

## mildly ugly hack to squeeze in the legend:
##  * remove title when plotting legend
##  * adjust x, y, text size so that it _just_ fits in the corner_
##  * use cowplot::ggdraw() to draw the title of the legend separately
## Otherwise, it's almost impossible to make the breaks large enough to see,
## without having the corner of the legend box obscure part of one of the tiles
## * Changing the legend background to transparent is ugly
## * Changing the drawing order of components so that the
##   legend box gets drawn before the tile (overlaid by it) but
##   after grid lines etc. (overlays them) seems nearly impossible
p1_targeted <- ((p1_temp %+% df_targeted)
  + ggtitle(TeX(sprintf("Targeted testing, w_S=%.1f, w_I=1",W_S_targeted)))
  ## REMOVE legend title here ... add it back in with ggdraw
  + theme(
        ##axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.14, 0.30),
        legend.text = element_text(size = 7),
        legend.background = element_rect(
            fill=background,
            ## adjustcolor("gray",alpha.f=0.1),
            ## fill=NA,
            ## skinny edge
            colour = "black",
            size=0.2))
)

p1_targeted_hacked <- (ggdraw(p1_targeted)
  + draw_label(TeX("$\\Delta$"), ## \u0394",  ## Unicode Delta?
               x=0.15,
               y=0.6,
               size=20)
)

ggsave(p1_targeted_hacked,
       filename = "R0contour_TTI.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 3. Plots for high rho Scenario:
# #################################
bins <- cut(sort(unique(c(df_random_h$Delta,df_targeted_h$Delta))),
            n_contour)
brks <- levels(bins)
brks_vec <- unique(as.numeric(unlist(lapply(strsplit(brks, ","), function(x) gsub("\\(|]", "", x)))))
brks_vec <- sort(ifelse(brks_vec[]< 0,0,brks_vec[]), decreasing = TRUE) ## replace neg with 0 and reorder

df_temp2 <- df_random_h

p12 <- (ggplot(df_temp2,aes(x=1/omega,y=rho,z=Delta))
       + theme_bw()
       + xlab(TeX('$\\1/omega$, mean test return time (day)'))
       + ylab(TeX('$\\rho$, testing intensity (1/day per 1000)'))
)

p1_temp2 <- (p12
            # + geom_contour_filled(breaks=brks)
            + geom_contour_filled(breaks=brks_vec)
            ##+ geom_contour(breaks=threshold,alpha=0.5,colour="black",lty="11",lwd=1)
            + facet_grid(theta_w~theta_c, labeller=label_special)
            + scale_x_continuous(expand=expansion(c(0,0)), breaks=inv_omega_fast_brks,limits=c(inv_omega_fast_min,inv_omega_fast_max))
            ## scale by 1000, no digits after decimal
            # + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
            + scale_y_continuous(expand=expansion(c(0,0)), breaks=rho_fast_brks,
                                 labels=scales::number_format(scale=1000,accuracy=1))
            + scale_fill_viridis_d(name=parse(text="Delta"),drop=FALSE,
                                   labels=hack_breaks)
            + geom_rect(data=df_temp, fill=ifelse(df_temp$theta_c < df_temp$theta_w,background,"NA" ),
                        color= NA,
                        ymin=-1,
                        ymax=10,
                        xmin=-Inf,
                        xmax=Inf)
            + theme(panel.spacing=grid::unit(0,"lines"))
)
print(p1_temp2)
# #################################
# 4. Plot the Random Testing Scenario:
# #################################

ggsave(p1_temp2
       + ggtitle(TeX('Random testing, $w_S=w_I=1$'))
       + theme(legend.position = "none"),
       filename = "R0contour_random2.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 5. Plot the Targeted Testing Scenario:
# #################################

p1_targeted_h <- ((p1_temp2 %+% df_targeted_h)
                + ggtitle(TeX(sprintf("Targeted testing, w_S=%.1f, w_I=1",W_S_targeted)))
                ## REMOVE legend title here ... add it back in with ggdraw
                + theme(
                  ##axis.title.y=element_blank(),
                  legend.title = element_blank(),
                  legend.position = c(0.14, 0.30),
                  legend.text = element_text(size = 7),
                  legend.background = element_rect(
                    fill=background,
                    ## adjustcolor("gray",alpha.f=0.1),
                    ## fill=NA,
                    ## skinny edge
                    colour = "black",
                    size=0.2))
)

p1_targeted_hacked <- (ggdraw(p1_targeted_h)
                       + draw_label(TeX("$\\Delta$"), ## \u0394",  ## Unicode Delta?
                                    x=0.15,
                                    y=0.6,
                                    size=20)
)

ggsave(p1_targeted_hacked,
       filename = "R0contour_TTI2.pdf" ,
       width = 12, height = 12, units = "cm")

# #################################
# 6. Check: in random testing when rho is high, $\theta_w=0$/$\theta_c=0.75$ panel
# #################################

df_random_h_check <- df_random_h %>%
  filter(theta_w==0 & theta_c==.75)
p_temp <- (ggplot(df_random_h_check,aes(x=1/omega,y=rho,z=Delta))
        + theme_bw()
        + xlab(TeX('$\\1/omega$, mean test return time (day)'))
        + ylab(TeX('$\\rho$, testing intensity (1/day per 1000)'))
)

p_check <- (p_temp
             + geom_contour_filled(breaks=brks_vec)
             ## + facet_grid(theta_w~theta_c, labeller=label_special)
             + scale_x_continuous(expand=expansion(c(0,0)), breaks=inv_omega_fast_brks,limits=c(inv_omega_fast_min,inv_omega_fast_max))
             ## scale by 1000, no digits after decimal
             # + scale_y_continuous(expand=expansion(c(0,0)), n.breaks=3)
             + scale_y_continuous(expand=expansion(c(0,0)), breaks=rho_fast_brks,
                                  labels=scales::number_format(scale=1000,accuracy=1))
             + scale_fill_viridis_d(name=parse(text="Delta"),drop=FALSE,
                                    labels=hack_breaks)
)
print(p_check)
