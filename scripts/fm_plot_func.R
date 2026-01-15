# Helper functions for plotting, like my_theme, mylegend



# Load functions

# Helper for Box
n_fun <- function(x) {
  return(data.frame(y = median(x) * 1.3, label = paste0("~italic(n)", " == ", length(x))))
}

# without x labels
my_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 16),
    legend.position = "none"
  )


my_theme1 <- my_theme +
  theme(
    axis.text.x = element_text(size = 14, face = "bold")
  )



my_theme_pred <- my_theme + theme(
  axis.text.x = element_text(size = 16, face = "bold.italic"),
  axis.text.y = element_text(size = 16, face = "bold"),
  legend.position = "none"
)

my_theme_pred_quant <- my_theme + theme(
  axis.text.x = element_text(size = 16, face = "bold"),
  axis.text.y = element_blank(), # removes the y‐axis labels
  axis.ticks.y = element_blank(), # removes the y‐axis tick marks
  legend.position = "none"
)




# Create the significances legend
library(cowplot)


levels <- c("A", "S1", "B", "C")
num <- c(5, 10, 15, 20)
ymin <- c(0, 0, 0, 0)
ymax <- c(1, 2, 3, 4)

Legend_DF <- data.frame(levels, num, ymin, ymax)

# Define significance-based colors and labels
sig_levels <- c("A" = "NS", "S1" = "*", "B" = "**", "C" = "***")
sig_colors <- c("A" = "black", "S1" = "cyan", "B" = "#ee9a00", "C" = "red")

plot <- ggplot(Legend_DF, aes(x = levels, y = num, colour = levels)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), size = 1.1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm"),
    legend.text = element_text(size = 20)
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_color_manual(
    values = sig_colors,
    breaks = names(sig_levels),
    labels = paste0(sig_levels, "     ") # Adding spacing for alignment
  )

mylegend <- get_legend(plot)


levels <- c("A", "B", "C")
num <- c(5, 15, 20)
ymin <- c(0, 0, 0)
ymax <- c(1, 3, 4)
Legend_DF <- data.frame(levels, num, ymin, ymax)
plot2 <- ggplot(Legend_DF, aes(x = levels, y = num, colour = levels)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), size = 1.1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm"),
    legend.text = element_text(size = 20)
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_color_manual(
    values = c("black", "#ee9a00", "red"),
    breaks = c("A", "B", "C"),
    labels = c("NS     ", "**     ", "***")
  )
mylegend23 <- get_legend(plot2)

# Delete unnecessary objects
rm(Legend_DF)
rm(plot)
rm(plot2)

# Final Model Functions ####



library(ggeffects)

# model = Stress.PP.fit.full

get_model_var_levels <- function(model) {
  # Extract the data frame used by the model
  mf <- model.frame(model)

  # Identify which columns are factors
  factor_vars <- names(mf)[vapply(mf, is.factor, logical(1))]

  numeric_vars <- names(mf)[vapply(mf, is.numeric, logical(1))]
  sprintf("Factor variables in model: %s", paste(factor_vars, collapse = ", "))
  sprintf("Numeric variables in model: %s", paste(numeric_vars, collapse = ", "))
  # For each factor, get its levels
  levels_list <- lapply(mf[factor_vars], levels)
  return(levels_list)
}


# var_levels <- get_model_var_levels(Stress.PP.fit.full)
# names(var_levels)



get_term_colors <- function(model, term, default = "gray") {
  # 1) model frame & coefficient table
  mf <- model.frame(model)
  cf <- as.data.frame(summary(model)$coefficients)
  # ensure rownames(cf) are the coef names e.g. "QTypeA", "QTypeW", etc.

  # 2) add a color column based on p‐value
  cf$color <- cut(
    cf[, "Pr(>|z|)"],
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("red", "orange", "cyan", "black")
  )

  # 3) if term is not a factor, just return its single color
  if (!term %in% names(mf) || !is.factor(mf[[term]])) {
    if (term %in% rownames(cf)) {
      return(setNames(as.character(cf[term, "color"]), term))
    } else {
      #return(setNames("default", term))
      return(setNames("black", term))
    }
  }

  # 4) otherwise, loop over each level
  lvls <- levels(mf[[term]])
  cols <- sapply(lvls, function(lvl) {
    nm <- paste0(term, lvl)
    if (nm %in% rownames(cf)) {
      as.character(cf[nm, "color"])
    } else {
      default
    }
  }, USE.NAMES = TRUE)

  return(cols)
}


# get_term_colors(Stress.PP.fit.full, "QType")
# get_term_colors(Stress.PP.fit.full, "Gender")
# get_term_colors(Stress.PP.fit.full, "Grade")
# get_term_colors(Stress.PP.fit.full, "SAI")



# finalModel=Stress.PP.fit.full

fm_plot_func <- function(finalModel, term = "QType", ymin = 0.0, ymax = 1.0) {
  
  
  # Allows to plot the model with different variables, QType, Gender, Grade, SAI
  # Use the function `fm_plot_func(Stress.PP.fit.full, term = "QType", ymin = 0.25, ymax = 0.7)` 
  # to plot the model with QType variables
  
  
  # Get the coefficient matrix
  coeff_matrix <- as.data.frame(coef(summary(finalModel)))


  # Round all numeric columns to 3 significant digits
  coeff_matrix[] <- lapply(coeff_matrix, function(x) {
    if (is.numeric(x)) signif(x, 3) else x
  })

  names(coeff_matrix)
  # Add  'Significance' column based on p-values
  coeff_matrix$Significance <- ifelse(
    coeff_matrix$`Pr(>|z|)` < 0.001, "***",
    ifelse(coeff_matrix$`Pr(>|z|)` < 0.01, "**",
      ifelse(coeff_matrix$`Pr(>|z|)` < 0.05, "*", "")
    )
  )

  coeff_matrix$color <- ifelse(
    coeff_matrix$`Pr(>|z|)` < 0.001, "red",
    ifelse(coeff_matrix$`Pr(>|z|)` < 0.01, "orange",
      ifelse(coeff_matrix$`Pr(>|z|)` < 0.05, "cyan", "black")
    )
  )


  # Suppose your model is named finalModel
  var_levels <- get_model_var_levels(finalModel)
  names(var_levels)


  ylims <- c(ymin, ymax) # y-axis limits for the plots


  if (term == "QType") {
    print("QType")


    group_levels <- var_levels[[term]]
    group_colors <- get_term_colors(finalModel, term)


    # strip off the names, leaving c("gray","orange","red")
    group_colors <- unname(group_colors)


    # Generate ggpredict output
    final_qt_plot <- ggpredict(finalModel, terms = c("QType")) %>%
      mutate(group = as.factor(c("1", "2", "3"))) %>%
      drop_na() %>%
      gather(key = "key", value = "value", conf.low, conf.high) %>%
      mutate(yminn = min(value), ymaxx = max(value)) %>%
      ggplot() +
      geom_line(aes(x = x, y = value, color = group), size = 2) +
      geom_point(aes(x = x, y = predicted, color = group), size = 3) +
      # scale_y_continuous(limits = c(0.3, 0.7)) +
      # scale_y_continuous(limits = c(yminn[1], ymaxx[1])) +
      scale_y_continuous(limits = ylims,
                         breaks = seq(ymin, ymax, by = 0.2),
                         labels = scales::label_number(accuracy = 0.1) 
                         ) + 
      
      my_theme_pred +
      theme(plot.title = element_text(size = 15)) +
      scale_color_manual(values = c(group_colors)) + # Use dynamically assigned colors
      labs(
        x = "",
        y = expression(bold(paste(bolditalic(P(S)), "")))
      )

    return(final_qt_plot)
  } else if (term == "Gender") {
    term <- "Gender"
    print("Gender")

    group_levels <- var_levels[[term]]
    group_colors <- get_term_colors(finalModel, term)


    # strip off the names, leaving c("gray","orange","red")
    group_colors <- unname(group_colors)

    pp.gender.pm <- ggpredict(finalModel,
      terms = c("Gender")
    ) %>%
      mutate(group = as.factor(c("1", "2"))) %>%
      drop_na() %>% # check # of levels you've
      #  gather data using conf.low and conf.high
      gather(key = "key", value = "value", conf.low, conf.high) %>%
      ggplot() +
      geom_line(aes(x = x, y = value, color = group), size = 2) +
      geom_point(aes(x = x, y = predicted, color = group, size = 3)) +
      theme_bw() +
      # scale_color_manual(values = c("gray", "#ee9a00")) +
      scale_color_manual(values = c(group_colors)) + # Use dynamically assigned colors
      scale_y_continuous(limits = ylims, guide = "none",  # "none" removes the y-axis ticks and labels
                         breaks = seq(ymin, ymax, by = 0.2),
                         labels = scales::label_number(accuracy = 0.1) ) + 
      my_theme_pred +
      labs(x = "", y = "", title = "Gender")

    return(pp.gender.pm)
  } else if (term == "Grade") {
    print("Grade")

    group_levels <- var_levels[[term]]
    group_colors <- get_term_colors(finalModel, term)


    # strip off the names, leaving c("gray","orange","red")
    group_colors <- unname(group_colors)

    # Correctness
    pp_grade.pm <- ggpredict(finalModel, "Grade") %>%
      mutate(group = as.factor(c("1", "2"))) %>%
      drop_na() %>% # check # of levels you've
      #  gather data using conf.low and conf.high
      gather(key = "key", value = "value", conf.low, conf.high) %>%
      ggplot() +
      geom_line(aes(x = x, y = value, color = group), size = 2) +
      geom_point(aes(x = x, y = predicted, color = group, size = 3)) +
      theme_bw() +
      # scale_color_manual(values = c("gray", "black")) +
      scale_color_manual(values = c(group_colors)) + # Use dynamically assigned colors
      #scale_y_continuous(limits = ylims, guide = "none") + # "none" removes the y-axis ticks and labels
      scale_y_continuous(limits = ylims, guide = "none",  # "none" removes the y-axis ticks and labels
                         breaks = seq(ymin, ymax, by = 0.2),
                         labels = scales::label_number(accuracy = 0.1) ) + 
      my_theme_pred_quant +
      labs(x = "", y = "", title = "Grade")

    return(pp_grade.pm)
  } else if (term == "SAI") {
    print("SAI")

    group_levels <- var_levels[[term]]
    group_colors <- get_term_colors(finalModel, term)

    group_colors <- unname(group_colors)

    line_color <- group_colors[1] # Use the first color for the line

    sai.pm <- ggpredict(finalModel, "SAI") %>%
      ggplot(aes(x, predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
      geom_line(colour = line_color, size = 2) + # use dynamically assigned color
      #scale_y_continuous(limits = ylims, guide = "none") + # "none" removes the y-axis ticks and labels
      scale_y_continuous(limits = ylims, guide = "none",  # "none" removes the y-axis ticks and labels
                         breaks = seq(ymin, ymax, by = 0.2),
                         labels = scales::label_number(accuracy = 0.1) ) + 
      geom_vline(
        xintercept = mean(Qlevel$SAI.Score),
        linetype = "dashed",
        color = "gray",
        size = 1
      ) +
      my_theme_pred +
      theme(axis.text.x = element_text(size = 14, face = "bold")) +
      labs(x = "", y = "", title = "SAI Score")

    return(sai.pm)
  } else {
    print("Unknown term")
  }




  print("ZAMTING ZONG")
}
