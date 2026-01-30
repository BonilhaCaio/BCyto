# csbSTATS ----
#' Analyse and visualize biostatistics data
#'
#' @description
#' The csbSTATS R package is an open-source project that provides a
#' user-friendly, interactive UI for exploratory data analysis,
#' visualization, and statistical testing.
#'
#' @docType package
#' @name csbSTATS
NULL

#' Run the csbSTATS application
#'
#' @description
#' runCsbSTATS() launches the csbSTATS Shiny application.
#'
#' @details
#' Use runCsbSTATS() at the R console to initiate the csbSTATS UI.
#'
#' @export
#'
runCsbSTATS <- function() {
  shinyApp(
    ui = ui,
    server = server,
    onStart = function() {
      shinyjs::useShinyjs()
    }
  )
}

# Imports ----------------------------------------------------------------
#' @import shiny
#' @import shinyjs
#' @import ggplot2
#' @import DT
#' @import jsonlite
#' @import colourpicker
#' @import fontawesome

#' @importFrom stats aov
#' @importFrom stats kruskal.test
#' @importFrom stats pairwise.wilcox.test
#' @importFrom stats shapiro.test
#' @importFrom stats sd
#' @importFrom stats t.test

#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom utils write.table
#' @importFrom stats TukeyHSD wilcox.test setNames
#' @importFrom utils combn

#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom base64enc dataURI

NULL

css_tabs <- "
body {
  background-color: #ffffff;
}

.nav-tabs {
  border-bottom: 1px solid #000000;
}

.nav-tabs > li > a {
  color: #000000 !important;
  border: 1px solid transparent;
}

.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover,
.nav-tabs > li.active > a:focus {
  color: #000000 !important;
  border: 1px solid #000000;
  border-bottom-color: transparent;
}
"

ui <- fluidPage(
  useShinyjs(),

  tags$head(
    tags$style(HTML("
  .tab-title-text {
    font-size: 14px;
    font-weight: 600;
  }
")),
    tags$style("#about {border-color:white; font-size:12px}"),
    tags$style(HTML("
    .instruction { font-size: 15px; margin-bottom: 10px; }
    .group-box {
      background-color: #ffffff;
      border: 1px solid #ddd;
      padding: 10px;
      margin-top: 10px;
    }
  ")),
    tags$style(HTML(css_tabs)),
    tags$style(HTML("
  .btn-primary {
    background-color: #000000 !important;
    border-color: #000000 !important;
    color: #ffffff !important;
  }

  .btn-primary:hover,
  .btn-primary:focus,
  .btn-primary:active {
    background-color: #000000 !important;
    border-color: #000000 !important;
    color: #ffffff !important;
  }
")),
    tags$style(HTML("
  .nav-tabs > li > a.disabled {
    pointer-events: none;
    color: #bdbdbd !important;
    cursor: not-allowed;
  }
"))
  ),

  headerPanel(
    fluidRow(
      column(
        width = 3,
        tags$div(
          style = "display:flex; align-items:center; gap:8px;",
          tags$img(
            src = base64enc::dataURI(
              file = system.file("icons", "Logo.png", package = "csbSTATS"),
              mime = "image/png"
            ),
            height = "50px"
          ),
          tags$span(
            "csbSTATS",
            style = "font-weight:600; font-size:30px;"
          )
        )
      ),
      column(
        width=2, offset=7, align="right",
        actionButton(inputId="about", label="About csbSTATS")
      )
    ),
  ),

  tabsetPanel(
    id = "main_tabs",

    ## TAB 1
    tabPanel(
      title = tagList(icon("cloud-upload-alt"), "File"),
      value = "file_tab",
      div(
        style = "
      display: flex;
      flex-direction: column;
      justify-content: center;
      align-items: center;
      height: 70vh;
    ",
        fileInput(
          "csv_file",
          label = NULL,
          buttonLabel = tagList(icon("file-upload"), "Load excel file"),
          placeholder = "Select a .csv file",
          accept = ".csv"
        ),
        br(),
        textOutput("file_status")
      )
    ),

    ## TAB 2
    tabPanel(
      title = tagList(icon("table"), "Data selection"),
      value = "selection_tab",
      br(),
      fluidRow(
        column(4, uiOutput("group_name_ui")),
        column(4, actionButton("add_group", "Assign group")),
        column(4, actionButton("finalize_groups", "Proceed", class = "btn-primary"))
      ),
      br(),
      div(class = "group-box",
          div(
            icon("list-ul"),
            "Defined groups",
            class = "tab-title-text"
          ),
          uiOutput("group_list")
      ),
      br(),
      DTOutput("data_table"),
    ),

    ## TAB 3
    tabPanel(
      title = tagList(icon("chart-bar"), "Plot"),
      value = "plot_tab",
      br(),
      fluidRow(

        ## LEFT COLUMN
        column(
          3,
          div(
            style = "border: 1px solid #000000; padding: 10px;",

            selectInput(
              "plot_type",
              tagList(icon("chart-column"), "Type"),
              c("Bar plot", "Box plot")
            ),
            selectInput(
              "error_type",
              tagList(icon("ruler-vertical"), "Error bars"),
              c("SD", "SEM")
            ),
            checkboxInput(
              "show_points",
              tagList("Show individual points"),
              TRUE
            ),
            br(),
            selectInput(
              "x_label_angle",
              tagList("Group name orientation"),
              choices = c("Horizontal" = 0, "Inclined" = 45, "Vertical" = 90),
              selected = 0
            ),
            uiOutput("group_order_selector"),
            br(),
            uiOutput("group_colour_selectors"),
            br(),
            ## Y-axis label input + OK button
            fluidRow(
              column(
                8,
                textInput("y_axis_label", "Y-axis label", "")
              ),
              column(
                4,
                conditionalPanel(
                  condition = "input.y_axis_label && input.y_axis_label.length > 0",
                  actionButton("apply_y_label", "OK")
                ),
                conditionalPanel(
                  condition = "!input.y_axis_label || input.y_axis_label.length === 0",
                  tags$button(
                    "OK",
                    class = "btn btn-default",
                    disabled = "disabled",
                    style = "margin-top: 25px;"
                  )
                )
              )
            ),
            br(),
            actionButton("go_stats", "Proceed to Stats", class = "btn-primary")
          )
        ),
        ## MIDDLE COLUMN
        column(
          6,
          div(
            style = "
              display: flex;
              justify-content: center;
              align-items: center;
              height: 100%;
            ",
            div(
              style = "
                width: 100%;
              ",
              plotOutput("main_plot"),
              tags$style(HTML("
                .shiny-plot-output .xlab,
                .shiny-plot-output .legend-title {
                  display: none;
                }
              "))
            )
          )
        ),
        ## RIGHT COLUMN
        column(
          3,
          div(
            style = "border: 1px solid #000000; padding: 10px;",
            sliderInput(
              "plot_width",
              tagList(icon("arrows-left-right")),
              250, 600, 350
            ),
            sliderInput(
              "plot_height",
              tagList(icon("arrows-up-down")),
              250, 600, 400
            ),
            sliderInput(
              "line_width",
              tagList(icon("pen-ruler")),
              0.5, 2, 1
            ),
            sliderInput(
              "bar_width",
              tagList(icon("grip-lines-vertical")),
              0.6, 1, 0.7
            ),
            sliderInput(
              "dot_size",
              tagList(icon("circle")),
              3, 8, 4
            ),
            sliderInput(
              "text_size",
              tagList(icon("text-height")),
              10, 30, 16
            )
          )
        )
      )
    ),

    ## TAB 4
    tabPanel(
      title = tagList(icon("calculator"), "Stats"),
      value = "stats_tab",
      br(),
      fluidRow(
        ## LEFT COLUMN
        column(
          3,
          div(
            style = "border: 1px solid #000000; padding: 10px;",
            uiOutput("comparison_selector"),
            br(),
            actionButton(
              "back_to_selection",
              tagList(icon("table"), "Return to Data")
            ),
            br(),
            br(),
            actionButton(
              "back_to_plot",
              tagList(icon("chart-bar"), "Return to Plot")
            ),
            br(),
            br(),
            downloadButton("save_plot", "Save plot", class = "btn-primary")
          )
        ),
        ## MIDDLE COLUMN
        column(
          6,
          div(
            style = "
              display: flex;
              justify-content: center;
              align-items: center;
              height: 100%;
            ",
            div(
              style = "width: 100%;",
              plotOutput("stats_plot")
            )
          )
        ),
        ## RIGHT COLUMN
        column(
          3,
          div(
            style = "border: 1px solid #000000; padding: 10px;",

            sliderInput(
              "sig_bar_offset",
              tagList(icon("arrows-up-down"), "Plot to bar"),
              min = 0, max = 0.2, value = 0.06, step = 0.01
            ),
            sliderInput(
              "sig_text_offset",
              tagList(icon("arrows-up-down"), "Bar to symbol"),
              min = 0, max = 0.1, value = 0.03, step = 0.005
            ),
            sliderInput(
              "sig_text_size",
              "Significance symbol size",
              min = 15, max = 50, value = 20
            )
          )
        )
      ),
      br(),
      fluidRow(
        column(
          12,
          verbatimTextOutput("stats_text")
        )
      )
    )
  )
)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "group",
    "value",
    "x1",
    "x2",
    "y",
    "label"
  ))
}

server <- function(input, output, session) {

  observeEvent(input$about, {
    showModal(
      modalDialog(
        title = "csbSTATS",
        "This is an open-source project that provides an user-friendly UI for
        automatd statistical analysis in R.",
        br(),
        br(),
        strong("Project website:"),
        "https://github.com/BonilhaCaio/csbSTATS",
        footer = NULL,
        size = "m",
        easyClose = TRUE
      )
    )
  })

  rv <- reactiveValues(
    data = NULL,
    groups = list(),
    locked_cells = NULL,
    y_axis_label = ""
  )

  disable("selection_tab")
  disable("plot_tab")
  disable("stats_tab")
  disable("finalize_groups")

  observeEvent(input$csv_file, {
    req(input$csv_file)
    if (tolower(tools::file_ext(input$csv_file$name)) != "csv") return()
    rv$data <- read.csv(input$csv_file$datapath, check.names = FALSE)
    output$file_status <- renderText("File loaded successfully.")
    enable("selection_tab")
    updateTabsetPanel(session, "main_tabs", selected = "selection_tab")
  })

  observeEvent(input$apply_y_label, {
    rv$y_axis_label <- input$y_axis_label
  })

  output$data_table <- renderDT({
    req(rv$data)
    locked_js <- if (is.null(rv$locked_cells)) "[]" else
      toJSON(lapply(seq_len(nrow(rv$locked_cells)), function(i)
        as.numeric(rv$locked_cells[i, ])), keep_vec_names = FALSE)

    datatable(
      rv$data,
      selection = list(mode = "multiple", target = "cell"),
      options = list(
        pageLength = nrow(rv$data),
        dom = "t",
        rowCallback = JS(sprintf(
          "function(row, data, displayIndex) {
     var locked = %s;
     locked.forEach(function(cell) {
       if (cell[0] === displayIndex + 1) {
         $('td:eq(' + cell[1] + ')', row).css('color', '#bdbdbd');
       }
     });
   }", locked_js))
      )
    )
  })

  observeEvent(input$add_group, {
    req(input$data_table_cells_selected, input$group_name)
    rv$groups[[input$group_name]] <- input$data_table_cells_selected
    rv$locked_cells <- rbind(rv$locked_cells, input$data_table_cells_selected)
    updateTextInput(session, "group_name", value = "")
  })

  output$group_list <- renderUI({
    tagList(lapply(names(rv$groups), function(g) {
      cells <- rv$groups[[g]]
      values <- as.numeric(apply(cells, 1, function(x)
        rv$data[x[1], x[2]]))

      div(paste0(
        "- ",
        g,
        ". N=",
        length(values),
        " (",
        paste(values, collapse = ", "),
        ")"
      ))
    }))
  })

  observeEvent(input$finalize_groups, {
    enable("plot_tab")
    updateTabsetPanel(session, "main_tabs", selected = "plot_tab")
  })

  observeEvent(input$back_to_selection, {
    rv$groups <- list()
    rv$locked_cells <- NULL
    updateTextInput(session, "group_name", value = "")
    updateTabsetPanel(session, "main_tabs", selected = "selection_tab")
    disable("plot_tab")
    disable("stats_tab")

  })

  observeEvent(input$back_to_plot, {
    disable("stats_tab")
    updateTabsetPanel(session, "main_tabs", selected = "plot_tab")
  })

  plot_data <- reactive({
    req(length(rv$groups) > 0)

    do.call(rbind, lapply(names(rv$groups), function(g) {
      cells <- rv$groups[[g]]
      data.frame(
        group = g,
        value = as.numeric(apply(cells, 1, function(x)
          rv$data[x[1], x[2]]))
      )
    }))
  })

  output$group_order_selector <- renderUI({
    selectInput(
      "group_order",
      tagList(icon("exchange-alt")),
      choices = names(rv$groups),
      selected = names(rv$groups),
      multiple = TRUE
    )
  })

  output$group_colour_selectors <- renderUI({
    fluidRow(lapply(names(rv$groups), function(g)
      column(3, colourpicker::colourInput(
        paste0("col_", g), g, "#BDBDBD"))))
  })

  observeEvent(input$plot_type, {
    if (input$plot_type == "Bar plot") enable("error_type") else disable("error_type")
  }, ignoreInit = TRUE)

  output$group_name_ui <- renderUI({
    textInput(
      "group_name",
      paste0(
        "Group name (click on cells for Group ",
        length(rv$groups) + 1,
        ")"
      )
    )
  })

  output$main_plot <- renderPlot({
    df <- plot_data()
    df$group <- factor(df$group, levels = input$group_order)
    fills <- sapply(levels(df$group), function(g) input[[paste0("col_", g)]])
    names(fills) <- levels(df$group)
    err_fun <- if (input$error_type == "SD") sd else function(x) sd(x)/sqrt(length(x))
    p <- ggplot(df, aes(group, value, fill = group)) +
      scale_fill_manual(values = fills) +
      labs(x = NULL, fill = NULL, y = rv$y_axis_label) +
      theme_classic(base_size = input$text_size) +
      theme(
        axis.text.x = element_text(
          angle = as.numeric(input$x_label_angle),
          hjust = if (as.numeric(input$x_label_angle) == 0) 0.5 else 1
        )
      )
    if (input$plot_type == "Bar plot") {
      p <- p +
        stat_summary(fun = mean, geom = "bar", colour = "black",
                     linewidth = input$line_width, width = input$bar_width) +
        stat_summary(fun.data = function(x) {
          m <- mean(x); e <- err_fun(x)
          data.frame(y = m, ymin = m - e, ymax = m + e)
        }, geom = "errorbar", colour = "black",
        linewidth = input$line_width, width = input$bar_width / 2)
    } else {
      p <- p + geom_boxplot(colour = "black", linewidth = input$line_width,
                            width = input$bar_width, outlier.shape = NA)
    }
    if (input$show_points) {
      p <- p + geom_jitter(
        aes(fill = group),
        shape = 21,
        colour = "black",
        size = input$dot_size,
        stroke = input$line_width,
        width = 0.1,
        show.legend = FALSE
      )
    }
    p
  }, width = function() input$plot_width,
  height = function() input$plot_height)

  observeEvent(input$go_stats, {
    enable("stats_tab")
    updateTabsetPanel(session, "main_tabs", selected = "stats_tab")
  })

  observe({
    if (length(rv$groups) < 2) {
      disable("finalize_groups")
    } else {
      enable("finalize_groups")
    }
  })

  observe({
    if (is.null(input$data_table_cells_selected) ||
        nrow(input$data_table_cells_selected) < 2) {
      disable("add_group")
    } else {
      enable("add_group")
    }
  })

  observe({
    tabs <- c("file_tab", "selection_tab", "plot_tab", "stats_tab")

    lapply(tabs, function(tb) {
      shinyjs::runjs(sprintf(
        "$('a[data-value=\"%s\"]').addClass('disabled');", tb
      ))
    })

    shinyjs::runjs(sprintf(
      "$('a[data-value=\"%s\"]').removeClass('disabled');",
      input$main_tabs
    ))
  })

  output$comparison_selector <- renderUI({
    req(length(rv$groups) >= 2)

    cmb <- combn(unique(plot_data()$group), 2, simplify = FALSE)

    checkboxGroupInput(
      "shown_comparisons",
      tagList(icon("grip-lines-vertical"), "Comparisons"),
      choices = sapply(cmb, paste, collapse = " vs "),
      selected = sapply(cmb, paste, collapse = " vs ")
    )
  })

  output$stats_text <- renderPrint({
    req(length(rv$groups) >= 2)
    df <- plot_data()
    sh <- by(df$value, df$group, function(x) {
      if (length(x) < 3) return(NA)
      if (sd(x, na.rm = TRUE) == 0) return(NA)
      shapiro.test(x)$p.value
    })
    sh_vals <- unlist(sh)
    testable <- sh_vals[!is.na(sh_vals)]
    normal <- if (length(testable) == 0) {
      TRUE
    } else {
      all(testable > 0.05)
    }
    cat("Normality (Shapiro-Wilk): ",
        ifelse(normal, "normal distribution\n", "non-normal distribution\n"))
    if (length(unique(df$group)) == 2) {
      test <- if (normal)
        t.test(value ~ group, df)
      else
        wilcox.test(value ~ group, df, exact = FALSE)
      cat("\nOmnibus test: ",
          ifelse(normal, "t-test", "Wilcoxon test"),
          "\np-value: ", signif(test$p.value, 3),
          "\nSignificant: ", ifelse(test$p.value < 0.05, "YES", "NO"), "\n")
    } else {
      omni <- if (normal) aov(value ~ group, df) else kruskal.test(value ~ group, df)
      p.omni <- if (normal) summary(omni)[[1]][["Pr(>F)"]][1] else omni$p.value
      cat("\nOmnibus test: ",
          ifelse(normal, "One-way ANOVA", "Kruskal-Wallis"),
          "\np-value: ", signif(p.omni, 3),
          "\nSignificant: ", ifelse(p.omni < 0.05, "YES", "NO"), "\n")
      if (p.omni < 0.05) {
        cat("\nMultiple comparisons:\n")
        pw <- if (normal)
          TukeyHSD(omni)$group[, "p adj"] else
            pairwise.wilcox.test(df$value, df$group, p.adjust.method = "BH")$p.value
        if (normal) {
          for (nm in names(pw)) {
            comps <- strsplit(nm, "-")[[1]]
            p <- pw[nm]
            sym <- if (p < 0.0001) "****"
            else if (p < 0.001) "***"
            else if (p < 0.01) "**"
            else if (p < 0.05) "*"
            else "ns"
            cat(comps[1], "vs", comps[2],
                ": p =", ifelse(p < 0.0001,
                                "< 0.0001",
                                formatC(p, format = "f", digits = 4)),
                "(", sym, ")\n")
          }
        } else {
          for (i in rownames(pw)) for (j in colnames(pw)) {
            if (!is.na(pw[i, j])) {
              p <- pw[i, j]
              sym <- if (p < 0.0001) "****"
              else if (p < 0.001) "***"
              else if (p < 0.01) "**"
              else if (p < 0.05) "*"
              else "ns"
              cat(i, "vs", j,
                  ": p =", ifelse(p < 0.0001,
                                  "< 0.0001",
                                  formatC(p, format = "f", digits = 4)),
                  "(", sym, ")\n")
            }
          }
        }
      }
    }
  })

  output$stats_plot <- renderPlot({
    df <- plot_data()
    df$group <- factor(df$group, levels = input$group_order)
    x_map <- setNames(seq_along(levels(df$group)), levels(df$group))
    fills <- sapply(levels(df$group), function(g) input[[paste0("col_", g)]])
    names(fills) <- levels(df$group)
    err_fun <- if (input$error_type == "SD") sd else function(x) sd(x)/sqrt(length(x))
    p <- ggplot(df, aes(group, value, fill = group)) +
      scale_fill_manual(values = fills) +
      labs(x = NULL, fill = NULL, y = rv$y_axis_label) +
      theme_classic(base_size = input$text_size) +
      theme(
        axis.text.x = element_text(
          angle = as.numeric(input$x_label_angle),
          hjust = if (as.numeric(input$x_label_angle) == 0) 0.5 else 1
        )
      )
    ymax <- tapply(df$value, df$group, max)
    range_y <- diff(range(df$value))
    if (input$plot_type == "Bar plot") {
      p <- p +
        stat_summary(fun = mean, geom = "bar", colour = "black",
                     linewidth = input$line_width, width = input$bar_width) +
        stat_summary(fun.data = function(x) {
          m <- mean(x); e <- err_fun(x)
          data.frame(y = m, ymin = m - e, ymax = m + e)
        }, geom = "errorbar", colour = "black",
        linewidth = input$line_width, width = input$bar_width / 2)
    } else {
      p <- p + geom_boxplot(colour = "black", linewidth = input$line_width,
                            width = input$bar_width, outlier.shape = NA)
    }
    if (input$show_points) {
      p <- p + geom_jitter(
        aes(fill = group),
        shape = 21,
        colour = "black",
        size = input$dot_size,
        stroke = input$line_width,
        width = 0.1,
        show.legend = FALSE
      )
    }
    if (length(input$shown_comparisons)) {
      comps <- strsplit(input$shown_comparisons, " vs ")
      bars <- do.call(rbind, lapply(seq_along(comps), function(i) {
        data.frame(
          x1 = x_map[comps[[i]][1]],
          x2 = x_map[comps[[i]][2]],
          base_y = max(ymax[comps[[i]]])
        )
      }))
      bars$span <- abs(bars$x2 - bars$x1)
      bars <- bars[order(bars$span), ]
      bars$y <- bars$base_y +
        range_y * input$sig_bar_offset * seq_len(nrow(bars))
      sh <- by(df$value, df$group, function(x) {
        if (length(x) < 3) return(NA)
        if (sd(x, na.rm = TRUE) == 0) return(NA)
        shapiro.test(x)$p.value
      })
      sh_vals <- unlist(sh)
      testable <- sh_vals[!is.na(sh_vals)]
      normal <- if (length(testable) == 0) {
        TRUE
      } else {
        all(testable > 0.05)
      }
      if (length(levels(df$group)) == 2) {
        test <- if (normal)
          t.test(value ~ group, df)
        else
          wilcox.test(value ~ group, df, exact = FALSE)
        pvals <- setNames(test$p.value,
                          paste(levels(df$group)[1],
                                levels(df$group)[2],
                                sep = " vs "))
      } else {
        if (normal) {
          omni <- aov(value ~ group, df)
          pw <- TukeyHSD(omni)$group[, "p adj"]
        } else {
          pw <- pairwise.wilcox.test(
            df$value,
            df$group,
            p.adjust.method = "BH",
            exact = FALSE
          )$p.value
        }
        pvals <- c()
        if (normal) {
          for (nm in names(pw)) {
            comps <- strsplit(nm, "-")[[1]]
            pvals[paste(comps[1], comps[2], sep = " vs ")] <- pw[nm]
            pvals[paste(comps[2], comps[1], sep = " vs ")] <- pw[nm]
          }
        } else {
          for (i in rownames(pw)) for (j in colnames(pw)) {
            if (!is.na(pw[i, j])) {
              pvals[paste(i, j, sep = " vs ")] <- pw[i, j]
            }
          }
        }
      }
      bars$label <- if (length(input$shown_comparisons)) {
        sapply(input$shown_comparisons, function(comp) {
          p <- pvals[comp]

          if (length(p) == 0 || is.na(p)) return("ns")

          if (p < 0.0001) "****"
          else if (p < 0.001) "***"
          else if (p < 0.01) "**"
          else if (p < 0.05) "*"
          else "ns"
        })
      } else {
        character(0)
      }
      p <- p +
        geom_segment(data = bars,
                     aes(x = x1, xend = x2, y = y, yend = y),
                     inherit.aes = FALSE,
                     linewidth = input$line_width) +
        geom_text(data = bars,
                  aes(x = (x1 + x2) / 2,
                      y = y + range_y * input$sig_text_offset,
                      label = label),
                  inherit.aes = FALSE,
                  fontface = "bold",
                  vjust = 0,
                  size = input$sig_text_size / 3) +
        expand_limits(y = max(bars$y + range_y * (input$sig_text_offset + 0.05)))
    }
    p
  }, width = function() input$plot_width,
  height = function() input$plot_height)

  output$save_plot <- downloadHandler(
    filename = function() "stats_plot.png",
    content = function(file) {
      ggsave(file, plot = last_plot(),
             width = input$plot_width / 72,
             height = input$plot_height / 72)
    }
  )
}
