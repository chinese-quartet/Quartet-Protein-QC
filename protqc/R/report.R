# ---------------------------------------------------------------------------- #
#' @title Generate Quartet Proteomics report
#'
#' @description Use calculated Met result to generate report
#'
#' @param qc_result list
#' @param report_template character
#' @param report_dir character
#' @param report_name character
#'
#' @importFrom dplyr %>%
#' @importFrom flextable flextable
#' @importFrom flextable theme_vanilla
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom officer body_add_par
#' @importFrom flextable body_add_flextable
#' @importFrom officer body_add_gg
#' @importFrom officer body_add_break
#' @importFrom officer read_docx
#'
#' @examples
#' # 加载示例 Protein_result 对象
#' protein_qc_result_path <- system.file("extdata", "protein_qc_result_example", package = "protqc")
#' load(protein_qc_result_path)
#'
#' # 指定包内文档的路径
#' report_template <- system.file("extdata", "quartet_template.docx", package = "protqc")
#'
#' # 运行函数
#' generate_protein_report(qc_result = qc_result, report_template = report_template)
#'
#' @export
generate_protein_report <- function(qc_result,
                                    report_template,
                                    report_dir = NULL,
                                    report_name = NULL) {
  if (is.null(qc_result) || is.null(report_template)) {
    stop("All arguments (qc_result, report_template) are required.")
  }

  if (is.null(report_dir)) {
    path <- getwd()
    report_dir <- file.path(path, "output")
    dir.create(report_dir, showWarnings = FALSE)
  }

  ### 读取quarter报告模板并生成报告
  if (is.null(report_name)) {
    report_name <- "Quartet_protein_report.docx"
  }
  output_file <- file.path(report_dir, report_name)

  ### 创建Evaluate Metrics 表格

  ft1 <- flextable(qc_result$conclusion)
  ft1 <- ft1 %>%
    color(~ Performance == "Bad", color = "#B80D0D", ~Performance) %>%
    color(~ Performance == "Fair", color = "#D97C11", ~Performance) %>%
    color(~ Performance == "Good", color = "#70C404", ~Performance) %>%
    color(~ Performance == "Great", color = "#0F9115", ~Performance) %>%
    width(width = 1.25) %>%
    align(align = "center", part = "all") %>%
    bold(i = 7, part = "body")

  ### 绘制Total score 历史分数排名散点图
  qc_result$rank_table$Total_norm <- as.numeric(qc_result$rank_table$Total_norm)

  p_rank_scatter_plot <- ggplot(data = qc_result$rank_table) +
    # 添加四个区域
    geom_rect(aes(xmin = 1, xmax = 3.47, ymin = -Inf, ymax = Inf), fill = "#B80D0D", alpha = 0.08) +
    geom_rect(aes(xmin = 3.47, xmax = 4.19, ymin = -Inf, ymax = Inf), fill = "#D97C11", alpha = 0.08) +
    geom_rect(aes(xmin = 4.19, xmax = 6.12, ymin = -Inf, ymax = Inf), fill = "#70C404", alpha = 0.08) +
    geom_rect(aes(xmin = 6.12, xmax = 10, ymin = -Inf, ymax = Inf), fill = "#0F9115", alpha = 0.08) +
    # 添加基础点图层
    geom_point(aes(x = Total_norm, y = reorder(Batch, Total_norm))) +
    # 突出显示 "QUERIED DATA" 对应的点
    geom_point(
      data = subset(qc_result$rank_table, Batch == "QUERIED DATA"),
      aes(x = Total_norm, y = reorder(Batch, Total_norm)),
      color = "orange", size = 3
    ) +
    # 自定义x轴刻度
    scale_x_continuous(breaks = c(1, 3.47, 4.19, 6.12, 10)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      x = " ",
      y = " ",
      title = "Total Score"
    )

  #### 设置输出文本

  ###### 第一部分
  ### Assessment Summary
  text_1 <- "The performance of the submitted data will be graded as Bad, Fair, Good, or Great based on the ranking by comparing the total score with the historical datasets. The total score is the harmonic mean of the scaled values of the number of features, missing percentage, absolute correlation, coefficient of variation, Signal-to-Noise Ratio (SNR), and relative correlation with reference datasets (RC)."
  ### Four levels of performance
  text_1_sup_1 <- "Based on the scaled total score, the submitted data will be ranked together with all Quartet historical datasets. The higher the score, the higher the ranking. After this, the performance levels will be assigned based on their ranking ranges."
  text_1_sup_2 <- "· Bad - the bottom 20%."
  text_1_sup_3 <- "· Fair - between bottom 20% and median 50%."
  text_1_sup_4 <- "· Good - between median 50% and top 20%."
  text_1_sup_5 <- "· Great - the top 20%."

  #### 第二部分 Quality control metric

  ### Performance Score
  text_2 <- "Scores of evaluation metrics for the current batch and all historical batches assessed. Please note that the results shown here are scaled values for all batches in each metric. The name of your data is Queried_Data."
  ### Signal-to-Noise Ratio
  text_3 <- "SNR is established to characterize the power in discriminating multiple groups. The PCA plot is used to visualize the metric."
  ### Correlation with Reference Datasets
  text_4 <- "Relative correlation with reference datasets metric which was representing the numerical consistency of the relative expression profiles."

  ### Method
  supplementary_info_1 <- "The QC pipeline starts from the expression profiles at peptide/protein levels, and enables to calculate 6 metrics. A Total score is the geometric mean of the linearly normalized values of these metrics."
  supplementary_info_1_1 <- "We expect as many proteins (mapped to gene symbols) as possible for downstreaming analyses."
  supplementary_info_1_2 <- "Too many missing values interfere with comparability. This metric is calculated globally."
  supplementary_info_1_3 <- "A CV value is calculated to indicate the dispersion within replicates feature by feature."
  supplementary_info_1_4 <- "Pearson correlation reflects overall reproducibility within replicates. We calculate correlation coefficients between each two replicates within each biological sample (D5, D6, F7, M8), and take the median as the final value for absolute correlation."
  supplementary_info_1_5 <- "SNR is established to characterize the ability of a platform or lab or batch, which is able to distinguish intrinsic differences among distinct biological sample groups ('signal') from variations in technical replicates of the same sample group ('noise')."
  supplementary_info_1_6 <- "RC is used for assessment of quantitative consistency with the reference dataset at relative levels. For shotgun proteomics, quantitation at peptide levels is theoretically more reliable. Therefore, the reference dataset is established by benchmarking the relative expression values (log2FCs), for each peptide sequence of each sample pair (D5/D6, F7/D6, M8/D6), in historical datasets at peptide levels. We calculate relatively qualified (satisfied with thresholds of p < 0.05) log2FCs of the queried data, for overlapped peptides with the reference dataset, as the input for the assessment of quantitative consistency. Then RC value is Pearson correlation coefficient between the test dataset and the reference dataset."


  ### Reference
  supplementary_info_ref1 <- "1. Zheng, Y. et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology 1–17 (2023)."
  supplementary_info_ref2 <- "2. Tian, S. et al. Quartet protein reference materials and datasets for multi-platform assessment of label-free proteomics. Genome Biology 24, 202 (2023)."

  ### Contact us
  supplementary_info_2_1 <- "Fudan University Pharmacogenomics Research Center"
  supplementary_info_2_2 <- "Project manager: Quartet Team"
  supplementary_info_2_3 <- "Email: quartet@fudan.edu.cn"

  ### Disclaimer
  supplementary_info_3 <- "This quality control report is only for this specific test data set and doesn’t represent an evaluation of the business level of the sequencing company. This report is only used for scientific research, not for clinical or commercial use. We don’t bear any economic and legal liabilities for any benefits or losses (direct or indirect) from using the results of this report."

  read_docx(report_template) %>%
    ## 添加报告标题
    body_add_par(value = "Quartet Report for Proteomics", style = "heading 1") %>%
    ## 第一部分，Assessment Summary
    body_add_par(value = "Assessment Summary", style = "heading 2") %>%
    body_add_flextable(ft1) %>%
    body_add_break() %>%
    ### 第二部分 Quality control metric
    body_add_par(value = "Quality Control Metric", style = "heading 2") %>%
    body_add_par(value = supplementary_info_1, style = "Normal") %>%
    body_add_par(value = "Number of features:", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_1, style = "Normal") %>%
    body_add_par(value = "Missing percentage (%):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_2, style = "Normal") %>%
    body_add_par(value = "Coefficient of variantion (CV %):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_3, style = "Normal") %>%
    body_add_par(value = "Absolute Correlation:", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_4, style = "Normal") %>%
    body_add_par(value = "Signal-to-Noise Ratio (SNR):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_5, style = "Normal") %>%
    body_add_par(value = "Relative Correlation with Reference Datasets (RC):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_6, style = "Normal") %>%
    body_add_par(value = "Total Score:", style = "heading 3") %>%
    body_add_par(value = text_1, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    body_add_par(value = "Performance Category:", style = "heading 3") %>%
    body_add_par(value = text_1_sup_1, style = "Normal") %>%
    body_add_par(value = text_1_sup_2, style = "Normal") %>%
    body_add_par(value = text_1_sup_3, style = "Normal") %>%
    body_add_par(value = text_1_sup_4, style = "Normal") %>%
    body_add_par(value = text_1_sup_5, style = "Normal") %>%
    ### 排名散点图
    body_add_par(value = "Performance Score", style = "heading 2") %>%
    body_add_gg(value = p_rank_scatter_plot, style = "centered") %>%
    body_add_par(value = text_2, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ## 信噪比
    body_add_par(value = "Signal-to-Noise Ratio", style = "heading 2") %>%
    body_add_gg(qc_result$results$snr_results$snr_plot, style = "centered") %>%
    body_add_par(value = text_3, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ## RC
    body_add_par(value = "Correlation with Reference Datasets", style = "heading 2") %>%
    body_add_gg(qc_result$results$cor_results$cor_plot, style = "centered") %>%
    body_add_par(value = text_4, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ### 附加信息
    body_add_par(value = "Supplementary Information", style = "heading 2") %>%
    # body_add_par(value = "Method", style = "heading 3") %>%
    # body_add_par(value = supplementary_info_1_1, style = "Normal") %>%
    # body_add_par(value = supplementary_info_1_2, style = "Normal") %>%
    # body_add_par(value = supplementary_info_1_3, style = "Normal") %>%

    body_add_par(value = "Reference", style = "heading 3") %>%
    body_add_par(value = supplementary_info_ref1, style = "Normal") %>%
    body_add_par(value = supplementary_info_ref2, style = "Normal") %>%
    body_add_par(value = "Contact us", style = "heading 3") %>%
    body_add_par(value = supplementary_info_2_1, style = "Normal") %>%
    body_add_par(value = supplementary_info_2_2, style = "Normal") %>%
    body_add_par(value = supplementary_info_2_3, style = "Normal") %>%
    body_add_par(value = "Disclaimer", style = "heading 3") %>%
    body_add_par(value = supplementary_info_3, style = "Normal") %>%
    ## 输出文件
    print(target = output_file)
}
