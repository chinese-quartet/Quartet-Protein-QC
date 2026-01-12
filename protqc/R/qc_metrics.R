#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data --------------------------------
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample

  # Replace zero by NA -----------------------
  d[d == 0] <- NA

  # Statistics: number of features -----------
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)

  # Statistics: missing percentage -----------
  d_all_num <- nrow(d) * (ncol(d) - 1)
  d_missing_num <- length(which(is.na(d)))
  prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)

  # Check if replicates available ------------
  samples <- table(s)
  rep_samples <- samples[samples > 1]
  rep_num <- length(rep_samples)
  if (rep_num == 0) {
    stop("No replicates are available.")
  } else {
    # Calculating: absolute correlation ------
    d_mtx <- d[, 2:ncol(d)]
    d_cortest <- corr.test(d_mtx, method = "pearson", adjust = "fdr")
    d_pmtx <- d_cortest$p
    d_cormtx <- d_cortest$r
    d_cormtx[d_pmtx > 0.05] <- 0
    d_cordf <- melt(d_cormtx)
    d_cordf <- d_cordf[d_cordf$Var2 != d_cordf$Var1, ]
    d_cordf <- merge(d_cordf, m, by.x = "Var1", by.y = "library")
    d_cordf <- merge(d_cordf, m, by.x = "Var2", by.y = "library")
    d_cordf <- d_cordf[d_cordf$sample.x == d_cordf$sample.y, ]
    cor_value <- median(d_cordf$value)
    stat_acor <- round(cor_value, 3)
  }

  # Calculating: CV --------------------------
  d_long <- melt(d)
  d_long <- na.omit(d_long)
  if (length(d_long$value[d_long$value < 0])) {
    d_long$value <- 2^(d_long$value)
    message("There are negative values. Default to log2-transformed values.")
  }
  d_long <- merge(d_long, m, by.x = "variable", by.y = "library")
  colnames(d_long) <- c("library", "feature", "value", "sample")
  d_cv <- aggregate(
    value ~ feature + sample,
    data = d_long,
    FUN = function(x) sd(x) / mean(x)
  )
  stat_cv <- round(median(d_cv$value, na.rm = T) * 100, 3)

  # Output -----------------------------------
  stat_all <- c(stat_num, stat_missing, stat_acor, stat_cv)

  return(stat_all)
}

#' Get color mapping for samples
#' @param samples Vector of sample names
#' @return Named vector of colors
#' @export
get_sample_colors <- function(samples) {
  # Define fixed color palette
  color_palette <- c(
    "D5" = "#4CC3D9", # Blue
    "D6" = "#7BC8A4", # Green
    "F7" = "#FFC65D", # Yellow
    "M8" = "#F16745" # Red
  )

  # Return colors only for existing samples
  available_colors <- color_palette[samples]
  return(available_colors)
}

#' Calculating SNR value; Plotting a PCA panel
#' @param expr_dt A expression profile (at protein level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @import stats
#' @import utils
#' @importFrom rlang :=
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 ggsave
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # Load data --------------------------------------
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])

  # Replace NA by zero -----------------------------
  expr_df[is.na(expr_df)] <- 0

  # Label the grouping info ------------------------
  ids <- colnames(expr_df)
  group <- meta_dt$sample
  ids_group_mat <- data.table(id = ids, group = group)

  # PCA --------------------------------------------
  expr_df_t <- t(expr_df)

  # Remove constant/zero variance features before PCA
  col_vars <- apply(expr_df_t, 2, var, na.rm = TRUE)
  zero_var_cols <- which(col_vars == 0 | is.na(col_vars))

  # 记录原始特征数和过滤后特征数
  n_features_original <- ncol(expr_df_t)
  n_features_removed <- length(zero_var_cols)

  if (n_features_removed > 0) {
    message(sprintf(
      "Removing %d features with zero variance before PCA analysis.",
      n_features_removed
    ))
    expr_df_t <- expr_df_t[, -zero_var_cols, drop = FALSE]
  }

  n_features_used <- ncol(expr_df_t)

  # Check if enough features remain
  if (n_features_used < 2) {
    stop("Not enough features with non-zero variance for PCA analysis. At least 2 features required.")
  }

  pca_prcomp <- prcomp(expr_df_t, retx = T, scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$sample_id <- rownames(pcs)
  pcs$sample <- meta_dt$sample

  # Calculating: SNR -------------------------------
  dt_perc_pcs <- data.table(
    PCX = 1:nrow(pcs),
    Percent = summary(pca_prcomp)$importance[2, ],
    AccumPercent = summary(pca_prcomp)$importance[3, ]
  )

  dt_dist <- data.table(
    id_a = rep(ids, each = length(ids)),
    id_b = rep(ids, time = length(ids))
  )

  dt_dist$group_a <- ids_group_mat[match(dt_dist$id_a, ids_group_mat$id)]$group
  dt_dist$group_b <- ids_group_mat[match(dt_dist$id_b, ids_group_mat$id)]$group

  dt_dist[, type := ifelse(id_a == id_b, "Same",
    ifelse(group_a == group_b, "Intra", "Inter")
  )]

  dt_dist[, dist := (dt_perc_pcs[1]$Percent * (pcs[id_a, 1] - pcs[id_b, 1])^2 +
    dt_perc_pcs[2]$Percent * (pcs[id_a, 2] - pcs[id_b, 2])^2)]

  dt_dist_stats <- dt_dist[, list(avg_dist = mean(dist)), by = list(type)]
  setkey(dt_dist_stats, type)
  signoise <- dt_dist_stats["Inter"]$avg_dist / dt_dist_stats["Intra"]$avg_dist
  signoise_db <- round(10 * log10(signoise), 3)

  # Plot -------------------------------------------
  if (plot) {
    # Get unique samples from metadata
    unique_samples <- unique(meta_dt$sample)
    colors_custom <- get_sample_colors(unique_samples)

    text_custom_theme <- element_text(
      size = 16,
      face = "plain",
      color = "black",
      hjust = 0.5
    )
    scale_axis_x <- c(min(pcs$PC1), max(pcs$PC1))
    scale_axis_y <- c(min(pcs$PC2), max(pcs$PC2))

    pc1_prop <- summary(pca_prcomp)$importance[2, 1]
    pc2_prop <- summary(pca_prcomp)$importance[2, 2]
    text_axis_x <- sprintf("PC1(%.2f%%)", pc1_prop * 100)
    text_axis_y <- sprintf("PC2(%.2f%%)", pc2_prop * 100)
    limit_x <- c(1.1 * scale_axis_x[1], 1.1 * scale_axis_x[2])
    limit_y <- c(1.1 * scale_axis_y[1], 1.1 * scale_axis_y[2])


    # 修改图表标题，显示实际使用的特征数
    p_title <- paste("SNR = ", signoise_db, sep = "")
    p_subtitle <- paste("(Proteins used in PCA = ", n_features_used,
      "/", n_features_original, ")",
      sep = ""
    )
    # p_title <- paste("SNR = ", signoise_db, sep = "")
    # p_subtitle <- paste("(Number of proteins = ", nrow(expr_dt), ")", sep = "")
    p <- ggplot(pcs, aes(x = .data$PC1, y = .data$PC2)) +
      geom_point(aes(color = sample), size = 8) +
      theme_few() +
      theme(
        plot.title = text_custom_theme,
        plot.subtitle = text_custom_theme,
        axis.title = text_custom_theme,
        axis.text = text_custom_theme,
        legend.title = text_custom_theme,
        legend.text = element_text(size = 16, color = "gray40")
      ) +
      labs(
        x = text_axis_x,
        y = text_axis_y,
        title = p_title,
        subtitle = p_subtitle
      ) +
      scale_color_manual(values = colors_custom) +
      scale_x_continuous(limits = limit_x) +
      scale_y_continuous(limits = limit_y) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }

  pc_num <- ncol(pcs)
  output <- data.table(pcs[, c((pc_num - 1):pc_num, 1:(pc_num - 2))])

  # Save & Output ----------------------------------
  if (!is.null(output_dir)) {
    if (plot) {
      output_dir_final1 <- file.path(output_dir, "pca_plot.png")
      ggsave(output_dir_final1, p, width = 6, height = 5.5)
    }
    output_dir_final2 <- file.path(output_dir, "pca_table.tsv")
    write.table(output, output_dir_final2, sep = "\t", row.names = F)
  }

  # return(list(table = output, SNR = signoise_db)
  return(list(table = output, SNR = signoise_db, snr_plot = p))
}

#' Analysis: differential expression
#' @param expr A expression table file (at peptide level)
#' @param group The grouping info
#' @import stats
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export
dep_analysis <- function(expr, group) {
  dge <- DGEList(counts = expr)
  design <- model.matrix(~group)

  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot = F)
  fit <- lmFit(v, design)

  fit <- eBayes(fit)
  result <- topTable(fit, coef = ncol(design), sort.by = "logFC", number = Inf)
  result$Sequence <- rownames(result)
  result$Sequence.Number <- nrow(result)
  result$Sample1 <- levels(group)[1]
  result$Sample2 <- levels(group)[2]
  result$Sample.Pair <- paste(levels(group)[2], levels(group)[1], sep = "/")

  return(result)
}

#' Calculating RC value; Plotting a scatterplot
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @param show_sample_pairs if True, samples in plot will be labeled.
#' @import stats
#' @import utils
#' @importFrom rlang .data
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggthemes theme_few
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggsave
#' @export
qc_cor <- function(expr_dt, meta_dt,
                   output_dir = NULL, plot = FALSE, show_sample_pairs = TRUE) {
  # Load data ------------------------------------------------------

  # # 1. 加载内置的定量数据集 (Quantitative / RC)
  # utils::data("reference_dataset_quant", package = "protqc", envir = environment())
  # ref_raw <- reference_dataset_quant
  # 
  # # 2. 列名映射与提取
  # # 你的 CSV 列名: group, peptide_sequence, value
  # # 我们需要映射为: Sample.Pair, Sequence, logFC 以适配原有逻辑
  # if (!all(c("group", "peptide_sequence", "value") %in% colnames(ref_raw))) {
  #   stop("Built-in reference_dataset_quant is missing: group, peptide_sequence, value")
  # }
  # 
  # ref_dt <- data.frame(
  #   Sequence = ref_raw$peptide_sequence,
  #   Sample.Pair = ref_raw$group,
  #   logFC = ref_raw$value,
  #   stringsAsFactors = FALSE
  # )
  
  # 1. 回滚：使用原来的 reference_dataset
  utils::data("reference_dataset", package = "protqc", envir = environment())
  ref_dt <- reference_dataset
  
  # 3. 数据预处理 (保持原有逻辑)
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])
  expr_matrix <- as.matrix(expr_df)
  expr_matrix <- expr_matrix[, meta_dt$library]

  # Replace NA by zero ---------------------------------------------
  expr_matrix[is.na(expr_matrix)] <- 0

  # Check if negative values exist ----------------------------------
  if (length(expr_matrix[expr_matrix < 0])) {
    expr_matrix <- 2^(expr_matrix)
    message("There are negative values. Default to log2-transformed values.")
  }

  # Check the grouping info ----------------------------------------
  samples <- as.character(unique(meta_dt$sample))

  # Modified: Flexibly identify which samples exist
  if (length(samples) < 2) {
    stop("At least 2 different sample types are required for correlation analysis.")
  }

  # Check if D6 exists (reference sample)
  check_d6 <- "D6" %in% samples
  if (!check_d6) {
    warning("D6 (reference sample) is not available. Using the first sample as reference.")
    reference_sample <- samples[1]
    samples <- c(reference_sample, samples[!samples %in% reference_sample])
  } else {
    reference_sample <- "D6"
    samples <- c("D6", samples[!samples %in% "D6"])
  }

  pair_num <- length(samples) - 1

  # Analysis: Differential expression ------------------------------
  result_final <- c()
  for (j in 2:(pair_num + 1)) {
    sample_pair <- paste(samples[j], reference_sample, sep = "/")
    ref_tmp <- ref_dt[ref_dt$Sample.Pair %in% sample_pair, ]

    col1 <- which(meta_dt$sample %in% samples[j])
    col2 <- which(meta_dt$sample %in% reference_sample)

    e_tmp <- expr_matrix[, c(col1, col2)]
    e_tmp <- e_tmp[apply(e_tmp, 1, function(x) length(which(x == 0)) < min(length(col1), length(col2))), ]
    expr_grouped <- e_tmp[rownames(e_tmp) %in% ref_tmp$Sequence, ]

    # sample_pairs <- factor(
    #   x = rep(c(samples[j], reference_sample), each = c(length(col1), length(col2))),
    #   levels = c(reference_sample, samples[j]),
    #   ordered = T
    # )
    sample_pairs <- factor(
      # 修改：将 each 改为 times，以支持不同长度的重复
      x = rep(c(samples[j], reference_sample), times = c(length(col1), length(col2))),
      levels = c(reference_sample, samples[j]),
      ordered = T
    )
    
    if (nrow(expr_grouped) > 0) {
      result_tmp <- dep_analysis(expr = expr_grouped, group = sample_pairs)
      result_tmp <- result_tmp[result_tmp$adj.P.Val < 0.05, ]
      result_final <- rbind(result_final, result_tmp)
    }
    
  }

  # Calculating: RC -----------------------------------------------
  if (is.null(result_final) || nrow(result_final) == 0) {
    return(list(DEPs = NULL, logfc = NULL, COR = NA, cor_plot = NULL))
  }
  
  result_final <- as.data.table(result_final)
  # 提取测试结果的关键列
  # 注意：dep_analysis 返回结果中 Sequence 是第一列/行名
  result_trim <- data.frame(
    Sequence = result_final$Sequence,
    Sample.Pair = result_final$Sample.Pair,
    logFC.Test = result_final$logFC
  )
  
  # 构建合并键
  result_trim$name <- paste(result_trim$Sequence, result_trim$Sample.Pair)
  ref_dt$name <- paste(ref_dt$Sequence, ref_dt$Sample.Pair)
  # 合并 (Inner Join)
  result_withref <- merge(result_trim, ref_dt, by = "name")

  # 整理最终表格用于绘图和计算
  df_test <- data.frame(
    Name = result_withref$name,
    Sequence = result_withref$Sequence.x,
    Sample.Pair = result_withref$Sample.Pair.x,
    logFC.Test = result_withref$logFC.Test,
    logFC.Reference = result_withref$log2FC # 来自参考集的 value
  )

  cor_value <- cor(x = df_test$logFC.Test, y = df_test$logFC.Reference)
  cor_value <- round(cor_value, 3)

  # # Plot ----------------------------------------------------------
  # if (plot) {
  #   text_custom_theme <- element_text(size = 16,
  #                                     face = "plain",
  #                                     color = "black",
  #                                     hjust = 0.5)
  #
  #   scale_axis_r <- c(min(df_test$logFC.Reference),
  #                     max(df_test$logFC.Reference))
  #   scale_axis_t <- c(min(df_test$logFC.Test),
  #                     max(df_test$logFC.Test))
  #   limit <- max(abs(c(scale_axis_r, scale_axis_t)))
  #   limit_axis <- c(- limit, limit)
  #
  #   plot_title <- paste("RC = ", cor_value, sep = "")
  #   plot_subtitle <- paste("(Number of peptides = ", nrow(df_test), ")", sep = "")
  #
  #   p <- ggplot(df_test, aes(x = .data$logFC.Reference, y = .data$logFC.Test)) +
  #     theme_few() +
  #     theme(plot.title = text_custom_theme,
  #           plot.subtitle = text_custom_theme,
  #           axis.title = text_custom_theme,
  #           axis.text = text_custom_theme,
  #           legend.title = text_custom_theme,
  #           legend.text = element_text(size = 16, color = "gray40")) +
  #     labs(y = "log2FC (Test Datasets)",
  #         x = "log2FC (Reference Datasets)",
  #         title = plot_title,
  #         subtitle = plot_subtitle) +
  #     coord_fixed(xlim = limit_axis, ylim = limit_axis)
  #
  #   if (show_sample_pairs == T) {
  #     # Create color mapping based on numerator (first sample in pair)
  #     # Extract numerator from Sample.Pair
  #     df_test$Numerator <- sapply(strsplit(as.character(df_test$Sample.Pair), "/"), function(x) x[1])
  #     unique_numerators <- unique(df_test$Numerator)
  #     pair_colors <- get_sample_colors(unique_numerators)
  #
  #     # Create full pair color mapping
  #     colors_custom <- c()
  #     for (pair in unique(df_test$Sample.Pair)) {
  #       numerator <- strsplit(as.character(pair), "/")[[1]][1]
  #       colors_custom[pair] <- pair_colors[numerator]
  #     }
  #
  #     p <- p +
  #       geom_point(aes(color = .data$Sample.Pair), size = 2.5, alpha = .5) +
  #       scale_color_manual(values = colors_custom)
  #   } else {
  #     p <- p + geom_point(color = "steelblue4", size = 2.5, alpha = .1)
  #   }
  #
  # }

  # Plot ----------------------------------------------------------
  p <- NULL
  if (plot) {
    text_custom_theme <- element_text(
      size = 16,
      face = "plain",
      color = "black",
      hjust = 0.5
    )

    scale_axis_r <- c(
      min(df_test$logFC.Reference),
      max(df_test$logFC.Reference)
    )
    scale_axis_t <- c(
      min(df_test$logFC.Test),
      max(df_test$logFC.Test)
    )
    limit <- max(abs(c(scale_axis_r, scale_axis_t)))
    limit_axis <- c(-limit, limit)

    plot_title <- paste("RC = ", cor_value, sep = "")
    plot_subtitle <- paste("(Number of peptides = ", nrow(df_test), ")", sep = "")

    # 动态生成样本对的颜色映射（基于分子的颜色）
    unique_comps <- unique(df_test$Sample.Pair)
    pair_colors <- sapply(unique_comps, function(comp_name) {
      # 提取样本对的分子（numerator），格式为 "Sample1/Sample2"
      numerator_sample <- strsplit(as.character(comp_name), "/")[[1]][1]
      color_palette <- c(
        "D5" = "#4CC3D9",
        "D6" = "#7BC8A4",
        "F7" = "#FFC65D",
        "M8" = "#F16745"
      )
      if (numerator_sample %in% names(color_palette)) {
        return(color_palette[[numerator_sample]])
      } else {
        return("gray") # 默认颜色
      }
    })
    names(pair_colors) <- unique_comps

    p <- ggplot(df_test, aes(x = .data$logFC.Reference, y = .data$logFC.Test)) +
      theme_few() +
      theme(
        plot.title = text_custom_theme,
        plot.subtitle = text_custom_theme,
        axis.title = text_custom_theme,
        axis.text = text_custom_theme,
        legend.title = text_custom_theme,
        legend.text = element_text(size = 16, color = "gray40")
      ) +
      labs(
        y = "log2FC (Test Datasets)",
        x = "log2FC (Reference Datasets)",
        title = plot_title,
        subtitle = plot_subtitle
      ) +
      coord_fixed(xlim = limit_axis, ylim = limit_axis)

    if (show_sample_pairs == T) {
      p <- p +
        geom_point(aes(color = .data$Sample.Pair), size = 2.5, alpha = .5) +
        scale_color_manual(values = pair_colors, name = "Sample Pair")
    } else {
      p <- p + geom_point(color = "steelblue4", size = 2.5, alpha = .1)
    }
  }

  # # Save & Output -------------------------------------------------
  # if (!is.null(output_dir)) {
  #   if (plot) {
  #     output_dir_final1 <- file.path(output_dir, "corr_plot.png")
  #     ggsave(output_dir_final1, p, height = 5.5, width = 5.5)
  #   }
  #   output_dir_final2 <- file.path(output_dir, "deps_table.tsv")
  #   output_dir_final3 <- file.path(output_dir, "corr_table.tsv")
  #   write.table(result_final, output_dir_final2, sep = "\t", row.names = F)
  #   write.table(df_test, output_dir_final3, sep = "\t", row.names = F)
  # }

  output_list <- list(
    DEPs = result_final,
    logfc = df_test,
    COR = cor_value,
    cor_plot = p
  )

  return(output_list)
}

#' Calculating Recall of nominal characteristics
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file (to map library to sample type)
#' @export
qc_recall <- function(expr_dt, meta_dt) {
  
  # 1. 加载内置的定性数据集 (Qualitative / Recall)
  # 使用 utils::data 加载，确保不受用户环境影响
  utils::data("reference_dataset_quali", package = "protqc", envir = environment())
  ref_dt <- reference_dataset_quali
  
  # 2. 数据清洗
  # A. 检查列名 (根据你的 CSV: sample, peptide_sequence)
  if (!all(c("sample", "peptide_sequence") %in% colnames(ref_dt))) {
    stop("Built-in reference_dataset_quali is missing required columns: sample, peptide_sequence")
  }
  
  # B. 样本名标准化: 去除 "Quartet " 前缀 (例如 "Quartet D5" -> "D5")
  # 这样才能和 metadata 中的 sample (D5, D6...) 匹配
  ref_dt$sample_clean <- gsub("Quartet ", "", ref_dt$sample)
  
  # 3. 确定共有样本类型
  sample_types <- unique(meta_dt$sample)
  ref_sample_types <- unique(ref_dt$sample_clean)
  common_types <- intersect(sample_types, ref_sample_types)
  
  if (length(common_types) == 0) {
    warning("No common sample types found between metadata and reference dataset.")
    return(NA)
  }
  
  group_recalls <- c()
  
  # 4. 遍历每一类样本 (D5, D6, F7, M8) 进行计算
  for (st in common_types) {
    # 分母: 该类样本的参考肽段集合
    ref_seqs <- unique(ref_dt$peptide_sequence[ref_dt$sample_clean == st])
    n_ref <- length(ref_seqs)
    
    if (n_ref == 0) next
    
    # 获取该类样本对应的 Library ID
    libs <- meta_dt$library[meta_dt$sample == st]
    valid_libs <- libs[libs %in% colnames(expr_dt)]
    
    if (length(valid_libs) == 0) next
    
    # 对每一个样本(Library)单独计算 Recall
    lib_recalls <- c()
    for (lib in valid_libs) {
      # 提取检测到的肽段 (假设非 NA 且 > 0 为检测到)
      vals <- expr_dt[[lib]]
      detected_idx <- which(!is.na(vals) & vals > 0)
      
      # expr_dt 第一列通常是 Peptide Sequence
      det_seqs <- expr_dt[[1]][detected_idx]
      
      # 分子: 交集数量
      n_det_ref <- length(intersect(det_seqs, ref_seqs))
      
      rec <- n_det_ref / n_ref
      lib_recalls <- c(lib_recalls, rec)
    }
    
    # 取该类样本的平均 Recall
    if (length(lib_recalls) > 0) {
      group_recalls <- c(group_recalls, mean(lib_recalls))
    }
  }
  
  if (length(group_recalls) == 0) return(NA)
  
  # 5. 最终结果: 所有类别的平均 Recall
  final_recall <- mean(group_recalls)
  return(final_recall)
}