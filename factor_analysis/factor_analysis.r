# Factor Analysis on Pima Indians Diabetes Dataset


helpers <- c("corrplot", "devtools", "heplots", "klaR", "psych")
lapply(helpers, library, character.only = TRUE) # Import all helper libraries

library(mlbench)
data(PimaIndiansDiabetes)


generate_factor_analysis <- function(
    p_partition = "pos",
    p_drop_zero = TRUE
) {
    #' Performs simple factor analysis on the Pima Indians Diabetes dataset.
    #' 
    #' @description Blah
    #'
    #' @param p_partition A string to characterize how to partition the samples
    #'  based on the `diabetes` column (either "pos" or "neg"). If neither is
    #'  specified, it uses all data (default is "pos").
    #' @param p_drop_zero A boolean indicating whether to drop rows with zero
    #'  values (default is TRUE).
    #' @returns A list containing the partitioned data, PCA results, covariance
    #'  matrix, covariance eigenvalues, covariance eigenvectors, and the factor
    #'  loading matrix.
    #' @examples
    #'  \dontrun{
    #'  generate_factor_analysis(p_partition = "pos", p_drop_zero = TRUE)
    #'  generate_factor_analysis(p_partition = "neg", p_drop_zero = FALSE)
    #'  }
    
    # Partitioning data
    if (p_partition == "pos") {
        part <- PimaIndiansDiabetes[
            PimaIndiansDiabetes$diabetes == "pos",
            !(names(PimaIndiansDiabetes) %in% "diabetes")
        ]
    } else if (p_partition == "neg") {
        part <- PimaIndiansDiabetes[
            PimaIndiansDiabetes$diabetes == "neg",
            !(names(PimaIndiansDiabetes) %in% "diabetes")
        ]
    } else {
        print("Argument `p_partition` is invalid; using all data")
        part <- PimaIndiansDiabetes[
            !(names(PimaIndiansDiabetes) %in% "diabetes")
        ]
        p_partition <- "all data"
    }

    # Remove rows with zero values, if specified
    if (p_drop_zero) {
        part <- part[apply(part != 0, 1, all), ]
    }

    # Compute the spectral decomposition
    ev <- eigen((cov(part)))$values
    eve <- eigen(cov(part))$vectors
    plot(
        ev,
        type = "b",
        main = sprintf("Eigenvalues of Covariance Matrix (%s)", p_partition),
        xlab = "Index",
        ylab = "Eigenvalue"
    )

    # Compute the covariance matrix
    s_cov_matrix <- cor(part)
    s_eigenvals <- eigen(s_cov_matrix)$values
    s_eigenvecs <- eigen(s_cov_matrix)$vectors

    # Scree plots
    par(mfrow = c(1, 2))
    plot(s_eigenvals, type = "b", main = "Scree Plot")
    screeplot(
        princomp(s_cov_matrix),
        main = "Scree Plot"
    )

    # Cumulative percentage plot
    cumsum(s_eigenvals) / length(s_eigenvals)
    par(mfrow = c(1, 1))
    plot(
        cumsum(s_eigenvals) / length(s_eigenvals),
        type = "b",
        main = "Cumulative Percentage Plot",
        xlab = "Index",
        ylab = "Cumulative Explained Variance"
    )

    # Perform Principal Component Analysis (PCA)
    data_pca <- princomp(part, cor = TRUE)

    # Factor Loadings
    load_matrix <- as.matrix(loadings(data_pca))
    png(
        filename = sprintf("./factors_%s.png", p_partition),
        width = 800,
        height = 600
    )
    corrplot(
        load_matrix,
        is.corr = FALSE,
        method = "circle",
        tl.col = "black",
        tl.cex = 2,
        number.cex = 2
    )
    dev.off()

    # Return values
    result <- list(
        data_partition = part,
        data_pca = data_pca,
        covariance_matrix = s_cov_matrix,
        covariance_eigenvalues = s_eigenvals,
        covariance_eigenvectors = s_eigenvecs,
        loadings = load_matrix
    )

    return (result)
}


create_factor_ranking <- function(
    data_matrix,
    eigenvalues,
    eigenvectors,
    num_components = 4
) {
    #' Ranking of the samples based on latent factor analysis.
    #' 
    #' @description Creates a ranking based on the eigenvalues and eigenvectors
    #' of the covariance matrix using the principal component method of factor
    #' analysis (see Section 13.3.1 of Rencher, Methods of Multivariate
    #' Analysis).
    #' 
    #' @param data_matrix A numeric matrix of the data.
    #' @param eigenvalues A numeric vector of eigenvalues of the covariance
    #'  matrix.
    #' @param eigenvectors A numeric matrix of eigenvectors of the covariance
    #'  matrix.
    #' @param num_components An integer specifying the number of components to
    #'  extract (default is 4).
    #' @returns A numeric matrix of the transformed data based on the selected
    #'  principal components. Can be interpreted as a ranking of the samples
    #'  based on the latent factors.


    # Extract the first `num_components` values
    pcs <- (eigenvectors[, 1:num_components])
    sqrt_eigs <- sqrt(eigenvalues[1:num_components])

    return (as.matrix(data_matrix) %*% pcs %*% as.matrix(sqrt_eigs))
}
