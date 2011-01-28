plotGenMap <- function (map, dense = FALSE, nMarker = TRUE, ...)
{
    if (class(map) == "gpData") map <- map$map
    chr <- unique(map$chr)
    chr <- chr[!is.na(chr)]

    # add legend to the left side
    if (dense)  layout(matrix(2:1, ncol = 2), width = c(0.82, 0.18))

    if (dense) {
    # colors from RColorBrewer red - green
    cols <- c("#A50026", "#D73027", "#F46D43", "#FDAE61",
                "#FEE08B", "#FFFFBF", "#D9EF8B", "#A6D96A", "#66BD63",
                "#1A9850", "#006837")[11:1]
        par(mar = c(5, 0, 4, 3.8) + 0.1)
        image(seq(-0.4, 0.4, length = 20), seq(from = 0, to = 1,
            length = 11), matrix(rep(seq(from = 0, to = 1, length = 11),
            20), nrow = 20, byrow = TRUE), col = cols, axes = FALSE,
            xlab = "")
        axis(side = 4, at = round(seq(from = 0, to = 1, length = 11),
            4), las = 1)
        par(mar = c(5, 4, 4, 1) + 0.1)
    }

   # make an empty plot 
    plot(map, type = "n", xaxt = "n", xlim = c(0.5, length(chr) +
        0.5), ylim = c(min(map$pos, na.rm = TRUE), max(map$pos,
        na.rm = TRUE) * 1.1), ...)
        
   # x-axis     
    axis(side = 1, at = seq(along = chr), labels = chr)


   # plot each chromosome
    for (i in seq(along = chr)) {
        n <- sum(map$chr == chr[i], na.rm = TRUE)
        start <- min(map$pos[map$chr == chr[i]], na.rm = TRUE)
        end <- max(map$pos[map$chr == chr[i]], na.rm = TRUE)
        if (dense) {
            bw <- (end - start)/n
            densEst <- density(map$pos[map$chr == chr[i]], kernel = "rectangular",
                from = min(map$pos[map$chr == chr[i]], na.rm = TRUE),
                to = max(map$pos[map$chr == chr[i]], na.rm = TRUE),
                cut = TRUE, bw = bw, na.rm=TRUE)
            densEst$y <- densEst$y/max(densEst$y)

            image(seq(i - 0.4, i + 0.4, length = 20), densEst$x,
                matrix(rep(densEst$y, 20), nrow = 20, byrow = TRUE),
                col = cols, add = TRUE)
        }
        else {
            lines(x = c(i, i), y = c(start, end))
            for (j in 1:n) {
                lines(x = c(i - 0.4, i + 0.4), y = rep(map$pos[map$chr ==
                  chr[i]][j], 2))
            }
        }
        # add nr. of markers
        if (nMarker)
            text(i, max(map$pos) * 1.05, sum(map$chr == chr[i],na.rm=TRUE))
    }

}

