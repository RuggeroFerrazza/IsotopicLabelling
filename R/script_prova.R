## 25.02.2016
## Script prova per modificare IsotopicLabeling according to PietroFranceschi


# Load the package
library(IsotopicLabeling)

# Load data
data("xcms_obj")

# Obtain table of peaks
peak_table <- table_xcms(xcms_obj)

# Fitting ....
fitted_abundances <- main_labeling(peak_table, compound="X40H77NO8P", labeling="C",
                                   mass_shift=0.05, RT=285, RT_shift=20, chrom_width=7, initial_abundance=NA)

# Summary
summary(fitted_abundances)

# Plot of patterns
plot(fitted_abundances, type="patterns", saveplots=T)

# Plot of residuals
plot(fitted_abundances, type="residuals", saveplots=F)

# Plot of summary
plot(fitted_abundances, type="summary", saveplots=F)

# Save the results
save_labeling(fitted_abundances)
