library(dplyr)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = FALSE)

## DATA LOADING ##
# Load data files from defined path
path <- '~/advRassignment_EDann/'
geneX.table <- read.delim(paste0(path,'data.txt'), row.names = 1)
slow.growth.profile <- read.delim(paste0(path,'slowgrowthprofile.txt'), row.names = 1, col.names = c('gene','slow.growth.profile')) 
doubling.time <- read.delim(paste0(path, 'reldoublingtime.txt'), col.names = c('rel.growth.rate', 'genename'))

## DATA RESHAPING ##
# Match information about growth rate and resemblance to slow growth profile for each deletion strain
strain.growth <- doubling.time[match(colnames(geneX.table), doubling.time$genename),] %>% 
  mutate(slow.growth.cov = sapply(genename,   # Add column for covariance of strain expression and slow growth profile
                                  function(strain) 
                                    cov(geneX.table[,strain],y = slow.growth.profile)))

## LINEAR MODEL ## 
strain.growth.LM <- function(strain.growth){
  # Compute linear model and relative statistics
  with(data = strain.growth, {
    correlation <- cor(rel.growth.rate, slow.growth.cov)
    model <- lm(rel.growth.rate ~ slow.growth.cov)
    coeffs <- sapply(coef(model), round, digits=2) # Round coefficients to 2 significant digits
    s <- summary(model)
    pval <- round(s$coefficients[2,'Pr(>|t|)'], digits=4)
    return(list(pcc = round(correlation, digits = 2),
                lm = model, 
                regression.coefficients = coeffs, 
                regression.pvalue = pval))
  })
}

print.lm.stats <- function(pcc, coeffs, pval){
  # Print to screen formatted statistics for linear regression
  print( noquote(paste('Pearson Correlation Coefficient =', round(pcc, digits = 2))) )
  print( noquote(paste('Regression line equation:', 'y =', coeffs[1], '+', coeffs[2], 'x')) )
  print( noquote(paste('P-value:', pval)) )
  }

yeast.LM <- strain.growth.LM(strain.growth) 
print.lm.stats(yeast.LM$pcc, yeast.LM$regression.coefficients, yeast.LM$regression.pvalue )

## PLOT ##
# Plot scatterplot
plot <- ggplot(strain.growth, aes(x=rel.growth.rate, y=slow.growth.cov)) +
  geom_point() # Scatter plot

# Add regression line and confidence intervals
plot <- plot + geom_smooth(fill='grey', 
              method = "lm",
              formula = y ~ x,
              color = "red") 

# Add labels and text
plot <- plot +
  ylab('covariance expression profile - slow growth profile') +
  xlab('relative doubling time') +
  geom_text_repel(aes(label = genename)) + # Add labels
  geom_label(x=0.08,y=0.6, label=paste0('y = ', yeast.LM$regression.coefficients[1], ' + ', yeast.LM$regression.coefficients[2], 'x',
                                       ' (p = ', round(yeast.LM$regression.pvalue, digits = 2), ')\n',
                                       'PCC = ', yeast.LM$pcc
                                      ),
             vjust=1, hjust=0) + # define vertical and horizontal positioning
  ggtitle("Effect of slow growth expression profile") 

# Show and save
plot + ggsave(paste0(path, "yeastGrowthVSslowGrowthProfile.pdf"))


