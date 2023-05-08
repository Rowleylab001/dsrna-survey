
ntlenKHR.df = read.csv("~/Downloads/ntlenKHR.txt")
ntlenKHS.df = read.csv("~/Downloads/ntlenKHS.txt")

ggplot(data=ntlenKHR.df, aes(x=ntlen)) +
  geom_histogram(bins = 30) +
  labs(title="KHR lengths for all BLAST hits")

ggplot(data=ntlenKHS.df, aes(x=ntlen)) +
  geom_histogram(bins = 30) +
  labs(title="KHS lengths for all BLAST hits")