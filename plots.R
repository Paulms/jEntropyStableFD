library(rio)
setwd("~/code/jEntropyStableDF")
data = import('test_2_400.txt')
names(data) = c("x", "MS", "ESC-0.2", "ESNC-0.2")
library(ggplot2)
library(ggthemes)
ggplot(data=data, aes(x=x, y=MS)) +
  geom_line()+geom_point()
library(reshape2)
melted <- melt(data, id="x")

ref = import('test_2_reference.txt')
names(ref) = c("x","REF") 
library(rbokeh)
figure(data= melted) %>%
  ly_lines(x, value, color = variable) %>%
  ly_lines(x, REF, data = ref, legend = "REF")

library(plotly)
plot_ly(melted, x = ~x, y = ~value, color = ~variable, mode='lines')



ggplot() +
  geom_line(data=melted, aes(x=x, y=value, color=variable)) +
  geom_line(data=ref, aes(x=x, y=REF)) +
  scale_color_ptol("cyl") +
  theme_minimal()
