#author: "Markus Viljanen"
#date: "2024-03-19"
library(ggplot2)
library(dplyr)
library(tidyr)

fn <- "results_simulation.csv"
df <- read.table(fn, header=T, sep=",", row.names=NULL)
cov.true <- 1.0
df$true <- cov.true

# Results from first simulation
df %>% filter((Simulation == 1) & Confounded & Control.z2) %>% 
  select(Spatial.Misalignment, Measurement.Error, 
         Model, mean, sd, rmse.y, rmse.z1)

# Summary of simulations: coefficient and predictive accuracy
df %>% filter(Confounded & Control.z2) %>% #%>% filter(Simulation <= 100) 
  group_by(Spatial.Misalignment, Measurement.Error, Model) %>% 
  summarize(n=n(), 
            mean.coef=round(mean(mean,na.rm=T),2), 
            mean.coef.sd=round(mean(sd,na.rm=T),2),#sd.coef=round(sd(mean),2), 
            #coverage=round(mean((true >= X0.025quant) & (true <= X0.975quant)),2),
            rmse.coef=round(mean(sqrt((true-mean)**2), na.rm=T),2),  
            rmse.y=round(mean(rmse.y),3), 
            rmse.z1=round(mean(rmse.z1),3)) %>%
  write.csv(file="results_paper/results_simulation.csv", row.names=F)

# Summary of simulations: estimated coefficient of X1 in two models
temp <- df %>% filter(Confounded & Control.z2 & Spatial.Misalignment & Measurement.Error) %>%
  select(Simulation, Model, mean) %>% pivot_wider(names_from=Model, values_from=mean)
temp %>% ggplot(aes(x=`Two-Stage`, y=Joint)) + geom_point() + geom_abline() +
  geom_smooth(method='lm', formula= y~x)

# Plot estimated coefficients from joint vs. two-stage model
df.subset <- df %>% filter(Confounded & Control.z2 & Spatial.Misalignment & Measurement.Error) %>% 
  mutate(Model = factor(Model, levels=c("Two-Stage", "Joint")))
plot <- ggplot(df.subset, aes(x=Simulation, y=mean, col=Model, group=Model)) + 
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(position=position_dodge(width=0.5), aes(ymin=X0.025quant, ymax=X0.975quant)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(name='Estimated Coefficient') + geom_hline(yintercept=cov.true, linetype='dashed') #+ 
  #facet_wrap( ~ Spatial, ncol = 2, labeller = labeller(Spatial = c("no" = "No Spatial SDM",
  #                                                                 "yes" = "Spatial SDM")))
plot
ggsave('results_paper/simulation_coef.png', width=3000, heigh=1200, units='px')


