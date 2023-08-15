args <- commandArgs(trailingOnly = TRUE)

library(INLA)
library(dplyr)
library(stringr)
library(tibble)

sf.plots <- readRDS("sf.plots.RDS")
sf.grid <- readRDS("sf.grid.RDS")

xy.plots <- sf.plots %>% as.data.frame %>% select(X,Y)
xy.grid <- sf.grid %>% as.data.frame %>% select(X,Y)

step <- 10000 / 1000
domain <- inla.nonconvex.hull(as.matrix(xy.plots), convex=-0.04)
mesh <- inla.mesh.2d(boundary=domain, max.edge=c(step*2, step*4), cutoff=step)

spde <- inla.spde2.matern(mesh=mesh,alpha=2)##print(spde$n.spde)
s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 

A_y <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.plots))
A_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.grid))
variables <- c("ORG_STOF.z", "CN.z", "N_totaal.z", "P_totaal.z", "K_CaCl2.z")



terms.abiotics <- list()
i <- 1
for (variable in variables) {
  print(variable)

  stack.est <- inla.stack(data=list(y = sf.plots[[variable]]),
                          A=list(A_y), effects=list(c(s.index, list(intercept=1))), 
                          tag="est")
  stack.pred <- inla.stack(data=list(y=NA), 
                           A=list(A_pred), effects=list(c(s.index, list(intercept=1))),
                           tag="pred") 
  join.stack <- inla.stack(stack.est, stack.pred)
  
  formula <- y ~ -1 + intercept + f(spatial.field, model=spde)
  
  output <- inla(formula, family="gaussian", num.threads=1,
                      data=inla.stack.data(join.stack, spde=spde),
                      control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE, link=1),
                control.inla = list(strategy = "adaptive", int.strategy = "eb"))
  
  # predictions at field visits
  index.est <- inla.stack.index(join.stack, tag="est")$data
  variable.est <- output$summary.fitted.values[index.est,]
  colnames(variable.est) <- sprintf("%s.%s", variable, colnames(variable.est))
  sf.plots[colnames(variable.est)] <- variable.est

  # predictions at grid
  index.pred <- inla.stack.index(join.stack, tag="pred")$data
  variable.pred <- output$summary.fitted.values[index.pred,]
  colnames(variable.pred) <- sprintf("%s.%s", variable, colnames(variable.pred))
  sf.grid[colnames(variable.pred)] <- variable.pred

  # transform posterior marginal of precision into variance and calculate statistics
  i_field <- inla.spde2.result(output, name="spatial.field", spde)
  var.marginal <- data.frame(inla.zmarginal(inla.tmarginal(function(x) log(1/x), output$marginals.hyperpar[[1]]), silent=T)) %>%
    rename(mean=mean, sd=sd, `0.025quant`=quant0.025, `0.25quant`=quant0.25, `0.5quant`=quant0.5, `0.75quant`=quant0.75, `0.975quant`=quant0.975)
  # predictions
  cols <- sprintf(c("%s.mean", "%s.sd", "%s.0.025quant", "%s.0.5quant", "%s.0.975quant"), variable)
  i.pred <- sf.grid %>% as.data.frame %>% select(all_of(cols))
  colnames(i.pred) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  
  # Abiotic factor regression terms and predictions
  terms.abiotics[[variable]] <- rbind(
  	# Intercept
  	output$summary.fixed %>% 
  		mutate(ID=row.names(output$summary.fixed)) %>% remove_rownames %>% mutate(term="fixed") %>% 
  		select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Precision
  	var.marginal %>% 
  		mutate(ID="log.variance.eps") %>% remove_rownames %>% mutate(term="spatial") %>% 
  		select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
  	# Spatial component
  	bind_rows(list(log.variance=i_field$summary.log.variance.nominal, 
  	               log.range=i_field$summary.log.range.nominal), .id="ID") %>% 
  	  remove_rownames %>% mutate(term="spatial") %>% 
  		select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Predictions
    cbind(sf.grid %>% as.data.frame %>% mutate(term="Predict", ID=Plot), i.pred) %>% 
    select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  	)
  i <- i + 1
}

terms.abiotics <- bind_rows(terms.abiotics, .id="model") %>% 
  mutate(species="", meta="Certain") %>%
  select(species, meta, model, term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)

terms <- terms.abiotics
fn <- "/mnt/scratch_dir/viljanem/results_certain/abiotics.csv"
write.table(terms, fn, sep = ",", col.names = F, row.names=F, append=F)
