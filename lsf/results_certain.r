args <- commandArgs(trailingOnly = TRUE)
species.take <- args[[1]]
print(species.take)

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


stk_y <- inla.stack(data = list(y=sf.plots[[species.take]]),
              A = list(A_y, 1, 1, 1, 1, 1, 1, 1, 1), 
              effects = list(c(list(u.field=1:spde$n.spde), list(intercept=1)), 
                              list(terrain=sf.plots$terrain), 
                              list(Dagnummer.z=sf.plots$Dagnummer.z), 
                              list(pH.z=sf.plots$pH.z), 
                              list(ORG_STOF.z=sf.plots$ORG_STOF.z.mean),
                              list(CN.z=sf.plots$CN.z.mean), 
                              list(N_totaal.z=sf.plots$N_totaal.z.mean), 
                              list(P_totaal.z=sf.plots$P_totaal.z.mean), 
                              list(K_CaCl2.z=sf.plots$K_CaCl2.z.mean)
                              ),
              tag="est.y")

n <- nrow(sf.grid)
stk_pred <- inla.stack(data=list(y=rep(NA, n)),  
                        effects=list(c(list(u.field=1:spde$n.spde), list(intercept=1)),
                                      list(terrain=sf.grid$terrain), 
                                      list(Dagnummer.z=sf.grid$Dagnummer.z), 
                                      list(pH.z=sf.grid$pH.z), 
                                      list(ORG_STOF.z=sf.grid$ORG_STOF.z.mean),
                                      list(CN.z=sf.grid$CN.z.mean), 
                                      list(N_totaal.z=sf.grid$N_totaal.z.mean), 
                                      list(P_totaal.z=sf.grid$P_totaal.z.mean), 
                                      list(K_CaCl2.z=sf.grid$K_CaCl2.z.mean)
                              ),
                        A=list(A_pred, 1, 1, 1, 1, 1, 1, 1, 1),
                        tag="est.pred")

stk <- inla.stack(stk_y, stk_pred)

formula <- y ~ -1 + intercept + f(inla.group(terrain), model = "rw2") + 
  f(inla.group(Dagnummer.z), model = "rw2") +  f(inla.group(pH.z), model = "rw2") +
  CN.z + ORG_STOF.z + N_totaal.z + P_totaal.z + K_CaCl2.z + f(u.field, model=spde)

output <-  tryCatch(
      {inla(formula, num.threads=1,
                data=inla.stack.data(stk, spde=spde), family="binomial",
                control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), #
                control.inla = list(strategy = "adaptive", int.strategy = "eb")) #int.strategy = "eb"
      }, error = function(cond) 
          {message(cond)
          return(NULL)}
  )

# Model results
if (!is.null(output)) {
  u_field <- inla.spde2.result(output, name="u.field", spde)
  index.pred <- inla.stack.index(stk, tag="est.pred")$data
  
  # SDM regression terms and predictions
  terms.sdm <- rbind(
    # Splines for day number, pH, terrain
    bind_rows(output$summary.random, .id="term") %>% 
        filter(term %in% c("inla.group(Dagnummer.z)", "inla.group(pH.z)", "inla.group(terrain)")) %>% 
        mutate(term=str_replace_all(term, c("inla.group\\("="", "\\)"=""))) %>% 
        select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Linear terms for CN, ORG_STOF, N_totaal, P_totaal, K_CaCl2
    output$summary.fixed %>% 
      mutate(ID=row.names(output$summary.fixed)) %>% remove_rownames %>% mutate(term="fixed") %>% 
      select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Spatial component
    bind_rows(list(log.variance=u_field$summary.log.variance.nominal, 
                    log.range=u_field$summary.log.range.nominal), .id="ID") %>% 
      remove_rownames %>% mutate(term="spatial") %>% 
      select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Predictions
    cbind(sf.grid %>% as.data.frame %>% mutate(term="Predict", ID=Plot), 
                      output$summary.fitted.values[index.pred,]) %>% 
    select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Predictions (log-odds scale)
    cbind(sf.grid %>% as.data.frame %>% mutate(term="APredict", ID=Plot), 
                      output$summary.linear.predictor[index.pred,]) %>% 
    select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  ) %>% 
    mutate(species=species.take, meta="Certain", model = "SDM") %>%
    select(species, meta, model, term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  
  # Save
  terms <- terms.sdm#rbind(, terms.abiotics %>% mutate(species=species.take))
  fn <- sprintf("/mnt/scratch_dir/viljanem/results_certain/%s.csv", species.take)
  write.table(terms, fn, sep = ",", col.names = F, row.names=F, append=F)
}

