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

i <- 1
stacks <- list()
for (var in variables) {
  # Observations of the variable
  sf.var <- sf.plots %>% filter(!is.na(get(var)))
  A_var <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sf.var %>% as.data.frame %>% select(X,Y)))
  # One output for every estimated variable (+ species occurrence)
  y <- matrix(nrow=nrow(sf.var), ncol=length(variables) + 1)
  y.0 <- matrix(nrow=nrow(sf.grid), ncol=length(variables) + 1)
  y[,i] <- sf.var[[var]]
  # Stack
  effects <- list()
  effects[[sprintf("%s.i.field", var)]] = 1:spde$n.spde
  effects[[sprintf("%s.intercept", var)]] = 1 #asdf
  stacks[[var]] <- inla.stack(data=list(y=y), effects=list(effects), A=list(A_var), tag=sprintf("est.%s", var))
  stacks[[sprintf("%s.pred", var)]] <- inla.stack(data=list(y=y.0), effects=list(effects), A=list(A_pred), tag=sprintf("est.pred.%s", var))
  i <- i + 1
}


# Species occurrence
y <- matrix(nrow=nrow(sf.plots), ncol=length(variables) + 1)
y[,i] <- sf.plots[[species.take]]
# Stack
stk_y <- inla.stack(data=list(y=y),  
                        effects=list(c(list(ORG_STOF.z.field=1:spde$n.spde, 
                                            CN.z.field=1:spde$n.spde,
                                            N_totaal.z.field=1:spde$n.spde, 
                                            P_totaal.z.field=1:spde$n.spde, 
                                            K_CaCl2.z.field=1:spde$n.spde,
                                            u.field=1:spde$n.spde), list(intercept=1)),
                                      list(terrain=sf.plots$terrain), 
                                      list(Dagnummer.z=sf.plots$Dagnummer.z), 
                                      list(pH.z=sf.plots$pH.z)),
                        A=list(A_y, 1, 1, 1),
                        tag="est.y")
# Predicted over grid
y <- matrix(nrow=nrow(sf.grid), ncol=length(variables) + 1)
stk_pred <- inla.stack(data=list(y=y),  
                        effects=list(c(list(ORG_STOF.z.field=1:spde$n.spde, 
                                            CN.z.field=1:spde$n.spde,
                                            N_totaal.z.field=1:spde$n.spde, 
                                            P_totaal.z.field=1:spde$n.spde, 
                                            K_CaCl2.z.field=1:spde$n.spde,
                                            u.field=1:spde$n.spde), list(intercept=1)),
                                      list(terrain=sf.grid$terrain), 
                                      list(Dagnummer.z=sf.grid$Dagnummer.z), 
                                      list(pH.z=sf.grid$pH.z)),
                        A=list(A_pred, 1, 1, 1),
                        tag="est.pred")

stk <- inla.stack(stacks[["ORG_STOF.z"]], stacks[["ORG_STOF.z.pred"]],
                  stacks[["CN.z"]], stacks[["CN.z.pred"]], 
                  stacks[["N_totaal.z"]], stacks[["N_totaal.z.pred"]],
                  stacks[["P_totaal.z"]], stacks[["P_totaal.z.pred"]],
                  stacks[["K_CaCl2.z"]], stacks[["K_CaCl2.z.pred"]],
                  stk_y, stk_pred)

formula <- y ~ -1 + intercept + ORG_STOF.z.intercept + CN.z.intercept + N_totaal.z.intercept + P_totaal.z.intercept + K_CaCl2.z.intercept +
f(inla.group(terrain), model = "rw2") + f(inla.group(Dagnummer.z), model = "rw2") +  f(inla.group(pH.z), model = "rw2") +
  f(ORG_STOF.z.i.field, model=spde) + f(ORG_STOF.z.field, copy="ORG_STOF.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) +
  f(CN.z.i.field, model=spde) + f(CN.z.field, copy="CN.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) +
  f(N_totaal.z.i.field, model=spde) + f(N_totaal.z.field, copy="N_totaal.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001))))  +
  f(P_totaal.z.i.field, model=spde) + f(P_totaal.z.field, copy="P_totaal.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001))))  +
  f(K_CaCl2.z.i.field, model=spde) + f(K_CaCl2.z.field, copy="K_CaCl2.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
  f(u.field, model=spde)

families <- c(rep("gaussian", length(variables)), "binomial")

output <-  tryCatch(
      {inla(formula, family=families, num.threads=1,
                data=inla.stack.data(stk), 
                control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=i),
                control.inla = list(strategy = "adaptive", int.strategy = "eb"))#int.strategy = "eb"
      }, error = function(cond) 
          {message(cond)
          return(NULL)}
  )

# Model results
if (!is.null(output)) {
  # random field and index in stack
  u_field <- inla.spde2.result(output, name="u.field", spde)
  index.pred <- inla.stack.index(stk, tag="est.pred")$data
  
  # SDM regression terms and predictions
  terms.sdm <- rbind(
    # Splines for day number, pH, terrain
    bind_rows(output$summary.random, .id="term") %>% 
        filter(term %in% c("inla.group(Dagnummer.z)", "inla.group(pH.z)", "inla.group(terrain)")) %>% 
        mutate(term=str_replace_all(term, c("inla.group\\("="", "\\)"=""))) %>% 
        select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Intercept
    output$summary.fixed %>% 
      mutate(ID=row.names(output$summary.fixed)) %>% remove_rownames %>% mutate(term="fixed") %>% 
      filter(ID  == "intercept") %>%
      select(term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`),
    # Linear terms for CN, ORG_STOF, N_totaal, P_totaal, K_CaCl2
    output$summary.hyperpar %>% 
      mutate(ID=row.names(output$summary.hyperpar)) %>% remove_rownames %>% mutate(term="fixed") %>% 
      filter(ID %>% startsWith("Beta")) %>%
      mutate(ID=str_replace_all(ID, c("Beta for "="", ".field"=""))) %>%
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
    mutate(species=species.take, meta="Uncertain", model = "SDM") %>%
    select(species, meta, model, term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  
  j <- 1
  terms.abiotics <- list()
  for (var in variables) {
    ## transform posterior marginal of precision into variance and calculate statistics
    i_field <- inla.spde2.result(output, name=sprintf("%s.i.field", var), spde)
    var.marginal <- data.frame(inla.zmarginal(inla.tmarginal(function(x) log(1/x), output$marginals.hyperpar[[j]]), silent=T)) %>%
      rename(mean=mean, sd=sd, `0.025quant`=quant0.025, `0.25quant`=quant0.25, `0.5quant`=quant0.5, `0.75quant`=quant0.75, `0.975quant`=quant0.975)
    # predictions
    index.pred <- inla.stack.index(stk, tag=sprintf("est.pred.%s", var))$data
    i.pred <- output$summary.linear.predictor[index.pred,] %>% select(mean, sd, `0.025quant`, `0.5quant`,`0.975quant`)
    ##  TODO: the standard deviation is not correct
    #i.mean <- (A_pred %*% i_field$summary.values[,2])[,1]
    #i.sd <- (A_pred %*% i_field$summary.values[,3])[,1]
    #i.pred <- data.frame(mean=i.mean, sd=i.sd, `0.025quant`=i.mean+qnorm(0.025)*i.sd, `0.5quant`=i.mean, `0.975quant`=i.mean+qnorm(0.975)*i.sd)%>%
    #  rename(mean=mean, sd=sd, `0.025quant`=X0.025quant, `0.5quant`=X0.5quant, `0.975quant`=X0.975quant)
    
    # Abiotic factor regression terms and predictions
    terms.abiotics[[var]] <- rbind(
      # Intercept
      output$summary.fixed %>% 
        mutate(ID=row.names(output$summary.fixed)) %>% remove_rownames %>% mutate(term="fixed") %>% 
        filter(ID  == sprintf("%s.intercept", var)) %>% mutate(ID = "intercept") %>% 
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
    j <- j + 1
  }
  terms.abiotics <- bind_rows(terms.abiotics, .id="model") %>% 
    mutate(species=species.take, meta="Uncertain") %>%
    select(species, meta, model, term, ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`)
  
  # Save
  terms <- rbind(terms.sdm, terms.abiotics)
  fn <- sprintf("/mnt/scratch_dir/viljanem/results_uncertain/%s.csv", species.take)
  write.table(terms, fn, sep = ",", col.names = F, row.names=F, append=F)
}

