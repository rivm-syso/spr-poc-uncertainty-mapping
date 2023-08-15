args <- commandArgs(trailingOnly = TRUE)
species.take <- args[[1]]
provincie <- args[[2]]
print(species.take)
print(provincie)

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

sf.train = sf.plots %>% filter(naam != provincie) #%>% filter(bodemgebruik %in% c("Bos", "Droog natuurlijk terrein", "Nat natuurlijk terrein"))
sf.test = sf.plots %>% filter(naam == provincie) #%>% filter(bodemgebruik %in% c("Bos", "Droog natuurlijk terrein", "Nat natuurlijk terrein"))
A_train <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sf.train %>% as.data.frame %>% select(X,Y)))
A_test <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sf.test %>% as.data.frame %>% select(X,Y)))

### Uncertain abiotics ###

# Species occurrence
y <- matrix(nrow=nrow(sf.train), ncol=length(variables) + 1)
y[,i] <- sf.train[[species.take]]
# Stack
stk_y <- inla.stack(data=list(y=y),  
                        effects=list(c(list(ORG_STOF.z.field=1:spde$n.spde, 
                                            CN.z.field=1:spde$n.spde,
                                            N_totaal.z.field=1:spde$n.spde, 
                                            P_totaal.z.field=1:spde$n.spde, 
                                            K_CaCl2.z.field=1:spde$n.spde,
                                            u.field=1:spde$n.spde), list(intercept=1)),
                                        list(terrain=sf.train$terrain), 
                                        list(Dagnummer.z=sf.train$Dagnummer.z), 
                                        list(pH.z=sf.train$pH.z)),
                        A=list(A_train, 1, 1, 1),
                        tag="est.y")
# Predicted over grid
y <- matrix(nrow=nrow(sf.test), ncol=length(variables) + 1)
stk_pred <- inla.stack(data=list(y=y),  
                        effects=list(c(list(ORG_STOF.z.field=1:spde$n.spde, 
                                            CN.z.field=1:spde$n.spde,
                                            N_totaal.z.field=1:spde$n.spde, 
                                            P_totaal.z.field=1:spde$n.spde, 
                                            K_CaCl2.z.field=1:spde$n.spde,
                                            u.field=1:spde$n.spde), list(intercept=1)),
                                        list(terrain=sf.test$terrain), 
                                        list(Dagnummer.z=sf.test$Dagnummer.z), 
                                        list(pH.z=sf.test$pH.z)),
                        A=list(A_test, 1, 1, 1),
                        tag="est.pred")

stk <- inla.stack(stacks[["ORG_STOF.z"]], stacks[["CN.z"]], stacks[["N_totaal.z"]],
                    stacks[["P_totaal.z"]], stacks[["K_CaCl2.z"]], stk_y, stk_pred)

formula <- y ~ -1 + intercept + ORG_STOF.z.intercept + CN.z.intercept + N_totaal.z.intercept + P_totaal.z.intercept + K_CaCl2.z.intercept +
f(inla.group(terrain), model = "rw2") + f(inla.group(Dagnummer.z), model = "rw2") +  f(inla.group(pH.z), model = "rw2") +
  f(ORG_STOF.z.i.field, model=spde) + f(ORG_STOF.z.field, copy="ORG_STOF.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) +
  f(CN.z.i.field, model=spde) + f(CN.z.field, copy="CN.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) +
  f(N_totaal.z.i.field, model=spde) + f(N_totaal.z.field, copy="N_totaal.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001))))  +
  f(P_totaal.z.i.field, model=spde) + f(P_totaal.z.field, copy="P_totaal.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001))))  +
  f(K_CaCl2.z.i.field, model=spde) + f(K_CaCl2.z.field, copy="K_CaCl2.z.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
  f(u.field, model=spde)

index.pred <- inla.stack.index(stk, tag="est.pred")$data

families <- c(rep("gaussian", length(variables)), "binomial")
output <-  tryCatch(
    {inla(formula, family=families, num.threads=1, #verbose=T,
            data=inla.stack.data(stk), 
            control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=i),
            control.inla = list(strategy = "adaptive", int.strategy = "eb"), control.compute=list(config = TRUE)) #int.strategy = "eb"
    }, error = function(cond) 
        {message(cond)
        return(NULL)}
)

# Second model results
if (!is.null(output)) {
    # get inla posterior sample
    s1 <- inla.posterior.sample(400, output) 
    linear_predictor <- inla.posterior.sample.eval("APredictor", s1)[index.pred,]
    P = 1/(1+exp(-linear_predictor))
    sample_mean <- function(p) mean(rbinom(length(p), 1, p))
    area_prevalences <- apply(P, 2, sample_mean) #colMeans(1/(1+exp(-linear_predictor)))

    # Predicted vs observed
    predicted <- output$summary.fitted.values[index.pred,c("mean", "sd", "0.025quant", "0.975quant")]
    observed <- sf.test[[species.take]]
    
    predictions <- cbind(predicted, observed) %>% mutate(meta = "Uncertain", model="SDM", species=species.take, naam = provincie)
    prevalences <- data.frame(area_prevalences, meta = "Uncertain", model="SDM", species=species.take, naam = provincie)

    fn.1 <- sprintf("/mnt/scratch_dir/viljanem/predictions_uncertain/%s_%s.csv", species.take, provincie)
    fn.2 <- sprintf("/mnt/scratch_dir/viljanem/prevalences_uncertain/%s_%s.csv", species.take, provincie)
    write.table(predictions, fn.1, sep = ",", col.names = F, row.names=F, append=F)
    write.table(prevalences, fn.2, sep = ",", col.names = F, row.names=F, append=F)
}
