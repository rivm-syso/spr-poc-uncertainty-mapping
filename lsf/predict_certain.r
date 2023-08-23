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


sf.train = sf.plots %>% filter(naam != provincie) #%>% filter(bodemgebruik %in% c("Bos", "Droog natuurlijk terrein", "Nat natuurlijk terrein"))
sf.test = sf.plots %>% filter(naam == provincie) #%>% filter(bodemgebruik %in% c("Bos", "Droog natuurlijk terrein", "Nat natuurlijk terrein"))
A_train <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sf.train %>% as.data.frame %>% select(X,Y)))
A_test <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sf.test %>% as.data.frame %>% select(X,Y)))

### Certain abiotics ###

# Fit
stk_y <- inla.stack(data = list(y=sf.train[[species.take]]),
                A = list(A_train, 1, 1, 1, 1, 1, 1, 1, 1), 
                effects = list(c(list(u.field=1:spde$n.spde), list(intercept=1)), 
                                list(terrain=sf.train$terrain), 
                                list(Dagnummer.z=sf.train$Dagnummer.z), 
                                list(pH.z=sf.train$pH.z), 
                                list(ORG_STOF.z=sf.train$ORG_STOF.z.mean),
                                list(CN.z=sf.train$CN.z.mean), 
                                list(N_totaal.z=sf.train$N_totaal.z.mean), 
                                list(P_totaal.z=sf.train$P_totaal.z.mean), 
                                list(K_CaCl2.z=sf.train$K_CaCl2.z.mean)
                                ),
                tag="est.y")

stk_pred <- inla.stack(data=list(y=rep(NA, nrow(sf.test))),  
                        effects=list(c(list(u.field=1:spde$n.spde), list(intercept=1)),
                                        list(terrain=sf.test$terrain), 
                                        list(Dagnummer.z=sf.test$Dagnummer.z), 
                                        list(pH.z=sf.test$pH.z), 
                                        list(ORG_STOF.z=sf.test$ORG_STOF.z.mean),
                                        list(CN.z=sf.test$CN.z.mean), 
                                        list(N_totaal.z=sf.test$N_totaal.z.mean), 
                                        list(P_totaal.z=sf.test$P_totaal.z.mean), 
                                        list(K_CaCl2.z=sf.test$K_CaCl2.z.mean)
                                ),
                        A=list(A_test, 1, 1, 1, 1, 1, 1, 1, 1),
                        tag="est.pred")

stk <- inla.stack(stk_y, stk_pred)

formula <- y ~ -1 + intercept + f(inla.group(terrain), model = "rw2") + 
    f(inla.group(Dagnummer.z), model = "rw2") +  f(inla.group(pH.z), model = "rw2") +
    CN.z + ORG_STOF.z + N_totaal.z + P_totaal.z + K_CaCl2.z + f(u.field, model=spde)

index.pred <- inla.stack.index(stk, tag="est.pred")$data

output <-  tryCatch(
    {inla(formula, num.threads=1,
                data=inla.stack.data(stk, spde=spde), family="binomial",
                control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), 
                control.inla = list(strategy = "adaptive", int.strategy = "eb"), control.compute=list(config = TRUE)) #int.strategy = "eb"
    }, error = function(cond) 
        {message(cond)
        return(NULL)}
)

# First model results
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
    
    predictions <- cbind(predicted, observed) %>% mutate(meta = "Certain", model="SDM", species=species.take, naam = provincie)
    prevalences <- data.frame(area_prevalences, meta = "Certain", model="SDM", species=species.take, naam = provincie)

    fn <- sprintf("/mnt/scratch_dir/viljanem/predictions_certain/%s_%s.csv", species.take, provincie)
    write.table(predictions, fn, sep = ",", col.names = F, row.names=F, append=F)
    fn <- sprintf("/mnt/scratch_dir/viljanem/prevalences_certain/%s_%s.csv", species.take, provincie)
    write.table(prevalences, fn, sep = ",", col.names = F, row.names=F, append=F)
}
