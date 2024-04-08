#author: "Markus Viljanen"
#date: "2024-03-19"

library(INLA)
library(sf)
library(ggplot2)
library(colorspace)
library(cowplot)
library(dplyr)
library(tidyr)

fn <- "results_simulation.csv"
if (file.exists(fn)) file.remove(fn)

# Simulation parameters
for (i in seq(1,200)) {
  print(sprintf("   Simulation %d", i))
  #i <- 1
  
  # Simulation parameters
  sd.obs <- sqrt(0.3)
  sd.spde <- sqrt(0.5)
  N1 <- 100
  N2 <- 1000
  
  # True effects
  b.0 <- -2.0
  b.Z1 <- 1.0
  b.Z2 <- 1.0
  
  # Create mesh
  limits <- cbind(c(0,1,1,0,0), c(0,0,1,1,0))
  mesh <- inla.mesh.2d(loc.domain=limits, cutoff=0.03, max.edge=c(0.03,0.12), offset=-0.03) 
  
  # Create grid
  sf.boundary <- st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
  sf.rectangles <- st_make_grid(sf.boundary, cellsize = .02, square = T)
  sf.grid <- st_sf(geometry=sf.rectangles)
  sf.grid[c("X", "Y")] <- sf.grid %>% st_geometry %>% st_centroid() %>% st_coordinates()
  
  # Base R plot of grid and mesh
  plot(sf.boundary, col='blue')
  plot(sf.rectangles, add=T)
  plot(mesh, add=T)
  
  # Define random field (Z)
  range0 <- 0.3
  sigma0 <- sqrt(sd.spde)
  kappa0 <- sqrt(8)/range0
  tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
  
  spde <- inla.spde2.matern(mesh,
                            B.tau=matrix(c(log(tau0), -1, +1),nrow=1,ncol=3),
                            B.kappa=matrix(c(log(kappa0), 0, -1),nrow=1,ncol=3), 
                            theta.prior.mean=c(0, 0),
                            theta.prior.prec=c(0.1, 0.1))
  mesh.index <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  
  # Realization of the random field
  Q <- inla.spde2.precision(spde=spde, theta=c(0,0))
  z <- as.vector(inla.qsample(n=1, Q=Q))#, seed=1434
  z1 <- as.vector(inla.qsample(n=1, Q=Q))#, seed=1434
  z2 <- as.vector(inla.qsample(n=1, Q=Q))#, seed=1434
  z0 <- as.vector(inla.qsample(n=1, Q=Q))#, seed=1434
  w <- as.vector(inla.qsample(n=1, Q=Q))
  
  # True values of covariate Z at grid points
  A.grid <- inla.spde.make.A(mesh, loc=as.matrix(subset(as.data.frame(sf.grid), select=c(X,Y))))
  sf.grid$Z.true <- as.vector(A.grid %*%  z)
  sf.grid$W.true <- as.vector(A.grid %*%  w)
  sf.grid$Z1.true <- sf.grid$Z.true + as.vector(A.grid %*%  z1)
  sf.grid$Z2.true <- sf.grid$Z.true + as.vector(A.grid %*%  z2)
  sf.grid$Z2.0cor <- as.vector(A.grid %*%  z0) + as.vector(A.grid %*%  z2)
  c(var(sf.grid$Z.true), 
    var(sf.grid$Z1.true), 
    var(sf.grid$Z2.true), 
    cov(sf.grid$Z1.true, sf.grid$Z2.true))
  
  # Create sampling locations
  N <- N1 + N2
  xy.obs <- data.frame(X=sample(1:N / N), Y=sample(1:N / N))
  
  # Run repeated samplings of the given simulation
  for (spatial.misalignment in c(F,T)) {
    for (measurement.error in c(F,T)) {
      for (confounded in c(T)) {#c(F,T)
        options.z2 <- c(T)#if (confounded) c(F,T) else c(F)
        for (control.z2 in options.z2) {
          for (spatial in c(T)) { #c(F,T)
            print(sprintf("Spatial Misalignment %s, Measurement error %s - Counfounded %s, Control %s, Spatial %s",
                          spatial.misalignment, measurement.error, confounded, control.z2, spatial))
            start <- Sys.time()
            
            # All locations where species or abiotics were sampled
            sf.obs <- st_sf(geometry=st_as_sf(xy.obs, coords = c("X","Y")))
            sf.obs[c("X", "Y")] <- sf.obs %>% st_geometry %>% st_centroid() %>% st_coordinates()
            # True values of covariate Z at at all locations
            A.obs <- inla.spde.make.A(mesh, loc=as.matrix(subset(as.data.frame(sf.obs), select=c(X,Y))))
            sf.obs$Z.true <- as.vector(A.obs %*%  z)
            sf.obs$W.true <- as.vector(A.obs %*%  w)
            sf.obs$Z1.true <- sf.obs$Z.true + as.vector(A.obs %*%  z1)
            sf.obs$Z2.true <- sf.obs$Z.true + as.vector(A.obs %*%  z2)
            sf.obs$Z2.0cor <- as.vector(A.obs %*%  z0) + as.vector(A.obs %*%  z2)
            # Observed values of covariate Z at all locations
            sf.obs$Z1.obs <- sf.obs$Z1.true 
            sf.obs$Z2.obs <- sf.obs$Z2.true 
            if (measurement.error) {
              sf.obs$Z1.obs <- sf.obs$Z1.obs + rnorm(N,mean=0,sd=sd.obs)
              sf.obs$Z2.obs <- sf.obs$Z2.obs + rnorm(N,mean=0,sd=sd.obs)
            }
            # Censor values of covariate Z at species locations
            sf.obs$abiotic <- sample(c(rep(T, N1), rep(F, N2)))
            if (spatial.misalignment) {
              sf.obs$Z1.obs <- ifelse(sf.obs$abiotic, sf.obs$Z1.obs, NA)
              sf.obs$Z2.obs <- ifelse(sf.obs$abiotic, sf.obs$Z2.obs, NA)
            }
            
            # True species suitability at grid points
            if (confounded) {
              theta <- b.0 + b.Z1*sf.grid$Z1.true + b.Z2*sf.grid$Z2.true + sf.grid$W.true
            } else {
              theta <- b.0 + b.Z1*sf.grid$Z1.true + b.Z2*sf.grid$Z2.0cor + sf.grid$W.true
            }
            sf.grid$y <- 1/(1+exp(-theta))
            
            # Observed values of species occurrence y at sampling locations
            if (confounded) {
              theta <- b.0 + b.Z1*sf.obs$Z1.true + b.Z2*sf.obs$Z2.true + sf.obs$W.true
            } else {
              theta <- b.0 + b.Z1*sf.obs$Z1.true + b.Z2*sf.obs$Z2.0cor + sf.obs$W.true
            }
            y <- rbinom(N, 1, prob=1/(1+exp(-theta)))
            sf.obs$y <- ifelse(!sf.obs$abiotic, y, NA) 
            
            # Estimate effect (direct covariate value)
            if (!spatial.misalignment) {
              print("    Model 0")
              stack.est <- inla.stack(data=list(y=sf.obs$y),
                                      A=list(A.obs, 1),
                                      effects=list(list(field=1:spde$n.spde, intercept=1), 
                                                   list(z1=sf.obs$Z1.obs, z2=sf.obs$Z2.obs)),
                                      tag="est")
              stk <- inla.stack(stack.est)
              if (spatial) {
                if(control.z2) {
                  formula <- y ~ -1 + intercept + z1 + z2 + f(field,model=spde) 
                } else {
                  formula <- y ~ -1 + intercept + z1 + f(field,model=spde) 
                }
              } else {
                if(control.z2) {
                  formula <- y ~ -1 + intercept + z1 + z2
                } else {
                  formula <- y ~ -1 + intercept + z1
                }
              }
              output <- inla(formula, data=inla.stack.data(stk,spde=spde),
                             family="binomial", control.predictor=list(A=inla.stack.A(stk),compute=TRUE,link=1))
              results.0 <- data.frame(subset(output$summary.fixed[2,], select=c(mean, sd, `0.025quant`, `0.5quant`, `0.975quant`, mode)),
                                      rmse.y=NA, sd.obs=sd.obs, sd.z1=NA, rmse.z1=NA, 
                                      N1=N1, N2=N2, Simulation=i, 
                                      Spatial.Misalignment=spatial.misalignment, Measurement.Error=measurement.error,
                                      Spatial=spatial, Confounded=confounded, Control.z2=control.z2, Model="Direct")
              write.table(results.0, fn, sep = ",", col.names = !file.exists(fn), row.names=F, append=T)
              
            }
            
            # Kriging of abiotic variable
            print("    Kriging")
            sd <- list()
            for (var in c("Z1", "Z2")) {
              print(sprintf("    %s", var))
              var.obs <- sprintf("%s.obs", var)
              var.pred <- sprintf("%s.pred", var)
              stack.est <- inla.stack(data=list(y = sf.obs[[var.obs]]),
                                      A=list(A.obs), effects=list(c(mesh.index, list(intercept=1))), 
                                      tag="est")
              stack.pred <- inla.stack(data=list(y=NA), 
                                       A=list(A.grid), effects=list(c(mesh.index, list(intercept=1))),
                                       tag="pred") 
              join.stack <- inla.stack(stack.est, stack.pred)
              
              formula <- y ~ -1 + intercept + f(field, model=spde)
              
              output <- inla(formula, family="gaussian", 
                             data=inla.stack.data(join.stack, spde=spde),
                             control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE, link=1))
             # spde.result <- inla.spde2.result(inla=output,name="field",spde=spde)
              
              # predictions at field visits
              index.est <- inla.stack.index(join.stack, tag="est")$data
              sf.obs[[var.pred]] <- output$summary.fitted.values[index.est,]$mean
              
              # predictions at grid
              index.pred <- inla.stack.index(join.stack, tag="pred")$data
              sf.grid[[var.pred]] <- output$summary.fitted.values[index.pred,]$mean
              
              sd[[var.obs]] <- mean(output$summary.fitted.values[index.pred,]$sd)
            }
            
            # Estimate effect (predicted covariate)
            print("    Model 1")
            stack.est <- inla.stack(data=list(y=sf.obs$y),
                                    A=list(A.obs, 1),
                                    effects=list(list(field=1:spde$n.spde, intercept=1), 
                                                 list(z1=sf.obs$Z1.pred, z2=sf.obs$Z2.pred)),
                                    tag="est")
            stack.pred <- inla.stack(data=list(y=NA),
                                     A=list(A.grid, 1),
                                     effects=list(list(field=1:spde$n.spde, intercept=1), 
                                                  list(z1=sf.grid$Z1.pred, z2=sf.grid$Z2.pred)),
                                     tag="pred")
            stk <- inla.stack(stack.est, stack.pred)
            if (spatial) {
              if(control.z2) {
                formula <- y ~ -1 + intercept + z1 + z2 + f(field,model=spde) 
              } else {
                formula <- y ~ -1 + intercept + z1 + f(field,model=spde) 
              }
            } else {
              if(control.z2) {
                formula <- y ~ -1 + intercept + z1 + z2
              } else {
                formula <- y ~ -1 + intercept + z1
              }
            }
            output <- inla(formula, data=inla.stack.data(stk,spde=spde),
                           family="binomial", control.predictor=list(A=inla.stack.A(stk),compute=TRUE,link=1))
            y.pred <- output$summary.fitted.values[inla.stack.index(stk, tag="pred")$data,]$mean
            rmse.z1 <- sqrt(mean((sf.grid$Z1.true-sf.grid$Z1.pred)^2))
            rmse.y <- sqrt(mean((sf.grid$y-y.pred)^2))
            results.1 <- data.frame(subset(output$summary.fixed[2,], select=c(mean, sd, `0.025quant`, `0.5quant`, `0.975quant`, mode)),
                                    rmse.y=rmse.y, sd.obs=sd.obs, sd.z1=sd$Z1.obs, rmse.z1=rmse.z1, 
                                    N1=N1, N2=N2, Simulation=i, 
                                    Spatial.Misalignment=spatial.misalignment, Measurement.Error=measurement.error,
                                    Spatial=spatial, Confounded=confounded, Control.z2=control.z2, Model="Two-Stage")
            write.table(results.1, fn, sep = ",", col.names = !file.exists(fn), row.names=F, append=T)
            
            # Estimate effect (uncertain covariate)
            print("    Model 2")
            ncol <- if(control.z2) 3 else 2
            # z1 covariate
            y <- matrix(nrow=nrow(sf.obs), ncol=ncol)
            y[,1] <- sf.obs$Z1.obs
            stack.est.z1 <- inla.stack(data=list(y = y),
                                       A=list(A.obs), effects=list(list(z1.i.field=1:spde$n.spde, z1.intercept=1)), 
                                       tag="est.z1")
            y.0 <- matrix(nrow=nrow(sf.grid), ncol=ncol)
            stack.pred.z1 <- inla.stack(data=list(y = y.0), 
                                        A=list(A.grid), effects=list(list(z1.i.field=1:spde$n.spde, z1.intercept=1)),
                                        tag="pred.z1") 
            # species
            y <- matrix(nrow=nrow(sf.obs), ncol=ncol)
            y[,ncol] <- sf.obs$y
            stack.est <- inla.stack(data=list(y=y),
                                    A=list(A.obs),
                                    effects=list(list(field=1:spde$n.spde, z1.field=1:spde$n.spde, z2.field=1:spde$n.spde, intercept=1)),
                                    tag="est")
            y <- matrix(nrow=nrow(sf.grid), ncol=ncol)
            stack.pred <- inla.stack(data=list(y=y),
                                     A=list(A.grid),
                                     effects=list(list(field=1:spde$n.spde, z1.field=1:spde$n.spde, z2.field=1:spde$n.spde, intercept=1)),
                                     tag="pred")
            if (control.z2) {
              # z2 covariate
              y <- matrix(nrow=nrow(sf.obs), ncol=ncol)
              y[,2] <- sf.obs$Z2.obs
              stack.est.z2 <- inla.stack(data=list(y = y),
                                         A=list(A.obs), effects=list(list(z2.i.field=1:spde$n.spde, z2.intercept=1)), 
                                         tag="est.z2")
              y.0 <- matrix(nrow=nrow(sf.grid), ncol=ncol)
              stack.pred.z2 <- inla.stack(data=list(y = y.0), 
                                          A=list(A.grid), effects=list(list(z2.i.field=1:spde$n.spde, z2.intercept=1)),
                                          tag="pred.z2") 
              stk <- inla.stack(stack.est.z1, stack.pred.z1, stack.est.z2, stack.pred.z2, stack.est, stack.pred)
              families <- c("gaussian", "gaussian", "binomial")
            } else {
              stk <- inla.stack(stack.est.z1, stack.pred.z1, stack.est, stack.pred)
              families <- c("gaussian", "binomial")
            }
            
            if (spatial) { 
              if(control.z2) {
                formula <- y ~ -1 + intercept + 
                  z1.intercept + f(z1.i.field, model=spde) + f(z1.field, copy="z1.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
                  z2.intercept + f(z2.i.field, model=spde) + f(z2.field, copy="z2.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
                  f(field,model=spde)
              } else {
                formula <- y ~ -1 + intercept + 
                  z1.intercept + f(z1.i.field, model=spde) + f(z1.field, copy="z1.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
                  f(field,model=spde)
              }
            } else {
              if(control.z2) {
                formula <- y ~ -1 + intercept + 
                  z1.intercept + f(z1.i.field, model=spde) + f(z1.field, copy="z1.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) + 
                  z2.intercept + f(z2.i.field, model=spde) + f(z2.field, copy="z2.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001)))) 
              } else {
                formula <- y ~ -1 + intercept + 
                  z1.intercept + f(z1.i.field, model=spde) + f(z1.field, copy="z1.i.field", fixed=FALSE, hyper=list(theta=list(param=c(0, 0.001))))
              }
            }
            
            output <- inla(formula, family=families, data=inla.stack.data(stk),
                           control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=ncol))
            y.pred <- output$summary.fitted.values[inla.stack.index(stk, tag="pred")$data,]$mean
            z1.pred <- output$summary.linear.predictor[inla.stack.index(stk, tag="pred.z1")$data,]$mean
            z1.sd <- mean(output$summary.linear.predictor[inla.stack.index(stk, tag="pred.z1")$data,]$sd)
            rmse.y <- sqrt(mean((sf.grid$y-y.pred)^2))
            rmse.z1 <- sqrt(mean((sf.grid$Z1.true-z1.pred)^2))
            if (control.z2) {
              idx <- ifelse(spatial, 9, 7)
            } else {
              idx <- ifelse(spatial, 6, 4)
            }
            results.2 <- data.frame(subset(output$summary.hyperpar[idx,], select=c(mean, sd, `0.025quant`, `0.5quant`, `0.975quant`, mode)),
                                    rmse.y=rmse.y, sd.obs=sd.obs, sd.z=z1.sd, rmse.z1=rmse.z1, 
                                    N1=N1, N2=N2, Simulation=i, 
                                    Spatial.Misalignment=spatial.misalignment, Measurement.Error=measurement.error,
                                    Spatial=spatial, Confounded=confounded, Control.z2=control.z2, Model="Joint")
            write.table(results.2, fn, sep = ",", col.names = !file.exists(fn), row.names=F, append=T)
            
            end <- Sys.time()
            print(end-start)
            
          }
        }
      }
    }
  }
}

# Plot illustration of simulation

#Plot latent random field and observations of Z1
plot1a <- ggplot() +
  geom_sf(aes(geometry=geometry, fill=Z1.true), data=sf.grid) +
  geom_sf(aes(geometry=geometry, fill=Z1.obs), color='black', pch=21,size=2,
          data=sf.obs %>% filter(!is.na(Z1.obs)))+
  scale_fill_gradient2(low='blue', midpoint=0, high='red', name=expression(X[1])) +
  coord_sf(xlim=c(0,1), ylim=c(0,1)) + # boundary of 0.05?
  theme(legend.position="bottom", legend.key.width = unit(0.1, "npc"))
#Plot latent random field and observations of Z1
plot1b <- ggplot() +
  geom_sf(aes(geometry=geometry, fill=Z2.true), data=sf.grid) +
  geom_sf(aes(geometry=geometry, fill=Z2.obs), color='black', pch=21,size=2,
          data=sf.obs %>% filter(!is.na(Z2.obs)))+
  scale_fill_gradient2(low='blue', midpoint=0, high='red', name=expression(X[2])) +
  coord_sf(xlim=c(0,1), ylim=c(0,1)) + # boundary of 0.05?
  theme(legend.position="bottom", legend.key.width = unit(0.1, "npc"))
# Plot species observations
sf.grid$log.odds <- log(sf.grid$y/(1-sf.grid$y))
plot2 <- ggplot() +
  geom_sf(aes(geometry=geometry, fill=log.odds), data=sf.grid) +
  geom_sf(aes(geometry=geometry, color=as.factor(y)), pch=19, size=2,
          data=sf.obs )+#%>% filter(!is.na(y))
  scale_fill_continuous_divergingx(palette = 'PiYG', mid = median(sf.grid$log.odds), name=expression(eta))+
  scale_colour_manual(values=c("darkmagenta", "darkgreen"),name=expression(Y)) +
  coord_sf(xlim=c(0,1), ylim=c(0,1)) + # boundary of 0.05?
  theme(legend.position="bottom", legend.key.width = unit(0.04, "npc"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
# Plot side-by-side
plot.all <- ggdraw() +
  draw_plot(plot1a, x = 0, y = 0, width = 0.333, height = 0.95) +
  draw_plot(plot1b, x = 0.333, y = 0, width = 0.333, height = 0.95) +
  draw_plot(plot2, x = 0.666, y = 0, width = 0.333, height = 0.95) +
  draw_plot_label(label = c("A", "B", "C"), size = 15, x = c(0.01, 0.34, 0.68), y = c(0.99,0.99,0.99)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
plot.all
ggsave('results_paper/simulation.png', width=3000, heigh=1200, units='px')


# Plot illustration of Kriging: Berkson & Classical type errors

# True values of covariate Z at grid points
Y.choice <- 0.80 #0.75
sf.line <- data.frame(X=seq(0,1,0.001), Y=Y.choice)
A.line <- inla.spde.make.A(mesh, loc=as.matrix(sf.line))
# True Z1 and observation standard error
sf.line$Z1.true <-  as.vector(A.line %*%  z) + as.vector(A.line %*%  z1)
sf.line$Z1.eps.sd <- sd.obs
# Predict covariate at a given Y latitude line
stack.est <- inla.stack(data=list(y = sf.obs$Z1.obs),
                        A=list(A.obs), effects=list(c(mesh.index, list(intercept=1))), 
                        tag="est")
stack.pred <- inla.stack(data=list(y=NA), 
                         A=list(A.line), effects=list(c(mesh.index, list(intercept=1))),
                         tag="pred") 
join.stack <- inla.stack(stack.est, stack.pred)
formula <- y ~ -1 + intercept + f(field, model=spde)
output <- inla(formula, family="gaussian", 
               data=inla.stack.data(join.stack, spde=spde),
               control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE, link=1))
index.pred <- inla.stack.index(join.stack, tag="pred")$data
# Predicted Z1 and prediction standard error
sf.line$Z1.pred <- output$summary.fitted.values[index.pred,]$mean
sf.line$Z1.pred.sd <- output$summary.fitted.values[index.pred,]$sd
# Data frame with true vs. predicted lines
sf.slice <- sf.line %>% pivot_longer(c(Z1.true, Z1.pred),
                                     names_to="Covariate",
                                     values_to='Value') %>%
  mutate(SD = ifelse(Covariate == "Z1.true", Z1.eps.sd, Z1.pred.sd)) %>% 
  mutate(Covariate=recode(Covariate, Z1.true="X1.true", Z1.pred="X1.pred"))
# Data frame with observations
sf.obs.slice <- sf.obs %>% filter((Y > Y.choice-0.05) & (Y < Y.choice+0.05)) %>% 
  select(X, Z1.obs) %>% filter(!is.na(Z1.obs)) %>%
  pivot_longer(Z1.obs, names_to="Covariate", values_to='Value') %>% 
  mutate(Covariate=recode(Covariate, Z1.obs="X1.obs"))
# Plot
plot.line <- ggplot(sf.slice) + geom_line(aes(x=X,y=Value, col=Covariate)) + 
  geom_ribbon(aes(x=X, ymax=Value+1.96*SD, ymin=Value-1.96*SD, fill=Covariate), alpha=0.3) +
  geom_point(aes(x=X, y=Value, col=Covariate), data = sf.obs.slice) +
  geom_hline(yintercept=0, linetype = 'dotted') + theme(legend.position = "top") +
  scale_color_manual(name='Covariate',
                     breaks=c('X1.obs', 'X1.pred', 'X1.true'),
                     labels=c(expression(paste(X[1],'*',(r))),expression(paste("E(",X[1](r),"|", D[1],")")),expression(X[1](r))),
                     values=c('X1.obs'='blue', 
                              'X1.pred'='orange',
                              'X1.true'='blue')) +
  scale_fill_manual(name='Covariate',
                    breaks=c('X1.obs', 'X1.pred', 'X1.true'),
                    labels=c(expression(paste(X[1],'*',(r))),expression(paste("E(",X[1](r),"|", D[1],")")),expression(X[1](r))),
                    values=c('X1.obs'='blue', 
                             'X1.pred'='orange',
                             'X1.true'='blue'), guide="none") + 
  guides(color = guide_legend(override.aes = list(linetype = c(0, 1, 1), 
                                                  shape=c(16,NA,NA)))) + 
  xlab(expression(r))
plot.line
ggsave('results_paper/simulation_kriging.png', width=1500, heigh=900, units='px')


# Load saved results
fn <- "results_simulation.csv"
df <- read.table(fn, header=T, sep=",", row.names=NULL)
df$true <- b.Z1

# Results from first simulation
df %>% filter((Simulation == 1) & Confounded & Control.z2) %>% 
  select(Spatial.Misalignment, Measurement.Error, 
         Model, mean, sd, rmse.y, rmse.z1)

# Summary of simulations: coefficient and predictive accuracy
df %>% filter(Simulation <= 100) %>% filter(Confounded & Control.z2) %>%
  group_by(Spatial.Misalignment, Measurement.Error, Model) %>% 
  summarize(n=n(), 
            mean.coef=round(mean(mean,na.rm=T),2), 
            mean.coef.sd=round(mean(sd,na.rm=T),2),#sd.coef=round(sd(mean),2), 
            rmse.coef=round(mean(sqrt((true-mean)**2), na.rm=T),2),  
            rmse.y=round(mean(rmse.y),3), 
            rmse.z1=round(mean(rmse.z1),3)) %>%
  write.csv(file="results_paper/results_simulation.csv", row.names=F)

# Plot estimated coefficients from joint vs. two-stage model
df.subset <- df %>% filter(Confounded & Control.z2 & Spatial.Misalignment & Measurement.Error) %>% 
  mutate(Model = factor(Model, levels=c("Two-Stage", "Joint")))
plot <- ggplot(df.subset, aes(x=Simulation, y=mean, col=Model, group=Model)) + 
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(position=position_dodge(width=0.5), aes(ymin=X0.025quant, ymax=X0.975quant)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(name='Estimated Coefficient') + geom_hline(yintercept=b.Z1, linetype='dashed') #+ 
plot
ggsave('results_paper/simulation_coef.png', width=3000, heigh=1200, units='px')



