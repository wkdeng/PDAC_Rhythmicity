#' CHIRAL with Dual Penalty System: Sampling Window + Reference Phase Constraints
#'
#' Benchmark-specific copy of scripts/CHIRAL_dual_penalty.R.
#' The legacy script is intentionally left unchanged. This copy fixes the
#' reference peak/trough scoring used during benchmarking so candidate sample
#' phases are scored against the observed reference-gene expression pattern.
#'
#' Enhanced version of CHIRAL that combines sampling window constraints with reference phase penalties
#' to improve absolute phase inference accuracy and consistency.
#'
#' @param E Required matrix of gene expression. Samples should be on the columns, genes on the rows.
#' @param iterations Number of maximum iterations, default 500.
#' @param clockgenes Set of clock genes to be used for the inference. Default NULL.
#' @param tau2 Tau parameter for the prior on the gene coefficient, default NULL (auto-calculated).
#' @param u u parameter for the prior on the gene means, default NULL (auto-calculated).
#' @param sigma2 Standard deviation of the data points, default NULL (auto-calculated).
#' @param TSM Switches the two state model on/off, default TRUE.
#' @param mean.centre.E Center around empirical mean of the data, default TRUE.
#' @param pbar Shows progress bar, default TRUE.
#' @param phi.start Initial guess for phases, default NULL.
#' @param standardize Standardize the matrix, default FALSE.
#' @param GTEX_names Convert GTEx names format, default FALSE.
#' @param constrain_12h_window Enable 12h sampling window constraint, default FALSE.
#' @param window_start Start of sampling window in radians [0, 2π), default 0.
#' @param window_width Width of sampling window in radians, default 2π.
#' @param window_prior_strength Strength of window constraint, default 2.0.
#' @param post_hoc_hard_constrain Enable hard constrain post hoc to force all samples fall in designated time window, default FALSE
#' @param reference_genes Vector of gene names to use as phase references, default NULL.
#' @param reference_phases Expected peak phases for reference genes in radians, default NULL.
#' @param reference_troughs Expected trough phases for reference genes in radians, default NULL (auto-calculated as peak + π).
#' @param reference_strength Strength of reference phase constraint, default 0.1.
#' @param use_trough_constraint Whether to use trough time constraints in addition to peak, default TRUE.
#' @param auto_adjust_anti_phase Automatically detect and correct π-shift (anti-phase) after convergence, default TRUE.
#' @param reference_sample_indices Indices of samples corresponding to reference genes (if applicable), default NULL.
#' @param n_cores Number of CPU cores to use for parallel processing, default NULL (auto-detect).
#' @return List containing inferred phases, parameters, and convergence information.
#' @examples
#' # Basic usage with window constraint only
#' result1 <- CHIRAL_dual_penalty(data, constrain_12h_window=TRUE, window_start=pi/4, window_width=pi)
#'
#' # With both window and reference constraints (peak + trough)
#' ref_genes <- c("PER1", "CLOCK", "BMAL1")
#' ref_phases <- c(pi/2, pi, 3*pi/2)  # Peaks: 6AM, 12PM, 6PM
#' ref_troughs <- c(3*pi/2, 0, pi/2)  # Troughs: 6PM, 12AM, 6AM (optional, auto-calculated if NULL)
#' result2 <- CHIRAL_dual_penalty(data,
#'                               constrain_12h_window=TRUE, window_start=0, window_width=pi,
#'                               reference_genes=ref_genes, reference_phases=ref_phases,
#'                               reference_troughs=ref_troughs,  # Optional
#'                               window_prior_strength=2.0, reference_strength=1.5,
#'                               use_trough_constraint=TRUE)
#' @export
if (requireNamespace("tictoc", quietly = TRUE)) library(tictoc)
library(foreach)
library(doParallel)

CHIRAL_dual_penalty <- function(E, iterations=500, clockgenes=NULL, tau2=NULL, u=NULL,
                                sigma2=NULL, TSM=TRUE, mean.centre.E=TRUE, q=0.1, update.q=FALSE, pbar=TRUE,
                                phi.start=NULL, standardize=FALSE, GTEx_names=FALSE,
                                constrain_12h_window=FALSE, window_start=0, window_width=2*pi, window_prior_strength=2.0,
                                post_hoc_hard_constrain=FALSE,
                                reference_genes=NULL, reference_phases=NULL, reference_troughs=NULL,
                                reference_strength=0.1, use_trough_constraint=TRUE, auto_adjust_anti_phase=TRUE,
                                reference_sample_indices=NULL, n_cores=NULL, 
                                reference_decay_enabled=TRUE, decay_midpoint=0.05, decay_sharpness=100) {

  # Setup parallel backend
  if(is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1  # Leave one core free
  }
  if(n_cores <= 1) {
    foreach::registerDoSEQ()
  } else {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))  # Ensure cluster is stopped when function exits
  }

  # Load required functions from original CHIRAL script

  # Validate inputs
  if(!is.null(reference_genes) && !is.null(reference_phases)) {
    if(length(reference_genes) != length(reference_phases)) {
      stop("reference_genes and reference_phases must have same length")
    }

    # Auto-calculate trough phases if not provided (peak + π)
    if(is.null(reference_troughs) && use_trough_constraint) {
      reference_troughs <- (reference_phases + pi) %% (2*pi)
      message("Auto-calculated trough phases as peak + π")
    }

    # Validate trough phases if provided
    if(!is.null(reference_troughs) && use_trough_constraint) {
      if(length(reference_troughs) != length(reference_genes)) {
        stop("reference_troughs must have same length as reference_genes")
      }
    }
  }

  # E is the data matrix - genes on rows, samples on columns
  id=as.vector(c(1,1,1,1))
  E=as.matrix(E)

  # Validate window parameters
  if(constrain_12h_window && window_width > pi) {
    warning("12h window constraint requested but window_width > π. Setting window_width = π")
    window_width = pi
  }

  # Normalize window parameters
  window_start = window_start %% (2*pi)
  window_end = (window_start + window_width) %% (2*pi)

  E=E[ , colSums(is.na(E)) == 0]
  if(standardize==TRUE){E=sweep(E,1,apply(E,1,sd),FUN="/")}
  E.full=E

  ### Clock gene selection ###
  rownames(E)=toupper(rownames(E))
  if(is.null(clockgenes)){
    if(is.null(rownames(E))){
      rownames(E)=1:nrow(E)
    }
    clockgenes=rownames(E)
  }
  clockgenes=toupper(clockgenes)
  gene.list=rownames(E)

  E.full=E

  if(GTEx_names){
    gene.list=gsub("^.*_", "",gene.list)
    gene.list=gsub("\\|.*$","", gene.list)
  }
  clock.coord=match(clockgenes,gene.list)
  clock.coord=clock.coord[!is.na(clock.coord)]
  E=E[clock.coord,]
  geni=gene.list[clock.coord]

  # Prepare reference gene mapping
  reference_gene_indices <- NULL
  if(!is.null(reference_genes)) {
    reference_genes <- toupper(reference_genes)
    reference_gene_indices <- match(reference_genes, geni)
    reference_gene_indices <- reference_gene_indices[!is.na(reference_gene_indices)]

    if(length(reference_gene_indices) == 0) {
      warning("No reference genes found in the dataset. Reference penalty will be ignored.")
      reference_genes <- NULL
      reference_phases <- NULL
    } else {
      # Adjust reference_phases to match found genes
      valid_ref_mask <- !is.na(match(reference_genes, geni))
      if(sum(valid_ref_mask) < length(reference_genes)) {
        warning(paste("Some reference genes not found:",
                      paste(reference_genes[!valid_ref_mask], collapse=", ")))
        reference_genes <- reference_genes[valid_ref_mask]
        reference_phases <- reference_phases[valid_ref_mask]
        if(!is.null(reference_troughs)) {
          reference_troughs <- reference_troughs[valid_ref_mask]
        }
      }
    }
  }

  if(is.matrix(E)==FALSE){return()}
  if(mean.centre.E==TRUE){E=sweep(E,1,rowMeans(E),FUN="-")}

  ### Initialization of parameters ###
  if(is.null(sigma2)){sigma2 = mean(apply(E, 1, var))}
  if(is.null(u)){u=0.2}
  if(is.null(tau2)){
    if(constrain_12h_window) {
      tau2 = 2/(12 + ncol(E)) * window_prior_strength
    } else {
      tau2 = 4/(24 + ncol(E))
    }
  }
  if(nrow(E)>500){tau2=tau2/10}
  phi=phi.start

  ### Phase initialization with dual constraints ###
  if(is.null(phi.start)){
    # Choose initialization strategy based on reference gene availability
    if(!is.null(reference_genes) && !is.null(reference_phases)) {
      if(constrain_12h_window) {
        message("### Using reference-based initialization WITH window constraint ###")
      } else {
        message("### Using reference-based initialization (no window constraint) ###")
      }

      # First, get fallback initialization using spin-glass
      beta=1000
      J=J.tilde(E)
      Zeta = Zeta.mf.ordered(J, beta, ncol(E))
      fallback_init <- Zeta[,2] + runif(ncol(E),-0.5,0.5)

      # Use reference gene expression to initialize phases
      # IMPORTANT: Pass window parameters to initialization function
      # This ensures candidates are ONLY tested within the valid sampling window
      phi_init <- initialize_phases_by_reference(
        E, reference_gene_indices, reference_phases, reference_troughs,
        constrain_window = constrain_12h_window,
        window_start = window_start,
        window_width = window_width,
        fallback_init = fallback_init
      )

      # No need for additional window constraint - already handled in initialization
      # But apply it anyway as a safety check
      if(constrain_12h_window) {
        phi_init <- constrain_phases_to_window(phi_init, window_start, window_width)
      }
    } else {
      message("### Using standard spin-glass initialization ###")
      # Standard initialization without reference genes
      beta=1000
      J=J.tilde(E)
      Zeta = Zeta.mf.ordered(J, beta, ncol(E))
      phi_init <- Zeta[,2] + runif(ncol(E),-0.5,0.5)

      # Apply window constraint if enabled
      if(constrain_12h_window) {
        phi_init <- constrain_phases_to_window(phi_init, window_start, window_width)
      }
    }

    phi <- phi_init %% (2*pi)
  }

  if(pbar){pb = txtProgressBar(min = 0, max = iterations, style=3)}

  Ng=length(E[,1])
  N=length(E[1,])
  Ns=N
  sigma2.0=sigma2
  T=diag(3)*tau2
  T[1,1]=u^2
  dTinv=1/det(T)
  sigma2.m0=apply(E,1,var)
  sigma2.m1=sigma2.m0

  phi.0=phi
  c=cos(phi)
  s=sin(phi)
  one=rep(1,Ns)
  X=cbind(one,c,s)
  curvature=matrix(0, nrow = 2, ncol = Ns)
  S=t(E)%*%E/Ng

  S<-list()
  for(l in 1:Ng){
    S[[l]]=E[l,]%o%E[l,]
  }

  W<-rep(1,nrow(E))
  W.0<-W

  # Initialize decay tracking
  decay_history <- numeric(iterations)
  current_decay_factor <- 1.0

  print('Start enhanced EM procedure')
  print(paste('Total iteration:',iterations))
  if(reference_decay_enabled && !is.null(reference_genes)) {
    print(paste('Reference penalty decay enabled: midpoint =', decay_midpoint, ', sharpness =', decay_sharpness))
  }
  ### Start of Enhanced EM procedure ###
  for(i in 1:iterations){
    phiold=phi
    Xold=X
    sigma2old=sigma2
    W.old=W

    Nn=(t(X)%*%X)
    Nninv<-solve(Nn)
    Tinv=solve(T)
    M=Nn+sigma2*Tinv
    Minv=solve(M)

    alpha=Minv%*%t(X)%*%t(E)
    alpha<-t(alpha)
    colnames(alpha)=c("mu","a","b")
    alphat=as.matrix(alpha)
    alpha=as.data.frame(alpha)
    if(TSM==FALSE){W<-W.0}
    Tot<-sum(W)

    A=sum(W*(alpha$a*alpha$a/sigma2+Minv[2,2]))/Tot
    B=sum(W*(alpha$a*alpha$b/sigma2+Minv[2,3]))/Tot
    C=sum(W*(alpha$b*alpha$a/sigma2+Minv[3,2]))/Tot
    D=sum(W*(alpha$b*alpha$b/sigma2+Minv[3,3]))/Tot

    K=matrix(c(A,B,C,D),nrow = 2,ncol = 2,byrow = T)
    if (any(is.nan(K))){return(list(alpha=alpha,weights=W,iteration=i))}

    al=apply(E,2, function(x){
      return(sum(W*(alpha$a*(x-alpha$mu)/sigma2-Minv[1,2]))/Tot)
    })
    be=apply(E,2, function(x){
      return(sum(W*(alpha$b*(x-alpha$mu)/sigma2-Minv[1,3]))/Tot)
    })
    O=as.matrix(rbind(al,be))
    if (any(is.nan(O))){return(list(alpha=alpha,weights=W,iteration=i))}

    ### Calculate sigmoid decay factor based on convergence (delta_phi) ###
    if(i > 1) {
      delta_phi <- max(abs(phi - phiold))
    } else {
      delta_phi <- 1.0  # First iteration: no delta_phi yet, use full penalty
    }

    # Sigmoid decay: factor = 1 / (1 + exp(sharpness * (midpoint - delta_phi)))
    # When delta_phi >> midpoint: factor ≈ 1 (full penalty, still exploring)
    # When delta_phi ≈ midpoint: factor ≈ 0.5 (transitioning)
    # When delta_phi << midpoint: factor ≈ 0 (no penalty, converged)
    if(reference_decay_enabled && !is.null(reference_genes)) {
      current_decay_factor <- 1 / (1 + exp(decay_sharpness * (decay_midpoint - delta_phi)))
    } else {
      current_decay_factor <- 1.0  # No decay if disabled
    }

    # Store decay history
    decay_history[i] <- current_decay_factor

    ### Enhanced root finding with dual penalties ###
    rooted=apply(O,2,function(x){
      zero=B*B*C*C+A*A*D*D-x[1]*x[1]*(D*D+C*C)-x[2]*x[2]*(A*A+B*B)+2*x[1]*x[2]*(A*B+C*D)-2*A*B*C*D
      one=2*((A+D)*(A*D-B*C)-x[1]*x[1]*D-x[2]*x[2]*A+x[1]*x[2]*(B+C))
      two=A*A+D*D+4*A*D-x[1]*x[1]-x[2]*x[2]-2*B*C
      three=2*(A+D)
      four=1
      roots=polyroot(c(zero,one,two,three,four))
      return(roots)
    })
    ### Test each root with dual penalty system (parallelized) ###
    ze=foreach(j=1:Ns, .combine='cbind', .packages=c('base'),
                .export=c('apply_window_constraint', 'calculate_window_penalty',
                          'calculate_reference_penalty_sample', 'reference_expression_score_sample',
                          'current_decay_factor')) %dopar% {
      possible=lapply(rooted[,j],function(y){
        if(abs(Im(y)) > 1.0e-8){return(c(0,0,100000,100000))}
        K.lambda=K+diag(2)*Re(y)
        zet=solve(K.lambda,O[,j])

        phit=atan2(zet[2],zet[1])

        # Apply soft window constraint if enabled
        if(constrain_12h_window) {
          phit <- apply_window_constraint(phit, window_start, window_width, window_prior_strength)
        }

        Xt=c(1,cos(phit),sin(phit))
        Mt=Xt%o%Xt

        Xt.old=c(1,cos(phiold[j]),sin(phiold[j]))
        Mt.old=Xt.old%o%Xt.old

        Qs=sapply(1:Ng,function(k){
          return(alphat[k,] %*% Mt %*% alphat[k,] - 2*alphat[k,]%*%Xt*E[k,j]+sigma2*sum(diag(Minv%*%Mt)))
        })

        Qs.old=sapply(1:Ng,function(k){
          return(alphat[k,] %*% Mt.old %*% alphat[k,] - 2*alphat[k,]%*%Xt.old*E[k,j]+sigma2*sum(diag(Minv%*%Mt.old)))
        })

        Q=sum(Qs*W/sigma2)/Tot
        Q.old=sum(Qs.old*W/sigma2)/Tot
        # Add window penalty if enabled
        if(constrain_12h_window) {
          window_penalty <- calculate_window_penalty(phit, window_start, window_width, window_prior_strength)
          Q <- Q + window_penalty
        }
        # Add reference phase penalty if enabled (with decay)
        if(!is.null(reference_genes) && !is.null(reference_phases)) {
          # Apply decayed penalty strength
          effective_strength <- reference_strength * current_decay_factor
          reference_penalty <- calculate_reference_penalty_sample(phit, j, reference_gene_indices,
                                                                  reference_phases, reference_troughs,
                                                                  E, alpha, effective_strength, use_trough_constraint)
          Q <- Q + reference_penalty
        }
        return(c(zet, Q, Q.old))
      })
      pox=matrix(unlist(possible),4,4,byrow = T)
      mpox=which.min(pox[,3])
      if(pox[mpox,3]==100000){
        cat("\n",pox,"not any solution on the circle at iteration",i,"\n")
        stop()
        return(list(phi=phiold,Qhist=Qhist,sigma=sigma2old, alpha=alpha, weights=W))
      }
      return(c(pox[mpox,],pox[,3]))
    }
    ### Update phases with dual constraints ###
    # ze is already a matrix from foreach with .combine='cbind'
    rownames(ze)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")

    phi_new = atan2(ze[2,],ze[1,]) %% (2*pi)

    # Apply hard constraints post-update
    if(constrain_12h_window & post_hoc_hard_constrain) {
      phi_new <- constrain_phases_to_window(phi_new, window_start, window_width)
    }

    phi <- phi_new

    # Store iteration history with decay factor
    if(i==1){
      Qhist=as.data.frame(t(ze))
      colnames(Qhist)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
      Qhist$iteration=i
      Qhist$sample=c(1:Ns)
      Qhist$decay_factor=current_decay_factor
    }
    else{
      Qtemp=data.frame(t(ze))
      colnames(Qtemp)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
      Qtemp$iteration=i
      Qtemp$sample=c(1:Ns)
      Qtemp$decay_factor=current_decay_factor
      Qhist=rbind(Qhist,Qtemp)}

    # Update weights
    P1.0=sapply(1:Ng, function(p){
      gamm=sqrt(dTinv*sigma2^3/det(M))
      exponent=E[p,]%*%X%*%Minv%*%t(X)%*%E[p,]/(2*sigma2)
      return(exp(exponent)*gamm)
    })
    W<-q*P1.0/(1-q+q*P1.0)
    W[is.nan(W)]<-1

    ### Check convergence ###
    if(max(abs(phi-phiold))<(0.01)){
      if(pbar){cat("\n algorithm has converged \n")}

      # Auto-adjust for potential anti-phase problem if reference genes are provided
      phi_adjusted <- phi
      if(auto_adjust_anti_phase && !is.null(reference_genes) && !is.null(reference_phases)) {
        if(constrain_12h_window) {
          cat("\n### Skipping post-hoc π-shift: 12h sampling window already constrains anti-phase ambiguity ###\n")
        } else {
          cat("\n### Post-hoc phase adjustment check ###\n")
          phi_adjusted <- adjust_phase_by_reference(
            phi, alpha, E, reference_gene_indices, reference_phases, reference_troughs, use_trough_constraint
          )
        }
      }

      return(list(phi=phi_adjusted, phi_raw=phi, sigma=sigma2, alpha=alpha, weights=W, iter=i, sigma.m1=sigma2.m1,
                  E=E.full, Qhist=Qhist, geni=geni, clock=E,
                  constrain_12h_window=constrain_12h_window, window_start=window_start, window_width=window_width,
                  reference_genes=reference_genes, reference_phases=reference_phases, reference_troughs=reference_troughs,
                  use_trough_constraint=use_trough_constraint, reference_strength=reference_strength,
                  auto_adjusted=!identical(phi, phi_adjusted),
                  decay_history=decay_history[1:i], reference_decay_enabled=reference_decay_enabled,
                  decay_midpoint=decay_midpoint, decay_sharpness=decay_sharpness))
    }

    ### Update parameters for next iteration ###
    c=cos(phi)
    s=sin(phi)
    one=rep(1,length(phi))
    X=cbind(one,c,s)

    Mold=Nn+sigma2*Tinv
    Moldinv=solve(Mold)

    sigma2.m1=sapply(1:Ng, function(s){
      return(sum(diag(S[[s]]-S[[s]]%*%Xold%*%Moldinv%*%t(X)))/Ns+0.01)
    })
    sigma2.m0=apply(E,1,var)
    sigma2<-mean(sigma2.m1*W+sigma2.m0*(1-W))
    sigma2.m1=sum(sigma2.m1*W)/sum(W)

    if(update.q){
      q=mean(W)
      if(q<0.05) q=0.05
      if(q>0.3) q=0.3}
    if(pbar){Sys.sleep(0.1)
      setTxtProgressBar(pb, i)}
  }

  if(pbar){close(pb); cat("\n")}

  # Auto-adjust for potential anti-phase problem if reference genes are provided (max iterations case)
  phi_adjusted <- phi
  if(auto_adjust_anti_phase && !is.null(reference_genes) && !is.null(reference_phases)) {
    if(constrain_12h_window) {
      cat("\n### Skipping post-hoc π-shift: 12h sampling window already constrains anti-phase ambiguity ###\n")
    } else {
      cat("\n### Post-hoc phase adjustment check (max iterations reached) ###\n")
      phi_adjusted <- adjust_phase_by_reference(
        phi, alpha, E, reference_gene_indices, reference_phases, reference_troughs, use_trough_constraint
      )
    }
  }

  return(list(phi=phi_adjusted, phi_raw=phi, Qhist=Qhist, sigma=sigma2, alpha=alpha, weights=W, iter=i,
              sigma.m1=sigma2.m1, E=E.full,
              constrain_12h_window=constrain_12h_window, window_start=window_start, window_width=window_width,
              reference_genes=reference_genes, reference_phases=reference_phases, reference_troughs=reference_troughs,
              use_trough_constraint=use_trough_constraint, reference_strength=reference_strength,
              auto_adjusted=!identical(phi, phi_adjusted),
              decay_history=decay_history, reference_decay_enabled=reference_decay_enabled,
              decay_midpoint=decay_midpoint, decay_sharpness=decay_sharpness))
}

#' Calculate penalty for deviation from reference gene phases (peaks and troughs)
#' @param phi_sample Current phase estimate for a sample
#' @param sample_idx Index of current sample
#' @param ref_gene_indices Indices of reference genes in the gene list
#' @param ref_phases Expected peak phases for reference genes
#' @param ref_troughs Expected trough phases for reference genes (can be NULL)
#' @param E Expression matrix
#' @param alpha Current gene parameters
#' @param strength Penalty strength
#' @param use_trough Whether to use trough constraint
calculate_reference_penalty_sample <- function(phi_sample, sample_idx, ref_gene_indices, ref_phases, ref_troughs,
                                                E, alpha, strength, use_trough=TRUE) {
  if(is.null(ref_gene_indices) || length(ref_gene_indices) == 0) {
    return(0)
  }

  return(strength * reference_expression_score_sample(
    phi_sample = phi_sample,
    sample_idx = sample_idx,
    ref_gene_indices = ref_gene_indices,
    ref_phases = ref_phases,
    ref_troughs = ref_troughs,
    E = E,
    use_trough = use_trough
  ))
}

#' Score reference-gene expression against a candidate sample phase.
#'
#' Lower scores indicate that high/low observed reference-gene expression is
#' consistent with the candidate phase and the supplied peak/trough references.
reference_expression_score_sample <- function(phi_sample, sample_idx, ref_gene_indices, ref_phases,
                                              ref_troughs = NULL, E, use_trough = TRUE) {
  total_score <- 0

  for(i in 1:length(ref_gene_indices)) {
    gene_idx <- ref_gene_indices[i]
    gene_expr_all_samples <- as.numeric(E[gene_idx, ])
    actual_expr <- as.numeric(E[gene_idx, sample_idx])

    if(length(unique(gene_expr_all_samples)) <= 1 || is.na(actual_expr)) {
      next
    }

    actual_percentile <- mean(gene_expr_all_samples <= actual_expr, na.rm = TRUE)

    peak_distance <- abs(((phi_sample - ref_phases[i] + pi) %% (2*pi)) - pi)
    expected_peak_expression <- (cos(peak_distance) + 1) / 2
    total_score <- total_score + (actual_percentile - expected_peak_expression)^2

    if(use_trough && !is.null(ref_troughs)) {
      trough_distance <- abs(((phi_sample - ref_troughs[i] + pi) %% (2*pi)) - pi)
      expected_trough_state <- (cos(trough_distance) + 1) / 2
      expected_low_expression <- 1 - expected_trough_state
      total_score <- total_score + 0.7 * (actual_percentile - expected_low_expression)^2
    }
  }

  return(total_score)
}

#' Apply reference gene bias during initialization
#' @param phi_init Initial phase estimates
#' @param ref_gene_indices Indices of reference genes
#' @param ref_phases Expected phases for reference genes
#' @param E Expression matrix
#' @param bias_strength Strength of bias (0-1)
apply_reference_bias <- function(phi_init, ref_gene_indices, ref_phases, E, bias_strength=0.3) {
  if(is.null(ref_gene_indices) || length(ref_gene_indices) == 0) {
    return(phi_init)
  }

  # For each sample, bias toward phases where reference genes would be expressed as expected
  for(s in 1:length(phi_init)) {
    reference_score <- 0
    total_weight <- 0

    for(i in 1:length(ref_gene_indices)) {
      gene_idx <- ref_gene_indices[i]
      expected_phase <- ref_phases[i]

      # Score how well current phase aligns with expected reference gene expression
      phase_diff <- abs(((phi_init[s] - expected_phase + pi) %% (2*pi)) - pi)
      alignment_score <- cos(phase_diff)  # Higher when phases align

      expr_level <- abs(E[gene_idx, s]) + 0.1
      reference_score <- reference_score + alignment_score * expr_level
      total_weight <- total_weight + expr_level
    }

    if(total_weight > 0) {
      avg_reference_score <- reference_score / total_weight

      # If reference score is low, bias toward nearest reference phase
      if(avg_reference_score < 0.5) {
        # Find nearest reference phase
        distances <- sapply(ref_phases, function(rp) {
          abs(((phi_init[s] - rp + pi) %% (2*pi)) - pi)
        })
        nearest_ref_phase <- ref_phases[which.min(distances)]

        # Bias toward nearest reference phase
        phi_init[s] <- phi_init[s] * (1 - bias_strength) + nearest_ref_phase * bias_strength
      }
    }
  }

  return(phi_init %% (2*pi))
}

#' Utility function to convert hours to phases
#' @param hours Time in hours [0, 24)
#' @return Phases in radians [0, 2π]
hours_to_phases <- function(hours) {
  phi <- (hours / 24) * 2 * pi
  return(phi %% (2*pi))
}

#' Utility function to convert phases to hours
#' @param phi Phases in radians [0, 2π]
#' @return Time in hours [0, 24)
phases_to_hours <- function(phi) {
  hours <- (phi / (2*pi)) * 24
  return(hours)
}


# Backward compatibility functions - original CHIRAL functions kept unchanged

#' Original CHIRAL function for backward compatibility
CHIRAL<- function(E, iterations=500, clockgenes=NULL,tau2=NULL, u=NULL,
                  sigma2=NULL, TSM=TRUE, mean.centre.E=TRUE, q=0.1, update.q=FALSE, pbar=TRUE,
                  phi.start=NULL, standardize=FALSE, GTEx_names=FALSE){

  # Call the new function with no window constraints for backward compatibility
  return(CHIRAL_dual_penalty(E, iterations=iterations, clockgenes=clockgenes, tau2=tau2, u=u,
                    sigma2=sigma2, TSM=TSM, mean.centre.E=mean.centre.E, q=q,
                    update.q=update.q, pbar=pbar, phi.start=phi.start,
                    standardize=standardize, GTEx_names=GTEx_names,
                    constrain_12h_window=FALSE))
}

cor.c<- function(x,y=NULL){
  if(is.null(y)){return(cor.c.g(x))}
  if(length(x)!=length(y)){
    cat("the two vectors have different length")
    return(NULL)
  }
  n=length(x)
  U1=0
  U2=0
  U3=0
  for(i in 1:(n-1)){
    for(j in ((i+1):n)){
      U1=U1+sin(x[i]-x[j])*sin(y[i]-y[j])
      U2=U2+sin(x[i]-x[j])^2
      U3=U3+sin(y[i]-y[j])^2
    }
  }
  corr=U1/sqrt(U2*U3)
  return(corr)
}

cor.c.g<-function(x){
  off=ncol(x)
  cr=matrix(0, ncol=off, nrow=off)
  dimnames(cr)=list(colnames(x), colnames(x))
  for(l in 1:off){
    for(m in l:off){
      cr[l,m]=cr[m,l]=cor.c(x[,l], x[,m])
    }
  }
  return(cr)
}

delta.phi<-function(phi.0, phi,period=2*pi, mode="forgotten", N=200, median_scale=1){
  phi.0=phi.0%%(period)
  phi=phi%%(period)
  isP=abs(period-2*pi)>0.001
  if(!isP){
    obj=calc.delta(phi.0, phi, N)
    bestphi=obj$phi
    mad=obj$median*median_scale
  }
  else{
    obj=calc.delta(phi.0/period*2*pi, phi/period*2*pi, N)
    bestphi=obj$phi*period/2/pi
    mad=obj$median*period/2/pi*median_scale
  }
  if(mode=="say"){
    cat("median:", mad, "\n")
    return(bestphi)}
  else if(mode=="return"){
    return(list(phi=bestphi, median=mad))
  }
  else if(mode=="no_median"){
    return(bestphi)
  }
  else{
    cat("median:", mad, "\n")
    return(list(phi=bestphi, median=(mad)))
  }
}

calc.delta<-function(phi.0, phi, N=200){
  mad<-12
  sdel=NULL
  for(j in 1:N){
    offset<-j/N*2*pi
    theta<-(phi-offset)%%(2*pi)
    del<-abs(theta-phi.0)%%(2*pi)
    delta<-del
    for (i in 1:length(phi)) {
      delta[i]=min(del[i], 2*pi-del[i])
    }
    if(median(delta)<mad){
      mad<-median(delta)
      bestphi<-theta
      j_temp=j
      sdel=delta}
    phi<-(-phi)%%(2*pi)
    theta<-(phi-offset)%%(2*pi)
    del<-abs(theta-phi.0)%%(2*pi)
    delta<-del
    for (i in 1:length(phi)) {
      delta[i]=min(del[i], 2*pi-del[i])
    }
    if(median(delta)<mad){
      mad<-median(delta)
      bestphi<-theta
      j_temp=-j
      sdel=delta}
  }
  return(list(phi=bestphi, median=(mad)))
}

adjust.phases<-function(realphi, infphi, period=2*pi){
  if(abs(period-2*pi)>0.001){return(adjust.phases(realphi*2*pi/period, infphi*2*pi/period)*period/pi/2)}
  for(i in 1:length(realphi)){
    if(realphi[i]-infphi[i]>pi){infphi[i]=2*pi+infphi[i]}
    if(realphi[i]-infphi[i]<(-pi)){infphi[i]=infphi[i]-2*pi}
  }
  return(infphi)
}

### Spin glass approximation to have initial condition for EM ###

Zeta.mf.ordered<-function(J, beta, samples, A.0=0.1){
  iterations<-1000
  time_symmetry<-0
  A<-rep(A.0, samples)
  Theta<-runif(samples, 0, 2*pi)
  for (time in 1: iterations) {
    A.c=A*cos(Theta)
    A.s=A*sin(Theta)
    for (k in 1: samples) {
      u=beta*sum(A.c*J[k,])
      v=beta*sum(A.s*J[k,])

      modulo<-sqrt(u*u+v*v)
      Zeta_k<-c(u/modulo, v/modulo)
      Theta[k]<-atan2(Zeta_k[2],Zeta_k[1])
      A[k]<-(besselI(modulo,1)/(besselI(modulo,0)))
      if(is.na(modulo)){A[k] <- 1}
      else if(modulo>20){A[k] <- 1}
      if(is.nan(u) || is.nan(v)){Theta[k]<-runif(1,0,2*pi)}
      else if(u==0){
        if(v==0){
          Theta[k]<-runif(1,0,2*pi)
        }
      }
    }
  }
  for (k in 1:(samples-1)) {
    if (Theta[k] < Theta[k+1]){
      time_symmetry<-time_symmetry+1
    }
    if (Theta[k] > Theta[k+1]){
      time_symmetry<-time_symmetry-1
    }
  }

  time_symmetry<-sign(time_symmetry)
  if(time_symmetry==-1){
    for (k in 1:samples) {
      Theta[k]<-2*pi-Theta[k]
    }
  }
  Theta<-Theta%%(2*pi)
  Zeta<-cbind(A, Theta)
  return(Zeta)
}

### Calculation of spin glass interaction matrix ###

J.tilde<-function(E,n.genes=0,n.samples=0){
  sda<-as.matrix(E)
  if(n.genes==0){n.genes=nrow(E)}
  if(n.samples==0){n.samples=ncol(E)}
  Jtilde<-matrix(0,ncol = n.samples, nrow = n.samples)
  for(i in 1:n.samples){
    for (j in 1:n.samples) {
      Jtilde[i,j]<-sum(sda[,i]*sda[,j])/(n.samples*n.genes)
    }
  }
  diag(Jtilde)<-0
  rownames(Jtilde)<-colnames(E)
  colnames(Jtilde)<-colnames(E)
  return(Jtilde)
}

#' Apply soft window constraint during optimization
apply_window_constraint <- function(phi, window_start, window_width, strength) {
  # Soft constraint - don't force hard boundaries during optimization
  # Let the optimization explore but bias toward window

  window_center <- (window_start + window_width/2) %% (2*pi)

  # Calculate distance from window center
  dist_from_center <- abs(((phi - window_center + pi) %% (2*pi)) - pi)

  # If outside acceptable range, bias toward window center
  max_allowed_dist <- window_width/2 * 1.2  # Allow some flexibility

  if(dist_from_center > max_allowed_dist) {
    # Bias toward window center
    bias_factor <- 0.8
    phi <- phi * (1 - bias_factor) + window_center * bias_factor
  }

  return(phi %% (2*pi))
}

#' Calculate penalty for being outside sampling window
calculate_window_penalty <- function(phi, window_start, window_width, strength) {
  window_end <- (window_start + window_width) %% (2*pi)

  # Calculate distance outside window using circular distance
  if(window_end < window_start) {
    # Wrap-around case (e.g., 10PM to 10AM)
    in_window <- (phi >= window_start) | (phi <= window_end)
    if(in_window) return(0)

    # Calculate circular distance to nearest window boundary
    dist_to_start <- pmin(abs(phi - window_start),
                          abs(phi - window_start + 2*pi),
                          abs(phi - window_start - 2*pi))
    dist_to_end <- pmin(abs(phi - window_end),
                        abs(phi - window_end + 2*pi),
                        abs(phi - window_end - 2*pi))
    dist_out <- pmin(dist_to_start, dist_to_end)
  } else {
    # Normal case (e.g., 6AM to 6PM)
    # Check if phi is in window using circular logic
    if(phi >= window_start && phi <= window_end) return(0)

    # Calculate circular distance to nearest window boundary
    dist_to_start <- pmin(abs(phi - window_start),
                          abs(phi - window_start + 2*pi),
                          abs(phi - window_start - 2*pi))
    dist_to_end <- pmin(abs(phi - window_end),
                        abs(phi - window_end + 2*pi),
                        abs(phi - window_end - 2*pi))
    dist_out <- pmin(dist_to_start, dist_to_end)
  }

  # Quadratic penalty
  penalty <- strength * dist_out^2
  return(penalty)
}

#' Post-hoc phase adjustment to fix potential π-shift (anti-phase) problem
#'
#' This function checks EACH SAMPLE independently to see if it should be shifted by π
#' (12h) to better align with reference gene expectations.
#'
#' @param phi Inferred phases in radians [0, 2π]
#' @param alpha Gene coefficient matrix (Ng x 3, columns: mu, a, b)
#' @param E Expression matrix (genes x samples)
#' @param ref_gene_indices Indices of reference genes
#' @param ref_phases Expected peak phases for reference genes
#' @param ref_troughs Expected trough phases for reference genes (can be NULL)
#' @param use_trough Whether to use trough constraint
#' @return Adjusted phases that better align with reference genes (per-sample decision)
#' @examples
#' # After running CHIRAL_dual_penalty:
#' result <- CHIRAL_dual_penalty(data, reference_genes=..., reference_phases=...)
#' adjusted_phi <- adjust_phase_by_reference(
#'   result$phi, result$alpha, result$clock,
#'   ref_gene_indices, ref_phases, ref_troughs
#' )
adjust_phase_by_reference <- function(phi, alpha, E, ref_gene_indices, ref_phases,
                                      ref_troughs = NULL, use_trough = TRUE) {

  if(is.null(ref_gene_indices) || length(ref_gene_indices) == 0) {
    message("No reference genes provided. Returning original phases.")
    return(phi)
  }

  n_samples <- length(phi)
  phi_adjusted <- phi
  n_adjusted <- 0

  # Evaluate each sample independently
  for(s in 1:n_samples) {
    phi_original <- phi[s]
    phi_shifted <- (phi[s] + pi) %% (2*pi)

    # Calculate alignment scores for this sample
    score_original <- calculate_sample_reference_score(
      phi_original, s, alpha, E, ref_gene_indices, ref_phases, ref_troughs, use_trough
    )

    score_shifted <- calculate_sample_reference_score(
      phi_shifted, s, alpha, E, ref_gene_indices, ref_phases, ref_troughs, use_trough
    )

    # Choose the version with better alignment for THIS sample
    if(score_shifted < score_original) {
      phi_adjusted[s] <- phi_shifted
      n_adjusted <- n_adjusted + 1
    }
  }

  # Report per-sample results
  if(n_adjusted > 0) {
    pct_adjusted <- (n_adjusted / n_samples) * 100
    message(sprintf(
      "[Per-sample] Adjusted %d/%d samples (%.1f%%) by π-shift",
      n_adjusted, n_samples, pct_adjusted
    ))
  } else {
    message("[Per-sample] No individual samples required adjustment")
  }

  # NOW: Check if shifting ALL samples together is even better (global adjustment)
  message("\n[Global check] Testing if shifting all samples by π improves overall alignment...")

  phi_after_individual <- phi_adjusted
  phi_global_shifted <- (phi_adjusted + pi) %% (2*pi)

  # Calculate total scores for both options
  score_individual <- 0
  score_global_shift <- 0

  for(s in 1:n_samples) {
    score_individual <- score_individual + calculate_sample_reference_score(
      phi_after_individual[s], s, alpha, E, ref_gene_indices, ref_phases, ref_troughs, use_trough
    )
    score_global_shift <- score_global_shift + calculate_sample_reference_score(
      phi_global_shifted[s], s, alpha, E, ref_gene_indices, ref_phases, ref_troughs, use_trough
    )
  }

  # Choose the better option
  if(score_global_shift < score_individual) {
    improvement <- ((score_individual - score_global_shift) / score_individual) * 100
    message(sprintf(
      "[Global] Shifting ALL samples by π improves alignment by %.1f%%. Applying global shift.",
      improvement
    ))
    message(sprintf(
      "[Global] Score after individual adjustment: %.3f, Score after global shift: %.3f",
      score_individual, score_global_shift
    ))
    return(phi_global_shifted)
  } else {
    message(sprintf(
      "[Global] No global shift needed. Individual adjustments are optimal."
    ))
    message(sprintf(
      "[Global] Score after individual adjustment: %.3f, Score if globally shifted: %.3f",
      score_individual, score_global_shift
    ))
    return(phi_after_individual)
  }
}

#' Calculate reference alignment score for a SINGLE sample
#'
#' Helper function to score how well a sample's reference gene expression
#' aligns with expected peak/trough patterns at a given phase.
#'
#' @param phi_sample Phase to evaluate for this sample
#' @param sample_idx Index of the sample
#' @param alpha Gene coefficients
#' @param E Expression matrix
#' @param ref_gene_indices Indices of reference genes
#' @param ref_phases Expected peak phases
#' @param ref_troughs Expected trough phases (can be NULL)
#' @param use_trough Whether to use trough constraint
#' @return Misalignment score for this sample (lower = better)
calculate_sample_reference_score <- function(phi_sample, sample_idx, alpha, E,
                                              ref_gene_indices, ref_phases,
                                              ref_troughs, use_trough) {
  return(reference_expression_score_sample(
    phi_sample = phi_sample,
    sample_idx = sample_idx,
    ref_gene_indices = ref_gene_indices,
    ref_phases = ref_phases,
    ref_troughs = ref_troughs,
    E = E,
    use_trough = use_trough
  ))
}

#' Initialize phases based on reference gene expression patterns
#'
#' This function uses reference gene expression levels to directly infer
#' initial phase estimates, avoiding the anti-phase trap. When window constraints
#' are provided, it restricts candidate phases to the valid sampling window.
#'
#' @param E Expression matrix (genes x samples)
#' @param ref_gene_indices Indices of reference genes
#' @param ref_phases Expected peak phases for reference genes
#' @param ref_troughs Expected trough phases (optional)
#' @param constrain_window Whether to constrain to sampling window (default FALSE)
#' @param window_start Start of sampling window in radians (default 0)
#' @param window_width Width of sampling window in radians (default 2π)
#' @param fallback_init Fallback initialization if reference-based fails
#' @return Vector of initial phase estimates [0, 2π]
#' @examples
#' # Without window constraint
#' phi_init <- initialize_phases_by_reference(
#'   E, ref_gene_indices, ref_phases, ref_troughs
#' )
#'
#' # With window constraint (e.g., daytime sampling only)
#' phi_init <- initialize_phases_by_reference(
#'   E, ref_gene_indices, ref_phases, ref_troughs,
#'   constrain_window = TRUE, window_start = 0, window_width = pi
#' )
initialize_phases_by_reference <- function(E, ref_gene_indices, ref_phases,
                                          ref_troughs = NULL,
                                          constrain_window = FALSE,
                                          window_start = 0, window_width = 2*pi,
                                          fallback_init = NULL) {

  n_samples <- ncol(E)
  phi_init <- numeric(n_samples)

  # Generate candidate phases
  if(constrain_window && window_width < 2*pi) {
    # IMPORTANT: If window constraint is enabled, only test phases WITHIN the window
    # This represents the biological reality: samples were collected during a limited time
    # E.g., patient samples collected only during clinic hours (8am-8pm)

    message(sprintf(
      "Constraining initialization to window [%.2f, %.2f] radians (%.1f-%.1f hours)",
      window_start, (window_start + window_width) %% (2*pi),
      window_start * 12 / pi, ((window_start + window_width) %% (2*pi)) * 12 / pi
    ))

    # Generate candidates within window only (every 15 degrees)
    n_candidates <- max(8, ceiling(window_width / (pi/12)))  # At least 8 candidates
    candidate_phases <- seq(window_start, window_start + window_width, length.out = n_candidates)
    # NOTE: Do NOT apply modulo here - keep phases in their original range
    # This preserves proper ordering for windows that may extend beyond 2π
  } else {
    # Test full 24-hour range (every 15 degrees = π/12 = 1 hour)
    candidate_phases <- seq(0, 2*pi - 0.01, by = pi/12)
  }

  # Strategy: For each sample, find the phase where reference genes
  # would have expression patterns most similar to observed

  for(s in 1:n_samples) {
    # Get expression levels of reference genes in this sample
    ref_expr <- E[ref_gene_indices, s]

    best_score <- Inf
    best_phase <- 0

    for(phi_test in candidate_phases) {
      score <- 0

      # For each reference gene, calculate expected vs observed expression similarity
      for(i in 1:length(ref_gene_indices)) {
        expected_peak <- ref_phases[i]
        actual_expr <- ref_expr[i]

        # Expected expression pattern: high near peak, low near trough
        # Use cosine similarity to expected peak time
        peak_distance <- abs(((phi_test - expected_peak + pi) %% (2*pi)) - pi)

        # Expected relative expression (0 to 1 scale)
        # When φ = expected_peak: expected_rel_expr = 1
        # When φ = expected_peak ± π: expected_rel_expr = 0
        expected_rel_expr <- (cos(peak_distance) + 1) / 2

        # Normalize actual expression to [0, 1] using percentile rank
        # This handles different gene expression levels
        gene_expr_all_samples <- E[ref_gene_indices[i], ]
        actual_percentile <- sum(gene_expr_all_samples <= actual_expr) / n_samples

        # Score is squared difference between expected and observed
        score <- score + (expected_rel_expr - actual_percentile)^2

        # Add trough constraint if provided
        if(!is.null(ref_troughs)) {
          expected_trough <- ref_troughs[i]
          trough_distance <- abs(((phi_test - expected_trough + pi) %% (2*pi)) - pi)
          expected_rel_expr_trough <- (cos(trough_distance) + 1) / 2

          # At trough, expression should be LOW (percentile should be low)
          # So we want: when expected_rel_expr_trough is high → actual should be low
          score <- score + 0.5 * (expected_rel_expr_trough - (1 - actual_percentile))^2
        }
      }

      # Choose phase with best (lowest) score
      if(score < best_score) {
        best_score <- score
        best_phase <- phi_test
      }
    }

    phi_init[s] <- best_phase
  }

  # Add small random noise to avoid exact duplicates
  phi_init <- phi_init + runif(n_samples, -0.1, 0.1)
  phi_init <- phi_init %% (2*pi)

  # Validate initialization
  if(any(is.na(phi_init)) || any(is.nan(phi_init))) {
    warning("Reference-based initialization failed. Using fallback.")
    if(!is.null(fallback_init)) {
      return(fallback_init)
    } else {
      return(runif(n_samples, 0, 2*pi))
    }
  }

  message(sprintf(
    "Initialized phases using reference gene expression patterns (range: %.2f to %.2f radians)",
    min(phi_init), max(phi_init)
  ))

  return(phi_init)
}

#' Helper function to constrain phases to sampling window
#' @param phi Vector of phases in [0, 2π]
#' @param window_start Start of window in radians
#' @param window_width Width of window in radians
#' @return Constrained phases
constrain_phases_to_window <- function(phi, window_start, window_width) {
  window_end <- (window_start + window_width) %% (2*pi)

  # Handle wrap-around case (e.g., 10PM to 10AM)
  if(window_end < window_start) {
    # Window wraps around midnight
    in_window <- (phi >= window_start) | (phi <= window_end)
    phi[!in_window] <- find_nearest_window_phase(phi[!in_window], window_start, window_end, wrap=TRUE)
  } else {
    # Normal case (e.g., 6AM to 6PM)
    # BUT: Special handling when window_start is near 0 (< 0.1 radians)
    # Phases near 2π should be treated as equivalent to 0
    if(window_start < 0.1) {
      # For phases very close to 2π, check if they're circularly within the window
      # E.g., if window is [0, π] and φ=6.2, check if (φ - 2π) would be in window
      phi_normalized <- phi %% (2*pi)
      alternative_phi <- ifelse(phi_normalized > 1.5*pi, phi_normalized - 2*pi, phi_normalized)
      in_window <- ((phi_normalized >= window_start) & (phi_normalized <= window_end)) |
                    ((alternative_phi >= window_start) & (alternative_phi <= window_end))
    } else {
      in_window <- (phi >= window_start) & (phi <= window_end)
    }
    phi[!in_window] <- find_nearest_window_phase(phi[!in_window], window_start, window_end, wrap=FALSE)
  }

  return(phi)
}


#' Find nearest phase within the sampling window
find_nearest_window_phase <- function(phi, window_start, window_end, wrap=FALSE) {
  # ALWAYS use circular distance, even for non-wrapping windows
  # This is critical when window_start=0, as phases near 2π are very close to 0

  # Circular distance calculation
  dist_to_start <- pmin(abs(phi - window_start), abs(phi - window_start + 2*pi), abs(phi - window_start - 2*pi))
  dist_to_end <- pmin(abs(phi - window_end), abs(phi - window_end + 2*pi), abs(phi - window_end - 2*pi))

  nearest <- ifelse(dist_to_start < dist_to_end, window_start, window_end)

  return(nearest)
}
