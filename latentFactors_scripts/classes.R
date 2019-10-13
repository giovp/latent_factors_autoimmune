setClass("latentNeighbors", slots=list(
  factors.df = "data.table",
  loadings.df = "data.table",
  diss = "data.table",
  clust.memb = "numeric",
  latent.cons.df = "data.table",
  kCluster = "integer",
  kLatent = "integer",
  samp = "integer",
  amari.dist = "numeric",
  entropy = "numeric",
  mean.svcca = "numeric",
  sil="numeric",
  loss="numeric",
  corrToMeanPatt="numeric"
),
prototype = list(
  diss = data.table::data.table(),
  kLatent = integer(),
  latent.cons.df = data.table::data.table()
)
)

classShow <- function(object){
  allSlots <- slotNames(object)
  for (slt in allSlots){
    clss <- class(slot(object, slt))
    cat(paste0("Slot: ", slt),"\n")
    if(clss=="data.table"){
      vals <- dplyr::trunc_mat(slot(object, slt), n = 3, n_extra = 2)
    } else if(clss=="integer"){
      vals <- sort(unique(slot(object, slt)))
    } else if(clss=="numeric"){
      vals <- length(slot(object, slt))
    }
    print(vals)
    cat("\n")
  }
}

setMethod("show","latentNeighbors",classShow)
