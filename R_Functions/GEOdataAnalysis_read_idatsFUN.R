
read_idats <- function(idat_files,quiet=FALSE){
  
J = length(idat_files)
zipped = !file.exists(paste0(idat_files, "_Grn.idat"))
suffix = rep(".idat", times = J)
suffix[zipped] = ".idat.gz"
ex = file.exists(paste0(idat_files, "_Grn", suffix)) & file.exists(paste0(idat_files, 
                                                                          "_Red", suffix))
if (!all(ex)) 
  stop("Some .idat files are missing")
P = illuminaio::readIDAT(paste0(idat_files[1], "_Grn", suffix[1]))$nSNPsRead
print(P)
if (P %in% c(1051815, 1051943, 1052641)) {
  platform = "EPIC"
  manifest = data.table::copy(ewastools:::manifest_epic)
  controls = data.table::copy(ewastools:::controls_epic)
}
else if (P == 622399) {
  platform = "450K"
  manifest = data.table::copy(ewastools:::manifest_450K)
  controls = data.table::copy(ewastools:::controls_450K)
}
else {
  stop("Unknown platform")
  platform = "UNK"
}
setkeyv(manifest, c("chr", "mapinfo"))
manifest[, `:=`(index, 1L:.N)]
controls[, `:=`(index, 1L:.N)]
manifest[channel == "Grn", `:=`(OOBi, 1:.N)]
manifest[channel == "Red", `:=`(OOBi, 1:.N)]
M = U = matrix(NA_real_, nrow = nrow(manifest), ncol = J)
S = T = matrix(NA_integer_, nrow = nrow(manifest), ncol = J)
N = V = matrix(NA_integer_, nrow = nrow(manifest), ncol = J)
ctrlG = ctrlR = matrix(NA_real_, nrow = nrow(controls), ncol = J)
ctrlN = matrix(NA_integer_, nrow = nrow(controls), ncol = J)
oobG = list(M = matrix(NA_real_, nrow = manifest[channel == 
                                                   "Red", .N], ncol = J), U = matrix(NA_real_, nrow = manifest[channel == 
                                                                                                                 "Red", .N], ncol = J))
oobR = list(M = matrix(NA_real_, nrow = manifest[channel == 
                                                   "Grn", .N], ncol = J), U = matrix(NA_real_, nrow = manifest[channel == 
                                                                                                                 "Grn", .N], ncol = J))
if (!quiet) 
  pb <- txtProgressBar(min = 0, max = J, style = 3)
barcodes = rep(NA_character_, J)
dates = rep(NA_character_, J)
for (j in 1:J) {
  red = illuminaio::readIDAT(paste0(idat_files[j], "_Red", 
                                    suffix[j]))
  grn = illuminaio::readIDAT(paste0(idat_files[j], "_Grn", 
                                    suffix[j]))
  idat_order = red$MidBlock
  if (!identical(idat_order, grn$MidBlock)) 
    stop("Red and green .idat files do not agree!")
  barcodes[j] = red$Barcode
  if (nrow(red$RunInfo) > 1) 
    dates[j] = red$RunInfo[2, 1]
  manifest[, `:=`(Ui, match(addressU, idat_order))]
  manifest[, `:=`(Mi, match(addressM, idat_order))]
  controls[, `:=`(i, match(address, idat_order))]
  setindexv(manifest, "channel")
  manifest["Both", `:=`(Mi, Ui), on = "channel"]
  i = manifest["Red", on = "channel"]
  U[i$index, j] = red$Quants[i$Ui, 1]
  T[i$index, j] = red$Quants[i$Ui, 2]
  V[i$index, j] = red$Quants[i$Ui, 3]
  M[i$index, j] = red$Quants[i$Mi, 1]
  S[i$index, j] = red$Quants[i$Mi, 2]
  N[i$index, j] = red$Quants[i$Mi, 3]
  oobG$U[i$OOBi, j] = grn$Quants[i$Ui, 1]
  oobG$M[i$OOBi, j] = grn$Quants[i$Mi, 1]
  i = manifest["Grn", on = "channel"]
  U[i$index, j] = grn$Quants[i$Ui, 1]
  T[i$index, j] = grn$Quants[i$Ui, 2]
  V[i$index, j] = grn$Quants[i$Ui, 3]
  M[i$index, j] = grn$Quants[i$Mi, 1]
  S[i$index, j] = grn$Quants[i$Mi, 2]
  N[i$index, j] = grn$Quants[i$Mi, 3]
  oobR$U[i$OOBi, j] = red$Quants[i$Ui, 1]
  oobR$M[i$OOBi, j] = red$Quants[i$Mi, 1]
  i = manifest["Both", on = "channel"]
  U[i$index, j] = red$Quants[i$Ui, 1]
  T[i$index, j] = red$Quants[i$Ui, 2]
  V[i$index, j] = red$Quants[i$Ui, 3]
  M[i$index, j] = grn$Quants[i$Mi, 1]
  S[i$index, j] = grn$Quants[i$Mi, 2]
  N[i$index, j] = grn$Quants[i$Mi, 3]
  ctrlR[controls$index, j] = red$Quants[controls$i, 1]
  ctrlG[controls$index, j] = grn$Quants[controls$i, 1]
  ctrlN[controls$index, j] = red$Quants[controls$i, 3]
  if (!quiet) 
    setTxtProgressBar(pb, j)
}
if (!quiet) 
  close(pb)
M[N == 0] = NA
U[V == 0] = NA
S[N == 0 | N == 1] = NA
T[V == 0 | V == 1] = NA
ctrlG[ctrlN == 0] = NA
ctrlR[ctrlN == 0] = NA
sample_ids = strsplit(x = idat_files, split = "/")
sample_ids = sapply(sample_ids, tail, n = 1L)
meta = data.table(sample_id = sample_ids, date = as.IDate(dates, 
                                                          "%m/%d/%Y %r"), time = as.ITime(dates, "%m/%d/%Y %r"))
raw = list(platform = platform, manifest = manifest, U = U, 
           T = T, V = V, M = M, S = S, N = N, controls = controls, 
           ctrlG = ctrlG, ctrlR = ctrlR, ctrlN = ctrlN, oobG = oobG, 
           oobR = oobR, meta = meta)
return(raw)
}
