
AMINOACIDSYMBOLS <- c(
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
    "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U", "O",
    "B", "J", "Z", "X"
)


getBaseSequence<-function(sequence, remove_dots=TRUE) {
  ans = sequence;
  if (remove_dots) {
    ans = gsub("\\.","", ans);
  }
  ans = gsub("\\(p","", ans);
  ans = gsub("\\)","", ans);
  return(ans);
}

pvalue.to.zscore <-function(p, one.sided=TRUE, log.p=FALSE) {
  ans = qnorm(p, lower.tail = FALSE, log.p = log.p)
  if (one.sided) {
    ans[ans < 0] = 0;
  }
  return(ans);
}



replacePhosWithE<-function(sequences) {

  ans = sequences
  for (seq_idx in 1:length(sequences)) {
    sequence = sequences[seq_idx]

    sequence = gsub("\\(pS\\)", "E", sequence)
    sequence = gsub("\\(pT\\)", "E", sequence)
    sequence = gsub("\\(pY\\)", "E", sequence)
    base_seq = getBaseSequence(sequence);
    ans[seq_idx] = sequence;
  }
  return(ans)
}

#' Get Alakazam peptide property features
#'
#' @param Sequences peptide aa sequence
#'
#' @return matrix with alakazam values in column and sequence in rows
#' @export
#'
#' @importFrom alakazam gravy
#' @importFrom alakazam bulk
#' @importFrom alakazam polar
#' @importFrom alakazam aliphatic
#' @importFrom alakazam charge
#' @importFrom alakazam countPatterns
#' @examples
#' getPeptidePropertyFeatures("LDTLRETQE")
getPeptidePropertyFeatures<-function(Sequences) {
  ans = data.frame(
    sequence_length = nchar(gsub("X","",Sequences)), #Length of sequence
    gravy = gravy(Sequences), # Grand average of hydrophobicity
    bulk  = bulk(Sequences), # Average bulkiness
    polar = polar(Sequences), # Average polarity
    aliphatic = aliphatic(Sequences), # Normalized aliphatic index
    charge = charge(Sequences), # Normalized Net charge;
    stringsAsFactors=FALSE
  )

  ans = cbind(
    ans,
    countPatterns(
        Sequences,
        nt=FALSE, c(ACIDIC="[DE]", BASIC="[RHK]", PHOS="[STY]",
                    AROMATIC="[FWHY]"))) #ACIDIC and Phosphorl Sites aas
  return(ans)
}


getPosMemeFeatures<-function(Sequences) {


  motif_features = NULL;

  meme_Ba = rep(0, length(Sequences));

  balpha_w8_path = system.file("extdata", "balpha_w8.xml", package="PP2A.B56.SLiMs");
  balpha_w7_path = system.file("extdata", "balpha_w7.xml", package="PP2A.B56.SLiMs");
  balpha_w6_path = system.file("extdata", "balpha_w6.xml", package="PP2A.B56.SLiMs");

  for (idx in 1:length(Sequences)) {
    ba_fimo = meme$fimo(seq=Sequences[idx], meme_xml=balpha_w6_path, thresh=1.0, no_qvalue = TRUE);


    if (nrow(ba_fimo) > 0) {
           ba_fimo = ba_fimo[which(ba_fimo$p.value == min(ba_fimo$p.value))[1],];

      aa_pos = strsplit(ba_fimo$motif.sequence,"")[[1]];

      motif_features = rbind(motif_features,aa_pos);

    }
  }
  motif_features = as.data.frame(motif_features, stringsAsFactors=TRUE);
  #rownames(motif_features) = Sequences;
  colnames(motif_features) = paste0("Pos",1:ncol(motif_features));


  return(motif_features);

}

pad_sequence<-function(seq, width) {
  n = nchar(seq);
  if (n < width) {
    add = width - n;
    padding = paste(rep("X",add), collapse="", sep="")
    seq = paste0(seq, padding);
  }
  return(seq);
}

pad_sequences<-function(seqs, width) {
  ans = seqs;
  for (idx in 1:length(seqs)) {
    ans[idx] = pad_sequence(seqs[idx], width)
  }
  return(ans);
}





getPhosPositions<-function(Sequences) {
  ans = list();
  #Sequences = gsub("\\.","", Sequences);
  for (row_idx in 1:length(Sequences)) {
    pos = unlist(gregexpr("\\(p", Sequences[row_idx]));
    if (length(pos) > 1) {
      for (idx in 2:length(pos)) {
        pos[idx] = pos[idx] - (idx-1) * 3
      }
    }
    ans[[Sequences[row_idx]]] = pos;
  }
  return(ans);
}


getPosMemeFeatures2 <-function (Sequences, max.pos = 9, debug=FALSE)
{

    aa_levels = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U", "O", "B", "J", "Z", "X")

    motif_features = NULL
    meme_Ba = rep(0, length(Sequences))
    balpha_w8_path = system.file("extdata", "balpha_w8.xml",
        package = "PP2A.B56.SLiMs")

    balpha_w7_path = system.file("extdata", "balpha_w7.xml",
      package = "PP2A.B56.SLiMs")

    balpha_w6_path = system.file("extdata", "balpha_w6.xml", package="PP2A.B56.SLiMs");

    for (idx in 1:length(Sequences)) {
      #if (debug) {cat(idx, " ",Sequences[idx],"\n")}
        ba_fimo = try(meme$fimo(seq = Sequences[idx], meme_xml = balpha_w6_path,
            thresh = 1.0, no_qvalue = TRUE));

        begin.idx = 1;
        if (!("try-error" %in% attributes(ba_fimo)) && nrow(ba_fimo) > 0) {
            ba_fimo = ba_fimo[which(ba_fimo$p.value == min(ba_fimo$p.value))[1],
                ]
            begin.idx = ba_fimo$motif.start;
        } else {
          begin.idx = 1;
        }
        current_seq = rep("X", max.pos)
        end.idx = min(begin.idx+max.pos-1,nchar(Sequences[idx]));
        clength = end.idx - begin.idx+1
        #if (debug) {cat("*clength:",clength); cat("current seq:",current_seq,"\n");}
        for (cidx in 1:clength) {
          current_seq[cidx] = substr(Sequences[idx],begin.idx+cidx-1,begin.idx+cidx-1)
        }
        #if (debug) {cat("*current seq:",current_seq,"\n");}

#            return(current_seq);
            aa_pos = current_seq
            motif_features = rbind(motif_features, aa_pos)
    }
    motif_features = as.data.frame(motif_features, stringsAsFactors = TRUE)

    for (idx in 1:ncol(motif_features)) {
      motif_features[,idx] = factor(motif_features[,idx], levels = aa_levels)
    }


    colnames(motif_features) = paste0("Pos", 1:ncol(motif_features))
    return(motif_features)
}

#' Get Feature data set for passed in amino acid sequences
#'
#' @param sequences peptide aa sequence
#'
#' @return feature matrix of features used to predict ITC values
#' @export
#'
#' @examples
#' getFS("LDTLRETQE")
getFS<-function(sequences) {
  Sequence = replacePhosWithED(sequences);
  pep_fs = getPeptidePropertyFeatures(Sequence);
  meme_fs = getMemeFeatures(Sequence);
  meme_pos_fs = getPosMemeFeatures2(Sequence);
  charge_pos_fs = getPosChargeFeatures(Sequence);
  mer_fs = getMerFeatures(
      pad_sequences(Sequence,9),
      min_mer=1, max_mer=1,
      filter.redundant=FALSE,
      filter.sd=FALSE)/pep_fs$sequence_length
  aa_levels = c("A", "R", "N", "D",
                "C", "Q", "E", "G",
                "H", "I", "L", "K",
                "M", "F", "P", "S",
                "T", "W", "Y", "V",
                "U", "O", "B", "J", "Z", "X"
  )
  for (aa in aa_levels) {
    if (!aa %in% colnames(mer_fs)) {
      mer_fs[,aa] = 0;
    }
  }
  mer_fs = mer_fs[,aa_levels];

  #Rule 1 - E or D at 2nd or 7th location, 10x stronger
  POS2_ED = meme_pos_fs$Pos2 == "E" | meme_pos_fs$Pos2 == "D"
  POS7_ED = meme_pos_fs$Pos7 == "E" | meme_pos_fs$Pos7 == "D"
  Rule1 = POS2_ED | POS7_ED

  POS8_ED = meme_pos_fs$Pos8 == "E" | meme_pos_fs$Pos8 == "D"
  POS9_ED = meme_pos_fs$Pos9 == "E" | meme_pos_fs$Pos9 == "D"


  #Rule 2 - E or D at 8th or 9th location, 2-3x stronger
  Rule2 = POS8_ED | POS9_ED

  #Rule 3 - K or R at 2nd or 7th location, 10x weaker
  POS2_KR = meme_pos_fs$Pos2 == "K" | meme_pos_fs$Pos2 == "R"
  POS7_KR = meme_pos_fs$Pos7 == "K" | meme_pos_fs$Pos7 == "R"
  Rule3 = POS2_KR | POS7_KR

  #Rule 4 - K or R at 8th or 9th location, 2-3x weaker
  POS8_KR = meme_pos_fs$Pos8 == "K" | meme_pos_fs$Pos8 == "R"
  POS9_KR = meme_pos_fs$Pos9 == "K" | meme_pos_fs$Pos9 == "R"
  Rule4 = POS8_KR | POS9_KR

  #5) L is standard for the 1st locaion, M or V or F or Y at the 1st, 2-4x weaker, I at the 1st, 2x weaker, other residues in the above list at the 1st, >4x weaker

  POS1_L = meme_pos_fs$Pos1 == "L"

  POS1_MVFY = meme_pos_fs$Pos1 %in% c("M","V","F","Y")

  POS1_I = meme_pos_fs$Pos1 == "I"

  POS1_Other = !meme_pos_fs$Pos1 %in% c("L","M","V","F","Y","I")

  Rule5 = rep("Unknown", nrow(meme_pos_fs))
  Rule5[POS1_L] = "L";
  Rule5[POS1_MVFY] = "MVFY";
  Rule5[POS1_I] = "I"
  Rule5[POS1_Other] = "Other"

  Rule5 = factor(Rule5, levels=c("L","I","MVFY","Other"))

  #6) I is standard for the  4th locaion, M or V or F or Y at the 4th, 2-4x weaker, L at the 4th, 2x weaker, other residues in the above list at the 4th, >4x weaker

  POS4_I = meme_pos_fs$Pos4 == "I"
  POS4_MVFY = meme_pos_fs$Pos4 %in% c("M","V","F","Y")
  POS4_L = meme_pos_fs$Pos4 == "L"
  POS4_Other = !meme_pos_fs$Pos4 %in% c("I", "M", "V", "F", "Y", "L")

  Rule6 = rep("Unknown", nrow(meme_pos_fs))
  Rule6[POS4_I] = "I";
  Rule6[POS4_MVFY] = "MVFY";
  Rule6[POS4_L] = "L"
  Rule6[POS4_Other] = "Other"

  Rule6 = factor(Rule6, levels = c("L","I","MVFY","Other"))

  #Rule 7 - 6th location could only be E. Even D would lose binding.
  Rule7 = meme_pos_fs$Pos6 != "E" #Calculate coefficient when Pos != "E", should be much weaker.


  data_df = data.frame(
  POS2_ED = POS2_ED,
  POS7_ED = POS7_ED,
  Rule1 = Rule1,
  POS8_ED = POS8_ED,
  POS9_ED = POS9_ED,
  Rule2 = Rule2,
  Rule3 = Rule3,
  Rule4 = Rule4,
  Rule5 = Rule5,
  Rule6 = Rule6,
  Rule7 = Rule7
  );


  new_fs = cbind(Sequence, pep_fs, meme_fs, mer_fs, meme_pos_fs, charge_pos_fs, data_df);
  return(new_fs);

}

getITCByPath<-function(sequences, itc.path) {
  new_fs = getFS(sequences);
  load(itc.path);

  library(randomForest)
  new_litc = rForestTest(litc_model$model, new_fs, litc_features)

  return(2^new_litc$test.predict);

}



getTestFxn<-function(regression_algorithm) {
  if (regression_algorithm == "rforest") {
    return(rForestTest);
  }
  if (regression_algorithm == "rfsrc")
  {
    return (rfsrcTestReg);
  }
  if (regression_algorithm == "rfsrc_w1")
  {
    return(rfsrcTestReg);
  }
  if (regression_algorithm == "rfsrc_w2")
  {
    return(rfsrcTestReg);
  }
  if (regression_algorithm == "rfsrc_wp") {
    return(rfsrcTestReg);
  }
  if (regression_algorithm == "M5P") {
    return(M5PTest);
  }
  if (regression_algorithm == "cforest") {
    return(function(train_data, test_data, class_name, feature_names, ...) {
      return(cForestTest(train_data, test_data, class_name, feature_names, ntree=500, regression=TRUE,...))
    });
  }
  if (regression_algorithm == "glmnet") {
    return(glmNetTestReg);
  }

  if (regression_algorithm == "glmnet_w1") {
    return(glmNetTestReg);
  }

  if (regression_algorithm == "glmnet_w2") {
    return(glmNetTestReg);
  }
  if (regression_algorithm == "glm") {
    return(glmTestReg);
  }
  if (regression_algorithm == "lm_wp") {
    return(lmTestWPReg);
  }
  if (regression_algorithm == "lmi_wp") {
    return(lmTestWPReg);
  }

  stop("Unknown regression:", regression_algorithm);

}

getITCByPath2<-function(sequences, itc.path, enforce_6E = FALSE, debug = FALSE) {
  vars = load(itc.path);

  if (debug){print(vars);}

  #Use the getFeatureSet function if provided in model file.
  if ("itc_get_fs_fxn" %in% vars) {
      if (debug) {cat("getting FS using custom function\n");}
      new_fs = itc_get_fs_fxn(sequences);
  } else {
      if (debug) {cat("getting FS using standard route\n");}
      new_fs = getFS(sequences);
  }
  if (!("itc_model" %in% vars)) {stop("missing ITC model\n")}

  test_fxn = getTestFxn(itc_regressor)

  p = test_fxn(itc_model, new_fs, itc_features);
  if (debug) {print(p);}
  new_itc = undoTransform(p$test.predict, itc_transform);

  #Apply caps to prevent really bad predictions.
  if ("itc_top_cap" %in% vars) {
    new_itc[new_itc > itc_top_cap] = itc_top_cap;
  }
  if ("itc_bottom_cap" %in% vars) {
    new_itc[new_itc < itc_bottom_cap] = itc_bottom_cap;
  }

  if (enforce_6E) {
    #message("Enforce 6E")
    new_itc[new_fs$Clip_Pos6_E == "No"] <- 95
  }

  if (debug) {
    if ("itc_label" %in% vars) {
      attr(new_itc, "itc_label") = itc_label;
    } else {
      attr(new_itc, "itc_label") = basename(itc_path);
    }
  }

  if (debug) {print(new_itc);}
  if (debug) {attr(new_itc, "fs") <- new_fs}
  return(new_itc);
}



getITC<-function(sequences, model2=TRUE) {
  new_fs <- getFS(sequences)

  if (model2) {
    itc.path = system.file("extdata", "litc.model2.RData", package="PP2A.B56.SLiMs");

  } else {
    itc.path = system.file("extdata", "litc.model.RData", package="PP2A.B56.SLiMs");
  }
  return(getITCByPath(sequences, itc.path))

}

#' Predict ITC using random forest model
#'
#' @param sequences amino acid peptide sequences
#' @param debug return more debugging info
#' @return vector of itc values
#' @export
#'
#' @examples
#' getITC_2024_06_27_cv("LDTLRETQE")
getITC_2024_06_27_cv <- function(sequences, enforce_6E = FALSE, debug=FALSE) {
    itc.path <- system.file(
        "extdata",
        "res_tr.m107_rem1_95_all_rfsrc_log.2024_06_27_cv.model.RData",
        package="PP2A.B56.SLiMs"
    )
    return(getITCByPath2(sequences = sequences, itc.path = itc.path, debug = debug, enforce_6E = enforce_6E))
}


getITCHasValue<-function(sequences, itc.path = system.file("extdata", "litc.hasvalue.model.RData", package="PP2A.B56.SLiMs")) {
  new_fs = getFS(sequences);
  load(itc.path);
  library(randomForest)
  new_litc = rForestTest(litc_hasvalue_model$model,
                                          new_fs, litc_hasvalue_fn
                                          )

  ans = new_litc$test.prob[,"TRUE"];
  attr(ans, "res") = new_litc;
  return(ans);
}

#'Gets the full sequence, replacing the phosphorylations with the E/D equivalents
getFullSequence<-function(sequences, remove_dots=FALSE) {
  ans = gsub(" ","", sequences); #Remove blank spaces
  if (remove_dots) { #Remove the dots if requested
    ans = gsub("\\.", "", ans)
  }
  ans = replacePhosWithED(ans)
  return(ans);
}



getClipSequence<-function(sequences) {
  full_sequences = getFullSequence(sequences, remove_dots=FALSE);
  ans = unlist(
    lapply(
      strsplit(full_sequences, "\\."),
      function(l){
        if (length(l)==1) {return(l)}
        else {return(l[2])}
      }
    )
  )
  return(ans);
}



# (pS) => E, (Sp) => D
replacePhosWithED<-function(sequences, debug=FALSE) {

  ans = sequences;
  for (seq_idx in 1:length(sequences)) {

    sequence = sequences[seq_idx];
    if (debug){cat(seq_idx, " ", sequence);}

    sequence = gsub("\\(Sp\\)", "D", sequence)
    sequence = gsub("\\(Tp\\)", "D", sequence)
    sequence = gsub("\\(Yp\\)", "D", sequence)

    sequence = gsub("\\(pS\\)", "E", sequence)
    sequence = gsub("\\(pT\\)", "E", sequence)
    sequence = gsub("\\(pY\\)", "E", sequence)
    if (debug) {cat(" => ", sequence,"\n");}
    ans[seq_idx] = sequence;
  }

  return(ans);
}

getPhosPositionsE = getPhosPositions;
getPhosPositionsD <- function(Sequences) {
  ans = list()
  for (row_idx in 1:length(Sequences)) {
    posS = unlist(gregexpr("\\(S", Sequences[row_idx]))
    if (length(posS) > 1) {
      for (idx in 2:length(posS)) {
        posS[idx] = posS[idx] - (idx - 1) * 3
      }
    }
    posT = unlist(gregexpr("\\(T", Sequences[row_idx]))
    if (length(posT) > 1) {
      for (idx in 2:length(posT)) {
        posT[idx] = posT[idx] - (idx-1) * 3
      }
    }
    posY = unlist(gregexpr("\\(Y", Sequences[row_idx]))
    if (length(posY) > 1) {
      for (idx in 2:length(posY)) {
        posY[idx] = posY[idx] - (idx-1) * 3
      }
    }
    pos = unique(c(posS, posT, posY))

    pos = pos[order(pos, decreasing = FALSE)]
    ans[[Sequences[row_idx]]] = pos
  }
  return(ans)
}

getBaseSequence2<-function(sequence, remove_dots = TRUE) {
  ans = sequence
  if (remove_dots) {
    ans = gsub("\\.", "", ans)
  }
  ans = gsub(" ", "", ans)
  ans = gsub("\\(p", "", ans)
  ans = gsub("\\(Sp", "S", ans)
  ans = gsub("\\(Tp", "T", ans)
  ans = gsub("\\(Yp", "Y", ans)
  ans = gsub("\\)", "", ans)
  return(ans)
}


clipSequence<-function(base_seq, full_seq) {

  ans = full_seq;

  dots = stringr::str_locate_all(base_seq, "\\.")

  for (idx in 1:length(dots)) {
    current = dots[[idx]]
    if (nrow(current) == 1) {
      dot_loc = current[1,1]
      ans[idx] = substr(ans[idx], dot_loc, nchar(ans[idx]))
    } else if (nrow(current) > 1) {
      print(current)
      stop("Implement!")
    }
  }

  return(ans);

}




#' get feature set for clip and full sequence
#'
#' @param sequences amino acid sequences
#' This function returns a set of features for the clipped sequence (between the .'s)
#' And the full sequence (removing the dots).  We also replace the (pSTY) => E and the
#' (STYp) with D.
#'
#'
#' @return data.frame with features for both the clipped and full sequence
#' @export
#'
#' @examples
#' getFS2b("LDTLRETQE")
getFS2b<-function(sequences) {


  full_sequences = getFullSequence(sequences, remove_dots=TRUE)

  fs_full <- PP2A.B56.SLiMs:::getFS(full_sequences)
  colnames(fs_full) = paste0("Full_", colnames(fs_full));

  clip_sequences = getClipSequence(sequences);


  fs_clip = PP2A.B56.SLiMs::getFS(clip_sequences);
  colnames(fs_clip) = paste0("Clip_", colnames(fs_clip));

  fs = cbind(fs_clip,
fs_full[,c("Full_Sequence","Full_sequence_length","Full_gravy","Full_bulk","Full_polar",
"Full_aliphatic", "Full_charge")]);

  clip_pos_indices = grep("Clip_Pos", colnames(fs))

  #Create features for each position for a particular amino acid.
  # I.e. Clip_Pos1_E?

  aas = AMINOACIDSYMBOLS

  for (clip_pos_idx in clip_pos_indices) {
    for (aa in aas) {
      label = paste0(colnames(fs)[clip_pos_idx], "_", aa);
      fs[,label] = fs[,clip_pos_idx] == aa;
    }
  }

  for (fn in extractFN(fs)) {
      if (is.logical(fs[,fn])) {
          bf <- rep("No", nrow(fs))
          bf[fs[,fn]] <- "Yes"
          fs[,fn] <- factor(bf, levels = c("No", "Yes"))
      }
  }
  return(fs);
}


