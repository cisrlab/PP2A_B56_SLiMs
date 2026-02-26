#' Title
#'
#' @param Sequences
#' @param min.width
#' @param debug
#'
#' @return
#' @export
#'
#' @importFrom alakazam charge
#'
#' @examples
getPosChargeFeatures<-function(Sequences, min.width = 9, debug=FALSE) {
    charge_features = NULL;

    meme_Ba = rep(0, length(Sequences));

    balpha_w8_path = system.file("extdata", "balpha_w8.xml", package="PP2A.B56.SLiMs");
    balpha_w7_path = system.file("extdata", "balpha_w7.xml", package="PP2A.B56.SLiMs");
    balpha_w6_path = system.file("extdata", "balpha_w6.xml", package="PP2A.B56.SLiMs");



    for (idx in 1:length(Sequences)) {
        #if (debug) {cat(idx, " ", Sequences[idx],"\n")}

        #    ba_fimo = try(meme$fimo(seq=Sequences[idx], meme_xml=balpha_w6_path, thresh=1.0, no_qvalue = TRUE));



        #    if (!("try-error" %in% attributes(ba_fimo)) && nrow(ba_fimo) > 0) {
        #      ba_fimo = ba_fimo[which(ba_fimo$p.value == min(ba_fimo$p.value))[1],];
        #      aa_pos = strsplit(ba_fimo$motif.sequence,"")[[1]];
        #    } else {
        begin.idx = 1;
        end.idx = min(nchar(Sequences[idx]), min.width)
        sseq = substr(Sequences[idx], begin.idx, end.idx);
        sseq = pad_sequence(sseq, min.width);
        aa_pos = strsplit(sseq,"")[[1]];
        #    }
        if (debug) {cat(idx," ",aa_pos,"\n")}
        if (debug) {cat(idx," ",charge(aa_pos),"\n")}
        charge_features = rbind(charge_features, charge(aa_pos));
    }
    charge_features = as.data.frame(charge_features, stringsAsFactors=TRUE);
    #rownames(charge_features) = Sequences;
    colnames(charge_features) = paste0("Charge",1:ncol(charge_features));


    return(charge_features);

}
