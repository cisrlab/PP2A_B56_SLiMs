getMemeFeatures<-function(Sequences, use_internal=!have_ext_fimo(), debug=FALSE) {
    SeqLength = nchar(Sequences);
    meme_Ba = rep(0, length(Sequences))
    meme_Ba_idx = rep(-1, length(Sequences))

    balpha_w8_path = system.file("extdata", "balpha_w8.xml", package="PP2A.B56.SLiMs");
    balpha_w7_path = system.file("extdata", "balpha_w7.xml", package="PP2A.B56.SLiMs");
    balpha_w6_path = system.file("extdata", "balpha_w6.xml", package="PP2A.B56.SLiMs");

    for (idx in 1:length(Sequences)) {
        ba_fimo = try(
          meme$fimo(
            seq=Sequences[idx],
            meme_xml=balpha_w6_path,
            thresh=1.0,
            no_qvalue = TRUE,
            use_internal = use_internal,
            debug = debug
          )
        )
        if (!"try-error" %in% attributes(ba_fimo) && nrow(ba_fimo) > 0) {
            meme_Ba[idx] = pvalue.to.zscore(min(ba_fimo$p.value, na.rm=TRUE));
            meme_Ba_idx[idx] = ba_fimo$motif.start[ba_fimo$p.value == min(ba_fimo$p.value)][1]
        }
    }
    ans = data.frame(
        meme_Ba   = round(meme_Ba, digits=6),
        meme_Ba_idx = meme_Ba_idx,
        stringsAsFactors=FALSE
    );
    #rownames(ans) = Sequences;
    NAs = is.na(ans);
    ans[NAs] = 0;
    return(ans);
}
