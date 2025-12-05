test_that("2024_06_27_cv works", {


    #Test to see if prediction are still the same and what was found before
    seq <- c(
        "WTSFFSG.CSPIEEEAH",
        "K.L(pS)PIIED(pS)",
        "LDTLRETQE",
        "LSIKK.LSPIIEDSREA",
        "WTSFFSG.LSPIE(Sp)EAH",
        "L(pS)PIIEDDREADH"
    )
    orig_predict <- c(
        49.0743950,
        2.2218900,
        5.9678875,
        6.7478491,
        57.6646175,
        1.4119304
    )

    new_predict <- as.numeric(getITC_2024_06_27_cv(seq))
    for (idx in seq_len(length(seq))) {
        expect_equal(
            orig_predict[idx], new_predict[idx],
            tolerance = 0.001
        )
    }
})
