test_that("2024_06_27_cv works", {

    # Load the original training predictions
    itc.path <- system.file(
        "extdata",
        "res_tr.m107_rem1_95_all_rfsrc_log.2024_06_27_cv.model.RData",
        package="PP2A.B56.SLiMs"
    )
    vars <- load(itc.path)

    seq <- names(itc_train_predict)
    orig_predict <- as.numeric(itc_train_predict)

    #Test to see if prediction are still the same and what was found before
    new_predict <- as.numeric(getITC_2024_06_27_cv(seq))
    print(new_predict)
    for (idx in seq_len(length(seq))) {
        expect_equal(
            new_predict[idx], orig_predict[idx],
            tolerance = 0.001
        )
    }
})
