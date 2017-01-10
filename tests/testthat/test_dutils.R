context("dutils")

test_that("dsort", {
    data(hubble)
    expect_equivalent(order(dsort(hubble, ~sigma)$sigma),
                      seq_len(nrow(hubble)))
})
