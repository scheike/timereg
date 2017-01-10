context("dutils")

test_that("dsort", {
    data("hubble",package="lava")
    expect_equivalent(order(dsort(hubble, ~sigma)$sigma),
                      seq_len(nrow(hubble)))
})

test_that("dsort", {
    data("hubble",package="lava")
    h1 <- hubble
    drename(h1,fun=toupper) <- ~.    
    expect_equivalent(colnames(h1),toupper(names(hubble)))
})
