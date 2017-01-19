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


test_that("daggregate", {    
    dd <- data.frame(a=1:20,b=20:1,
                    g1=rep(0:1,10),
                    g2=rep(0:1,each=10),
                    g3=rbinom(20,1,0.5))
    dd$g1[1:2] <- NA
    dd$g2[2:3] <- NA
    dd$a[3:7] <- NA
    ##dcor(dd,.~g1+g2)
    
    dcor(dd, "[acg][12]*$", regex=TRUE, use="pairwise")
    dcor(dd, "[ab]",subset=is.na(g1), regex=TRUE, use="pairwise")
    dcor(dd, ~.|is.na(g1)|is.na(g2))
    dcor(dd, use="pairwise")
})


