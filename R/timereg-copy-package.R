

#' The Bone Marrow Transplant Data
#' 
#' Bone marrow transplant data with 408 rows and 5 columns.
#' 
#' 
#' @format The data has 408 rows and 5 columns. \describe{ \item{cause}{a
#' numeric vector code.  Survival status. 1: dead from treatment related
#' causes, 2: relapse , 0: censored.} \item{time}{ a numeric vector. Survival
#' time.  } \item{platelet}{a numeric vector code. Plalelet 1: more than 100 x
#' \eqn{10^9} per L, 0: less.} \item{tcell}{a numeric vector. T-cell depleted
#' BMT 1:yes, 0:no.} \item{age}{a numeric vector code. Age of patient, scaled
#' and centered ((age-35)/15).} }
#' @references NN
#' @name bmt
#' @docType data
#' @source Simulated data
#' @keywords package 
#' @examples
#' 
#' data(bmt)
#' names(bmt)
#' 
NULL





#' The multicenter AIDS cohort study
#' 
#' CD4 counts collected over time.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{obs}{a numeric vector. Number of observations.} \item{id}{a numeric
#' vector. Id of subject.} \item{visit}{ a numeric vector. Timings of the
#' visits in years.} \item{smoke}{a numeric vector code. 0: non-smoker, 1:
#' smoker.} \item{age}{a numeric vector. Age of the patient at the start of the
#' trial.} \item{cd4}{a numeric vector. CD4 percentage at the current visit.}
#' \item{cd4.prev}{a numeric vector. CD4 level at the preceding visit.}
#' \item{precd4}{a numeric vector. Post-infection CD4 percentage.} \item{lt}{a
#' numeric vector. Gives the starting time for the time-intervals.} \item{rt}{a
#' numeric vector. Gives the stopping time for the time-interval.} }
#' @references Kaslow et al. (1987), The multicenter AIDS cohort study:
#' rational, organisation and selected characteristics of the participants.
#' Am. J. Epidemiology 126, 310--318.
#' @source MACS Public Use Data Set Release PO4 (1984-1991). See reference.
#' @name cd4 
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(cd4)
#' names(cd4)
#' 
NULL





#' CSL liver chirrosis data
#' 
#' Survival status for the liver chirrosis patients of Schlichting et al.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{id}{ a numeric vector. Id of subject.  }
#' 
#' \item{time}{ a numeric vector. Time of measurement.  } \item{prot}{ a
#' numeric vector. Prothrombin level at measurement time.  } \item{dc}{ a
#' numeric vector code. 0: censored observation, 1: died at eventT.  }
#' \item{eventT}{ a numeric vector. Time of event (death).  } \item{treat}{ a
#' numeric vector code. 0: active treatment of prednisone, 1: placebo
#' treatment.  } \item{sex}{ a numeric vector code. 0: female, 1: male.  }
#' \item{age}{ a numeric vector. Age of subject at inclusion time subtracted
#' 60.  } \item{prot.base}{ a numeric vector. Prothrombin base level before
#' entering the study.  } \item{prot.prev}{ a numeric vector. Level of
#' prothrombin at previous measurement time.  } \item{lt}{ a numeric vector.
#' Gives the starting time for the time-intervals.  } \item{rt}{ a numeric
#' vector. Gives the stopping time for the time-intervals.  } }
#' @references Schlichting, P., Christensen, E., Andersen, P., Fauerholds, L.,
#' Juhl, E., Poulsen, H. and Tygstrup, N. (1983), The Copenhagen Study Group
#' for Liver Diseases, Hepatology 3, 889--895
#' @name  csl
#' @docType data
#' @keywords package 
#' @source P.K. Andersen
#' @examples
#' 
#' data(csl)
#' names(csl)
#' 
NULL





#' The Diabetic Retinopathy Data
#' 
#' The data was colleceted to test a laser treatment for delaying blindness in
#' patients with dibetic retinopathy. The subset of 197 patiens given in Huster
#' et al. (1989) is used.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{id}{a numeric vector. Patient code.} \item{agedx}{a numeric vector.
#' Age of patient at diagnosis.} \item{time}{a numeric vector. Survival time:
#' time to blindness or censoring.} \item{status}{ a numeric vector code.
#' Survival status. 1: blindness, 0: censored.} \item{trteye}{a numeric vector
#' code. Random eye selected for treatment. 1: left eye 2: right eye.}
#' \item{treat}{a numeric vector. 1: treatment 0: untreated.} \item{adult}{a
#' numeric vector code. 1: younger than 20, 2: older than 20.} }
#' @source Huster W.J. and Brookmeyer, R. and Self. S. (1989) MOdelling paired
#' survival data with covariates, Biometrics 45, 145-56.
#' @name  diabetes
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(diabetes)
#' names(diabetes)
#' 
NULL





#' Melanoma data and Danish population mortality by age and sex
#' 
#' Melanoma data with background mortality of Danish population.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{id}{ a numeric vector. Gives patient id.  } \item{sex}{ a numeric
#' vector. Gives sex of patient.  } \item{start}{ a numeric vector.  Gives the
#' starting time for the time-interval for which the covariate rate is
#' representative.  } \item{stop}{ a numeric vector. Gives the stopping time
#' for the time-interval for which the covariate rate is representative.  }
#' \item{status}{ a numeric vector code. Survival status. 1: dead from
#' melanoma, 0: alive or dead from other cause. } \item{age}{ a numeric vector.
#' Gives the age of the patient at removal of tumor.  } \item{rate}{ a numeric
#' vector. Gives the population mortality for the given sex and age.  Based on
#' Table A.2 in Andersen et al. (1993).  } }
#' @source Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993),
#' \emph{Statistical Models Based on Counting Processes}, Springer-Verlag.
#' @name  mela.pop 
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(mela.pop)
#' names(mela.pop)
#' 
NULL





#' The Melanoma Survival Data
#' 
#' The melanoma data frame has 205 rows and 7 columns.  It contains data
#' relating to survival of patients after operation for malignant melanoma
#' collected at Odense University Hospital by K.T.  Drzewiecki.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{no}{ a numeric vector. Patient code. } \item{status}{ a numeric vector
#' code. Survival status. 1: dead from melanoma, 2: alive, 3: dead from other
#' cause. } \item{days}{ a numeric vector. Survival time. } \item{ulc}{ a
#' numeric vector code. Ulceration, 1: present, 0: absent. } \item{thick}{ a
#' numeric vector. Tumour thickness (1/100 mm). } \item{sex}{ a numeric vector
#' code. 0: female, 1: male. } }
#' @source Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993),
#' \emph{Statistical Models Based on Counting Processes}, Springer-Verlag.
#' 
#' Drzewiecki, K.T., Ladefoged, C., and Christensen, H.E. (1980), Biopsy and
#' prognosis for cutaneous malignant melanoma in clinical stage I. Scand. J.
#' Plast. Reconstru. Surg. 14, 141-144.
#' @name  melanoma 
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(melanoma)
#' names(melanoma)
#' 
NULL





#' The TRACE study group of myocardial infarction
#' 
#' The TRACE data frame contains 1877 patients and is a subset of a data set
#' consisting of approximately 6000 patients.  It contains data relating
#' survival of patients after myocardial infarction to various risk factors.
#' 
#' sTRACE is a subsample consisting of 300 patients.
#' 
#' tTRACE is a subsample consisting of 1000 patients.
#' 
#' 
#' @aliases TRACE sTRACE tTRACE
#' @format This data frame contains the following columns: \describe{
#' \item{id}{a numeric vector. Patient code. } \item{status}{ a numeric vector
#' code. Survival status. 9: dead from myocardial infarction, 0: alive, 7: dead
#' from other causes.  } \item{time}{ a numeric vector. Survival time in years.
#' } \item{chf}{ a numeric vector code. Clinical heart pump failure, 1:
#' present, 0: absent. } \item{diabetes}{ a numeric vector code. Diabetes, 1:
#' present, 0: absent. } \item{vf}{ a numeric vector code. Ventricular
#' fibrillation, 1: present, 0: absent. } \item{wmi}{ a numeric vector.
#' Measure of heart pumping effect based on ultrasound measurements where 2 is
#' normal and 0 is worst. } \item{sex}{ a numeric vector code. 1: female, 0:
#' male. } \item{age}{ a numeric vector code. Age of patient. } }
#' @source The TRACE study group.
#' 
#' Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L., F. E. Nielsen,
#' Melchior, T., Joen, T. and P. K. Andersen (1997), Does in-hospital
#' ventricular fibrillation affect prognosis after myocardial infarction?,
#' European Heart Journal 18, 919--924.
#' @name  TRACE
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(TRACE)
#' names(TRACE)
#' 
NULL



