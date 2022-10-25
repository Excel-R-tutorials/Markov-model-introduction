
#
input_tab <- 
  tibble::tribble(
    ~Name, ~Value, ~Description,
    "cAsymp", cAsymp, "Cost of one cycle in the asymptomatic disease state",
    "cDeath", cDeath, "Cost associated with transition to the dead state",
    "cDrug", cDrug, "Cost of drug for one cycle",
    "tpDcm", tpDcm, "Probability of dying from the disease in a single cycle",
    "cProg", cProg, "Cost of one cycle in the progressive disease state",
    "effect", effect, "Effectiveness of drug in terms of reducing disease progression",
    "tpProg", 0.01, "Coefficent of increase for probability of entering the progressive disease state",
    "uAsymp", uAsymp, "Quality of life weight for one cycle in the asymptomatic disease state",
    "uProg", uProg, "Quality of life weight for one cycle in the progressive disease state",
    "oDr", oDr, "Discount rate for outcomes",
    "cDr", cDr, "Discount rate for costs",
    "cycle", cycle,	"Length in years of one cycle",
    "ini_age", Initial_age, "The initial age at which patients are deeemed to start the model",
    "nD35", 0.0017, "Natural death risk for over 35's (from standard life-tables)",
    "nD45", 0.0044, "Natural death risk for over 45's (from standard life-tables)",
    "nD55",	0.0138, "Natural death risk for over 55's (from standard life-tables)",
    "nD65",	0.0379, "Natural death risk for over 65's (from standard life-tables)",
    "nD75", 0.0912, "Natural death risk for over 75's (from standard life-tables)",
    "nD85",	0.1958, "Natural death risk for over 85's (from standard life-tables)"
  )

write.csv(input_tab, file = "input_tab.csv")


#
output_tab <-
  rbind(
    c("Strategy", "Eff", paste0("Cost (", enc2utf8("\u00A3"),")"), paste0(enc2utf8("\u0394"), "Eff"),
      paste0(enc2utf8("\u0394"), "Cost (", enc2utf8("\u00A3"), ")"),  "ICER"),
    c("No drug", 7.76, format(9265, big.mark=","), "" , "", ""),
    c("Drug", 8.62, format(16155, big.mark=","), 0.87, 6891, format(7931, big.mark=",")))
  
write.csv(output_tab, file = "output_tab.csv")
