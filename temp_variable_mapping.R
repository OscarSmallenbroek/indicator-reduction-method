# Create a mapping of abbreviated variable codes to full names
variable_mapping <- data.frame(
  code = c("GovEff", "ResPPop", "LogPerf", "FinanceSSGEMplus", "ResTalBus", "EcComplexity", "EntMedia", "RegQua", "RuleOL", "PISA", "RDExp", "KnowIntsEmp", "GERDFinBus", "PublicPrivCopubs", "THEunirank"),
  full_name = c("Government Effectiveness", "Researchers, Full-Time Equivalent per 1,000 Labour Force", "Logistics Performance Index", "Financial Support to SMEs", "Researchers, Tertiary Sector", "Economic Complexity", "Entertainment Media", "Regulatory Quality", "Rule of Law", "PISA Score", "Research and Development Expenditure as Percentage of GDP", "Knowledge-intensive Employment", "GERD Financing from Business", "Public-Private Publication Co-authorships", "THE University Ranking Score"),
  domain = c("Wellbeing Today", "Resources for Future", "Connectivity", "Resources for Future", "Societal Resilience", "Societal Resilience", "Societal Resilience", "Institutions", "Institutions", "Resources for Future", "Resources for Future", "Societal Resilience", "Resources for Future", "Societal Resilience", "Resources for Future")
)

# Save mapping for reference
write.csv(variable_mapping, "data/variable_mapping.csv", row.names = FALSE)