#Visit data

fit <- glm(formula= ARW_NetComplexity ~   factor(Diagnosis_group) + factor(Clinic)
          + Action_count + Role_count + Workstation_count, data=visits)



# Fixed Effects model for Clinic Days
library(plm)
fixed <- plm(ARW_NetComplexity ~ Action_count + Role_count + Workstation_count + NumUniqueProcedures+
               NumUniqueDiagnosisGroups + NumPhysicians + NumVisits,
             data=clinic_days, index = c("Clinic", "YMD_date"), model="within")

summary(fixed)
fixef(fixed)
pFtest(fixed, fit)
