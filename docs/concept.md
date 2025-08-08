# Concept Note

## Objective

Develop the first effective therapy for late‑stage (symptomatic) rabies by combining in‑silico screening of BBB‑penetrant small‑molecule inhibitors targeting the rabies virus (RABV) L and P proteins with rapid in‑vitro validation.

## Rationale

Current rabies infection is nearly always fatal once neurological symptoms appear; the virus quickly reaches the central nervous system and hides behind the blood–brain barrier (BBB). Traditional antiviral agents and monoclonal antibodies struggle to cross the BBB and reach infected neurons. Small molecules capable of crossing the BBB offer a complementary approach to antibody‑based therapies.

## Key Tasks and Timeline (Year 1)

- **Target Identification:** Collect available structures of the RABV L polymerase and P phosphoprotein and build homology models if necessary.
- **Library Preparation:** Compile drug‑like compound libraries (e.g. from ChEMBL, DrugBank) and filter for BBB permeability using metrics such as CNS MPO and logBB.
- **Docking and Scoring:** Run high‑throughput docking using AutoDock Vina and machine‑learning scoring functions (DeepChem) to rank compounds.
- **Selection of Hits:** Perform clustering and prioritization to select ~50 top‑scoring compounds for further testing.
- **In Vitro Validation:** Evaluate selected compounds in neuronal cell culture infected with rabies virus to determine EC₅₀ values.
- **Patenting and Reporting:** Assess the patent landscape and prepare filings for promising chemical scaffolds.

## Repository Structure

This repository contains code and documentation to support the computational part of the project. Laboratory protocols will be kept separately.
