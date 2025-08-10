# Workflow

```mermaid
graph LR
  A[Prepare inputs] --> B[Docking]
  B --> C[Rescore]
  C --> D[Report]
  subgraph Outputs
    D --> E(results.parquet)
    D --> F(report.html)
  end
```
