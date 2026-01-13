# AWS HealthOmics Quick Reference Guide

## 🚀 Quick Start Commands

### Initial Setup
```bash
# Configure AWS CLI
aws configure

# Check HealthOmics access
aws omics list-sequence-stores --region us-east-1

# Run workshop setup
./healthomics_setup.sh
source ~/healthomics-workshop/env_vars.sh
```

---

## 📦 Store Management

### Create Stores
```python
import boto3
omics = boto3.client('omics', region_name='us-east-1')

# Sequence Store
response = omics.create_sequence_store(name='my-seq-store')
seq_store_id = response['id']

# Reference Store
response = omics.create_reference_store(name='my-ref-store')
ref_store_id = response['id']

# Variant Store (requires reference)
response = omics.create_variant_store(
    name='my-var-store',
    reference={'referenceArn': reference_arn}
)

# Annotation Store
response = omics.create_annotation_store(
    name='my-anno-store',
    storeFormat='TSV'
)
```

### List Stores
```bash
# List all sequence stores
aws omics list-sequence-stores --region us-east-1

# List all reference stores
aws omics list-reference-stores --region us-east-1

# List all variant stores
aws omics list-variant-stores --region us-east-1

# List all annotation stores
aws omics list-annotation-stores --region us-east-1
```

---

## 📤 Data Upload

### Upload Reference Genome
```python
response = omics.start_reference_import_job(
    referenceStoreId='ref-store-id',
    roleArn='arn:aws:iam::123456789012:role/HealthOmicsRole',
    sources=[{
        'sourceFile': 's3://bucket/reference.fasta',
        'name': 'GRCh38'
    }]
)
```

### Upload Sequencing Data
```python
response = omics.start_read_set_import_job(
    sequenceStoreId='seq-store-id',
    roleArn='arn:aws:iam::123456789012:role/HealthOmicsRole',
    sources=[{
        'subjectId': 'PATIENT001',
        'sampleId': 'SAMPLE001',
        'sourceFiles': {
            'source1': 's3://bucket/sample_R1.fastq.gz',
            'source2': 's3://bucket/sample_R2.fastq.gz'
        }
    }]
)
```

### Check Import Status
```bash
# Check reference import
aws omics get-reference-import-job \
  --id job-id \
  --reference-store-id ref-store-id \
  --region us-east-1

# Check read set import
aws omics get-read-set-import-job \
  --id job-id \
  --sequence-store-id seq-store-id \
  --region us-east-1
```

---

## 🔄 Workflow Management

### Start a Workflow Run
```python
response = omics.start_run(
    workflowId='workflow-id',
    roleArn='arn:aws:iam::123456789012:role/HealthOmicsRole',
    name='my-workflow-run',
    parameters={
        'reference': 'GRCh38',
        'input_fastq_1': 's3://bucket/sample_R1.fastq.gz',
        'input_fastq_2': 's3://bucket/sample_R2.fastq.gz'
    },
    outputUri='s3://bucket/outputs/'
)
run_id = response['id']
```

### Monitor Workflow Run
```bash
# Get run status
aws omics get-run --id run-id --region us-east-1

# List all runs
aws omics list-runs --region us-east-1

# List runs by status
aws omics list-runs --status RUNNING --region us-east-1
```

### Cancel a Run
```bash
aws omics cancel-run --id run-id --region us-east-1
```

---

## 🧬 Variant Management

### Create Variant Store with Reference
```python
response = omics.create_variant_store(
    name='my-variant-store',
    reference={
        'referenceArn': 'arn:aws:omics:us-east-1:123456789012:referenceStore/ref-store/reference/ref-id'
    }
)
```

### Import VCF to Variant Store
```python
response = omics.start_variant_import_job(
    destinationName='my-variant-store',
    roleArn='arn:aws:iam::123456789012:role/HealthOmicsRole',
    items=[{
        'source': 's3://bucket/variants.vcf.gz'
    }],
    runLeftNormalization=True
)
```

### Import Annotations
```python
response = omics.start_annotation_import_job(
    destinationName='my-annotation-store',
    roleArn='arn:aws:iam::123456789012:role/HealthOmicsRole',
    items=[{
        'source': 's3://bucket/annotations.tsv'
    }]
)
```

---

## 🔍 Athena Queries

### Connect to Athena
```python
athena = boto3.client('athena', region_name='us-east-1')

def run_query(query):
    response = athena.start_query_execution(
        QueryString=query,
        ResultConfiguration={
            'OutputLocation': 's3://bucket/athena-results/'
        }
    )
    return response['QueryExecutionId']
```

### Common Queries

#### Find Variants by Gene
```sql
SELECT contigname, start, referenceallele, alternatealleles
FROM "healthomics"."variant_store"
WHERE gene = 'BRCA1'
LIMIT 100;
```

#### Pathogenic Variants
```sql
SELECT v.contigname, v.start, v.referenceallele, v.alternatealleles,
       a.gene, a.clinical_significance
FROM "healthomics"."variant_store" v
JOIN "healthomics"."annotation_store" a
  ON v.contigname = a.chr AND v.start = a.pos
WHERE a.clinical_significance IN ('Pathogenic', 'Likely_pathogenic');
```

#### Rare Variants
```sql
SELECT v.contigname, v.start, a.gene, a.population_af
FROM "healthomics"."variant_store" v
JOIN "healthomics"."annotation_store" a
  ON v.contigname = a.chr AND v.start = a.pos
WHERE CAST(a.population_af AS DOUBLE) < 0.01
ORDER BY CAST(a.population_af AS DOUBLE);
```

#### Variants in Region
```sql
SELECT contigname, start, referenceallele, alternatealleles, qual
FROM "healthomics"."variant_store"
WHERE contigname = 'chr17'
  AND start BETWEEN 43044295 AND 43170245
ORDER BY start;
```

#### Variant Count by Chromosome
```sql
SELECT contigname, COUNT(*) as count
FROM "healthomics"."variant_store"
GROUP BY contigname
ORDER BY contigname;
```

---

## 🛠️ Common Operations

### Delete Resources

#### Delete Read Set
```bash
aws omics delete-read-set \
  --id read-set-id \
  --sequence-store-id seq-store-id \
  --region us-east-1
```

#### Delete Reference
```bash
aws omics delete-reference \
  --id reference-id \
  --reference-store-id ref-store-id \
  --region us-east-1
```

#### Delete Store
```bash
# Delete sequence store
aws omics delete-sequence-store \
  --id seq-store-id \
  --region us-east-1

# Delete reference store
aws omics delete-reference-store \
  --id ref-store-id \
  --region us-east-1
```

### Tag Resources
```python
omics.tag_resource(
    resourceArn='arn:aws:omics:us-east-1:123456789012:sequenceStore/store-id',
    tags={
        'Project': 'Cancer-Genomics',
        'Environment': 'Production',
        'CostCenter': 'Research'
    }
)
```

---

## 📊 Monitoring & Troubleshooting

### Check Store Status
```bash
# Get sequence store details
aws omics get-sequence-store --id store-id --region us-east-1

# Get reference store details
aws omics get-reference-store --id store-id --region us-east-1

# Get variant store details
aws omics get-variant-store --name store-name --region us-east-1
```

### View CloudWatch Logs
```bash
# List log groups
aws logs describe-log-groups \
  --log-group-name-prefix /aws/omics \
  --region us-east-1

# Get log events
aws logs get-log-events \
  --log-group-name /aws/omics/WorkflowLog \
  --log-stream-name run-id \
  --region us-east-1
```

### Get Run Details
```python
response = omics.get_run(id='run-id')

print(f"Status: {response['status']}")
print(f"Started: {response.get('startTime')}")
print(f"Stopped: {response.get('stopTime')}")

# Check task statuses
for task in response.get('tasks', []):
    print(f"{task['name']}: {task['status']}")
```

---

## 💰 Cost Optimization

### Storage Pricing
```
Sequence Store: $0.01/GB/month (after compression)
Reference Store: $0.02/GB/month
Variant Store: $0.40/GB/month
Annotation Store: $0.40/GB/month
```

### Compute Pricing
```
Workflows: Based on vCPU-hours and GB-hours
Check current pricing: https://aws.amazon.com/omics/pricing/
```

### Best Practices
- Delete intermediate files after processing
- Use S3 lifecycle policies
- Right-size workflow resources
- Enable workflow caching
- Monitor with Cost Explorer

---

## 🔐 Security Best Practices

### IAM Policy Example
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "omics:CreateSequenceStore",
        "omics:GetSequenceStore",
        "omics:ListSequenceStores",
        "omics:StartReadSetImportJob",
        "omics:GetReadSet"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject"
      ],
      "Resource": "arn:aws:s3:::your-bucket/*"
    }
  ]
}
```

### Encryption
- All data encrypted at rest by default
- Uses AWS KMS for key management
- Option to use customer-managed keys (CMK)

---

## 📞 Support & Resources

### Documentation
- Main docs: https://docs.aws.amazon.com/omics/
- API Reference: https://docs.aws.amazon.com/omics/latest/api/
- Best Practices: https://docs.aws.amazon.com/omics/latest/dev/best-practices.html

### Sample Code
- GitHub: https://github.com/aws-samples/amazon-omics-tutorials
- Workflow Examples: https://github.com/aws-samples/amazon-omics-workflows

### Getting Help
- AWS Support Console
- AWS re:Post: https://repost.aws/
- Service Quotas: Check and request increases in AWS Console

---

## ⚡ Performance Tips

1. **Use partitioning** in Athena queries (by chromosome)
2. **Enable caching** for repeated workflow runs
3. **Batch operations** when processing multiple samples
4. **Monitor CloudWatch** metrics for bottlenecks
5. **Optimize queries** with EXPLAIN and proper indexing
6. **Use run groups** to organize and limit concurrent workflows

---

## 🔄 Common Workflows

### End-to-End Variant Calling Pipeline
```
1. Upload FASTQ → Sequence Store
2. Upload Reference → Reference Store
3. Run BWA-MEM alignment workflow
4. Run GATK variant calling workflow
5. Import VCF → Variant Store
6. Load annotations → Annotation Store
7. Query with Athena for analysis
```

### RNA-Seq Analysis Pipeline
```
1. Upload FASTQ → Sequence Store
2. Upload transcriptome reference → Reference Store
3. Run STAR alignment workflow
4. Run Salmon quantification workflow
5. Export results for downstream analysis
```

---

## 📝 Quick Troubleshooting

| Issue | Solution |
|-------|----------|
| Access Denied | Check IAM permissions |
| Import Failed | Verify file format and S3 permissions |
| Workflow Timeout | Increase resource allocation |
| Query Too Expensive | Add partitioning filters |
| No Results in Query | Check table names and schemas |

---

**Version**: 1.0
**Last Updated**: October 2025