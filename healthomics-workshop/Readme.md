# AWS HealthOmics Workshop - Complete Training Guide

## 🧬 Welcome to AWS HealthOmics!

This comprehensive workshop will take you from beginner to advanced level in using AWS HealthOmics for genomics data storage, analysis, and querying.

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Workshop Structure](#workshop-structure)
4. [Getting Started](#getting-started)
5. [Module Descriptions](#module-descriptions)
6. [Quick Start Guide](#quick-start-guide)
7. [Code Explanations](#code-explanations)
8. [Advanced Topics](#advanced-topics)
9. [Troubleshooting](#troubleshooting)
10. [Resources](#resources)

---

## 🎯 Overview

### What You'll Learn

By completing this workshop, you will:

1. ✅ **Understand AWS HealthOmics architecture** - Learn how the service works
2. ✅ **Create and manage stores** - Set up storage for genomic data
3. ✅ **Upload and process data** - Load sequencing data and references
4. ✅ **Run bioinformatics workflows** - Execute analysis pipelines at scale
5. ✅ **Work with variants** - Load and manage variant calls (VCFs)
6. ✅ **Query data with SQL** - Use Athena for powerful variant queries

### Time Commitment

- **Quick Start**: 30 minutes (understanding concepts only)
- **Core Workshop**: 1-2 hours (hands-on with provided data)
- **Advanced**: 3-4 hours (complete implementation with your data)

---

## 🔧 Prerequisites

### Required

- **AWS Account** with HealthOmics access
- **AWS CLI** installed and configured
- **Basic genomics knowledge**:
  - FASTQ files (sequencing reads)
  - FASTA files (reference genomes)
  - VCF files (variant calls)
  - BAM/CRAM files (aligned reads)

### Optional but Recommended

- Python 3.8+ with boto3
- Understanding of bioinformatics workflows
- SQL knowledge for querying
- S3 bucket for data storage

### AWS Regions

HealthOmics is available in:
- us-east-1 (N. Virginia)
- us-west-2 (Oregon)
- eu-west-1 (Ireland)
- eu-west-2 (London)
- eu-central-1 (Frankfurt)
- ap-southeast-1 (Singapore)

---

## 📚 Workshop Structure

### File Organization

```
healthomics-workshop/
├── README.md                     # This file
├── workshop_guide.md             # Detailed conceptual guide
├── healthomics_setup.sh          # Initial setup script
│
├── Module 1: Store Creation
│   ├── 01_create_stores.sh       # Bash version
│   └── 01_create_stores.py       # Python version (recommended)
│
├── Module 2: Data Upload
│   └── 02_upload_data.py         # Upload sequencing data and references
│
├── Module 3: Workflows
│   └── 03_workflow_management.py # Create and run workflows
│
├── Module 4: Variants
│   └── 04_variant_management.py  # Load and manage VCF data
│
└── Module 5: Queries
    └── 05_athena_queries.py      # Query variants with Athena
```

---

## 🚀 Getting Started

### Step 1: Initial Setup

```bash
# Make setup script executable
chmod +x healthomics_setup.sh

# Run setup
./healthomics_setup.sh

# Load environment variables
source ~/healthomics-workshop/env_vars.sh
```

### Step 2: Verify AWS Configuration

```bash
# Check AWS credentials
aws sts get-caller-identity

# Verify HealthOmics access
aws omics list-sequence-stores --region us-east-1
```

### Step 3: Start with Module 1

```bash
# Make Python scripts executable
chmod +x *.py

# Run Module 1
python3 01_create_stores.py
```

---

## 📖 Module Descriptions

### Module 1: Create HealthOmics Stores

**Purpose**: Set up the foundation for storing genomic data

**What You'll Create**:
- Sequence Store (for FASTQ/BAM/CRAM files)
- Reference Store (for reference genomes)
- Annotation Store (for variant annotations)

**Files**: `01_create_stores.py` or `01_create_stores.sh`

**Duration**: 5-10 minutes

**Key Concepts**:
```python
# Creating a Sequence Store
response = omics_client.create_sequence_store(
    name="my-sequence-store",
    description="Store for raw sequencing data"
)

# Why? Sequence stores:
# - Automatically compress data (50-70% savings)
# - Encrypt data at rest
# - Provide fast, efficient access
```

---

### Module 2: Upload Data

**Purpose**: Load sequencing data and reference genomes

**What You'll Upload**:
- Reference genome (FASTA format)
- Sequencing reads (FASTQ files)

**Files**: `02_upload_data.py`

**Duration**: 15-30 minutes (depends on file size)

**Key Concepts**:
```python
# Uploading a reference genome
response = omics_client.start_reference_import_job(
    referenceStoreId=store_id,
    sources=[{
        'sourceFile': 's3://my-bucket/references/hg38.fasta',
        'name': 'GRCh38'
    }]
)

# Why upload to HealthOmics?
# - Optimized storage format
# - Built-in compression
# - Fast retrieval for workflows
# - Automatic indexing
```

---

### Module 3: Run Workflows

**Purpose**: Execute bioinformatics analysis pipelines

**What You'll Learn**:
- Creating workflow definitions
- Configuring workflow parameters
- Starting workflow runs
- Monitoring progress

**Files**: `03_workflow_management.py`

**Duration**: 30-60 minutes

**Key Concepts**:
```python
# Starting a workflow run
response = omics_client.start_run(
    workflowId=workflow_id,
    parameters={
        'reference': 'GRCh38',
        'sample_fastq_1': 's3://bucket/sample_R1.fastq.gz',
        'sample_fastq_2': 's3://bucket/sample_R2.fastq.gz'
    },
    outputUri='s3://bucket/outputs/'
)

# Common workflows:
# - Alignment (BWA-MEM, Bowtie2)
# - Variant calling (GATK, FreeBayes)
# - RNA-Seq (STAR, Salmon)
# - Quality control (FastQC)
```

---

### Module 4: Manage Variants

**Purpose**: Load and organize variant data

**What You'll Do**:
- Create Variant Store
- Import VCF files
- Load annotations
- Manage variant metadata

**Files**: `04_variant_management.py`

**Duration**: 15-30 minutes

**Key Concepts**:
```python
# Importing VCF to Variant Store
response = omics_client.start_variant_import_job(
    destinationName=variant_store_name,
    items=[{
        'source': 's3://bucket/variants/sample.vcf.gz'
    }],
    runLeftNormalization=True  # Standardize variant representation
)

# Benefits:
# - Queryable variant database
# - Join with annotations
# - Fast filtering
# - Cohort analysis
```

---

### Module 5: Query with Athena

**Purpose**: Analyze variants using SQL

**What You'll Query**:
- Variants by gene
- Pathogenic variants
- Rare variants
- Genomic regions

**Files**: `05_athena_queries.py`

**Duration**: 30-45 minutes

**Key Concepts**:
```sql
-- Find pathogenic variants in BRCA1
SELECT 
    v.contigname as chromosome,
    v.start as position,
    v.referenceallele as ref,
    v.alternatealleles as alt,
    a.clinical_significance
FROM variant_store v
JOIN annotation_store a 
ON v.contigname = a.chr AND v.start = a.pos
WHERE a.gene = 'BRCA1' 
AND a.clinical_significance = 'Pathogenic';

-- Why Athena?
-- - Standard SQL interface
-- - Serverless (no infrastructure)
-- - Integrated with HealthOmics
-- - Cost-effective for large datasets
```

---

## 🎓 Code Explanations

### Understanding the Store Creation Process

#### What is a Store?

A **store** in HealthOmics is a specialized database optimized for genomic data types:

1. **Sequence Store**: 
   - Stores raw sequencing files (FASTQ, BAM, CRAM)
   - Automatically applies compression algorithms
   - Typically achieves 50-70% size reduction
   - Example: 100GB FASTQ → 30-50GB in Sequence Store

2. **Reference Store**:
   - Stores reference genomes (FASTA files)
   - Used by workflows and variant calling
   - Examples: GRCh38, GRCh37, mouse genome (mm10)

3. **Variant Store**:
   - Stores variant calls (VCF files)
   - Creates queryable database
   - Enables SQL-based analysis
   - Supports joins with annotations

4. **Annotation Store**:
   - Stores variant annotations
   - Sources: ClinVar, dbSNP, gnomAD, etc.
   - Enriches variants with additional information

#### Code Walkthrough: Creating a Sequence Store

```python
def create_sequence_store(self):
    """
    Create a Sequence Store for raw sequencing data
    """
    # Step 1: Call AWS HealthOmics API
    response = self.omics_client.create_sequence_store(
        name=f"{self.workshop_name}-sequence-store",  # Unique name
        description="Workshop sequence store",         # Optional description
        tags={                                         # Tags for organization
            'Project': self.workshop_name,
            'Purpose': 'Workshop'
        }
    )
    
    # Step 2: Extract the store ID from response
    store_id = response['id']
    # Example ID: "1234567890"
    
    # Step 3: Wait for store to become active
    # Stores go through: CREATING → ACTIVE
    self._wait_for_sequence_store(store_id)
    
    # Step 4: Store is now ready to accept data
    return store_id
```

**Why these steps?**
- Creating a store is asynchronous (takes a few seconds)
- We need to wait for `ACTIVE` status before using it
- The store ID is used for all subsequent operations

---

### Understanding Data Upload

#### How Data Upload Works

```python
def upload_sequencing_data(self, s3_uri_r1, subject_id, sample_id):
    """
    Upload FASTQ files to Sequence Store
    """
    # Step 1: Start import job
    response = self.omics_client.start_read_set_import_job(
        sequenceStoreId=self.store_id,  # Where to store
        roleArn=self.role_arn,           # Permissions
        sources=[{
            'subjectId': subject_id,      # Patient/sample identifier
            'sampleId': sample_id,
            'sourceFiles': {
                'source1': s3_uri_r1      # S3 location of FASTQ
            }
        }]
    )
    
    # Step 2: HealthOmics processes the file
    # - Downloads from S3
    # - Validates format
    # - Applies compression
    # - Stores in optimized format
    
    # Step 3: Monitor progress
    job_id = response['id']
    self._wait_for_import_completion(job_id)
    
    # Step 4: Get read set ID
    read_set_id = self._get_read_set_id(job_id)
    return read_set_id
```

**What happens during import?**

1. **Validation**: Checks file format (valid FASTQ/BAM/CRAM)
2. **Compression**: Applies specialized genomic compression
3. **Indexing**: Creates internal indexes for fast access
4. **Storage**: Saves in HealthOmics optimized format

**Example compression:**
```
Original FASTQ: 10 GB
↓ HealthOmics compression
Stored size: 3-5 GB (50-70% savings)
```

---

### Understanding Workflows

#### Workflow Lifecycle

```python
def run_workflow(self, workflow_id, parameters):
    """
    Execute a bioinformatics workflow
    """
    # Step 1: Start workflow run
    response = self.omics_client.start_run(
        workflowId=workflow_id,
        parameters=parameters,      # Input files, settings
        outputUri='s3://outputs/',  # Where to save results
        roleArn=self.role_arn
    )
    
    run_id = response['id']
    
    # Step 2: Workflow engine processes tasks
    # - Provisions compute resources
    # - Downloads input files
    # - Runs each workflow step
    # - Uploads outputs to S3
    
    # Step 3: Monitor execution
    while True:
        status = self._get_run_status(run_id)
        if status == 'COMPLETED':
            break
        elif status == 'FAILED':
            raise Exception("Workflow failed")
        time.sleep(30)
    
    return run_id
```

**Workflow execution flow:**

```
User starts run
      ↓
Engine allocates resources
      ↓
Downloads inputs from S3/Stores
      ↓
Executes workflow steps:
  - Task 1: Alignment
  - Task 2: Sort BAM
  - Task 3: Mark duplicates
  - Task 4: Call variants
      ↓
Uploads outputs to S3
      ↓
Workflow completes
```

---

### Understanding Variant Queries

#### How Athena Queries Work

When you import a VCF into a Variant Store, HealthOmics:

1. **Parses the VCF file**
2. **Creates an Athena table** automatically
3. **Stores data** in queryable format
4. **Enables SQL queries** on the data

```python
def query_pathogenic_variants(self):
    """
    Find clinically significant variants
    """
    query = """
    SELECT 
        v.contigname as chromosome,
        v.start as position,
        v.referenceallele as ref,
        v.alternatealleles as alt,
        a.clinical_significance
    FROM variant_store v
    JOIN annotation_store a 
    ON v.contigname = a.chr 
    AND v.start = a.pos
    WHERE a.clinical_significance = 'Pathogenic'
    LIMIT 100;
    """
    
    # Execute query
    results = self.athena_client.start_query_execution(
        QueryString=query,
        ResultConfiguration={
            'OutputLocation': 's3://results/'
        }
    )
    
    return results
```

**Query execution:**

```
1. User submits SQL query
      ↓
2. Athena parses query
      ↓
3. Reads data from Variant/Annotation Stores
      ↓
4. Filters and joins data
      ↓
5. Returns results
      ↓
6. User analyzes results
```

---

## 🔬 Advanced Topics

### 1. Workflow Optimization

#### Resource Allocation

```python
# Right-size your workflow for cost efficiency
workflow_params = {
    'cpu': 16,           # CPUs per task
    'memory': 32768,     # Memory in MB
    'storage': 100       # Storage in GB
}

# Tips:
# - Start small, scale up if needed
# - Monitor resource usage in CloudWatch
# - Use spot instances for cost savings (automatic)
```

#### Parallel Processing

```python
# Process multiple samples simultaneously
run_group = omics_client.create_run_group(
    name='batch-processing',
    maxCpus=100000,      # Total CPU limit
    maxRuns=100          # Concurrent runs
)

# Submit multiple runs to the group
for sample in samples:
    omics_client.start_run(
        workflowId=workflow_id,
        runGroupId=run_group_id,
        parameters=sample_params
    )
```

---

### 2. Cost Optimization Strategies

#### Storage Costs

```python
# Sequence Store: $0.01/GB/month (after compression)
# Example: 1TB of FASTQ → 300-500GB stored
# Monthly cost: $3-$5

# Optimize:
# 1. Delete intermediate files
# 2. Use lifecycle policies on S3
# 3. Archive old data to Glacier
```

#### Compute Costs

```python
# Workflow costs based on:
# - vCPU-hours used
# - Memory-hours used
# - Storage accessed

# Optimization strategies:
# 1. Right-size resources
# 2. Use efficient algorithms
# 3. Cache reference data
# 4. Monitor and optimize slow steps
```

---

### 3. Integration Patterns

#### With AWS Step Functions

```python
# Orchestrate complex pipelines
import boto3

stepfunctions = boto3.client('stepfunctions')

# Define state machine
state_machine = {
    'StartAt': 'ImportData',
    'States': {
        'ImportData': {
            'Type': 'Task',
            'Resource': 'arn:aws:lambda:...',
            'Next': 'RunWorkflow'
        },
        'RunWorkflow': {
            'Type': 'Task',
            'Resource': 'arn:aws:omics:...',
            'Next': 'LoadVariants'
        },
        'LoadVariants': {
            'Type': 'Task',
            'Resource': 'arn:aws:lambda:...',
            'End': True
        }
    }
}
```

#### With Lambda for Automation

```python
# Trigger workflows on S3 uploads
def lambda_handler(event, context):
    # Get S3 event
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = event['Records'][0]['s3']['object']['key']
    
    # Start HealthOmics workflow
    omics = boto3.client('omics')
    omics.start_run(
        workflowId='workflow-123',
        parameters={
            'input_fastq': f's3://{bucket}/{key}'
        }
    )
```

---

## 🔍 Troubleshooting

### Common Issues and Solutions

#### 1. "Access Denied" Errors

**Problem**: Cannot create stores or import data

**Solution**:
```bash
# Check IAM permissions
aws iam get-user-policy --user-name your-user --policy-name HealthOmicsAccess

# Required permissions:
# - omics:CreateSequenceStore
# - omics:StartReadSetImportJob
# - s3:GetObject (for import source)
# - s3:PutObject (for workflow outputs)
```

#### 2. Import Job Failures

**Problem**: Data import fails or hangs

**Common causes**:
- Invalid file format
- Corrupted files
- S3 permissions
- Network issues

**Solution**:
```python
# Check file integrity
import gzip
with gzip.open('sample.fastq.gz', 'rt') as f:
    # Try reading first few lines
    for i, line in enumerate(f):
        if i >= 4:
            break
        print(line)

# Verify S3 permissions
aws s3 ls s3://your-bucket/your-file.fastq.gz
```

#### 3. Workflow Failures

**Problem**: Workflow runs fail

**Debugging steps**:
```python
# Get run details
response = omics_client.get_run(id=run_id)

# Check task logs
for task in response['tasks']:
    if task['status'] == 'FAILED':
        print(f"Failed task: {task['name']}")
        print(f"Error: {task.get('statusMessage')}")

# View CloudWatch logs
logs_client = boto3.client('logs')
logs_client.get_log_events(
    logGroupName=f"/aws/omics/WorkflowLog",
    logStreamName=run_id
)
```

#### 4. Athena Query Errors

**Problem**: Queries fail or return no results

**Solutions**:
```sql
-- Check if table exists
SHOW TABLES IN healthomics_database;

-- Verify column names
DESCRIBE variant_store;

-- Test simple query first
SELECT COUNT(*) FROM variant_store;

-- Common fixes:
-- 1. Check table permissions
-- 2. Verify partition metadata
-- 3. Ensure correct data types in WHERE clauses
```

---

## 📊 Performance Optimization

### Query Performance

```sql
-- Bad: Full table scan
SELECT * FROM variant_store 
WHERE gene = 'BRCA1';

-- Good: Use partitioning
SELECT * FROM variant_store 
WHERE contigname = 'chr17'  -- Partition key
AND start BETWEEN 43044295 AND 43170245  -- Index
AND gene = 'BRCA1';

-- Best: Specific filters, limited results
SELECT contigname, start, referenceallele, alternatealleles
FROM variant_store 
WHERE contigname = 'chr17'
AND start BETWEEN 43044295 AND 43170245
LIMIT 1000;
```

### Workflow Performance

```python
# Enable caching for repeated runs
workflow_params = {
    'enable_caching': True,  # Reuse intermediate results
    'cache_location': 's3://bucket/cache/'
}

# Use appropriate instance types
# - Memory-intensive: r5 instances
# - Compute-intensive: c5 instances
# - General purpose: m5 instances
```

---

## 📚 Additional Resources

### Official Documentation
- [AWS HealthOmics Docs](https://docs.aws.amazon.com/omics/)
- [HealthOmics API Reference](https://docs.aws.amazon.com/omics/latest/api/)
- [Best Practices Guide](https://docs.aws.amazon.com/omics/latest/dev/best-practices.html)

### Sample Workflows
- [GitHub: AWS HealthOmics Examples](https://github.com/aws-samples/amazon-omics-tutorials)
- [Nextflow Workflows](https://nf-co.re/)
- [WDL Workflows](https://github.com/gatk-workflows/)

### Community
- [AWS re:Post](https://repost.aws/)
- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Nextflow Community](https://nextflow.io/community.html)

---

## ✅ Workshop Checklist

Use this to track your progress:

- [ ] Completed environment setup
- [ ] Created all four store types
- [ ] Uploaded a reference genome
- [ ] Uploaded sequencing data
- [ ] Created/imported a workflow
- [ ] Successfully ran a workflow
- [ ] Imported VCF to Variant Store
- [ ] Loaded annotations
- [ ] Executed Athena queries
- [ ] Joined variants with annotations
- [ ] Filtered by clinical significance
- [ ] Exported query results

---

## 🎉 Congratulations!

You've completed the AWS HealthOmics workshop! You now know how to:

✅ Store and manage genomic data at scale
✅ Run bioinformatics workflows efficiently
✅ Query variants using SQL
✅ Integrate HealthOmics into your analysis pipelines

### Next Steps

1. **Apply to Your Data**: Use these scripts with your own genomic datasets
2. **Customize Workflows**: Adapt existing workflows for your use case
3. **Automate**: Build automated pipelines with Step Functions
4. **Scale**: Process cohorts and population-scale data
5. **Share**: Help others by contributing workflows and best practices

---

## 📧 Support

For issues or questions:
- AWS Support: Open a case in AWS Console
- Workshop issues: Check the troubleshooting section
- Feature requests: Contact your AWS account team

---

**Last Updated**: October 2025
**Version**: 1.0
**Maintained by**: AWS HealthOmics Workshop Team