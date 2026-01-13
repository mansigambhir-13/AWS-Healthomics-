# AWS HealthOmics Workflows - Architecture & Implementation Guide

## 📊 Understanding the Diagram

Your uploaded diagram shows the complete HealthOmics Workflows lifecycle in 3 steps:

---

## STEP 1: Container Setup 🐳

**What Happens:**
```
Tooling Container → Amazon ECR Private Registry
```

**Explanation:**
- Bioinformatics tools (BWA, GATK, STAR, etc.) are packaged in Docker containers
- Containers are pushed to Amazon ECR (Elastic Container Registry)
- Workflows reference these containers to run analyses

**Why Containers?**
- ✅ Reproducibility: Same environment every time
- ✅ Portability: Works anywhere
- ✅ Version control: Pin specific tool versions
- ✅ Isolation: No dependency conflicts

**Example Containers:**
- `quay.io/biocontainers/bwa:0.7.17`
- `broadinstitute/gatk:4.4.0.0`
- `biocontainers/samtools:1.15`

---

## STEP 2: Workflow Definition 📝

**What Happens:**
```
WDL/Nextflow/CWL Source Code
    ↓
Zip Bundle + parameters.json
    ↓
CreateWorkflow API
    ↓
Workflow ID
```

**Explanation:**

### 2.1 Workflow Languages

**WDL (Workflow Definition Language)**
- Developed by Broad Institute
- Used for GATK pipelines
- Human-readable syntax
- Strong type system

**Nextflow**
- Popular for complex pipelines
- Groovy-based DSL
- Excellent parallelization
- Active community

**CWL (Common Workflow Language)**
- Standard format
- Tool-agnostic
- YAML/JSON based
- Portable across platforms

### 2.2 Workflow Structure

A workflow contains:
1. **Tasks/Processes**: Individual analysis steps
2. **Inputs**: Data and parameters needed
3. **Outputs**: Results produced
4. **Dependencies**: Task execution order

### 2.3 The Zip Bundle

Contains:
- Main workflow file (e.g., `main.wdl`, `main.nf`)
- Sub-workflows (if any)
- Task definitions
- Configuration files

### 2.4 parameters.json

Defines workflow parameters:
```json
{
  "reference_genome": {
    "description": "Reference genome for alignment",
    "optional": false
  },
  "threads": {
    "description": "Number of CPU threads",
    "optional": true
  }
}
```

---

## STEP 3: Workflow Execution 🚀

**What Happens:**
```
Workflow ID + Execution Role + inputs.json + Input Data
    ↓
Start Run
    ↓
Run Group (optional)
    ↓
Workflow Run
    ↓
Run Task 1 → Run Task 2 → ... → Run Task N
    ↓
Intermediate Results (shared file system)
    ↓
Final Results (outputUri)
```

**Explanation:**

### 3.1 Required Components

**Workflow ID**: 
- Created in Step 2
- Unique identifier for your workflow definition

**Execution Role**:
- IAM role with permissions
- Read from S3 (inputs)
- Write to S3 (outputs)
- Access HealthOmics stores
- Pull from ECR

**inputs.json**:
```json
{
  "sample_fastq_1": "s3://bucket/sample_R1.fastq.gz",
  "sample_fastq_2": "s3://bucket/sample_R2.fastq.gz",
  "reference": "GRCh38",
  "threads": 8
}
```

**Input Data**:
- Raw files in S3
- Data from HealthOmics stores
- Reference genomes

### 3.2 Workflow Execution Flow

```
1. Start Run
   ↓
2. HealthOmics provisions compute resources
   ↓
3. Downloads input data
   ↓
4. Runs Task 1 (e.g., Alignment)
   ↓
5. Stores intermediate results in shared file system
   ↓
6. Runs Task 2 (e.g., Sort BAM)
   ↓
7. More intermediate results
   ↓
8. Runs Task N (e.g., Call Variants)
   ↓
9. Uploads final results to outputUri
   ↓
10. Workflow completes
```

### 3.3 Run Groups

**Purpose**: Organize and manage related workflow runs

**Benefits**:
- Set resource limits (max CPUs, max runs)
- Group by project or sample cohort
- Track costs together
- Manage quotas

---

## 🔄 Complete Workflow Lifecycle

### Detailed Step-by-Step

#### Phase 1: Preparation

1. **Identify Analysis Need**
   - What analysis? (alignment, variant calling, RNA-Seq)
   - Which tools? (BWA, GATK, STAR)
   - What inputs needed?

2. **Prepare Containers**
   - Use existing containers from registries OR
   - Build custom containers
   - Push to ECR

3. **Write Workflow Definition**
   - Choose language (WDL/Nextflow)
   - Define tasks
   - Specify inputs/outputs
   - Set dependencies

#### Phase 2: Workflow Creation

4. **Bundle Workflow**
   - Create zip file with workflow files
   - Upload to S3
   - Create parameters.json

5. **Create Workflow in HealthOmics**
   - Call CreateWorkflow API
   - Provide S3 URI of zip bundle
   - Get Workflow ID back

6. **Verify Workflow**
   - Check workflow is ACTIVE
   - Review parameters
   - Test with sample data

#### Phase 3: Execution

7. **Prepare Input Data**
   - Upload data to S3 OR
   - Use data from HealthOmics stores
   - Create inputs.json

8. **Configure Run**
   - Set workflow parameters
   - Specify output location
   - Set compute resources (optional)
   - Assign to run group (optional)

9. **Start Run**
   - Call StartRun API
   - Provide all inputs
   - Get Run ID back

10. **Monitor Execution**
    - Check run status
    - View task progress
    - Monitor resource usage
    - Check CloudWatch logs

#### Phase 4: Results

11. **Retrieve Outputs**
    - Outputs saved to specified S3 location
    - Organized by run ID
    - Includes task outputs
    - Logs and metrics available

12. **Analysis & Iteration**
    - Review results
    - Optimize if needed
    - Run additional samples
    - Scale to cohorts

---

## 🎯 Key Concepts Explained

### Workflow vs Run

**Workflow**:
- The **definition** or **blueprint**
- Created once
- Reusable
- Contains logic and structure

**Run**:
- An **execution** of a workflow
- Created each time you run
- Uses specific input data
- Produces specific outputs

**Analogy**: 
- Workflow = Recipe
- Run = Actually cooking the meal

### Task Execution

**How Tasks Run:**

1. **Task Isolation**
   - Each task runs in its own container
   - Clean environment
   - No interference between tasks

2. **Parallel Execution**
   - Independent tasks run in parallel
   - Automatic by HealthOmics
   - Maximizes throughput

3. **Sequential Execution**
   - Dependent tasks wait for prerequisites
   - Automatic dependency management
   - Data passed through shared file system

4. **Resource Allocation**
   - CPU, memory, storage per task
   - Can be specified in workflow
   - Auto-scaled by HealthOmics

### Shared File System

**Purpose**: Efficient data sharing between tasks

**How it Works:**
```
Task 1 Output → Shared FS → Task 2 Input
                  ↓
            No S3 upload/download
            Faster execution
```

**Benefits:**
- No intermediate S3 uploads
- Faster task transitions
- Lower costs
- Automatic cleanup

### Output Management

**outputUri**: S3 location for final results

**Structure:**
```
s3://bucket/outputs/run-id/
├── task1/
│   ├── output1.bam
│   └── output1.bam.bai
├── task2/
│   ├── output2.vcf
│   └── output2.vcf.idx
└── final/
    └── results.vcf.gz
```

---

## 💡 Real-World Example: Variant Calling Pipeline

### The Analysis Pipeline

```
FASTQ files (raw reads)
    ↓
[Task 1] BWA alignment → BAM file
    ↓
[Task 2] Sort BAM → Sorted BAM
    ↓
[Task 3] Mark duplicates → Deduped BAM
    ↓
[Task 4] GATK HaplotypeCaller → VCF file
    ↓
Final VCF output
```

### Workflow Components

**Containers Needed:**
- BWA: `biocontainers/bwa:0.7.17`
- Samtools: `biocontainers/samtools:1.15`
- GATK: `broadinstitute/gatk:4.4.0.0`

**Inputs:**
- FASTQ R1 file
- FASTQ R2 file
- Reference genome
- Known variants (for GATK)

**Outputs:**
- Aligned BAM file
- Variant calls (VCF)
- QC metrics

**Parameters:**
- Number of threads
- Memory allocation
- Quality thresholds

### Execution Flow

1. Upload FASTQ files to S3
2. Reference genome already in HealthOmics Reference Store
3. Start workflow run with inputs
4. HealthOmics:
   - Provisions compute
   - Runs BWA alignment (parallel if multiple samples)
   - Sorts BAM files
   - Marks duplicates
   - Calls variants
   - Uploads results to S3
5. Results ready in ~2-4 hours (depending on data size)

---

## 📊 Resource Management

### Compute Resources

**Automatic Scaling:**
- HealthOmics automatically provisions resources
- Based on task requirements
- Scales up for parallel tasks
- Scales down when idle

**Resource Limits:**
- Set via run groups
- Max CPUs across all runs
- Max concurrent runs
- Prevents runaway costs

### Cost Optimization

**Tips:**
1. Use run groups to set limits
2. Right-size task resources
3. Use spot-like pricing (automatic)
4. Delete intermediate files
5. Optimize workflow logic

**Pricing Model:**
- Pay per vCPU-hour
- Pay per GB-hour (memory)
- Pay per GB (storage accessed)
- No charges when idle

---

## 🔍 Monitoring & Debugging

### Run Status

**States:**
- PENDING: Queued, waiting for resources
- STARTING: Initializing
- RUNNING: Executing tasks
- COMPLETED: Successfully finished
- FAILED: Encountered error
- CANCELLED: User cancelled

### Task Monitoring

**Per-Task Information:**
- Status (PENDING, RUNNING, COMPLETED, FAILED)
- Start time
- End time
- Resources used (CPU, memory)
- Logs

### CloudWatch Logs

**Log Groups:**
- `/aws/omics/WorkflowLog`
- One log stream per run
- Contains task outputs
- Includes error messages

**Viewing Logs:**
```bash
aws logs get-log-events \
  --log-group-name /aws/omics/WorkflowLog \
  --log-stream-name run-id
```

### Debugging Failed Runs

**Steps:**
1. Check run status and error message
2. View CloudWatch logs
3. Identify failed task
4. Review task logs
5. Check input data validity
6. Verify container availability
7. Test with smaller dataset

---

## 🎓 Best Practices

### Workflow Design

1. **Modular Tasks**: Break into small, reusable tasks
2. **Clear Dependencies**: Explicit task ordering
3. **Error Handling**: Retry logic for transient failures
4. **Validation**: Check inputs before processing
5. **Documentation**: Comment your workflow code

### Container Management

1. **Pin Versions**: Use specific tags (not `latest`)
2. **Test Locally**: Verify containers work before deployment
3. **Optimize Size**: Smaller containers = faster startup
4. **Use Registries**: Public registries for common tools
5. **Private ECR**: For custom or proprietary tools

### Data Management

1. **Organize Inputs**: Clear naming conventions
2. **Stage Data**: Upload to S3 before running
3. **Clean Outputs**: Delete unnecessary intermediates
4. **Compression**: Use .gz for large files
5. **Metadata**: Tag data with sample IDs

### Cost Management

1. **Test Small**: Use subset of data for testing
2. **Monitor Spending**: CloudWatch billing alarms
3. **Resource Limits**: Set run group limits
4. **Optimize Workflows**: Remove unnecessary steps
5. **Batch Processing**: Group similar samples

---

## 🚀 Next Steps

1. **Learn Workflow Languages**
   - WDL tutorial: https://openwdl.org/
   - Nextflow tutorial: https://nextflow.io/

2. **Explore Example Workflows**
   - AWS samples: https://github.com/aws-samples/amazon-omics-tutorials
   - GATK workflows: https://github.com/gatk-workflows
   - nf-core pipelines: https://nf-co.re/

3. **Start Simple**
   - Begin with single-task workflow
   - Test with small data
   - Gradually increase complexity

4. **Scale Up**
   - Add more tasks
   - Process more samples
   - Optimize performance

---

## 📚 Additional Resources

- **HealthOmics Workflows Docs**: https://docs.aws.amazon.com/omics/latest/dev/workflows.html
- **WDL Specification**: https://github.com/openwdl/wdl
- **Nextflow Docs**: https://www.nextflow.io/docs/latest/
- **Container Registries**:
  - BioContainers: https://biocontainers.pro/
  - Quay.io: https://quay.io/
  - Docker Hub: https://hub.docker.com/

---

This guide covers the complete workflow lifecycle from container preparation to result analysis. Ready to implement? Let's create the actual code!