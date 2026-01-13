# AWS HealthOmics Workflows - Quick Reference

## 📊 The 3-Step Architecture

```
STEP 1: Container → ECR
STEP 2: Workflow Definition → CreateWorkflow → Workflow ID
STEP 3: Workflow ID + Inputs → StartRun → Results
```

---

## 🚀 Quick Commands

### List Available Workflows
```bash
aws omics list-workflows --region us-east-1
```

### Create Workflow
```bash
aws omics create-workflow \
  --name "MyWorkflow" \
  --engine WDL \
  --definition-uri s3://bucket/workflow.zip \
  --region us-east-1
```

### Start Workflow Run
```bash
aws omics start-run \
  --workflow-id 1234567 \
  --role-arn arn:aws:iam::123456789012:role/OmicsRole \
  --parameters file://inputs.json \
  --output-uri s3://bucket/outputs/ \
  --region us-east-1
```

### Check Run Status
```bash
aws omics get-run --id run-123456 --region us-east-1
```

### List Runs
```bash
aws omics list-runs --status RUNNING --region us-east-1
```

### Get Run Logs
```bash
aws logs get-log-events \
  --log-group-name /aws/omics/WorkflowLog \
  --log-stream-name run-123456
```

---

## 💻 Python Quick Start

### Initialize
```python
import boto3
omics = boto3.client('omics', region_name='us-east-1')
```

### Create Workflow
```python
response = omics.create_workflow(
    name='MyWorkflow',
    engine='WDL',
    definitionUri='s3://bucket/workflow.zip'
)
workflow_id = response['id']
```

### Start Run
```python
response = omics.start_run(
    workflowId=workflow_id,
    roleArn='arn:aws:iam::123456789012:role/OmicsRole',
    parameters={
        'input_file': 's3://bucket/input.fastq',
        'reference': 'GRCh38'
    },
    outputUri='s3://bucket/outputs/'
)
run_id = response['id']
```

### Monitor Run
```python
while True:
    response = omics.get_run(id=run_id)
    status = response['status']
    print(f"Status: {status}")
    
    if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
        break
    
    time.sleep(30)
```

---

## 📝 Workflow File Structure

### WDL Workflow
```wdl
version 1.0

workflow MyWorkflow {
    input {
        File input_file
        String param
    }
    
    call MyTask {
        input:
            file = input_file,
            parameter = param
    }
    
    output {
        File result = MyTask.output_file
    }
}

task MyTask {
    input {
        File file
        String parameter
    }
    
    command <<<
        tool --input ~{file} --param ~{parameter}
    >>>
    
    output {
        File output_file = "result.txt"
    }
    
    runtime {
        docker: "container-image:tag"
        memory: "8GB"
        cpu: 4
    }
}
```

### Nextflow Workflow
```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = null
params.param = 'default'

process MyTask {
    cpus 4
    memory '8 GB'
    container 'container-image:tag'
    
    input:
    path input_file
    val param
    
    output:
    path 'result.txt'
    
    script:
    """
    tool --input ${input_file} --param ${param}
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)
    param_ch = Channel.of(params.param)
    
    MyTask(input_ch, param_ch)
}
```

---

## 🐳 Container Quick Reference

### Popular Bioinformatics Containers

**Alignment:**
- BWA: `quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8`
- Bowtie2: `quay.io/biocontainers/bowtie2:2.4.5--py39hd2f7db1_3`

**Variant Calling:**
- GATK: `broadinstitute/gatk:4.4.0.0`
- FreeBayes: `quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2`

**RNA-Seq:**
- STAR: `quay.io/biocontainers/star:2.7.10b--h6b7c446_1`
- Salmon: `quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1`

**Utilities:**
- Samtools: `quay.io/biocontainers/samtools:1.15--h3843a85_0`
- BCFtools: `quay.io/biocontainers/bcftools:1.15--h0ea216a_0`

---

## 🔧 Common Workflow Patterns

### Sequential Tasks
```wdl
workflow Sequential {
    call Task1
    call Task2 { input: data = Task1.output }
    call Task3 { input: data = Task2.output }
}
```

### Parallel Tasks
```wdl
workflow Parallel {
    scatter (sample in samples) {
        call ProcessSample { input: sample = sample }
    }
}
```

### Conditional Execution
```wdl
workflow Conditional {
    if (run_qc) {
        call QualityControl
    }
    call MainAnalysis
}
```

---

## 📊 Input/Output Patterns

### inputs.json Example
```json
{
    "WorkflowName.input_file": "s3://bucket/input.fastq",
    "WorkflowName.reference": "s3://bucket/reference.fa",
    "WorkflowName.sample_name": "SAMPLE001",
    "WorkflowName.threads": 8
}
```

### Output Structure
```
s3://bucket/outputs/run-id/
├── task1/
│   └── output1.bam
├── task2/
│   └── output2.vcf
└── logs/
    └── task_logs.txt
```

---

## 🎯 Run Group Management

### Create Run Group
```python
response = omics.create_run_group(
    name='ProjectA',
    maxCpus=100000,
    maxRuns=100
)
run_group_id = response['id']
```

### Start Run in Group
```python
omics.start_run(
    workflowId=workflow_id,
    runGroupId=run_group_id,
    parameters=params,
    outputUri='s3://bucket/outputs/'
)
```

---

## 🔍 Monitoring & Debugging

### Check Run Status
```python
response = omics.get_run(id=run_id)
print(f"Status: {response['status']}")
print(f"Tasks: {len(response.get('tasks', []))}")
```

### View Failed Tasks
```python
response = omics.get_run(id=run_id)
for task in response.get('tasks', []):
    if task['status'] == 'FAILED':
        print(f"Failed: {task['name']}")
        print(f"Error: {task.get('statusMessage')}")
```

### Get Logs
```python
logs = boto3.client('logs')
response = logs.get_log_events(
    logGroupName='/aws/omics/WorkflowLog',
    logStreamName=run_id
)
for event in response['events']:
    print(event['message'])
```

---

## 💰 Cost Optimization

### Tips:
1. Right-size CPU/memory in runtime blocks
2. Use spot-like pricing (automatic)
3. Set run group limits
4. Delete unnecessary intermediate files
5. Optimize workflow logic

### Monitor Costs:
```python
response = omics.get_run(id=run_id)
if 'resourceDigests' in response:
    print("Resource Usage:")
    for key, value in response['resourceDigests'].items():
        print(f"  {key}: {value}")
```

---

## ⚠️ Troubleshooting

### Common Errors

**"Container not found"**
```
Solution: Verify container URI and version tag
```

**"Out of memory"**
```
Solution: Increase memory in runtime block
runtime { memory: "32GB" }
```

**"Access Denied"**
```
Solution: Check IAM role permissions
- S3 read/write
- ECR pull
- HealthOmics run
```

**"File not found"**
```
Solution: Verify S3 URIs
aws s3 ls s3://bucket/file
```

---

## 🎓 Best Practices

1. **Pin Container Versions**: Use specific tags, not `latest`
2. **Test Locally First**: Use Cromwell (WDL) or Nextflow locally
3. **Start Small**: Test with subset of data
4. **Add Retries**: Handle transient failures
5. **Document Workflows**: Add comments and README
6. **Use Run Groups**: Organize and limit resources
7. **Monitor Costs**: Set up billing alerts
8. **Clean Up**: Delete old runs and intermediate files

---

## 📚 Resources

- **Workflows Docs**: https://docs.aws.amazon.com/omics/latest/dev/workflows.html
- **WDL Specification**: https://openwdl.org/
- **Nextflow Docs**: https://www.nextflow.io/
- **Example Workflows**: https://github.com/aws-samples/amazon-omics-tutorials
- **GATK Workflows**: https://github.com/gatk-workflows/
- **nf-core**: https://nf-co.re/

---

## 🚀 Quick Start Checklist

- [ ] Create workflow definition (WDL/Nextflow)
- [ ] Bundle workflow into zip file
- [ ] Upload bundle to S3
- [ ] Create workflow in HealthOmics
- [ ] Prepare input data in S3
- [ ] Create inputs.json file
- [ ] Create IAM execution role
- [ ] Start workflow run
- [ ] Monitor progress
- [ ] Retrieve results from S3

---

## 💡 Pro Tips

- **Use existing workflows**: Start with GATK or nf-core pipelines
- **Leverage run groups**: Better cost control and organization
- **Cache intermediate results**: Reuse computation when possible
- **Batch similar samples**: Process cohorts together
- **Monitor CloudWatch**: Set up alarms for failures
- **Version your workflows**: Track changes over time
- **Test incrementally**: Add complexity gradually

---

**Ready to run workflows?** Use this reference alongside the implementation code!