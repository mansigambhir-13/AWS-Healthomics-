#!/usr/bin/env python3
"""
AWS HealthOmics Workflows - Quick Start Implementation
Complete working example based on the 3-step architecture diagram
"""

import boto3
import json
import time
from datetime import datetime

# ============================================================================
# STEP 1: CONTAINER SETUP (ECR)
# ============================================================================

def step1_explore_containers():
    """
    STEP 1 from diagram: Tooling Container → Amazon ECR
    
    Explores available containers for bioinformatics workflows
    """
    print("="*70)
    print("STEP 1: CONTAINER SETUP")
    print("="*70)
    
    ecr = boto3.client('ecr')
    
    print("\n1. Checking your ECR repositories...")
    try:
        repos = ecr.describe_repositories()
        if repos['repositories']:
            print(f"✓ Found {len(repos['repositories'])} private repositories")
            for repo in repos['repositories'][:3]:
                print(f"  • {repo['repositoryName']}")
        else:
            print("  No private repositories found")
    except Exception as e:
        print(f"  Cannot access ECR: {e}")
    
    print("\n2. Public containers you can use:")
    containers = {
        "Alignment": {
            "BWA": "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8",
            "Bowtie2": "quay.io/biocontainers/bowtie2:2.4.5--py39hd2f7db1_3"
        },
        "Variant Calling": {
            "GATK": "broadinstitute/gatk:4.4.0.0",
            "FreeBayes": "quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2"
        },
        "SAM/BAM Tools": {
            "Samtools": "quay.io/biocontainers/samtools:1.15--h3843a85_0",
            "Picard": "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"
        }
    }
    
    for category, tools in containers.items():
        print(f"\n  {category}:")
        for tool, image in tools.items():
            print(f"    • {tool}: {image}")
    
    print("\n✓ STEP 1 Complete: Containers identified")
    print("="*70)
    return containers


# ============================================================================
# STEP 2: WORKFLOW CREATION
# ============================================================================

def step2_create_workflow():
    """
    STEP 2 from diagram: WDL/Nextflow → zip bundle → CreateWorkflow
    
    Creates a simple alignment workflow and uploads it to HealthOmics
    """
    print("\n" + "="*70)
    print("STEP 2: WORKFLOW CREATION")
    print("="*70)
    
    # 2.1: Create WDL workflow definition
    print("\n1. Creating WDL workflow definition...")
    
    wdl_content = '''version 1.0

# Simple BWA Alignment Workflow
workflow AlignmentWorkflow {
    input {
        File reference_fasta
        File input_fastq_r1
        File input_fastq_r2
        String sample_name
        Int threads = 4
    }
    
    call BwaAlign {
        input:
            reference = reference_fasta,
            fastq_r1 = input_fastq_r1,
            fastq_r2 = input_fastq_r2,
            sample = sample_name,
            cpu = threads
    }
    
    call SortBam {
        input:
            input_bam = BwaAlign.aligned_bam,
            sample = sample_name
    }
    
    output {
        File sorted_bam = SortBam.sorted_bam
        File bam_index = SortBam.bam_index
    }
}

task BwaAlign {
    input {
        File reference
        File fastq_r1
        File fastq_r2
        String sample
        Int cpu
    }
    
    command <<<
        bwa mem -t ~{cpu} -R '@RG\\tID:~{sample}\\tSM:~{sample}' \\
            ~{reference} ~{fastq_r1} ~{fastq_r2} | \\
            samtools view -bS - > ~{sample}.bam
    >>>
    
    output {
        File aligned_bam = "~{sample}.bam"
    }
    
    runtime {
        docker: "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
        memory: "16GB"
        cpu: cpu
    }
}

task SortBam {
    input {
        File input_bam
        String sample
    }
    
    command <<<
        samtools sort -o ~{sample}.sorted.bam ~{input_bam}
        samtools index ~{sample}.sorted.bam
    >>>
    
    output {
        File sorted_bam = "~{sample}.sorted.bam"
        File bam_index = "~{sample}.sorted.bam.bai"
    }
    
    runtime {
        docker: "quay.io/biocontainers/samtools:1.15--h3843a85_0"
        memory: "8GB"
        cpu: 2
    }
}
'''
    
    # Save workflow file
    import os
    os.makedirs('workflow_bundle', exist_ok=True)
    
    with open('workflow_bundle/alignment.wdl', 'w') as f:
        f.write(wdl_content)
    
    print("  ✓ Created: workflow_bundle/alignment.wdl")
    
    # 2.2: Create parameters.json
    print("\n2. Creating parameters.json...")
    
    parameters = {
        "reference_fasta": {
            "description": "Reference genome FASTA file",
            "optional": False
        },
        "input_fastq_r1": {
            "description": "FASTQ R1 reads",
            "optional": False
        },
        "input_fastq_r2": {
            "description": "FASTQ R2 reads",
            "optional": False
        },
        "sample_name": {
            "description": "Sample identifier",
            "optional": False
        },
        "threads": {
            "description": "Number of CPU threads",
            "optional": True
        }
    }
    
    with open('workflow_bundle/parameters.json', 'w') as f:
        json.dump(parameters, f, indent=2)
    
    print("  ✓ Created: workflow_bundle/parameters.json")
    
    # 2.3: Create zip bundle
    print("\n3. Creating zip bundle...")
    
    import zipfile
    bundle_name = 'alignment_workflow.zip'
    
    with zipfile.ZipFile(bundle_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write('workflow_bundle/alignment.wdl', 'alignment.wdl')
        zipf.write('workflow_bundle/parameters.json', 'parameters.json')
    
    print(f"  ✓ Created: {bundle_name}")
    print(f"  Size: {os.path.getsize(bundle_name)} bytes")
    
    # 2.4: Upload to S3 (you need to provide bucket name)
    print("\n4. Ready to upload to S3 and create workflow")
    print("\nTo complete STEP 2, run:")
    print("""
    # Upload bundle to S3
    aws s3 cp alignment_workflow.zip s3://YOUR-BUCKET/workflows/
    
    # Create workflow in HealthOmics
    aws omics create-workflow \\
      --name "AlignmentWorkflow" \\
      --engine WDL \\
      --definition-uri s3://YOUR-BUCKET/workflows/alignment_workflow.zip \\
      --region us-east-1
    
    # This returns: Workflow ID
    """)
    
    print("\n✓ STEP 2 Complete: Workflow defined and bundled")
    print("="*70)
    
    return bundle_name


def step2_create_workflow_in_omics(bucket_name):
    """
    Actually create the workflow in HealthOmics (if you have S3 bucket)
    """
    print("\n" + "="*70)
    print("STEP 2: CREATING WORKFLOW IN HEALTHOMICS")
    print("="*70)
    
    s3 = boto3.client('s3')
    omics = boto3.client('omics')
    
    bundle_name = 'alignment_workflow.zip'
    
    # Upload to S3
    print(f"\n1. Uploading workflow bundle to s3://{bucket_name}/workflows/")
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    s3_key = f"workflows/{timestamp}/{bundle_name}"
    
    try:
        s3.upload_file(bundle_name, bucket_name, s3_key)
        s3_uri = f"s3://{bucket_name}/{s3_key}"
        print(f"  ✓ Uploaded to: {s3_uri}")
    except Exception as e:
        print(f"  ✗ Upload failed: {e}")
        print(f"\n  Please upload manually:")
        print(f"  aws s3 cp {bundle_name} s3://{bucket_name}/workflows/")
        return None
    
    # Create workflow
    print("\n2. Creating workflow in HealthOmics...")
    
    try:
        response = omics.create_workflow(
            name='AlignmentWorkflow',
            description='Simple BWA alignment workflow',
            engine='WDL',
            definitionUri=s3_uri,
            tags={
                'Project': 'Tutorial',
                'CreatedBy': 'QuickStart'
            }
        )
        
        workflow_id = response['id']
        print(f"  ✓ Workflow created!")
        print(f"  Workflow ID: {workflow_id}")
        print(f"  ARN: {response['arn']}")
        print(f"  Status: {response.get('status', 'CREATING')}")
        
        # Wait for workflow to be active
        print("\n3. Waiting for workflow to become ACTIVE...")
        for i in range(30):
            time.sleep(5)
            status_response = omics.get_workflow(id=workflow_id)
            status = status_response.get('status', 'UNKNOWN')
            print(f"  [{i+1}/30] Status: {status}", end='\r')
            
            if status == 'ACTIVE':
                print(f"\n  ✓ Workflow is ACTIVE and ready!")
                break
            elif status in ['FAILED', 'DELETED']:
                print(f"\n  ✗ Workflow entered {status} state")
                return None
        
        print("\n✓ STEP 2 Complete: Workflow ready to use!")
        print("="*70)
        
        return workflow_id
        
    except Exception as e:
        print(f"  ✗ Failed to create workflow: {e}")
        return None


# ============================================================================
# STEP 3: WORKFLOW EXECUTION
# ============================================================================

def step3_run_workflow(workflow_id, role_arn, output_bucket):
    """
    STEP 3 from diagram: Workflow ID + Role + Inputs → Start Run
    
    Starts a workflow run and monitors it
    """
    print("\n" + "="*70)
    print("STEP 3: WORKFLOW EXECUTION")
    print("="*70)
    
    omics = boto3.client('omics')
    
    # 3.1: Prepare inputs
    print("\n1. Preparing inputs.json...")
    
    inputs = {
        "AlignmentWorkflow.reference_fasta": "s3://YOUR-BUCKET/reference/hg38.fa",
        "AlignmentWorkflow.input_fastq_r1": "s3://YOUR-BUCKET/samples/sample_R1.fastq.gz",
        "AlignmentWorkflow.input_fastq_r2": "s3://YOUR-BUCKET/samples/sample_R2.fastq.gz",
        "AlignmentWorkflow.sample_name": "SAMPLE001",
        "AlignmentWorkflow.threads": 8
    }
    
    print("  Sample inputs:")
    for key, value in inputs.items():
        print(f"    {key}: {value}")
    
    # 3.2: Start run
    print("\n2. Starting workflow run...")
    
    run_name = f"alignment-run-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    output_uri = f"s3://{output_bucket}/outputs/{run_name}/"
    
    print(f"  Run name: {run_name}")
    print(f"  Output: {output_uri}")
    
    try:
        # Note: You need real S3 URIs for this to work
        print("\n  To start the run, execute:")
        print(f"""
    omics = boto3.client('omics')
    
    response = omics.start_run(
        workflowId='{workflow_id}',
        roleArn='{role_arn}',
        name='{run_name}',
        parameters={{
            'AlignmentWorkflow.reference_fasta': 's3://bucket/reference.fa',
            'AlignmentWorkflow.input_fastq_r1': 's3://bucket/sample_R1.fastq.gz',
            'AlignmentWorkflow.input_fastq_r2': 's3://bucket/sample_R2.fastq.gz',
            'AlignmentWorkflow.sample_name': 'SAMPLE001',
            'AlignmentWorkflow.threads': 8
        }},
        outputUri='{output_uri}'
    )
    
    run_id = response['id']
    print(f"Run started: {{run_id}}")
        """)
        
        print("\n3. Run monitoring:")
        print("""
    # Monitor run status
    while True:
        response = omics.get_run(id=run_id)
        status = response['status']
        print(f"Status: {status}")
        
        if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
            break
        
        time.sleep(30)
        """)
        
        print("\n4. Task inspection:")
        print("""
    # View individual tasks
    response = omics.get_run(id=run_id)
    for task in response.get('tasks', []):
        print(f"Task: {task['name']}")
        print(f"  Status: {task['status']}")
        print(f"  CPUs: {task.get('cpus', 'N/A')}")
        print(f"  Memory: {task.get('memory', 'N/A')} GB")
        """)
        
        print("\n5. View logs:")
        print("""
    # Get CloudWatch logs
    logs = boto3.client('logs')
    
    response = logs.get_log_events(
        logGroupName='/aws/omics/WorkflowLog',
        logStreamName=run_id,
        limit=100
    )
    
    for event in response['events']:
        print(event['message'])
        """)
        
        print("\n✓ STEP 3 Complete: Workflow execution demonstrated")
        print("="*70)
        
    except Exception as e:
        print(f"  Note: {e}")


# ============================================================================
# COMPLETE EXAMPLE WITH REAL CODE
# ============================================================================

def complete_example_with_monitoring(workflow_id, role_arn, parameters, output_uri):
    """
    Complete working example with real monitoring
    """
    print("\n" + "="*70)
    print("COMPLETE WORKFLOW EXECUTION WITH MONITORING")
    print("="*70)
    
    omics = boto3.client('omics')
    logs = boto3.client('logs')
    
    # Start run
    print("\n1. Starting workflow run...")
    
    try:
        response = omics.start_run(
            workflowId=workflow_id,
            roleArn=role_arn,
            name=f"run-{datetime.now().strftime('%Y%m%d%H%M%S')}",
            parameters=parameters,
            outputUri=output_uri
        )
        
        run_id = response['id']
        print(f"✓ Run started: {run_id}")
        
        # Monitor run
        print("\n2. Monitoring run progress...")
        
        attempt = 0
        while True:
            attempt += 1
            response = omics.get_run(id=run_id)
            status = response['status']
            
            print(f"\n[Check #{attempt}] Status: {status}")
            
            # Show task progress
            tasks = response.get('tasks', [])
            if tasks:
                completed = sum(1 for t in tasks if t.get('status') == 'COMPLETED')
                running = sum(1 for t in tasks if t.get('status') == 'RUNNING')
                total = len(tasks)
                
                print(f"  Tasks: {completed}/{total} completed")
                if running > 0:
                    print(f"  Currently running: {running} tasks")
                
                # Show running tasks
                for task in tasks:
                    if task.get('status') == 'RUNNING':
                        print(f"    → {task.get('name', 'Unknown')}")
            
            # Check if done
            if status == 'COMPLETED':
                print(f"\n✓ Workflow COMPLETED successfully!")
                
                # Show results
                print(f"\nResults location: {response.get('outputUri', 'N/A')}")
                
                if 'resourceDigests' in response:
                    print("\nResource Usage:")
                    for key, value in response['resourceDigests'].items():
                        print(f"  {key}: {value}")
                
                break
                
            elif status == 'FAILED':
                print(f"\n✗ Workflow FAILED")
                error = response.get('statusMessage', 'No error message')
                print(f"Error: {error}")
                
                # Show failed tasks
                for task in tasks:
                    if task.get('status') == 'FAILED':
                        print(f"\nFailed task: {task.get('name', 'Unknown')}")
                        print(f"  Error: {task.get('statusMessage', 'N/A')}")
                
                break
                
            elif status == 'CANCELLED':
                print(f"\n⚠ Workflow was cancelled")
                break
            
            # Wait before checking again
            time.sleep(30)
        
        # Get logs
        print("\n3. Retrieving logs...")
        try:
            log_response = logs.get_log_events(
                logGroupName='/aws/omics/WorkflowLog',
                logStreamName=run_id,
                limit=50
            )
            
            events = log_response.get('events', [])
            if events:
                print(f"\nShowing last {len(events)} log entries:\n")
                for event in events[-10:]:  # Last 10
                    timestamp = datetime.fromtimestamp(event['timestamp']/1000)
                    print(f"[{timestamp}] {event['message']}")
            else:
                print("  No logs available yet")
                
        except Exception as e:
            print(f"  Could not retrieve logs: {e}")
        
        print("\n" + "="*70)
        print("✓ Complete example finished!")
        print("="*70)
        
        return run_id
        
    except Exception as e:
        print(f"✗ Error: {e}")
        return None


# ============================================================================
# MAIN DEMONSTRATION
# ============================================================================

def main():
    """
    Main demonstration of all 3 steps
    """
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║          AWS HealthOmics Workflows - Complete Implementation         ║
║                     Based on 3-Step Architecture                     ║
╚══════════════════════════════════════════════════════════════════════╝

This script demonstrates the complete workflow lifecycle:

STEP 1: Container Setup (ECR)
  - List available containers
  - Show public bioinformatics containers
  
STEP 2: Workflow Creation
  - Create WDL workflow definition
  - Bundle with parameters.json
  - Upload to S3
  - Create workflow in HealthOmics
  
STEP 3: Workflow Execution
  - Start workflow run
  - Monitor progress
  - Inspect tasks
  - View logs
  - Retrieve results

Let's begin!
""")
    
    input("Press Enter to start STEP 1...")
    
    # STEP 1: Explore containers
    containers = step1_explore_containers()
    
    input("\nPress Enter to continue to STEP 2...")
    
    # STEP 2: Create workflow
    bundle_name = step2_create_workflow()
    
    print(f"""
╔══════════════════════════════════════════════════════════════════════╗
║                       NEXT STEPS TO COMPLETE                         ║
╚══════════════════════════════════════════════════════════════════════╝

You now have:
  ✓ Workflow definition (WDL)
  ✓ Parameters file (JSON)
  ✓ Bundle file: {bundle_name}

To continue:

1. Upload to S3:
   aws s3 cp {bundle_name} s3://YOUR-BUCKET/workflows/

2. Create workflow:
   python3
   >>> from {__file__.split('/')[-1].replace('.py', '')} import step2_create_workflow_in_omics
   >>> workflow_id = step2_create_workflow_in_omics('YOUR-BUCKET')

3. Run workflow:
   >>> from {__file__.split('/')[-1].replace('.py', '')} import complete_example_with_monitoring
   >>> run_id = complete_example_with_monitoring(
   ...     workflow_id='YOUR-WORKFLOW-ID',
   ...     role_arn='arn:aws:iam::123456789012:role/OmicsRole',
   ...     parameters={{...}},
   ...     output_uri='s3://bucket/outputs/'
   ... )

Or use the full implementation in workflows_implementation.py!
""")


if __name__ == '__main__':
    main()