#!/usr/bin/env python3
"""
AWS HealthOmics Workflows - Complete Implementation
This module covers: Creating workflows, running them, monitoring tasks, and viewing logs
"""

import boto3
import json
import time
import zipfile
import os
from datetime import datetime
from botocore.exceptions import ClientError

class HealthOmicsWorkflowComplete:
    """
    Complete implementation of HealthOmics Workflows
    
    Based on the 3-step architecture:
    1. Container setup (ECR)
    2. Workflow creation (CreateWorkflow)
    3. Workflow execution (StartRun)
    """
    
    def __init__(self, region='us-east-1'):
        """Initialize AWS clients"""
        self.region = region
        self.omics = boto3.client('omics', region_name=region)
        self.s3 = boto3.client('s3', region_name=region)
        self.logs = boto3.client('logs', region_name=region)
        self.ecr = boto3.client('ecr', region_name=region)
        
        # Get account info
        self.account_id = boto3.client('sts').get_caller_identity()['Account']
        
        print(f"Initialized HealthOmics Workflows Manager")
        print(f"Account: {self.account_id}")
        print(f"Region: {region}\n")
    
    # ========================================================================
    # STEP 1: CONTAINER MANAGEMENT
    # ========================================================================
    
    def list_ecr_repositories(self):
        """
        List ECR repositories that contain workflow containers
        
        Containers are Docker images containing bioinformatics tools
        Examples: BWA, GATK, STAR, Samtools
        """
        print("="*70)
        print("ECR Repositories (Container Registry)")
        print("="*70)
        
        try:
            response = self.ecr.describe_repositories()
            repos = response.get('repositories', [])
            
            if repos:
                print(f"\nFound {len(repos)} repositories:\n")
                for repo in repos:
                    print(f"Repository: {repo['repositoryName']}")
                    print(f"  URI: {repo['repositoryUri']}")
                    print(f"  Created: {repo.get('createdAt', 'N/A')}")
                    print()
            else:
                print("\nNo ECR repositories found.")
                print("You can use public registries like:")
                print("  - quay.io/biocontainers/")
                print("  - docker.io/biocontainers/")
                print("  - broadinstitute/gatk")
        
        except ClientError as e:
            print(f"Error listing ECR repositories: {e}")
    
    def get_public_container_examples(self):
        """Show examples of public bioinformatics containers"""
        print("\n" + "="*70)
        print("Public Container Examples")
        print("="*70)
        
        containers = {
            "Alignment": [
                "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8",
                "quay.io/biocontainers/bowtie2:2.4.5--py39hd2f7db1_3"
            ],
            "SAM/BAM Processing": [
                "quay.io/biocontainers/samtools:1.15--h3843a85_0",
                "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"
            ],
            "Variant Calling": [
                "broadinstitute/gatk:4.4.0.0",
                "quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2"
            ],
            "RNA-Seq": [
                "quay.io/biocontainers/star:2.7.10b--h6b7c446_1",
                "quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1"
            ],
            "Quality Control": [
                "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1",
                "quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0"
            ]
        }
        
        for category, images in containers.items():
            print(f"\n{category}:")
            for image in images:
                print(f"  • {image}")
        
        print("\n" + "="*70)
    
    # ========================================================================
    # STEP 2: WORKFLOW CREATION
    # ========================================================================
    
    def create_simple_wdl_workflow(self):
        """
        Create a simple WDL workflow example
        
        This workflow demonstrates:
        - Basic WDL syntax
        - Task definition
        - Workflow structure
        - Input/output specification
        """
        print("="*70)
        print("Creating Example WDL Workflow")
        print("="*70)
        
        # Simple alignment workflow
        wdl_content = '''version 1.0

# Simple BWA Alignment Workflow
# Aligns FASTQ reads to a reference genome

workflow SimpleAlignment {
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
            sample = sample_name,
            cpu = threads
    }
    
    output {
        File sorted_bam = SortBam.sorted_bam
        File sorted_bam_index = SortBam.bam_index
    }
}

task BwaAlign {
    input {
        File reference
        File fastq_r1
        File fastq_r2
        String sample
        Int cpu
        Int memory_gb = 16
    }
    
    command <<<
        set -e
        
        # Index reference if needed
        bwa index ~{reference}
        
        # Perform alignment
        bwa mem -t ~{cpu} -R '@RG\\tID:~{sample}\\tSM:~{sample}\\tPL:ILLUMINA' \\
            ~{reference} ~{fastq_r1} ~{fastq_r2} | \\
            samtools view -bS - > ~{sample}.bam
    >>>
    
    output {
        File aligned_bam = "~{sample}.bam"
    }
    
    runtime {
        docker: "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
        memory: "~{memory_gb}GB"
        cpu: cpu
    }
}

task SortBam {
    input {
        File input_bam
        String sample
        Int cpu
        Int memory_gb = 8
    }
    
    command <<<
        samtools sort -@ ~{cpu} -o ~{sample}.sorted.bam ~{input_bam}
        samtools index ~{sample}.sorted.bam
    >>>
    
    output {
        File sorted_bam = "~{sample}.sorted.bam"
        File bam_index = "~{sample}.sorted.bam.bai"
    }
    
    runtime {
        docker: "quay.io/biocontainers/samtools:1.15--h3843a85_0"
        memory: "~{memory_gb}GB"
        cpu: cpu
    }
}
'''
        
        # Save workflow
        workflow_dir = 'healthomics_workflows'
        os.makedirs(workflow_dir, exist_ok=True)
        
        wdl_file = f'{workflow_dir}/simple_alignment.wdl'
        with open(wdl_file, 'w') as f:
            f.write(wdl_content)
        
        print(f"✓ WDL workflow created: {wdl_file}")
        
        # Create parameters file
        parameters = {
            "reference_fasta": {
                "description": "Reference genome in FASTA format",
                "optional": False
            },
            "input_fastq_r1": {
                "description": "FASTQ file with forward reads",
                "optional": False
            },
            "input_fastq_r2": {
                "description": "FASTQ file with reverse reads",
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
        
        params_file = f'{workflow_dir}/parameters.json'
        with open(params_file, 'w') as f:
            json.dump(parameters, f, indent=2)
        
        print(f"✓ Parameters file created: {params_file}")
        
        return wdl_file, params_file
    
    def create_workflow_bundle(self, workflow_file, params_file, output_zip='workflow_bundle.zip'):
        """
        Create a zip bundle for workflow upload
        
        The bundle contains:
        - Main workflow file (WDL/Nextflow)
        - Parameters file
        - Any sub-workflows or config files
        """
        print("\n" + "="*70)
        print("Creating Workflow Bundle")
        print("="*70)
        
        try:
            with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
                # Add workflow file
                zipf.write(workflow_file, os.path.basename(workflow_file))
                print(f"  Added: {workflow_file}")
                
                # Add parameters file
                if params_file:
                    zipf.write(params_file, os.path.basename(params_file))
                    print(f"  Added: {params_file}")
            
            print(f"\n✓ Workflow bundle created: {output_zip}")
            print(f"  Size: {os.path.getsize(output_zip)} bytes")
            
            return output_zip
            
        except Exception as e:
            print(f"✗ Error creating bundle: {e}")
            raise
    
    def upload_workflow_to_s3(self, bundle_file, bucket, key=None):
        """Upload workflow bundle to S3"""
        if key is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            key = f"workflows/{timestamp}/{bundle_file}"
        
        print(f"\nUploading to S3...")
        print(f"  Bucket: {bucket}")
        print(f"  Key: {key}")
        
        try:
            self.s3.upload_file(bundle_file, bucket, key)
            s3_uri = f"s3://{bucket}/{key}"
            print(f"✓ Uploaded: {s3_uri}")
            return s3_uri
        except ClientError as e:
            print(f"✗ Error uploading: {e}")
            raise
    
    def create_workflow(self, name, definition_uri, engine='WDL', description=""):
        """
        Create a workflow in HealthOmics
        
        This implements STEP 2 of the diagram:
        Workflow definition → CreateWorkflow → Workflow ID
        
        Args:
            name: Workflow name
            definition_uri: S3 URI of workflow bundle
            engine: 'WDL' or 'NEXTFLOW'
            description: Workflow description
        
        Returns:
            Workflow ID
        """
        print("\n" + "="*70)
        print("Creating Workflow in HealthOmics")
        print("="*70)
        print(f"Name: {name}")
        print(f"Definition: {definition_uri}")
        print(f"Engine: {engine}\n")
        
        try:
            response = self.omics.create_workflow(
                name=name,
                description=description or f"Workflow created at {datetime.now()}",
                engine=engine,
                definitionUri=definition_uri,
                tags={
                    'Project': 'HealthOmics-Workshop',
                    'CreatedBy': 'Workshop-Script'
                }
            )
            
            workflow_id = response['id']
            workflow_arn = response['arn']
            
            print(f"✓ Workflow created successfully!")
            print(f"  Workflow ID: {workflow_id}")
            print(f"  ARN: {workflow_arn}")
            print(f"  Status: {response.get('status', 'CREATING')}")
            
            # Wait for workflow to become active
            print("\nWaiting for workflow to become ACTIVE...")
            self._wait_for_workflow_active(workflow_id)
            
            return workflow_id
            
        except ClientError as e:
            print(f"✗ Error creating workflow: {e}")
            raise
    
    def _wait_for_workflow_active(self, workflow_id, max_attempts=30):
        """Wait for workflow to become ACTIVE"""
        for attempt in range(max_attempts):
            try:
                response = self.omics.get_workflow(id=workflow_id)
                status = response.get('status', 'UNKNOWN')
                
                print(f"  [{attempt+1}/{max_attempts}] Status: {status}", end='\r')
                
                if status == 'ACTIVE':
                    print("\n✓ Workflow is ACTIVE and ready to use")
                    return
                elif status in ['FAILED', 'DELETED']:
                    raise Exception(f"Workflow entered {status} state")
                
                time.sleep(5)
            except ClientError as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(5)
        
        raise Exception("Workflow creation timed out")
    
    # ========================================================================
    # STEP 3: WORKFLOW EXECUTION
    # ========================================================================
    
    def start_workflow_run(self,
                          workflow_id,
                          role_arn,
                          parameters,
                          output_uri,
                          name=None,
                          run_group_id=None,
                          priority=None):
        """
        Start a workflow run
        
        This implements STEP 3 of the diagram:
        Workflow ID + Inputs + Role → Start Run → Workflow Run
        
        Args:
            workflow_id: ID of the workflow to run
            role_arn: IAM role ARN for execution
            parameters: Dictionary of input parameters
            output_uri: S3 URI for outputs
            name: Run name (optional)
            run_group_id: Run group ID (optional)
            priority: Run priority 0-100000 (optional)
        
        Returns:
            Run ID
        """
        print("="*70)
        print("Starting Workflow Run")
        print("="*70)
        print(f"Workflow ID: {workflow_id}")
        print(f"Output URI: {output_uri}")
        
        if name is None:
            name = f"run-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
        
        print(f"Run Name: {name}\n")
        
        try:
            run_params = {
                'workflowId': workflow_id,
                'roleArn': role_arn,
                'name': name,
                'parameters': parameters,
                'outputUri': output_uri
            }
            
            if run_group_id:
                run_params['runGroupId'] = run_group_id
            if priority is not None:
                run_params['priority'] = priority
            
            response = self.omics.start_run(**run_params)
            
            run_id = response['id']
            
            print(f"✓ Workflow run started!")
            print(f"  Run ID: {run_id}")
            print(f"  ARN: {response['arn']}")
            print(f"  Status: {response.get('status', 'PENDING')}")
            
            return run_id
            
        except ClientError as e:
            print(f"✗ Error starting run: {e}")
            raise
    
    def get_run_status(self, run_id):
        """Get current status of a workflow run"""
        try:
            response = self.omics.get_run(id=run_id)
            return response
        except ClientError as e:
            print(f"Error getting run status: {e}")
            raise
    
    def monitor_run(self, run_id, poll_interval=30, show_tasks=True):
        """
        Monitor a workflow run until completion
        
        Shows:
        - Run status
        - Task progress
        - Resource usage
        - Logs (optional)
        """
        print("\n" + "="*70)
        print(f"Monitoring Run: {run_id}")
        print("="*70)
        
        attempt = 0
        while True:
            try:
                response = self.omics.get_run(id=run_id)
                status = response['status']
                attempt += 1
                
                # Print status update
                print(f"\n[Check #{attempt}] Status: {status}")
                print(f"  Started: {response.get('startTime', 'N/A')}")
                
                # Show task progress if available
                if show_tasks and 'tasks' in response:
                    self._print_task_status(response['tasks'])
                
                # Check for terminal states
                if status == 'COMPLETED':
                    print(f"\n✓ Workflow run COMPLETED successfully!")
                    self._print_run_summary(response)
                    return response
                    
                elif status == 'FAILED':
                    print(f"\n✗ Workflow run FAILED")
                    error_msg = response.get('statusMessage', 'No error message')
                    print(f"  Error: {error_msg}")
                    self._print_failed_tasks(response.get('tasks', []))
                    return response
                    
                elif status == 'CANCELLED':
                    print(f"\n⚠ Workflow run was CANCELLED")
                    return response
                
                # Continue monitoring
                time.sleep(poll_interval)
                
            except ClientError as e:
                print(f"Error monitoring run: {e}")
                time.sleep(poll_interval)
    
    def _print_task_status(self, tasks):
        """Print status of individual tasks"""
        if not tasks:
            return
        
        total = len(tasks)
        completed = sum(1 for t in tasks if t.get('status') == 'COMPLETED')
        running = sum(1 for t in tasks if t.get('status') == 'RUNNING')
        failed = sum(1 for t in tasks if t.get('status') == 'FAILED')
        
        print(f"\n  Task Progress: {completed}/{total} completed")
        if running > 0:
            print(f"    Running: {running}")
        if failed > 0:
            print(f"    Failed: {failed}")
        
        # Show running tasks
        for task in tasks:
            if task.get('status') == 'RUNNING':
                print(f"    • {task.get('name', 'Unknown')}: RUNNING")
    
    def _print_run_summary(self, run_info):
        """Print summary of completed run"""
        print("\nRun Summary:")
        print(f"  Name: {run_info.get('name', 'N/A')}")
        print(f"  Output: {run_info.get('outputUri', 'N/A')}")
        
        if 'startTime' in run_info and 'stopTime' in run_info:
            duration = run_info['stopTime'] - run_info['startTime']
            print(f"  Duration: {duration}")
        
        if 'resourceDigests' in run_info:
            print(f"\nResource Usage:")
            for key, value in run_info['resourceDigests'].items():
                print(f"  {key}: {value}")
    
    def _print_failed_tasks(self, tasks):
        """Print information about failed tasks"""
        failed = [t for t in tasks if t.get('status') == 'FAILED']
        if failed:
            print("\nFailed Tasks:")
            for task in failed:
                print(f"  • {task.get('name', 'Unknown')}")
                if 'statusMessage' in task:
                    print(f"    Error: {task['statusMessage']}")
    
    # ========================================================================
    # TASK INSPECTION
    # ========================================================================
    
    def get_run_tasks(self, run_id):
        """
        Get detailed information about all tasks in a run
        
        Shows:
        - Task names
        - Status
        - Start/end times
        - Resource usage
        - Logs
        """
        print("="*70)
        print(f"Tasks for Run: {run_id}")
        print("="*70)
        
        try:
            response = self.omics.get_run(id=run_id)
            tasks = response.get('tasks', [])
            
            if not tasks:
                print("\nNo tasks found (run may still be initializing)")
                return []
            
            print(f"\nFound {len(tasks)} tasks:\n")
            
            for i, task in enumerate(tasks, 1):
                print(f"Task {i}: {task.get('name', 'Unknown')}")
                print(f"  Status: {task.get('status', 'N/A')}")
                print(f"  Task ID: {task.get('taskId', 'N/A')}")
                
                if 'startTime' in task:
                    print(f"  Started: {task['startTime']}")
                if 'stopTime' in task:
                    print(f"  Stopped: {task['stopTime']}")
                
                if 'cpus' in task:
                    print(f"  CPUs: {task['cpus']}")
                if 'memory' in task:
                    print(f"  Memory: {task['memory']} GB")
                
                if task.get('status') == 'FAILED' and 'statusMessage' in task:
                    print(f"  Error: {task['statusMessage']}")
                
                print()
            
            return tasks
            
        except ClientError as e:
            print(f"Error getting tasks: {e}")
            raise
    
    # ========================================================================
    # LOG VIEWING
    # ========================================================================
    
    def get_run_logs(self, run_id, task_name=None, limit=100):
        """
        Retrieve CloudWatch logs for a workflow run
        
        Args:
            run_id: Run ID
            task_name: Specific task name (optional, shows all if None)
            limit: Maximum number of log events to retrieve
        
        Returns:
            List of log events
        """
        print("="*70)
        print(f"CloudWatch Logs for Run: {run_id}")
        if task_name:
            print(f"Task: {task_name}")
        print("="*70)
        
        log_group = '/aws/omics/WorkflowLog'
        log_stream = run_id
        
        try:
            response = self.logs.get_log_events(
                logGroupName=log_group,
                logStreamName=log_stream,
                limit=limit
            )
            
            events = response.get('events', [])
            
            if not events:
                print("\nNo logs found yet (run may still be starting)")
                return []
            
            print(f"\nShowing {len(events)} log entries:\n")
            
            for event in events:
                timestamp = datetime.fromtimestamp(event['timestamp']/1000)
                message = event['message']
                
                # Filter by task if specified
                if task_name and task_name.lower() not in message.lower():
                    continue
                
                print(f"[{timestamp}] {message}")
            
            return events
            
        except ClientError as e:
            if e.response['Error']['Code'] == 'ResourceNotFoundException':
                print("\nLog stream not found (run may not have started yet)")
            else:
                print(f"Error retrieving logs: {e}")
            return []
    
    # ========================================================================
    # RUN MANAGEMENT
    # ========================================================================
    
    def list_runs(self, status=None, max_results=20):
        """List workflow runs with optional status filter"""
        print("="*70)
        print("Workflow Runs")
        if status:
            print(f"Filter: {status}")
        print("="*70)
        
        try:
            params = {'maxResults': max_results}
            if status:
                params['status'] = status
            
            response = self.omics.list_runs(**params)
            runs = response.get('items', [])
            
            if not runs:
                print("\nNo runs found")
                return []
            
            print(f"\nFound {len(runs)} runs:\n")
            
            for run in runs:
                print(f"Run: {run['name']}")
                print(f"  ID: {run['id']}")
                print(f"  Status: {run.get('status', 'N/A')}")
                print(f"  Started: {run.get('startTime', 'N/A')}")
                if run.get('stopTime'):
                    print(f"  Stopped: {run['stopTime']}")
                print()
            
            return runs
            
        except ClientError as e:
            print(f"Error listing runs: {e}")
            raise
    
    def cancel_run(self, run_id):
        """Cancel a running workflow"""
        print(f"Cancelling run: {run_id}")
        
        try:
            self.omics.cancel_run(id=run_id)
            print(f"✓ Cancellation request submitted")
        except ClientError as e:
            print(f"✗ Error cancelling run: {e}")
            raise
    
    def list_workflows(self):
        """List all workflows"""
        print("="*70)
        print("Available Workflows")
        print("="*70)
        
        try:
            response = self.omics.list_workflows()
            workflows = response.get('items', [])
            
            if not workflows:
                print("\nNo workflows found")
                return []
            
            print(f"\nFound {len(workflows)} workflows:\n")
            
            for wf in workflows:
                print(f"Workflow: {wf['name']}")
                print(f"  ID: {wf['id']}")
                print(f"  Type: {wf.get('type', 'N/A')}")
                print(f"  Status: {wf.get('status', 'N/A')}")
                print(f"  Created: {wf.get('creationTime', 'N/A')}")
                print()
            
            return workflows
            
        except ClientError as e:
            print(f"Error listing workflows: {e}")
            raise


def demonstrate_complete_workflow():
    """
    Demonstrate the complete workflow lifecycle
    """
    print("="*70)
    print("AWS HealthOmics Workflows - Complete Implementation")
    print("="*70)
    print("""
This demonstration shows the complete workflow lifecycle:

STEP 1: Container Setup (ECR)
  • List available containers
  • Show public container examples
  
STEP 2: Workflow Creation
  • Create WDL workflow definition
  • Bundle workflow files
  • Upload to S3
  • Create workflow in HealthOmics
  
STEP 3: Workflow Execution
  • Start workflow run
  • Monitor progress
  • View tasks
  • Check logs
  • Retrieve results

Let's begin!
""")
    
    # Initialize manager
    manager = HealthOmicsWorkflowComplete()
    
    # Step 1: Containers
    print("\n" + "="*70)
    print("STEP 1: CONTAINER SETUP")
    print("="*70)
    manager.list_ecr_repositories()
    manager.get_public_container_examples()
    
    # Step 2: Workflow Creation
    print("\n" + "="*70)
    print("STEP 2: WORKFLOW CREATION")
    print("="*70)
    wdl_file, params_file = manager.create_simple_wdl_workflow()
    
    print("\nTo complete workflow creation, you need to:")
    print("  1. Create workflow bundle: bundle_file = manager.create_workflow_bundle(wdl_file, params_file)")
    print("  2. Upload to S3: s3_uri = manager.upload_workflow_to_s3(bundle_file, 'your-bucket')")
    print("  3. Create workflow: workflow_id = manager.create_workflow('MyWorkflow', s3_uri)")
    
    # Step 3: Workflow Execution
    print("\n" + "="*70)
    print("STEP 3: WORKFLOW EXECUTION")
    print("="*70)
    manager.list_workflows()
    manager.list_runs()
    
    print("\nTo start a workflow run:")
    print("""
    parameters = {
        'reference_fasta': 's3://bucket/reference.fa',
        'input_fastq_r1': 's3://bucket/sample_R1.fastq.gz',
        'input_fastq_r2': 's3://bucket/sample_R2.fastq.gz',
        'sample_name': 'SAMPLE001',
        'threads': 8
    }
    
    run_id = manager.start_workflow_run(
        workflow_id='your-workflow-id',
        role_arn='arn:aws:iam::123456789012:role/OmicsUnifiedJobRole',
        parameters=parameters,
        output_uri='s3://bucket/outputs/'
    )
    
    # Monitor the run
    manager.monitor_run(run_id)
    
    # View tasks
    manager.get_run_tasks(run_id)
    
    # Check logs
    manager.get_run_logs(run_id)
    """)
    
    print("\n" + "="*70)
    print("Workflow Implementation Complete!")
    print("="*70)
    print("""
You now have all the code to:
  ✓ Manage containers in ECR
  ✓ Create workflow definitions
  ✓ Upload workflows to HealthOmics
  ✓ Start and monitor workflow runs
  ✓ Inspect tasks and logs
  ✓ Manage workflow lifecycle

Next steps:
  1. Create your own workflow definition (WDL or Nextflow)
  2. Test with sample data
  3. Scale to production workloads
""")


if __name__ == '__main__':
    demonstrate_complete_workflow()